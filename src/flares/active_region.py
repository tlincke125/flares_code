from pathlib import Path
import sys
path_root = Path(__file__).parents[2]
sys.path.append(str(path_root))

from flares.utils import *
from flares.data import *

import torch
import torch.nn.functional as F
from skimage.morphology import square, binary_dilation
from skimage.measure import label
import numpy as np
import networkx as nx
import warnings

from scipy.stats import skew
from skimage.filters import threshold_local
from skimage.morphology import square
from skimage.measure import label


NUM_FEATURES = 55
GRAPH_EMB_SIZE = 128


# Feature Library
class ActiveRegion:
    def __init__(self, hnum: int, date: str, root: str):
        """An Active Region is an entry point into physical features for 
        parameterization, segmentation, and graph methods. If the specified
        (hnum, date) active region contains a nan, self.valid turns to false
        and subsequent calls to ActiveRegion methods may fail given the existence
        of nans

        Args:
            hnum (int): The specified harpnumber - file must exist in root/magnetogram/sharp_{hnum} 
            date (datetime): The specified active region date and time - this date must exist in the specified harpnumber data folder 
            root (asdkln): The path to the data. Root must be a directory that holds both root/magnetogram and root/continuum. Inside both
            of these subfolders, there must be a series of folders labeled sharp_{hnum} that contain the sequence of fits files for extraction
        """
    
        # Generate xyz components of magnetic field and continuum
        self.Bz, self.Bx, self.By, self.cont = get_data(hnum, date, root)
        self.shape = self.Bz.shape
        self.valid = True # Valid is false

        if np.count_nonzero(np.isnan(self.Bz)) / self.Bz.size > 0.0:
            self.valid = False
            warnings.warn(f"Hnum {hnum} date {date} has a nan, skipping")
            return

        
        # To prevent nan's popping up do to small numbers,
        # cut everything below minimum val to minimum val
        minimum_val = 0.01
        self.Bz[np.abs(self.Bz) < minimum_val] = minimum_val
        self.Bx[np.abs(self.Bx) < minimum_val] = minimum_val
        self.By[np.abs(self.By) < minimum_val] = minimum_val


        # All the fields / features
        self.Bh = None
        self.gamma = None
        self.B = None
        self.grad_B_x = None
        self.grad_B_y = None
        self.grad_Bh_x = None
        self.grad_Bh_y = None
        self.grad_Bz_x = None
        self.grad_Bz_y = None
        self.grad_Bx_x = None
        self.grad_Bx_y = None
        self.grad_By_x = None
        self.grad_By_y = None
        self.grad_Bm = None
        self.J = None
        self.Jh = None
        self.hc = None
        self.shear = None
        self.rho = None
        self.Bpx = None
        self.Bpy = None
        self.Bpz = None

        # The three datasets
        self.baseline = self.physical_features(np.ones(self.shape, dtype = bool))
        self.segmented = np.zeros(4 * NUM_FEATURES)
        self.graph_ds = np.zeros(GRAPH_EMB_SIZE)

        # List of node masks
        # Each index of node_masks corresponds to its feature vector in feature_vec
        # Ditto w/ G
        self.node_masks = np.zeros((0, self.shape[0], self.shape[1]))
        self.node_feature_vecs = np.zeros((0, NUM_FEATURES))
        self.G = nx.Graph() 

        self.umbra = None
        self.pumbra = None
        self.nl = None
        self.background = None
        self.graph_edges_declared = False


    # Main Two Functions for accessing data
    def physical_features(self, mask):
        """Extracts the physical fetures from the active region.

        Args:
            mask (Numpy boolean Array of the same shape as self): A mask (subset) of self to extract physical 
            features on. ie, which pixels of self should this function compute

        Returns:
            1 dimensional Numpy array: an array with all of the physical features computed on the subset provided by mask
        """
        mask = torch.from_numpy(mask).to(dev)
        self.switch_to_gpu()

        moment_features = torch.tensor([
            *self.z_moments(mask), \
            *self.h_moments(mask), \
            *self.gamma_moments(mask), \
            *self.B_grad_moments(mask), \
            *self.Bz_grad_moments(mask), \
            *self.Bh_grad_moments(mask), \
            *self.J_moments(mask), \
            *self.Jh_moments(mask), \
            *self.twist_moments(mask), \
            *self.hc_moments(mask), \
            *self.shear_moments(mask), \
            *self.rho_moments(mask)
            ])
        val = torch.tensor([
            self.phitot(mask), \
            self.phitotabs(mask), \
            self.itot(mask), \
            self.itotabs(mask), \
            self.hctot(mask), \
            self.hctotabs(mask), \
            self.totrho(mask), \
            *moment_features
            ])
        self.come_back_from_gpu()
        mask = mask.detach().cpu().numpy()
        return val

    def switch_to_gpu(self):
        self.Bz = torch.from_numpy(self.Bz).to(dev)
        self.Bx = torch.from_numpy(self.Bx).to(dev)
        self.By = torch.from_numpy(self.By).to(dev)
        self.cont = torch.from_numpy(self.cont).to(dev)

    def come_back_from_gpu(self):
        self.Bz = self.Bz.detach().cpu().numpy()
        self.Bx = self.Bx.detach().cpu().numpy()
        self.By = self.By.detach().cpu().numpy()
        self.cont = self.cont.detach().cpu().numpy()

    def assert_masks(self):
        self.assert_neutral_lines()
        self.assert_umbra_pumbra()
        self.assert_background()

    def assert_neutral_lines(self, kernel_radius = 1):
        if self.nl is None:
            # Find neutral Lines
            thresh = 150
            nl_mask = binary_dilation(self.Bz < -thresh, square(3)) & binary_dilation(self.Bz > thresh, square(3))


            labeled, labels, sizes = self.group_pixels(nl_mask)
            labels, sizes = self.remove_small_groups(labeled, labels, sizes, 10)
            labels, sizes = self.remove_percentage_max(labeled, labels, sizes)
            labels, sizes = self.largest_n_clusters(labels, sizes)

            # Add all the graph nodes
            nl_mask = np.zeros(self.shape, dtype = bool)
            cur_node = len(self.node_masks)

            if len(sizes) == 0:
                self.segmented[0:NUM_FEATURES] = 0
                self.nl = nl_mask
                return

            for i in labels:
                mask = labeled == i
                self.node_masks = np.concatenate((self.node_masks, mask[None,...]), axis = 0)
                self.node_feature_vecs = np.vstack((self.node_feature_vecs, self.physical_features(mask)))
                self.G.add_node(cur_node)
                cur_node += 1
                nl_mask |= mask

            # Compute the segmented data set
            self.segmented[0:NUM_FEATURES] = self.physical_features(nl_mask)
            self.nl = nl_mask


    def assert_umbra_pumbra(self):
        if self.pumbra is None or self.umbra is None:
            cont_bounded = (255 * (self.cont - np.min(self.cont)) / np.ptp(self.cont)).astype(np.uint8)
 
            block_size = np.min(self.shape)
            if block_size % 2 == 0:
                block_size -= 1
                
            offset = 10
            binary_adaptive = cont_bounded < (threshold_local(cont_bounded, block_size, offset = offset) - offset)

            labeled_0, labels, sizes = self.group_pixels(binary_adaptive)
            labels, sizes = self.remove_small_groups(labeled_0, labels, sizes)
            labels, sizes = self.remove_bordering_pixels(labeled_0, labels, sizes)
            labels, sizes = self.remove_percentage_max(labeled_0, labels, sizes)
            labels, sizes = self.largest_n_clusters(labels, sizes)

            um_mask = np.zeros(self.shape, dtype = bool)
            pu_mask = np.zeros(self.shape, dtype = bool)
            cur_node = 0

            if len(sizes) == 0:
                self.segmented[NUM_FEATURES:2*NUM_FEATURES] = 0
                self.segmented[2*NUM_FEATURES:3*NUM_FEATURES] = 0
                self.umbra = um_mask
                self.pumbra = pu_mask
                return

            for i in labels:
                mask = labeled_0 == i
                mx = np.max(self.cont[mask])
                mn = np.min(self.cont[mask])
                t = (mx - mn) / 2 + mn
                if mx - mn > 21000:
                    # Both umbra and penumbra
                    um = mask & (self.cont <= t)
                    pu = mask & (self.cont > t)


                    # Further segment the umbra node again
                    labeled, labels, sizes = self.group_pixels(pu)
                    labels, sizes = self.remove_small_groups(labeled, labels, sizes, 10)
                    labels, sizes = self.remove_percentage_max(labeled, labels, sizes)
                    labels, sizes = self.largest_n_clusters(labels, sizes)

                    for i in labels:
                        mask = labeled == i
                        self.node_masks = np.concatenate((self.node_masks, mask[None,...]), axis = 0)
                        self.node_feature_vecs = np.vstack((self.node_feature_vecs, self.physical_features(mask)))
                        self.G.add_node(cur_node)
                        cur_node += 1
                        pu_mask |= mask

                    # Further segment the umbra node again
                    labeled, labels, sizes = self.group_pixels(um)
                    labels, sizes = self.remove_small_groups(labeled, labels, sizes, 10)
                    labels, sizes = self.remove_percentage_max(labeled, labels, sizes)
                    labels, sizes = self.largest_n_clusters(labels, sizes)

                    for i in labels:
                        mask = labeled == i
                        self.node_masks = np.concatenate((self.node_masks, mask[None,...]), axis = 0)
                        self.node_feature_vecs = np.vstack((self.node_feature_vecs, self.physical_features(mask)))
                        self.G.add_node(cur_node)
                        cur_node += 1
                        um_mask |= mask

                else:
                    # Only umbra
                    um = mask & (self.cont <= t)
                    self.node_masks = np.concatenate((self.node_masks, um[None,...]), axis = 0)
                    self.node_feature_vecs = np.vstack((self.node_feature_vecs, self.physical_features(um)))
                    self.G.add_node(cur_node)
                    cur_node += 1
                    um_mask |= um
            
            # Compute the segmented data set
            self.segmented[NUM_FEATURES:2*NUM_FEATURES] = self.physical_features(um_mask)
            self.segmented[2*NUM_FEATURES:3*NUM_FEATURES] = self.physical_features(pu_mask)
            self.umbra = um_mask
            self.pumbra = pu_mask

    # Grouping and filtering operations
    def group_pixels(self, mask):
        labeled = label(mask, connectivity = 2)
        labels = np.unique(labeled)[1:]
        sizes = np.array([np.count_nonzero(labeled==i) for i in labels])
        return labeled, labels, sizes
    def remove_small_groups(self, labeled, labels, sizes, p = 500):
        if len(sizes) == 0:
            return labels, sizes
        filt = np.argwhere((sizes < p))
        return np.delete(labels, filt), np.delete(sizes, filt)
    def remove_bordering_pixels(self, labeled, labels, sizes):
        if len(sizes) == 0:
            return labels, sizes
        bordered = []
        for i in range(len(labels)):
            rows, cols = np.where(labeled == labels[i])
            if min(rows) == 0 or min(cols) == 0:
                bordered.append(i)
            if max(cols) == self.shape[1] - 1 or max(rows) == self.shape[0] - 1:
                bordered.append(i)
        return np.delete(labels, bordered), np.delete(sizes, bordered)
    def remove_percentage_max(self, labeled, labels, sizes, p = 0.1):
        if len(sizes) == 0:
            return labels, sizes
        filt = np.argwhere(sizes < p * np.max(sizes))
        return np.delete(labels, filt), np.delete(sizes, filt)
    def largest_n_clusters(self, labels, sizes, n = 6):
        if len(sizes) == 0:
            return labels, sizes
        n = min(n, len(labels))
        a = np.partition(sizes, -n)[-n]
        return labels[sizes >= a], sizes[sizes >= a]

    def assert_background(self):
        if self.background is None:
            self.background = torch.zeros(self.shape, dtype=bool)
            self.assert_neutral_lines()
            self.assert_umbra_pumbra()

            #  node = {"COM" : (0, 0), "MASK" : ~(self.nl | self.umbra | self.pumbra)}
            background = ~(self.nl | self.umbra | self.pumbra)
            self.segmented[3*NUM_FEATURES:4*NUM_FEATURES] = self.physical_features(background)
            self.background = background

    def connect_edges(self):
        self.assert_masks()
        if not self.graph_edges_declared:
            for i in range(len(self.node_masks) - 1):
                for j in range(i + 1, len(self.node_masks)):
                    if np.count_nonzero(binary_dilation(self.node_masks[i], square(3)) & binary_dilation(self.node_masks[j], square(3))) > 0:
                        self.G.add_edge(i, j)
            self.graph_edges_declared = True
            

    # Tools
    def moment(self, data):
        avg = torch.mean(data)
        std = torch.sqrt(torch.mean((data - avg)**2))
        skw = torch.mean(((data - avg)/std)**3)
        krt = torch.mean(((data - avg)/std)**4) - 3.0
        return torch.tensor([avg, std, skw, krt])

    def norm(self, data):
        n = 0
        for i in data:
            n += torch.abs(i)
        return n


    # Physical Properties
    def z_moments(self, mask):
        return self.moment(self.Bz[mask])
    def phitot(self, mask):
        return torch.sum(torch.abs(self.Bz[mask]))
    def phitotabs(self, mask):
        return torch.abs(torch.sum(self.Bz[mask]))
    def h_moments(self, mask):
        self.assert_Bh()
        return self.moment(self.Bh[mask])
    def gamma_moments(self, mask):
        self.assert_gamma()
        return self.moment(self.gamma[mask])
    def B_grad_moments(self, mask):
        self.assert_grad_B()
        return self.moment(self.norm((self.grad_B_x[mask], self.grad_B_y[mask])))
    def Bz_grad_moments(self, mask):
        self.assert_grad_Bz()
        return self.moment(self.norm((self.grad_Bz_x[mask], self.grad_Bz_y[mask])))
    def Bh_grad_moments(self, mask):
        self.assert_grad_Bh()
        return self.moment(self.norm((self.grad_Bh_x[mask], self.grad_Bh_y[mask])))
    def J_moments(self, mask):
        self.assert_J()
        return self.moment(self.J[mask])
    def itot(self, mask):
        self.assert_J()
        return torch.sum(torch.abs(self.J[mask]))
    def itotabs(self, mask):
        self.assert_J()
        return torch.abs(torch.sum(self.J[mask]))
    def itot_polarity(self, mask):
        self.assert_J()
        return torch.abs(torch.sum(self.J[self.Bz > 0 & mask])) + torch.abs(torch.sum(self.J[self.Bz < 0 & mask]))
    def Jh_moments(self, mask):
        self.assert_Jh()
        return self.moment(self.Jh[mask])
    def ihtot(self, mask):
        self.assert_Jh()
        return torch.sum(torch.abs(self.Jh[mask]))
    def ihtotabs(self, mask):
        self.assert_Jh()
        return torch.abs(torch.sum(self.Jh[mask]))
    def twist_moments(self, mask):
        self.assert_J()
        return self.moment(self.J[mask] / self.Bz[mask])
    # Too computationally intensive
    #def twist_force_free_m(self):
    #    pass
    def hc_moments(self, mask):
        self.assert_hc()
        return self.moment(self.hc[mask])
    def hctot(self, mask):
        self.assert_hc()
        return torch.sum(torch.abs(self.hc[mask]))
    def hctotabs(self, mask):
        self.assert_hc()
        return torch.abs(torch.sum(self.hc[mask]))
    def shear_moments(self, mask):
        self.assert_shear() # Pretty big function call right here
        return self.moment(self.shear[mask])
    def rho_moments(self, mask):
        self.assert_rho()
        return self.moment(self.rho[mask])
    def totrho(self, mask):
        self.assert_rho()
        return torch.sum(self.rho[mask])


    # ASSERTIONS - designed so that each component is only calculated ONCE - example: if I want Bh,
    # simply assert_Bh. If Bh already exists, then nothing is done
    def assert_Bh(self):
        if self.Bh is None:
            self.Bh = self.norm((self.Bx, self.By))
    def assert_gamma(self):
        if self.gamma is None:
            self.assert_Bh()
            self.gamma = torch.arctan(self.Bz / self.Bh)
    def assert_B(self):
        if self.B is None:
            self.B = self.norm((self.Bx, self.By, self.Bz))
    def assert_grad_B(self):
        if self.grad_B_x is None or self.grad_B_y is None:
            self.assert_B()
            self.grad_B_x, self.grad_B_y = gradient(self.B)
    def assert_grad_Bh(self):
        if self.grad_Bh_x is None or self.grad_Bh_y is None:
            self.assert_Bh()
            self.grad_Bh_x, self.grad_Bh_y = gradient(self.Bh)
    def assert_grad_Bz(self):
        if self.grad_Bz_x is None or self.grad_Bz_y is None:
            self.grad_Bz_x, self.grad_Bz_y = gradient(self.Bz)
    def assert_grad_Bx(self):
        if self.grad_Bx_x is None or self.grad_Bx_y is None:
            self.grad_Bx_x, self.grad_Bx_y = gradient(self.Bx)
            self.grad_Bx_x = -self.grad_Bx_x
    def assert_grad_By(self):
        if self.grad_By_x is None or self.grad_By_y is None:
            self.grad_By_x, self.grad_By_y = gradient(self.By)
            self.grad_By_x = -self.grad_By_x
    def assert_grad_Bm(self):
        if self.grad_Bm is None:
            self.assert_grad_B()
            self.grad_Bm = self.norm((self.grad_B_x, self.grad_B_y))
    def assert_J(self):
        if self.J is None:
            self.assert_grad_Bx()
            self.assert_grad_By()
            self.J = self.grad_By_x - self.grad_Bx_y
    def assert_Jh(self):
        if self.Jh is None:
            self.assert_grad_Bx()
            self.assert_grad_By()
            self.assert_B()
            self.Jh = (self.By * self.grad_Bx_y - self.Bx * self.grad_By_x) / self.B
    def assert_hc(self):
        if self.hc is None:
            self.assert_J()
            self.hc = self.Bz * self.J
    def assert_shear(self):
        if self.shear is None:
            self.assert_Bp()
            self.assert_B()
            dot = self.Bx * self.Bpx + self.By * self.Bpy + self.Bz * self.Bpz
            magp = self.norm((self.Bpx, self.Bpy, self.Bpz))
            self.shear = torch.arccos(dot / (self.B * magp))
    def assert_rho(self):
        if self.rho is None:
            self.assert_Bp()
            self.rho = (self.B - self.norm((self.Bpx, self.Bpy, self.Bpz)))**2
    def assert_Bp(self):
        if self.Bpx is None or self.Bpy is None or self.Bpz is None:

            Bz = F.pad(self.Bz, (radius, radius, radius, radius)).float()

            dz = 0.0001
            Xmx = Bz.shape[0]
            Ymx = Bz.shape[1]

            pot = torch.zeros(self.shape, device = dev).float()

            # Distance kernel - a kernel with values filled in the "circle" (by def of norm) as the distance from
            # the center multiplied by dz (for integration)

            kern = dist_kern.to(dev)

            # Convolution -- integrate over each pixel
            pot = F.conv2d(Bz[None, None, ...], kern[None, None, ...])[0][0]

            # Save potential
            self.potential = pot

            # Get Potential Fields
            self.Bpz = self.Bz
            grad = gradient(self.potential)
            self.Bpx, self.Bpy = -grad[0], -grad[1]
