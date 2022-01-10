import torch
import numpy as np


###### CONSTANTS
dev = torch.device("cpu")

###### UTILITIES
def norm(data):
    n = 0
    for i in data:
        n += np.abs(i)
    return n


radius = 10
dz = 0.0001
dist_kern = torch.zeros((2*radius + 1, 2*radius + 1))
for x0 in range(radius + 1):
    for y0 in range(radius + 1):
        if norm((x0, y0)) <= radius:
            v = dz / norm((x0, y0, dz))
            dist_kern[radius + x0][radius + y0] = v
            dist_kern[radius - x0][radius + y0] = v
            dist_kern[radius + x0][radius - y0] = v
            dist_kern[radius - x0][radius - y0] = v


def gradient(data):
    retrows = torch.zeros(data.shape, device = dev)
    retcols = torch.zeros(data.shape, device = dev)

    retrows[1:-1,:] = data[2:,:] - data[:-2,:]

    retrows[0] = (-3 * data[0] + 4*data[1] - data[2])
    retrows[-1] = -(-3 * data[-1] + 4*data[-2] - data[-3])

    retcols[:,1:-1] = data[:,2:] - data[:,:-2]
    retcols[:,0] = (-3 * data[:,0] + 4*data[:,1] - data[:,2])
    retcols[:,-1] = -(-3 * data[:,-1] + 4*data[:,-2] - data[:,-3])

    return 0.5 * retrows, 0.5 * retcols

def cov(m, rowvar=False):
    if m.dim() > 2:
        raise ValueError('m has more than 2 dimensions')
    if m.dim() < 2:
        m = m.view(1, -1)
    if not rowvar and m.size(0) != 1:
        m = m.t()
    # m = m.type(torch.double)  # uncomment this line if desired
    fact = 1.0 / (m.size(1) - 1)
    m -= torch.mean(m, dim=1, keepdim=True)
    mt = m.t()  # if complex: mt = m.t().conj()
    return fact * m.matmul(mt).squeeze()

def fractal_dimension(Z, threshold = 0.9):
    # Only for 2d image
    assert(len(Z.shape) == 2)

    # From https://github.com/rougier/numpy-100 (#87)
    def boxcount(Z, k):
        S = np.add.reduceat(
            np.add.reduceat(Z, np.arange(0, Z.shape[0], k), axis=0),
                               np.arange(0, Z.shape[1], k), axis=1)

        # We count non-empty (0) and non-full boxes (k*k)
        return len(np.where((S > 0) & (S < k*k))[0])


    # Transform Z into a binary array
    Z = (Z < threshold)

    # Minimal dimension of image
    p = min(Z.shape)

    # Greatest power of 2 less than or equal to p
    n = 2**np.floor(np.log(p)/np.log(2))

    # Extract the exponent
    n = int(np.log(n)/np.log(2))

    # Build successive box sizes (from 2**n down to 2**1)
    sizes = 2**np.arange(n, 1, -1)

    # Actual box counting with decreasing size
    counts = []
    for size in sizes:
        counts.append(boxcount(Z, size))

    # Fit the successive log(sizes) with log (counts)
    coeffs = np.polyfit(np.log(sizes), np.log(counts), 1)
    return -coeffs[0]