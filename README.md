---
description: |
    API documentation for modules: flares, flares.active_region, flares.data, flares.utils.

lang: en

classoption: oneside
geometry: margin=1in
papersize: a4

linkcolor: blue
links-as-notes: true
...


    
# Module `flares` {#id}

# flares_code

This is the start of my package


    
## Sub-modules

* [flares.active_region](#flares.active_region)
* [flares.data](#flares.data)
* [flares.utils](#flares.utils)






    
# Module `flares.active_region` {#id}







    
## Classes


    
### Class `ActiveRegion` {#id}




>     class ActiveRegion(
>         hnum: int,
>         date: str,
>         root: str
>     )


An Active Region is an entry point into physical features for 
parameterization, segmentation, and graph methods. If the specified
(hnum, date) active region contains a nan, self.valid turns to false
and subsequent calls to ActiveRegion methods may fail given the existence
of nans


Args
-----=
**```hnum```** :&ensp;<code>int</code>
:   The specified harpnumber - file must exist in root/magnetogram/sharp_{hnum} 


**```date```** :&ensp;<code>datetime</code>
:   The specified active region date and time - this date must exist in the specified harpnumber data folder 


**```root```** :&ensp;<code>asdkln</code>
:   The path to the data. Root must be a directory that holds both root/magnetogram and root/continuum. Inside both


of these subfolders, there must be a series of folders labeled sharp_{hnum} that contain the sequence of fits files for extraction







    
#### Methods


    
##### Method `B_grad_moments` {#id}




>     def B_grad_moments(
>         self,
>         mask
>     )




    
##### Method `Bh_grad_moments` {#id}




>     def Bh_grad_moments(
>         self,
>         mask
>     )




    
##### Method `Bz_grad_moments` {#id}




>     def Bz_grad_moments(
>         self,
>         mask
>     )




    
##### Method `J_moments` {#id}




>     def J_moments(
>         self,
>         mask
>     )




    
##### Method `Jh_moments` {#id}




>     def Jh_moments(
>         self,
>         mask
>     )




    
##### Method `assert_B` {#id}




>     def assert_B(
>         self
>     )




    
##### Method `assert_Bh` {#id}




>     def assert_Bh(
>         self
>     )




    
##### Method `assert_Bp` {#id}




>     def assert_Bp(
>         self
>     )




    
##### Method `assert_J` {#id}




>     def assert_J(
>         self
>     )




    
##### Method `assert_Jh` {#id}




>     def assert_Jh(
>         self
>     )




    
##### Method `assert_background` {#id}




>     def assert_background(
>         self
>     )




    
##### Method `assert_gamma` {#id}




>     def assert_gamma(
>         self
>     )




    
##### Method `assert_grad_B` {#id}




>     def assert_grad_B(
>         self
>     )




    
##### Method `assert_grad_Bh` {#id}




>     def assert_grad_Bh(
>         self
>     )




    
##### Method `assert_grad_Bm` {#id}




>     def assert_grad_Bm(
>         self
>     )




    
##### Method `assert_grad_Bx` {#id}




>     def assert_grad_Bx(
>         self
>     )




    
##### Method `assert_grad_By` {#id}




>     def assert_grad_By(
>         self
>     )




    
##### Method `assert_grad_Bz` {#id}




>     def assert_grad_Bz(
>         self
>     )




    
##### Method `assert_hc` {#id}




>     def assert_hc(
>         self
>     )




    
##### Method `assert_masks` {#id}




>     def assert_masks(
>         self
>     )




    
##### Method `assert_neutral_lines` {#id}




>     def assert_neutral_lines(
>         self,
>         kernel_radius=1
>     )




    
##### Method `assert_rho` {#id}




>     def assert_rho(
>         self
>     )




    
##### Method `assert_shear` {#id}




>     def assert_shear(
>         self
>     )




    
##### Method `assert_umbra_pumbra` {#id}




>     def assert_umbra_pumbra(
>         self
>     )




    
##### Method `come_back_from_gpu` {#id}




>     def come_back_from_gpu(
>         self
>     )




    
##### Method `connect_edges` {#id}




>     def connect_edges(
>         self
>     )




    
##### Method `gamma_moments` {#id}




>     def gamma_moments(
>         self,
>         mask
>     )




    
##### Method `group_pixels` {#id}




>     def group_pixels(
>         self,
>         mask
>     )




    
##### Method `h_moments` {#id}




>     def h_moments(
>         self,
>         mask
>     )




    
##### Method `hc_moments` {#id}




>     def hc_moments(
>         self,
>         mask
>     )




    
##### Method `hctot` {#id}




>     def hctot(
>         self,
>         mask
>     )




    
##### Method `hctotabs` {#id}




>     def hctotabs(
>         self,
>         mask
>     )




    
##### Method `ihtot` {#id}




>     def ihtot(
>         self,
>         mask
>     )




    
##### Method `ihtotabs` {#id}




>     def ihtotabs(
>         self,
>         mask
>     )




    
##### Method `itot` {#id}




>     def itot(
>         self,
>         mask
>     )




    
##### Method `itot_polarity` {#id}




>     def itot_polarity(
>         self,
>         mask
>     )




    
##### Method `itotabs` {#id}




>     def itotabs(
>         self,
>         mask
>     )




    
##### Method `largest_n_clusters` {#id}




>     def largest_n_clusters(
>         self,
>         labels,
>         sizes,
>         n=6
>     )




    
##### Method `moment` {#id}




>     def moment(
>         self,
>         data
>     )




    
##### Method `norm` {#id}




>     def norm(
>         self,
>         data
>     )




    
##### Method `phitot` {#id}




>     def phitot(
>         self,
>         mask
>     )




    
##### Method `phitotabs` {#id}




>     def phitotabs(
>         self,
>         mask
>     )




    
##### Method `physical_features` {#id}




>     def physical_features(
>         self,
>         mask
>     )


Extracts the physical fetures from the active region.


Args
-----=
**```mask```** :&ensp;<code>Numpy boolean Array</code> of <code>the same shape as self</code>
:   A mask (subset) of self to extract physical 


features on. ie, which pixels of self should this function compute

Returns
-----=
<code>1 dimensional Numpy array</code>
:   an array with all of the physical features computed on the subset provided by mask



    
##### Method `remove_bordering_pixels` {#id}




>     def remove_bordering_pixels(
>         self,
>         labeled,
>         labels,
>         sizes
>     )




    
##### Method `remove_percentage_max` {#id}




>     def remove_percentage_max(
>         self,
>         labeled,
>         labels,
>         sizes,
>         p=0.1
>     )




    
##### Method `remove_small_groups` {#id}




>     def remove_small_groups(
>         self,
>         labeled,
>         labels,
>         sizes,
>         p=500
>     )




    
##### Method `rho_moments` {#id}




>     def rho_moments(
>         self,
>         mask
>     )




    
##### Method `shear_moments` {#id}




>     def shear_moments(
>         self,
>         mask
>     )




    
##### Method `switch_to_gpu` {#id}




>     def switch_to_gpu(
>         self
>     )




    
##### Method `totrho` {#id}




>     def totrho(
>         self,
>         mask
>     )




    
##### Method `twist_moments` {#id}




>     def twist_moments(
>         self,
>         mask
>     )




    
##### Method `z_moments` {#id}




>     def z_moments(
>         self,
>         mask
>     )






    
# Module `flares.data` {#id}






    
## Functions


    
### Function `get_data` {#id}




>     def get_data(
>         harpnum,
>         date,
>         root,
>         nantozero=False
>     )




    
### Function `get_dates` {#id}




>     def get_dates(
>         harpnum,
>         root,
>         sort=False
>     )


$$\theta = 5x + \frac{1}{\exp{53}}$$




    
# Module `flares.utils` {#id}






    
## Functions


    
### Function `cov` {#id}




>     def cov(
>         m,
>         rowvar=False
>     )




    
### Function `fractal_dimension` {#id}




>     def fractal_dimension(
>         Z,
>         threshold=0.9
>     )




    
### Function `gradient` {#id}




>     def gradient(
>         data
>     )




    
### Function `norm` {#id}




>     def norm(
>         data
>     )






-----
Generated by *pdoc* 0.10.0 (<https://pdoc3.github.io>).
