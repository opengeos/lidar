import os
import pkg_resources
import richdem as rd
from scipy import ndimage
import numpy as np
import time

# convert numpy array to rdarray
def np2rdarray(in_array, no_data, projection, geotransform):
    out_array = rd.rdarray(in_array, no_data=no_data)
    out_array.projection = projection
    out_array.geotransform = geotransform
    return out_array


# mean filter
def MeanFilter(in_dem, kernel_size=3, out_file=None):

    print("Mean filtering ...")
    start_time = time.time()
    dem = rd.LoadGDAL(in_dem)
    no_data = dem.no_data
    projection = dem.projection
    geotransform = dem.geotransform

    weights = np.full((kernel_size, kernel_size), 1.0/(kernel_size * kernel_size))
    mean = ndimage.filters.convolve(dem, weights)
    mean = np2rdarray(mean, no_data, projection, geotransform)
    print("Run time: {:.4f} seconds".format(time.time() - start_time))

    if out_file is not None:
        print("Saving dem ...")
        rd.SaveGDAL(out_file, mean)
        return out_file

    return mean


# median filter
def MedianFilter(in_dem, kernel_size=3, out_file=None):

    print("Median filtering ...")
    start_time = time.time()
    dem = rd.LoadGDAL(in_dem)
    no_data = dem.no_data
    projection = dem.projection
    geotransform = dem.geotransform

    med = ndimage.median_filter(dem, size=kernel_size)
    med = np2rdarray(med, no_data, projection, geotransform)
    print("Run time: {:.4f} seconds".format(time.time() - start_time))

    if out_file is not None:
        print("Saving dem ...")
        rd.SaveGDAL(out_file, med)
        return out_file

    return med


# Gaussian filter
def GaussianFilter(in_dem, sigma=1, out_file=None):

    print("Gaussian filtering ...")
    start_time = time.time()
    dem = rd.LoadGDAL(in_dem)
    no_data = dem.no_data
    projection = dem.projection
    geotransform = dem.geotransform

    gau = ndimage.gaussian_filter(dem, sigma=sigma)
    gau = np2rdarray(gau, no_data, projection, geotransform)
    print("Run time: {:.4f} seconds".format(time.time() - start_time))

    if out_file is not None:
        print("Saving dem ...")
        rd.SaveGDAL(out_file, gau)
        return out_file

    return gau


# #####################################  main script
if __name__ == '__main__':

    # identify the sample data directory of the package
    package_name = 'lidar'
    data_dir = pkg_resources.resource_filename(package_name, 'data/')
    print("Sample data directory: {}".format(data_dir))

    # use the sample dem. Change it to your own dem if needed
    in_dem = os.path.join(data_dir, 'dem.tif')
    out_dir = os.path.join(os.path.expanduser("~"), "temp")

    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    median_dem = os.path.join(out_dir, "median.tif")
    gaussian_dem = os.path.join(out_dir, "gaussian.tif")
    mean_dem = os.path.join(out_dir, "mean.tif")

    dem = rd.LoadGDAL(in_dem)                   # original dem

    mean = MeanFilter(in_dem, kernel_size=3)
    mean_diff = mean - dem
    med = MedianFilter(in_dem, kernel_size=3)
    med_diff = med - dem
    gau = GaussianFilter(in_dem, sigma=1)
    gau_diff = gau - dem

    # plotting data
    dem_fig = rd.rdShow(dem, ignore_colours=[0], axes=False, cmap='jet', figsize=(6, 5.5))
    mean_fig = rd.rdShow(mean, ignore_colours=[0], axes=False, cmap='jet', figsize=(6, 5.5))
    med_fig = rd.rdShow(med, ignore_colours=[0], axes=False, cmap='jet', figsize=(6, 5.5))
    gau_fig = rd.rdShow(gau, ignore_colours=[0], axes=False, cmap='jet', figsize=(6, 5.5))
    mean_diff_fig = rd.rdShow(mean_diff, ignore_colours=[0], axes=False, cmap='jet', figsize=(6, 5.5))
    med_diff_fig = rd.rdShow(med_diff, ignore_colours=[0], axes=False, cmap='jet', figsize=(6, 5.5))
    gau_diff_fig = rd.rdShow(gau_diff, ignore_colours=[0], axes=False, cmap='jet', figsize=(6, 5.5))

    # save results
    rd.SaveGDAL(mean_dem, mean)
    rd.SaveGDAL(median_dem, med)
    rd.SaveGDAL(gaussian_dem, gau)
