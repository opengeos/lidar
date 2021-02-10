"""Module for applying filters to image.

"""
import os
import pkg_resources
import richdem as rd
from scipy import ndimage
import numpy as np
import time


def np2rdarray(in_array, no_data, projection, geotransform):
    """Converts an numpy array to rdarray.

    Args:
        in_array (np.array): The input numpy array.
        no_data (float): The no_data value of the array.
        projection (str): The projection of the image.
        geotransform (str): The geotransform of the image.

    Returns:
        object: The richDEM array.
    """
    out_array = rd.rdarray(in_array, no_data=no_data)
    out_array.projection = projection
    out_array.geotransform = geotransform
    return out_array


def MeanFilter(in_dem, kernel_size=3, out_file=None):
    """Applies a mean filter to an image.

    Args:
        in_dem (str): File path to the input image.
        kernel_size (int, optional): The size of the moving window. Defaults to 3.
        out_file (str, optional): File path to the output image. Defaults to None.

    Returns:
        np.array: The numpy array containing the filtered image.
    """
    print("Mean filtering ...")
    start_time = time.time()
    dem = rd.LoadGDAL(in_dem)
    no_data = dem.no_data
    projection = dem.projection
    geotransform = dem.geotransform

    weights = np.full((kernel_size, kernel_size), 1.0 / (kernel_size * kernel_size))
    mean = ndimage.filters.convolve(dem, weights)
    mean = np2rdarray(mean, no_data, projection, geotransform)
    print("Run time: {:.4f} seconds".format(time.time() - start_time))

    if out_file is not None:
        print("Saving dem ...")
        rd.SaveGDAL(out_file, mean)
        return out_file

    return mean


def MedianFilter(in_dem, kernel_size=3, out_file=None):
    """Applies a median filter to an image.

    Args:
        in_dem (str): File path to the input image.
        kernel_size (int, optional): The size of the moving window. Defaults to 3.
        out_file (str, optional): File path to the output image. Defaults to None.

    Returns:
        np.array: The numpy array containing the filtered image.
    """
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


def GaussianFilter(in_dem, sigma=1, out_file=None):
    """Applies a Gaussian filter to an image.

    Args:
        in_dem (str): File path to the input image.
        sigma (int, optional): Standard deviation. Defaults to 1.
        out_file (str, optional): File path to the output image. Defaults to None.

    Returns:
        np.array: The numpy array containing the filtered image.
    """
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
