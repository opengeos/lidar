from __future__ import print_function
from __future__ import division
from __future__ import print_function
from __future__ import division
import sys
from skimage.external.tifffile import TiffFile
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import os
from scipy import ndimage
import arcpy
import matplotlib as mpl
from matplotlib.colors import rgb_to_hsv, hsv_to_rgb

DEF_AZIMUTH = 135  # degrees
DEF_ELEVATION = 45  # degrees
DEF_AMBIENT_WEIGHT = 1
DEF_LAMP_WEIGHT = 5
DO_SANITY_CHECKS = (
    True  # If True intermediate results will be checked for boundary values.
)
INTENSITY_CMAP = plt.cm.get_cmap("gray")
INTENSITY_CMAP.set_bad("red")
INTENSITY_CMAP.set_over("blue")  # to check that no intensity is above 1
INTENSITY_CMAP.set_under("yellow")  # to check that no intensity is below 0
DEF_CMAP = plt.cm.get_cmap("gist_earth")


def weighted_intensity(
    terrain,
    azimuth=DEF_AZIMUTH,
    elevation=DEF_ELEVATION,
    ambient_weight=DEF_AMBIENT_WEIGHT,
    lamp_weight=DEF_LAMP_WEIGHT,
):
    """Calculates weighted average of the ambient illumination and the that of one or more lamps.

    The azimuth and elevation parameters can be scalars or lists. Use the latter for multiple
    lamps. They should be of equal length.

    The lamp_weight can be given per lamp or one value can be specified, which is then used for
    all lamps sources.

    See also the hill_shade doc string.
    """
    # Make sure input is in the correct shape
    azimuths = enforce_list(azimuth)
    elevations = enforce_list(elevation)
    assert_same_length(azimuths, elevations, "azimuths", "elevations")

    lamp_weights = enforce_list(lamp_weight)
    if len(lamp_weights) == 1:
        lamp_weights = lamp_weights * len(azimuths)
    assert_same_length(azimuths, lamp_weights, "azimuths", "lamp_weights")

    # Create weights and rel_intensities arrays
    rel_intensities = [np.ones_like(terrain)]
    weights = [ambient_weight]
    for azim, elev, lmpw in zip(azimuths, elevations, lamp_weights):
        rel_int = relative_surface_intensity(terrain, azimuth=azim, elevation=elev)
        rel_intensities.append(rel_int)
        weights.append(lmpw)

    rel_intensities = np.dstack(rel_intensities)
    weights = np.array(weights)

    # The actual weighted-average calculation
    unit_weights = weights / np.sum(weights)
    surface_intensity = np.average(rel_intensities, axis=2, weights=unit_weights)
    return surface_intensity


def relative_surface_intensity(terrain, azimuth=DEF_AZIMUTH, elevation=DEF_ELEVATION):
    """Calculates the intensity that falls on the surface for light of intensity 1.
    This equals cosine(theta) where theta is the angle between the direction of the light
    source and the surface normal. When the cosine is negative, the angle is > 90 degrees.
    In that case the surface receives no light so we clip to 0. The result of this function is
    therefore always between 0 and 1.
    """
    # cosine(theta) is the dot-product of the normal vector and the vector that contains the
    # direction of the light source. Both vectors must be unit vectors (have length 1).
    normals = surface_unit_normals(terrain)
    light = polar_to_cart3d(azimuth, elevation)
    intensity = np.dot(normals, light)

    if DO_SANITY_CHECKS:
        np.testing.assert_approx_equal(
            np.linalg.norm(light),
            1.0,
            err_msg="sanity check: light vector should have length 1",
        )
        assert np.all(intensity >= -1.0), "sanity check: cos(theta) should be >= -1"
        assert np.all(intensity <= 1.0), "sanity check: cos(theta) should be <= 1"

    # Where the dot product is smaller than 0 the angle between the light source and the surface
    # is larger than 90 degrees. These pixels receive no light so we clip the intensity to 0.
    intensity = np.clip(intensity, 0.0, 1.0)
    return intensity


def surface_unit_normals(terrain):
    """Returns an array of shape (n_rows, n_cols, 3) with unit surface normals.
    That is, each result[r,c,:] contains the vector of length 1, perpendicular to the surface.
    """
    dr, dc = np.gradient(terrain)

    # Vectors that do a step of 1 in the row direction, 0 in the column direction and dr upwards
    vr = np.dstack(
        (dr, np.ones_like(dr), np.zeros_like(dr))
    )  # shape = (n_rows, n_cols, 3)

    # Vectors that do a step of 0 in the row direction, 1 in the column direction and dc upwards
    vc = np.dstack(
        (dc, np.zeros_like(dc), np.ones_like(dc))
    )  # shape = (n_rows, n_cols, 3)

    # The surface normals can be calculated as the cross product of those vector pairs
    surface_normals = np.cross(vr, vc)  # surface_normals.shape = (n_rows, n_cols, 3)

    # Divide the normals by their magnitude to get unit vectors.
    # (Add artificial dimension of length 1 so that we can use broadcasting)
    normal_magnitudes = np.linalg.norm(surface_normals, axis=2)
    return surface_normals / np.expand_dims(normal_magnitudes, axis=2)


def polar_to_cart3d(azimuth, elevation):
    """Converts the polar (azimuth, elevation) unit vector to (height, row, col) coordinates."""
    azimuth_rad = azimuth * np.pi / 180.0
    elevation_rad = elevation * np.pi / 180.0

    height = np.sin(elevation_rad)
    row = np.cos(elevation_rad) * np.sin(azimuth_rad)
    col = np.cos(elevation_rad) * np.cos(azimuth_rad)

    return np.array([height, row, col])


def assert_same_length(s0, s1, label0, label1):
    """Asserts list s1 and s2 have the same lengths"""
    if len(s0) != len(s1):
        raise AssertionError(
            "size mismatch between {} (len={}) and {} (len={})".format(
                label0, s0, label1, s1
            )
        )


def enforce_list(var):
    """Runs the list() constructor on the var parameter."""
    try:
        return list(var)  # iterable
    except TypeError:
        return [var]


# The mpl_surface_intensity is included to compare the matplotlib implementation with
# the relative_surface_intensity() results.
#
def mpl_surface_intensity(
    terrain, azimuth=165, elevation=DEF_ELEVATION, azim0_is_east=False, normalize=False
):
    """Calculates the intensity that falls on the surface when illuminated with intensity 1

    This is the implementation as is used in matplotlib.
    Forked from Ran Novitsky's blog (no license found).
    The original source is the LightSource.shade_rgb function of the matplotlib.colors module.
    See:
        http://rnovitsky.blogspot.nl/2010/04/using-hillshade-image-as-intensity.html
        https://github.com/matplotlib/matplotlib/blob/master/lib/matplotlib/colors.py

    input:
        terrain - a 2-d array of the terrain
        azimuth - where the light comes from: 0 south ; 90 east ; 180 north ;
                    270 west
        elevation - where the light comes from: 0 horizon ; 90 zenith

    output:
        a 2-d array of normalized hillshade
    """
    from numpy import pi, cos, sin, gradient, arctan, hypot, arctan2

    # Convert alt, az to radians
    az = azimuth * pi / 180.0
    alt = elevation * pi / 180.0

    # gradient in x and y directions
    dx, dy = gradient(terrain)

    slope = 0.5 * pi - arctan(hypot(dx, dy))
    if azim0_is_east:
        # The arctan docs specify that the parameters are (y, x), in that order.
        # This makes an azimuth of 0 correspond to east.
        aspect = arctan2(dy, dx)
    else:
        aspect = arctan2(dx, dy)
    intensity = sin(alt) * sin(slope) + cos(alt) * cos(slope) * cos(
        -az - aspect - 0.5 * pi
    )

    if DO_SANITY_CHECKS:
        assert np.all(intensity >= -1.0), "sanity check: cos(theta) should be >= -1"
        assert np.all(intensity <= 1.0), "sanity check: cos(theta) should be <= 1"

    # The matplotlib source just normalizes the intensities. However, I believe that their
    # intensities are the same as mine so that, where they are < 0 the angle between the light
    # source and the surface is larger than 90 degrees. These pixels receive no light so
    # they should be clipped. This is done when the normalize parameter is set to False.

    if normalize:
        intensity = (intensity - intensity.min()) / (intensity.max() - intensity.min())
    else:
        intensity = np.clip(intensity, 0.0, 1.0)

    return intensity


def is_non_finite_mask(array):
    "Returns mask with ones where the data is infite or Nan"
    np.logical_not(np.isfinite(array))


def replace_nans(array, array_nan_value, mask=None):
    """Returns a copy of the array with the NaNs replaced by nan_value"""
    finite_array = np.copy(array)
    if mask is None:
        mask = is_non_finite_mask(array)
    if np.any(mask):
        finite_array[mask] = array_nan_value
    return finite_array


def normalize(values, vmin=None, vmax=None, norm=None):
    """Normalize values between using a mpl.colors.Normalize object or (vmin, vmax) interval.
    If norm is specified, vmin and vmax are ignored.
    If norm is None and vmin and vmax are None, the values are autoscaled.
    """
    if norm is None:
        norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)

    return norm(values)


def color_data(data, cmap, vmin=None, vmax=None, norm=None):
    """Auxiliary function that colors the data."""
    norm_data = normalize(data, vmin=vmin, vmax=vmax, norm=norm)
    rgba = cmap(norm_data)
    return rgba


def no_blending(rgba, norm_intensities):
    """Just returns the intensities. Use in hill_shade to just view the calculated intensities"""
    assert norm_intensities.ndim == 2, "norm_intensities must be 2 dimensional"
    return norm_intensities


def rgb_blending(rgba, norm_intensities):
    """Calculates image colors by multiplying the rgb value with the normalized intensities

    :param rgba: [nrows, ncols, 3|4] RGB or RGBA array. The alpha layer will be ignored.
    :param norm_intensities: normalized intensities

    Returns 3D array that can be plotted with matplotlib.imshow(). The last dimension is RGB.
    """
    assert rgba.ndim == 3, "rgb must be 3 dimensional"
    assert norm_intensities.ndim == 2, "norm_intensities must be 2 dimensional"

    # Add artificial dimension of length 1 at the end of norm_intensities so that it can be
    # multiplied with the rgb array using numpy broad casting
    expanded_intensities = np.expand_dims(norm_intensities, axis=2)
    rgb = rgba[:, :, :3]

    return rgb * expanded_intensities


def hsv_blending(rgba, norm_intensities):
    """Calculates image colors by placing the normalized intensities in the Value layer of the
    HSV color of the normalized data.

    IMPORTANT: may give incorrect results for color maps that include colors close to black
        (e.g. cubehelix or hot).

    :param rgba: [nrows, ncols, 3|4] RGB or RGBA array. The alpha layer will be ignored.
    :param norm_intensities: normalized intensities

    Returns 3D array that can be plotted with matplotlib.imshow(). The last dimension is RGB.
    """
    rgb = rgba[:, :, :3]
    hsv = rgb_to_hsv(rgb)
    hsv[:, :, 2] = norm_intensities
    return hsv_to_rgb(hsv)


def pegtop_blending(rgba, norm_intensities):
    """Calculates image colors with the Pegtop Light shading of ImageMagick

    See:
        http://www.imagemagick.org/Usage/compose/#pegtoplight

    Forked from Ran Novitsky's blog (no license found)
        http://rnovitsky.blogspot.nl/2010/04/using-hillshade-image-as-intensity.html

    :param rgba: [nrows, ncols, 3|4] RGB or RGBA array. The alpha layer will be ignored.
    :param norm_intensities: normalized intensities

    Returns 3D array that can be plotted with matplotlib.imshow(). The last dimension is RGB.
    """
    # get rgb of normalized data based on cmap
    rgb = rgba[:, :, :3]

    # form an rgb eqvivalent of intensity
    d = norm_intensities.repeat(3).reshape(rgb.shape)

    # simulate illumination based on pegtop algorithm.
    return 2 * d * rgb + (rgb ** 2) * (1 - 2 * d)


def hill_shade(
    data,
    terrain=None,
    azimuth=DEF_AZIMUTH,
    elevation=DEF_ELEVATION,
    ambient_weight=DEF_AMBIENT_WEIGHT,
    lamp_weight=DEF_LAMP_WEIGHT,
    cmap=DEF_CMAP,
    vmin=None,
    vmax=None,
    norm=None,
    blend_function=rgb_blending,
):
    """Calculates a shaded relief given a 2D array of surface heights.

    You can specify data properties and terrain height in separate parameters. The data array
    determines the (unshaded) color, the terrain is used to calculate the shading component.
    If the terrain is left to None, the data array will be used as terrain as well. The terrain
    parameter can also be used to scale the surface heights. E.g. use: terrain = data * 10

    The relief is calculated from one or more artificial light sources which positions are
    specified by their azimuth and elevation angle. Each lamp can have a weight, which
    corresponds to its strength. Ambient light is light that reaches the surface via indirect
    illumination. If the ambient weight set to 0, pixels that are not illuminated by any lamp
    are rendered completely black, which is usually undesirable.

    The algorithm uses a color map to color the relief. The minimum and maximum values of the
    color scale can be specified by vmin and vmax (or by giving a normalization function). If
    these are all None (the default), the color bar will be auto-scaled.

    The blend_function is the function that merges the color and shade components into the
    final result. It was found that rbg_blending (the default) gives the best results. If set
    to no_blending, only the intensities of the shade component are returned. This is useful
    for debugging.

    :param data: 2D array with terrain properties
    :param terrain: 2D array with terrain heights
    :param azimuth: azimuth angle [degrees] of the lamp direction(s). Can be scalar or list.
    :param elevation: elevation angle [degrees] of the lamp direction(s). Can be scalar or list.
    :param ambient_weight: the relative strength of the ambient illumination (default = 1)
    :param lamp_weight: the relative strength of the lamp or lamps (default = 5)
    :param cmap: matplotlib color map to color the data (default: 'gist_earth')
    :param vmin: use to set a minimum value of the color scale
    :param vmax: use to set a maximum value of the color scale
    :param norm: colorbar normalization function. E.g.: mpl.colors.Normalize(vmin=0.0, vmax=1.0)
    :param blend_function: function that blends shading and color (default = rbg_blending)

    :returns: 3D array (n_rows, n_cols, 3) with for each pixel an RGB color.
        If blend_function=no_blending the result is a 2D array with only shading intensities.
    """
    if terrain is None:
        terrain = data

    assert data.ndim == 2, "data must be 2 dimensional"
    assert terrain.shape == data.shape, "{} != {}".format(terrain.shape, data.shape)

    surface_intensity = weighted_intensity(
        terrain,
        azimuth=azimuth,
        elevation=elevation,
        ambient_weight=ambient_weight,
        lamp_weight=lamp_weight,
    )

    rgba = color_data(data, cmap=cmap, vmin=vmin, vmax=vmax, norm=norm)
    return blend_function(rgba, surface_intensity)


# derive distance image from an array
def calculateDist(array):
    distance = ndimage.distance_transform_edt(array)
    distance[distance != 1] = 0
    dist = np.ma.masked_where(distance == 0, distance)
    return dist


# maximize plot window
def maxPlotWindow():
    figManager = plt.get_current_fig_manager()
    figManager.window.showMaximized()


# display image and legend
def display_image(img, title, legend="", max_plot=False):
    if max_plot == True:
        maxPlotWindow()
    im = plt.imshow(img)
    plt.suptitle(title, fontsize=24, fontweight="bold", color="black")
    if legend != "":
        values = np.unique(img.ravel())[1:]
        colors = [im.cmap(im.norm(value)) for value in values]
        # create a patch (proxy artist) for every color
        patches = [
            mpatches.Patch(
                color=colors[i], label=legend + " {l}".format(l=int(values[i]))
            )
            for i in range(len(values))
        ]
        # put those patched as legend-handles into the legend
        plt.legend(handles=patches, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.0)
    # plt.show()
    return True


# animate the images in a selected folder
def visual(background, img_path, interval, iterations):
    rgb = hill_shade(background)  # create hillshade image
    rows = rgb.shape[0]
    cols = rgb.shape[1]

    display_image(
        rgb, "Inundation simulation with rainfall intensity 5 cm/h", max_plot=False
    )
    plt.pause(3)
    layer = None
    for iter in range(1, iterations + 1):
        img_files = os.listdir(img_path)
        if iter % 2 == 0:  # simulate water level decreasing
            img_files = reversed(img_files)
            layer.remove()
        elif iter > 1 and iter % 2 == 1:
            plt.imshow(rgb)
            plt.pause(1)
        text = None
        for i, img in enumerate(img_files):
            if i > 0:
                layer.remove()
            img_file = os.path.join(img_path, img)
            with TiffFile(img_file) as tif:
                img_arr = tif.asarray()
            img_arr = np.ma.masked_where(img_arr == 0, img_arr)
            layer = plt.imshow(img_arr, alpha=0.5)
            if iter % 2 != 0:
                text = plt.text(
                    cols / 2,
                    -5,
                    "Water level increasing: Time = {} hours".format(i + 1),
                    fontsize=20,
                    ha="center",
                )
            else:
                text = plt.text(
                    cols / 2, -5, "Water level decreasing ...", fontsize=20, ha="center"
                )
            plt.pause(interval)
            text.remove()
    return True


if __name__ == "__main__":
    # in_dem = "../../data/dem_full.tif"
    # # in_dem = r"C:\temp\ppr\case\case_dem.tif"
    # images = r"C:\temp\ppr\flood"

    in_dem = arcpy.GetParameterAsText(0)
    iterations = int(arcpy.GetParameterAsText(1))
    images = arcpy.GetParameterAsText(2)

    desc = arcpy.Describe(in_dem)
    in_dem = desc.catalogPath

    with TiffFile(in_dem) as tif:
        bk_img = tif.asarray()

    interval = 0.0001
    # iterations = 3

    visual(bk_img, images, interval, iterations)
    plt.show()
    quit()
