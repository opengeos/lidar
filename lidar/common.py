import csv
import os
import shutil
import subprocess
import tarfile
import urllib.request
import zipfile


def in_colab_shell():
    """Tests if the code is being executed within Google Colab."""
    try:
        import google.colab  # pylint: disable=unused-variable

        return True
    except ImportError:
        return False


def is_drive_mounted():
    """Checks whether Google Drive is mounted in Google Colab.

    Returns:
        bool: Returns True if Google Drive is mounted, False otherwise.
    """
    drive_path = "/content/drive/My Drive"
    if os.path.exists(drive_path):
        return True
    else:
        return False


def download_file(
    url=None,
    output=None,
    quiet=False,
    proxy=None,
    speed=None,
    use_cookies=True,
    verify=True,
    id=None,
    fuzzy=False,
    resume=False,
    unzip=True,
    overwrite=False,
):
    """Download a file from URL, including Google Drive shared URL.

    Args:
        url (str, optional): Google Drive URL is also supported. Defaults to None.
        output (str, optional): Output filename. Default is basename of URL.
        quiet (bool, optional): Suppress terminal output. Default is False.
        proxy (str, optional): Proxy. Defaults to None.
        speed (float, optional): Download byte size per second (e.g., 256KB/s = 256 * 1024). Defaults to None.
        use_cookies (bool, optional): Flag to use cookies. Defaults to True.
        verify (bool | str, optional): Either a bool, in which case it controls whether the server's TLS certificate is verified, or a string, in which case it must be a path to a CA bundle to use. Default is True.. Defaults to True.
        id (str, optional): Google Drive's file ID. Defaults to None.
        fuzzy (bool, optional): Fuzzy extraction of Google Drive's file Id. Defaults to False.
        resume (bool, optional): Resume the download from existing tmp file if possible. Defaults to False.
        unzip (bool, optional): Unzip the file. Defaults to True.
        overwrite (bool, optional): Overwrite the file if it already exists. Defaults to False.

    Returns:
        str: The output file path.
    """

    import gdown

    if output is None:
        if isinstance(url, str) and url.startswith("http"):
            output = os.path.basename(url)

    if isinstance(url, str):
        if os.path.exists(os.path.abspath(output)) and (not overwrite):
            print(
                f"{output} already exists. Skip downloading. Set overwrite=True to overwrite."
            )
            return os.path.abspath(output)
        else:
            url = github_raw_url(url)

    if "https://drive.google.com/file/d/" in url:
        fuzzy = True

    output = gdown.download(
        url, output, quiet, proxy, speed, use_cookies, verify, id, fuzzy, resume
    )

    if unzip and output.endswith(".zip"):

        with zipfile.ZipFile(output, "r") as zip_ref:
            if not quiet:
                print("Extracting files...")
            zip_ref.extractall(os.path.dirname(output))

    return os.path.abspath(output)


def download_folder(
    url=None,
    id=None,
    output=None,
    quiet=False,
    proxy=None,
    speed=None,
    use_cookies=True,
    remaining_ok=False,
):
    """Downloads the entire folder from URL.

    Args:
        url (str, optional): URL of the Google Drive folder. Must be of the format 'https://drive.google.com/drive/folders/{url}'. Defaults to None.
        id (str, optional): Google Drive's folder ID. Defaults to None.
        output (str, optional):  String containing the path of the output folder. Defaults to current working directory.
        quiet (bool, optional): Suppress terminal output. Defaults to False.
        proxy (str, optional): Proxy. Defaults to None.
        speed (float, optional): Download byte size per second (e.g., 256KB/s = 256 * 1024). Defaults to None.
        use_cookies (bool, optional): Flag to use cookies. Defaults to True.
        resume (bool, optional): Resume the download from existing tmp file if possible. Defaults to False.

    Returns:
        list: List of files downloaded, or None if failed.
    """
    import gdown

    files = gdown.download_folder(
        url, id, output, quiet, proxy, speed, use_cookies, remaining_ok
    )
    return files


def download_from_url(url, out_file_name=None, out_dir=".", unzip=True, verbose=True):
    """Download a file from a URL (e.g., https://github.com/giswqs/whitebox/raw/master/examples/testdata.zip)

    Args:
        url (str): The HTTP URL to download.
        out_file_name (str, optional): The output file name to use. Defaults to None.
        out_dir (str, optional): The output directory to use. Defaults to '.'.
        unzip (bool, optional): Whether to unzip the downloaded file if it is a zip file. Defaults to True.
        verbose (bool, optional): Whether to display or not the output of the function
    """
    in_file_name = os.path.basename(url)

    if out_file_name is None:
        out_file_name = in_file_name
    out_file_path = os.path.join(os.path.abspath(out_dir), out_file_name)

    if verbose:
        print("Downloading {} ...".format(url))

    try:
        urllib.request.urlretrieve(url, out_file_path)
    except Exception:
        raise Exception("The URL is invalid. Please double check the URL.")

    final_path = out_file_path

    if unzip:
        # if it is a zip file
        if ".zip" in out_file_name:
            if verbose:
                print("Unzipping {} ...".format(out_file_name))
            with zipfile.ZipFile(out_file_path, "r") as zip_ref:
                zip_ref.extractall(out_dir)
            final_path = os.path.join(
                os.path.abspath(out_dir), out_file_name.replace(".zip", "")
            )

        # if it is a tar file
        if ".tar" in out_file_name:
            if verbose:
                print("Unzipping {} ...".format(out_file_name))
            with tarfile.open(out_file_path, "r") as tar_ref:
                def is_within_directory(directory, target):
                    
                    abs_directory = os.path.abspath(directory)
                    abs_target = os.path.abspath(target)
                
                    prefix = os.path.commonprefix([abs_directory, abs_target])
                    
                    return prefix == abs_directory
                
                def safe_extract(tar, path=".", members=None, *, numeric_owner=False):
                
                    for member in tar.getmembers():
                        member_path = os.path.join(path, member.name)
                        if not is_within_directory(path, member_path):
                            raise Exception("Attempted Path Traversal in Tar File")
                
                    tar.extractall(path, members, numeric_owner=numeric_owner) 
                    
                
                safe_extract(tar_ref, out_dir)
            final_path = os.path.join(
                os.path.abspath(out_dir), out_file_name.replace(".tart", "")
            )

    if verbose:
        print("Data downloaded to: {}".format(final_path))

    return


def download_from_gdrive(gfile_url, file_name, out_dir=".", unzip=True, verbose=True):
    """Download a file shared via Google Drive
       (e.g., https://drive.google.com/file/d/18SUo_HcDGltuWYZs1s7PpOmOq_FvFn04/view?usp=sharing)

    Args:
        gfile_url (str): The Google Drive shared file URL
        file_name (str): The output file name to use.
        out_dir (str, optional): The output directory. Defaults to '.'.
        unzip (bool, optional): Whether to unzip the output file if it is a zip file. Defaults to True.
        verbose (bool, optional): Whether to display or not the output of the function
    """
    try:
        from google_drive_downloader import GoogleDriveDownloader as gdd
    except ImportError:
        print("GoogleDriveDownloader package not installed. Installing ...")
        subprocess.check_call(
            ["python", "-m", "pip", "install", "googledrivedownloader"]
        )
        from google_drive_downloader import GoogleDriveDownloader as gdd

    file_id = gfile_url.split("/")[5]
    if verbose:
        print("Google Drive file id: {}".format(file_id))

    dest_path = os.path.join(out_dir, file_name)
    gdd.download_file_from_google_drive(file_id, dest_path, True, unzip)

    return


def csv_points_to_shp(in_csv, out_shp, latitude="latitude", longitude="longitude"):
    """Converts a csv file containing points (latitude, longitude) into a shapefile.

    Args:
        in_csv (str): File path or HTTP URL to the input csv file. For example, https://raw.githubusercontent.com/giswqs/data/main/world/world_cities.csv
        out_shp (str): File path to the output shapefile.
        latitude (str, optional): Column name for the latitude column. Defaults to 'latitude'.
        longitude (str, optional): Column name for the longitude column. Defaults to 'longitude'.

    """
    import whitebox

    if in_csv.startswith("http") and in_csv.endswith(".csv"):
        out_dir = os.path.join(os.path.expanduser("~"), "Downloads")
        out_name = os.path.basename(in_csv)

        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        download_from_url(in_csv, out_dir=out_dir)
        in_csv = os.path.join(out_dir, out_name)

    wbt = whitebox.WhiteboxTools()
    in_csv = os.path.abspath(in_csv)
    out_shp = os.path.abspath(out_shp)

    if not os.path.exists(in_csv):
        raise Exception("The provided csv file does not exist.")

    with open(in_csv, encoding="utf-8") as csv_file:
        reader = csv.DictReader(csv_file)
        fields = reader.fieldnames
        xfield = fields.index(longitude)
        yfield = fields.index(latitude)

    wbt.csv_points_to_vector(in_csv, out_shp, xfield=xfield, yfield=yfield, epsg=4326)


def csv_to_shp(in_csv, out_shp, latitude="latitude", longitude="longitude"):
    """Converts a csv file with latlon info to a point shapefile.

    Args:
        in_csv (str): The input csv file containing longitude and latitude columns.
        out_shp (str): The file path to the output shapefile.
        latitude (str, optional): The column name of the latitude column. Defaults to 'latitude'.
        longitude (str, optional): The column name of the longitude column. Defaults to 'longitude'.
    """
    import csv
    import shapefile as shp

    if in_csv.startswith("http") and in_csv.endswith(".csv"):
        out_dir = os.path.join(os.path.expanduser("~"), "Downloads")
        out_name = os.path.basename(in_csv)

        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        download_from_url(in_csv, out_dir=out_dir)
        in_csv = os.path.join(out_dir, out_name)

    out_dir = os.path.dirname(out_shp)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    try:
        points = shp.Writer(out_shp, shapeType=shp.POINT)
        with open(in_csv, encoding="utf-8") as csvfile:
            csvreader = csv.DictReader(csvfile)
            header = csvreader.fieldnames
            [points.field(field) for field in header]
            for row in csvreader:
                points.point((float(row[longitude])), (float(row[latitude])))
                points.record(*tuple([row[f] for f in header]))

        out_prj = out_shp.replace(".shp", ".prj")
        with open(out_prj, "w") as f:
            prj_str = 'GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["Degree",0.0174532925199433]] '
            f.write(prj_str)

    except Exception as e:
        print(e)


def clone_repo(out_dir=".", unzip=True):
    """Clones the lidar GitHub repository.

    Args:
        out_dir (str, optional): Output folder for the repo. Defaults to '.'.
        unzip (bool, optional): Whether to unzip the repository. Defaults to True.
    """
    url = "https://github.com/giswqs/lidar/archive/master.zip"
    filename = "lidar-master.zip"
    download_from_url(url, out_file_name=filename, out_dir=out_dir, unzip=unzip)


def check_install(package):
    """Checks whether a package is installed. If not, it will install the package.

    Args:
        package (str): The name of the package to check.
    """
    import subprocess

    try:
        __import__(package)
        # print('{} is already installed.'.format(package))
    except ImportError:
        print("{} is not installed. Installing ...".format(package))
        try:
            subprocess.check_call(["python", "-m", "pip", "install", package])
        except Exception as e:
            print("Failed to install {}".format(package))
            print(e)
        print("{} has been installed successfully.".format(package))


def update_package():
    """Updates the lidar package from the lidar GitHub repository without the need to use pip or conda.
    In this way, I don't have to keep updating pypi and conda-forge with every minor update of the package.

    """
    import shutil

    try:
        download_dir = os.path.join(os.path.expanduser("~"), "Downloads")
        if not os.path.exists(download_dir):
            os.makedirs(download_dir)
        clone_repo(out_dir=download_dir)

        pkg_dir = os.path.join(download_dir, "lidar-master")
        work_dir = os.getcwd()
        os.chdir(pkg_dir)

        if shutil.which("pip") is None:
            cmd = "pip3 install ."
        else:
            cmd = "pip install ."

        os.system(cmd)
        os.chdir(work_dir)

        print(
            "\nPlease comment out 'lidar.update_package()' and restart the kernel to take effect:\nJupyter menu -> Kernel -> Restart & Clear Output"
        )

    except Exception as e:
        raise Exception(e)


def check_package(name, URL=""):

    try:
        __import__(name.lower())
    except Exception:
        raise ImportError(
            f"{name} is not installed. Please install it before proceeding. {URL}"
        )


def is_tool(name):
    """Check whether `name` is on PATH and marked as executable."""

    from shutil import which

    return which(name) is not None


def random_string(string_length=3):
    """Generates a random string of fixed length.

    Args:
        string_length (int, optional): Fixed length. Defaults to 3.

    Returns:
        str: A random string
    """
    import random
    import string

    # random.seed(1001)
    letters = string.ascii_lowercase
    return "".join(random.choice(letters) for i in range(string_length))


def github_raw_url(url):
    """Get the raw URL for a GitHub file.

    Args:
        url (str): The GitHub URL.
    Returns:
        str: The raw URL.
    """
    if isinstance(url, str) and url.startswith("https://github.com/") and "blob" in url:
        url = url.replace("github.com", "raw.githubusercontent.com").replace(
            "blob/", ""
        )
    return url


def mosaic(images, output, ext=".tif", merge_args={}, verbose=True, **kwargs):
    """Mosaics a list of images into a single image. Inspired by https://bit.ly/3A6roDK.

    Args:
        images (str | list): An input directory containing images or a list of images.
        output (str): The output image filepath.
        ext (str, optional): The image file extension. Defaults to '.tif'.
        merge_args (dict, optional): A dictionary of arguments to pass to the rasterio.merge function. Defaults to {}.
        verbose (bool, optional): Whether to print progress. Defaults to True.

    """
    from rasterio.merge import merge
    import rasterio as rio
    from pathlib import Path
    import shutil

    output = os.path.abspath(output)

    if isinstance(images, str):
        path = Path(images)
        raster_files = list(path.iterdir())
        raster_files = [f for f in raster_files if f.suffix == ext]
    elif isinstance(images, list):
        raster_files = images
    else:
        raise ValueError("images must be a list of raster files.")

    if len(raster_files) == 0:
        print("No raster files found.")
        return
    elif len(raster_files) == 1:
        shutil.copyfile(raster_files[0], output)
        return

    raster_to_mosiac = []

    if not os.path.exists(os.path.dirname(output)):
        os.makedirs(os.path.dirname(output))

    for index, p in enumerate(raster_files):
        if verbose:
            print(f"Reading {index+1}/{len(raster_files)}: {os.path.basename(p)}")
        raster = rio.open(p, **kwargs)
        raster_to_mosiac.append(raster)

    if verbose:
        print("Merging rasters...")
    arr, transform = merge(raster_to_mosiac, **merge_args)

    output_meta = raster.meta.copy()
    output_meta.update(
        {
            "driver": "GTiff",
            "height": arr.shape[1],
            "width": arr.shape[2],
            "transform": transform,
        }
    )

    with rio.open(output, "w", **output_meta) as m:
        m.write(arr)


def geometry_bounds(geometry, decimals=4):
    """Returns the bounds of a geometry.

    Args:
        geometry (dict): A GeoJSON geometry.
        decimals (int, optional): The number of decimal places to round the bounds to. Defaults to 4.

    Returns:
        list: A list of bounds in the form of [minx, miny, maxx, maxy].
    """
    if isinstance(geometry, dict):
        if "geometry" in geometry:
            coords = geometry["geometry"]["coordinates"][0]
        else:
            coords = geometry["coordinates"][0]

    else:
        raise ValueError("geometry must be a GeoJSON-like dictionary.")

    x = [p[0] for p in coords]
    y = [p[1] for p in coords]
    west = round(min(x), decimals)
    east = round(max(x), decimals)
    south = round(min(y), decimals)
    north = round(max(y), decimals)
    return [west, south, east, north]


def reproject_image(image, output, dst_crs="EPSG:4326", resampling="nearest", **kwargs):
    """Reprojects an image.

    Args:
        image (str): The input image filepath.
        output (str): The output image filepath.
        dst_crs (str, optional): The destination CRS. Defaults to "EPSG:4326".
        resampling (Resampling, optional): The resampling method. Defaults to "nearest".
        **kwargs: Additional keyword arguments to pass to rasterio.open.

    """
    import rasterio as rio
    from rasterio.warp import calculate_default_transform, reproject, Resampling

    if isinstance(resampling, str):
        resampling = getattr(Resampling, resampling)

    image = os.path.abspath(image)
    output = os.path.abspath(output)

    if not os.path.exists(os.path.dirname(output)):
        os.makedirs(os.path.dirname(output))

    with rio.open(image, **kwargs) as src:
        transform, width, height = calculate_default_transform(
            src.crs, dst_crs, src.width, src.height, *src.bounds
        )
        kwargs = src.meta.copy()
        kwargs.update(
            {
                "crs": dst_crs,
                "transform": transform,
                "width": width,
                "height": height,
            }
        )

        with rio.open(output, "w", **kwargs) as dst:
            for i in range(1, src.count + 1):
                reproject(
                    source=rio.band(src, i),
                    destination=rio.band(dst, i),
                    src_transform=src.transform,
                    src_crs=src.crs,
                    dst_transform=transform,
                    dst_crs=dst_crs,
                    resampling=resampling,
                    dst_resolution=(10, 10),
                    **kwargs,
                )


def check_file_path(file_path, make_dirs=True):

    """Gets the absolute file path.

    Args:
        file_path (str): The path to the file.
        make_dirs (bool, optional): Whether to create the directory if it does not exist. Defaults to True.

    Raises:
        FileNotFoundError: If the directory could not be found.
        TypeError: If the input directory path is not a string.

    Returns:
        str: The absolute path to the file.
    """
    if isinstance(file_path, str):
        if file_path.startswith("~"):
            file_path = os.path.expanduser(file_path)
        else:
            file_path = os.path.abspath(file_path)

        file_dir = os.path.dirname(file_path)
        if not os.path.exists(file_dir) and make_dirs:
            os.makedirs(file_dir)

        return file_path

    else:
        raise TypeError("The provided file path must be a string.")


def temp_file_path(extension):
    """Returns a temporary file path.

    Args:
        extension (str): The file extension.

    Returns:
        str: The temporary file path.
    """

    import tempfile
    import uuid

    if not extension.startswith("."):
        extension = "." + extension
    file_id = str(uuid.uuid4())
    file_path = os.path.join(tempfile.gettempdir(), f"{file_id}{extension}")

    return file_path


def clip_image(image, mask, output):
    """Clip an image by mask.

    Args:
        image (str): Path to the image file in GeoTIFF format.
        mask (str | list | dict): The mask used to extract the image. It can be a path to vector datasets (e.g., GeoJSON, Shapefile), a list of coordinates, or m.user_roi.
        output (str): Path to the output file.

    Raises:
        ImportError: If the fiona or rasterio package is not installed.
        FileNotFoundError: If the image is not found.
        ValueError: If the mask is not a valid GeoJSON or raster file.
        FileNotFoundError: If the mask file is not found.
    """
    import json

    try:
        import fiona
        import rasterio
        import rasterio.mask
    except ImportError as e:
        raise ImportError(e)

    if not os.path.exists(image):
        raise FileNotFoundError(f"{image} does not exist.")

    if not output.endswith(".tif"):
        raise ValueError("Output must be a tif file.")

    output = check_file_path(output)

    if isinstance(mask, str):
        if mask.startswith("http"):
            mask = download_file(mask, output)
        if not os.path.exists(mask):
            raise FileNotFoundError(f"{mask} does not exist.")
    elif isinstance(mask, list) or isinstance(mask, dict):

        if isinstance(mask, list):
            geojson = {
                "type": "FeatureCollection",
                "features": [
                    {
                        "type": "Feature",
                        "properties": {},
                        "geometry": {"type": "Polygon", "coordinates": [mask]},
                    }
                ],
            }
        else:
            geojson = {
                "type": "FeatureCollection",
                "features": [mask],
            }
        mask = temp_file_path(".geojson")
        with open(mask, "w") as f:
            json.dump(geojson, f)

    with fiona.open(mask, "r") as shapefile:
        shapes = [feature["geometry"] for feature in shapefile]

    with rasterio.open(image) as src:
        out_image, out_transform = rasterio.mask.mask(src, shapes, crop=True)
        out_meta = src.meta

    out_meta.update(
        {
            "driver": "GTiff",
            "height": out_image.shape[1],
            "width": out_image.shape[2],
            "transform": out_transform,
        }
    )

    with rasterio.open(output, "w", **out_meta) as dest:
        dest.write(out_image)


def join_tables(in_shp, in_csv, out_shp):
    """Joins a CSV table to a shapefile.

    Args:
        in_shp (str): Path to the input shapefile.
        in_csv (str): Path to the input CSV file.
        out_shp (str): Path to the output shapefile.
    """
    import geopandas as gpd
    import pandas as pd

    dep_df = gpd.read_file(in_shp)
    info_df = pd.read_csv(in_csv)
    if len(info_df) > 0:
        info_df.columns = [col.replace("-", "_")[:10] for col in info_df.columns]
        info_df["id"] = info_df["region_id"]
        info_df.drop("region_id", axis=1, inplace=True)
        df = pd.merge(dep_df, info_df, on="id")
        df.to_file(out_shp)
    else:
        print("No data to join")


def download_ned_by_huc(
    huc_id,
    huc_type="huc8",
    datasets=None,
    out_dir=None,
    return_url=False,
    download_args={},
    **kwargs,
):
    """Download the US National Elevation Datasets (NED) for a Hydrologic Unit region. See https://apps.nationalmap.gov/tnmaccess/#/ for more information.

    Args:
        huc_id (str): The HUC ID, for example, "01010002"
        huc_type (str, optional): The HUC type, e.g., huc2, huc4, huc8. Defaults to "huc8".
        datasets (str, optional): Comma-delimited list of valid dataset tag names. The commonly used datasets include:
            Digital Elevation Model (DEM) 1 meter
            National Elevation Dataset (NED) 1/3 arc-second Current
            National Elevation Dataset (NED) 1/9 arc-second Current
            National Elevation Dataset (NED) 1 arc-second Current
            For more information, see https://apps.nationalmap.gov/tnmaccess/#/product
            Defaults to None, which will be the (NED) 1/3 arc-second
        out_dir (str, optional): The output directory. Defaults to None, which will use the current working directory.
        return_url (bool, optional): If True, the URL will be returned instead of downloading the data. Defaults to False.
        download_args (dict, optional): The download arguments to be passed to the download_file function. Defaults to {}.

    Returns:
        list: The list of downloaded files.
    """

    import requests

    endpoint = "https://tnmaccess.nationalmap.gov/api/v1/products?"

    if datasets is None:
        datasets = "National Elevation Dataset (NED) 1/3 arc-second Current"

    if out_dir is None:
        out_dir = os.getcwd()

    kwargs["datasets"] = datasets
    kwargs["polyType"] = huc_type
    kwargs["polyCode"] = huc_id

    result = requests.get(endpoint, params=kwargs).json()
    if "errorMessage" in result:
        raise ValueError(result["errorMessage"])
    else:
        links = [x["downloadURL"] for x in result["items"]]
        for index, link in enumerate(links):
            if "historical" in link:
                link = link.replace("historical", "current")[:-13] + ".tif"
                links[index] = link

    if return_url:
        return links
    else:
        for index, link in enumerate(links):

            r = requests.head(link)
            if r.status_code == 200:
                filepath = os.path.join(out_dir, os.path.basename(link))
                print(
                    f"Downloading {index + 1} of {len(links)}: {os.path.basename(link)}"
                )
                download_file(link, filepath, **download_args)
            else:
                print(f"{link} does not exist.")


def download_ned_by_bbox(
    bbox,
    datasets=None,
    out_dir=None,
    return_url=False,
    download_args={},
    **kwargs,
):
    """Download the US National Elevation Datasets (NED) for a bounding box. See https://apps.nationalmap.gov/tnmaccess/#/ for more information.

    Args:
        bbox (list): The bounding box in the form [xmin, ymin, xmax, ymax].
        huc_type (str, optional): The HUC type, e.g., huc2, huc4, huc8. Defaults to "huc8".
        datasets (str, optional): Comma-delimited list of valid dataset tag names. The commonly used datasets include:
            Digital Elevation Model (DEM) 1 meter
            National Elevation Dataset (NED) 1/3 arc-second Current
            National Elevation Dataset (NED) 1/9 arc-second Current
            National Elevation Dataset (NED) 1 arc-second Current
            For more information, see https://apps.nationalmap.gov/tnmaccess/#/product
            Defaults to None, which will be the (NED) 1/3 arc-second
        out_dir (str, optional): The output directory. Defaults to None, which will use the current working directory.
        return_url (bool, optional): If True, the URL will be returned instead of downloading the data. Defaults to False.
        download_args (dict, optional): The download arguments to be passed to the download_file function. Defaults to {}.

    Returns:
        list: The list of downloaded files.
    """

    import requests

    endpoint = "https://tnmaccess.nationalmap.gov/api/v1/products?"

    if datasets is None:
        datasets = "National Elevation Dataset (NED) 1/3 arc-second Current"

    if out_dir is None:
        out_dir = os.getcwd()

    if isinstance(bbox, list):
        bbox = ",".join([str(x) for x in bbox])

    kwargs["datasets"] = datasets
    kwargs["bbox"] = bbox

    result = requests.get(endpoint, params=kwargs).json()
    if "errorMessage" in result:
        raise ValueError(result["errorMessage"])
    else:
        links = [x["downloadURL"] for x in result["items"]]
        for index, link in enumerate(links):
            if "historical" in link:
                link = link.replace("historical", "current")[:-13] + ".tif"
                links[index] = link

    if return_url:
        return links
    else:
        for index, link in enumerate(links):

            r = requests.head(link)
            if r.status_code == 200:
                filepath = os.path.join(out_dir, os.path.basename(link))
                print(
                    f"Downloading {index + 1} of {len(links)}: {os.path.basename(link)}"
                )
                download_file(link, filepath, **download_args)
            else:
                print(f"{link} does not exist.")


def resample(src, dst, resolution, **kwargs):
    """Resample a raster to a new resolution.

    Args:
        src (str): The source raster.
        dst (str): The destination raster.
        resolution (float): The new resolution.
    """
    from osgeo import gdal

    gdal.Warp(dst, src, xRes=resolution, yRes=resolution, **kwargs)
