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
                tar_ref.extractall(out_dir)
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
