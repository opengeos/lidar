# Installation

**lidar** supports a variety of platforms, including Microsoft Windows,
macOS, and Linux operating systems. Note that you will need to have
**Python 3.x** (&lt; 3.9) installed. Python 2.x is not supported.
**lidar** is available on both [PyPI](https://pypi.python.org/pypi/lidar) and [conda-forge](https://anaconda.org/conda-forge/lidar).
lidar has a [GDAL](https://gdal.org/) dependency, which can be challenging to install using pip on Windows.
Therefore, it is highly recommended to install lidar from the conda-forge channel.
If you encounter any errors, please check the [Dependencies](#dependencies) section below.

## Install from PyPI

To install **lidar** from PyPI, run this command in your terminal:

```console
pip install lidar
```

## Install from conda-forage

If you have [Anaconda](https://www.anaconda.com/distribution/#download-section) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html)
installed on your computer, you can create a fresh conda environment to install lidar:

```console
conda create -n geo python=3.11
conda activate geo
conda install -c conda-forge mamba
mamba install -c conda-forge lidar
```

## Upgrade lidar

If you have installed lidar before and want to upgrade to the latest version, you can run the following command in your terminal:

```console
pip install -U lidar
```

If you use conda, you can update lidar to the latest version by running the following command in your terminal:

```console
mamba update lidar -c conda-forge
```

To install the development version from GitHub directly using Git, run the following code:

```console
pip install git+https://github.com/opengeos/lidar
```

## Dependencies

lidar's Python dependencies are listed in its [requirements.txt](https://github.com/opengeos/lidar/blob/master/requirements.txt) file. In
addition, lidar has a C library dependency: GDAL &gt;=1.11.2. How to
install GDAL in different operating systems will be explained below.
More information about GDAL can be found [here](https://trac.osgeo.org/gdal/wiki/DownloadingGdalBinaries).

### Linux

#### Debian-based Linux

The following commands can be used to install GDAL for Debian-based
Linux distributions (e.g., Ubuntu, Linux Mint).

```console
sudo add-apt-repository ppa:ubuntugis/ppa
sudo apt-get update
sudo apt-get install gdal-bin libgdal-dev
```

If you encounter any compiling errors, try the following commands.

```console
sudo apt-get install --reinstall build-essential
sudo apt-get install python3-dev
pip install wheel
```

#### Pacman-based Linux

The following commands can be used to install GDAL for Pacman-based
Linux distributions (e.g., Arch Linux, Manjaro). You might need to use
**sudo** if you encounter permission errors.

```console
sudo pacman -S yaourt --noconfirm
yaourt -S gdal --noconfirm
yaourt -S python-gdal --noconfirm
```

### macOS

For a Homebrew based Python environment, do the following.

```console
brew update
brew install gdal
```

Alternatively, you can install GDAL binaries from [kyngchaos](http://www.kyngchaos.com/software/frameworks#gdal_complete). You will
then need to add the installed location
`/Library/Frameworks/GDAL.framework/Programs` to your system path.

### Windows

The instruction below assumes that you have installed [Anaconda](https://www.anaconda.com/download). Open
**Anaconda Prompt** and enter the following commands to create a conda
environment and install required packages

```console
conda create -n geo python=3.11
conda activate geo
conda install -c conda-forge mamba
mamba install -c conda-forge lidar
```

When installing the **lidar** package, if you encounter an error
saying `Microsoft Visual C++ 14.0 is required`, please follow the steps
below to fix the error and reinstall **lidar**. More information can
be found at this link [Fix Python 3 on Windows error - Microsoft Visual C++ 14.0 is required](https://www.scivision.co/python-windows-visual-c++-14-required/).

-   Download [Microsoft Build Tools for Visual Studio 2017](https://visualstudio.microsoft.com/thank-you-downloading-visual-studio/?sku=BuildTools&rel=15)
-   Double click to install the downloaded installer - **Microsoft Build Tools for Visual Studio 2017**.
-   Open **Microsoft Build Tools for Visual Studio 2017**
-   Select **Workloads --> Visual C++ build tools** and click the install button
