# Installation

## conda-forge

**lidar** supports a variety of platforms, including Microsoft Windows,
macOS, and Linux operating systems. Note that you will need to have
**Python 3.x** (&lt; 3.9) installed. Python 2.x is not supported. The
**lidar** Python package can be installed using the following command.
If you encounter any errors, please check the [Dependencies](#dependencies) section
below. The instruction below assumes that you have installed
[Anaconda](https://www.anaconda.com/download). Open **Anaconda Prompt** and enter the
following commands to create a conda environment and install required
packages.

```console
conda create -n py38 python=3.8
conda activate py38
conda install -c conda-forge mamba
mamba install -c conda-forge lidar 
```

## pip

If you have installed **lidar** before and want to upgrade to the latest
version, you can use the following command:

```console
pip install lidar -U
```

## Dependencies

lidar's Python dependencies are listed in its requirements.txt file. In
addition, lidar has a C library dependency: GDAL &gt;=1.11.2. How to
install GDAL in different operating systems will be explained below.
More informaton about GDAL can be found [here](https://trac.osgeo.org/gdal/wiki/DownloadingGdalBinaries).

It is highly recommended that you use a Python virtual environment
(e.g., conda) to test the lidar package. Please follow the [conda user
guide] to install conda if necessary. Once you have conda installed, you
can use Terminal or an Anaconda Prompt to create a Python virtual
environment. Check [managing Python environment] for more information.

Once GDAL has been installed, you can then proceed to install the
**lidar** Python package using the following command:

```console
conda create -n py38 python=3.8
conda activate py38
conda install -c conda-forge gdal 
pip install lidar
```    

### Linux

#### Debian-based Linux

The following commands can be used to install GDAL for Debian-based
Linux distributions (e.g., Ubuntu, Linux Mint).

```console
sudo add-apt-repository ppa:ubuntugis/ppa
sudo apt-get update
sudo apt-get install gdal-bin libgdal-dev
pip install lidar
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
pip install lidar
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
conda create -n py38 python=3.8
conda activate py38
conda install -c conda-forge gdal 
pip install richdem
pip install lidar
```  

When installing the **richdem** package, if you encounter an error
saying 'Microsoft Visual C++ 14.0 is required', please follow the steps
below to fix the error and reinstall **richdem**. More infomration can
be found at this link Fix [Python 3 on Windows error - Microsoft Visual C++ 14.0 is required](https://www.scivision.co/python-windows-visual-c++-14-required/).

* Download [Microsoft Build Tools for Visual Studio 2017](https://visualstudio.microsoft.com/thank-you-downloading-visual-studio/?sku=BuildTools&rel=15)
* Double click to install the downloaded installer - **Microsoft Build Tools for Visual Studio 2017**.
* Open **Microsoft Build Tools for Visual Studio 2017**
* Select **Workloads --> Visual C++ build tools** and click the install button