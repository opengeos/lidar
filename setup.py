#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""
import io
import json
import os
import shutil
from os import path as op
import urllib.request
from setuptools import setup, find_packages

# Find available package versions


def pkg_versions(package_name):
    url = "https://pypi.python.org/pypi/%s/json" % (package_name,)
    text = urllib.request.urlopen(url).read()
    data = json.loads(text)
    versions = data["releases"].keys()
    return list(versions)


# Find a matching version
def find_version(version, version_list):
    match_version = None
    for v in version_list:
        if v.startswith(version):
            match_version = v
            return match_version

    return match_version


if shutil.which("gdal-config") is None:
    print("GDAL is not installed. Installing GDAL ...")
    cmd = "pip install --find-links=https://girder.github.io/large_image_wheels --no-cache GDAL"
    os.system(cmd)

# check GDAL version installed in the system
# GDAL_VERSION = os.popen("gdal-config --version").read().rstrip()
GDAL_INFO = os.popen("gdalinfo --version").read().rstrip()
GDAL_VERSION = GDAL_INFO.split(",")[0].replace("GDAL", "").lstrip()
GDAL_VERSION_NUM = str(GDAL_VERSION.replace(".", ""))
PYGDAL_VERSION = find_version(GDAL_VERSION, pkg_versions("pygdal"))

# if PYGDAL_VERSION is None:
#     print(
#         "GDAL version not found in PyPI. Please install GDAL version %s or higher."
#         % (GDAL_VERSION,)
#     )
#     exit(1)

print("GDAL version: %s" % (GDAL_VERSION,))


with open("README.md", mode="rb") as readme_file:
    readme = readme_file.read().decode("utf-8")

here = op.abspath(op.dirname(__file__))

# get the dependencies and installs
with io.open(op.join(here, "requirements.txt"), encoding="utf-8") as f:
    all_reqs = f.read().split("\n")

install_requires = [x.strip() for x in all_reqs if "git+" not in x]

# install_requires.append('pygdal==' + PYGDAL_VERSION)

dependency_links = [x.strip().replace("git+", "") for x in all_reqs if "git+" not in x]


extras_requires = {
    "all": ["geopandas", "rasterio"],
}

requirements = [
    "Click>=6.0",
]

setup_requirements = []

test_requirements = []

setup(
    author="Qiusheng Wu",
    author_email="giswqs@gmail.com",
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
    ],
    description="A Python package for delineating nested surface depressions in digital elevation data",
    entry_points={
        "console_scripts": [
            "lidar=lidar.cli:main",
        ],
    },
    install_requires=install_requires,
    extras_require=extras_requires,
    dependency_links=dependency_links,
    license="MIT license",
    long_description=readme,
    long_description_content_type="text/markdown",
    include_package_data=True,
    keywords="lidar",
    name="lidar",
    packages=find_packages(include=["lidar"]),
    setup_requires=setup_requirements,
    test_suite="tests",
    tests_require=test_requirements,
    url="https://github.com/giswqs/lidar",
    version="0.7.3",
    zip_safe=False,
)
