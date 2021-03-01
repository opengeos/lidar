#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""
import os, platform
from os import path as op
import io
from setuptools import setup, find_packages

# check GDAL version installed in the system
# GDAL_VERSION = os.popen("gdal-config --version").read().rstrip()
GDAL_INFO = os.popen("gdalinfo --version").read().rstrip()
GDAL_VERSION = GDAL_INFO.split(',')[0].replace('GDAL', '').lstrip()
GDAL_VERSION_NUM = str(GDAL_VERSION.replace(".", ""))
PYGDAL_VERSION = '2.2.2.3' # default pygdal version to install


# pygdal version to install based on the GDAL version
# GDAL version history: https://trac.osgeo.org/gdal/wiki/DownloadSource
# pygdal version history: https://pypi.org/project/pygdal/#history
if not GDAL_VERSION_NUM.isdigit():
    print("GDAL cannot be detected in your system. Please install GDAL first!")
# elif GDAL_VERSION[:3] == '2.3':
#     PYGDAL_VERSION = GDAL_VERSION + '.5'
# elif GDAL_VERSION[:3] == '2.4':
#     PYGDAL_VERSION = GDAL_VERSION + '.5'    
# elif GDAL_VERSION[:3] == '3.0':
#     PYGDAL_VERSION = GDAL_VERSION + '.5' 
else:
    PYGDAL_VERSION = GDAL_VERSION + '.5'

with open('README.md', mode = 'rb') as readme_file:
    readme = readme_file.read().decode('utf-8')

# with open('changelog.md', mode = 'rb') as history_file:
#     history = history_file.read().decode('utf-8')

here = op.abspath(op.dirname(__file__))

# get the dependencies and installs
with io.open(op.join(here, 'requirements.txt'), encoding='utf-8') as f:
    all_reqs = f.read().split('\n')

install_requires = [x.strip() for x in all_reqs if 'git+' not in x]
# if platform.system() != "Windows":
#     install_requires.append('pygdal==' + PYGDAL_VERSION)

# install pygdal for travis CI
if GDAL_VERSION == '2.2.2':
    install_requires.append('pygdal==' + PYGDAL_VERSION)

dependency_links = [x.strip().replace('git+', '') for x in all_reqs if 'git+' not in x]

requirements = ['Click>=6.0', ]

setup_requirements = [ ]

test_requirements = [ ]

setup(
    author="Qiusheng Wu",
    author_email='giswqs@gmail.com',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    description="A Python package for delineating nested surface depressions in digital elevation data",
    entry_points={
        'console_scripts': [
            'lidar=lidar.cli:main',
        ],
    },
    install_requires=install_requires,
    dependency_links=dependency_links,
    license="MIT license",
    long_description=readme,
    long_description_content_type='text/markdown',
    include_package_data=True,
    keywords='lidar',
    name='lidar',
    packages=find_packages(include=['lidar']),
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/giswqs/lidar',
    version='0.6.1',
    zip_safe=False,
)
