#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""
from os import path as op
import io
from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

here = op.abspath(op.dirname(__file__))

# get the dependencies and installs
with io.open(op.join(here, 'requirements.txt'), encoding='utf-8') as f:
    all_reqs = f.read().split('\n')

install_requires = [x.strip() for x in all_reqs if 'git+' not in x]
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
        "Programming Language :: Python :: 2",
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
    ],
    description="Terrain and hydrological analysis based on LiDAR-derived digital elevation models (DEM)",
    entry_points={
        'console_scripts': [
            'lidar=lidar.cli:main',
        ],
    },
    install_requires=install_requires,
    dependency_links=dependency_links,
    license="MIT license",
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    keywords='lidar',
    name='lidar',
    packages=find_packages(include=['lidar']),
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/giswqs/lidar',
    version='0.1.6',
    zip_safe=False,
)
