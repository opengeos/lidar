# -*- coding: utf-8 -*-

"""Top-level package for lidar."""

__author__ = """Qiusheng Wu"""
__email__ = 'giswqs@gmail.com'
__version__ = '0.1.6'

from .filling import ExtractSinks
from .slicing import DelineateDepressions
from .filtering import MeanFilter, MedianFilter, GaussianFilter
