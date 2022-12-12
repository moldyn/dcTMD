# -*- coding: utf-8 -*-
""".. include:: ../../README.md"""
__all__ = ["work", "force", "io", "utils"]

########## 
# Crappy fix that unsures that the dcTMD module is in sys.path
import sys, os
sys.path.append(os.getcwd() + '/../')
##########


from dcTMD import work
from dcTMD import force
from dcTMD import io
from dcTMD import utils

__version__ = '0.2.0'
