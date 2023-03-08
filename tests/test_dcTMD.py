# -*- coding: utf-8 -*-
"""
Tests for the clustering module.

MIT License
Copyright (c) 2021-2022, Victor Tänzel
All rights reserved.
"""

import numpy as np
import pytest
import os
import dcTMD


VELOCITY = 0.001
RESOLUTION = 1
VERBOSE = True
TEMPERATURE = 300
INDICES = np.array([1, 3, 5])
HERE = os.path.dirname(__file__)

