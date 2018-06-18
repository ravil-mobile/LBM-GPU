#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 12 21:43:20 2018

@author: hasan
"""
import os
import numpy as np
array = 255 * np.ones ( (630, 930), dtype=int)
np.savetxt( os.path.abspath(os.path.curdir) + "/empty_mesh.txt"
           , array, fmt="%0.f")