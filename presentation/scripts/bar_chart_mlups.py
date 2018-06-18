#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  8 17:46:32 2018

@author: hasan
"""

import numpy as np
import matplotlib.pyplot as plt

mlups = [60.33, 116.502]
x_indices = ['float', 'double']

plt.ylim(0,180)
bar = plt.bar(x_indices, mlups)
plt.title('Floats vs Doubles on GTX 950m')
plt.ylabel('MLUPS')
