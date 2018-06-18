#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  8 17:46:32 2018

@author: hasan
"""

import matplotlib.pyplot as plt

mlups = [42, 120, 364, 650]
x_indices = ['GeForce 610M','GTX 950M', 'Tesla 2090', 'GTX 980Ti']

fig, ax = plt.subplots()
ax.bar(x_indices, mlups)
#plt.xticks(['GTX 950M', 'Tesla 2090', 'GTX 980Ti'], (1,2,3) )
ax.set_ylim(0,700)
ax.set_title('Performance', fontsize=24)
ax.set_ylabel('MLUPS',fontsize=20)
ax.set_xlabel ('GPU', fontsize=20)
ax.tick_params(labelsize=20)