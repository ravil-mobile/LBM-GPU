#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 11 02:10:15 2018

@author: hasan
"""
import matplotlib.pyplot as plt 

# GB/s
freq_gtx_950 = 914
freq_m_2090 = 1850
freq_gtx_980 = 1000
freq_geforce_610m = 672

frequency = [freq_geforce_610m, freq_gtx_950, freq_m_2090, freq_gtx_980]
# with floats
mlups_geforce_610m = 42
mlups_gtx_950 = 120
mlups_m_2090 = 364
mlups_gtx_980 = 650


mlups = [mlups_geforce_610m, mlups_gtx_950, mlups_m_2090, mlups_gtx_980]

text = ["Geforce 610M","GTX 950M", "Tesla 2090", "GTX 980 Ti" ]

fig, ax = plt.subplots()
ax.scatter(frequency, mlups, marker="D", s=100, color=["r", "g", "b", "black"])
#ax.plot(frequency, mlups, '--')

ax.set_xlim(0,2200)
ax.set_ylim(0,700)
ax.set_xlabel("Frequency MHz", fontsize=15)
ax.set_ylabel("MLUPS",fontsize=15)
ax.set_title("Frequency-Performance correlation", fontsize = 24)
ax.tick_params(labelsize=18)
    
for index, label in enumerate(text):
    ax.annotate(label, (frequency[index], mlups[index]),
                xytext=(frequency[index] + 20 , mlups[index] - 30 ),
                size = 20)
    
fig.show()
