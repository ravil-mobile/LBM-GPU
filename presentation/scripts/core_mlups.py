#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 11 15:47:42 2018

@author: hasan
"""

#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 11 02:10:15 2018

@author: hasan
"""
import matplotlib.pyplot as plt 

# GB/s
core_gtx_950 = 640
core_m_2090 = 512
core_gtx_980 = 2816
core_geforce_610m = 48

cores = [core_geforce_610m, core_gtx_950, core_m_2090, core_gtx_980]
# with floats
mlups_geforce_610m = 42
mlups_gtx_950 = 120
mlups_m_2090 = 364
mlups_gtx_980 = 650


mlups = [mlups_geforce_610m, mlups_gtx_950, mlups_m_2090, mlups_gtx_980]

text = ["Geforce 610M","GTX 950M", "Tesla 2090", "GTX 980 Ti" ]

fig, ax = plt.subplots()
ax.scatter(cores, mlups, marker="D", s=100, color=["r", "g", "b", "black"])
#ax.plot(cores, mlups, '--')

ax.set_xlim(0,3000)
ax.set_ylim(0,700)
ax.set_xlabel("Cores", fontsize=20)
ax.set_ylabel("MLUPS",fontsize=20)
ax.set_title("Core-Performance correlation", fontsize = 24)
ax.tick_params(labelsize=18)

for index, label in enumerate(text):
    # GTX 980 is outside boundaries.
    # Handle it separately
    if index == 3: 
        ax.annotate(label, (cores[index], mlups[index]),
                xytext=(cores[index] - 200 , mlups[index] - 80 ),
                size = 20)
    else:
        
        ax.annotate(label, (cores[index], mlups[index]),
                    xytext=(cores[index] + 30 , mlups[index] - 30 ),
                    size = 20)
        
fig.show()
