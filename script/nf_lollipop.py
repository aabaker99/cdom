#!usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PatchCollection
from matplotlib import gridspec
from matplotlib.patches import Rectangle

gs = gridspec.GridSpec(2, 1, height_ratios=[4,1])
ax1 = plt.subplot(gs[0])
ax2 = plt.subplot(gs[1])

# NF1
length = 2839
domains = [{
  'id': 'IPR001936',
  'start': 1187,
  'stop': 1557,
  'label': 'Ras GTPase-activating'
}, {
  'id': 'IPR001251',
  'start': 1580,
  'stop': 1738,
  'label': 'CAL-TRIO lipid binding'
}]
 
# create data
values=np.random.uniform(size=40)

# plot with no marker
# values is a list with size equal to size of gene in bp
ax1.stem(values)
ax1_xticks = ax1.get_xticks()
 
# change color and shape and size and edges
#(markers, stemlines, baseline) = ax1.stem(values)
#plt.setp(markers, marker='D', markersize=10, markeredgecolor="orange", markeredgewidth=2)

TRACK_HEIGHT = 5
rect = Rectangle((0,0), 5, TRACK_HEIGHT)
cx = 2.5
cy = 2.5

pc = PatchCollection([rect]) # facecolor=facecolor, alpha=alpha, edgecolor=edgecolor)
ax2.set_ylim(0, TRACK_HEIGHT)
ax2.add_collection(pc)

ax2.set_xticks(ax1.get_xticks())
ax2.set_xbound(ax1.get_xbound())

ax2.tick_params(axis='y', which='both', left=False, labelleft=False)
ax2.annotate('hello', (cx, cy), color='w', weight='bold', ha='center', va='center')

ax1.set_title('NF1 Lollipop Plot')
ax1.set_ylabel('Number of mutations')
ax2.set_xlabel('Position')

plt.savefig('out.png')
