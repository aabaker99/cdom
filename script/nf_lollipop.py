# library
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
 
# change color and shape and size and edges
#(markers, stemlines, baseline) = ax1.stem(values)
#plt.setp(markers, marker='D', markersize=10, markeredgecolor="orange", markeredgewidth=2)

rect = Rectangle((0,0), 5, 5)
cx = 2.5
cy = 2.5

pc = PatchCollection([rect]) # facecolor=facecolor, alpha=alpha, edgecolor=edgecolor)
ax2.set_xlim(0, 10)
ax2.set_ylim(0, 10)
ax2.add_collection(pc)
ax2.annotate('hello', (cx, cy), color='w', weight='bold', ha='center', va='center')

plt.show()
