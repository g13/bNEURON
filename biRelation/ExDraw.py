import numpy as np
import sys, os, getopt
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot, gridspec
fmt='eps'

def drawExample(ax1,directory,typeEI):
    k = np.load(directory+'/k'+typeEI+'.npy')
    rs = np.load(directory+'/rs'+typeEI+'.npy')
    ax1.plot(k

fign = 'exampleK-RS.'+fmt
fig = pyplot.figure(fign,figsize = (15,15)

gs = gridspec.GridSpec(2, 2)
ax1 = fig.add_subplot(gs[0])

ax2 = 
drawExample(ax1,'cK-62-0-238660','EI')
'
