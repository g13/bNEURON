import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot#, ticker
import numpy as np

def collect(dV, vRange, ax):
    dv = data['dV']
    vrange = data['vRange']
    ax.plot(vrange,dv,lw=.5)
    dV = np.append(dV,dv)
    vRange = np.append(vRange,vrange)
    return vRange, dV

fmt = 'png'
directory = './resting_potential_shift/'
fn = 'dV-nonNormVrest'
sel = '0'
vRange = np.empty([],dtype='double')
dV = np.empty([],dtype='double')
fign = 'resting_potential_shift'
fig = pyplot.figure(fign,figsize=(8,4))
ax = fig.add_subplot(1,1,1)

i = 3
data = np.load(directory+fn+sel+'-'+str(i)+'.npy')
vRange, dV = collect(dV, vRange, ax)

i = 1
data = np.load(directory+fn+sel+'-'+str(i)+'.npy')
vRange, dV = collect(dV, vRange, ax)

i = 0
data = np.load(directory+fn+sel+'-'+str(i)+'.npy')
vRange, dV = collect(dV, vRange, ax)

i = 2
data = np.load(directory+fn+sel+'-'+str(i)+'.npy')
vRange, dV = collect(dV, vRange, ax)
print vRange

pyplot.savefig(fign+'.'+fmt,format=fmt,bbox_inches='tight',dpi=600)
