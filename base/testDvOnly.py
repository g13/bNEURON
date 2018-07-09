from neuron import h
h.load_file('stdlib.hoc')
h.load_file('stdrun.hoc')
from neuroAlter import prepCell
import multiprocessing as mp
import numpy as np
from n128 import *
import matplotlib
matplotlib.use('Agg')
#matplotlib.use('Qt5Agg')
from matplotlib import pyplot#, ticker
from sv import bproceed, gather_and_distribute_results

def test(vi,vRange,cell,vThres,gList,loc,pos,sL,vSL,n,tstep,run_t,alphaR,normVrest, spikeTrain, sel, sender):
    v0 = vRange[vi]
    #fign = 'singlet'+str(v0)
    fign = ''
    untilTol = True
    if normVrest:
        cell.set_volt(v0)
    V, _, _ = bproceed(cell, v0, vThres, sL, gList, vSL, spikeTrain, n, sel, run_t, tstep, fign, pos, loc, alphaR, normVrest, untilTol)
    if fired:
        print 'fired'
    print 'test', vi,
    sender.send(np.array([vi, V[-1]-V[0]],dtype='object'))
    print 'sent'

def plotTest(vRange,cell,vSL,sL,pos,loc,gList,vThres,n,tstep,run_t,alphaR,normVrest,spikeTrain,sel,fign,fmt,datafn):
    nv = vRange.size
    dV = np.empty((nv))
    jobs = []
    receivers = []
    for vi in xrange(nv):
        (receiver, sender) = mp.Pipe(False)
        job = mp.Process(target=test, args = (vi,vRange,cell,vThres,gList,loc,pos,sL,vSL,n,tstep,run_t,alphaR,normVrest,spikeTrain,sel,sender))
        jobs.append(job)
        receivers.append(receiver)
        job.start()
    gather_and_distribute_results(receivers, jobs, test.__name__, nv, 0, V, dendV)
    fig = pyplot.figure(fign,figsize=(8,4))
    ax = fig.add_subplot(1,1,1)
    ax.plot(vRange,dV,'*',ms=2.0)
    pyplot.savefig(fign+'.'+fmt,format=fmt,bbox_inches='tight',dpi=600)
    with open(datafn+'.npz','w') as datafile:
        np.savez(datafile, dV=dV)

if __name__ == '__main__':
    fmt = 'png'
    datafn0 = 'dVonly-nonNormVrest'
    g0 = 32.0*2e-4
    maxE = 5
    tstep = 1.0/10.0
    run_t = 500.0
    
    #locE = np.array([],dtype='int')
    #locI = np.array([14],dtype='int')
    locE = np.array([60, 72, 78, 84, 90, 98],dtype='int')
    locI = np.array([14, 28, 30],dtype='int')
    #gE = np.array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0]) * g0
    #gI = -g0*np.array([10.0, 10.0, 10.0])
    
    #locE = np.array([79, 82, 83, 98, 120, 124],dtype='int')
    #locI = np.array([14, 28, 40],dtype='int')
    gE = np.array([0.6, 0.6, 0.2, 0.6, 0.15, 0.6]) * g0
    gI = -g0*np.array([6.0, 10.0, 8.0])
    
    posE = np.array([0.3,0.3,0.9,0.6,0.4,0.2])
    posI = np.array([0.7,0.2,0.5])
    gE = gE[:locE.size]
    gI = gI[:locI.size]
    posE = posE[:locE.size]
    posI = posI[:locI.size]
    pos = np.concatenate((posE, posI))
    loc = np.concatenate((locE, locI))
    gList = np.concatenate((gE, gI))
    vrest = -70.0
    vThres = -54.0
    alphaR = True 
    #vRange = np.array([-62],dtype='double')
    #vRange = np.array([-63,-62,-61,-60,-59,-58],dtype='double')
    #vRange = np.array([-62.0,-61.5,-61.0,-60.5,-60.0],dtype='double')
    vRange = np.arange(-61.40,-60.95,0.005)
    normVrest = False
    run_nt = int(round(run_t/tstep))+1
    nv = vRange.size
    n = loc.size
    cell, vSL, sL, _ = prepCell(gList, loc, pos, n, vrest, alphaR)

    selected = np.arange(n)
    for i in selected:
        sel = [i]
        spikeTrain = [np.array([0])]
        datafn = datafn0 + str(i)
        plotTest(vRange,cell,vSL,sL,pos,loc,gList,vThres,n,tstep,run_t,alphaR,normVrest,spikeTrain,sel,datafn,fmt,datafn)
