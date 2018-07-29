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
from matplotlib import pyplot, ticker

def run(cell, v0, vThres, cpi):
    f = open('pylog','a')
    if cpi:
        cell.cp_init()
    else:
        cell.init(v0)
    print '    cell initiated'
    print >>f, '    cell initiated'
    tsps = np.empty(int((h.tstop)/h.dt/2))
    vold = v0
    nc = 0
    tsp = 0
    firing = 0 
    h.t = 0
    while (int(round(h.t/h.dt)) < int(round(h.tstop/h.dt))):
        h.fadvance()
        if cell.soma(0.5).v > vThres+15 and cell.soma(0.5).v < vold and not firing:
            tsp = h.t
            tsps[nc] = tsp
            firing = 1
            nc = nc + 1
            print >>f, '    fired at', tsp, ', ', nc, 'spike(s) in total'
            print '    fired at', tsp, ', ', nc, 'spike(s) in total'
        vold = cell.soma(0.5).v
        if (cell.soma(0.5).v <= vThres+7.5):
            firing = 0

    print  'stopping with v', cell.soma(0.5).v
    print >>f, 'stopping with v', cell.soma(0.5).v
    f.close()
    return nc
#
def bproceed(cell, v0, synList, gList, vecStimList, spikeTrain, n, sel, tend, tstep, name, pos, loc, alphaR = True, cpi = False):
    f = open('pylog','a')
    h.tstop = tend
    h.dt = tstep
    print "tend = ", tend, ', v0 = ', v0
    print >> f, "tend = ", tend, ', v0 = ', v0

    print 'gList:'
    print gList
    v = h.Vector()
    v.record(cell.soma(0.5)._ref_v)
    t = h.Vector()
    t.record(h._ref_t)
    dendv = []
    for i in xrange(n):
        dendv.append(h.Vector())
        dendv[i].record(cell.dend[loc[i]](pos[i])._ref_v)
    for i in xrange(len(sel)):
        if spikeTrain[i].size > 0:
            print ' assigning spike events to', sel[i], spikeTrain[i]
            print >>f, ' assigning spike events to', sel[i], spikeTrain[i]
            vecStimList[sel[i]].play(h.Vector(spikeTrain[i]))
        else:
            vecStimList[sel[i]].play(h.Vector([]))
        vecStimList[sel[i]].dt = h.dt
        if alphaR:
            synList[sel[i]].f = abs(gList[sel[i]])
            synList[sel[i]].deltat = h.dt
        else:
            synList[sel[i]].gmax = abs(gList[sel[i]])
    for i in xrange(n):
        if i not in sel:
            vecStimList[i].play(h.Vector([]))
            vecStimList[i].dt = h.dt
            if alphaR:
                synList[i].deltat = h.dt
                synList[i].f = 0
            else:
                synList[i].gmax = 0

    fired = run(cell, v0, vThres, cpi)

    v1 = v.as_numpy()
    dendv1 = np.empty((run_nt,n))
    for i in xrange(n):
        dvtmp = dendv[i].as_numpy()
        dendv1[:,i] = dvtmp
    
    if name:
        t1 = t.as_numpy()
        fig = pyplot.figure(name,figsize=(8,4))
        ax = fig.add_subplot(1,1,1)
        ax.plot(t1,v1)
        ax.plot(t1,dendv1)
        pyplot.savefig(name+'.png',format='png',bbox_inches='tight',dpi=900)
        pyplot.close()
    f.close()
    return v1.copy(), fired, dendv1

def getLeaky(vi,vRange,cell,gList,loc,pos,sL,vSL,n,tstep,run_t,alphaR,sender):
    spikeTrain = [np.array([])]
    sel = []
    v0 = vRange[vi]
    leakyV, fired, leakyDendv = bproceed(cell, v0, sL, gList, vSL, spikeTrain, n, sel, run_t, tstep, '', pos, loc, alphaR)
    if fired:
        print "spontaneously fired"
        assert(not fired)
    sender.send(np.array([vi, leakyV, leakyDendv]))

def getSinglets(i,v0,cell,gList,loc,pos,sL,vSL,n,tstep,run_t,alphaR,include0,leakyV0,leakyDendV0,sender):
    sel = [i]
    spikeTrain = [np.array([0])]

    V, fired, dendV = bproceed(cell, v0, synList, gList, vSL, spikeTrain, n, sel, run_t, tstep, '', pos, loc, alphaR)
    if fired:
        print "single input leads to fire"
        assert(not fired)
    V = V - v0
    dendV = dendV - v0
    tMax = np.abs(V).argmax()
    returnArray = [i,V,dendV,tMax]
    if include0:
        cpi = True
        V0, fired, dendV0 = bproceed(cell, v0, synList, gList, vSL, spikeTrain, n, sel, run_t, tstep, fign, pos, loc, alphaR, cpi)
        if fired:
            print "single input leads to fire"
            assert(not fired)
        V0 = V0 - leakyV0
        dendV0 = dendV0 - leakyDendV0
        tMax0 = np.abs(V0).argmax()
        returnArray.append(V0)
        returnArray.append(dendV0)
        returnArray.append(tMax0)
    sender.send(returnArray)

def getSingletsV(vi,vRange,cell,gList,loc,pos,sL,vSL,n,tstep,run_t,run_nt,alphaR,include0,leakyV0,leakDendV0,sender):
    cell.set_volt(v0)
    V = np.empty((run_nt,n))
    dendV = np.empty((run_nt,n,n))
    tMax = np.empty(n)
    if include0:
        V0 = np.empty((run_nt,n))
        dendV0 = np.empty((run_nt,n,n))
        tMax0 = np.empty(n)

    jobs = []
    receivers = []
    for i in xrange(n):
        receiver, sender = mp.Pipe(False)
        job = mp.Process(target=getSinglets, args = (i,v0,cell,gList,loc,pos,sL,vSL,n,tstep,run_t,alphaR,include0,leakyV0,leakyDendV0,sender))
         jobs.append(job)
         receivers.append(receiver)
         job.start()
    result, argi = gather_and_distribute_results(receivers, jobs, getSinglets.__name__, n, sortingInd = 0):
    for i in xrange(n):
        iSorted = argi[i]
        V[:,i] = result[i][1]
        dendV[:,:,i] = result[i][2]
        tMax[i] = result[i][3]
        if include0:
            V0[:,i] = result[i][4]
            dendV0[:,:,i] = result[i][5]
            tMax0[i] = result[i][6]
    returnArray = [vi,V,dendV,tMax]
    if include0:
        returnArray.append(V0)
        returnArray.append(dendV0)
        returnArray.append(tMax0)
    sender.send(returnArray)

def gather_and_distribute_results(receivers, jobs, targetName, n, sortingInd = 0):
    print "gather " + targetName + " results"
    # block job until data recevied
    result = np.array([receiver.recv() for receiver in receivers])
    # wait for all jobs to finish
    for job in jobs:
        job.join()
    # sort received result indices
    ind = [result[i][sortingInd] for i in xrange(n)]
    argi = np.argsort(ind)
    return result, argi

if __name__ == '__main__':
    compare = False
    include0 = True
    g0 = 32.0*5e-4
    fmt = 'png'
    loadSingle = True 
    tstep = 1.0/10.0
    run_t = 200.0
    run_nt = int(round(run_t/tstep))+1

    locE = np.array([60, 72],dtype='int')
    locI = np.array([14, 28],dtype='int')
    #locE = np.array([60, 72, 78, 84, 90, 98],dtype='int')
    #locI = np.array([14, 28, 30],dtype='int')
    gE = np.array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0]) * g0
    gI = -g0*np.array([10.0, 10.0, 10.0])

    #locE = np.array([79, 82, 83, 98, 120, 124],dtype='int')
    #locI = np.array([14, 28, 40],dtype='int')
    #gE = np.array([0.6, 0.6, 0.2, 0.6, 0.15, 0.6]) * g0
    #gI = -g0*np.array([6.0, 10.0, 8.0])

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
    vThres = -60.0
    n = loc.size
    alphaR = True 
    cell, vSL, synList, _ = prepCell(gList, loc, pos, n, vrest, alphaR)
    
    #vRange = np.array([-74,-70,-66,-62],dtype='double')
    vRange = np.array([-74,-70,-66],dtype='double')
    dt = 10.0
    theme = 'test'
    datafn = 'leakySingle-' + theme + '.npy'
    # linear and leaky
    plotLeakySingle = False 
    t = np.arange(run_nt)*tstep
    cpi = True
    nv = vRange.size
    if not loadSingle:
        V = np.empty((run_nt,n,nv))
        dendV = np.empty((run_nt,n,n,nv))
        tMax = np.empty((n,nv))
        if include0:
            V0 = np.empty((run_nt,n,nv))
            dendV0 = np.empty((run_nt,n,n,nv))
            tMax0 = np.empty((n,nv))
            leakyV0 = np.empty((run_nt,nv))
            leakyDendV0 = np.empty((run_nt,n,nv))
            print "get Leaky"
            jobs = []
            receivers = []
            for vi in xrange(nv):
                receiver, sender = mp.Pipe(False)
                job = mp.Process(target=getLeaky, args = (vi,vRange,cell,gList,loc,pos,sL,vSL,n,v0,tstep,run_t,alphaR,sender))
                jobs.append(job)
                receivers.append(receiver)
                job.start()
            result, argvi = gather_and_distribute_results(receivers, jobs, getLeaky.__name__, nv)
            for vi in xrange(n):
                viSorted = argvi[vi]
                leaky0[:,vi] = result[viSorted][1]
                leakyDend0[:,vi] = result[viSorted][2]

        print "get singlets"
        jobs = []
        receivers = []
        for vi in xrange(nv):
            receiver, sender = mp.Pipe(False)
            job = mp.Process(target=getSingletsV, args = (vi,vRange,cell,gList,loc,pos,sL,vSL,n,tstep,run_t,run_nt,alphaR,include0,leakyV0,leakDendV0,sender))
            jobs.append(job)
            receivers.append(receiver)
            job.start()
            result, argvi = gather_and_distribute_results(receivers, jobs, getSingletsV.__name__, nv)
        for vi in xrange(nv):
            viSorted = argvi[vi]
            V = [:,:,vi] = result[viSorted][1]
            dendV[:,:,:,vi] = result[viSorted][2]
            tMax[:,vi] = result[viSorted][3]
            if include0:
                V0 = [:,:,vi] = result[viSorted][4]
                dendV0[:,:,:,vi] = result[viSorted][5]
                tMax0[:,vi] = result[viSorted][6]

        if include0:
            with open(datafn, 'w') as leakySingleData:
                np.savez(leakySingleData, V=V, dendV=dendV,tMax=tMax, V0=V0, dendV0=dendV0,tMax0=tMax0, leakyV0=leakyV0, leakyDendV0=leakyDendV0)
        else:
            with open(datafn, 'w') as leakySingleData:
                np.savez(leakySingleData, V=V, dendV=dendV,tMax=tMax)

    else:
        lsData = np.load(datafn)
        V = lsData['V']
        dendV = lsData['dendV']
        tMax = lsData['tMax']
        print V.shape
        print dendV.shape
        print tMax.shape
        if include0:
            V0 = lsData['V0']
            dendV0 = lsData['dendV0']
            tMax0 = lsData['tMax0']
            leakyV0 = lsData['leakyV0']
            leakyDendV0 = lsData['leakyDendV0']
            print V0.shape
            print dendV0.shape
            print tMax0.shape
            print leakyV0.shape
            print leakyDendV0.shape
        #leakyDend = lsData['leakyDend']
        #leakyDend0 = lsData['leakyDend0']
        #sDendv = lsData['sDendv']
        #sDendv0 = lsData['sDendv0']
    if plotLeakySingle:
        fign = 'leaky-single'
        print "plot " + fign
        prop_cycle = pyplot.rcParams['axes.prop_cycle']
        colors = prop_cycle.by_key()['color']
        fig = pyplot.figure(fign, figsize=(8,4))
        ax1 = fig.add_subplot(2,2,1)
        ax1.plot(t,leaky0)
        ax1.set_title('leaky')

        ax3 = fig.add_subplot(2,2,3)
        for vi in xrange(nv):
            ax3.set_prop_cycle('color',colors)
            ax3.plot(t,sv[:,:,vi],lw=(vi+1)*0.5)
        ax3.set_title('sv_from_rest')

        ax4 = fig.add_subplot(2,2,4)
        for vi in xrange(nv):
            ax4.set_prop_cycle('color',colors)
            ax4.plot(t,sv0[:,:,vi],lw=(vi+1)*0.5)
        ax4.set_title('sv')

        ax2 = fig.add_subplot(2,2,2)
        for vi in xrange(nv):
            ax2.set_prop_cycle('color',colors)
            ax2.plot(t,sv[:,:,vi]-sv0[:,:,vi],lw=(vi+1)*0.5)
        ax2.set_title('dv')
        pyplot.savefig(fign+'.'+fmt,format=fmt,bbox_inches='tight',dpi=900)
        print "saved "+fign +'.'+fmt
