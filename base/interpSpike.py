from neuroAlter import proceed, vASproceed, prepCell
from n128 import *
import numpy as np
import matplotlib
matplotlib.use('Agg')
#matplotlib.use('Qt5Agg')
from matplotlib import pyplot

if __name__ == '__main__':
    theme = 'activeT'
    tstep = 1.0/10.0
    run_t = 50
    run_nt = int(round(run_t/tstep))
    tol_t = min([300,run_t])
    tol_nt = int(round(tol_t/tstep))+1
    seed = 231271 
    np.random.seed(seed)
    locE = np.array([60, 72, 78, 84, 90, 98],dtype='int')
    locI = np.array([14, 28, 30],dtype='int')
    g0 = 32.0*5e-4
    gE = np.array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0]) * g0
    gI = -g0*np.array([10.0, 10.0, 10.0])
    posE = np.array([0.3,0.3,0.9,0.6,0.4,0.2])
    posI = np.array([0.7,0.2,0.5])
    pos = np.concatenate((posE, posI))
    loc = np.concatenate((locE, locI))
    gList = np.concatenate((gE, gI))
    v0 = -70
    n = loc.size
    alphaR = True 
    cell, vecStimList, synList, distance = prepCell(gList, loc, pos, n, v0, alphaR)

    sel = range(n)

    vecTuple = [np.array([run_t+1]),np.array([run_t+1]),np.array([run_t+1]),np.array([run_t+1]),np.array([run_t+1]),np.array([run_t+1]),np.array([run_t+1]),np.array([run_t+1]),np.array([run_t+1])]
    #vecTuple = (np.array([0,330,run_t+1]),np.array([run_t+1]),np.array([run_t+1]),np.array([run_t+1]),np.array([run_t+1]),np.array([run_t+1]),np.array([run_t+1]),np.array([run_t+1]),np.array([0,660,run_t+1]))
    #RList = np.array([[0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0]],dtype='double')
    RList0 = np.zeros((9,2))
    #dendVclamp = np.array([1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000])
    dendVclamp = np.array([-50, -50, -50, -50, -50, -50, 1000, 1000, 1000])
    t0 = 0.0
    tref = 13
    vThres = -65.0
    #v_pre = -70
    v_pre = vThres
    vBack = -61.0
    printR = True
    getDendV = False 
    oneGo = True 
    clampDend = False
    monitorDend = False 
    trans0 = 11
    pas = False 
    printSlice = False
    vASt = 5
    nvAS = 6
    vASrange = np.linspace(-65,-55,nvAS)
    print vASrange
    vASnt = int(round(vASt/tstep))
    vAS = np.empty((vASnt,nvAS))
    i = 0
    dendVclamp[:locE.size] = -50
    pyplot.figure('spikeShape', figsize=(8,4))
    RList = RList0.copy();
    RList[:6,0] = 1e-3
    trans = trans0
    spike = np.empty((nvAS,run_nt+1))
    spikeSize = np.empty(nvAS,dtype='int')
    vA = np.empty((run_nt+1,nvAS))
    vSrange = np.empty(nvAS)
    for v_ in vASrange:
        v_pre = v_
        trans = trans - 1
        ntrans = round(trans/tstep)
        dendVclamp[6:9] = v_
        v, fired, tsp, ntrans = proceed(cell, v_pre, synList, RList, vecStimList, vecTuple, n, trans, trans+run_t, vBack, tref, vThres, oneGo, t0, tstep, loc, pos, dendVclamp, alphaR, getDendV, False, pas, printSlice)
        if printSlice:
            pyplot.figure('slice', figsize=(8,4))
            pyplot.savefig('vAS0'+str(v_)+'.png',format='png',bbox_inches='tight',dpi=900)
        if fired:
            v = v - cell.Vrest
            vA[:,i] = v[ntrans:]
            b = v.argmin()
            vAS[:,i] = v[b:b+vASnt]
            spikeShape = v[ntrans:b]
            print 'size of spikeShape', b-ntrans
            spikeSize[i] = b-ntrans
            spike[i,:spikeSize[i]] = spikeShape
            vSrange[i] = v_
            i = i + 1
            pyplot.figure('spikeShape', figsize=(8,4))
            pyplot.plot(spikeShape)
    pyplot.close('slice')
    spike.dump('spikes.npy')

    vi = 3
    vj = 4
    vr = 0.5
    ns = np.min([spikeSize[vi],spikeSize[vj]])
    spikeInterp = spike[vi][:ns] + vr*(spike[vj][:ns]-spike[vi][:ns])
    pyplot.figure('spikeShape', figsize=(8,4))
    pyplot.plot(spikeInterp,':')
    pyplot.savefig('spikeShape.png',format='png',bbox_inches='tight',dpi=900)
