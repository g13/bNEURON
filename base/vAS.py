from neuroAlter import proceed, prepCell
from n128 import *
import numpy as np
import sys
import matplotlib
matplotlib.use('Agg')
#matplotlib.use('Qt5Agg')
from matplotlib import pyplot

def write_one(data_filename,data,mode='ab'):
    with open(data_filename,mode) as data_file:
        data.tofile(data_file)

if __name__ == '__main__':
    pas = False 
    theme = 'act-BT'
    tstep = 1.0/10.0
    run_t = 300
    transL = 180
    #run_t = 200
    #transL = 10
    run_nt = int(round(run_t/tstep))
    tol_t = min([300,run_t])
    tol_nt = int(round(tol_t/tstep))+1
    seed = 231271 
    np.random.seed(seed)
    #locE = np.array([79, 82, 83, 108, 124, 129],dtype='int')
    locE = np.array([60, 72, 78, 84, 90, 98],dtype='int')
    locI = np.array([14, 28, 30],dtype='int')
    #locE = np.array([72,79],dtype='int')
    #locI = np.array([14,40],dtype='int')
    #locE = np.array([79, 82, 83, 98, 120, 124],dtype='int')
    #locI = np.array([14, 28, 40],dtype='int')
    g0 = 32.0*1e-3
    #gE = np.array([0.6, 0.6, 0.2, 0.6, 0.15, 0.6]) * g0
    #gI = -g0*np.array([6.0, 10.0, 8.0])
    gE = np.array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0]) * g0
    gI = -g0*np.array([10.0, 10.0, 10.0])
    posE = np.array([0.3,0.3,0.9,0.6,0.4,0.2])
    posI = np.array([0.7,0.2,0.5])

    if pas:
        vNS = np.zeros((2,2))
        vAS = vNS
        norm_vAS = vAS
        i = 2
        write_one('../library/'+theme+'/'+theme+'-vNSrange.bin',vNS[0,:i],'wb')
        write_one('../library/'+theme+'/'+theme+'-vNS.bin',vNS[:,:i],'wb')
        write_one('../library/'+theme+'/'+theme+'-vASrange.bin',vAS[0,:i],'wb')
        write_one('../library/'+theme+'/'+theme+'-vAS.bin',norm_vAS,'wb')
        sys.exit(0)
    posE = posE[:locE.size]
    posI = posI[:locI.size]
    gE = gE[:locE.size]
    gI = gI[:locI.size]
    pos = np.concatenate((posE, posI))
    loc = np.concatenate((locE, locI))
    gList = np.concatenate((gE, gI))
    vrest = -70
    n = loc.size
    alphaR = True 
    cell, vecStimList, synList, distance = prepCell(gList, loc, pos, n, vrest, alphaR)

    sel = range(n)

    vecTuple = [np.array([run_t+1]),np.array([run_t+1]),np.array([run_t+1]),np.array([run_t+1]),np.array([run_t+1]),np.array([run_t+1]),np.array([run_t+1]),np.array([run_t+1]),np.array([run_t+1])]
    #vecTuple = (np.array([0,330,run_t+1]),np.array([run_t+1]),np.array([run_t+1]),np.array([run_t+1]),np.array([run_t+1]),np.array([run_t+1]),np.array([run_t+1]),np.array([run_t+1]),np.array([0,660,run_t+1]))
    #RList = np.array([[0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0]],dtype='double')
    RList0 = np.zeros((9,2))
    dendVclamp0 = np.array([1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000])
    t0 = 0.0
    tref = 13
    vThres = -65.0
    #v_pre = -70
    v_pre = vThres
    vBack = -61.0
    printR = True
    getDendV = False 
    spike = True 
    nonSpike = True
    oneGo = True 
    clampDend = False
    monitorDend = False 
    trans0 = 11
    ntransL = int(round(transL/tstep))
    printSlice = True 
    #vASt = 100 
    #nvAS = 3
    vASt = 150 
    nvAS = 5
    passive = True
    if passive:
        ASdendEv = -50
        vASrange = np.linspace(-60,-54,nvAS)
    else:
        ASdendEv = -50
        vASrange = np.linspace(-65,-55,nvAS)
    
    #vNSt = 100
    #nvNS = 3
    vNSt = 200
    nvNS = 5
    if passive:
        vNSrange = np.linspace(-65,-60,nvNS)
        NSdendEv = 5
    else:
        vNSrange = np.linspace(-66,-58,nvNS)
        NSdendEv = 5
    vBuffer = -1;

    if nonSpike:
        print vNSrange
        i = 0
        i2 = 0
        dendVclamp = dendVclamp0.copy()
        pyplot.figure('superShape', figsize=(8,4))
        vNSnt = int(round(vNSt/tstep))+1
        vNS = np.empty((vNSnt,nvNS))
        dvNS = np.empty((vNSnt,nvNS))
        RList = RList0.copy();
        trans = trans0
        superShape = np.empty((nvNS,run_nt+1))
        superSize = np.empty(nvNS,dtype='int')
        for v_ in vNSrange:
            v_pre = v_
            trans = trans - 1
            ntrans = round(trans/tstep)
            dendVclamp[:locE.size] = NSdendEv + v_pre
            dendVclamp[6:9] = v_
            fign = 'vNS0'+str(v_)
            v, fired, tsp, ntrans = proceed(cell, v_pre, synList, RList, vecStimList, vecTuple, n, trans, trans+run_t, vBack, tref, vThres, oneGo, t0, tstep, loc, pos, dendVclamp, alphaR, getDendV, False, pas, printSlice,fign)
            if not fired:
                v = v[ntrans:]
                b = 0
                for iv in np.arange(v.size):
                    if v[iv] <= v_pre+vBuffer:
                        b = iv
                        break
                assert(b>0)
                vNS[:,i] = v[b:b+vNSnt]
                superSize[i] = b
                superShape[i,:superSize[i]] = v[:b]
                print 'size of superShape', b
                pyplot.figure('superShape', figsize=(8,4))
                pyplot.plot(superShape[i,:superSize[i]],':')
                fign = 'vNS2'+str(v_)
                v0, fired0, _, _ = proceed(cell, v_pre, synList, RList, vecStimList, vecTuple, n, transL, transL+run_t, vBack, tref, vThres, oneGo, t0, tstep, loc, pos, dendVclamp0, alphaR, getDendV, False, pas, printSlice, fign)
                if not fired0:
                    v0 = v0[ntransL:]
                    b0 = 0
                    for iv in np.arange(v0.size):
                        if v0[iv] <= v_pre+vBuffer:
                            b0 = iv
                            break
                    assert(b0>0)
                    dvNS[:,i2] = vNS[:,i] - v0[b0:b0+vNSnt]
                    i2 = i2 + 1
                else:
                    print "leaky fire, why?"
                i = i + 1
        pyplot.savefig('superShape.png',format='png',bbox_inches='tight',dpi=900)
        pyplot.close('slice')
        write_one('../library/'+theme+'/'+theme+'-vNSrange.bin',vNS[0,:i],'wb')
        print vNS[0,:i]
        write_one('../library/'+theme+'/'+theme+'-vNS.bin',vNS[:,:i],'wb')
        pyplot.figure('dvNS', figsize=(8,4))
        pyplot.plot(dvNS[:,:i2])
        pyplot.savefig('dvNS.png',format='png',bbox_inches='tight',dpi=900)
        pyplot.figure('vNS', figsize=(8,4))
        pyplot.plot(vNS[:,:i])
        pyplot.savefig('vNS.png',format='png',bbox_inches='tight',dpi=900)
         

    if spike:
        print vASrange
        i = 0
        dendVclamp = dendVclamp0.copy()
        dendVclamp[:locE.size] = ASdendEv
        vASnt = int(round(vASt/tstep))+1
        vAS = np.empty((vASnt,nvAS))
        RList = RList0.copy();
        trans = trans0
        spikeShape = np.empty((nvAS,run_nt+1))
        spikeSize = np.empty(nvAS,dtype='int')
        for v_ in vASrange:
            v_pre = v_
            trans = trans - 1
            ntrans = round(trans/tstep)
            dendVclamp[6:9] = v_
            fign = 'vAS0'+str(v_)
            v, fired, tsp, ntrans = proceed(cell, v_pre, synList, RList, vecStimList, vecTuple, n, trans, trans+run_t, vBack, tref, vThres, oneGo, t0, tstep, loc, pos, dendVclamp, alphaR, getDendV, False, pas, printSlice,fign)
            if fired:
                v = v - cell.Vrest
                v = v[ntrans:]
                b = v.argmin()
                vAS[:,i] = v[b:b+vASnt]
                spikeSize[i] = b
                spikeShape[i,:spikeSize[i]] = v[:b]
                print 'size of spikeShape', b
                pyplot.figure('spikeShape', figsize=(8,4))
                pyplot.plot(spikeShape[i,:spikeSize[i]],':')
                i = i + 1
        pyplot.savefig('spikeShape.png',format='png',bbox_inches='tight',dpi=900)
        pyplot.close('slice')
        norm_vAS = vAS[:,:i]/np.abs(vAS[:,:i]).max(0).reshape(1,i)

        pyplot.figure('vAS',figsize=(8,4))
        t = np.arange(vASnt)*tstep
        t0 = np.arange(run_nt+1)*tstep
        pyplot.plot(t,vAS[:,:i],':')
        pyplot.xlim(0,vASt)
        pyplot.ylim(-7,0)
        pyplot.savefig('vAS.png',format='png',bbox_inches='tight',dpi=900)
        pyplot.close()
        pyplot.figure('norm_vAS',figsize=(8,4))
        pyplot.plot(t,norm_vAS,':')
        pyplot.xlim(0,vASt)
        pyplot.savefig('norm_vAS.png',format='png',bbox_inches='tight',dpi=900)
        vAS = vAS + cell.Vrest
        write_one('../library/'+theme+'/'+theme+'-vASrange.bin',vAS[0,:i],'wb')
        print vAS[0,:i]
        write_one('../library/'+theme+'/'+theme+'-vAS.bin',norm_vAS,'wb')
