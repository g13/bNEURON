import numpy as np
import os, time, sys
from neuroAlter import *
import multiprocessing as mp
from getPSP import write_one
from shutil import copy
from datetime import datetime
from matplotlib import pyplot
import matplotlib
matplotlib.use('Agg')

def shift(v, rdt):
    if rdt > 0:
        for i in xrange(v.size):
            if i < v.size-1:
                v[i] = v[i] + (v[i+1] - v[i]) * r
            else:
                v[i] = v[i] + (0-v[i]) * r
    return v

def shiftAdd(v, vlet, shift, i0, nt):
    asize = min(v.size,i0+nt)
    if shift == 1:
        v[i0:asize] = v[i0:asize] + vlet[1:asize-i0+1]
    else:
        for i in xrange(i0,asize):
            j = i - i0
            v[i] = v[i] + vlet[j] + (vlet[j+1] - vlet[j]) * shift

taue0 = 0.098814229249
taue1 = 8.333333333333
taui0 = 0.98814229249  
taui1 = 50.0          
def getGH(dt,f,g=0,h=0):
    if f<0:
        f = -f
        t0 = taui0
        t1 = taui1
    else:
        t0 = taue0
        t1 = taue1
    h = h + f/t0
    etr = np.exp(-dt/t0)
    etd = np.exp(-dt/t1)
    c = t0/(t1-t0)*(etd-etr)
    g = g*etd + c*h 
    h = h*etr
    return np.array([g, h])

def subOnly(gList,loc,pos,vrest,v0,presetSpike,trans,run_t,tol_t,tstep,spList,seed,name,rE,rI,linear,bilinear,bilinear0,noMoreInput_t,s,did,rdpi):
    # only available on Unix
    f = open('sOlog-'+str(did),'w')
    nE = sum([x>0 for x in gList])
    nI = sum([x<=0 for x in gList])
    n = nI + nE
    
    cell, vecStimList, synList = prepCell(gList, loc, pos, n, vrest)
    
    tol_nt = int(tol_t/tstep)
    ntrans = int(trans/tstep)
    run_nt = int(run_t/tstep)
    
    
    if presetSpike:
        spikeTrain = list([])
        for i in xrange(n):
            tmpTrain = np.empty(spList[i].size)
            k = 0
            for j in xrange(spList[i].size):
                ts = round(spList[i][j]/tstep)*tstep
                if ts <= noMoreInput_t:
                    tmpTrain[k] = ts
                    k = k+1
                else:
                    break
            spikeTrain.append(np.resize(tmpTrain,k))
    else:
        rateE = np.ones(nE)*rE/1000
        rateI = np.ones(nI)*rI/1000
        rate = np.concatenate((rateE, rateI))
        spikeTrain = list([])
        for i in xrange(n):
            tmpTrain = np.empty(int(10*rate[i]*noMoreInput_t))
            if tmpTrain.size > 0:
                tmpT = 0
                lastT = 0
                j = 0
                while tmpT < noMoreInput_t and j < tmpTrain.size:
                    dt = - np.log(np.random.random_sample(1))/rate[i]
                    tmpT = tmpT + dt
                    if dt > tstep or lastT == 0:
                        tmpTrain[j] = round(tmpT/tstep)*tstep
                        j = j + 1
                    lastT = tmpT
                if tmpTrain[j-1] > noMoreInput_t:
                    spikeTrain.append(np.resize(tmpTrain,(j-1)))
                else:
                    spikeTrain.append(np.resize(tmpTrain,(j)))
            else:
                spikeTrain.append(np.empty(0))
            print i, "th spikeTrain:", spikeTrain[i]
            print >> f, i, "th spikeTrain:", spikeTrain[i]
    
    totSpikes = sum([x.size for x in spikeTrain])
    print " total inputs = ", totSpikes
    print >> f, " total inputs = ", totSpikes
    print spikeTrain
    print >> f, spikeTrain
    
    singleSpikeTrain = np.empty(totSpikes,dtype=float)
    singleSpikeID = np.empty(totSpikes,dtype=int)
    j = 0
    for i in xrange(n):
        if spikeTrain[i].size > 0:
            singleSpikeTrain[j:j+spikeTrain[i].size] = spikeTrain[i]
            singleSpikeID[j:j+spikeTrain[i].size] = i
            j = j + spikeTrain[i].size
    args = np.argsort(singleSpikeTrain, kind='mergesort')
    print singleSpikeTrain
    print >> f, singleSpikeTrain
    print singleSpikeID
    print >> f, singleSpikeID
    print " "
    print >> f, " "
    singleSpikeTrain = singleSpikeTrain[args]
    singleSpikeID = singleSpikeID[args]
    
    print singleSpikeTrain
    print >> f, singleSpikeTrain
    print singleSpikeID
    print >> f, singleSpikeID

    assert(not singleSpikeTrain[-1] == 0)
    
    print " original sim:"
    print >> f, " original sim:"
    RList = np.zeros((n,2))
    rel = np.zeros(n,dtype=int)
    sel = range(0,n)
    simV,fired,_,_ = proceed(cell, v0, synList, RList, vecStimList, spikeTrain, n, trans, trans + run_t, 0, 0, 0, 1, 0, tstep, loc, pos)
    print simV.size
    print >> f, simV.size
    #assert(fired==0)
    
    v_pre = v0
    liV = np.zeros(max(run_nt,int((singleSpikeTrain[-1]+tol_t)/tstep))+1) + vrest
    leakyV, _= leaky(cell, v0, synList, RList, vecStimList, n, trans, trans + run_t, tstep, loc, pos)
    leakyV = leakyV - vrest
    leakyV = leakyV[ntrans:ntrans+run_nt+1]
    liV[:run_nt+1] = liV[:run_nt+1] + leakyV
    #leakyV, _, _ = proceed(cell, v0, synList, RList, vecStimList, np.empty((n,0)), n, trans, tol_t + trans, 0, 0, 0, 1, 0, tstep)
    #leakyV = leakyV - vrest
    #leakyV = leakyV[ntrans:ntrans+tol_nt+1]
    #liV[:tol_nt+1] = liV[:tol_nt+1] + leakyV
    biV = liV.copy()
    biV0 = liV.copy()
    pyplot.figure(name+'/'+name+'-'+'leakyV',figsize=(8,4))
    pyplot.plot(np.linspace(0,run_t,run_nt+1),liV[:run_nt+1])
    pyplot.savefig(name+'/'+name+'-leakyV.png',format='png',bbox_inches='tight',dpi=rdpi)
    pyplot.close()
   
    if linear:
        print " singlets for linear:"
        print >> f, " singlets for linear:"
        for i in xrange(totSpikes):
            if gList[singleSpikeID[i]] < 0.0:
                label = 'I'+str(singleSpikeID[i])+'-'
            else:
                label = 'E'+str(singleSpikeID[i])+'-'

            t0 = singleSpikeTrain[i]
            if t0 > run_t:
                break
            i0 = int(round(t0/tstep))
            i1 = i0 + 1
            v_pre = liV[i0]
            sel = singleSpikeID[i]
        
            v,_ = sproceed(cell, v_pre, synList, gList, RList, vecStimList, 0.0, n, sel, trans, trans + tol_t, tstep, '')
            leakyV, _ = leaky(cell, v_pre, synList, RList, vecStimList, n, trans, trans + tol_t, tstep, loc, pos)
            v = v[ntrans:ntrans+tol_nt+1] - leakyV[ntrans:ntrans+tol_nt+1]
        
            pyplot.figure(name+'/'+name+'-s'+label+str(i),figsize=(8,4))
            pyplot.title(label+str(i))
            pyplot.plot(np.linspace(0,tol_t,tol_nt+1),v,'g',lw = 1.5)
            pyplot.savefig(name+'/'+name+'-s'+label+str(i)+'.png',format='png',bbox_inches='tight',dpi=rdpi)
            pyplot.close()
            shiftAdd(liV,v,1,i1,tol_nt)
                
        pyplot.figure(name+'/'+name+'-liV',figsize=(8,4))
        pyplot.plot(np.linspace(0,run_t,run_nt+1),simV[ntrans:ntrans+run_nt+1], lw = 2)
        pyplot.plot(np.linspace(0,run_t,run_nt+1),liV[:run_nt+1], lw = 1)
        pyplot.savefig(name+'/'+name+'-liV.png',format='png',bbox_inches='tight',dpi=rdpi)
        pyplot.close()
    
    rel = np.zeros(n,dtype=int)
    if bilinear:
        for i in xrange(totSpikes):
            if gList[singleSpikeID[i]] < 0:
                label1 = 'I'+str(singleSpikeID[i])+'-'
            else:
                label1 = 'E'+str(singleSpikeID[i])+'-'
            print ''
            print >> f, ''
            t0 = singleSpikeTrain[i]
            if t0 > run_t:
                break
            i0 = int(round(t0/tstep))
            i1 = i0 + 1
            v_pre = biV[i0]
            sel = singleSpikeID[i]
            RList = np.zeros((n,2))
            #print 'v_pre = ', v_pre
            #print >> f, 'v_pre = ', v_pre
            v2, _ = sproceed(cell, v_pre, synList, gList, RList, vecStimList, 0.0, n, sel, trans, trans + tol_t, tstep, '')
            leakyV, _ = leaky(cell, v_pre, synList, RList, vecStimList, n, trans, trans + tol_t, tstep, loc, pos)
            v2 = v2[ntrans:ntrans+tol_nt+1] - leakyV[ntrans:ntrans+tol_nt+1]
            shiftAdd(biV,v2,1,i1,tol_nt)
            pyplot.figure(name+'/'+name+'-s'+label1+str(i),figsize=(8,4))
            pyplot.title(label1+str(i))
            pyplot.plot(np.linspace(0,tol_t,tol_nt+1),v2,'b',lw = 1.0)
            pyplot.savefig(name+'/'+name+'-s'+label1+str(i)+'.png',format='png',bbox_inches='tight',dpi=rdpi)
        
            for j in xrange(0,i):
                print i, '-', j
                print >> f, i, '-', j
                if gList[singleSpikeID[j]] < 0:
                    label2 = 'I'+str(singleSpikeID[j])+'-'
                else:
                    label2 = 'E'+str(singleSpikeID[j])+'-'
                label = label2 + label1
                t00 = singleSpikeTrain[j]
                dt = t0 - t00

                if dt < tol_t:
                    idt = int(round(dt/tstep))
                    relit = tol_nt + 1 - idt
                    RList = np.zeros((n,2))
                    if singleSpikeID[j] == singleSpikeID[i]:
                        sel = [singleSpikeID[j]]
                        spikes = np.array([[]])
                        rl0 = getGH(dt,gList[sel[0]])
                        RList[sel[0],:] = getGH(0,gList[sel[0]],rl0[0],rl0[1])
                    else:
                        sel = [singleSpikeID[j],singleSpikeID[i]]
                        spikes = np.array([np.array([]),np.array([])])
                        RList[sel[0],:] = getGH(dt,gList[sel[0]])
                        RList[sel[1],:] = getGH(0,gList[sel[1]])
                    print 'synapse:', sel
                    print >> f, 'synapse:', sel
                    v,_ = bproceed(cell, v_pre, synList, gList, RList, vecStimList, spikes, n, sel, trans, trans + tol_t - dt, 0.0, tstep, '')
                    v = v[ntrans:ntrans + relit]
                    v = v - leakyV[ntrans:ntrans+relit]
                
                    RList = np.zeros((n,2))
                    sel = [singleSpikeID[j]]
                    RList[sel[0],:] = getGH(dt,gList[sel[0]])
                    spikes = np.array([[]])
                    v1,_ = bproceed(cell, v_pre, synList, gList, RList, vecStimList, spikes, n, sel, trans, trans + tol_t - dt, 0.0, tstep, '')
                    v1 = v1[ntrans:ntrans+relit]
                    v1 = v1 - leakyV[ntrans:ntrans+relit]
        
                    addv = v2[:relit] + v1
                    bv = v - addv
                    shiftAdd(biV,bv,1,i1,relit-1)
                    pyplot.figure(name+'/'+name+'-b'+label+str(i)+'-'+str(j),figsize=(8,4))
                    pyplot.title(label+str(i)+'-'+str(j))
                    t = np.linspace(t0,t0+relit*tstep,relit)
                    pyplot.plot(t,v,'k',lw=2)
                    pyplot.plot(t,v1,'r', lw=1)
                    pyplot.plot(t,v2[:relit],'b',lw = 0.5)
                    pyplot.plot(t,addv,':g',lw = 2)
                    pyplot.savefig(name+'/'+name+'-b'+label+str(i)+'-'+str(j)+'.png',format='png',bbox_inches='tight',dpi=rdpi)
                    pyplot.close()
        
    if bilinear0:
        v_pre = vrest
        RList = np.zeros((n,2))
        for i in xrange(totSpikes):
            if gList[singleSpikeID[i]] < 0:
                label1 = 'I'+str(singleSpikeID[i])+'-'
            else:
                label1 = 'E'+str(singleSpikeID[i])+'-'
            print ''
            print >> f, ''
            t0 = singleSpikeTrain[i]
            if t0 > run_t:
                break
            i0 = int(round(t0/tstep))
            i1 = i0 + 1
            sel = singleSpikeID[i]
            #print 'v_pre = ', v_pre
            #print >> f, 'v_pre = ', v_pre
            v2, _ = sproceed(cell, v_pre, synList, gList, RList, vecStimList, 0.0, n, sel, trans, trans + tol_t, tstep, '')
            v2 = v2[ntrans:ntrans+tol_nt+1] - vrest
            shiftAdd(biV0,v2,1,i1,tol_nt)
            pyplot.figure(name+'/'+name+'-s'+label1+str(i),figsize=(8,4))
            pyplot.title(label1+str(i))
            pyplot.plot(np.linspace(0,tol_t,tol_nt+1),v2, 'r', lw = 0.5)
            pyplot.savefig(name+'/'+name+'-s'+label1+str(i)+'.png',format='png',bbox_inches='tight',dpi=rdpi)
        
            for j in xrange(0,i):
                print i, '-', j
                print >> f, i, '-', j
                if gList[singleSpikeID[j]] < 0:
                    label2 = 'I'+str(singleSpikeID[j])+'-'
                else:
                    label2 = 'E'+str(singleSpikeID[j])+'-'
                label = label2 + label1
                t00 = singleSpikeTrain[j]
                dt = t0 - t00
        
                if dt < tol_t:
                    idt = int(round(dt/tstep))
                    if singleSpikeID[j] == singleSpikeID[i]:
                        sel = [singleSpikeID[j]]
                        spikes = np.array([[0.0,dt]])
                    else:
                        sel = [singleSpikeID[j],singleSpikeID[i]]
                        spikes = np.array([[0.0],[dt]])
                    print 'synapse:', sel
                    print >> f, 'synapse:', sel
                    v,_ = bproceed0(cell, v_pre, synList, gList, RList, vecStimList, spikes, n, sel, trans, 0.0, trans + tol_t, 0.0, tstep, '')
                    v = v[ntrans+idt:ntrans+tol_nt+1]
                    relit = v.size
                    v = v - vrest
                    assert(relit == tol_nt+1-idt)
                
                    sel = singleSpikeID[j]
                    v1,_ = sproceed(cell, v_pre, synList, gList, RList, vecStimList, 0.0, n, sel, trans, trans + tol_t, tstep, '')
                    v1 = v1[ntrans+idt:ntrans+tol_nt+1]
                    v1 = v1 - vrest
        
                    addv = v2[:relit] + v1
                    bv = v - addv
                    shiftAdd(biV0,bv,1,i1,relit-1)
                    pyplot.figure(name+'/'+name+'-b0'+label+str(i)+'-'+str(j),figsize=(8,4))
                    pyplot.title(label+str(i)+'-'+str(j))
                    t = np.linspace(t0,t0+relit*tstep,relit)
                    pyplot.plot(t,v,'k',lw=2)
                    pyplot.plot(t,v1,'r',lw=1)
                    pyplot.plot(t,v2[:relit],'b',lw=0.5)
                    pyplot.plot(t,addv,':g',lw=2)
                    pyplot.savefig(name+'/'+name+'-b0'+label+str(i)+'-'+str(j)+'.png',format='png',bbox_inches='tight',dpi=rdpi)
                    pyplot.close()

    if bilinear:
        fig = pyplot.figure(name+'/'+name+'-biV',figsize=(8,4))
        ax = fig.add_subplot(211)
        minv = np.amin([np.amin(biV),np.amin(liV),np.amin(simV)])
        maxv = np.amax([np.amax(biV),np.amax(liV),np.amax(simV)])
        t = np.linspace(0,run_t,run_nt+1)
        if linear:
            li, = ax.plot(t,liV[:run_nt+1], lw=0.5, color = '0.75',label = 'linear')
            write_one(name+'/'+name+'-figdata.bin',liV,mode='w')
        bi, = ax.plot(t,biV[:run_nt+1],'g',lw=1.5,label='bilinear')
        if bilinear0:
            bi0, = ax.plot(t,biV0[:run_nt+1],'c',lw=1.0,label='bilinear0')
        si, = ax.plot(t,simV[ntrans:ntrans+run_nt+1],':b',label = 'sim')
        ymin, ymax = ax.get_ylim()
        ax.plot(np.array([noMoreInput_t, noMoreInput_t]),np.array([ymin,ymax]),':k')
        for i in xrange(n):
            if gList[i] > 0:
                for j in xrange(spikeTrain[i].size):
                    ts = spikeTrain[i][j]
                    if ts > run_t:
                        break
                    it = int(np.floor(ts/tstep))
                    vtar = simV[ntrans+it]
                    ax.plot((ts,ts),(minv,vtar),':r')
                    ax.text(ts, minv+(vtar-minv)*0.1*i, str(i)+'-'+str(j), color='red',fontsize=6)
            else:
                for j in xrange(spikeTrain[i].size):
                    ts = spikeTrain[i][j]
                    if ts > run_t:
                        break
                    it = int(np.floor(ts/tstep))
                    vtar = simV[ntrans+it]
                    ax.plot((ts,ts),(minv,vtar),':b')
                    ax.text(ts, minv+(vtar-minv)*0.1*i, str(i)+'-'+str(j), color='blue',fontsize=6)
        write_one(name+'/'+name+'-figdata.bin',biV)
        write_one(name+'/'+name+'-figdata.bin',simV)
        write_one(name+'/'+name+'-figdata.bin',biV0)
        if bilinear0 and linear:
            ax.legend(handles=[si,li,bi,bi0])
        elif linear:
            ax.legend(handles=[si,li,bi])
        else: 
            ax.legend(handles=[si,bi,bi0])
        pyplot.xlabel('t (ms)')
        pyplot.ylabel('V (mV)')

        ax = fig.add_subplot(212)
        if linear:
            li, = ax.plot(t,liV[:run_nt+1]-simV[ntrans:ntrans+run_nt+1], lw=0.5, color = '0.75',label = r'$\Delta$linear')
        if bilinear0:
            bi0, = ax.plot(t,biV0[:run_nt+1]-simV[ntrans:ntrans+run_nt+1],'c', lw=1.0, label=r'$\Delta$bilinear0')
        bi, = ax.plot(t,biV[:run_nt+1]-simV[ntrans:ntrans+run_nt+1],'g', lw=1.5, label=r'$\Delta$bilinear')
        ax.plot(t,np.zeros(run_nt+1),':k')
        if bilinear0 and linear:
            ax.legend(handles=[si,li,bi,bi0])
        elif linear:
            ax.legend(handles=[si,li,bi])
        else: 
            ax.legend(handles=[si,bi,bi0])

        pyplot.savefig(name+'/'+name+'-biV.png',format='png',bbox_inches='tight',dpi=rdpi)
        pyplot.close()
    dbv = simV[ntrans:ntrans+run_nt+1] - biV[:run_nt+1]
    dlv = simV[ntrans:ntrans+run_nt+1] - liV[:run_nt+1]
    db0v = simV[ntrans:ntrans+run_nt+1] - biV0[:run_nt+1]
    f.close()
    if s is not None:
        s.send(np.array([dbv,dlv,db0v,did]))
    else:
        return dbv, dlv, db0v

if __name__ == '__main__':
    #stratum_radiatum 0-78
    #apical dends 0-134
    #basal dends 135-191
    #axon secs 192-199
    #apical trunk [0, 2, 14, 28, 30, 32, 40, 44, 52, 60, 66, 72, 74, 78, 84, 90, 92, 98]
    #name = 't'
    #seed = 14
    f=open('sOlog0','w')
    start = time.time()

    input_start = 0
    run_t = 80 + input_start
    noMoreInput_t = run_t/4
    tol_t = min(240,run_t)
    trans = 110
    rErange = np.array([50,100,150,200,250,300,350,400,450,500])
    rIrange = rErange
    #rIrange = [50,100,150,200,250,300,350,400,450,500]
    nsp = np.array([1,2,3,4,5,6,7,8,9,10])
    tstep = 1.0/10.0
    v0 = -69
    vrest = -70
    rdpi = 160
    presetSpike = 0
    name0 = 'freqLimit4'
    orig = 0
    seed = 1505530645 
    #seed = int(time.mktime(datetime.now().timetuple())) 
    print "seed = ", seed
    np.random.seed(seed)
    linear = 1
    bilinear = 1
    bilinear0 = 1
    loadData = False 
    #locE = np.array([36, 74, 83, 97, 117, 132],dtype='int')
    if presetSpike:
        nE = 2
        nI = 2
        locE = np.random.randint(53,134,nE)
        locI = np.random.randint(0,53,nI)
    else:
        nE = 6
        nI = 2
        locE = np.random.randint(53,134,nE)
        locI = np.random.randint(0,53,nI)
    #gE = 3.0e-4 + np.random.random_sample(locE.size)*0
    gE = 4.e-3*np.ones(locE.size)
    posE = np.random.random_sample(locE.size)
    #locI = np.array([7, 28, 137],dtype='int')
    #gI = 1.5e-4 + np.random.random_sample(locI.size) * 0) * (-1)
    gI = -1.0e-2*np.ones(locI.size)
    posI = np.random.random_sample(locI.size)
    print name0
    print >>f, name0

    if presetSpike:
        ni = nsp.size
        rErange = np.empty(ni)
        rIrange = np.empty(ni)
        freq = nsp.astype('float')/(noMoreInput_t/1000.0)
    else:
        ni = rErange.size
        assert(rErange.size == rIrange.size)
        freq = rErange
    dbv = np.empty((2,ni))
    db0v = np.empty((2,ni))
    dlv = np.empty((2,ni))
    pos = np.concatenate((posE, posI))
    loc = np.concatenate((locE, locI))
    gList = np.concatenate((gE, gI))
    spList = np.empty(ni,dtype=list)
    if not loadData:
        jobs = []
        recv = []
        for i in xrange(ni):
            name = str(i)
            r, send = mp.Pipe(False)
            if presetSpike:
                assert(noMoreInput_t < run_t)
                sp1 = np.linspace(input_start, noMoreInput_t, nsp[i]+1)
                sp1 = sp1[:-1]
                if sp1.size == 1:
                    dsp = noMoreInput_t - input_start
                else:
                    dsp = sp1[1]-sp1[0]
                sp2 = sp1 + dsp/2
                sp3 = sp1 + dsp/2
                sp4 = sp1 + dsp
                spList[i] = [sp1,sp2,sp3,sp4]

            if not os.path.isdir(name):
                os.mkdir('./'+name)
            p = mp.Process(target = subOnly, args = (gList,loc,pos,vrest,v0,presetSpike,trans,run_t,tol_t,tstep,spList[i],seed,name,rErange[i],rIrange[i],linear,bilinear,bilinear0,noMoreInput_t,send,i,rdpi))
            jobs.append(p)
            recv.append(r)
            p.start()
        for p in jobs:
            p.join()
        dbs = [x.recv() for x in recv]
        ind = [dbs[i][3] for i in xrange(ni)]
        argi = np.argsort(ind)
        print "indexes", ind
        for i in xrange(ni):
            dbv[0,i] = np.mean(np.absolute(dbs[argi[i]][0]))
            dlv[0,i] =np.mean(np.absolute(dbs[argi[i]][1]))
            db0v[0,i] = np.mean(np.absolute(dbs[argi[i]][2]))
            dbv[1,i] =  np.std(np.absolute(dbs[argi[i]][0]))
            dlv[1,i] = np.std(np.absolute(dbs[argi[i]][1]))
            db0v[1,i] =  np.std(np.absolute(dbs[argi[i]][2]))
            #dbv[0,i] = np.mean(np.absolute(dbvl))
            #db0v[0,i] = np.mean(np.absolute(db0vl))
            #dlv[0,i] = np.mean(np.absolute(dlvl))
            #dbv[1,i] = np.std(np.absolute(dbvl))
            #db0v[1,i] = np.std(np.absolute(db0vl))
            #dlv[1,i] = np.std(np.absolute(dlvl))

        dbv.dump(name0+'dbv.npy')
        dlv.dump(name0+'dlv.npy')
        db0v.dump(name0+'db0v.npy')
    else:
        dbv = np.load(name0+'-dbv.npy')
        dlv = np.load(name0+'-dlv.npy')
        db0v = np.load(name0+'-db0v.npy')
    fig = pyplot.figure(name0+'-delta', figsize=(8,4))
    ax = fig.add_subplot(111)
    ax.errorbar(freq+0.1,dbv[0,:],dbv[1,:],color='b',ls=':', label=r'$\Delta$bilinear',marker='*',ms=2)
    ax.errorbar(freq+0.2,db0v[0,:],db0v[1,:],color='c',ls=':', label=r'$\Delta$bilinear0',marker='*',ms=2)
    ax.errorbar(freq+0.3,dlv[0,:], dlv[1,:],color='g',ls=':',label=r'$\Delta$linear',marker='*',ms=2)
    pyplot.xlabel('Freq (Hz)')
    pyplot.ylabel(r'$\Delta$V (mV)')
    ax.legend()
    ax.set_yscale('log')
    pyplot.savefig(name0+'-delta.png',format='png',bbox_inches='tight',dpi=900)
    pyplot.close()
    end = time.time()
    print (end - start)/3600.0, ' hrs'
    f.close()
