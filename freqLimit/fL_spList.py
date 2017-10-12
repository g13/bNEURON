import numpy as np
import os, time, sys
from neuroAlter import *
from getPSP import write_one
from shutil import copy
from datetime import datetime
from matplotlib import pyplot
import matplotlib
matplotlib.use('Agg')

f=open('sOlog','w')

class receptor(object):
    def __init__(self, a, b, c, d, cm):
        self.a = a
        self.b = b
        self.c = c
        self.d = d
        self.cm = cm
        self.Rtau = 1 / (cm*a + b)
        self.Rinf = self.Rtau * cm*a 
        self.expr = np.exp(-c/self.Rtau)
        print "decayTau = ", 1/self.b
        print >> f, "decayTau = ", 1/self.b
        print "Rtau = ", self.Rtau
        print >> f, "Rtau = ", self.Rtau
        print "Rinf = ", self.Rinf
        print >> f, "Rinf = ", self.Rinf

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

def subOnly(gList,loc,pos,vrest,v0,presetSpike,trans,run_t,tol_t,tstep,spikeTrain,seed,rE,rI,name,linear,bilinear,bilinear0,noMoreInput_t,rdpi):
    
    nE = sum([x>0 for x in gList])
    nI = sum([x<=0 for x in gList])
    n = nI + nE
    
    cell, vecStimList, synList = prepCell(gList, loc, pos, n, vrest)
    
    tol_nt = int(tol_t/tstep)
    ntrans = int(trans/tstep)
    run_nt = int(run_t/tstep)
    
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
    pyplot.figure(name+'-'+'leakyV',figsize=(8,4))
    pyplot.plot(np.linspace(0,run_t,run_nt+1),liV[:run_nt+1])
    pyplot.savefig(name+'-leakyV.png',format='png',bbox_inches='tight',dpi=rdpi)
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
        
            pyplot.figure(name+'-s'+label+str(i),figsize=(8,4))
            pyplot.title(label+str(i))
            pyplot.plot(np.linspace(0,tol_t,tol_nt+1),v,'g',lw = 1.5)
            pyplot.savefig(name+'-s'+label+str(i)+'.png',format='png',bbox_inches='tight',dpi=rdpi)
            pyplot.close()
            shiftAdd(liV,v,1,i1,tol_nt)
                
        pyplot.figure(name+'-liV',figsize=(8,4))
        pyplot.plot(np.linspace(0,run_t,run_nt+1),simV[ntrans:ntrans+run_nt+1], lw = 2)
        pyplot.plot(np.linspace(0,run_t,run_nt+1),liV[:run_nt+1], lw = 1)
        pyplot.savefig(name+'-liV.png',format='png',bbox_inches='tight',dpi=rdpi)
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
            pyplot.figure(name+'-s'+label1+str(i),figsize=(8,4))
            pyplot.title(label1+str(i))
            pyplot.plot(np.linspace(0,tol_t,tol_nt+1),v2,'b',lw = 1.0)
            pyplot.savefig(name+'-s'+label1+str(i)+'.png',format='png',bbox_inches='tight',dpi=rdpi)
        
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
                    pyplot.figure(name+'-b'+label+str(i)+'-'+str(j),figsize=(8,4))
                    pyplot.title(label+str(i)+'-'+str(j))
                    t = np.linspace(t0,t0+relit*tstep,relit)
                    pyplot.plot(t,v,'k',lw=2)
                    pyplot.plot(t,v1,'r', lw=1)
                    pyplot.plot(t,v2[:relit],'b',lw = 0.5)
                    pyplot.plot(t,addv,':g',lw = 2)
                    pyplot.savefig(name+'-b'+label+str(i)+'-'+str(j)+'.png',format='png',bbox_inches='tight',dpi=rdpi)
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
            pyplot.figure(name+'-s'+label1+str(i),figsize=(8,4))
            pyplot.title(label1+str(i))
            pyplot.plot(np.linspace(0,tol_t,tol_nt+1),v2, 'r', lw = 0.5)
            pyplot.savefig(name+'-s'+label1+str(i)+'.png',format='png',bbox_inches='tight',dpi=rdpi)
        
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
                    pyplot.figure(name+'-b0'+label+str(i)+'-'+str(j),figsize=(8,4))
                    pyplot.title(label+str(i)+'-'+str(j))
                    t = np.linspace(t0,t0+relit*tstep,relit)
                    pyplot.plot(t,v,'k',lw=2)
                    pyplot.plot(t,v1,'r',lw=1)
                    pyplot.plot(t,v2[:relit],'b',lw=0.5)
                    pyplot.plot(t,addv,':g',lw=2)
                    pyplot.savefig(name+'-b0'+label+str(i)+'-'+str(j)+'.png',format='png',bbox_inches='tight',dpi=rdpi)
                    pyplot.close()

    simV = simV[ntrans:ntrans+run_nt+1]
    biV = biV[:run_nt+1]
    liV = liV[:run_nt+1]
    biV0 = biV0[:run_nt+1]
    if bilinear:
        fig = pyplot.figure(name+'-biV',figsize=(8,4))
        ax = fig.add_subplot(211)
        minv = np.amin([np.amin(biV),np.amin(liV),np.amin(simV)])
        maxv = np.amax([np.amax(biV),np.amax(liV),np.amax(simV)])
        t = np.linspace(0,run_t,run_nt+1)
        if linear:
            li, = ax.plot(t,liV, lw=0.5, color = '0.75',label = 'linear')
            write_one(name+'-figdata.bin',liV,mode='w')
        bi, = ax.plot(t,biV,'g',lw=1.5,label='bilinear')
        if bilinear0:
            bi0, = ax.plot(t,biV0,'c',lw=1.0,label='bilinear0')
        si, = ax.plot(t,simV,':b',label = 'sim')
        ymin, ymax = ax.get_ylim()
        ax.plot(np.array([noMoreInput_t, noMoreInput_t]),np.array([ymin,ymax]),':k')
        for i in xrange(n):
            if gList[i] > 0:
                for j in xrange(spikeTrain[i].size):
                    ts = spikeTrain[i][j]
                    if ts > run_t:
                        break
                    it = int(np.floor(ts/tstep))
                    vtar = simV[it]
                    ax.plot((ts,ts),(minv,vtar),':r')
                    ax.text(ts, minv+(vtar-minv)*0.1*i, str(i)+'-'+str(j), color='red',fontsize=6)
            else:
                for j in xrange(spikeTrain[i].size):
                    ts = spikeTrain[i][j]
                    if ts > run_t:
                        break
                    it = int(np.floor(ts/tstep))
                    vtar = simV[it]
                    ax.plot((ts,ts),(minv,vtar),':b')
                    ax.text(ts, minv+(vtar-minv)*0.1*i, str(i)+'-'+str(j), color='blue',fontsize=6)
        write_one(name+'-figdata.bin',biV)
        write_one(name+'-figdata.bin',simV)
        write_one(name+'-figdata.bin',biV0)
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
            li, = ax.plot(t,liV-simV, lw=0.5, color = '0.75',label = r'$\Delta$linear')
        if bilinear0:
            bi0, = ax.plot(t,biV0-simV,'c', lw=1.0, label=r'$\Delta$bilinear0')
        bi, = ax.plot(t,biV-simV,'g', lw=1.5, label=r'$\Delta$bilinear')
        ax.plot(t,np.zeros(run_nt+1),':k')
        if bilinear0 and linear:
            ax.legend(handles=[si,li,bi,bi0])
        elif linear:
            ax.legend(handles=[si,li,bi])
        else: 
            ax.legend(handles=[si,bi,bi0])

        pyplot.savefig(name+'-biV.png',format='png',bbox_inches='tight',dpi=rdpi)
        pyplot.close()
    return simV, biV, liV, biV0

if __name__ == '__main__':
    #stratum_radiatum 0-78
    #apical dends 0-134
    #basal dends 135-191
    #axon secs 192-199
    #apical trunk [0, 2, 14, 28, 30, 32, 40, 44, 52, 60, 66, 72, 74, 78, 84, 90, 92, 98]

    start = time.time()

    run_t = 80
    noMoreInput_t = 25
    tol_t = min(240,run_t)
    trans = 110
    tstep = 1.0/10.0
    v0 = -67
    vrest = -70
    rdpi = 300 
    name0 = 'fL_spList'
    seed = 231278
    #seed = int(time.mktime(datetime.now().timetuple())) 
    print "seed = ", seed
    np.random.seed(seed)
    locE = np.array([52, 60, 72, 78, 84, 98],dtype='int')
    locI = np.array([14, 28, 40],dtype='int')
    gE = (1e-5 + np.random.random_sample(locE.size) * (1e-1-1e-5)) * (1.0/2.0)
    gI = (1e-5 + np.random.random_sample(locI.size) * (1e-1-1e-5)) * (-1.0/2.0)
    posE = np.random.random_sample(locE.size)
    posI = np.random.random_sample(locI.size)

    pos = np.concatenate((posE, posI))
    loc = np.concatenate((locE, locI))
    spList = np.array([np.array([ 49.7]), np.array([  0.5,  20.3,  39.8]), np.array([  2.1,  16.9,  52.1]), np.array([  3.8,  14.7,  15.9,  51.2]), np.array([  5.7,  13.4,  18.4,  23.8,  40.1]), np.array([  7.8,  12.9,  22.8,  41.6]), np.array([ 41.]), np.array([ 52.4]), np.array([ 66.2])])
    for i in xrange(loc.size):
        spList[i] = spList[i][:-1]
    linear = 1
    bilinear = 1
    bilinear0 = 1
    gList = np.concatenate((gE, gI))
    simV, biV, liV, biV0 = subOnly(gList,loc,pos,vrest,v0,1,trans,run_t,tol_t,tstep,spList,seed,0,0,name0,linear,bilinear,bilinear0,noMoreInput_t,rdpi)

    simV.dump(name0+'simV.npy')
    biV.dump(name0+'biV.npy')
    liV.dump(name0+'liV.npy')
    biV0.dump(name0+'biV0.npy')
    end = time.time()
    print (end - start)/3600.0, ' hrs'
f.close()
