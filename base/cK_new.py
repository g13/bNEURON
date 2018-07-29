#coding:utf-8
from neuron import h
h.load_file('stdlib.hoc')
h.load_file('stdrun.hoc')
from neuroAlter import prepCell
import multiprocessing as mp
import numpy as np
import sys, os, getopt
from n128 import *
import matplotlib
matplotlib.use('Agg')
#matplotlib.use('Qt5Agg')
from matplotlib import pyplot#, ticker
from sv import bproceed, write_one, gather_and_distribute_results
from neuroAlter import proceed

#### theoretical precision test (without interpolation error and tail cut off #########

def getKv(v1,v2,v,v2p,nstep,nS):
    dv = v-v1-v2 
    r = v2p/v2
    r[np.isnan(r)] = 0.0
    kV = dv*r
    return kV

def plotv(ax,data,i,color,marker):
    target = np.abs(data)
    sign = np.sign(np.average(data))
    if sign == 0.0:
        sign = 1.0
    ax.errorbar(i,sign*np.average(target),np.average(np.std(target,1)),ecolor=color,marker=marker)

def cK(cid, seed, fmt, run_t, tstep, rE, rI, testLinear, testBilinear, plotSindivid, plotBindivid, newSv, sender):
    directory = str(seed)
    g0 = 32.0*6e-4
    #testLinear = True

# initialize cell
    #locE = np.array([60],dtype='int')
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

    run_nt = int(round(run_t/tstep))+1
    n = loc.size

    cell, vSL, sL, _ = prepCell(gList, loc, pos, n, vrest, alphaR)
# population events
    
    rE = rE/1000.0
    rI = rI/1000.0
    rates = np.empty(n)
    rates[:locE.size] = np.random.randn(locE.size)*rE*0.1+rE
    rates[locE.size:] = np.random.randn(locI.size)*rI*0.1+rI
    print rates
    np.random.seed(seed)
    spikeTrain = np.empty(n,dtype=object)
    IDs = np.empty(n,dtype=object)
    activePeriod = run_t/3.0
    for i in xrange(n):
        if rates[i] <= 0.0:
            spikeTrain[i] = np.array([run_t+1])
            continue
        ts = 0
        events = []
        while ts < activePeriod:
            tt = np.random.poisson(1.0/rates[i],1)
            if tt > 0:
                ts += int(round(tt/tstep))*tstep
            events.append(ts)
        print str(i)+'th synapse:', len(events), 'spikes'
        spikeTrain[i] = np.array(events)
    RList = np.zeros((n,2))
    dendVclamp = np.zeros(n)+ 1000
    for i in xrange(n):
        IDs[i] = i + np.zeros(spikeTrain[i].size,dtype=int)
    IDs = np.hstack(IDs)
    spikes = np.hstack(spikeTrain)
    indSorted = np.argsort(spikes)
    sortedID = IDs[indSorted]
    sortedSpikes = spikes[indSorted]
    itsp = np.round(sortedSpikes/tstep).astype(int)

    if testLinear:
        plotLinear = True
    else:
        plotLinear = False
    if testBilinear:
        plotBilinear = True
    else:
        plotBilinear = False
    if not os.path.exists(directory):
        os.mkdir(directory)
    directory = directory + '/'
# run simulation
    v0 = vrest
    trans = 0
    t0 = 0
    oneGo = True
    getDendV = True
    monitorDend = False
    pas = False
    plotSlice = True
    tref = 10.0
    vBack = vThres - 2.0
    simv, _, _, _, dendv = proceed(cell, v0, sL, RList, vSL, spikeTrain, n, trans, run_t, vBack, tref, vThres, oneGo, t0, tstep, loc, pos, dendVclamp, alphaR, getDendV, monitorDend, pas, plotSlice, fign=directory+'simvOnly')

# run linear
    if testLinear:
        liv = np.zeros(run_nt,dtype=float) + vrest
        liv0 = np.zeros(run_nt,dtype=float) + vrest
        liv1 = np.zeros(run_nt,dtype=float) + vrest
        tol = True
        basic_run_t = 1
        for i in xrange(sortedSpikes.size):
            it = itsp[i]
            sel = [sortedID[i]]
            if sortedSpikes[i] >= run_t:
                break
            if i==0:
                v0 = vrest
            else:
                v0 = liv[it]
            ## liv
            tsp = np.array([[0]])
            cpi = False
            sv, _, dendv = bproceed(cell, v0, vThres, sL, gList, vSL, tsp, n, sel, basic_run_t, tstep, '', pos, loc, alphaR, cpi, tol)
            if it + sv.size >= run_nt:
                et = run_nt-it
            else:
                et = sv.size
            sv = sv[:et] - v0
            dendv = dendv[:,:et] - v0
            liv[it:it+et] += sv
            if plotSindivid:
                fign = directory+'liv'+str(i)
                fig = pyplot.figure(fign,figsize=(8,4))
                ax = fig.add_subplot(1,1,1)
                ax.plot(np.arange(sv.size)*tstep,sv)
                ax.plot(np.arange(sv.size)*tstep,dendv.T)
                pyplot.savefig(fign+'.'+fmt,format=fmt,dpi=900)
            ## liv0
            tsp = np.array([[0]])
            cpi = False
            v0 = liv0[it]
            sv, _, dendv = bproceed(cell, vrest, vThres, sL, gList, vSL, tsp, n, sel, basic_run_t, tstep, '', pos, loc, alphaR, cpi, tol)
            if it + sv.size >= run_nt:
                et = run_nt-it
            else:
                et = sv.size
            sv = sv[:et] - vrest
            dendv = dendv[:,:et] - vrest
            liv0[it:it+et] += sv
            if plotSindivid:
                fign = directory+'liv0'+str(i)
                fig = pyplot.figure(fign,figsize=(8,4))
                ax = fig.add_subplot(1,1,1)
                ax.plot(np.arange(sv.size)*tstep,sv)
                ax.plot(np.arange(sv.size)*tstep,dendv.T)
                pyplot.savefig(fign+'.'+fmt,format=fmt,dpi=900)
            ## liv1
            tsp = np.array([[0]])
            v0 = liv1[it]
            print 'liv1'+'-'+str(i), v0
            cpi = True
            cell.set_volt(v0)
            sv, _, dendv = bproceed(cell, vrest, vThres, sL, gList, vSL, tsp, n, sel, basic_run_t, tstep, '', pos, loc, alphaR, cpi, tol)
            if it + sv.size >= run_nt:
                et = run_nt-it
            else:
                et = sv.size
            sel = []
            tsp = np.array([[]])
            leaky_run_t = (et-1)*tstep
            leakyV, _, leakyDendV = bproceed(cell, vrest, vThres, sL, gList, vSL, tsp, n, sel, leaky_run_t, tstep, '', pos, loc, alphaR, cpi)
            sv = sv[:et] - leakyV
            dendv = dendv[:,:et] - leakyDendV
            liv1[it:it+et] += sv
            if plotSindivid:
                fign = directory+'liv1'+str(i)
                fig = pyplot.figure(fign,figsize=(8,4))
                ax = fig.add_subplot(1,1,1)
                ax.plot(np.arange(sv.size)*tstep,sv)
                ax.plot(np.arange(sv.size)*tstep,dendv.T)
                pyplot.savefig(fign+'.'+fmt,format=fmt,dpi=900)

# run bilinear
    if testBilinear:
        tol = True
        biv = np.zeros(run_nt,dtype=float) + vrest
        biv0 = np.zeros(run_nt,dtype=float) + vrest
        biv1 = np.zeros(run_nt,dtype=float) + vrest
        nonv = np.zeros(sortedSpikes.size,dtype=float) + vrest
        nonv1 = np.zeros(sortedSpikes.size,dtype=float) + vrest
        sv = np.empty(sortedSpikes.size,dtype = object)
        sv0 = np.empty(sortedSpikes.size,dtype = object)
        sv1 = np.empty(sortedSpikes.size,dtype = object)
        for i in xrange(sortedSpikes.size):
            if sortedSpikes[i] >= run_t:
                break
            it = itsp[i]
            sel = [sortedID[i]]
            if i==0:
                v0 = vrest
            else:
                v0 = biv[it]
            nonv[i] = v0
            cpi = False
            tsp = np.array([[0]])
            sv[i], _, dendv = bproceed(cell, v0, vThres, sL, gList, vSL, tsp, n, sel, run_t, tstep, '', pos, loc, alphaR, cpi, tol)
            if it + sv[i].size >= run_nt:
                et = run_nt-it
            else:
                et = sv[i].size
            sv[i] = sv[i][:et] - v0
            dendv = dendv[:,:et] - v0
            biv[it:it+et] += sv[i]
            if plotSindivid:
                fign = directory+'biv_'+str(i)
                fig = pyplot.figure(fign,figsize=(8,4))
                ax = fig.add_subplot(1,1,1)
                ax.plot(np.arange(sv[i].size)*tstep,sv[i])
                ax.plot(np.arange(sv[i].size)*tstep,dendv.T)
                pyplot.savefig(fign+'.'+fmt,format=fmt,dpi=900)

            for j in xrange(i):
                jt = itsp[j]
                if jt + sv[j].size < it:
                    continue
                if sortedID[j] == sortedID[i]:
                    sel = [sortedID[j]]
                    tsp = np.array([[0,sortedSpikes[i]-sortedSpikes[j]]])
                else:
                    sel = [sortedID[j],sortedID[i]]
                    tsp = np.array([[0],[sortedSpikes[i]-sortedSpikes[j]]])
                idt = it-jt
                b_run_t = sv[j].size*tstep
                v0 = nonv[j]
                bv, _, dendv = bproceed(cell, v0, vThres, sL, gList, vSL, tsp, n, sel, b_run_t, tstep, '', pos, loc, alphaR, cpi, tol) 
                if jt + bv.size >= run_nt:
                    et = run_nt-jt
                else:
                    et = bv.size
                bv = bv[:et] - v0
                dendv = dendv[:,:et] - v0
                if sv[j].size > et:
                    s1t = et
                else:
                    s1t = sv[j].size

                if newSv:
                    sel = [sortedID[i]]
                    tsp = np.array([[0]])
                    sv2, _, dendv = bproceed(cell, v0, vThres, sL, gList, vSL, tsp, n, sel, (et - idt)*tstep, tstep, '', pos, loc, alphaR, cpi) 
                    sv2 = sv2 - v0
                else:
                    sv2 = sv[i]

                if et - idt > sv2.size:
                    s2t = sv2.size
                else:
                    s2t = et - idt
                print 's2:', it,', ', s2t, '<=', sv2.size
                print 's1:', jt,', ', s1t, '<=', sv[j].size, 'kv start:', idt 
                kv = bv.copy()
                kv[:s1t] -= sv[j][:s1t]
                kv[idt:idt+s2t] -= sv2[:s2t]
                kv[idt+s2t:] = 0
                addv = np.zeros(bv.size)
                addv[:s1t] += sv[j][:s1t]
                if newSv:
                    if sv[i].size < s2t:
                        sv2p = np.append(sv[i],np.zeros(s2t-sv[i].size))
                    else:
                        sv2p = sv[i][:s2t]
                    print 'size: ', sv2.size, sv2p.size, s2t - sv[i].size, sv[i].size
                    kv[idt:idt+s2t] = kv[idt:idt+s2t]*sv2p/sv2[:s2t]
                    addv[idt:idt+s2t] += sv2p
                else:
                    addv[idt:idt+s2t] += sv2[:s2t]
                
                biv[jt:jt+et] += kv

                if plotBindivid:
                    fign = directory+'biv'+str(i)+'-'+str(j)
                    fig = pyplot.figure(fign,figsize=(8,4))
                    ax = fig.add_subplot(1,1,1)
                    t = np.arange(bv.size)*tstep
                    ax.plot(t,bv,'b')
                    #ax.plot(t,dendv.T)
                    ax.plot(t,kv,':k')
                    ax.plot(t,addv,':g')
                    ax.plot(t[:s1t],sv[j][:s1t],':r')
                    ax.plot(t[idt:idt+s2t],sv2[:s2t],':b')
                    pyplot.savefig(fign+'.'+fmt,format=fmt,dpi=900)

            v0 = biv0[it]
            sel = [sortedID[i]]
            tsp = np.array([[0]])
            sv0[i], _, dendv = bproceed(cell, vrest, vThres, sL, gList, vSL, tsp, n, sel, run_t, tstep, '', pos, loc, alphaR, cpi, tol)
            if it + sv0[i].size >= run_nt:
                et = run_nt-it
            else:
                et = sv0[i].size
            sv0[i] = sv0[i][:et] - vrest
            dendv = dendv[:,:et] - vrest
            biv0[it:it+et] += sv0[i]
            if plotSindivid:
                fign = directory+'biv0_'+str(i)
                fig = pyplot.figure(fign,figsize=(8,4))
                ax = fig.add_subplot(1,1,1)
                ax.plot(np.arange(sv0[i].size)*tstep,sv0[i])
                ax.plot(np.arange(sv0[i].size)*tstep,dendv.T)
                pyplot.savefig(fign+'.'+fmt,format=fmt,dpi=900)

            for j in xrange(i):
                jt = itsp[j]
                if jt + sv0[j].size < it:
                    continue
                if sortedID[j] == sortedID[i]:
                    sel = [sortedID[j]]
                    tsp = np.array([[0,sortedSpikes[i]-sortedSpikes[j]]])
                else:
                    sel = [sortedID[j],sortedID[i]]
                    tsp = np.array([[0],[sortedSpikes[i]-sortedSpikes[j]]])
                idt = it-jt
                b_run_t = sv0[j].size*tstep
                bv, _, dendv = bproceed(cell, vrest, vThres, sL, gList, vSL, tsp, n, sel, b_run_t, tstep, '', pos, loc, alphaR, cpi, tol) 
                if jt + bv.size >= run_nt:
                    et = run_nt-jt
                else:
                    et = bv.size
                bv = bv[:et] - vrest
                dendv = dendv[:,:et] - vrest
                if sv0[j].size > et:
                    s1t = et
                else:
                    s1t = sv0[j].size

                sv2 = sv0[i]

                if et - idt > sv2.size:
                    s2t = sv2.size
                else:
                    s2t = et - idt
                print 's2:', it,', ', s2t, '<=', sv2.size
                print 's1:', jt,', ', s1t, '<=', sv0[j].size, 'kv start:', idt 
                kv = bv.copy()
                kv[:s1t] -= sv0[j][:s1t]
                kv[idt:idt+s2t] -= sv2[:s2t]
                addv = np.zeros(bv.size)
                addv[:s1t] += sv0[j][:s1t]
                addv[idt:idt+s2t] += sv2[:s2t]
                
                biv0[jt:jt+et] += kv

                if plotBindivid:
                    fign = directory+'biv0_'+str(i)+'-'+str(j)
                    fig = pyplot.figure(fign,figsize=(8,4))
                    ax = fig.add_subplot(1,1,1)
                    t = np.arange(bv.size)*tstep
                    ax.plot(t,bv,'b')
                    #ax.plot(t,dendv.T)
                    ax.plot(t,kv,':k')
                    ax.plot(t,addv,':g')
                    ax.plot(t[:s1t],sv0[j][:s1t],':r')
                    ax.plot(t[idt:idt+s2t],sv2[:s2t],':b')
                    pyplot.savefig(fign+'.'+fmt,format=fmt,dpi=900)

            sel = [sortedID[i]]
            if i==0:
                v0 = vrest
            else:
                v0 = biv1[it]
            nonv1[i] = v0
            cpi = True
            cell.set_volt(v0)
            tsp = np.array([[0]])
            sv1[i], _, dendv = bproceed(cell, v0, vThres, sL, gList, vSL, tsp, n, sel, run_t, tstep, '', pos, loc, alphaR, cpi, tol)
            if it + sv1[i].size >= run_nt:
                et = run_nt-it
            else:
                et = sv1[i].size
            leaky_run_t = (et-1)*tstep
            sel = []
            tsp = np.array([[]])
            cell.set_volt(v0)
            leakyV, _, leakyDendV = bproceed(cell, v0, vThres, sL, gList, vSL, tsp, n, sel, leaky_run_t, tstep, '', pos, loc, alphaR, cpi)
            sv1[i] = sv1[i][:et] - leakyV
            dendv = dendv[:,:et] - leakyDendV
            biv1[it:it+et] += sv1[i]
            if plotSindivid:
                fign = directory+'biv1_'+str(i)
                fig = pyplot.figure(fign,figsize=(8,4))
                ax = fig.add_subplot(1,1,1)
                ax.plot(np.arange(sv1[i].size)*tstep,sv1[i])
                ax.plot(np.arange(sv1[i].size)*tstep,dendv.T)
                pyplot.savefig(fign+'.'+fmt,format=fmt,dpi=900)

            for j in xrange(i):
                jt = itsp[j]
                if jt + sv1[j].size < it:
                    continue
                if sortedID[j] == sortedID[i]:
                    sel = [sortedID[j]]
                    tsp = np.array([[0,sortedSpikes[i]-sortedSpikes[j]]])
                else:
                    sel = [sortedID[j],sortedID[i]]
                    tsp = np.array([[0],[sortedSpikes[i]-sortedSpikes[j]]])
                idt = it-jt
                b_run_t = sv1[j].size*tstep
                v0 = nonv1[j]
                cell.set_volt(v0)
                bv, _, dendv = bproceed(cell, v0, vThres, sL, gList, vSL, tsp, n, sel, b_run_t, tstep, '', pos, loc, alphaR, cpi, tol) 
                if jt + bv.size >= run_nt:
                    et = run_nt-jt
                else:
                    et = bv.size
                leaky_run_t = (et-1)*tstep
                sel = []
                tsp = np.array([[]])
                cell.set_volt(v0)
                leakyV, _, leakyDendV = bproceed(cell, v0, vThres, sL, gList, vSL, tsp, n, sel, leaky_run_t, tstep, '', pos, loc, alphaR, cpi)
                bv = bv[:et] - leakyV
                dendv = dendv[:,:et] - leakyDendV
                if sv1[j].size > et:
                    s1t = et
                else:
                    s1t = sv1[j].size

                if newSv:
                    sel = [sortedID[i]]
                    tsp = np.array([[0]])
                    v0 = bv[idt]
                    cell.set_volt(v0)
                    sv2, _, dendv = bproceed(cell, v0, vThres, sL, gList, vSL, tsp, n, sel, (et - idt)*tstep, tstep, '', pos, loc, alphaR, cpi) 
                    leaky_run_t = (sv2.size-1)*tstep
                    sel = []
                    tsp = np.array([[]])
                    cell.set_volt(v0)
                    leakyV, _, leakyDendV = bproceed(cell, v0, vThres, sL, gList, vSL, tsp, n, sel, leaky_run_t, tstep, '', pos, loc, alphaR, cpi)
                    sv2 = sv2 - leakyV
                else:
                    sv2 = sv1[i]

                if et - idt > sv2.size:
                    s2t = sv2.size
                else:
                    s2t = et - idt
                print 's2:', it,', ', s2t, '<=', sv2.size
                print 's1:', jt,', ', s1t, '<=', sv1[j].size, 'kv start:', idt 
                kv = bv.copy()
                kv[:s1t] -= sv1[j][:s1t]
                kv[idt:idt+s2t] -= sv2[:s2t]
                addv = np.zeros(bv.size)
                addv[:s1t] += sv1[j][:s1t]
                addv[idt:idt+s2t] += sv2[:s2t]
                
                biv1[jt:jt+et] += kv

                if plotBindivid:
                    fign = directory+'biv1_'+str(i)+'-'+str(j)
                    fig = pyplot.figure(fign,figsize=(8,4))
                    ax = fig.add_subplot(1,1,1)
                    t = np.arange(bv.size)*tstep
                    ax.plot(t,bv,'b')
                    #ax.plot(t,dendv.T)
                    ax.plot(t,kv,':k')
                    ax.plot(t,addv,':g')
                    ax.plot(t[:s1t],sv1[j][:s1t],':r')
                    ax.plot(t[idt:idt+s2t],sv2[:s2t],':b')
                    pyplot.savefig(fign+'.'+fmt,format=fmt,dpi=900)
# plot and compare
    t = np.arange(run_nt)*tstep
    fign = 'compare'+str(seed)
    fig = pyplot.figure(fign,figsize=(8,4))
    ax1 = fig.add_subplot(2,1,1)
    ax1.plot(t,simv,'k')
    if plotLinear:
        ax1.plot(t,liv,'r')
        ax1.plot(t,liv0,'g')
        ax1.plot(t,liv1,'m')
    if plotBilinear:
        ax1.plot(t,biv,'b')
        ax1.plot(t,biv0,'c')
        ax1.plot(t,biv1,'y')

    ax2 = fig.add_subplot(2,1,2)
    if plotLinear:
        ax2.plot(t,simv-liv,':r')
        ax2.plot(t,simv-liv0,':g')
        ax2.plot(t,simv-liv1,':m')
    if plotBilinear:
        ax2.plot(t,simv-biv,':b')
        ax2.plot(t,simv-biv0,':c')
        ax2.plot(t,simv-biv1,':y')
    
    pyplot.savefig(fign+'.'+fmt,format=fmt,dpi=900)
    datafn = directory+'data.bin'
    write_one(datafn,np.array([np.int(run_nt)]),'wb') 
    write_one(datafn,simv) 
    write_one(datafn,np.array([np.int(testLinear)])) 
    if testLinear:
        write_one(datafn,liv) 
        write_one(datafn,liv0) 
        write_one(datafn,liv1) 
    write_one(datafn,np.array([np.int(testBilinear)])) 
    if testBilinear:
        write_one(datafn,biv) 
        write_one(datafn,biv0) 
        write_one(datafn,biv1) 
    print cid, 'finished'
    sender.send(np.array([cid,datafn],dtype=object))

if __name__ == '__main__':
    fmt = 'png'
    theme = 'oldSv'
    generateData = True 
    seed0 = 6723412
    run_t = 80
    tstep = 1.0/10.0
    nproc = 2
    rE = 20
    rI = 10
    np.random.seed(seed0)
    seed = np.random.randint(234535,945678,nproc)

    testLinear = False
    plotSindivid = True 
    testBilinear = True
    newSv = True
    plotBindivid = True 

    run_nt = int(round(run_t/tstep))+1
    dataFn =np.empty(nproc,dtype=object)
    simv = np.empty((nproc,run_nt))
    liv = np.empty((nproc,run_nt))
    liv0 = np.empty((nproc,run_nt))
    liv1 = np.empty((nproc,run_nt))
    biv = np.empty((nproc,run_nt))
    biv0 = np.empty((nproc,run_nt))
    biv1 = np.empty((nproc,run_nt))
    if generateData:
        jobs = []
        receivers = []
        for i in xrange(nproc):
            receiver, sender = mp.Pipe(False)
            job = mp.Process(target=cK, args = (i,seed[i],fmt,run_t,tstep,rE,rI,testLinear,testBilinear,plotSindivid,plotBindivid,newSv,sender))
            jobs.append(job)
            receivers.append(receiver)
            job.start()
        gather_and_distribute_results(receivers, jobs, cK.__name__ + str(seed0), nproc, 0, dataFn)
        print dataFn
    else:
        for i in xrange(nproc):
            dataFn[i] = str(seed[i])+'/data.bin'
    for i in xrange(nproc):
        with open(dataFn[i],'rb') as datafile:
            tmp = np.fromfile(datafile, dtype=int,count=1)
            simv[i,:] = np.fromfile(datafile, dtype=float, count=run_nt)
            testLinear = np.fromfile(datafile, dtype=int, count=1)
            if testLinear:
                liv[i,:] = np.fromfile(datafile, dtype=float, count=run_nt)
                liv0[i,:] = np.fromfile(datafile, dtype=float, count=run_nt)
                liv1[i,:] = np.fromfile(datafile, dtype=float, count=run_nt)
            testBilinear = np.fromfile(datafile, dtype=int, count=1)
            if testBilinear:
                biv[i,:] = np.fromfile(datafile, dtype=float, count=run_nt)
                biv0[i,:] = np.fromfile(datafile, dtype=float, count=run_nt)
                biv1[i,:] = np.fromfile(datafile, dtype=float, count=run_nt)

    theme = '-'+theme+'-'+str(seed0)
    fign = 'StatsError'+theme
    fig = pyplot.figure(fign,figsize=(8,4))
    ax = fig.add_subplot(1,1,1)
    t = np.arange(run_nt)*tstep
    marker = '*'
    color = 'r'
    plotv(ax,liv-simv,1,color,marker)
    color = 'g'
    plotv(ax,liv0-simv,2,color,marker)
    color = 'm'
    plotv(ax,liv1-simv,3,color,marker)
    color = 'b'
    plotv(ax,biv-simv,4,color,marker)
    color = 'c'
    plotv(ax,biv0-simv,5,color,marker)
    color = 'y'
    plotv(ax,biv1-simv,6,color,marker)
    pyplot.savefig(fign+'.'+fmt,format=fmt,dpi=900)
