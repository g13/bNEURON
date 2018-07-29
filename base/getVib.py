#coding:utf-8
from neuron import h
#from neuroAlter import *
h.load_file('stdlib.hoc')
h.load_file('stdrun.hoc')
import sys, getopt, os, time
import numpy as np
import multiprocessing as mp
from sv import *
from bv import *

def getVib(argv):
    seed = 231278
    rdpi = 300
    tstep = 1.0/10.0
    #==================
    run_t = 340.0
    #==================
    #run_t = 100.0
    
    dtRange = np.array([0,4,8,15,20,25,30,50,70,140,210],dtype='float')
    #===========================================================
    #dtRangeE = np.array([6,8,10,12,15,17,19,21,23,26,30],dtype='float')
    #dtRangeI = np.array([0,0.5,1,1.5,2,4,40,50,60,70,80,90,110,125,150,175],dtype='float')
    #dtRange = np.sort(np.hstack([dtRangeE,dtRangeI]))
    #===========================================================
    #dtRange = np.array([0,30,60],dtype='float')

    assert(dtRange[0] == 0)
    assert(dtRange[-1] < run_t)
    np.random.seed(seed)
    #===========================================================
    locE = np.array([60, 72, 78, 84, 90, 98],dtype='int')
    locI = np.array([14, 28, 30],dtype='int')
    #locE = np.array([79, 82, 83, 98, 120, 124],dtype='int')
    #locI = np.array([14, 28, 40],dtype='int')
    vThres = -54.0
    vRange = np.array([-74,-70,-67,-65,-63,-62,-61,-60,-59,-58],dtype='double')
    #===========================================================
    #vRange = np.array([-74,-70,-65,-61],dtype='double')
    #locE = np.array([72,79],dtype='int')
    #locI = np.array([14,40],dtype='int')

    g0 = 32.0*1e-3
    gE = np.array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0]) * g0
    gI = -g0*np.array([10.0, 10.0, 10.0])
    #gE = np.array([0.6, 0.6, 0.2, 0.6, 0.15, 0.6]) * g0
    #gI = -g0*np.array([6.0, 10.0, 8.0])
    posE = np.array([0.3,0.3,0.9,0.6,0.4,0.2])
    posI = np.array([0.7,0.2,0.5])
    gE = gE[:locE.size]
    gI = gI[:locI.size]
    posE = posE[:locE.size]
    posI = posI[:locI.size]
    # defaults
    plotData = True 
    singleOnly = False 
    vrest = -70
    fmt = 'png'
    assert(vrest in vRange)
    usage = 'getVib.py --theme=<theme name> --vid=<vid> -p<plotData> -s<single only>// all the items in lists should be separated by comma without space'
    try:
        opts, args = getopt.getopt(argv,'psh',['theme=','vid=','trans=','fmt='])
    except getopt.GetoptError:
        print "error reading args"
        print argv
        print usage
        sys.exit(2)
    print opts
    opt0 = [opt[0] for opt in opts]
    if '-h' in opt0:
        print usage
        sys.exit()
    if '--vid' not in opt0 and '-s' not in opt0:
        print 'vid must be specified or -s (singleOnly) need to be used'
        sys.exit(2)
    if '--theme' not in opt0:
        print 'theme name must be specified'
        sys.exit(2)
    for opt, arg in opts:
        if opt == '--vid':
            vid = int(arg)
        elif opt == '--theme':
            theme = arg
        elif opt == '-p':
            plotData = True
        elif opt == '-s':
            singleOnly = True
        elif opt == '--fmt':
            fmt = arg

    alphaR = True
    sim_t = run_t
    tol_t = dtRange[-1]
    tol_nt = int(round(tol_t/tstep)) + 1
    run_nt = int(round(run_t/tstep)) + 1
    nv = vRange.size
    ndt = dtRange.size
    idtRange = np.round(dtRange/tstep).astype('int')
    
    t = np.arange(run_nt)*tstep
    pos = np.concatenate((posE, posI))
    loc = np.concatenate((locE, locI))
    gList = np.concatenate((gE, gI))
    nE = locE.size
    nI = locI.size
    n = nE + nI
    RList = np.zeros((n,2))
    directory = '.' 
    
    # sv
    if singleOnly:
        write_one(directory+'/'+theme+'-vRange.bin',vRange,'wb')
        write_one(directory+'/'+theme+'-dtRange.bin',dtRange,'wb')
        p = np.array([tstep,run_t,n])
        write_one(directory+'/'+theme+'-p.bin',p,'wb')
        write_one(directory+'/'+theme+'-R.bin',loc,'wb')
        write_one(directory+'/'+theme+'-R.bin',pos,'ab')
        write_one(directory+'/'+theme+'-R.bin',gList,'ab')

        plotSv(vRange,gList,loc,pos,vrest,vThres,run_t,tstep,fmt,theme,True)
        #plotSv(gList,loc,pos,vrest,vThres,run_t,tstep,fmt,theme,False)
        return

    datafn = 'singlets-' + theme
    tmp = float(1.0)
    floatSize = tmp.__sizeof__()
    with open(datafn+'.bin','rb') as datafile:
        datafile.seek(vid*run_nt*n*floatSize)
        sv = np.fromfile(datafile,dtype=float,count=n*run_nt)
        sv = sv.reshape(n,run_nt)
    fign = 'test_sv.png'
    fig = pyplot.figure('test_sv',figsize=(8,4))
    ax = fig.add_subplot(1,1,1)
    t = np.arange(run_nt)*tstep
    ax.plot(t,sv.T)
    pyplot.savefig(fign, format='png',bbox_inches='tight',dpi=900)

    theme = theme + '-V' + str(vid)
    v0 = vRange[vid]
    # bv
    bfire = np.zeros((ndt,n,n),dtype=int)
    for idt in xrange(ndt):
        print " idt = ", idt
        dt = dtRange[idt]
        kv = np.zeros((n,n,run_nt))
        jobs = []
        receivers = []
        for i in xrange(n):
            receiver, sender = mp.Pipe(False)
            job = mp.Process(target = getI, args = (i,dt,idt,dtRange,idtRange,iextraRange,ndt,vid,gList,pos,v0,cell,synList,vecStimList,n,trans,ntrans,run_t,tol_t,tstep,run_nt,tol_nt,leakyV,extralet,v1,send,savePlot,plotData,directory,rdpi,alphaR,nc))
            jobs.append(job)
            receivers.append(receiver)
            job.start()
        gather_and_distribute_results(receivers, jobs, getI.__name__ + str(vid), n, 0, kV, k, tMax, sFire)
        print 'gathering kv idt ', idt
        result = np.array([x.recv() for x in recv])
        for p in jobs:
            p.join()
        print 'joined kv idt ', idt
        ind = [result[i][2] for i in xrange(n)]
        argi = np.argsort(ind)
        print "distribute kv idt ", idt
        for i in xrange(n):
            ii = argi[i]
            kv[i,:,:,:] = result[ii][0]
            bfire[idt,i,:,:] = result[ii][1]
                    
        write_one(directory+'/'+theme+'.bin',kv,'ab')
    write_one(directory+'/'+theme+'.bin',bfire,'ab')

if __name__ == '__main__':
    getVib(sys.argv[1:])
