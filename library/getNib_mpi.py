from neuroAlter import *
import struct, sys, getopt, os, time
import numpy as np
import multiprocessing as mp
from neuron import h
from shutil import copy
h.load_file('stdlib.hoc')
h.load_file('stdrun.hoc')

figfmt = 'png'
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
def write_one(data_filename,data,mode='ab'):
    with open(data_filename,mode) as data_file:
        data.tofile(data_file)

def getK(v1,v2,v,nstep,nS):
    K = np.empty(nstep)
    vp = np.empty((nstep,nS))
    for i in xrange(nstep):
        dv = v[i,:]-v1[i,:]-v2[i,:] 
        A = v1[i,:]*v2[i,:]
        At = np.transpose(A)
        K[i] = np.dot(At,dv)/np.dot(At,A) 
        if np.isnan(K[i]):
            K[i] = 0
        vp[i,:] = v1[i,:]+v2[i,:]+A*K[i]
    return K, vp

def getNear(seq, target):
    n = seq.size 
    if target < seq[0]:
        istart = 1
        jnext = 0
        ratio = (target-seq[1])/(seq[0]-seq[1])
    else:
        if target >= seq[n-1]:
            istart = n-2 
            jnext = n-1
            ratio = (target-seq[istart])/(seq[jnext]-seq[istart]);
        else:
            for i in xrange(n-1):
                if target >= seq[i] and target < seq[i+1]:
                    istart = i; jnext = i+1;
                    ratio = (target-seq[istart])/(seq[jnext]-seq[istart]);
                    break;
    return istart, jnext, ratio

def getJ(i,j,dt,idt,dtRange,idtRange,iextraRange,ndt,vid,gList,pos,v0,cell,sL,vSL,n,trans,ntrans,run_t,tol_t,tstep,run_nt,tol_nt,leakyV,extralet,v1,s,savePlot,plotData,directory,rdpi):
    print "j = ", j
    bfire = np.empty(ndt)
    kv = np.zeros((ndt,run_nt))
    if i == j:
        if idt == 0:
            s.send(np.array([kv,bfire,j]))
            return
        else:
            spikes = np.array([[]])
            sel = [i]
    else:
        spikes = np.array([np.array([]),np.array([])])
        sel = [i, j]
    if plotData:
        dtfigname = 'kv-v'+str(vid)+'-i'+str(i)+'-j'+str(j)+'-dt'+str(idt)
        dtfig = pyplot.figure(dtfigname,figsize = (8,4))
    for jdt in xrange(idt,ndt):
        print " jdt = ", jdt
        dtt = dtRange[jdt]
        RList = np.zeros((n,2))
        if i==j:
            rl0 = getGH(dt,gList[i])
            RList[i,:] = getGH(dtt-dt,gList[i],rl0[0],rl0[1])
        else:
            RList[i,:] = getGH(dtt,gList[i])
            RList[j,:] = getGH(dtt-dt,gList[j])
        #if jdt == ndt-1:
        relit = run_nt - idtRange[jdt]
        relt = run_t-dtt
        v,fired = bproceed(cell, v0, sL, gList, RList, vSL, spikes, n, sel, trans, trans + relt , 0, tstep, '')
        v = v[ntrans:ntrans+relit]
        #else:
        #    relit = tol_nt - idtRange[jdt]
        #    relt = tol_t-dtt
        #    v,fired = bproceed(cell, v0, sL, gList, RList, vSL, spikes, n, sel, trans, trans + relt, 0, tstep, '')
        #    v = v[ntrans:ntrans+relit]
        bfire[jdt] = fired
        assert(dtt>=dt)
        print relit, ' == ', v.size
        assert(v.size == relit)
        v = v - leakyV[:relit]
        print "max V = ", np.amax(v)
        print "max it = ", np.argmax(v)
        dtij = idtRange[jdt] - idtRange[idt]
        idtij = np.argwhere(idtRange - dtij == 0)
        assert(idtij.size==1 or idtij.size==0)
        if idtij.size == 1:
            print idtij, " I skipped a singlet simulation"
            v2 = v1[idtij[0,0],j,dtij:dtij+relit]
        else:
            idtij = np.argwhere(iextraRange - dtij == 0)
            print idtij, " I skipped a extralet simulation"
            assert(idtij.size==1)
            v2 = extralet[idtij[0,0],j,dtij:dtij+relit]

        assert(v2.size == relit)
        addv = v1[jdt,i,idtRange[jdt]:idtRange[jdt]+relit] + v2
        kvtmp = v-addv
        kv[jdt,idtRange[jdt]:idtRange[jdt]+relit] = kvtmp
        print "max kv = ", np.amax(np.absolute(kvtmp))
        if plotData:
            pyplot.figure(dtfigname)
            t = np.arange(run_nt)*tstep
            tSel = np.arange(idtRange[jdt],idtRange[jdt]+relit)
            ax1 = dtfig.add_subplot(111)         
            ax1.plot(t[tSel],v,'c',lw=3)
            ax1.plot(t[tSel],v1[jdt,i,tSel],'r',lw=2)
            ax1.plot(t[tSel],v2,'b',lw=1)
            ax1.plot(t[tSel],addv,':g',lw=1)
            ax1.plot(t[tSel],kvtmp,':k')
    if savePlot and plotData:
        pyplot.figure(dtfigname)
        pyplot.savefig(directory+'/'+dtfigname+'.png',format='png',bbox_inches='tight',dpi=rdpi);
        pyplot.close(dtfigname)
    s.send(np.array([kv,bfire,j]))
    print 'j',j, '\'s out'

def getI(i,dt,idt,dtRange,idtRange,iextraRange,ndt,vid,gList,pos,v0,cell,sL,vSL,n,trans,ntrans,run_t,tol_t,tstep,run_nt,tol_nt,leakyV,extralet,v1,s,savePlot,plotData,directory,rdpi):
    print " i = ", i
    jobs = []
    recv = []
    kv = np.zeros((n,ndt,run_nt))
    bfire = np.empty((n,ndt))
    for j in xrange(n):
        r, send = mp.Pipe(False)
        p = mp.Process(target = getJ, args = (i,j,dt,idt,dtRange,idtRange,iextraRange,ndt,vid,gList,pos,v0,cell,sL,vSL,n,trans,ntrans,run_t,tol_t,tstep,run_nt,tol_nt,leakyV,extralet,v1,send,savePlot,plotData,directory,rdpi))
        jobs.append(p)
        recv.append(r)
        p.start()
    print 'gathering i ', i
    result = np.array([x.recv() for x in recv])
    for p in jobs:
        p.join()
    print 'joined i ', i
    ind = [result[j][2] for j in xrange(n)]
    argi = np.argsort(ind)
    for j in xrange(n):
        jj = argi[j]
        kv[j,:,:] = result[jj][0]
        bfire[j,:] = result[jj][1]
    s.send(np.array([kv,bfire,i]))
    print 'i',i,'\'s out'
        
def getJ0(i,j,dt,idt,idtRange,ndt,vid,gList,v0,cell,sL,vSL,n,trans,ntrans,run_t,tol_t,tstep,run_nt,tol_nt,vrest,v1,s,savePlot,plotData,directory,rdpi):
    kv0 = np.zeros((run_nt))
    print " j = ", j
    if i == j:
        if idt == 0:
            s.send(np.array([kv0,j]))
            return
        else:
            spikes = np.array([[0.0, dt]])
            sel = [i]
    else:
        spikes = np.array([np.array([0.0]),np.array([dt])])
        sel = [i, j]
    RList = np.zeros((n,2))
    if idt == ndt-1:
        relit = run_nt - idtRange[idt]
        relt = run_t-dt
        v, _ = bproceed(cell, v0, sL, gList, RList, vSL, spikes, n, sel, trans, trans + run_t , 0, tstep, '')
        v = v[ntrans+idtRange[idt]:ntrans+run_nt]
    else:
        relit = tol_nt - idtRange[idt]
        relt = tol_t-dt
        v, _ = bproceed(cell, v0, sL, gList, RList, vSL, spikes, n, sel, trans, trans + tol_t, 0, tstep, '')
        v = v[ntrans+idtRange[idt]:ntrans+tol_nt]
    print relit, ' == ', v.size
    assert(v.size == relit)
    v = v - vrest
    print "max V = ", np.amax(v)
    print "max it = ", np.argmax(v)
    v10 = v1[0,i,idtRange[idt]:idtRange[idt]+relit]
    v2 = v1[0,j,:relit]
    addv =  v10 + v2 
    kvtmp = v-addv
    kv0[idtRange[idt]:idtRange[idt]+relit] = kvtmp
    print "max kv0 = ", np.amax(np.absolute(kvtmp))
    t = np.arange(run_nt)*tstep
    if plotData:
        dtfign0 = 'kv0-v'+str(vid)+'-i'+str(i)+'-j'+str(j)+'-dt'+str(idt)
        dtfig0 = pyplot.figure(dtfign0,figsize = (8,4))
        ax1 = dtfig0.add_subplot(111)         
        ax1.plot(t[idtRange[idt]:idtRange[idt]+relit],v,'c',lw=2)
        ax1.plot(t[idtRange[idt]:idtRange[idt]+relit],v10,'r',lw=2)
        ax1.plot(t[idtRange[idt]:idtRange[idt]+relit],v2,'b',lw=1)
        ax1.plot(t[idtRange[idt]:idtRange[idt]+relit],addv,'g',lw=1)
        ax1.plot(t[idtRange[idt]:idtRange[idt]+relit],kvtmp,':k')
        if savePlot:
            pyplot.savefig(directory+'/'+dtfign0+'.png',format='png',bbox_inches='tight',dpi=rdpi);
            pyplot.close(dtfign0)
    s.send(np.array([kv0,j]))
    print 'j',j, '\'s out'

def getI0(i,dt,idt,idtRange,ndt,vid,gList,v0,cell,sL,vSL,n,trans,ntrans,run_t,tol_t,tstep,run_nt,tol_nt,vrest,v1,s,savePlot,plotData,directory,rdpi):
    kv0 = np.zeros((n,run_nt))
    print " i0 = ", i
    jobs = []
    recv = []
    for j in xrange(n):
        r, send = mp.Pipe(False)
        p = mp.Process(target = getJ0, args = (i,j,dt,idt,idtRange,ndt,vid,gList,v0,cell,sL,vSL,n,trans,ntrans,run_t,tol_t,tstep,run_nt,tol_nt,vrest,v1,send,savePlot,plotData,directory,rdpi))
        jobs.append(p)
        recv.append(r)
        p.start()
    print 'gathering i0 ', i
    result = np.array([x.recv() for x in recv])
    for p in jobs:
        p.join()
    print 'joined i0 ', i
    ind = [result[ii][1] for ii in xrange(n)]
    argi = np.argsort(ind)
    for j in xrange(n):
        kv0[j,:] = result[argi[j]][0]
    s.send(np.array([kv0,i]))
    print 'i0 ',i,'s out'

def getSinglets(i,idt,dt,gList,pos,cell,v0,sL,vSL,n,trans,ntrans,relt,tstep,relit,run_nt,idtRange,leakyV,leakyDendV,spikes,s):
    print " i = ", i
    v1 = np.zeros((run_nt),dtype='double')
    dendDv = np.zeros((run_nt),dtype='double')
    sel = [i]
    RList = np.zeros((n,2))
    RList[i,:] = getGH(dt,gList[i])
    v1tmp, fired, dendv = bproceed(cell, v0, sL, gList, RList, vSL, spikes, n, sel, trans, trans + relt, 0, tstep, '', pos[i])
    dendv = dendv[ntrans:ntrans+relit] - leakyDendV[:relit]
    v1tmp = v1tmp[ntrans:ntrans+relit] - leakyV[:relit]
    tmax = idtRange[idt] + np.argmax(v1tmp)
    s.send(np.array([v1tmp,dendv,tmax,fired,i]))
    print "i", i, "\'s out"

def getExtralets(i,idt,iextraRange,gList,dt,v0,cell,sL,vSL,spikes,tstep,n,trans,ntrans,run_nt,relt,relit,leakyV,s):
    print " i = ", i
    extralet = np.zeros((run_nt),dtype='double')
    sel = [i]
    RList = np.zeros((n,2))
    RList[i,:] = getGH(dt,gList[i])
    evtmp, _ = bproceed(cell, v0, sL, gList, RList, vSL, spikes, n, sel, trans, trans + relt, 0, tstep, '')
    extralet[iextraRange[idt]:] = evtmp[ntrans:ntrans+relit] - leakyV[:relit]
    s.send(np.array([extralet,i]))
    print "i", i, "\'s out"

def getNib(argv):
    seed = 231278
    tstep = 1.0/10.0
    run_t = 225.0
    trans = 110.0
    rdpi = 300
    #dtRange = np.array([0,30,60],dtype='float')
    
    #dtRange = np.array([0,70,140],dtype='float')
    #dtRange = np.array([0,2,4,6,8,10,12,15,20,25,30,50,70,140,210],dtype='float')
    dtRangeE = np.array([6,8,10,12,15,17,19,21,23,26,30],dtype='float')
    dtRangeI = np.array([0,0.5,1,1.5,2,4,40,50,60,70,80,90,110,125,150,175],dtype='float')
    dtRange = np.sort(np.hstack([dtRangeE,dtRangeI]))

    #dtRange = np.array([0,12,24,60],dtype='float')
    #dtRange = np.array([0,8,20,130],dtype='float')
    #dtRange = np.array([0,2,4,8,12,16,20,30,40,70,100,130,210],dtype='float')
    assert(dtRange[0] == 0)
    np.random.seed(seed)
    #locE = np.array([36, 53, 74, 79, 83, 90, 97, 101, 117, 125, 132, 174],dtype='int')
   # locE = np.array([36, 74, 83, 97, 117, 132],dtype='int')
   # E/I split 1-74 | 75-134
    #locE = np.random.randint(75,134,6)
    locE = np.array([79, 82, 83, 108, 124, 129],dtype='int')
    #vRange = np.arange(-74,-59)
    #vRange = np.array([-74,-70,-65,-61],dtype='double')
    vRange = np.array([-74,-70,-67,-65,-63,-62,-61,-60,-59,-58],dtype='double')
    #locE = np.array([60, 72, 78, 84, 98],dtype='int')
    #locE = np.array([60, 72],dtype='int')
    #locE = np.array([52, 60, 72, 78, 84, 98],dtype='int')
    #locE = np.array([52, 60, 72],dtype='int')
    #locE = np.array([72],dtype='int')
    #locI = np.random.randint(1,75,3)
    #locI = np.array([41, 42, 43],dtype='int')
    #locI = np.array([40],dtype='int')
    locI = np.array([14, 28, 40],dtype='int')
    #locI = np.array([14],dtype='int')
    ### swap one
    #swap = locE[0]
    #locE[0] = locI[-1]
    #locI[-1] = swap

    #locE = np.array([32, 52, 66, 78, 98, 136],dtype='int')
    #locE = np.array([74],dtype='int')
    gE = (1e-5 + np.random.random_sample(locE.size) * (1e-1-1e-5)) * (2.0/2.0)
    #gE = np.array([2e-2])
    posE = np.random.random_sample(locE.size)
    #posE = np.array([0.6])
    #posE = np.ones(locE.size)*0.5
    #locI = np.array([7, 28, 137],dtype='int')
    #locI = np.array([2, 14, 28],dtype='int')
    #locI = np.array([28],dtype='int')
    gI = (1e-5 + np.random.random_sample(locI.size) * (1e-1-1e-5)) * (-1.0/2.0)
    posI = np.random.random_sample(locI.size)
    #posI = np.ones(locI.size)*0.5

    plotData = True 
    savePlot = True 
    vrestOnly = False
    vrest = -70
    assert(vrest in vRange)
    usage = 'getNib.py --theme=<theme name> --vid=<vid> --trans=<vClamp time> -p<plotData> -s<savePlot> -r<vrest only>// all the items in lists should be separated by comma without space'
    try:
        opts, args = getopt.getopt(argv,'psr',['theme=','vid=','trans='])
    except getopt.GetoptError:
        print "error reading args"
        print argv
        print usage
        sys.exit(2)
    print opts
    opt0 = [opt[0] for opt in opts]
    if '--vid' not in opt0:
        print 'vid must be specified'
        sys.exit(2)
    if '--theme' not in opt0:
        print 'theme name must be specified'
        sys.exit(2)
    for opt, arg in opts:
        if opt == '--vid':
            vid = int(arg)
        elif opt == '--theme':
            theme = arg
        elif opt == '--trans':
            trans = float(arg)
        elif opt == '-p':
            plotData = True
        elif opt == '-s':
            savePlot = True
            plotData = True
        elif opt == '-r':
            vrestOnly = True
        elif opt == '-h':
            print usage
            sys.exit()
    if vrestOnly:
        v0 = vrest
    else:
        v0 = vRange[vid]
    sim_t = trans + run_t
    tol_t = dtRange[-1]
    tol_nt = int(round(tol_t/tstep)) + 1
    run_nt = int(round(run_t/tstep)) + 1
    ntrans = int(round(trans/tstep))
    ndt = dtRange.size
    idtRange = np.round(dtRange/tstep).astype('int')
    
    t = np.arange(run_nt)*tstep
    print 'trans for ', ntrans, ' steps.', 'sim last', run_nt , 'steps'
    print 'v0 ', v0
    pos = np.concatenate((posE, posI))
    loc = np.concatenate((locE, locI))
    gList = np.concatenate((gE, gI))
    nE = locE.size
    nI = locI.size
    n = nE + nI
    RList = np.zeros((n,2))
    directory = '.' 
    if vid==0 and not vrestOnly:
        write_one(directory+'/'+theme+'-vRange.bin',vRange,'wb')
        write_one(directory+'/'+theme+'-dtRange.bin',dtRange,'wb')
        p = np.array([tstep,run_t,trans,n])
        write_one(directory+'/'+theme+'-p.bin',p,'wb')
        write_one(directory+'/'+theme+'-R.bin',loc,'wb')
        write_one(directory+'/'+theme+'-R.bin',pos,'ab')
        write_one(directory+'/'+theme+'-R.bin',gList,'ab')
    if not vrestOnly:
        theme = theme+'-V'+str(vid)

    cell, vecStimList, synList = prepCell(gList, loc, pos, n, vrest)
    #leakyV, _, _, _ = proceed(cell, v0, synList, RList, vecStimList, np.empty((n,0)), n, trans, run_t + trans, 0, 0, 0, 1, 0, tstep)
    leakyV, leakyDendV = leaky(cell, v0, synList, RList, vecStimList, n, trans, run_t + trans, tstep, loc, pos)
    leakyV = leakyV[ntrans:ntrans+run_nt]
    leakyDendV = leakyDendV[:,ntrans:ntrans+run_nt]
    print leakyV.size
    #if plotData:
    #    pyplot.figure('leak', figsize = (8,4))
    #    pyplot.plot(t,leakyV)
    #    pyplot.show()
    v1 = np.zeros((ndt,n,run_nt))
    sfire = np.zeros((ndt,n),dtype='int')
    tmax = np.zeros((ndt,n),dtype='int')
    dendDv = np.zeros((ndt,n,run_nt),dtype='double')
    for idt in xrange(ndt):
        #k = np.empty((n,n,run_nt))
        print " idt = ", idt
        dt = dtRange[idt]
        relit = run_nt - idtRange[idt]
        relt = run_t-dt
        spikes = np.array([[]])
        jobs = []
        recv = []
        for i in xrange(n):
            r, send = mp.Pipe(False)
            p = mp.Process(target = getSinglets, args = (i,idt,dt,gList,pos,cell,v0,synList,vecStimList,n,trans,ntrans,relt,tstep,relit,run_nt,idtRange,leakyV,leakyDendV[i,:],spikes,send))
            jobs.append(p)
            recv.append(r)
            p.start()
        print "gathering singlets idt", idt
        result = np.array([x.recv() for x in recv])
        for p in jobs:
            p.join()
        print "joined singlets idt", idt
        ind = [result[i][4] for i in xrange(n)]
        argi = np.argsort(ind)
        print "distribute singlets idt", idt
        for i in xrange(n):
            ii = argi[i]
            v1[idt,i,idtRange[idt]:] = result[ii][0]
            dendDv[idt,i,idtRange[idt]:] = result[ii][1]
            tmax[idt,i] = result[ii][2]
            sfire[idt,i] = result[ii][3]

    print "============ singlets done ==============="
    if not vrestOnly:
        if plotData:
            for i in xrange(n):
                sfigname = 'singlets'+str(i)+'-v'+str(vid)
                pyplot.figure(sfigname,figsize = (8,4))
                for idt in xrange(ndt):
                    pyplot.plot(t[idtRange[idt]:],v1[idt,i,idtRange[idt]:],'r')
                if savePlot:
                    pyplot.savefig(directory+'/'+sfigname+'.png',format='png',bbox_inches='tight',dpi=rdpi);
                    pyplot.close(sfigname)
        write_one(directory+'/'+theme+'.bin',leakyV,'wb')
        write_one(directory+'/'+theme+'.bin',leakyDendV,'ab')
        write_one(directory+'/'+theme+'.bin',v1,'ab')
        write_one(directory+'/'+theme+'.bin',tmax,'ab')
        write_one(directory+'/'+theme+'.bin',sfire,'ab')
        write_one(directory+'/'+theme+'.bin',dendDv,'ab')
        iextraRange = np.array([]).astype('int')
        for idt in xrange(ndt):
            for jdt in xrange(idt):
                ijdt = idtRange[idt] - idtRange[jdt]
                if ijdt not in idtRange and ijdt not in iextraRange:
                    iextraRange = np.append(iextraRange,ijdt)
        nedt = iextraRange.size 
        extraRange = iextraRange*tstep
        print "extra dtRange:"
        print extraRange
        extralet = np.zeros((nedt,n,run_nt))
        for idt in xrange(nedt):
            #k = np.empty((n,n,run_nt))
            print " idt = ", idt
            dt = extraRange[idt]
            relit = run_nt - iextraRange[idt]
            relt = run_t-dt
            spikes = np.array([[]])
            jobs = []
            recv = []
            for i in xrange(n):
                print " i = ", i
                r, send = mp.Pipe(False)
                p = mp.Process(target = getExtralets, args = (i,idt,iextraRange,gList,dt,v0,cell,synList,vecStimList,spikes,tstep,n,trans,ntrans,run_nt,relt,relit,leakyV,send))
                jobs.append(p)
                recv.append(r)
                p.start()
            print "gathering extralets idt", idt
            result = np.array([x.recv() for x in recv])
            for p in jobs:
                p.join()
            print "joined extralets idt", idt
            ind = [result[i][1] for i in xrange(n)]
            argi = np.argsort(ind)
            print "distribute extralets idt", idt
            for i in xrange(n):
                extralet[idt,i,:] = result[argi[i]][0]

        print "============ extralets done =============="
        bfire = np.zeros((ndt,n,n,ndt),dtype=int)
        for idt in xrange(ndt):
            print " idt = ", idt
            dt = dtRange[idt]
            kv = np.zeros((n,n,ndt,run_nt))
            jobs = []
            recv = []
            for i in xrange(n):
                r, send = mp.Pipe(False)
                p = mp.Process(target = getI, args = (i,dt,idt,dtRange,idtRange,iextraRange,ndt,vid,gList,pos,v0,cell,synList,vecStimList,n,trans,ntrans,run_t,tol_t,tstep,run_nt,tol_nt,leakyV,extralet,v1,send,savePlot,plotData,directory,rdpi))
                jobs.append(p)
                recv.append(r)
                p.start()
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
    # kv0
    else:
        kv0 = np.zeros((ndt,n,n,run_nt))
        for idt in xrange(ndt):
            print " idt = ", idt
            dt = dtRange[idt]
            jobs = []
            recv = []
            for i in xrange(n):
                r, send = mp.Pipe(False)
                p = mp.Process(target = getI0, args = (i,dt,idt,idtRange,ndt,vid,gList,v0,cell,synList,vecStimList,n,trans,ntrans,run_t,tol_t,tstep,run_nt,tol_nt,vrest,v1,send,savePlot,plotData,directory,rdpi))
                jobs.append(p)
                recv.append(r)
                p.start()
            print "gathering kv0 idt", idt
            result = np.array([x.recv() for x in recv])
            for p in jobs:
                p.join()
            print "joined kv0 idt", idt
            ind = [result[i][1] for i in xrange(n)]
            argi = np.argsort(ind)
            print "distribute kv0 idt", idt
            for i in xrange(n):
                kv0[idt,i,:,:] = result[argi[i]][0]
        write_one(directory+'/'+theme+'-kv0.bin',kv0,'ab')
if __name__ == '__main__':
    getNib(sys.argv[1:])
