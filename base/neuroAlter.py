from neuron import h
h.load_file('stdlib.hoc')
h.load_file('stdrun.hoc')
import numpy as np
from n128 import *
import sys
import matplotlib
matplotlib.use('Agg')
#matplotlib.use('Qt5Agg')
from matplotlib import pyplot, ticker
f = open('pylog','w')
f.close()

class receptor(object):
    def __init__(self, a, b, c, d, cm):
        f = open('pylog','a')
        self.a = a
        self.b = b
        self.c = c
        self.d = d
        self.cm = cm
        self.Rtau = 1 / (cm*a + b)
        self.Rinf = self.Rtau * cm*a 
        self.expr = np.exp(-c/self.Rtau)
        print 'decayTau = ', 1/self.b
        print >>f, 'decayTau = ', 1/self.b
        print 'Rtau = ', self.Rtau
        print >>f, 'Rtau = ', self.Rtau
        print 'Rinf = ', self.Rinf
        print >>f, 'Rinf = ', self.Rinf
        f.close()

def getR(RList, gList, spikeTrain, t0, rel, front, sel, ampa, gaba):
    for i in sel:
        R = 0
        if gList[i] < 0:
            recep = gaba
        else:
            recep = ampa
        lastRel = -1e23
        for j in xrange(rel[i],front[i]) :
            nextRel = spikeTrain[i][j]
            if nextRel - lastRel > recep.d + recep.c:
                R = recep.Rinf + (R - recep.Rinf) * recep.expr
                R = R * np.exp (- recep.b * (nextRel - (lastRel + recep.c)))
                lastRel = nextRel

        if lastRel > t0 - recep.c:
            RList[i] = recep.Rinf + (R - recep.Rinf) * exp (- (t0 - lastRel) / recep.Rtau)
        else:
            RList[i] = recep.Rinf + (R - recep.Rinf) * recep.expr
            RList[i] = R * np.exp (- recep.b * (t0 - (lastRel + recep.c)))

def run(cell, v0, vBack, tref, vThres, synList, RList, n, trans, oneGo, t0, alphaR = True, loc=np.array([]), pos=np.array([]), pas = False, v=h.Vector, t=h.Vector(), rdv=[h.Vector()], cpi = False, cpt = -1.0):
    f = open('pylog','a')
    nc = 0
    #if not oneGo:
    #    cell.init(v0)
    #else:
    if cpi:
        cell.cp_init()
        if trans > 0:
            print "starting with copied voltage (and gating states), trans should be 0 if not intentional"
        print '    cell cp_initiated'
        print >>f, '    cell cp_initiated'
    else:
        cell.init()
        print '    cell initiated'
        print >>f, '    cell initiated'
    steps = 0
    tsps = np.empty(int((h.tstop-t0)/h.dt/2))
    tsp = -tref
    #plusOne = 0
    vold = v0
    firing = 0 
    h.t = t0
    print  '    t0 = ', h.t
    print >>f, '    t0 = ', h.t
    neTrans = int(round((t0+trans)/h.dt))
    while int(round(h.t/h.dt)) < int(round((t0+trans)/h.dt)):
        h.fadvance()
        steps = steps + 1
        #print  cell.soma(0.5).v
        #print >>f,  steps, h.t
    for i in xrange(n):
        #print  str(i), 'th's Synapse = ', RList[i]
        #print >>f,  str(i), 'th's Synapse = ', RList[i]
        if alphaR:
            synList[i].g = RList[i][0]
            synList[i].h = RList[i][1]
        else:
            synList[i].R = RList[i]
    #oldProgress = 0
    print  '   ', trans, 'ms trans complete, used', steps , 'steps, t+trans =', h.t, 'v =', '%7.5f.' % cell.soma(0.5).v
    print >>f,  '   ', trans, 'ms trans complete, used', steps , 'steps, t+trans =', h.t, 'v =', '%7.5f.' % cell.soma(0.5).v
    if len(loc)>0:
        for i in xrange(n): 
            print >> f, '    dend[', loc[i], ']' ,'(',pos[i],').v =', cell.dend[loc[i]](pos[i]).v, ',',
            print '    dend[', loc[i], ']' ,'(',pos[i],').v =', cell.dend[loc[i]](pos[i]).v, ',',
            if alphaR:
                print >> f, synList[i].g, ',', synList[i].h
                print synList[i].g, ',', synList[i].h
            else:
                print >>f, synList[i].R
                print synList[i].R
    if cpt > -1 and int(round(cpt/h.dt)) == int(round(h.t/h.dt)):
        cell.copy_volt()
        if cell.include_states:
            cell.copy_state()
    if not pas:
        if oneGo:
            while (int(round(h.t/h.dt)) < int(round(h.tstop/h.dt))):
                h.fadvance()
                if cpt > -1 and int(round(cpt/h.dt)) == int(round(h.t/h.dt)):
                    cell.copy_volt()
                    if cell.include_states:
                        cell.copy_state()
                if cell.soma(0.5).v > vThres+15 and cell.soma(0.5).v < vold and not firing:
                    tsp = h.t-trans
                    print 'nc =', nc, '<', tsps.size 
                    tsps[nc] = tsp
                    firing = 1
                    nc = nc + 1
                    print  tsp, 'fired expecting finish', tsp + tref
                    print >>f,  tsp, 'fired expecting finish', tsp + tref
                vold = cell.soma(0.5).v
                if (cell.soma(0.5).v <= vThres+7.5):
                    firing = 0
        else:
            while ((cell.soma(0.5).v > vBack or firing) and int(round(h.t/h.dt)) < int(round(h.tstop/h.dt))):
                h.fadvance()
                if cell.soma(0.5).v > vThres+15 and cell.soma(0.5).v < vold and not firing:
                    tsp = h.t-trans
                    tsps[nc] = tsp
                    firing = 1
                    nc = nc + 1
                    print >>f,  tsp, 'fired'
                    print  tsp, 'fired'
                if cell.soma(0.5).v <= vBack and cell.soma(0.5).v > vold and firing:
                    firing = 0
                vold = cell.soma(0.5).v
            print >> f, '    ((', cell.soma(0.5).v > vBack, ' or ', firing > 0, ') and ', int(round(h.t/h.dt)) < int(round(h.tstop/h.dt)), ')'
            print '    ((', cell.soma(0.5).v > vBack, ' or ', firing > 0, ') and ', int(round(h.t/h.dt)) < int(round(h.tstop/h.dt)), ')'
    else:
        assert(oneGo==True)
        while (int(round(h.t/h.dt)) < int(round(h.tstop/h.dt))):
            if cpt > -1 and int(round(cpt/h.dt)) == int(round(h.t/h.dt)):
                cell.copy_volt()
                if cell.include_states:
                    cell.copy_state()
            if cell.soma(0.5).v >= vThres:
                tsp = h.t-trans
                print >>f, '    nc =', nc, '<', tsps.size 
                print '    nc =', nc, '<', tsps.size 
                tsps[nc] = tsp
                nc = nc + 1
                print >>f,  tsp, 'fired expecting finish', tsp + tref
                print  tsp, 'fired expecting finish', tsp + tref

            if h.t < trans + tsp + tref:
                for sec in cell.all: 
                    sec.v = cell.Vrest
                h.t = h.t + h.dt
                v.append(cell.Vrest)
                t.append(h.t)
                for i in xrange(n):
                    rdv[i].append(cell.Vrest)
            else:
                h.fadvance()

    print  'stopping with v', cell.soma(0.5).v, '<', vBack, 'firing =', firing, 'or t', h.t - trans, '==', h.tstop - trans
    print >>f, 'stopping with v', cell.soma(0.5).v, '<', vBack, 'firing =', firing, 'or t', h.t - trans, '==', h.tstop - trans
    f.close()
    return nc, tsps[:nc]

def prepCell(gList, loc, pos, n, v0, alphaR = True):
    f = open('pylog','a')
    #print '    location,  strength,   pos'
    #print >>f, '    location,  strength,   pos'
    #for i in xrange(n):
    #    print '#{:2d}   {:3d}      {: 1.1e}    {:0.1f}'.format(i,loc[i],gList[i],pos[i])
    #    print >>f, '#{:2d}   {:3d}      {: 1.1e}    {:0.1f}'.format(i,loc[i],gList[i],pos[i])
    print >> f, 'gList = np.array([',
    print 'gList = np.array([',
    for i in xrange(n):
        if i == n-1:
            print >> f, gList[i],
            print gList[i],
        else:
            print >> f, gList[i],', ',
            print gList[i],', ',
    print >> f, '])'
    print '])'

    print >> f, 'loc = np.array([',
    print 'loc = np.array([',
    for i in xrange(n):
        if i == n-1:
            print >> f, loc[i],
            print loc[i],
        else:
            print >> f, loc[i],', ',
            print loc[i],', ',
    print >> f, '])'
    print '])'

    print >> f, 'pos = np.array([',
    print 'pos = np.array([',
    for i in xrange(n):
        if i == n-1:
            print >> f, pos[i],
            print pos[i],
        else:
            print >> f, pos[i],', ',
            print pos[i],', ',
    print >> f, '])'
    print '])'

    h.cvode.active(0)
    cell = RealisticNeuron(v0)
    cell.soma.push()
    h.distance()
    h.pop_section()
    vecStimList = [h.VecStim() for i in xrange(n)]
    synList = []
    dist = np.zeros((n));
    print >> f, 'distance = np.array([',
    print 'distance = np.array([',
    for i in xrange(n):
        cell.dend[loc[i]].push()
        if gList[i] < 0:
            if alphaR:
                synList.append(h.ALPHAI())
            else:
                synList.append(h.GABAa())
        else:
            if alphaR:
                synList.append(h.ALPHAE())
            else:
                synList.append(h.AMPA())
        synList[i].loc(pos[i])
        h.setpointer(vecStimList[i]._ref_y,'pre',synList[i])
        h.pop_section()
        if gList[i] < 0:
            if alphaR:
                synList[i].f = -gList[i]
            else:
                synList[i].gmax = -gList[i]
        else:
            if alphaR:
                synList[i].f = gList[i]
            else:
                synList[i].gmax = gList[i]
        dist[i] = h.distance(cell.dend[loc[i]](pos[i]));
        if i==n-1:
            print >> f, dist[i],
            print dist[i],
        else:
            print >> f, dist[i], ', ',
            print dist[i], ', ',

    print >> f, '])'
    print '])'
    f.close()
    return cell, vecStimList, synList, dist

def set_volt(cell,v):
    cell.set_volt(v)

def proceed(cell, v0, synList, RList, vecStimList, spikeTrain, n, trans, tend, vBack, tref, vThres, oneGo, t0, tstep, loc=np.array([]), pos=np.array([]), dendVclamp=np.array([]), alphaR = True, getDendV = False, monitorDend = False, pas = False, plotSlice=False, fign='slice', cpi=False, cpt=-1.0):
    f = open('pylog','a')
    h.tstop = tend
    print '    ', t0, ' with ', trans, ' trans ', ' to ', tend , ', tref: ', tref
    print >>f, '    ', t0, ' with ', trans, ' trans ', ' to ', tend , ', tref: ', tref
    print '    ', 'vBack: ', vBack, ' vThres: ', vThres, ' vStart ', v0
    print >>f,   '    ', 'vBack: ', vBack, ' vThres: ', vThres, ' vStart: ', v0
    print '    RList:', RList
    print >>f,   '    RList:', RList
    print '    spikeTrain:', spikeTrain
    print >>f,   '    spikeTrain:', spikeTrain
    print '    oneGo: ', oneGo
    print >>f, '    oneGo: ', oneGo
    print '    loc: ', loc
    print >>f, '    loc: ', loc
    print '    pos: ', pos
    print >>f, '    pos: ', pos
    h.dt = tstep
    v = h.Vector()
    t = h.Vector()
    v.record(cell.soma(0.5)._ref_v)
    t.record(h._ref_t)
    if trans > 0:
        cell.soma.push()
        vHold = h.SEClamp(0.5)
        vHold.dur2 = 0.0
        vHold.dur3 = 0.0
        vHold.rs = 1e-9
        vHold.dur1 = t0 + trans
        vHold.amp1 = v0
        h.pop_section()
        print '    soma clamp'
        vHolds = []
        if not monitorDend:
            j = 0
            for i in xrange(len(loc)):
                if not abs(dendVclamp[i] -1000) < 1e-10:
                    cell.dend[loc[i]].push()
                    vHolds.append(h.SEClamp(pos[i]))
                    vHolds[j].dur2 = 0.0
                    vHolds[j].dur3 = 0.0
                    if oneGo:
                        vHolds[j].rs = cell.dend[loc[i]].Ra
                    elif dendVclamp[i] > 0:
                        vHolds[j].rs = 1e-9
                        dendVclamp[i] = -dendVclamp[i]
                    vHolds[j].dur1 = t0 + trans
                    vHolds[j].amp1 = dendVclamp[i]
                    j = j + 1
                    h.pop_section()
                    print '    dend[',loc[i],'] clamped'
        else:
            print cell.ndend, '    dendrites clampled to ', v0
            for i in xrange(cell.ndend):
                cell.dend[i].push()
                vHolds.append(h.SEClamp(0.5))
                vHolds[i].dur2 = 0
                vHolds[i].dur3 = 0
                vHolds[i].rs = 1e-9
                vHolds[i].dur1 = t0 + trans
                vHolds[i].amp1 = cell.Vrest + (v0 - cell.Vrest) * 1.0
                h.pop_section()
    dendv = []
    if monitorDend:
        for i in xrange(cell.ndend):
            dendv.append(h.Vector())
            dendv[3*i].record(cell.dend[i](0.16666)._ref_v)
            dendv.append(h.Vector())
            dendv[3*i+1].record(cell.dend[i](0.5)._ref_v)
            dendv.append(h.Vector())
            dendv[3*i+2].record(cell.dend[i](0.83333)._ref_v)
    else:
        for i in xrange(n):
            dendv.append(h.Vector())
            dendv[i].record(cell.dend[loc[i]](pos[i])._ref_v)
    for i in xrange(n):
        if len(spikeTrain) > 0:
            vecStimList[i].play(h.Vector(spikeTrain[i] + trans))
            vecStimList[i].dt = h.dt
            if alphaR:
                print '    played ', i, ' s: ', synList[i].g
                print >>f, '    played ', i, ' s: ', synList[i].g
            else:
                print '    played ', i, ' s: ', synList[i].gmax
                print >>f, '    played ', i, ' s: ', synList[i].gmax
            if alphaR:
                synList[i].deltat = h.dt
    print '    ready to run'
    print >>f, '    ready to run'
    fired, tsp = run(cell, v0, vBack, tref, vThres, synList, RList, n, trans, oneGo, t0, alphaR, loc, pos, pas, v, t, dendv, cpi, cpt)
    v1 = v.as_numpy()
    if oneGo or int(round(h.t/h.dt)) > int(round(h.tstop/h.dt)):
        ntotal = int(round((tend-t0)/tstep)+1)
    else:
        ntotal = int(round((h.t-t0)/tstep)+1)
    #if __name__ == '__main__':
    t1 = t.as_numpy()
    if v1.size!=ntotal:
        print >>f, '    v1 steps ', v1.size, ', ntotal ', ntotal
        print '    v1 steps ', v1.size, ', ntotal ', ntotal
        print >>f, '    remove extra step'
        print '    remove extra step'
        v1 = v1[:-1]
        #if __name__ == '__main__':
        t1 = t1[:-1]
    assert(t1.size == v1.size)
    print >>f, '    is this nan ', v1[-1], v1[0]
    print   '    is this nan ', v1[-1], v1[0]
    print >>f, '    v size ', v1.size
    print   '    v size ', v1.size
    ntrans = int(round(trans/tstep))
    print >>f, '    trans size ', ntrans
    print '    trans size ', ntrans
    print >>f, '    trans end, v0, v1', v1[ntrans-1], v1[ntrans], v1[-1]
    print '    trans end, v0, v1', v1[ntrans-1], v1[ntrans], v1[-1]
    print >>f, '    v w/o trans size ', v1.size - ntrans
    print '    v w/o trans size ', v1.size - ntrans
    if trans > 0:
        vHold = None

    if plotSlice: 
        fig = pyplot.figure('slice',figsize=(8,4))
        ax = fig.add_subplot(1,1,1)
        ax.minorticks_on()
        labels = ['dend['+str(loc[x])+']' for x in range(n)]
        #istart = ntrans
        istart = 1 
        if not monitorDend:
            #pyplot.plot(t1[ntrans:]-trans,v1[ntrans:],'k')
            ymax=-60
            if not oneGo:
                vline, = ax.plot(t1[istart:]-trans,v1[istart:],'k')
            else:
                vline, = ax.plot(t1[istart:]-trans,v1[istart:],'k')
            if np.amax(v1) > ymax:
                ymax = np.amax(v1)
            dline = []
            for i in xrange(n):
                dv = dendv[i].as_numpy()
                dv = dv[:v1.size]
                #pyplot.plot(t1[ntrans:]-trans,dv[ntrans:],':')
                if not oneGo:
                    dl, = ax.plot(t1[istart:]-trans,dv[istart:],':')
                else:
                    dl, = ax.plot(t1[istart:]-trans,dv[istart:],':')
                dline.append(dl)
                if np.amax(dv) > ymax:
                    ymax = np.amax(dv)
            ax.set_ylim(-80,np.amin([-40,ymax]))
            #pyplot.ylim(-69,-62)
            pyplot.legend(dline,labels)
            ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
            ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
        else:
            nt = ntrans
            dendVec = np.empty((3*cell.ndend,ntotal),dtype='float')
            for i in xrange(cell.ndend):
                dendVec[3*i,:] = dendv[3*i].as_numpy()
                dendVec[3*i+1,:] = dendv[3*i+1].as_numpy()
                dendVec[3*i+2,:] = dendv[3*i+2].as_numpy()
            v1 = v.as_numpy()
            t1 = t.as_numpy()
            dendVec = np.reshape(dendVec,(3,199,ntotal))
            mdend = np.average(dendVec, axis = 0)
            print '    mdend ',mdend.shape
            #qdendlow = np.percentile(mdend,10, axis = 0)
            qdendlow = np.amin(mdend, axis = 0)
            print '    qdendlow',qdendlow.shape
            #qdendhigh = np.percentile(mdend,90, axis = 0)
            qdendhigh = np.amax(mdend, axis = 0)
            print '    qdendhigh',qdendhigh.shape
            dendAverage = np.average(mdend, axis=0)
            print '    dendAverage',dendAverage.shape
            denderrbar = np.vstack((qdendlow,qdendhigh))
            print '    denderrbar',denderrbar.shape
            istart = 0
            iend = t1.size
            #istart = nt 
            t2 = t1[istart:iend]
            ax.plot(t2-istart*tstep,v1[istart:iend])
            dA = dendAverage[istart:iend]
            dE = denderrbar[:,istart:iend]
            select = np.arange(0,t2.size,10)
            t2 = t2[select] - istart*tstep
            print '    t2',t2.shape
            dA = dA[select]
            print '    dA',dA.shape
            dE = dE[:,select]
            print '    dE',dE.shape
            ax.plot(t2,dA,'--k')
            ax.plot(t2,dE[0,:],':k')
            ax.plot(t2,dE[1,:],':k')
            ax.errorbar(t2,dA,np.absolute(np.vstack((dA,dA))-dE))
            qdendlow = np.percentile(mdend,5, axis = 0)
            qdendhigh = np.percentile(mdend,95, axis = 0)
            denderrbar = np.vstack((qdendlow,qdendhigh))
            dE = denderrbar[:,istart:]
            dE = dE[:,select]
            ax.errorbar(t2+0.2,dA,np.absolute(np.vstack((dA,dA))-dE),ls = 'None')
            qdendlow = np.percentile(mdend,25, axis = 0)
            qdendhigh = np.percentile(mdend,75, axis = 0)
            denderrbar = np.vstack((qdendlow,qdendhigh))
            dE = denderrbar[:,istart:]
            dE = dE[:,select]
            ax.errorbar(t2+0.4,dA,np.absolute(np.vstack((dA,dA))-dE),ls = 'None')
            qdendlow = np.percentile(mdend,45, axis = 0)
            qdendhigh = np.percentile(mdend,55, axis = 0)
            denderrbar = np.vstack((qdendlow,qdendhigh))
            dE = denderrbar[:,istart:]
            dE = dE[:,select]
            ax.errorbar(t2+0.6,dA,np.absolute(np.vstack((dA,dA))-dE),ls = 'None')
            ax.set_ylim(v0-0.25,v0+0.25)
        pyplot.savefig(fign+'.png',format='png',bbox_inches='tight',dpi=900)
        pyplot.close()
        print fign,'.png saved'
        
    f.close()
    if getDendV:
        dendVecOut = np.array([])
        for i in xrange(n):
            if oneGo and monitorDend:
                dvtmp = dendv[loc[i]].as_numpy().copy()
            else:
                dvtmp = dendv[i].as_numpy().copy()
            if dvtmp.size != ntotal:
                print '    dendv steps ', dvtmp.size, ', ntotal ', ntotal
                print '    remove extra step'
                dvtmp = dvtmp[:-1]
            dendVecOut = np.concatenate((dendVecOut,dvtmp))
        print 'returning dendV', dendVecOut.size
        return v1.copy(), fired, tsp, ntrans, dendVecOut
    else:
        print 'returning somaV only'
        return v1.copy(), fired, tsp, ntrans

def leaky(cell, v0, synList, RList, vecStimList, n, trans, tend, tstep, loc, pos, alphaR = True, cpi = False, cpt = -1.0):
    f = open('pylog','a')
    h.tstop = tend
    print loc
    print >>f, loc
    print pos
    print >>f, pos

    h.dt = tstep
    v = h.Vector()
    t = h.Vector()
    v.record(cell.soma(0.5)._ref_v)
    t.record(h._ref_t)
    if trans > 0:
        cell.soma.push()
        vHold = h.SEClamp(0.5)
        vHold.dur2 = 0
        vHold.dur3 = 0
        vHold.rs = 1e-9
        vHold.dur1 = trans
        vHold.amp1 = v0
        h.pop_section()
        print ' soma clamp'
    dendv = []
    for i in xrange(n):
        dendv.append(h.Vector())
        dendv[i].record(cell.dend[loc[i]](pos[i])._ref_v)
        print 'record ', i
        print >>f, 'record ', i
        if alphaR:
            synList[i].deltat = h.dt
        vecStimList[i].play(h.Vector([]))
        vecStimList[i].dt = h.dt
    print   'ready to run'
    print >>f, 'ready to run'
    run(cell, v0, 0, 0, 0, synList, RList, n, trans, 1, 0, alphaR, loc, pos, False, v, t, dendv, cpi, cpt)
    #print   'i ran and i killed ', fired, ' bedbug(s)'
    #print >>f,   'i ran and i killed ', fired, ' bedbug(s)'
    v1 = v.as_numpy()
    ntotal = int(round(tend/tstep)+1)
    if __name__ == '__main__':
        t1 = t.as_numpy()
    if v1.size!=ntotal:
        print ' steps ', v1.size, ', ntotal ', ntotal
        print ' remove extra step'
        print >>f, ' remove extra step'
        v1 = v1[:-1]
        if __name__ == '__main__':
            t1 = t1[:-1]
    print   ' is this nan ', v1[-1], v1[-2], v1[-3]
    print >>f,   ' is this nan ', v1[-1], v1[-2], v1[-3]
    print   ' v size ', v1.size
    print >>f,   ' v size ', v1.size
    ntrans = int(round(trans/tstep))
    #if v1[ntrans] == v1[ntrans-1]:
    #    ntrans = ntrans+1
    #assert(v1[ntrans] != v1[ntrans-1])
    print ' trans size ', ntrans
    print   ' trans end, v0, v1', v1[ntrans-1], v1[ntrans], v1[ntrans+1]
    print >>f,   ' trans end, v0, v1', v1[ntrans-1], v1[ntrans], v1[ntrans+1]
    print   ' v w/o trans size ', v1.size - ntrans
    print >>f,   ' v w/o trans size ', v1.size - ntrans
    if trans > 0:
        vHold = None
    if __name__ == '__main__':
        nt = int((trans)/tstep)
        pyplot.figure('slice',figsize=(8,4))
        pyplot.plot(t1[nt:],v1[nt:])
        #pyplot.show()
        pyplot.savefig('slice.png',format='png',bbox_inches='tight',dpi=900)
    #print 'highest voltage ', np.amax(v1), ' in ', v1.size, ' steps'
    #print >>f, 'highest voltage ', np.amax(v1), ' in ', v1.size, ' steps'
    dendvArray = np.empty((n,ntotal))
    for i in xrange(n):
        dvtmp = dendv[i].as_numpy().copy()
        if dvtmp.size != ntotal:
            print '    dendv steps ', dvtmp.size, ', ntotal ', ntotal
            print '    remove extra step'
            dvtmp = dvtmp[:-1]
        dendvArray[i,:] = dvtmp
    f.close()
    return v1.copy(), dendvArray

def sproceed(cell, v0, synList, gList, RList, vecStimList, spike, n, sel, trans, tend, tstep, name, alphaR = True):
    t0 = 0.0
    f = open('pylog','a')
    h.tstop = tend
    print ' single synapse event '
    print >>f, ' single synapse event '
    print t0, ' with ', trans, ' trans ', ' to ', tend 
    print >>f, t0, ' with ', trans, ' trans ', ' to ', tend 
    h.dt = tstep
    v = h.Vector()
    t = h.Vector()
    v.record(cell.soma(0.5)._ref_v)
    t.record(h._ref_t)
    if trans > 0:
        vHold = h.SEClamp(0.5)
        vHold.dur2 = 0
        vHold.dur3 = 0
        vHold.rs = 1e-9
        vHold.dur1 = t0 + trans
        vHold.amp1 = v0
    for i in xrange(n):
        if i == sel:
            print ' assigning single spike event', i, spike
            print >>f, ' assigning single spike event', i, spike
            if alphaR:
                synList[i].f = abs(gList[i])
                synList[i].deltat = h.dt
            else:
                synList[i].gmax = abs(gList[i])
            vecStimList[i].play(h.Vector([spike+trans]))
            vecStimList[i].dt = h.dt
        else:
            if alphaR:
                synList[i].deltat = h.dt
                synList[i].f = 0
            else:
                synList[i].gmax = 0
            vecStimList[i].play(h.Vector([]))
            vecStimList[i].dt = h.dt

    fired, _ = run(cell, v0, 0, 0, 0, synList, RList, n, trans, 1, t0, alphaR)
    if fired > 0:
        print >>f, 'fired, but not supposed to'
        print 'fired, but not supposed to'
    #print 'i ran and i killed ', fired, ' bedbug(s)'
    #print >>f, 'i ran and i killed ', fired, ' bedbug(s)'
    #assert(fired == 0)
    v1 = v.as_numpy()
    if trans > 0:
        vHold = None
    if name:
        t1 = t.as_numpy()
        nt = int((trans)/tstep)
        pyplot.figure(name,figsize=(8,4))
        pyplot.plot(t1[nt:],v1[nt:])
        pyplot.savefig(name+'.png',format='png',bbox_inches='tight',dpi=900)
    #print 'highest voltage ', np.amax(v1), ' in ', v1.size, ' steps'
    #print >>f, 'highest voltage ', np.amax(v1), ' in ', v1.size, ' steps'
    f.close()
    return v1.copy(), fired

def bproceed0(cell, v0, synList, gList, RList, vecStimList, spikeTrain, n, sel, trans, dt, tend, t0, tstep, name, pos=-1.0, printR = False, alphaR = True):
    f = open('pylog','a')
    h.tstop = tend
    print t0, ' with ', trans, ' trans ', ' to ', tend, ' v0 = ', v0

    print >>f, t0, ' with ', trans, ' trans ', ' to ', tend, ' v0 = ', v0
    h.dt = tstep
    v = h.Vector()
    if pos != -1.0:
        dendv = h.Vector()
        dendv.record(cell.dend[sel[0]](pos)._ref_v)
    t = h.Vector()
    v.record(cell.soma(0.5)._ref_v)
    t.record(h._ref_t)
    if printR:
        if alphaR:
            g = h.Vector()
            hh = h.Vector()
            g.record(synList[sel[0]]._ref_g)
            hh.record(synList[sel[0]]._ref_h)
        else:
            R = h.Vector()
            R.record(synList[sel[0]]._ref_R)
    if trans + dt> 0:
        vHold = h.SEClamp(0.5)
        vHold.dur2 = 0
        vHold.dur3 = 0
        vHold.rs = 1e-11
        vHold.dur1 = t0 + trans + dt
        vHold.amp1 = v0
    for i in xrange(len(sel)):
        if spikeTrain[i].size > 0:
            print ' assigning 2 spike events', sel[i], spikeTrain[i]+trans
            print >>f, ' assigning 2 spike events', sel[i], spikeTrain[i]+trans
            vecStimList[sel[i]].play(h.Vector(spikeTrain[i]+trans))
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
            #vecStimList[i] = h.VecStim()
            #h.setpointer(vecStimList[i]._ref_y,'pre',synList[i])
            vecStimList[i].dt = h.dt
            if alphaR:
                synList[i].deltat = h.dt
                synList[i].f = 0
            else:
                synList[i].gmax = 0

    fired, _ = run(cell, v0, 0, 0, 0, synList, RList, n, trans + dt, 1, t0, alphaR)
    #print 'i ran and i killed ', fired, ' bedbug(s)'
    #print >>f, 'i ran and i killed ', fired, ' bedbug(s)'
    if fired > 0:
        print >>f, 'fired, but not supposed to'
        print 'fired, but not supposed to'
    #assert(fired == 0)
    v1 = v.as_numpy()
    if printR:
        if alphaR:
            g1 = g.as_numpy()
            h1 = hh.as_numpy()
        else:
            r1 = R.as_numpy()

    if trans + dt > 0:
        vHold = None
    if name:
        t1 = t.as_numpy()
        nt = int((trans)/tstep)
        pyplot.figure(name,figsize=(8,4))
        pyplot.plot(t1[nt:],v1[nt:])
        pyplot.savefig(name+'.png',format='png',bbox_inches='tight',dpi=900)
        pyplot.close()
    #print 'highest voltage ', np.amax(v1), ' in ', v1.size, ' steps'
    #print >>f, 'highest voltage ', np.amax(v1), ' in ', v1.size, ' steps'
    f.close()
    if printR:
        if alphaR:
            return v1.copy(), fired, g1.copy(), h1.copy() 
        else:
            return v1.copy(), fired, r1.copy()
    else:
        if pos == -1.0:
            return v1.copy(), fired
        else:
            return v1.copy(), fired, dendv.as_numpy().copy()

def bproceed(cell, v0, synList, gList, RList, vecStimList, spikeTrain, n, sel, trans, tend, t0, tstep, name, pos = -1.0, loc = 0, alphaR = True, cpi = False, cpt = -1.0):
    f = open('pylog','a')
    h.tstop = tend
    print t0, ' with ', trans, ' trans ', ' to ', tend, ' v0 = ', v0

    print >>f, t0, ' with ', trans, ' trans ', ' to ', tend, ' v0 = ', v0
    print 'gList:'
    print gList
    h.dt = tstep
    v = h.Vector()
    dendv = h.Vector()
    if pos != -1.0: # must be singlet
        assert(len(sel) == 1)
        dendv.record(cell.dend[loc](pos)._ref_v)
    t = h.Vector()
    v.record(cell.soma(0.5)._ref_v)
    t.record(h._ref_t)
    if trans> 0:
        vHold = h.SEClamp(0.5)
        vHold.dur2 = 0
        vHold.dur3 = 0
        vHold.rs = 1e-11
        vHold.dur1 = t0 + trans
        vHold.amp1 = v0
    for i in xrange(len(sel)):
        if spikeTrain[i].size > 0:
            print ' assigning 2 spike events', sel[i], spikeTrain[i]+trans
            print >>f, ' assigning 2 spike events', sel[i], spikeTrain[i]+trans
            vecStimList[sel[i]].play(h.Vector(spikeTrain[i]+trans))
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
            #vecStimList[i] = h.VecStim()
            #h.setpointer(vecStimList[i]._ref_y,'pre',synList[i])
            vecStimList[i].dt = h.dt
            if alphaR:
                synList[i].deltat = h.dt
                synList[i].f = 0
            else:
                synList[i].gmax = 0

    fired, _ = run(cell, v0, 0, 0, 0, synList, RList, n, trans, 1, t0, alphaR, np.array([]), np.array([]), False, v, t, dendv, cpi, cpt)
    #print 'i ran and i killed ', fired, ' bedbug(s)'
    #print >>f, 'i ran and i killed ', fired, ' bedbug(s)'
    if fired > 0:
        print >>f, 'fired, but not supposed to'
        print 'fired, but not supposed to'
    #assert(fired == 0)
    v1 = v.as_numpy()
    if trans > 0:
        vHold = None
    if name:
        t1 = t.as_numpy()
        nt = int((trans)/tstep)
        pyplot.figure(name,figsize=(8,4))
        pyplot.plot(t1[nt:],v1[nt:])
        pyplot.savefig(name+'.png',format='png',bbox_inches='tight',dpi=900)
        pyplot.close()
    #print 'highest voltage ', np.amax(v1), ' in ', v1.size, ' steps'
    #print >>f, 'highest voltage ', np.amax(v1), ' in ', v1.size, ' steps'
    f.close()
    if pos == -1.0:
        return v1.copy(), fired
    else:
        return v1.copy(), fired, dendv.as_numpy().copy()

def srun(cell, v0, trans, t0):
    f = open('pylog','a')
    cell.init()
    print ' cell initiated'
    print >>f, ' cell initiated'
    steps = 0
    h.t = t0
    print  't0 = ', h.t
    print >>f,  't0 = ', h.t
    while int(round(h.t/h.dt)) < int(round((t0+trans)/h.dt)):
        h.fadvance()
        steps = steps + 1
    print trans, 'ms trans complete, used ', steps , ' steps, t+trans = ', h.t, 'v = ', '%7.5f.' % cell.soma(0.5).v
    print' distant apical dend v = ', '%7.5f.' % cell.dend[121](0.0).v
    print' distant basal dend v = ', '%7.5f.' % cell.dend[153](0.0).v
    print >>f,  trans, 'ms trans complete, used ', steps , ' steps, t+trans = ', h.t, 'v = ', '%7.5f.' % cell.soma(0.5).v
    print >>f, ' distant apical dend v = ', '%7.5f.' % cell.dend[121](0.0).v
    print >>f, ' distant basal dend v = ', '%7.5f.' % cell.dend[153](0.0).v
    while (h.t < h.tstop):
        h.fadvance()
        steps = steps + 1
    print 'srun ended in ', steps , ' steps, t+trans = ', h.t, 'v = ', '%7.5f.' % cell.soma(0.5).v
    print ' distant apical dend v = ', '%7.5f.' % cell.dend[121](0.0).v
    print ' distant basal dend v = ', '%7.5f.' % cell.dend[153](0.0).v
    print >>f,  ' srun ended in ', steps , ' steps, t+trans = ', h.t, 'v = ', '%7.5f.' % cell.soma(0.5).v
    print >>f, ' distant apical dend v = ', '%7.5f.' % cell.dend[121](0.0).v
    print >>f, ' distant basal dend v = ', '%7.5f.' % cell.dend[153](0.0).v
    f.close()

def clampEffects(cell, v0, trans, tend, t0, tstep, clampDend):
    f = open('pylog','a')
    h.tstop = tend
    print   '     ', t0, ' with ', trans, ' trans ', ' to ', tend
    print >>f,   '     ', t0, ' with ', trans, ' trans ', ' to ', tend 

    h.dt = tstep
    v = h.Vector()
    t = h.Vector()
    v.record(cell.soma(0.5)._ref_v)
    t.record(h._ref_t)
    if trans > 0:
        cell.soma.push()
        vHold = h.SEClamp(0.5)
        vHold.dur2 = 0
        vHold.dur3 = 0
        vHold.rs = 1e-9
        vHold.dur1 = t0 + trans
        vHold.amp1 = v0
        h.pop_section()
        print ' soma clamp'
        if clampDend:
            vHolds = []
            print cell.ndend, ' dendrites'
            for i in xrange(cell.ndend):
                    cell.dend[i].push()
                    vHolds.append(h.SEClamp(0.5))
                    vHolds[i].dur2 = 0
                    vHolds[i].dur3 = 0
                    vHolds[i].rs = 1e-9
                    vHolds[i].dur1 = t0 + trans
                    vHolds[i].amp1 = v0
                    h.pop_section()
    dendv = []
    for i in xrange(cell.ndend):
        dendv.append(h.Vector())
        dendv[3*i].record(cell.dend[i](0.16666)._ref_v)
        dendv.append(h.Vector())
        dendv[3*i+1].record(cell.dend[i](0.5)._ref_v)
        dendv.append(h.Vector())
        dendv[3*i+2].record(cell.dend[i](0.83333)._ref_v)
    print   'ready to run'
    print >>f, 'ready to run'
    srun(cell, v0, trans, t0)
    #print   'i ran and i killed ', fired, ' bedbug(s)'
    #print >>f,   'i ran and i killed ', fired, ' bedbug(s)'
    ntotal = int(round((tend-t0)/tstep)+1)
    nt = int((trans)/tstep)
    dendVec = np.empty((3*cell.ndend,ntotal),dtype='float')
    for i in xrange(cell.ndend):
        dendVec[3*i,:] = dendv[3*i].as_numpy()
        dendVec[3*i+1,:] = dendv[3*i+1].as_numpy()
        dendVec[3*i+2,:] = dendv[3*i+2].as_numpy()
    v1 = v.as_numpy()
    t1 = t.as_numpy()
    dendVec = np.reshape(dendVec,(3,199,ntotal))
    mdend = np.average(dendVec, axis = 0)
    print 'mdend ',mdend.shape
    qdendlow = np.percentile(mdend,10, axis = 0)
    print 'qdendlow',qdendlow.shape
    qdendhigh = np.percentile(mdend,90, axis = 0)
    print 'qdendhigh',qdendhigh.shape
    dendAverage = np.average(mdend, axis=0)
    print 'dendAverage',dendAverage.shape
    denderrbar = np.vstack((qdendlow,qdendhigh))
    print 'denderrbar',denderrbar.shape
    pyplot.figure('slice',figsize=(8,4))
    #istart = nt
    istart = 0 
    pyplot.plot(t1[istart:],v1[istart:])
    t2 = t1[istart:]
    dA = dendAverage[istart:]
    dE = denderrbar[:,istart:]
    select = np.arange(0,t2.size,50)
    t2 = t2[select]
    print 't2',t2.shape
    dA = dA[select]
    print 'dA',dA.shape
    dE = dE[:,select]
    print 'dE',dE.shape
    pyplot.errorbar(t2,dA,np.absolute(np.vstack((dA,dA))-dE))
    pyplot.draw()
    #pyplot.savefig('slice.png',format='png',bbox_inches='tight',dpi=900)
    #print 'highest voltage ', np.amax(v1), ' in ', v1.size, ' steps'
    #print >>f, 'highest voltage ', np.amax(v1), ' in ', v1.size, ' steps'
    f.close()

if __name__ == '__main__':
    tstep = 1.0/10.0
    #dtv = [10,25,50,100,170,240]
    run_t = 200.0
    run_nt = int(round(run_t/tstep))+1
    tol_t = min([300,run_t])
    tol_nt = int(round(tol_t/tstep))+1
    seed = 231271 
    np.random.seed(seed)
    #locE = np.array([60, 72, 78, 84, 90, 98],dtype='int')
    #locI = np.array([14, 28, 30],dtype='int')
    locE = np.array([79, 82, 83, 98, 120, 124],dtype='int')
    locI = np.array([14, 28, 40],dtype='int')

    g0 = 32.0*1e-3
    #gE = np.array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0]) * g0
    #gI = -g0*np.array([10.0, 10.0, 10.0])
    gE = np.array([0.6, 0.6, 0.2, 0.6, 0.15, 0.6]) * g0
    gI = -g0*np.array([6.0, 10.0, 8.0])
    posE = np.array([0.3,0.3,0.9,0.6,0.4,0.2])
    posI = np.array([0.7,0.2,0.5])
    pos = np.concatenate((posE, posI))
    loc = np.concatenate((locE, locI))
    gList = np.concatenate((gE, gI))
    print gList
    #loc = np.array([28])
    #gList = np.array([-1e-2])    
    #pos = np.array([0.5])
    #loc = np.array([137,101])
    #gList = np.array([-6.3e-04,1.16e-03])    
    #pos = np.array([0.08,0.01])
    v0 = -70
    n = loc.size
    alphaR = True 
    cell, vecStimList, synList, distance = prepCell(gList, loc, pos, n, v0, alphaR)

    sel = range(n)

    #vecTuple = (np.array([0, run_t+1]),np.array([100,run_t+1]),np.array([200,run_t+1]),np.array([300,run_t+1]),np.array([400,run_t+1]),np.array([500,run_t+1]),np.array([600,run_t+1]),np.array([800,run_t+1]),np.array([1000,run_t+1]))
    #vecTuple = (np.array([0,20,40, run_t+1]),np.array([100,120,140,run_t+1]),np.array([200,220,240,run_t+1]),np.array([300,320,340,run_t+1]),np.array([400,420,440,run_t+1]),np.array([500,520,540,run_t+1]),np.array([600,620,640,run_t+1]),np.array([800,run_t+1]),np.array([1000,run_t+1]))
    vecTuple = [np.array([0,run_t+1]),np.array([run_t+1]),np.array([run_t+1]),np.array([run_t+1]),np.array([run_t+1]),np.array([run_t+1]),np.array([run_t+1]),np.array([run_t+1]),np.array([run_t+1])]
    vecTuple0 = [np.array([run_t+1]),np.array([run_t+1]),np.array([run_t+1]),np.array([run_t+1]),np.array([run_t+1]),np.array([run_t+1]),np.array([run_t+1]),np.array([run_t+1]),np.array([run_t+1])]
    #vecTuple = (np.array([0,330,run_t+1]),np.array([run_t+1]),np.array([run_t+1]),np.array([run_t+1]),np.array([run_t+1]),np.array([run_t+1]),np.array([run_t+1]),np.array([run_t+1]),np.array([0,660,run_t+1]))
    taue0 = 0.098814229249 
    RList0 = [[0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0]]
    RList = [[0, gE[0]/taue0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0]]
    dendVclamp = np.array([1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000])
    #dendVclamp = np.array([-50, -50, -50, -50, -50, -50, 1000, 1000, 1000])
    t0 = 0.0
    tref = 13
    vThres = -65.0
    v_pre = -70
    #v_pre = vThres
    vBack = -vThres - 1
    printR = False
    getDendV = False 
    #RList = np.zeros((n,2))
    oneGo = True 
    clampDend = False
    monitorDend = False 
    pas = False 
    plotSv = True 
    plotSleak = False
    vRange = np.array([-74,-70,-67,-65,-63,-62,-61,-60,-59,-58],dtype='double')
    trans = 0
    vtrace = np.empty((run_nt,vRange.size))
    leakyV = np.empty((run_nt,vRange.size))
    pyplot.figure('v',figsize=((8,4)))
    ax1 = pyplot.subplot(1,2,1)
    ax2 = pyplot.subplot(1,2,2)
    #for i in xrange(vRange.size):
    for i in xrange(1):
        v_pre = vRange[i]
        cell.set_volt(v_pre)
        vtrace[:,i], fired, tsp, ntrans = proceed(cell, v_pre, synList, RList0, vecStimList, vecTuple, n, trans, trans+run_t, vBack, tref, vThres, oneGo, t0, tstep, loc, pos, dendVclamp, alphaR, getDendV, monitorDend, pas, plotSv, 'v-'+str(-v_pre), True)
        leakyV[:,i], fired, tsp, ntrans = proceed(cell, v_pre, synList, RList0, vecStimList, vecTuple0, n, trans, trans+run_t, vBack, tref, vThres, oneGo, t0, tstep, loc, pos, dendVclamp, alphaR, getDendV, monitorDend, pas, plotSleak, 'leak-'+str(-v_pre), True)
        ax1.plot(vtrace[:,i],'b',lw = i*0.25)
        ax1.plot(leakyV[:,i],'k',lw=0.25)
        ax2.plot(vtrace[:,i]-leakyV[:,i],'r',lw = i*0.25)
    pyplot.savefig('v.png',format='png',bbox_inches='tight',dpi=900)
    pyplot.close()

    vtrace0 = np.empty((run_nt,vRange.size))
    leakyV0 = np.empty((run_nt,vRange.size))
    pyplot.figure('v0',figsize=((8,4)))
    ax1 = pyplot.subplot(1,2,1)
    ax2 = pyplot.subplot(1,2,2)
    trans = 110.0
    ntrans = int(round(trans/tstep))
    for i in xrange(vRange.size):
        v_pre = vRange[i]
        vtmp, fired, tsp, ntrans = proceed(cell, v_pre, synList, RList0, vecStimList, vecTuple, n, trans, trans+run_t, vBack, tref, vThres, oneGo, t0, tstep, loc, pos, dendVclamp, alphaR, getDendV, monitorDend, pas, plotSv, 'v0-'+str(-v_pre))
        vtrace0[:,i] = vtmp[ntrans:]
        vtmp, fired, tsp, ntrans = proceed(cell, v_pre, synList, RList0, vecStimList, vecTuple0, n, trans, trans+run_t, vBack, tref, vThres, oneGo, t0, tstep, loc, pos, dendVclamp, alphaR, getDendV, monitorDend, pas, plotSleak, 'leak0-'+str(-v_pre))
        leakyV0[:,i] = vtmp[ntrans:]
        ax1.plot(vtrace0[:,i],'b',lw = i*0.25)
        ax1.plot(leakyV0[:,i],'k',lw=0.25)
        ax2.plot(vtrace0[:,i]-leakyV0[:,i],'r',lw = i*0.25)
    pyplot.savefig('v0.png',format='png',bbox_inches='tight',dpi=900)
    pyplot.close()
    pyplot.figure('dv',figsize=((8,4)))
    ax1 = pyplot.subplot(1,2,1)
    ax2 = pyplot.subplot(1,2,2)
    for i in xrange(vRange.size):
        ax1.plot(leakyV[:,i]-leakyV0[:,i],'k',lw=0.25*i)
        ax2.plot(vtrace[:,i]-leakyV[:,i]-vtrace0[:,i]+leakyV0[:,i],'r',lw=0.25*i)
    pyplot.savefig('dv.png',format='png',bbox_inches='tight',dpi=900)
