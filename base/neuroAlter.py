from neuron import h
h.load_file('stdlib.hoc')
h.load_file('stdrun.hoc')
import numpy as np
from n128 import *
import sys
import matplotlib
matplotlib.use('Agg')
#matplotlib.use('Qt5Agg')
from matplotlib import pyplot
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
        print "decayTau = ", 1/self.b
        print >>f, "decayTau = ", 1/self.b
        print "Rtau = ", self.Rtau
        print >>f, "Rtau = ", self.Rtau
        print "Rinf = ", self.Rinf
        print >>f, "Rinf = ", self.Rinf
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

def run(cell, v0, vBack, tref, vThres, synList, RList, n, trans, oneGo, t0, printR = False, alphaR = True, loc=np.array([]), pos=np.array([]), pas = False, v=h.Vector, t=h.Vector()):
    f = open('pylog','a')
    nc = 0
    cell.init()
    print " cell initiated"
    print >>f, " cell initiated"
    steps = 0
    tsps = np.empty(int((h.tstop-t0)/h.dt/2))
    tsp = -tref
    #plusOne = 0
    vold = v0
    firing = 0 
    h.t = t0
    print  "t0 = ", h.t
    print >>f,  "t0 = ", h.t
    while round(h.t/h.dt) < round((t0+trans)/h.dt):
        h.fadvance()
        steps = steps + 1
        #print  cell.soma(0.5).v
        #print >>f,  steps, h.t
    for i in xrange(n):
        #print  str(i), "th's Synapse = ", RList[i]
        #print >>f,  str(i), "th's Synapse = ", RList[i]
        if alphaR:
            synList[i].g = RList[i][0]
            synList[i].h = RList[i][1]
        else:
            synList[i].R = RList[i]
    #oldProgress = 0
    print  trans, "ms trans complete, used ", steps , ' steps, t+trans = ', h.t, 'v = ', '%7.5f.' % cell.soma(0.5).v
    print >>f,  trans, "ms trans complete, used ", steps , ' steps, t+trans = ', h.t, 'v = ', '%7.5f.' % cell.soma(0.5).v
    if len(loc)>0:
        for i in xrange(n): 
            print " dend[", loc[i], "]" ,"(",pos[i],").v = ", cell.dend[loc[i]](pos[i]).v
            print >> f, " dend[", loc[i], "]" ,"(",pos[i],").v = ", cell.dend[loc[i]](pos[i]).v
    checkList = [];
    if not pas:
        if oneGo:
            while (round(h.t/h.dt) < round(h.tstop/h.dt)):
                h.fadvance()
                if cell.soma(0.5).v > vThres+15 and cell.soma(0.5).v < vold and not firing:
                    tsp = h.t-trans
                    print "nc1 = ", nc
                    print "tsps size = ", tsps.size 
                    tsps[nc] = tsp
                    firing = 1
                    nc = nc + 1
                    print  tsp + trans, "fired expecting finish ", tsp + tref + trans
                    print >>f,  tsp + trans, "fired expecting finish ", tsp + tref + trans
                vold = cell.soma(0.5).v
                if len(checkList) > 0 and len(loc) > 0:
                    for checkTime in checkList:
                        if abs(h.t-checkTime-trans) < 1e-10:
                            print " check at ", checkTime
                            for i in xrange(n): 
                                print "     dend[", loc[i], "]" ,"(",pos[i],").v = ", cell.dend[loc[i]](pos[i]).v
                                print >> f, "   dend[", loc[i], "]" ,"(",pos[i],").v = ", cell.dend[loc[i]](pos[i]).v
                if (cell.soma(0.5).v <= vThres+7.5):
                    firing = 0
        else:
            while ((cell.soma(0.5).v > vBack or h.t < trans + tsp + tref or firing) and round(h.t/h.dt) < round(h.tstop/h.dt)):
                h.fadvance()
                if cell.soma(0.5).v > vThres+15 and cell.soma(0.5).v < vold and not firing:
                    tsp = h.t-trans
                    print "nc = ", nc
                    print "tsps size = ", tsps.size 
                    tsps[nc] = tsp
                    firing = 1
                    nc = nc + 1
                    print  tsp + trans, "fired expecting finish ", tsp + tref + trans
                    print >>f,  tsp + trans, "fired expecting finish ", tsp + tref + trans
                vold = cell.soma(0.5).v
                if (cell.soma(0.5).v <= vThres+7.5):
                    firing = 0
    else:
        assert(oneGo==1)
        while (round(h.t/h.dt) < round(h.tstop/h.dt)):
            if cell.soma(0.5).v >= vThres:
                tsp = h.t-trans
                print "nc1 = ", nc
                print "tsps size = ", tsps.size 
                tsps[nc] = tsp
                nc = nc + 1
                print  tsp + trans, "fired expecting finish ", tsp + tref + trans
                print >>f,  tsp + trans, "fired expecting finish ", tsp + tref + trans

            if len(checkList) > 0 and len(loc) > 0:
                for checkTime in checkList:
                    if abs(h.t-checkTime-trans) < 1e-10:
                        print " check at ", checkTime
                        for i in xrange(n): 
                            print "     dend[", loc[i], "]" ,"(",pos[i],").v = ", cell.dend[loc[i]](pos[i]).v
                            print >> f, "   dend[", loc[i], "]" ,"(",pos[i],").v = ", cell.dend[loc[i]](pos[i]).v

            if h.t < trans + tsp + tref:
                for sec in cell.all: 
                    sec.v = cell.Vrest
                h.t = h.t + h.dt
                v.append(cell.soma(0.5).v)
                t.append(h.t)
            else:
                h.fadvance()

    print "((", cell.soma(0.5).v > vBack, " or ", h.t < trans + tsp + tref, " or ", firing > 0, ") and ", round(h.t/h.dt) < round(h.tstop/h.dt), ")"
    print >> f, "((", cell.soma(0.5).v > vBack, " or ", h.t < trans + tsp + tref, " or ", firing > 0, ") and ", round(h.t/h.dt) < round(h.tstop/h.dt), ")"

    print  "stopping with v ", cell.soma(0.5).v, " < ", vBack, " firing = ", firing, " or t ", h.t, " == ", h.tstop
    print  >> f, "stopping with v ", cell.soma(0.5).v, " < ", vBack, " firing = ", firing, " or t ", h.t, " == ", h.tstop
    f.close()
    return nc, tsps[:nc]

def prepCell(gList, loc, pos, n, v0, alphaR = True):
    f = open('pylog','a')
    print "    location,  strength,   pos"
    print >>f, "    location,  strength,   pos"
    for i in xrange(n):
        print '#{:2d}   {:3d}      {: 1.1e}    {:0.1f}'.format(i,loc[i],gList[i],pos[i])
        print >>f, '#{:2d}   {:3d}      {: 1.1e}    {:0.1f}'.format(i,loc[i],gList[i],pos[i])
    h.cvode.active(0)
    cell = RealisticNeuron(v0)
    cell.soma.push()
    h.distance()
    h.pop_section()
    vecStimList = [h.VecStim() for i in xrange(n)]
    synList = []
    for i in xrange(n):
        print i
        print >>f, i
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
        h.setpointer(vecStimList[i]._ref_y,"pre",synList[i])
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
    f.close()
    return cell, vecStimList, synList

def proceed(cell, v0, synList, RList, vecStimList, spikeTrain, n, trans, tend, vBack, tref, vThres, oneGo, t0, tstep, loc=np.array([]), pos=np.array([]), dendVclamp=np.array([]), printR = False, alphaR = True, getDendV = False, monitorDend = False, pas = False):
    f = open('pylog','a')
    h.tstop = tend
    print   "     ", t0, ' with ', trans, ' trans ', ' to ', tend , ', tref: ', tref
    print >>f,   "     ", t0, ' with ', trans, ' trans ', ' to ', tend , ', tref: ', tref
    print   "     ", "vBack: ", vBack, " vThres: ", vThres, " vStart ", v0
    print >>f,   "     ", "vBack: ", vBack, " vThres: ", vThres, " vStart ", v0
    print   "      RList:", RList
    print >>f,   "      RList:", RList
    print   "      spikeTrain:", spikeTrain
    print >>f,   "      spikeTrain:", spikeTrain
    print "oneGo: ", oneGo
    print >>f, "oneGo: ", oneGo
    print loc
    print >>f, loc
    print pos
    print >>f, pos
    print dendVclamp 
    print >>f, dendVclamp 
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
        print " soma clamp"
        vHolds = []
        if not oneGo:
            if not monitorDend:
                j = 0
                for i in xrange(len(loc)):
                    if not abs(dendVclamp[i] -1000) < 1e-10:
                        cell.dend[loc[i]].push()
                        vHolds.append(h.SEClamp(pos[i]))
                        vHolds[j].dur2 = 0
                        vHolds[j].dur3 = 0
                        vHolds[j].rs = 1e-9
                        vHolds[j].dur1 = t0 + trans
                        vHolds[j].amp1 = dendVclamp[i]
                        j = j + 1
                        h.pop_section()
                        print " dend[",loc[i],"] clamped"
            else:
                print cell.ndend, " dendrites clampled to ", v0
                for i in xrange(cell.ndend):
                    cell.dend[i].push()
                    vHolds.append(h.SEClamp(0.5))
                    vHolds[i].dur2 = 0
                    vHolds[i].dur3 = 0
                    vHolds[i].rs = 1e-9
                    vHolds[i].dur1 = t0 + trans
                    vHolds[i].amp1 = cell.Vrest + (v0 - cell.Vrest) * 1.0
                    h.pop_section()
                for i in xrange(len(loc)):
                    if not abs(dendVclamp[i] -1000) < 1e-10:
                        cell.dend[loc[i]].push()
                        vHolds.append(h.SEClamp(pos[i]))
                        vHolds[loc[i]].dur2 = 0
                        vHolds[loc[i]].dur3 = 0
                        vHolds[loc[i]].rs = 1e-9
                        vHolds[loc[i]].dur1 = t0 + trans
                        vHolds[loc[i]].amp1 = dendVclamp[i]
                        h.pop_section()
                        print " dend[",loc[i],"] specifically clamped to ", dendVclamp[i]
    dendv = []
    if oneGo and monitorDend and __name__ == '__main__':
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
                print "played ", i, " s: ", synList[i].g
                print >>f, "played ", i, " s: ", synList[i].g
            else:
                print "played ", i, " s: ", synList[i].gmax
                print >>f, "played ", i, " s: ", synList[i].gmax
            if alphaR:
                synList[i].deltat = h.dt
    print   "ready to run"
    print >>f, "ready to run"
    fired, tsp = run(cell, v0, vBack, tref, vThres, synList, RList, n, trans, oneGo, t0, printR, alphaR, loc, pos, pas, v, t)
    #print   "i ran and i killed ", fired, " bedbug(s)"
    #print >>f,   "i ran and i killed ", fired, " bedbug(s)"
    v1 = v.as_numpy()
    if oneGo or round(h.t/h.dt) > round(h.tstop/h.dt):
        ntotal = int(round((tend-t0)/tstep)+1)
    else:
        ntotal = int(round((h.t-t0)/tstep)+1)
    #if __name__ == '__main__':
    t1 = t.as_numpy()
    if v1.size!=ntotal:
        print " steps ", v1.size, ", ntotal ", ntotal
        print " remove extra step"
        print >>f, " remove extra step"
        v1 = v1[:-1]
        #if __name__ == '__main__':
        t1 = t1[:-1]
    assert(t1.size == v1.size)
    print   " is this nan ", v1[-1], v1[0]
    print >>f,   " is this nan ", v1[-1], v1[0]
    print   " v size ", v1.size
    print >>f,   " v size ", v1.size
    ntrans = int(round(trans/tstep))
    print " trans size ", ntrans
    print   " trans end, v0, v1", v1[ntrans-1], v1[ntrans], v1[-1]
    print >>f,   " trans end, v0, v1", v1[ntrans-1], v1[ntrans], v1[-1]
    print   " v w/o trans size ", v1.size - ntrans
    print >>f,   " v w/o trans size ", v1.size - ntrans
    if trans > 0:
        vHold = None

    if oneGo and __name__ == '__main__':
        pyplot.figure('slice',figsize=(8,4))
        if not monitorDend:
            pyplot.plot(t1[ntrans:]-trans,v1[ntrans:])
        else:
            nt = ntrans
            dendVec = np.empty((3*cell.ndend,ntotal),dtype='float')
            for i in xrange(cell.ndend):
                dendVec[3*i,:] = dendv[3*i].as_numpy()
                dendVec[3*i+1,:] = dendv[3*i+1].as_numpy()
                dendVec[3*i+2,:] = dendv[3*i+2].as_numpy()
            v1 = v.as_numpy()
            t1 = t.as_numpy()
            apicalDend = (dendv[3*121].as_numpy() + dendv[3*121+1].as_numpy() + dendv[3*121+2].as_numpy())/3.0
            basalDend = (dendv[3*153].as_numpy() + dendv[3*153+1].as_numpy() + dendv[3*153+2].as_numpy())/3.0
            dendVec = np.reshape(dendVec,(3,199,ntotal))
            mdend = np.average(dendVec, axis = 0)
            print 'mdend ',mdend.shape
            #qdendlow = np.percentile(mdend,10, axis = 0)
            qdendlow = np.amin(mdend, axis = 0)
            print 'qdendlow',qdendlow.shape
            #qdendhigh = np.percentile(mdend,90, axis = 0)
            qdendhigh = np.amax(mdend, axis = 0)
            print 'qdendhigh',qdendhigh.shape
            dendAverage = np.average(mdend, axis=0)
            print 'dendAverage',dendAverage.shape
            denderrbar = np.vstack((qdendlow,qdendhigh))
            print 'denderrbar',denderrbar.shape
            #istart = 0 
            istart = nt 
            t2 = t1[istart:]
            pyplot.plot(t2-istart*tstep,v1[istart:])
            dA = dendAverage[istart:]
            dE = denderrbar[:,istart:]
            select = np.arange(0,t2.size,20)
            t2 = t2[select] - istart*tstep
            print 't2',t2.shape
            dA = dA[select]
            print 'dA',dA.shape
            dE = dE[:,select]
            print 'dE',dE.shape
            pyplot.errorbar(t2,dA,np.absolute(np.vstack((dA,dA))-dE))
            qdendlow = np.percentile(mdend,5, axis = 0)
            qdendhigh = np.percentile(mdend,95, axis = 0)
            denderrbar = np.vstack((qdendlow,qdendhigh))
            dE = denderrbar[:,istart:]
            dE = dE[:,select]
            pyplot.errorbar(t2+0.2,dA,np.absolute(np.vstack((dA,dA))-dE),ls = 'None')
            qdendlow = np.percentile(mdend,25, axis = 0)
            qdendhigh = np.percentile(mdend,75, axis = 0)
            denderrbar = np.vstack((qdendlow,qdendhigh))
            dE = denderrbar[:,istart:]
            dE = dE[:,select]
            pyplot.errorbar(t2+0.4,dA,np.absolute(np.vstack((dA,dA))-dE),ls = 'None')
            qdendlow = np.percentile(mdend,45, axis = 0)
            qdendhigh = np.percentile(mdend,55, axis = 0)
            denderrbar = np.vstack((qdendlow,qdendhigh))
            dE = denderrbar[:,istart:]
            dE = dE[:,select]
            pyplot.errorbar(t2+0.6,dA,np.absolute(np.vstack((dA,dA))-dE),ls = 'None')
        pyplot.savefig('slice.png',format='png',bbox_inches='tight',dpi=900)
        
    #print "highest voltage ", np.amax(v1), " in ", v1.size, " steps"
    #print >>f, "highest voltage ", np.amax(v1), " in ", v1.size, " steps"
    f.close()
    if getDendV:
        dendVecOut = np.array([])
        if oneGo and monitorDend and __name__ == '__main__':
            for i in xrange(n):
                dendVecOut = np.concatenate((dendVecOut,dendv[loc[i]].as_numpy().copy()))
        else:
            for i in xrange(n):
                dendVecOut = np.concatenate((dendVecOut,dendv[i].as_numpy().copy()))
        print "returning dendV", dendVecOut.size
        return v1.copy(), fired, tsp, ntrans, dendVecOut
    else:
        print "returning somaV only"
        return v1.copy(), fired, tsp, ntrans

def leaky(cell, v0, synList, RList, vecStimList, n, trans, tend, tstep, loc, pos, printR = False, alphaR = True):
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
        print " soma clamp"
    dendv = []
    for i in xrange(n):
        dendv.append(h.Vector())
        dendv[i].record(cell.dend[loc[i]](pos[i])._ref_v)
        print "record ", i
        print >>f, "record ", i
        if alphaR:
            synList[i].deltat = h.dt
        vecStimList[i].play(h.Vector([]))
        vecStimList[i].dt = h.dt
    print   "ready to run"
    print >>f, "ready to run"
    run(cell, v0, 0, 0, 0, synList, RList, n, trans, 1, 0, printR, alphaR, loc, pos)
    #print   "i ran and i killed ", fired, " bedbug(s)"
    #print >>f,   "i ran and i killed ", fired, " bedbug(s)"
    v1 = v.as_numpy()
    ntotal = int(round(tend/tstep)+1)
    if __name__ == '__main__':
        t1 = t.as_numpy()
    if v1.size!=ntotal:
        print " steps ", v1.size, ", ntotal ", ntotal
        print " remove extra step"
        print >>f, " remove extra step"
        v1 = v1[:-1]
        if __name__ == '__main__':
            t1 = t1[:-1]
    print   " is this nan ", v1[-1], v1[-2], v1[-3]
    print >>f,   " is this nan ", v1[-1], v1[-2], v1[-3]
    print   " v size ", v1.size
    print >>f,   " v size ", v1.size
    ntrans = int(round(trans/tstep))
    #if v1[ntrans] == v1[ntrans-1]:
    #    ntrans = ntrans+1
    #assert(v1[ntrans] != v1[ntrans-1])
    print " trans size ", ntrans
    print   " trans end, v0, v1", v1[ntrans-1], v1[ntrans], v1[ntrans+1]
    print >>f,   " trans end, v0, v1", v1[ntrans-1], v1[ntrans], v1[ntrans+1]
    print   " v w/o trans size ", v1.size - ntrans
    print >>f,   " v w/o trans size ", v1.size - ntrans
    if trans > 0:
        vHold = None
    if __name__ == '__main__':
        nt = int((trans)/tstep)
        pyplot.figure('slice',figsize=(8,4))
        pyplot.plot(t1[nt:],v1[nt:])
        #pyplot.show()
        pyplot.savefig('slice.png',format='png',bbox_inches='tight',dpi=900)
    #print "highest voltage ", np.amax(v1), " in ", v1.size, " steps"
    #print >>f, "highest voltage ", np.amax(v1), " in ", v1.size, " steps"
    dendvArray = np.empty((n,ntotal))
    for i in xrange(n):
        dendvArray[i,:] = dendv[i].as_numpy().copy()[:ntotal]
    f.close()
    return v1.copy(), dendvArray

def sproceed(cell, v0, synList, gList, RList, vecStimList, spike, n, sel, trans, tend, tstep, name, printR = False, alphaR = True):
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

    fired, _ = run(cell, v0, 0, 0, 0, synList, RList, n, trans, 1, t0, printR, alphaR)
    if fired > 0:
        print "I fired, but I'm not supposed to"
    #print "i ran and i killed ", fired, " bedbug(s)"
    #print >>f, "i ran and i killed ", fired, " bedbug(s)"
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
    #print "highest voltage ", np.amax(v1), " in ", v1.size, " steps"
    #print >>f, "highest voltage ", np.amax(v1), " in ", v1.size, " steps"
    f.close()
    return v1.copy(), fired

def bproceed0(cell, v0, synList, gList, RList, vecStimList, spikeTrain, n, sel, trans, dt, tend, t0, tstep, name, pos=-1.0, printR = False, alphaR = True):
    f = open('pylog','a')
    h.tstop = tend
    print t0, ' with ', trans, ' trans ', ' to ', tend, ' v0 = ', v0

    print >>f, t0, ' with ', trans, ' trans ', ' to ', tend, ' v0 = ', v0
    print gList
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
            #h.setpointer(vecStimList[i]._ref_y,"pre",synList[i])
            vecStimList[i].dt = h.dt
            if alphaR:
                synList[i].deltat = h.dt
                synList[i].f = 0
            else:
                synList[i].gmax = 0

    fired, _ = run(cell, v0, 0, 0, 0, synList, RList, n, trans + dt, 1, t0, printR, alphaR)
    #print "i ran and i killed ", fired, " bedbug(s)"
    #print >>f, "i ran and i killed ", fired, " bedbug(s)"
    if fired > 0:
        print "I fired, but I'm not supposed to"
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
    #print "highest voltage ", np.amax(v1), " in ", v1.size, " steps"
    #print >>f, "highest voltage ", np.amax(v1), " in ", v1.size, " steps"
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

def bproceed(cell, v0, synList, gList, RList, vecStimList, spikeTrain, n, sel, trans, tend, t0, tstep, name, pos=-1.0, printR = False, alphaR = True):
    f = open('pylog','a')
    h.tstop = tend
    print t0, ' with ', trans, ' trans ', ' to ', tend, ' v0 = ', v0

    print >>f, t0, ' with ', trans, ' trans ', ' to ', tend, ' v0 = ', v0
    print gList
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
            #h.setpointer(vecStimList[i]._ref_y,"pre",synList[i])
            vecStimList[i].dt = h.dt
            if alphaR:
                synList[i].deltat = h.dt
                synList[i].f = 0
            else:
                synList[i].gmax = 0

    fired, _ = run(cell, v0, 0, 0, 0, synList, RList, n, trans, 1, t0, printR, alphaR)
    #print "i ran and i killed ", fired, " bedbug(s)"
    #print >>f, "i ran and i killed ", fired, " bedbug(s)"
    if fired > 0:
        print "I fired, but I'm not supposed to"
    #assert(fired == 0)
    v1 = v.as_numpy()
    if printR:
        if alphaR:
            g1 = g.as_numpy()
            h1 = hh.as_numpy()
        else:
            r1 = R.as_numpy()

    if trans > 0:
        vHold = None
    if name:
        t1 = t.as_numpy()
        nt = int((trans)/tstep)
        pyplot.figure(name,figsize=(8,4))
        pyplot.plot(t1[nt:],v1[nt:])
        pyplot.savefig(name+'.png',format='png',bbox_inches='tight',dpi=900)
        pyplot.close()
    #print "highest voltage ", np.amax(v1), " in ", v1.size, " steps"
    #print >>f, "highest voltage ", np.amax(v1), " in ", v1.size, " steps"
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

def srun(cell, v0, trans, t0):
    f = open('pylog','a')
    cell.init()
    print " cell initiated"
    print >>f, " cell initiated"
    steps = 0
    h.t = t0
    print  "t0 = ", h.t
    print >>f,  "t0 = ", h.t
    while round(h.t/h.dt) < round((t0+trans)/h.dt):
        h.fadvance()
        steps = steps + 1
    print trans, "ms trans complete, used ", steps , ' steps, t+trans = ', h.t, 'v = ', '%7.5f.' % cell.soma(0.5).v
    print' distant apical dend v = ', '%7.5f.' % cell.dend[121](0.0).v
    print' distant basal dend v = ', '%7.5f.' % cell.dend[153](0.0).v
    print >>f,  trans, "ms trans complete, used ", steps , ' steps, t+trans = ', h.t, 'v = ', '%7.5f.' % cell.soma(0.5).v
    print >>f, ' distant apical dend v = ', '%7.5f.' % cell.dend[121](0.0).v
    print >>f, ' distant basal dend v = ', '%7.5f.' % cell.dend[153](0.0).v
    while (h.t < h.tstop):
        h.fadvance()
        steps = steps + 1
    print "srun ended in ", steps , ' steps, t+trans = ', h.t, 'v = ', '%7.5f.' % cell.soma(0.5).v
    print ' distant apical dend v = ', '%7.5f.' % cell.dend[121](0.0).v
    print ' distant basal dend v = ', '%7.5f.' % cell.dend[153](0.0).v
    print >>f,  " srun ended in ", steps , ' steps, t+trans = ', h.t, 'v = ', '%7.5f.' % cell.soma(0.5).v
    print >>f, ' distant apical dend v = ', '%7.5f.' % cell.dend[121](0.0).v
    print >>f, ' distant basal dend v = ', '%7.5f.' % cell.dend[153](0.0).v
    f.close()

def clampEffects(cell, v0, trans, tend, t0, tstep, clampDend):
    f = open('pylog','a')
    h.tstop = tend
    print   "     ", t0, ' with ', trans, ' trans ', ' to ', tend
    print >>f,   "     ", t0, ' with ', trans, ' trans ', ' to ', tend 

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
        print " soma clamp"
        if clampDend:
            vHolds = []
            print cell.ndend, " dendrites"
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
    print   "ready to run"
    print >>f, "ready to run"
    srun(cell, v0, trans, t0)
    #print   "i ran and i killed ", fired, " bedbug(s)"
    #print >>f,   "i ran and i killed ", fired, " bedbug(s)"
    ntotal = int(round((tend-t0)/tstep)+1)
    nt = int((trans)/tstep)
    dendVec = np.empty((3*cell.ndend,ntotal),dtype='float')
    for i in xrange(cell.ndend):
        dendVec[3*i,:] = dendv[3*i].as_numpy()
        dendVec[3*i+1,:] = dendv[3*i+1].as_numpy()
        dendVec[3*i+2,:] = dendv[3*i+2].as_numpy()
    v1 = v.as_numpy()
    t1 = t.as_numpy()
    apicalDend = (dendv[3*121].as_numpy() + dendv[3*121+1].as_numpy() + dendv[3*121+2].as_numpy())/3.0
    basalDend = (dendv[3*153].as_numpy() + dendv[3*153+1].as_numpy() + dendv[3*153+2].as_numpy())/3.0
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
    #print "highest voltage ", np.amax(v1), " in ", v1.size, " steps"
    #print >>f, "highest voltage ", np.amax(v1), " in ", v1.size, " steps"
    f.close()

if __name__ == '__main__':
    tstep = 1.0/10.0
    run_t = 7500
    #dtv = [10,25,50,100,170,240]
    run_nt = int(round(run_t/tstep))
    tol_t = min([300,run_t])
    tol_nt = int(round(tol_t/tstep))+1
    seed = 231271 
    np.random.seed(seed)
    #locE = np.array([36, 53, 74, 79, 83, 90, 97, 101, 117, 125, 132, 174],dtype='int')
    #locE = np.array([36, 74, 83, 97, 117, 132],dtype='int')
    #locE = np.array([36, 73, 97, 117, 146, 180],dtype='int')
    #locE = np.array([79, 82, 83, 108, 124, 129],dtype='int')
    locE = np.array([60, 72, 78, 84, 90, 98],dtype='int')
    #locI = np.array([14, 28, 30],dtype='int')
    locI = np.array([14, 28, 30],dtype='int')
    #locE = np.array([32, 52, 66, 78, 98, 136],dtype='int')
    #locE = np.array([74],dtype='int')
    #gE = 1e-4 + np.random.random_sample(locE.size) * (1e-4-1e-4) + 2e-2 + 1e-3*3
    #gE = 3.2e-2 * (1 + np.random.random_sample(locE.size) * 0.2 - 0.1)
    g0 = 3.2e-3
    gE = np.repeat(g0,locE.size)
    posE = np.random.random_sample(locE.size)
    #locI = np.array([2, 14, 28],dtype='int')
    #locI = np.array([7, 28, 137],dtype='int')
    #locI = np.array([28],dtype='int')
    #gI = (1e-1 + np.random.random_sample(locI.size) * (1-1e-1)) * (-0.1/2.0)
    gI = -np.repeat(g0,locI.size)
    posI = np.random.random_sample(locI.size)
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
    pos = np.array([0.618941, 0.647016, 0.57058, 0.457597, 0.526949, 0.796418, 0.391948, 0.975584, 0.130033])
    gList = np.array([0.012822, 0.00512783, 0.00986785, 0.0205335, 0.00415574, 0.0232459, -0.0123107, -0.0158747, -0.0376081])
    alphaR = True 
    cell, vecStimList, synList = prepCell(gList, loc, pos, n, v0, alphaR)

    sel = range(n)
    vecTuple = (np.array([3206.47, 3208.91, 3234.75, 3254.99, 3264.1, 3267.28, 3271.93, 3277.8, 3280.03, 3281.66, 3287.34, 3292.85, 3300.56, 3311.13, 3313.26, 3313.86, 3319.93, 3320.16, 3325.66, 3325.76, 3329.03, 3329.21, 3341.12, 3345, 3347.3, 3359.31, 3388.4, 3400.51, 3410.81, 3413.87, 3415.52, 3421.92, 3426.51, 3436.44, 3449.42, 3467.87, 3473.44, 3476.58, 3477.23, 3483.63, 3487.05, 3490.06, 3493.53, 3511.94, 3519.2, 3521.27, 3522.13, 3523.88, 3533.84, 3553.75, 3554.87, 3557.49, 3568, 3578.19, 3580.62, 3581.47, 3584.2, 3590.15, 3610.06, 3620.66, 3628.68, 3635.27, 3644.41, 3646.38, 3654.82, 3655.1, 3678.17, 3685.52, 3686.65, 3689.86, 3700.55, 3709.94, 3720.94, 3721.07, 3727.73, 3731.63, 3733.2, 3736.88, 3749.35, 3754.1, 3765.39, 3791.34, 3791.48, 3794.25, 3795.23, 3796.38, 3810.26, 3819.35, 3830.03, 3835.11, 3851.17, 3854.69, 3865.37, 3883.09, 3899.48, 3908.29, 3914.82, 3920.71, 3942.6, 3964.13, 3964.38, 3968.71, 3973.45, 3978.43, 3982.71, 3989.35, 3991.15, 3996.98, 4008.14, 4021.46, 4025.62, 4046.91, 4057.94, 4070.77, 4082.96, 4084.01, 4084.52, 4086.87, 4103.97, 4116.22, 4151.28, 4151.97, 4178.21, 4190.87, 4191.13, 4200.34, 4207.81, 4218.03, 4219.86, 4240.32, 4248.77, 4249.66, 4263.84, 4263.85, 4277.74, 4284.06, 4286.32, 4295.73, 4296.17, 4299.36, 4312, 4318.89, 4326.79, 4333.73, 4356.05, 4373.01, 4375.65, 4378.01, 4422.3, 4431.1, 4438.66, 4438.73, 4450.57, 4469.62, 4473.6, 4481.22, 4483.48, 4491.93, 4493.17, 4499.16, 4503.04, 4526.16, 4547.56, 4552.89, 4569.87, 4578.68, 4582.58, 4598.05, 4618.21, 4623.17, 4634.88, 4653.97, 4658.73, 4664.1, 4664.4, 4668.57, 4669.15, 4670.32, 4674.86, 4689.54, 4692.46, 4716.17, 4730.74, 4734.07, 4735.95, 4740.44, 4755.4, 4760.95, 4773.88, 4792, 4859.22, 4866.43, 4867.47, 4875.85, 4879.99, 4890.3, 4898.4, 4910.87, 4934.98, 4942.77, 4954.64, 4964.72, 4967.53, 4968.87, 4975.44, 4976.42, 4978.37, 4982.05, 4985.45, 4998.35, 4999.26, 4999.49, 5000.85, 5001.55, 5020.92, 5048.53, 5072.08, 5108.69, 5116.96, 5135.58, 5143.43, 5145.15, 5160.73, 5162.99, 5166.78, 5179.12, 5221.1, 5229.57, 5233.77, 5236.84, 5247.34, 5265.92, 5267.06, 5268.95, 5283.01, 5292.66, 5294.77, 5303.69, 5304.19, 5328.25, 5329.19, 5330.32, 5332.82, 5354.82, 5368.12, 5369.76, 5372.12, 5377.44, 5380.84, 5389.5, 5395.43, 5396.75, 5398.5, 5409.7, 5413.35, 5417.67, 5439.2, 5449.52, 5464.25, 5469.22, 5480.18, 5489.93, 5490.99, 5491.89, 5493.97, 5497.86, 5505.82, 5511.6, 5520.07, 5533.3, 5536.52, 5539.76, 5550.14, 5558.01, 5560.74, 5573.34, 5581.43, 5590.14, 5602.58, 5622.5, 5622.72, 5630.52, 5636.45, 5639.88, 5644.18, 5662.26, 5667.74, 5672.37, 5673.27, 5679.59, 5682.51, 5684.7, 5691.69, 5693.26, 5697.59, 5706.45, 5732.69, 5734.88, 5745.76, 5746.76, 5752.45, 5757.33, 5769.03, 5781.3, 5781.65, 5794.12, 5794.19, 5796.22, 5801.41, 5811.28, 5824.25, 5827.22, 5855.05, 5859.63, 5879.08, 5885.28, 5885.57, 5891.92, 5896.23, 5901.89, 5903.94, 5908.73, 5913.52, 5918.44, 5925.53, 5931.66, 5937.32, 5939.52, 5939.64, 5940.67, 5947.35, 5973.95, 5986.28, 6009.21, 6014.87, 6017.14, 6025.06, 6025.4, 6028.86, 6031.84, 6032.21, 6037.67, 6046.67, 6090.38, 6091.25, 6094.29, 6097.76, 6107.76, 6130.33, 6148.56, 6155.12, 6156.79, 6207.01, 6207.73, 6207.81, 6212.43, 6213.61, 6215.39, 6240.42, 6242.43, 6267.62, 6274.8, 6276.05, 6285.61, 6289.83, 6290.77, 6297.72, 6301.82, 6308.47, 6317.83, 6324.31, 6327.85, 6331.05, 6337.4, 6340.3, 6355.88, 6365.91, 6387.53, 6387.63, 6397.55, 6416.23, 6430.42, 6438.62, 6460.53, 6461.83, 6463.99, 6465.59, 6472.36, 6480.05, 6483.52, 6492.11, 6502.15, 6503.86, 6513.32, 6518.45, 6520.33, 6525.73, 6533.21, 6535.39, 6545.11, 6547.45, 6552.46, 6557.96, 6567.26, 6582.93, 6610.57, 6611.07, 6613.43, 6615.17, 6615.71, 6628.93, 6638.8, 6650.09, 6651.41, 6652.22, 6655.87, 6660.89, 6678.1, 6687.43, 6713.02, 6730.99, 6735.33, 6744.83, 6746.79, 6748.43, 6751.63, 6764.65, 6767.13, 6788.93, 6805.98, 6818.24, 6823.73, 6826.2, 6830.56, 6832.59, 6835.91, 6853.11, 6870.65, 6879.04, 6884.72, 6897.02, 6897.83, 6901.43, 6901.51, 6903.61, 6905.89, 6908.58, 6914.9, 6915.39, 6917.47, 6918.84, 6933.78, 6936.14, 6948.25, 6966.13, 6967.83, 6976.44, 6984.62, 7012.63, 7018.06, 7021.77, 7038.38, 7048.64, 7062.78, 7082.19, 7083.55, 7084.14, 7094, 7094.96, 7097.15, 7101.57, 7111.04, 7118.81, 7149.7, 7152.05, 7170.19, 7173.75, 7182.18, 7189.82, 7196.39, 7198.15, 7225.36, 7236.57, 7252.74, 7259.06, 7267.85, 7269.59, 7279.43, 7315.98, 7317.61, 7320.6, 7320.87, 7325.08, 7326.37, 7332.84, 7334.11, 7339.15, 7341.32, 7343.79, 7347.1, 7350.36, 7352.23, 7368.34, 7386.89, 7415.77, 7426.86, 7441.36, 7458.8, 7459.4, 7462, 7466.18, 7487.57, 7491.84, 7503.5]), np.array([3198.74, 3205.41, 3221.67, 3235.06, 3236.58, 3250.3, 3250.49, 3252.2, 3253.06, 3261.54, 3289.14, 3289.42, 3293, 3294.31, 3294.61, 3303.83, 3305.99, 3318.9, 3330.12, 3347.33, 3356.86, 3365.39, 3365.45, 3379.94, 3382.54, 3388.43, 3399.24, 3401.02, 3409.07, 3459.11, 3462.18, 3469.99, 3486.78, 3508.89, 3512.62, 3515.4, 3531.14, 3548.8, 3551.91, 3567.93, 3570.99, 3573.36, 3590.26, 3609.6, 3611.23, 3612.81, 3616.91, 3636.49, 3645.36, 3647.71, 3666.38, 3671.92, 3677.51, 3682.63, 3697.02, 3697.63, 3716.28, 3723.62, 3726.77, 3746.08, 3753.71, 3764.2, 3777.79, 3782.52, 3801.42, 3805.37, 3819.99, 3820.74, 3823.49, 3840.22, 3847.79, 3851.24, 3878.24, 3896.95, 3897.72, 3911.64, 3919.73, 3941.34, 3952.94, 3955.72, 3962.71, 3963.29, 3965.3, 3965.39, 3993.31, 3995.39, 4002.47, 4018.93, 4028.81, 4040.22, 4041.12, 4045.45, 4052.43, 4057.69, 4064.84, 4070.84, 4073.93, 4089.51, 4092.45, 4094.5, 4103.81, 4111.62, 4116.04, 4118.01, 4118.28, 4120.68, 4126.85, 4139.12, 4174.4, 4175.05, 4196.62, 4197.79, 4210.67, 4218.81, 4251.57, 4255.39, 4261.55, 4273.06, 4277.64, 4278.81, 4287.73, 4292.65, 4309.36, 4316.59, 4348.42, 4353.44, 4360.84, 4367.92, 4372.67, 4390.53, 4393.37, 4400.83, 4405.36, 4422.8, 4422.95, 4423.55, 4430.94, 4446.09, 4449.9, 4452.22, 4452.52, 4454.7, 4458.23, 4460.26, 4466, 4469.55, 4484.52, 4489.26, 4512.19, 4519.75, 4526.48, 4536.09, 4537.03, 4550.82, 4553.59, 4554.54, 4556.32, 4559.15, 4565.81, 4568.71, 4572.52, 4586.86, 4591.98, 4592.1, 4592.67, 4595.51, 4595.69, 4599.18, 4600.34, 4603.17, 4604.79, 4613.31, 4614.66, 4615.34, 4623.13, 4626.22, 4633.39, 4641.11, 4651.6, 4660.96, 4667.68, 4672.25, 4672.89, 4677.03, 4677.03, 4704.77, 4704.85, 4709.11, 4723.87, 4726.41, 4763.18, 4764.89, 4766.21, 4767.64, 4767.8, 4788.61, 4798.24, 4798.37, 4798.62, 4811.26, 4818.09, 4820.97, 4825.92, 4828.2, 4839.08, 4847.09, 4847.86, 4848.15, 4853.72, 4858.91, 4865.63, 4866.29, 4867.5, 4870.68, 4874.81, 4883.66, 4885.94, 4887.64, 4897.03, 4906.26, 4910.5, 4913.98, 4914.46, 4922.53, 4937.71, 4948.45, 4962.13, 4967.87, 4976.59, 4986.93, 5004.94, 5006.35, 5010.98, 5023.25, 5044.62, 5049.07, 5049.4, 5051.38, 5058.39, 5066.5, 5072.48, 5080.05, 5093.11, 5144, 5160.66, 5165.57, 5171.63, 5173.23, 5177.94, 5186.92, 5199.12, 5212.18, 5217.2, 5222.32, 5232.64, 5237.11, 5240.35, 5266.81, 5277.36, 5294.73, 5299.58, 5314.66, 5315.9, 5321.08, 5334.86, 5348.13, 5354.16, 5374.52, 5378.86, 5383.3, 5383.34, 5388.42, 5389.43, 5391.29, 5391.92, 5403.12, 5416.08, 5432.47, 5436.59, 5442.91, 5446.91, 5446.94, 5476.42, 5484.54, 5500.73, 5504.34, 5510.22, 5521.55, 5524.63, 5530.43, 5545.82, 5548.92, 5565.48, 5573.88, 5575.81, 5598.52, 5604.82, 5610.38, 5635.89, 5642.9, 5643.29, 5648.7, 5653.84, 5655.87, 5662.84, 5665.83, 5666.03, 5674.93, 5681.27, 5710.7, 5714.9, 5723.13, 5724.41, 5742.1, 5744.25, 5748.23, 5769.84, 5777.87, 5780.67, 5790.3, 5795.45, 5798.3, 5803.17, 5814.59, 5819.41, 5823.4, 5842.27, 5846.69, 5870.15, 5870.45, 5888.98, 5895.45, 5901.59, 5902.22, 5908.14, 5916.27, 5917.42, 5926.18, 5933.55, 5937.3, 5941.48, 5957.98, 5976.28, 5992.9, 6000.68, 6026.9, 6027.34, 6028.66, 6042.69, 6048.83, 6051.21, 6081.67, 6092.93, 6116.21, 6122.3, 6128.63, 6129.66, 6131.67, 6135.37, 6153.9, 6163.01, 6179.82, 6194.74, 6194.75, 6194.84, 6201.5, 6207.38, 6207.85, 6211.45, 6214.24, 6216.04, 6234.88, 6242.63, 6249.18, 6261.36, 6266.76, 6271.06, 6279.85, 6283.75, 6293.76, 6293.97, 6294.97, 6330.36, 6336.79, 6349.4, 6363.25, 6367.78, 6375.92, 6389.63, 6390.42, 6423.22, 6434.92, 6437.49, 6448.44, 6457.83, 6470.73, 6472.72, 6478.65, 6487.58, 6501.15, 6506.75, 6510.18, 6511.31, 6517.17, 6521.26, 6537.22, 6553.62, 6556, 6564.95, 6565.58, 6568.38, 6570.11, 6580.91, 6583.05, 6583.52, 6585.87, 6598.01, 6601.09, 6603.19, 6604.7, 6613.14, 6613.49, 6616.78, 6624.38, 6643.92, 6645.77, 6659.72, 6671.38, 6673.89, 6681.8, 6688.51, 6691.97, 6702.66, 6704.91, 6709.7, 6720.23, 6730.04, 6734.48, 6736.54, 6736.7, 6744.89, 6762.1, 6765.72, 6781.75, 6782.83, 6790.23, 6795.27, 6809.7, 6835.3, 6856.26, 6861.39, 6885.19, 6892.16, 6892.94, 6893.18, 6908.55, 6908.8, 6909.39, 6911.14, 6912.41, 6928.61, 6930.45, 6930.52, 6935.59, 6945.67, 6947.75, 6949.18, 6949.76, 6961.58, 6969.19, 6970.82, 6978.22, 7000.24, 7008.23, 7017.64, 7024.55, 7026.1, 7040.17, 7048.08, 7074.08, 7076.14, 7080.53, 7086.93, 7096.62, 7099.65, 7119.84, 7125.23, 7144.46, 7158.48, 7158.53, 7213.94, 7215.58, 7226.89, 7260.62, 7268.74, 7288.43, 7289.55, 7291.36, 7293.46, 7293.89, 7313.39, 7313.42, 7314.13, 7316.08, 7327.48, 7330.77, 7330.81, 7346.14, 7349.89, 7358.04, 7358.79, 7363.66, 7372.59, 7377.33, 7390.19, 7392.72, 7400.91, 7401.8, 7406.36, 7414.08, 7416.64, 7424.13, 7426.44, 7440.74, 7442.44, 7454.2, 7455.18, 7482.74, 7485.42, 7487.6, 7494.94, 7504.54]), np.array([3204.07, 3207.97, 3209.05, 3215.13, 3215.51, 3233.35, 3247.93, 3258.15, 3260.81, 3276.59, 3290.64, 3293.59, 3308.19, 3311.48, 3313.83, 3320.39, 3329.05, 3329.31, 3333.04, 3338.12, 3340.88, 3355.05, 3355.33, 3358.02, 3366.49, 3374.44, 3383.38, 3388.16, 3396.93, 3403.39, 3410.89, 3423.75, 3430.13, 3447.36, 3457.35, 3461.72, 3472.81, 3477.37, 3478.91, 3487.08, 3514.7, 3522.17, 3529.05, 3533.3, 3542.56, 3545.31, 3554.46, 3559.22, 3581.55, 3587.25, 3601.97, 3612.48, 3621.06, 3624.1, 3635.64, 3640.34, 3644.62, 3646.37, 3652.16, 3670.92, 3676.47, 3679.26, 3679.47, 3679.52, 3687.53, 3691.21, 3693.15, 3728.5, 3731.36, 3735.31, 3744.47, 3752.33, 3755.84, 3757.28, 3764.74, 3769.4, 3771.53, 3772.42, 3779.84, 3782.65, 3784.11, 3819.86, 3843.91, 3853.09, 3861.33, 3882.25, 3885.52, 3895.66, 3904.48, 3912.52, 3914.78, 3918.56, 3919.27, 3920.56, 3926.58, 3926.59, 3927.37, 3955.34, 3955.86, 3960.38, 3978.69, 3984.13, 3986.91, 3989.01, 3990.26, 3995.29, 4007.12, 4046.19, 4052.66, 4053.03, 4053.46, 4065.54, 4077.65, 4086.44, 4094.09, 4094.76, 4129.9, 4144.22, 4148.41, 4149.41, 4150.48, 4171.66, 4172.91, 4180.36, 4184.88, 4186.6, 4186.6, 4186.96, 4187.58, 4203.86, 4204.87, 4207.78, 4218.64, 4230.89, 4261.51, 4294.32, 4298.95, 4309.07, 4309.18, 4317.68, 4319.3, 4320.32, 4336.98, 4375.01, 4381.57, 4383.27, 4385.35, 4386.5, 4389.32, 4402.16, 4416.73, 4417.62, 4420.78, 4438.96, 4441.17, 4441.39, 4454.57, 4467.6, 4483.94, 4486.52, 4494.28, 4494.45, 4514.14, 4518.79, 4521.94, 4524.86, 4526.86, 4527.21, 4527.37, 4529.46, 4535.98, 4554.4, 4561.71, 4570.76, 4580.02, 4586.46, 4587.29, 4589.45, 4597.31, 4600.66, 4613.01, 4617.45, 4624.55, 4629.9, 4639.58, 4640.67, 4641.76, 4646.72, 4647.35, 4647.7, 4662.22, 4677.76, 4679.67, 4681.43, 4684.75, 4693.41, 4694.8, 4699.26, 4716.19, 4717.86, 4723.74, 4724.64, 4729.36, 4749.48, 4750.51, 4758.21, 4760.96, 4762.22, 4781.1, 4782, 4783.55, 4785.61, 4811.4, 4819.57, 4827.82, 4845.8, 4858.72, 4876.59, 4906.48, 4909.39, 4926.98, 4930.14, 4932.29, 4932.41, 4935.82, 4942, 4946.18, 4948.51, 4954.7, 4956.59, 4965.92, 4973.16, 4974.84, 4976.47, 5009.38, 5011.33, 5018.78, 5022.77, 5043.13, 5045.27, 5048.45, 5057.17, 5066.21, 5069.17, 5078.2, 5082.81, 5115.03, 5128.34, 5138.63, 5144.48, 5145.35, 5147.37, 5158.28, 5159.47, 5167.28, 5178.87, 5182.74, 5189.76, 5196.45, 5197.44, 5200.09, 5210.85, 5214.56, 5215.65, 5220.78, 5235.02, 5242.9, 5244.47, 5254.99, 5257.44, 5257.87, 5272.71, 5275.43, 5303.83, 5316.04, 5316.36, 5325.96, 5337.72, 5344.69, 5350.97, 5352.34, 5356.15, 5356.59, 5358.24, 5378.58, 5380.59, 5389.85, 5409.73, 5412.23, 5412.33, 5413.22, 5419.31, 5437.58, 5448.48, 5450.88, 5453.04, 5453.47, 5456.05, 5464.39, 5477.77, 5478.97, 5500.41, 5504.01, 5504.5, 5508.66, 5511.17, 5527.77, 5536.39, 5538.69, 5555.34, 5556.29, 5557.48, 5560.35, 5575.7, 5580.36, 5581.29, 5590.73, 5614.97, 5616.97, 5631.26, 5637.28, 5648.17, 5657.93, 5662.43, 5667.22, 5667.47, 5670.22, 5670.44, 5677.23, 5678.23, 5686.54, 5713.14, 5718.71, 5724.88, 5728.47, 5738.11, 5738.87, 5754.67, 5757.41, 5773.36, 5780.52, 5781.41, 5803.87, 5812.09, 5812.98, 5820.9, 5826.23, 5830.85, 5836.16, 5848.23, 5852.23, 5852.54, 5853.37, 5877.63, 5883.6, 5885.08, 5891.29, 5893.77, 5914.06, 5919.59, 5920.68, 5925.97, 5945.66, 5993.92, 6007.19, 6013.31, 6039.44, 6042.18, 6043.99, 6050.57, 6055.23, 6061.05, 6066.85, 6069.9, 6070.36, 6087.35, 6088.44, 6106.52, 6108.17, 6115.14, 6130.03, 6152.77, 6153.94, 6164.32, 6166.12, 6178.13, 6185.49, 6194.08, 6215.69, 6229.19, 6236.32, 6236.32, 6241.71, 6241.89, 6241.89, 6249.11, 6253.9, 6260.34, 6269.08, 6271.05, 6283.18, 6284.26, 6285.09, 6297.04, 6313.57, 6336.97, 6354.72, 6360.96, 6377.02, 6377.23, 6379.79, 6394, 6412.37, 6415.97, 6423.05, 6430.58, 6432.87, 6456.42, 6464.31, 6475.37, 6476.12, 6476.36, 6476.6, 6482.5, 6487.18, 6488.79, 6490.49, 6516.69, 6528.45, 6542.39, 6551.39, 6552.02, 6552.34, 6562.38, 6568, 6574.69, 6575.49, 6577.4, 6587.48, 6590.13, 6599.37, 6604.62, 6618.51, 6620.85, 6645.47, 6647.53, 6650.43, 6652.4, 6671.89, 6676.94, 6696.94, 6699.04, 6700, 6702.87, 6708.8, 6716.38, 6728.05, 6735.32, 6736.65, 6745.85, 6751.41, 6751.88, 6771.59, 6773.12, 6791.07, 6800.62, 6808.76, 6817.39, 6826.15, 6838.4, 6840.59, 6857.61, 6862.71, 6864.73, 6865.01, 6876.1, 6912.89, 6918.41, 6940.24, 6941.99, 6949.41, 6951.38, 6966.58, 6969.22, 6972.65, 6979.15, 6983.71, 6987.76, 6995.7, 7007.56, 7010.81, 7016.51, 7019.05, 7021.4, 7028.68, 7033.11, 7035.04, 7042.86, 7043.13, 7045.18, 7048.95, 7053.9, 7054.11, 7056.12, 7067.82, 7072.53, 7080.68, 7080.8, 7081.09, 7093.64, 7100.12, 7108.1, 7135.33, 7144.25, 7161.73, 7164.32, 7166.59, 7219.47, 7237.61, 7241.22, 7241.66, 7250.39, 7253.75, 7266.73, 7275.59, 7277.04, 7282.56, 7289.26, 7290.72, 7299.49, 7308.83, 7322.58, 7345.48, 7374.93, 7384.17, 7384.86, 7429.2, 7429.59, 7430.77, 7438.75, 7446.47, 7465.35, 7481.08, 7484.02, 7495.55, 7499.84, 7500.09]), np.array([3210.85, 3225.13, 3231.01, 3231.1, 3238.73, 3242.25, 3242.62, 3250.8, 3252.26, 3274.64, 3287.2, 3304.37, 3304.65, 3307.08, 3319.73, 3331.59, 3334.62, 3338.98, 3342.29, 3343.53, 3353.58, 3373.97, 3375.56, 3381.88, 3385.18, 3386.79, 3389.32, 3409.23, 3412.91, 3421.65, 3433.96, 3438.86, 3447.94, 3455.84, 3458.61, 3464.81, 3465.52, 3471.96, 3473.11, 3477.89, 3483.54, 3483.84, 3503.46, 3507.4, 3514.39, 3524.13, 3526.77, 3533.21, 3561.85, 3572.3, 3579.26, 3583.25, 3599.89, 3601.98, 3604.78, 3610.25, 3634.32, 3634.5, 3644.63, 3649.62, 3661.17, 3661.24, 3679.97, 3688.46, 3693.84, 3701.19, 3721.33, 3729.29, 3733.41, 3750.62, 3760.38, 3777.33, 3788.16, 3788.36, 3800.07, 3802.1, 3812.9, 3828.41, 3836.63, 3844.64, 3846.07, 3858.8, 3864.05, 3869.64, 3878.5, 3880.77, 3882.18, 3884.49, 3888.44, 3890.28, 3897.17, 3904.77, 3905.55, 3909.2, 3923.91, 3943.95, 3943.95, 3948.38, 3951.67, 3958.12, 3964.39, 3970.02, 3976.04, 3985.22, 3993.91, 3996.24, 4006.6, 4019.32, 4022.02, 4047.35, 4054.43, 4056.31, 4060.41, 4077.69, 4085.32, 4087.64, 4112.85, 4118.51, 4119.92, 4125.76, 4136.87, 4139.97, 4142.74, 4145.61, 4186.47, 4205.92, 4222.03, 4237.46, 4238.83, 4245.86, 4248.82, 4263.71, 4302.61, 4316.15, 4339.21, 4342.95, 4344.39, 4352.67, 4356.32, 4360.24, 4389.39, 4396.1, 4398.28, 4404.42, 4424.87, 4429.95, 4431.14, 4434.31, 4464.31, 4497.61, 4500.25, 4501.52, 4523.44, 4529.15, 4532.93, 4535.42, 4540.15, 4546.88, 4548.75, 4553.36, 4554.43, 4554.67, 4555.03, 4557.49, 4563.58, 4565.24, 4566.48, 4568.29, 4589.27, 4590.82, 4597.26, 4603.67, 4612.29, 4618.34, 4618.73, 4626.1, 4634.97, 4656.7, 4656.75, 4662.79, 4665.52, 4666.01, 4666.42, 4667.9, 4670.04, 4686.21, 4689.23, 4695.39, 4726.31, 4746.88, 4748.49, 4754.23, 4755.89, 4759.3, 4765.67, 4770.9, 4776.34, 4777.98, 4785.32, 4796.02, 4802.8, 4838.45, 4848.72, 4858.49, 4869.3, 4869.51, 4870.33, 4885.16, 4936.76, 4945.58, 4949.19, 4968.83, 4972.8, 4974.54, 4974.62, 4998.32, 5003.04, 5020.75, 5021.17, 5042.6, 5042.96, 5049.29, 5051.48, 5072.04, 5077.83, 5085.4, 5092.34, 5098.26, 5124.62, 5129.93, 5133.48, 5139.32, 5145.02, 5165.6, 5170.17, 5187.99, 5199.96, 5201.66, 5205.75, 5207.1, 5220.45, 5239.91, 5240.75, 5253.89, 5254.03, 5259.09, 5272.93, 5273.83, 5276.83, 5286.05, 5290.31, 5293.46, 5297.23, 5298.62, 5319.56, 5322.99, 5325.28, 5350.05, 5350.18, 5353.98, 5393.73, 5415.6, 5427.74, 5440.79, 5442.26, 5452.74, 5460.87, 5470.23, 5475.08, 5476.55, 5483.74, 5491.42, 5498.5, 5502.1, 5503.6, 5504.08, 5521.86, 5525.65, 5527.65, 5540.07, 5546.1, 5551.78, 5555.6, 5564.82, 5565.39, 5569.39, 5584.25, 5587.96, 5595.08, 5612.36, 5613.74, 5617.2, 5649.73, 5683.81, 5685.29, 5687.57, 5702.12, 5705.63, 5717.99, 5720.83, 5740.77, 5748.1, 5751.28, 5753.97, 5756.18, 5765.95, 5767.23, 5770.66, 5771.88, 5787.53, 5790.13, 5793.09, 5830.93, 5835.78, 5836.63, 5872.03, 5902.32, 5937.61, 5939.41, 5965.24, 5966.08, 5966.63, 5967.31, 5972.25, 5976.46, 5987.24, 5987.86, 6014.68, 6019.4, 6028.26, 6029.79, 6034.67, 6040.28, 6053.51, 6056.61, 6058.47, 6062.55, 6066.38, 6068.52, 6091.05, 6094.03, 6098.82, 6104.18, 6108.62, 6121.55, 6166.09, 6210.92, 6215.43, 6228.09, 6230.86, 6231.19, 6235.75, 6238.95, 6256.02, 6258.31, 6271.18, 6288.93, 6296.72, 6302.53, 6303.73, 6315.28, 6344.9, 6350.18, 6350.6, 6371.12, 6371.31, 6387.58, 6403.02, 6403.09, 6410.11, 6444.74, 6446.07, 6454.23, 6454.29, 6454.85, 6467.19, 6467.94, 6470.48, 6470.66, 6483.52, 6488.26, 6526.68, 6530.35, 6535.99, 6557.87, 6567.04, 6583.14, 6585.81, 6608.26, 6622.51, 6637.29, 6643.62, 6644.62, 6647.72, 6648.39, 6655.09, 6660.59, 6664.82, 6667.57, 6683.86, 6699.99, 6718.89, 6721.81, 6724.41, 6724.61, 6726.75, 6730.48, 6735.05, 6760.68, 6761.04, 6762.54, 6768.87, 6769.35, 6794.84, 6795.16, 6799.86, 6801.14, 6807.01, 6812.76, 6815.14, 6818.41, 6818.92, 6832.41, 6836.22, 6838.34, 6846.29, 6851.4, 6855.84, 6855.97, 6857.37, 6864.36, 6868.25, 6873.65, 6875.01, 6888.11, 6900.2, 6903.07, 6921.31, 6954.3, 6972.49, 6974.35, 6983.09, 6998.44, 6998.96, 7014.58, 7014.91, 7056.68, 7057.47, 7064.65, 7079.43, 7087.07, 7088.19, 7094.34, 7098.07, 7099.47, 7106.96, 7128.53, 7130.73, 7140.14, 7163.7, 7167.67, 7182.19, 7190.69, 7199.38, 7210.42, 7231.04, 7233.48, 7236.71, 7242.31, 7245.94, 7246.48, 7257.36, 7259.05, 7285.26, 7289.11, 7301.48, 7307.94, 7315.45, 7332.73, 7333.58, 7338.45, 7340.57, 7349.99, 7352.78, 7356.59, 7361.42, 7361.84, 7367.73, 7371.86, 7380.11, 7406.01, 7410.79, 7414.2, 7417.24, 7420.11, 7432.98, 7447.98, 7450.58, 7450.94, 7451.49, 7454.12, 7480.78, 7488.27, 7501.95]), np.array([3199.85, 3208.85, 3211.39, 3246.42, 3246.48, 3264.4, 3265.76, 3270.07, 3277.3, 3284.33, 3288.33, 3302.02, 3304.44, 3305.8, 3315.3, 3342.97, 3352.9, 3360.65, 3363.53, 3370.78, 3375.01, 3376.64, 3397.91, 3413.43, 3419.18, 3428.99, 3437.58, 3460.82, 3461.57, 3471.02, 3474.95, 3474.95, 3477.6, 3478.66, 3536.41, 3560.22, 3561.16, 3564.96, 3566.44, 3568.37, 3587.59, 3604.2, 3616.14, 3616.18, 3628.65, 3630.69, 3636.91, 3640.95, 3648.81, 3650.92, 3662.28, 3667.81, 3671.18, 3672.8, 3673.15, 3681.18, 3709.02, 3710.63, 3724.41, 3729.7, 3746.52, 3760.66, 3769.39, 3770.44, 3773.48, 3780.39, 3781.47, 3792.08, 3802.2, 3804.72, 3834.88, 3840.44, 3852.54, 3856.87, 3888.06, 3889.99, 3896.64, 3897.9, 3901.05, 3922.66, 3927.18, 3938.26, 3940.53, 3943.83, 3949.69, 3949.92, 3973.89, 3976.39, 3991.35, 4017.32, 4020.14, 4043.32, 4061.19, 4065.79, 4066.22, 4076.06, 4082.08, 4092.86, 4097.28, 4105.14, 4108.61, 4123.28, 4132.36, 4140.87, 4144.46, 4150.61, 4153.33, 4154.67, 4156.21, 4174.46, 4180.63, 4186.19, 4199.02, 4199.55, 4202.94, 4223.12, 4230.97, 4251.31, 4254.48, 4265.72, 4289.79, 4293.98, 4295.51, 4324.6, 4336.33, 4339.39, 4340.84, 4342.06, 4353.92, 4368.24, 4369.87, 4376.27, 4388.78, 4398.44, 4400.29, 4403.49, 4404.08, 4408.94, 4444.17, 4446.96, 4453.25, 4463.9, 4476.19, 4504.07, 4509.46, 4530.41, 4533.67, 4533.75, 4542.83, 4551.08, 4553.79, 4557.33, 4560.45, 4567.43, 4568.54, 4568.59, 4573.29, 4579.92, 4590.16, 4595.44, 4605.22, 4605.89, 4629.43, 4639.78, 4643.04, 4653.45, 4672.18, 4680.99, 4684.48, 4689.16, 4701, 4717.34, 4720.64, 4731.11, 4737.3, 4743.6, 4757.45, 4771.7, 4781.02, 4783.06, 4793.88, 4801.43, 4805.42, 4809.67, 4830.71, 4832.25, 4843.91, 4848.49, 4855.92, 4859.98, 4870.3, 4878.15, 4882.96, 4896.14, 4897.08, 4905.52, 4907.36, 4907.68, 4924.28, 4953.08, 4954.4, 4966.89, 4974.84, 4985.46, 4988.43, 4990.77, 5006.5, 5011.67, 5038.25, 5040.75, 5050.22, 5057.86, 5072.95, 5074.46, 5092.75, 5112.73, 5113.7, 5114.06, 5114.53, 5120.3, 5128.13, 5134.03, 5138.53, 5149.03, 5152.06, 5161.49, 5184.54, 5185.67, 5187.3, 5188.82, 5217.01, 5222.79, 5223.21, 5225.48, 5232.49, 5233.99, 5264.16, 5264.82, 5267.96, 5269.93, 5272.16, 5272.99, 5279.89, 5280.18, 5281.16, 5287.24, 5301.59, 5316.76, 5321.94, 5327.7, 5337.8, 5346.71, 5362.55, 5367.16, 5374.31, 5379.34, 5394.66, 5401.13, 5403.64, 5405.06, 5409.62, 5409.71, 5435.06, 5435.38, 5445.6, 5449.22, 5452.34, 5453.17, 5465.21, 5471.23, 5473.04, 5473.94, 5476.06, 5486.32, 5489.92, 5493.14, 5506.53, 5507.38, 5517.04, 5522.6, 5536.82, 5538.48, 5544.04, 5551.58, 5567.62, 5596.09, 5598.28, 5598.95, 5606.94, 5607.23, 5608.23, 5611.39, 5615.08, 5619.88, 5621.87, 5622.94, 5636.85, 5638.09, 5640.81, 5647.14, 5652.78, 5654.16, 5656.88, 5664.13, 5679.65, 5688.14, 5692.39, 5700.86, 5708.83, 5715.48, 5742.93, 5743.45, 5750.75, 5764.94, 5766.93, 5784.5, 5785.58, 5796.49, 5807.19, 5809.45, 5817.39, 5834.62, 5857.29, 5862.55, 5883.77, 5893.74, 5894.07, 5894.48, 5895.94, 5910.8, 5915.08, 5916.22, 5923.72, 5926.62, 5929.44, 5939.85, 5941.43, 5963.58, 5975.82, 5980.71, 5984.06, 5984.85, 5989.73, 5992.07, 5993.62, 5993.76, 6001.2, 6005.24, 6006.68, 6010.51, 6015.71, 6018.46, 6024.74, 6031.1, 6034.81, 6043.07, 6044.57, 6051.6, 6051.97, 6054.53, 6069.51, 6097.32, 6105.4, 6112.42, 6114.06, 6116.16, 6121.7, 6129.5, 6130.64, 6130.76, 6136.85, 6142.11, 6143.44, 6151.86, 6167.08, 6171.26, 6172.88, 6177.17, 6178.86, 6183.76, 6190.23, 6209.06, 6209.45, 6212.7, 6216.96, 6217.29, 6231.57, 6253.84, 6255.23, 6263.27, 6263.37, 6272.55, 6284.62, 6289.14, 6294.66, 6297.75, 6299.82, 6311.29, 6315.44, 6324.69, 6327.99, 6334.53, 6338.92, 6347.25, 6348.58, 6348.71, 6349.43, 6350.45, 6359.99, 6361.66, 6363.49, 6367.95, 6373.09, 6376.44, 6379.17, 6414.36, 6416.48, 6416.95, 6437.02, 6441.39, 6457.47, 6473.78, 6483.51, 6490.13, 6492.93, 6494.79, 6495.88, 6500.35, 6500.48, 6516.49, 6516.77, 6516.86, 6518.76, 6541.46, 6549.56, 6564.6, 6566.3, 6589.52, 6596.49, 6599.14, 6601.3, 6608.2, 6612.08, 6614.72, 6617.85, 6633.38, 6637.47, 6642.49, 6646.44, 6660.75, 6661.15, 6676.12, 6688.61, 6695.88, 6710.2, 6713.32, 6723.16, 6726.89, 6739.95, 6740.59, 6744.54, 6746.32, 6748.35, 6751.77, 6757.2, 6759.24, 6781.45, 6781.83, 6791.35, 6799.68, 6800.73, 6811.33, 6814.78, 6829.44, 6840.78, 6846.23, 6847.97, 6852.17, 6854.34, 6856.33, 6861.62, 6891.14, 6901.05, 6910.31, 6917.12, 6936.92, 6944.94, 6945.85, 6949.77, 6952, 6954.44, 6957.19, 6966.75, 6973.98, 6983.25, 6985.75, 6986.75, 7045.58, 7051.5, 7052.06, 7059.65, 7063.56, 7066.41, 7074.29, 7076.6, 7087.17, 7095.7, 7102.01, 7103.63, 7106.9, 7113.16, 7123.43, 7137.14, 7152.99, 7154.92, 7167.25, 7173.05, 7188.04, 7195.58, 7197.93, 7200.19, 7208.2, 7224.1, 7228.81, 7232.65, 7233.07, 7236.07, 7251.42, 7256.45, 7258.83, 7260.04, 7266.64, 7271.71, 7284.35, 7293.21, 7299.02, 7305.28, 7323.22, 7327.58, 7328.2, 7332.28, 7342.75, 7350.44, 7358.8, 7378.77, 7385.15, 7406.44, 7407.05, 7407.46, 7426.18, 7429.72, 7442.95, 7461.15, 7463.41, 7468.65, 7481.1, 7481.63, 7495.65, 7505.57]), np.array([3202.09, 3205.57, 3211.11, 3218.36, 3248.75, 3267.63, 3270.91, 3272.1, 3316.08, 3316.73, 3319.92, 3340.5, 3358.54, 3367.95, 3375.92, 3378.26, 3381.24, 3396.62, 3412.48, 3421.43, 3438.83, 3443.67, 3445.02, 3446.44, 3446.93, 3452.72, 3455.76, 3476.08, 3490.38, 3494.66, 3496.53, 3500.55, 3505.98, 3510.02, 3514.16, 3519.99, 3522.65, 3541.56, 3543.28, 3545.09, 3555.33, 3557.31, 3561.04, 3585.5, 3597.19, 3598.05, 3600.22, 3605.02, 3609.9, 3614.28, 3618.64, 3623.74, 3624.05, 3627.37, 3629.82, 3635.23, 3640.56, 3685.74, 3697.03, 3717.1, 3726.18, 3731.16, 3736.13, 3738.77, 3753.24, 3784.56, 3792.35, 3794.69, 3805.42, 3810.77, 3813.39, 3826.06, 3827.96, 3828.99, 3829.35, 3832.81, 3836.04, 3837.28, 3850.85, 3860.62, 3863.45, 3867.63, 3881.99, 3896.64, 3900.85, 3923.67, 3926.4, 3931.29, 3953.17, 3960.46, 3962.37, 3962.98, 3967.46, 3976.55, 3979.66, 3981.22, 3998.03, 4002.69, 4010.51, 4022.12, 4024.99, 4036.81, 4052.33, 4055.05, 4055.41, 4057.32, 4057.86, 4058.22, 4059.72, 4100.3, 4102.53, 4109.21, 4110.87, 4120.82, 4122, 4126.49, 4138.3, 4147.62, 4152.29, 4155.11, 4161.35, 4162.83, 4169.56, 4196.36, 4196.99, 4197.86, 4205.61, 4209.78, 4214.86, 4234.12, 4238.55, 4240.4, 4255.22, 4256.54, 4265.93, 4284.03, 4284.68, 4298.82, 4322.14, 4326.16, 4352.17, 4357.16, 4366.47, 4375.01, 4375.83, 4381.76, 4382.88, 4395.75, 4396.93, 4405.6, 4409.7, 4412.12, 4412.34, 4419.57, 4425.26, 4441.47, 4447.08, 4452.41, 4462.24, 4462.72, 4469.96, 4473.63, 4479.33, 4492.76, 4502.66, 4512.29, 4521.06, 4521.8, 4522.83, 4543.15, 4549.62, 4560.72, 4568.95, 4574.44, 4591.26, 4594.93, 4597.15, 4607.32, 4609.74, 4612.43, 4612.91, 4640.22, 4645.12, 4645.41, 4645.5, 4661.26, 4663.95, 4664.36, 4673.97, 4683.5, 4695.91, 4699.8, 4700.6, 4702.94, 4712.47, 4715.68, 4717.12, 4721.22, 4723.88, 4723.9, 4727.01, 4734.95, 4735.1, 4739.42, 4751.73, 4754.73, 4778.23, 4778.46, 4778.99, 4794.12, 4833.39, 4834.84, 4841.14, 4861.24, 4869.61, 4874.43, 4884.42, 4896.19, 4897.8, 4944.84, 4949.07, 4955.44, 4963.45, 4966.43, 4973.1, 4980.9, 4985.01, 4986.52, 5007.49, 5018.46, 5029.78, 5031.27, 5031.55, 5033.8, 5035.92, 5039.94, 5042.8, 5055.08, 5061.68, 5075.46, 5079.66, 5083.53, 5086.14, 5115.31, 5116.2, 5149.89, 5158.39, 5165, 5167.84, 5172.08, 5175.57, 5180.18, 5183.63, 5190.56, 5195.64, 5197.85, 5200.56, 5210.21, 5210.49, 5211.48, 5240.29, 5243.54, 5253.24, 5255.43, 5258.57, 5264.18, 5277.62, 5280.48, 5283.19, 5284.85, 5292.92, 5312.59, 5314.08, 5314.86, 5319.19, 5321.76, 5322.28, 5337.57, 5354.5, 5375.18, 5378.23, 5381.12, 5386.51, 5408.29, 5412.67, 5430.06, 5438.38, 5441.96, 5454.98, 5457.64, 5472.01, 5479.42, 5502.2, 5507.76, 5510.14, 5511.17, 5513.25, 5522.58, 5531.55, 5543.05, 5554.98, 5555.55, 5558.35, 5560.97, 5562.96, 5596.63, 5607.5, 5635.44, 5636.13, 5649.79, 5651.26, 5655.06, 5655.37, 5662.53, 5665.14, 5671.28, 5692.9, 5704.61, 5708.09, 5717.61, 5726.13, 5726.26, 5730.01, 5732.34, 5745.77, 5747.01, 5751.79, 5754.46, 5756.56, 5764.44, 5769.13, 5774.72, 5786.12, 5787.5, 5789.37, 5791.08, 5811.43, 5816.72, 5820.04, 5822.95, 5837.78, 5844.95, 5845.69, 5846.41, 5852.6, 5855.6, 5857.82, 5858.48, 5858.56, 5859.26, 5861.9, 5863.8, 5865.34, 5866.98, 5888.64, 5901.67, 5915.22, 5927.53, 5939.61, 5941.83, 5952.68, 5963.75, 5983.61, 6003.1, 6008.11, 6024.2, 6024.65, 6028.47, 6047.12, 6052.5, 6055.37, 6068.62, 6069.45, 6071.62, 6071.78, 6085.03, 6092.68, 6097.52, 6133.72, 6137.54, 6150.27, 6156.89, 6165.62, 6167.8, 6173.03, 6178.59, 6188.34, 6194.52, 6197.67, 6200.55, 6211.74, 6212.3, 6239.53, 6244.08, 6246.88, 6252.05, 6253.72, 6254.43, 6265.65, 6269.31, 6272.62, 6273.1, 6288.1, 6292.35, 6298.78, 6300.48, 6300.6, 6301.67, 6306.87, 6310.14, 6322.82, 6323.25, 6329.41, 6343.79, 6352.27, 6359.55, 6364.01, 6366.96, 6387.79, 6392, 6429.84, 6435.11, 6456.76, 6462.06, 6489.54, 6503.5, 6509.49, 6542.23, 6557.38, 6560.08, 6564.82, 6584.01, 6588.03, 6589.11, 6598.38, 6607.06, 6614.94, 6618.82, 6628.82, 6638.52, 6639, 6641.26, 6645.67, 6654.64, 6655.92, 6659.84, 6660.37, 6669.44, 6682.09, 6694.73, 6701.04, 6713.25, 6743.82, 6746.33, 6770.06, 6772.81, 6775.7, 6786.19, 6804.23, 6805.23, 6810.11, 6819.45, 6823.51, 6835.69, 6841.56, 6842.19, 6843.83, 6845.1, 6853.09, 6856.86, 6857.68, 6877.06, 6881.29, 6883.71, 6890.16, 6924.53, 6932.96, 6933.68, 6935.84, 6936.79, 6940.43, 6968.99, 6973.82, 6995.75, 7004.89, 7012.05, 7027.71, 7029.3, 7033.26, 7038.74, 7043.89, 7058.81, 7062.11, 7081.39, 7126.43, 7127.16, 7129.22, 7132.59, 7134.15, 7143.95, 7153.65, 7162.24, 7162.89, 7163.83, 7170.29, 7175.25, 7182.08, 7186.98, 7192.98, 7194.91, 7195.87, 7209.09, 7217.95, 7221.62, 7225.18, 7230.47, 7239.18, 7241.54, 7248.57, 7250.34, 7282.37, 7285.43, 7286.1, 7288.99, 7291.57, 7315.13, 7327.98, 7341.74, 7343.83, 7347.27, 7359.61, 7360.42, 7367.55, 7373.52, 7375.7, 7385.26, 7389.31, 7398.37, 7417.22, 7420.38, 7422.18, 7428.43, 7433.2, 7439.47, 7447.54, 7468.73, 7471.9, 7480.43, 7490.99, 7491.34, 7495.46, 7510.16]), np.array([3201.42, 3240.23, 3304.45, 3313.93, 3541.3, 3572.54, 3660.99, 3719.92, 3792.88, 3801.63, 3896.1, 4273.95, 4289.34, 4323.04, 4347.41, 4410.98, 4514.34, 4552.42, 4720.55, 4872.21, 4884.95, 4898.71, 5025.9, 5389.71, 5442.38, 5445.14, 5475.28, 5536.65, 5544.08, 5551.62, 5556.88, 5577.18, 5653.4, 5709.3, 5763.63, 5852.35, 5862.59, 5958.59, 5967.39, 6099.12, 6181.72, 6198.62, 6218.47, 6292.2, 6428.5, 6532.82, 6610.61, 6649.27, 6654.22, 6668.97, 6714.78, 6741.11, 6922.77, 7001.57, 7016.64, 7101.01, 7131.08, 7164.46, 7269.78, 7292.92, 7350.84, 7448.98, 7461.03, 7624.72]), np.array([3327.69, 3347.47, 3363.43, 3396.12, 3450.3, 3536.57, 3576.08, 3600.72, 3697.63, 3923.07, 3995.82, 4044.19, 4121.19, 4157.02, 4172.17, 4337.21, 4368.82, 4569.13, 4640.06, 4717.94, 4757.08, 4793.51, 4868.19, 4883.36, 4982.94, 5036.97, 5063.5, 5264.38, 5416.88, 5456.74, 5571.87, 5686.35, 5713.89, 5755.29, 5805.03, 5816.16, 5834.14, 5876.14, 5912.57, 5924.46, 5981.16, 5999.19, 6066.04, 6090.43, 6091.9, 6182.55, 6198.38, 6201.92, 6221.91, 6312.47, 6351.72, 6397.37, 6409.26, 6466.26, 6723.4, 6750.86, 6756.48, 6849.13, 6877.9, 6884.26, 7000.73, 7244.34, 7396.5, 7438.6, 7588.67]), np.array([3300.24, 3381.81, 3525.52, 3526.27, 3632.52, 3726.18, 3747.43, 3778.01, 3958.59, 3981.57, 4085.95, 4171.89, 4256.03, 4279.8, 4499.82, 4526.64, 4544.01, 4560.08, 4663.82, 4836.15, 4885.28, 5018.88, 5019.99, 5045.59, 5299.03, 5329.73, 5336.03, 5411.82, 5628.38, 5646.6, 5813.17, 5850.61, 5938.28, 5940.32, 6051.03, 6125.88, 6226, 6253.53, 6440.89, 6447.29, 6477.5, 6500.8, 6512.83, 6614.31, 6712.31, 6916.76, 7001.26, 7136.88, 7153.13, 7213.62, 7314.59, 7346.85, 7452.71, 7760.08]))
    RList = [[0.00138853, 2.19322e-23], [0.00126339, 4.905e-10], [0.00150304, 1.11925e-24], [0.00171945, 8.93607e-35], [0.00105208, 1.17261e-13], [5.5671e-06, 1.75078e-229], [0.000306829, 2.08171e-07], [0.000344971, 6.91911e-13], [0.00130844, 0.000105671]]
    dendVclamp = np.array([1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000])

    #vecTuple = (np.array([0, run_t+1]),np.array([100,run_t+1]),np.array([200,run_t+1]),np.array([300,run_t+1]),np.array([400,run_t+1]),np.array([500,run_t+1]),np.array([600,run_t+1]),np.array([800,run_t+1]),np.array([1000,run_t+1]))
    #vecTuple = (np.array([0,20,40, run_t+1]),np.array([100,120,140,run_t+1]),np.array([200,220,240,run_t+1]),np.array([300,320,340,run_t+1]),np.array([400,420,440,run_t+1]),np.array([500,520,540,run_t+1]),np.array([600,620,640,run_t+1]),np.array([800,run_t+1]),np.array([1000,run_t+1]))
    #vecTuple0 = [np.array([run_t+1]),np.array([run_t+1]),np.array([run_t+1]),np.array([run_t+1]),np.array([run_t+1]),np.array([run_t+1]),np.array([run_t+1]),np.array([run_t+1]),np.array([run_t+1])]
    #vecTuple = (np.array([0,330,run_t+1]),np.array([run_t+1]),np.array([run_t+1]),np.array([run_t+1]),np.array([run_t+1]),np.array([run_t+1]),np.array([run_t+1]),np.array([run_t+1]),np.array([0,660,run_t+1]))
    #RList = [[0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0]]
    #dendVclamp = np.array([1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000])
    v_pre = -59.9757
    t0 = 3198.7
    tref = 13
    vThres = -60.0
    vBack = -61.0
    printR = False
    getDendV = False 
    #RList = np.zeros((n,2))
    oneGo = False 
    clampDend = False
    monitorDend = False 
    trans = 16 
    pas = False 
    #trans = 100
    #clampEffects(cell, v0, trans, trans+run_t, t0, tstep, clampDend)
    #trans = 200
    #clampEffects(cel, v0, trans, trans+run_t, t0, tstep, clampDend)
    #trans = 250
    #clampEffects(cell, v0, trans, trans+run_t, t0, tstep, clampDend)
    count = 0
    #for v_pre in range(-74,-58,4):
    #    for i in range(0,locE.size):
    #        count = count + 1
    #        vecTuple = list(vecTuple0)
    #        vecTuple[i] = np.array([0,3,run_t+1]);
    v, fired, tsp, ntrans = proceed(cell, v_pre, synList, RList, vecStimList, vecTuple, n, trans, trans+run_t, vBack, tref, vThres, oneGo, t0, tstep, loc, pos, dendVclamp, printR, alphaR, getDendV, monitorDend, pas)
    pyplot.savefig('slice.png',format='png',bbox_inches='tight',dpi=900)
    print count
