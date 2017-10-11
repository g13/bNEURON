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

def run(cell, v0, vBack, tref, vThres, synList, RList, n, trans, oneGo, t0, printR = False, alphaR = True, loc=np.array([]), pos=np.array([])):
    f = open('pylog','a')
    nc = 0
    cell.init()
    print " cell initiated"
    print >>f, " cell initiated"
    steps = 0
    if oneGo:
        tsps = np.empty(int((h.tstop-t0)/h.dt))
    else:
        tsps = np.empty(1)
    tsp = 0
    #plusOne = 0
    vold = v0
    fired = 0 
    h.t = t0
    print  "t0 = ", h.t
    print >>f,  "t0 = ", h.t
    while round(h.t/h.dt) < round((t0+trans)/h.dt):
        h.fadvance()
        steps = steps + 1
        #print  cell.soma(0.0).v
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
    print  trans, "ms trans complete, used ", steps , ' steps, t+trans = ', h.t, 'v = ', '%7.5f.' % cell.soma(0.0).v
    print >>f,  trans, "ms trans complete, used ", steps , ' steps, t+trans = ', h.t, 'v = ', '%7.5f.' % cell.soma(0.0).v
    if len(loc)>0:
        for i in xrange(n): 
            print " dend[", loc[i], "]" ,"(",pos[i],").v = ", cell.dend[loc[i]](pos[i]).v
            print >> f, " dend[", loc[i], "]" ,"(",pos[i],").v = ", cell.dend[loc[i]](pos[i]).v
    checkList = [];
    #while ( ((cell.soma(0.0).v > vBack and not fired ) or ( (h.t < trans + tsp + tref or cell.soma(0.0).v > vThres-2)  and fired ) or oneGo) and h.t < h.tstop):
    while ( ((cell.soma(0.0).v > vBack and not fired ) or ( (h.t < trans + tsp + tref)  and fired ) or oneGo) and h.t < h.tstop):
        h.fadvance()
        if cell.soma(0.0).v > vThres+15 and cell.soma(0.0).v < vold and not fired:
            tsp = h.t-trans
            tsps[nc] = tsp
            fired = 1
            nc = nc + 1
            print  tsp, "fired expecting finish ", tsp + tref
            print >>f,  tsp, "fired expecting finish ", tsp + tref
        vold = cell.soma(0.0).v
        if oneGo:
            if len(checkList) > 0 and len(loc) > 0:
                for checkTime in checkList:
                    if abs(h.t-checkTime-trans) < 1e-10:
                        print " check at ", checkTime
                        for i in xrange(n): 
                            print "     dend[", loc[i], "]" ,"(",pos[i],").v = ", cell.dend[loc[i]](pos[i]).v
                            print >> f, "   dend[", loc[i], "]" ,"(",pos[i],").v = ", cell.dend[loc[i]](pos[i]).v
            if h.t > tsp + trans + tref:
                fired = 0
            #progress = h.t/h.tstop*100 % 10
            #if progress > oldProgress + 1:
            #    print  progress
            #    print >>f,  progress
            #    oldProgress = progress
    print  "stopping with v ", cell.soma(0.0).v, " t ", h.t
    print >>f,  "stopping with v ", cell.soma(0.0).v, " t ", h.t
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

def proceed(cell, v0, synList, RList, vecStimList, spikeTrain, n, trans, tend, vBack, tref, vThres, oneGo, t0, tstep, loc=np.array([]), pos=np.array([]), dendVclamp=np.array([]), printR = False, alphaR = True, getDendV = False, monitorDend = False):
    f = open('pylog','a')
    ferr = open('pyerr','w')
    #sys.stderr = ferr
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
    v.record(cell.soma(0.0)._ref_v)
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
            print "played ", i
            print >>f, "played ", i
            if alphaR:
                synList[i].deltat = h.dt
    print   "ready to run"
    print >>f, "ready to run"
    fired, tsp = run(cell, v0, vBack, tref, vThres, synList, RList, n, trans, oneGo, t0, printR, alphaR, loc, pos)
    #print   "i ran and i killed ", fired, " bedbug(s)"
    #print >>f,   "i ran and i killed ", fired, " bedbug(s)"
    v1 = v.as_numpy()
    if oneGo:
        ntotal = int(round((tend-t0)/tstep)+1)
        if __name__ == '__main__':
            t1 = t.as_numpy()
        if v1.size!=ntotal:
            print " steps ", v1.size, ", ntotal ", ntotal
            print " remove extra step"
            print >>f, " remove extra step"
            v1 = v1[:-1]
            if __name__ == '__main__':
                t1 = t1[:-1]
    else:
        ntotal = int(round((h.t-t0)/tstep)+1)
    print   " is this nan ", v1[-1], v1[0]
    print >>f,   " is this nan ", v1[-1], v1[0]
    print   " v size ", v1.size
    print >>f,   " v size ", v1.size
    ntrans = int(round(trans/tstep))
    #if v1[ntrans] == v1[ntrans-1]:
    #    ntrans = ntrans+1
    #assert(v1[ntrans] != v1[ntrans-1])
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
        pyplot.draw()
        #pyplot.savefig('slice.png',format='png',bbox_inches='tight',dpi=900)
        
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
    ferr = open('pyerr','w')
    #sys.stderr = ferr
    h.tstop = tend
    print loc
    print >>f, loc
    print pos
    print >>f, pos

    h.dt = tstep
    v = h.Vector()
    t = h.Vector()
    v.record(cell.soma(0.0)._ref_v)
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
    v.record(cell.soma(0.0)._ref_v)
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
    h.dt = tstep
    v = h.Vector()
    if pos != -1.0:
        dendv = h.Vector()
        dendv.record(cell.dend[sel[0]](pos)._ref_v)
    t = h.Vector()
    v.record(cell.soma(0.0)._ref_v)
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
    h.dt = tstep
    v = h.Vector()
    if pos != -1.0:
        dendv = h.Vector()
        dendv.record(cell.dend[sel[0]](pos)._ref_v)
    t = h.Vector()
    v.record(cell.soma(0.0)._ref_v)
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
    print trans, "ms trans complete, used ", steps , ' steps, t+trans = ', h.t, 'v = ', '%7.5f.' % cell.soma(0.0).v
    print' distant apical dend v = ', '%7.5f.' % cell.dend[121](0.0).v
    print' distant basal dend v = ', '%7.5f.' % cell.dend[153](0.0).v
    print >>f,  trans, "ms trans complete, used ", steps , ' steps, t+trans = ', h.t, 'v = ', '%7.5f.' % cell.soma(0.0).v
    print >>f, ' distant apical dend v = ', '%7.5f.' % cell.dend[121](0.0).v
    print >>f, ' distant basal dend v = ', '%7.5f.' % cell.dend[153](0.0).v
    while (h.t < h.tstop):
        h.fadvance()
        steps = steps + 1
    print "srun ended in ", steps , ' steps, t+trans = ', h.t, 'v = ', '%7.5f.' % cell.soma(0.0).v
    print ' distant apical dend v = ', '%7.5f.' % cell.dend[121](0.0).v
    print ' distant basal dend v = ', '%7.5f.' % cell.dend[153](0.0).v
    print >>f,  " srun ended in ", steps , ' steps, t+trans = ', h.t, 'v = ', '%7.5f.' % cell.soma(0.0).v
    print >>f, ' distant apical dend v = ', '%7.5f.' % cell.dend[121](0.0).v
    print >>f, ' distant basal dend v = ', '%7.5f.' % cell.dend[153](0.0).v
    f.close()

def clampEffects(cell, v0, trans, tend, t0, tstep, clampDend):
    f = open('pylog','a')
    #ferr = open('pyerr','w')
    #sys.stderr = ferr
    h.tstop = tend
    print   "     ", t0, ' with ', trans, ' trans ', ' to ', tend
    print >>f,   "     ", t0, ' with ', trans, ' trans ', ' to ', tend 

    h.dt = tstep
    v = h.Vector()
    t = h.Vector()
    v.record(cell.soma(0.0)._ref_v)
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
    run_t = 200 
    #dtv = [10,25,50,100,170,240]
    run_nt = int(round(run_t/tstep))
    tol_t = min([300,run_t])
    tol_nt = int(round(tol_t/tstep))+1
    seed = 13
    np.random.seed(seed)
    #locE = np.array([36, 53, 74, 79, 83, 90, 97, 101, 117, 125, 132, 174],dtype='int')
    #locE = np.array([36, 74, 83, 97, 117, 132],dtype='int')
    #locE = np.array([36, 73, 97, 117, 146, 180],dtype='int')
    locE = np.array([32, 52, 66, 78, 98, 136],dtype='int')
    #locE = np.array([74],dtype='int')
    #gE = 1e-4 + np.random.random_sample(locE.size) * (1e-4-1e-4) + 2e-2 + 1e-3*3
    #gE = 3.2e-2 * (1 + np.random.random_sample(locE.size) * 0.2 - 0.1)
    gE = 1e-5 + np.random.random_sample(locE.size) * (1e-2-1e-5)
    posE = np.random.random_sample(locE.size)
    locI = np.array([2, 14, 28],dtype='int')
    #locI = np.array([7, 28, 137],dtype='int')
    #locI = np.array([28],dtype='int')
    gI = (1e-5 + np.random.random_sample(locI.size) * (2e-3-1e-5)) * (-1.0/2.0)
    #gI = (1e-5 + np.random.random_sample(locI.size) * (1e-4-1e-5)) * (-0) - 8e-3 - 1e-3*3
    #gI = -1e-3 * (1 + np.random.random_sample(locI.size) * 0.2 - 0.1)
    posI = np.random.random_sample(locI.size)
    pos = np.concatenate((posE, posI))
    loc = np.concatenate((locE, locI))
    gList = np.concatenate((gE, gI))
    #loc = np.array([28])
    #gList = np.array([-1e-2])    
    #pos = np.array([0.5])
    #loc = np.array([137,101])
    #gList = np.array([-6.3e-04,1.16e-03])    
    #pos = np.array([0.08,0.01])
    n = loc.size
    alphaR = True 
    cell, vecStimList, synList = prepCell(gList, loc, pos, n, -70, alphaR)

    sel = range(n)

    #vecTuple = (np.array([84.7, 85.4, 108.8, 113, 117.9, 120, 134.7, 151, 181.5, 187.7, 197.8, 210.5]), np.array([56, 57.1, 70.4, 87.4, 97.9, 120.9, 124.6, 142, 162.7, 175.9, 186.9, 189.8, 210.4]), np.array([72.7, 75.7, 85.7, 86.7, 109.3, 112.4, 114.8, 122.5, 134.4, 141.7, 150.4, 153.5, 164.5, 167.6, 174.8, 179.4, 192.6, 193.1, 203.9]), np.array([71.3, 79, 89.6, 108.2, 113.8, 125, 147.2, 149.8, 161.4, 162.9, 173.6, 178.1, 197.9, 210.5]), np.array([73.7, 78.4, 83.5, 85.8, 87.2, 109.1, 111.2, 111.6, 129.9, 139.5, 141.7, 142.9, 143.6, 155.5, 175.4, 184.8, 193.4, 202.5]), np.array([90.2, 102.1, 123.6, 125.2, 132.8, 139.6, 148.2, 148.4, 154.4, 162.2, 163.9, 174.9, 175.8, 187, 224.2]), np.array([97.1, 178.3, 208.7]), np.array([60.3, 63.1, 102.6, 103.4, 149.2, 200.1]), np.array([188.2, 252.6]))
    #RList = [[0.00198759, 8.09972e-13], [0.000117574, 3.95497e-51], [0.00158127, 4.12393e-14], [0.00143576, 2.9786e-05], [0.000711666, 1.24271e-20], [0.000134815, 1.92519e-57], [2.32916e-06, 3.49001e-09], [1.80668e-05, 1.5604e-17], [4.35576e-06, 1.38783e-16]]
    #dendVclamp = np.array([-69.5971, -66.8555, -66.0087, -67.1818, -65.2955, -67.8925, 1000, 1000, -69.4652])
    #vecTuple = (np.array([15.6, 37.7, 40.6, 44.5, 49.7, 50.9, 84.7, 85.4, 108.8, 113, 117.9, 120, 134.7, 151, 181.5, 187.7, 197.8, 210.5]), np.array([20.1, 36, 42.3, 56, 57.1, 70.4, 87.4, 97.9, 120.9, 124.6, 142, 162.7, 175.9, 186.9, 189.8, 210.4]), np.array([27, 38.8, 49.8, 50.6, 72.7, 75.7, 85.7, 86.7, 109.3, 112.4, 114.8, 122.5, 134.4, 141.7, 150.4, 153.5, 164.5, 167.6, 174.8, 179.4, 192.6, 193.1, 203.9]), np.array([43.8, 52.6, 71.3, 79, 89.6, 108.2, 113.8, 125, 147.2, 149.8, 161.4, 162.9, 173.6, 178.1, 197.9, 210.5]), np.array([0.7, 7, 49.1, 73.7, 78.4, 83.5, 85.8, 87.2, 109.1, 111.2, 111.6, 129.9, 139.5, 141.7, 142.9, 143.6, 155.5, 175.4, 184.8, 193.4, 202.5]), np.array([1.8, 6, 8.1, 11.8, 12.6, 40.8, 90.2, 102.1, 123.6, 125.2, 132.8, 139.6, 148.2, 148.4, 154.4, 162.2, 163.9, 174.9, 175.8, 187, 224.2]), np.array([12.7, 22.5, 43.7, 97.1, 178.3, 208.7]), np.array([18.5, 22.1, 60.3, 63.1, 102.6, 103.4, 149.2, 200.1]), np.array([25.1, 188.2, 252.6]))
    RList = [[0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0]]
    vecTuple = ()
    #RList = [[0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0]]
    dendVclamp = np.array([1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000])
    #v_pre = -60.9841
    #t0 = 53.4
    v_pre = -66
    #v_pre = -74
    t0 = 0.0
    tref = 3
    vThres = -62.5
    vBack = -63.5
    printR = False
    getDendV = True 
    #RList = np.zeros((n,2))
    oneGo = 1
    clampDend = False 
    monitorDend = True 
    v0 = -74
    trans = 100
    clampEffects(cell, v0, trans, trans+run_t, t0, tstep, clampDend)
    trans = 200
    clampEffects(cell, v0, trans, trans+run_t, t0, tstep, clampDend)
    trans = 250
    clampEffects(cell, v0, trans, trans+run_t, t0, tstep, clampDend)
    #v, fired, tsp, ntrans, _ = proceed(cell, v_pre, synList, RList, vecStimList, vecTuple, n, trans, trans+run_t, vBack, tref, vThres, oneGo, t0, tstep, loc, pos, dendVclamp, printR, alphaR, getDendV, monitorDend)
    #v, fired, tsp, ntrans, _ = proceed(cell, v_pre, synList, RList, vecStimList, vecTuple, n, trans, trans+run_t, vBack, tref, vThres, oneGo, t0, tstep, loc, pos, dendVclamp, printR, alphaR, getDendV, monitorDend)
    pyplot.show()
    #print v
    #print fired
    #print tsp 
    #print ntrans
