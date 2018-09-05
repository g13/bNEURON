from neuron import h, gui
from n128 import *
import re
import numpy as np

seed = 231276
mks = 2
bs = 8
np.random.seed(seed)
#locE = np.random.randint(75,134,6)
locE = np.array([80, 89, 91, 97, 121, 126])
#locI = np.random.randint(1,75,3)
locI = np.array([25,31,43])
posE = np.random.random_sample(locE.size)
posI = np.random.random_sample(locI.size)

def p(type, x, y, keystate):
    if type == 2:
        d = vc.nearest(x, y)
        arc = vc.push_selected()
        if arc >= 0:
            cs = h.cas()
            #vc.color(2)
            sn = str(cs)
            match = re.search('\.(dend|soma)\d*\Z',sn)
            if not match == None:
                name = match.group(0)[5:]
                print '%g from %s' % (d, name)
                print "at x=",x,", y=",y
                vc.label(x,y,name,1,1,0.5,0.,5,2)
                #vc.label(x+378.941,y+457.031,name,1,1,0.5,0.5,2)
                #vc.exec_menu("Change Text")
        h.pop_section()
        vc.flush()
cell = RealisticNeuron(-70)
vc = h.PlotShape()
vc.exec_menu('Show Diam')
vc.color(2)
vc.flush()
for i in xrange(locE.size):
    print i, locE[i]
    sec = cell.dend[locE[i]]
    n = int(sec.n3d())
    L = np.zeros(n)
    for iseg in xrange(1,n):
        #print sec.x3d(iseg), sec.y3d(iseg), sec.z3d(iseg)
        L[iseg] = L[iseg-1] + np.sqrt(np.sum(np.square([sec.x3d(iseg) - sec.x3d(iseg-1),sec.y3d(iseg)-sec.y3d(iseg-1),sec.z3d(iseg)-sec.z3d(iseg-1)])))
    rL = L/L[-1]
    iseg = np.nonzero(np.greater(rL,posE[i]))[0][0]
    print "target seg ", iseg
    rl = rL[iseg] - posE[i]
    x = sec.x3d(iseg) + (sec.x3d(iseg+1)-sec.x3d(iseg)) * rl
    y = sec.y3d(iseg) + (sec.y3d(iseg+1)-sec.y3d(iseg)) * rl
    vc.mark(x,y,"o",mks,2,bs)
    vc.label(x+1,y+1,str(locE[i]))
for i in xrange(locI.size):
    print i, locI[i]
    sec = cell.dend[locI[i]]
    n = int(sec.n3d())
    L = np.zeros(n)
    for iseg in xrange(1,n):
        #print sec.x3d(iseg), sec.y3d(iseg), sec.z3d(iseg)
        L[iseg] = L[iseg-1] + np.sqrt(np.sum(np.square([sec.x3d(iseg) - sec.x3d(iseg-1),sec.y3d(iseg)-sec.y3d(iseg-1),sec.z3d(iseg)-sec.z3d(iseg-1)])))
    rL = L/L[-1]
    iseg = np.nonzero(np.greater(rL,posI[i]))[0][0]
    rl = rL[iseg] - posI[i]
    x = sec.x3d(iseg) + (sec.x3d(iseg+1)-sec.x3d(iseg)) * rl
    y = sec.y3d(iseg) + (sec.y3d(iseg+1)-sec.y3d(iseg)) * rl
    vc.mark(x,y,"o",mks,3,bs)
for i in xrange(75):
    cell.dend[i].push()
    vc.color(3)
    h.pop_section()
for i in xrange(75,135):
    cell.dend[i].push()
    vc.color(2)
    h.pop_section()
apical_range = [0,2,14,28,30,32,40,44,52,60,66,72,74,78,84,90,92,98]
for i in apical_range:
    cell.dend[i].push()
    vc.color(4)
    h.pop_section()
    

#vc.menu_tool('test', p)
#vc.exec_menu('test')
done = False
colored = np.ones(cell.ndend)
while not done:
    secstr = input("enter section index or 0 to exit\n")
    if secstr == '':
        done = True
    else:
        secnum = int(secstr)
        if secnum == 0:
            done = True
        else:
            if colored[secnum] == 0:
                if secnum < 75:
                    cell.dend[secnum].push()
                    vc.color(3)
                    h.pop_section()
                else:
                    cell.dend[secnum].push()
                    vc.color(2)
                    h.pop_section()
                colored[secnum] = 1
            else:
                cell.dend[secnum].push()
                vc.color(1)
                h.pop_section()
                colored[secnum] = 0
