from neuroAlter import prepCell, bproceed, leaky
import numpy as np
import sys, os, getopt
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot, gridspec
from shutil import copy
from datetime import datetime
import time
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
    
def getRS(pred,b,rSS = False):
    vpred = b-pred
    SSres = np.dot(vpred,vpred)
    vtot = b-np.average(b)
    SStot = np.dot(vtot,vtot)
    if np.absolute(SStot) == 0:
        rs = 1
    else:
        rs = 1- SSres/SStot
    if rSS:
        return rs, SSres, SStot
    else:
        return rs
        
def getKfig(fign,t,k,rs,intercept,coef,b,plim,directory,savePlot,EItype):
    pyplot.figure(fign,figsize = (8,4))
    #gs = gridspec.GridSpec(3, 1, height_ratios=[0.5, 1.5, 1], hspace = 0.18)
    gs = gridspec.GridSpec(2, 1, height_ratios=[2.0, 1.0], hspace = 0.0)
    #ax = pyplot.subplot(gs[0])
    #ax.plot(t,k,lw=1.0)
    ymin = -5
    ymax = 5
    iexpand = 1
    while np.logical_and(k > ymin, k < ymax).astype('float').sum()/k.size < plim:
        iexpand = iexpand + 1
        ymin = ymin - 5 
        ymax = ymax + 5 
    #ax.set_ylim(2, 50)
    #ax.set_yticks([10,30,50])
    ax0_0 = pyplot.subplot(gs[0])
    line0 = ax0_0.plot(t,k,lw=1.0,c='#E74C3C',label='k'+EItype)
    #line1 = ax0_0.plot(t,intercept,lw=1.5,c='darkorange',label='intercept');
    ax0_0.set_ylim(ymin, ymax)
    ax0_0.set_ylabel('k'+EItype + r' $(mV^{-1})$')
    ax0_1 = ax0_0.twinx()
    lines0 = ax0_1.plot(t,coef.T,':g',lw=1.0)
    lines1 = ax0_1.plot(t,b.T,':c',lw=1.0)
    lines0[0].set_label(r'$v_'+EItype[0]+'v_'+EItype[1] + '$')
    lines1[0].set_label(r'$v_{'+EItype+'}-'+ 'v_{'+EItype[0]+'}v_{'+EItype[1] + '}$')
    ylen = ymax - ymin
    tick0 = ymin + np.ceil(ylen*0.2/5) * 5
    tick1 = ymax - np.ceil(ylen*0.2/5) * 5
    ax0_0.set_yticks([tick0, ymin+ylen/2.0, tick1])
    
    #ax.spines['bottom'].set_visible(False)
    #ax.spines['right'].set_visible(False)
    #ax.spines['top'].set_visible(False)
    #ax.tick_params(bottom='off',labelbottom='off')
    #ax0_0.spines['top'].set_visible(False)
    ax0_1.spines['top'].set_visible(False)
    #ax0_0.tick_params(labeltop='off', labelbottom='off', bottom = 'off')
    ax0_1.tick_params(labelbottom='off', bottom = 'off')
    #pyplot.setp(ax.get_xticklabels(), visible=False)
    pyplot.setp(ax0_0.get_xticklabels(), visible=False)
    pyplot.setp(ax0_1.get_xticklabels(), visible=False)
    #ax.xaxis.tick_top()
    #d = 0.015
    #kwargs = dict(transform=ax.transAxes, color='k', clip_on = False)
    #ax.plot((-d,d),(-5*d,+5*d), **kwargs)
    #ax.plot((1-d,1+d),(-d,+d), **kwargs)

    #kwargs.update(transform=ax0_0.transAxes)
    #ax0_0.plot((-d,d),(1-2*d,1+2*d), **kwargs)
    #ax0_0.plot((1-d,1+d),(1-d,1+d), **kwargs)

    ax0_1.set_ylabel(r'$v_'+EItype[0]+'v_'+EItype[1] + ' (mV^2)$')
    ax0_1.ticklabel_format(axis='y',style='sci',scilimits=(-2,2))
    #ax0_1.yaxis.set_major_formatter(mpl.ticker.ScalarFormatter(useMathText=True, useOffset=False))

    ax0_0.legend()
    ax0_1.legend()

    ax1_0 = pyplot.subplot(gs[1], sharex = ax0_0)
    ax1_1 = ax1_0.twinx()

    line1, = ax1_0.plot(t,intercept,lw=1,c='darkorange', label='intercept');
    line2, = ax1_1.plot(t,rs,lw=1,c='#34495E',label=r'$R^{2}$');
    ax1_0.set_ylim((-0.2, 1.2))
    ax1_1.set_ylim((-0.2, 1.2))
    #pyplot.setp(ax1_1.get_yticklabels(), visible=False)
    ax1_0.set_xlabel('time (ms)')
    ax1_0.set_ylabel('intercept (mV)')
    ax1_1.set_ylabel(r'$R^{2}$')
    ax1_0.spines['top'].set_visible(False)
    ax1_1.spines['top'].set_visible(False)

    icross = np.diff(np.sign(coef[0,:])).nonzero()[0] + 1
    jcross = np.diff(np.sign(b[0,:])).nonzero()[0] + 1
    dy = (ymax-ymin)*0.035
    dx = (t[-1]-t[0])*0.015
    for i in icross:
        ax0_0.arrow(t[i], ymin-1.1*dy, 0.0, +dy, color='g',width = 0.0, head_length=0.9*dy,head_width = dx, clip_on=False, length_includes_head = True)
    for j in jcross:
        ax0_0.arrow(t[j], ymin-1.1*dy, 0.0, +dy, color='c',width = 0.0, head_length=0.9*dy,head_width = dx, clip_on=False, length_includes_head = True)
    #ax1_0.legend()
    #ax1_1.legend()
    if savePlot:
        pyplot.savefig(directory+'/'+fign+'.png',format='png',bbox_inches='tight',dpi=rdpi);
        pyplot.close() 
    else:
        pyplot.show()
def getWfig(fign,k,coef,b,directory,savePlot,EItype,it=-1):
    pyplot.figure(fign,figsize = (8,4))
    if it == -1:
        it = np.argmax(np.absolute(coef).max(0))
    pred = coef[:,it]*k[it]
    pyplot.figure(wfigname,figsize=(8,4))
    pyplot.plot(coef[:,it],b[:,it],'*k',ls='None',ms=8.5)
    rs = getRS(pred,b[:,it])
    coef = np.hstack((0,coef[:,it]))
    pred = np.hstack((0,pred))
    sid = coef.argsort()
    pyplot.plot(coef[sid],pred[sid])
    pyplot.xlabel(r'$v_'+EItype[0]+'v_'+EItype[1] + ' (mV^2)$')
    pyplot.title(r'$R^2 = ' + '%.2f'%rs + '$')
    pyplot.ylabel('k'+EItype+ r' $(mV^{-1})$')

    if np.absolute(pred0).max()> np.absolute(coef).max()*3:
        if np.absolute(coef.max()) > np.absolute(coef.min()):
            xmin = coef.min()
            xmax = np.absolute(pred0).max()
        else:
            xmin = -np.absolute(pred0).max()
            xmax = coef.max()
    else:
        xmin = coef.min()
        xmax = coef.max()
    xlen = xmax-xmin
    pyplot.xlim(xmin-xlen*0.2,xmax+xlen*0.2)
    pyplot.ticklabel_format(axis='both',style='sci',scilimits=(-2,2))

    if savePlot:
        pyplot.savefig(directory+'/'+fign+'.png',format='png',bbox_inches='tight',dpi=rdpi);
        pyplot.close() 
    else:
        pyplot.show()
        
argv = sys.argv[1:]
try:
    opts, args = getopt.getopt(argv,'ps',['vid=','vmin=','dV=','dt=','seed='])
except getopt.GetoptError:
    print "err"
    sys.exit(2)
print opts
opt0 = [opt[0] for opt in opts]
if '--vmin' not in opt0:
    print "specify vmin"
if '--vid' not in opt0:
    print "specify vid"
if '--seed' not in opt0:
    print "specify seed"
if '--dt' not in opt0:
    print "specify dt"
if '--dV' not in opt0:
    print "specify dV"
for opt,arg in opts:
    if opt == '--vmin':
        vmin = float(arg)
    elif opt == '--vid':
        vid = int(arg)
    elif opt == '--dt':
        dt = float(arg)
    elif opt == '--seed':
        seed = int(arg)
    elif opt == '--dV':
        dV = float(arg)

v0 = vmin + dV*vid
print "v0 == ",v0
#v0 = -67
#dt = 2
#seed = int(time.mktime(datetime.now().timetuple())) 
#print "seed = ", seed
nloc = 4
run_t = 240
#run_t = 22
tstep = 0.1
trans = 110
#trans = 2 
pos0 = 0.5
gList0 = np.array([4,8,16,32]) * 1e-4
plim = 0.7
#gList0 = np.array([1,3]) * 1e-4
rI = 1e1
testEE = True
#testEE = False 
testII = True
testIE = True
#testIE = False 
testEI = True
#testEI = False 
singleOnly = False
directory0 = 'cK' + str(int(v0)) + '-' + str(0) + '-' + str(seed)
singletE0Filename = directory0+'/singletE0.npy'
singletI0Filename = directory0+'/singletI0.npy'
if not os.path.isfile(singletE0Filename):
    if not dt==0 and (testIE or testEE):
        dt = 0
        testEE = True
        singleOnly = True
        testII = False 
        testIE = False 
        testEI = False 
        print " need to get Esinglet of dt=0 first, please change dt and rerun after finish"
    
if not os.path.isfile(singletI0Filename):
    if not dt==0 and (testEI or testII):
        dt = 0
        testII = True
        singleOnly = True
        testEE = False 
        testIE = False 
        testEI = False 
        print " need to get Isinglet of dt=0 first, please change dt and rerun after finish"
directory = 'cK' + str(int(v0)) + '-' + str(int(dt)) + '-' + str(seed)
print directory0
print directory

if not os.path.isdir(directory):
    os.mkdir('./'+directory)
dir0 = os.path.dirname(__file__)
filename = os.path.join(dir0,'confirmK.py')
copydir = os.path.join(dir0,directory)
copy(filename,copydir)
savePlot = True
plotData = True
rdpi = 300
run_nt = np.round(run_t/tstep + 1).astype('int')
ntrans = np.round(trans/tstep).astype('int')
print "run steps: ", run_nt, tstep, "ms per step"
np.random.seed(seed)
locE = np.random.randint(52,134,nloc)
locI = np.random.randint(0,52,nloc)
ng = gList0.size
idt = np.round(dt/tstep).astype('int')
end_t = run_t-dt
iend = run_nt-idt
cell, vSL, sL = prepCell([gList0[0]], [locE[0]], [pos0], 1, -70)
RList = np.zeros((1,2))
t = np.arange(run_nt)*tstep
leakyV, _ = leaky(cell, v0, sL, RList, vSL, 1, trans, run_t + trans, tstep, [locE[0]], [pos0])
leakyV = leakyV[ntrans:ntrans+run_nt]
leakyV.dump(directory+'/'+'leakyV.pickle')
lfigname = 'leakyV'
pyplot.figure()
pyplot.plot(t[idt:],leakyV[:iend])
if savePlot:
    pyplot.savefig(directory+'/'+lfigname+'.png',format='png',bbox_inches='tight',dpi=rdpi);
    pyplot.close()
if testEE:
    #Singlet
    print "singleE"
    if plotData:
        sfigname = 'Esinglets'
        pyplot.figure(sfigname,figsize = (8,4))
        ax1 = pyplot.subplot(1,2,1)
        ax2 = pyplot.subplot(1,2,2)
    singletE = np.empty((nloc,ng,run_nt))
    sfired = np.empty((nloc,ng))
    for i in xrange(nloc):
        print "single ", i
        loc = np.repeat(locE[i],ng)
        pos = np.repeat(pos0,ng)
        cell, vSL, sL = prepCell(gList0, loc, pos, ng, -70)
        spikeTime = np.array([[]])
        for ig in xrange(ng):
            gList = gList0
            RList = np.zeros((ng,2))
            RList[ig,:] = getGH(dt,gList[ig])
            print RList[ig,:]
            vtmp, sfired[i,ig] = bproceed(cell, v0, sL, gList, RList, vSL, spikeTime, ng, [ig], trans, trans + end_t, 0, tstep, '')
            if sfired[i,ig]:
                quit()
            singletE[i,ig,idt:] = vtmp[ntrans:ntrans+iend] - leakyV[:iend]
            ax1.plot(t[idt:],singletE[i,ig,idt:],':r')
            ax2.plot(t[idt:],vtmp[ntrans:ntrans+iend],':r')
            ax2.plot(t[idt:],leakyV[:iend],':b')
    if savePlot:
        pyplot.savefig(directory+'/'+sfigname+'.png',format='png',bbox_inches='tight',dpi=rdpi);
        pyplot.close()
    if idt == 0:
        singletE.dump(directory0+'/singletE0.npy')
        singletE.dump(directory+'/singletE.npy')
    else:
        singletE.dump(directory+'/singletE.npy')
    if not singleOnly:
        print "doubletEE"
        #Doublet
        gList = np.hstack((np.array([np.repeat(gList0,ng)]).transpose(),np.tile(gList0,(1,ng)).transpose())).reshape(ng,ng,2)
        print gList
        spikeTime = np.array([np.array([]),np.array([0.0])])  
        sel = [0, 1]
        pos = np.repeat(pos0,2)
        singletE0 = np.load(singletE0Filename)
        singletE = np.load(directory+'/singletE.npy')
        doubletEE = np.empty((nloc,nloc,ng*ng,run_nt))
        doubletEE0 = np.empty((nloc,nloc,ng*ng,run_nt))
        dfiredEE = np.empty((nloc,nloc,ng,ng))
        kEE = np.empty((nloc,nloc,run_nt))
        rsEE = np.empty((nloc,nloc,run_nt))
        interceptEE = np.empty((nloc,nloc,run_nt))
        kEE0 = np.empty((nloc,nloc,run_nt))
        rsEE0 = np.empty((nloc,nloc,run_nt))
        interceptEE0 = np.empty((nloc,nloc,run_nt))
        for i in xrange(nloc):
            for j in xrange(nloc):
                loc = np.array([locE[i], locE[j]])
                for ig in xrange(ng):
                    v1 = singletE[i,ig,idt:]
                    if plotData:
                        bfigname = 'doublets-E'+str(i)+'-E'+str(j)+'-ig'+str(ig)
                        bfig = pyplot.figure(bfigname,figsize = (8,4))
                        ax1 = bfig.add_subplot(111)         
                        ax1.plot(t[idt:],v1,'r',lw=3)
                        if v0 == -70:
                            bfign0 = 'doublets0-E'+str(i)+'-E'+str(j)+'-ig'+str(ig)
                            bfig0 = pyplot.figure(bfign0,figsize = (8,4))
                            ax2 = bfig0.add_subplot(111)         
                            ax2.plot(t[idt:],singletE0[i,ig,idt:],'r',lw=3)
                    for jg in xrange(ng):
                        print gList[ig,jg,:]
                        RList = np.zeros((2,2))
                        RList[0,:] = getGH(dt,gList[ig,jg,0])
                        print RList
                        cell, vSL, sL = prepCell(gList[ig,jg,:], loc, pos, 2, -70)
                        vtmp, dfiredEE[i,j,ig,jg] = bproceed(cell, v0, sL, gList[ig,jg,:], RList, vSL, spikeTime, 2, sel, trans, trans + end_t, 0, tstep, '')
                        if dfiredEE[i,j,ig,jg]:
                            quit()
                        vtmp = vtmp[ntrans:ntrans+iend] - leakyV[:iend]
                        doubletEE[i,j,ig*ng+jg,idt:] = vtmp
                        addv = v1 + singletE0[j,jg,:iend]
                        vs = vtmp - addv
                        if plotData:
                            ax1.plot(t[idt:],vtmp,'k',lw=2)
                            ax1.plot(t[idt:],singletE0[j,jg,:iend],'b',lw=1)
                            ax1.plot(t[idt:],addv,':g', lw =1)
                        if v0 == -70:
                            spikeTime0 = np.array([np.array([0.0]),np.array([dt])])  
                            RList = np.zeros((2,2))
                            vtmp, _ = bproceed(cell, v0, sL, gList[ig,jg,:], RList, vSL, spikeTime0, 2, sel, trans, trans + run_t, 0, tstep, '')
                            vtmp = vtmp[ntrans+idt:ntrans+run_nt] + 70
                            doubletEE0[i,j,ig*ng+jg,idt:] = vtmp
                            addv = singletE0[j,jg,:iend] + singletE0[i,ig,idt:]
                            vs = vtmp - addv
                            if plotData:
                                ax2.plot(t[idt:],vtmp,'k',lw=2)
                                ax2.plot(t[idt:],singletE0[j,jg,:iend],'b',lw=1)
                                ax2.plot(t[idt:],addv,':g', lw =1)
                            
                    if savePlot and plotData:
                        pyplot.figure(bfigname)
                        pyplot.savefig(directory+'/'+bfigname+'.png',format='png',bbox_inches='tight',dpi=rdpi);
                        if v0 == -70:
                            pyplot.figure(bfign0)
                            pyplot.savefig(directory+'/'+bfign0+'.png',format='png',bbox_inches='tight',dpi=rdpi);
                        pyplot.close()
                #k
                v1 = np.repeat(singletE[i,:,idt:],ng,axis=0).reshape(ng*ng,iend)
                v2 = np.vstack([singletE0[j,:,:iend],]*ng)
                addv = v1+v2
                coef = v1*v2
                b = doubletEE[i,j,:,idt:] - addv
                arrived = False
                for it in xrange(iend):
                    coefmat0 = np.vstack((coef[:,it],np.zeros((ng*ng)))).T
                    coefmat1 = np.vstack((coef[:,it],np.ones((ng*ng)))).T
                    if np.absolute(coef[:,it]).sum() == 0:
                        print 'inf slope due to delayed arrival of EPSP at soma'
                        kEE[i,j,idt+it] = 0
                        c = [0,0]
                    else:
                        if not arrived and np.std(coef[:,it]) > 0:
                            arrived = True
                        kEE[i,j,idt+it] = np.linalg.lstsq(coefmat0,b[:,it])[0][0]
                        c = np.linalg.lstsq(coefmat1,b[:,it])[0]
                    interceptEE[i,j,idt+it] = c[1]
                    pred0 = kEE[i,j,idt+it]*coef[:,it]
                    pred1 = c[0]*coef[:,it]+c[1]
                    rsEE[i,j,idt+it] = getRS(pred0,b[:,it])

                if plotData:
                    kfigname = 'k-E'+str(i)+'-E'+str(j)
                    wfigname = 'midK-E'+str(i)+'-E'+str(j)
                    getKfig(kfigname,t[idt:],kEE[i,j,idt:],rsEE[i,j,idt:],interceptEE[i,j,idt:],coef,b,plim,directory,savePlot,'EE')
                    getWfig(wfigname,kEE[i,j,idt:],coef,b,directory,savePlot,'EE')
                    wfigname = 'worstK-E'+str(i)+'-E'+str(j)
                    it = np.argmin(rsEE[i,j,idt:])
                    getWfig(wfigname,kEE[i,j,idt:],coef,b,directory,savePlot,'EE',it)
                if v0 == -70:
                    v1 = np.repeat(singletE0[i,:,idt:],ng,axis=0).reshape(ng*ng,iend)
                    v2 = np.vstack([singletE0[j,:,:iend],]*ng)
                    addv = v1+v2
                    coef = v1*v2
                    b = doubletEE0[i,j,:,idt:] - addv
                    arrived = False
                    for it in xrange(iend):
                        coefmat0 = np.vstack((coef[:,it],np.zeros((ng*ng)))).T
                        coefmat1 = np.vstack((coef[:,it],np.ones((ng*ng)))).T
                        if np.absolute(coef[:,it]).sum() == 0:
                            print 'inf slope due to delayed arrival of EPSP at soma'
                            kEE0[i,j,idt+it] = 0
                            c = [0,0]
                        else:
                            if not arrived and np.std(coef[:,it]) > 0:
                                arrived = True
                            kEE0[i,j,idt+it] = np.linalg.lstsq(coefmat0,b[:,it])[0][0]
                            c = np.linalg.lstsq(coefmat1,b[:,it])[0]
                        interceptEE0[i,j,idt+it] = c[1]
                        pred0 = kEE0[i,j,idt+it]*coef[:,it]
                        pred1 = c[0]*coef[:,it]+c[1]
                        rsEE0[i,j,idt+it] = getRS(pred0,b[:,it])

                    if plotData:
                        kfign0 = 'k0-E'+str(i)+'-E'+str(j)
                        wfign0 = 'midK0-E'+str(i)+'-E'+str(j)
                        getKfig(kfign0,t[idt:],kEE0[i,j,idt:],rsEE0[i,j,idt:],interceptEE0[i,j,idt:],coef,b,plim,directory,savePlot,'EE')
                        getWfig(wfign0,kEE0[i,j,idt:],coef,b,directory,savePlot,'EE')
                        wfign0 = 'worstK0-E'+str(i)+'-E'+str(j)
                        it = np.argmin(rsEE0[i,j,idt:])
                        getWfig(wfign0,kEE0[i,j,idt:],coef,b,directory,savePlot,'EE',it)

        if v0 == -70:
            doubletEE0.dump(directory+'/doubletEE0.npy')
            kEE0.dump(directory+'/kEE0.npy')
            rsEE0.dump(directory+'/rsEE0.npy') 
            interceptEE0.dump(directory+'/interceptEE0.npy')

        doubletEE.dump(directory+'/doubletEE.npy')
        kEE.dump(directory+'/kEE.npy')
        rsEE.dump(directory+'/rsEE.npy') 
        interceptEE.dump(directory+'/interceptEE.npy')
        dfiredEE.dump(directory+'/dfiredEE.npy')

if testII:
    #Singlet
    print "singleI"
    if plotData:
        sfigname = 'Isinglets'
        pyplot.figure(sfigname,figsize = (8,4))
        ax1 = pyplot.subplot(1,2,1)
        ax2 = pyplot.subplot(1,2,2)
    singletI = np.empty((nloc,ng,run_nt))
    sfired = np.empty((nloc,ng))
    for i in xrange(nloc):
        print "single ", i
        loc = np.repeat(locI[i],ng)
        pos = np.repeat(pos0,ng)
        cell, vSL, sL = prepCell(-gList0*rI, loc, pos, ng, -70)
        spikeTime = np.array([[]])
        for ig in xrange(ng):
            gList = -gList0*rI
            RList = np.zeros((ng,2))
            RList[ig,:] = getGH(dt,gList[ig])
            print RList[ig,:]
            st = np.array([[]])
            vtmp, sfired[i,ig] = bproceed(cell, v0, sL, gList, RList, vSL, spikeTime, ng, [ig], trans, trans + end_t, 0, tstep, '')
            singletI[i,ig,idt:] = vtmp[ntrans:ntrans+iend] - leakyV[:iend]
            ax1.plot(t[idt:],singletI[i,ig,idt:],':r')
            ax2.plot(t[idt:],vtmp[ntrans:ntrans+iend],':r')
            ax2.plot(t[idt:],leakyV[:iend],':b')
    if savePlot:
        pyplot.savefig(directory+'/'+sfigname+'.png',format='png',bbox_inches='tight',dpi=rdpi);
        pyplot.close()
    else:
        pyplot.show()
    if idt == 0:
        singletI.dump(directory0+'/singletI0.npy')
        singletI.dump(directory+'/singletI.npy')
    else:
        singletI.dump(directory+'/singletI.npy')
        
    if not singleOnly:
        print "doubletII"
        #Doublet
        gList = np.hstack((np.array([np.repeat(-gList0,ng)]).transpose(),np.tile(-gList0,(1,ng)).transpose())).reshape(ng,ng,2)*rI
        print gList
        spikeTime = np.array([np.array([]),np.array([0.0])])  
        sel = [0, 1]
        pos = np.repeat(pos0,2)
        singletI0 = np.load(singletI0Filename)
        singletI = np.load(directory+'/singletI.npy')
        doubletII = np.empty((nloc,nloc,ng*ng,run_nt))
        doubletII0 = np.empty((nloc,nloc,ng*ng,run_nt))
        dfiredII = np.empty((nloc,nloc,ng,ng))
        kII = np.empty((nloc,nloc,run_nt))
        rsII = np.empty((nloc,nloc,run_nt))
        interceptII = np.empty((nloc,nloc,run_nt))
        kII0 = np.empty((nloc,nloc,run_nt))
        rsII0 = np.empty((nloc,nloc,run_nt))
        interceptII0 = np.empty((nloc,nloc,run_nt))
        for i in xrange(nloc):
            for j in xrange(nloc):
                loc = np.array([locI[i], locI[j]])
                for ig in xrange(ng):
                    v1 = singletI[i,ig,idt:]
                    if plotData:
                        bfigname = 'doublets-I'+str(i)+'-I'+str(j)+'-ig'+str(ig)
                        bfig = pyplot.figure(bfigname,figsize = (8,4))
                        ax1 = bfig.add_subplot(111)         
                        ax1.plot(t[idt:],v1,'r',lw =3)
                        if v0 == -70:
                            bfign0 = 'doublets0-I'+str(i)+'-I'+str(j)+'-ig'+str(ig)
                            bfig0 = pyplot.figure(bfign0,figsize = (8,4))
                            ax2 = bfig0.add_subplot(111)         
                            ax2.plot(t[idt:],singletI0[i,ig,idt:],'r', lw=3)
                    for jg in xrange(ng):
                        print gList[ig,jg,:]
                        RList = np.zeros((2,2))
                        RList[0,:] = getGH(dt,gList[ig,jg,0])
                        print RList
                        cell, vSL, sL = prepCell(gList[ig,jg,:], loc, pos, 2, -70)
                        vtmp, dfiredII[i,j,ig,jg] = bproceed(cell, v0, sL, gList[ig,jg,:], RList, vSL, spikeTime, 2, sel, trans, trans + end_t, 0, tstep, '')
                        vtmp = vtmp[ntrans:ntrans+iend] - leakyV[:iend]
                        doubletII[i,j,ig*ng+jg,idt:] = vtmp
                        addv = v1 + singletI0[j,jg,:iend]
                        vs = vtmp - addv
                        if plotData:
                            ax1.plot(t[idt:],vtmp,'k',lw=2)
                            ax1.plot(t[idt:],singletI0[j,jg,:iend],'b',lw=1)
                            ax1.plot(t[idt:],addv,':g',lw=1)
                        if v0 == -70:
                            spikeTime0 = np.array([np.array([0.0]),np.array([dt])])  
                            RList = np.zeros((2,2))
                            vtmp, _ = bproceed(cell, v0, sL, gList[ig,jg,:], RList, vSL, spikeTime0, 2, sel, trans, trans + run_t, 0, tstep, '')
                            vtmp = vtmp[ntrans+idt:ntrans+run_nt] + 70
                            doubletII0[i,j,ig*ng+jg,idt:] = vtmp
                            addv = singletI0[j,jg,:iend] + singletI0[i,ig,idt:]
                            vs = vtmp - addv
                            if plotData:
                                ax2.plot(t[idt:],vtmp,'k',lw=2)
                                ax2.plot(t[idt:],singletI0[j,jg,:iend],'b',lw=1)
                                ax2.plot(t[idt:],addv,':g', lw =1)

                    if savePlot and plotData:
                        pyplot.figure(bfigname)
                        pyplot.savefig(directory+'/'+bfigname+'.png',format='png',bbox_inches='tight',dpi=rdpi);
                        if v0 == -70:
                            pyplot.figure(bfign0)
                            pyplot.savefig(directory+'/'+bfign0+'.png',format='png',bbox_inches='tight',dpi=rdpi);
                        pyplot.close()
                #k
                v1 = np.repeat(singletI[i,:,idt:],ng,axis=0).reshape(ng*ng,iend)
                v2 = np.vstack([singletI0[j,:,:iend],]*ng)
                addv = v1+v2
                coef = v1*v2
                b = doubletII[i,j,:,idt:] - addv
                arrived = False
                for it in xrange(iend):
                    coefmat0 = np.vstack((coef[:,it],np.zeros((ng*ng)))).T
                    coefmat1 = np.vstack((coef[:,it],np.ones((ng*ng)))).T
                    if np.absolute(coef[:,it]).sum() == 0:
                        print 'inf slope due to delayed arrival of IPSP at soma'
                        kII[i,j,idt+it] = 0
                        c = [0,0]
                    else:
                        if not arrived and np.std(coef[:,it]) > 0:
                            arrived = True
                        kII[i,j,idt+it] = np.linalg.lstsq(coefmat0,b[:,it])[0][0]
                        c = np.linalg.lstsq(coefmat1,b[:,it])[0]
                    interceptII[i,j,idt+it] = c[1]
                    pred0 = kII[i,j,idt+it]*coef[:,it]
                    pred1 = c[0]*coef[:,it]+c[1]
                    rsII[i,j,idt+it] = getRS(pred0,b[:,it])

                if plotData:
                    kfigname = 'k-I'+str(i)+'-I'+str(j)
                    wfigname = 'midK-I'+str(i)+'-I'+str(j)
                    getKfig(kfigname,t[idt:],kII[i,j,idt:],rsII[i,j,idt:],interceptII[i,j,idt:],coef,b,plim,directory,savePlot,'II')
                    getWfig(wfigname,kII[i,j,idt:],coef,b,directory,savePlot,'II')
                    wfigname = 'worstK-I'+str(i)+'-I'+str(j)
                    it = np.argmin(rsII[i,j,idt:])
                    getWfig(wfigname,kII[i,j,idt:],coef,b,directory,savePlot,'II',it)
                if v0 == -70:
                    v1 = np.repeat(singletI0[i,:,idt:],ng,axis=0).reshape(ng*ng,iend)
                    v2 = np.vstack([singletI0[j,:,:iend],]*ng)
                    addv = v1+v2
                    coef = v1*v2
                    b = doubletII0[i,j,:,idt:] - addv
                    arrived = False
                    for it in xrange(iend):
                        coefmat0 = np.vstack((coef[:,it],np.zeros((ng*ng)))).T
                        coefmat1 = np.vstack((coef[:,it],np.ones((ng*ng)))).T
                        if np.absolute(coef[:,it]).sum() == 0:
                            print 'inf slope due to delayed arrival of IPSP at soma'
                            kII0[i,j,idt+it] = 0
                            c = [0,0]
                        else:
                            if not arrived and np.std(coef[:,it]) > 0:
                                arrived = True
                            kII0[i,j,idt+it] = np.linalg.lstsq(coefmat0,b[:,it])[0][0]
                            c = np.linalg.lstsq(coefmat1,b[:,it])[0]
                        interceptII0[i,j,idt+it] = c[1]
                        pred0 = kII0[i,j,idt+it]*coef[:,it]
                        pred1 = c[0]*coef[:,it]+c[1]
                        rsII0[i,j,idt+it] = getRS(pred0,b[:,it])

                    if plotData:
                        kfign0 = 'k0-I'+str(i)+'-I'+str(j)
                        wfign0 = 'midK0-I'+str(i)+'-I'+str(j)
                        getKfig(kfign0,t[idt:],kII0[i,j,idt:],rsII0[i,j,idt:],interceptII0[i,j,idt:],coef,b,plim,directory,savePlot,'II')
                        getWfig(wfign0,kII0[i,j,idt:],coef,b,directory,savePlot,'II')
                        wfign0 = 'worstK0-I'+str(i)+'-I'+str(j)
                        it = np.argmin(rsII0[i,j,idt:])
                        getWfig(wfign0,kII0[i,j,idt:],coef,b,directory,savePlot,'II',it)

        if v0 == -70:
            doubletII0.dump(directory+'/doubletII0.npy')
            kII0.dump(directory+'/kII0.npy')
            rsII0.dump(directory+'/rsII0.npy') 
            interceptII0.dump(directory+'/interceptII0.npy')

        doubletII.dump(directory+'/doubletII.npy')
        dfiredII.dump(directory+'/dfiredII.npy')
        kII.dump(directory+'/kII.npy')
        rsII.dump(directory+'/rsII.npy') 
        interceptII.dump(directory+'/interceptII.npy')
if testIE and not singleOnly:
    print "doubletIE"
    #Doublet
    gList = np.hstack((np.array([np.repeat(-gList0*rI,ng)]).transpose(),np.tile(gList0,(1,ng)).transpose())).reshape(ng,ng,2)
    print gList
    spikeTime = np.array([np.array([]),np.array([0.0])])  
    sel = [0, 1]
    pos = np.repeat(pos0,2)
    singletI0 = np.load(singletI0Filename)
    singletE0 = np.load(singletE0Filename)
    singletI = np.load(directory+'/singletI.npy')
    doubletIE = np.empty((nloc,nloc,ng*ng,run_nt))
    doubletIE0 = np.empty((nloc,nloc,ng*ng,run_nt))
    dfiredIE = np.empty((nloc,nloc,ng,ng))
    kIE = np.empty((nloc,nloc,run_nt))
    rsIE = np.empty((nloc,nloc,run_nt))
    interceptIE = np.empty((nloc,nloc,run_nt))
    kIE0 = np.empty((nloc,nloc,run_nt))
    rsIE0 = np.empty((nloc,nloc,run_nt))
    interceptIE0 = np.empty((nloc,nloc,run_nt))
    for i in xrange(nloc):
        for j in xrange(nloc):
            loc = np.array([locI[i], locE[j]])
            for ig in xrange(ng):
                v1 = singletI[i,ig,idt:]
                if plotData:
                    bfigname = 'doublets-I'+str(i)+'-E'+str(j)+'-ig'+str(ig)
                    bfig = pyplot.figure(bfigname,figsize = (8,4))
                    ax1 = bfig.add_subplot(111)         
                    ax1.plot(t[idt:],v1,'b', lw =3)
                    if v0 == -70:
                        bfign0 = 'doublets0-I'+str(i)+'-E'+str(j)+'-ig'+str(ig)
                        bfig0 = pyplot.figure(bfign0,figsize = (8,4))
                        ax2 = bfig0.add_subplot(111)         
                        ax2.plot(t[idt:],singletI0[i,ig,idt:],'b', lw=3)
                for jg in xrange(ng):
                    print gList[ig,jg,:]
                    RList = np.zeros((2,2))
                    RList[0,:] = getGH(dt,gList[ig,jg,0])
                    print RList
                    cell, vSL, sL = prepCell(gList[ig,jg,:], loc, pos, 2, -70)
                    vtmp, dfiredIE[i,j,ig,jg] = bproceed(cell, v0, sL, gList[ig,jg,:], RList, vSL, spikeTime, 2, sel, trans, trans + end_t, 0, tstep, '')
                    vtmp = vtmp[ntrans:ntrans+iend] - leakyV[:iend]
                    doubletIE[i,j,ig*ng+jg,idt:] = vtmp
                    addv = v1 + singletE0[j,jg,:iend]
                    vs = vtmp - addv
                    if plotData:
                        ax1.plot(t[idt:],vtmp,'k',lw=2)
                        ax1.plot(t[idt:],singletE0[j,jg,:iend],'r',lw=1)
                        ax1.plot(t[idt:],addv,':g', lw=1)
                    if v0 == -70:
                        spikeTime0 = np.array([np.array([0.0]),np.array([dt])])  
                        RList = np.zeros((2,2))
                        vtmp, _ = bproceed(cell, v0, sL, gList[ig,jg,:], RList, vSL, spikeTime0, 2, sel, trans, trans + run_t, 0, tstep, '')
                        vtmp = vtmp[ntrans+idt:ntrans+run_nt] + 70
                        doubletIE0[i,j,ig*ng+jg,idt:] = vtmp
                        addv = singletI0[i,ig,idt:] + singletE0[j,jg,:iend] 
                        vs = vtmp - addv
                        if plotData:
                            ax2.plot(t[idt:],vtmp,'k',lw=2)
                            ax2.plot(t[idt:],singletE0[j,jg,:iend],'r',lw=1)
                            ax2.plot(t[idt:],addv,':g', lw =1)

                if savePlot and plotData:
                    pyplot.figure(bfigname)
                    pyplot.savefig(directory+'/'+bfigname+'.png',format='png',bbox_inches='tight',dpi=rdpi);
                    if v0 == -70:
                        pyplot.figure(bfign0)
                        pyplot.savefig(directory+'/'+bfign0+'.png',format='png',bbox_inches='tight',dpi=rdpi);
                    pyplot.close()
            #k
            v1 = np.repeat(singletI[i,:,idt:],ng,axis=0).reshape(ng*ng,iend)
            v2 = np.vstack([singletE0[j,:,:iend],]*ng)
            addv = v1+v2
            coef = v1*v2
            b = doubletIE[i,j,:,idt:] - addv
            arrived = False
            for it in xrange(iend):
                coefmat0 = np.vstack((coef[:,it],np.zeros((ng*ng)))).T
                coefmat1 = np.vstack((coef[:,it],np.ones((ng*ng)))).T
                if np.absolute(coef[:,it]).sum() == 0:
                    print 'inf slope due to delayed arrival of IPSP at soma'
                    kIE[i,j,idt+it] = 0
                    c = [0,0]
                else:
                    if not arrived and np.std(coef[:,it]) > 0:
                        arrived = True
                    kIE[i,j,idt+it] = np.linalg.lstsq(coefmat0,b[:,it])[0][0]
                    c = np.linalg.lstsq(coefmat1,b[:,it])[0]
                interceptIE[i,j,idt+it] = c[1]
                pred0 = kIE[i,j,idt+it]*coef[:,it]
                pred1 = c[0]*coef[:,it]+c[1]
                rsIE[i,j,idt+it] = getRS(pred0,b[:,it])

            if plotData:
                kfigname = 'k-I'+str(i)+'-E'+str(j)
                wfigname = 'midK-I'+str(i)+'-E'+str(j)
                getKfig(kfigname,t[idt:],kIE[i,j,idt:],rsIE[i,j,idt:],interceptIE[i,j,idt:],coef,b,plim,directory,savePlot,'IE')
                getWfig(wfigname,kIE[i,j,idt:],coef,b,directory,savePlot,'IE')
                wfigname = 'worstK-I'+str(i)+'-E'+str(j)
                it = np.argmin(rsIE[i,j,idt:])
                getWfig(wfigname,kIE[i,j,idt:],coef,b,directory,savePlot,'IE',it)
            if v0 == -70:
                v1 = np.repeat(singletI0[i,:,idt:],ng,axis=0).reshape(ng*ng,iend)
                v2 = np.vstack([singletE0[j,:,:iend],]*ng)
                addv = v1+v2
                coef = v1*v2
                b = doubletIE0[i,j,:,idt:] - addv
                arrived = False
                for it in xrange(iend):
                    coefmat0 = np.vstack((coef[:,it],np.zeros((ng*ng)))).T
                    coefmat1 = np.vstack((coef[:,it],np.ones((ng*ng)))).T
                    if np.absolute(coef[:,it]).sum() == 0:
                        print 'inf slope due to delayed arrival of IPSP at soma'
                        kIE0[i,j,idt+it] = 0
                        c = [0,0]
                    else:
                        if not arrived and np.std(coef[:,it]) > 0:
                            arrived = True
                        kIE0[i,j,idt+it] = np.linalg.lstsq(coefmat0,b[:,it])[0][0]
                        c = np.linalg.lstsq(coefmat1,b[:,it])[0]
                    interceptIE0[i,j,idt+it] = c[1]
                    pred0 = kIE0[i,j,idt+it]*coef[:,it]
                    pred1 = c[0]*coef[:,it]+c[1]
                    rsIE0[i,j,idt+it] = getRS(pred0,b[:,it])

                if plotData:
                    kfign0 = 'k0-I'+str(i)+'-E'+str(j)
                    wfign0 = 'midK0-I'+str(i)+'-E'+str(j)
                    getKfig(kfign0,t[idt:],kIE0[i,j,idt:],rsIE0[i,j,idt:],interceptIE0[i,j,idt:],coef,b,plim,directory,savePlot,'IE')
                    getWfig(wfign0,kIE0[i,j,idt:],coef,b,directory,savePlot,'IE')
                    wfign0 = 'worstK0-I'+str(i)+'-E'+str(j)
                    it = np.argmin(rsIE0[i,j,idt:])
                    getWfig(wfign0,kIE0[i,j,idt:],coef,b,directory,savePlot,'IE',it)

    if v0 == -70:
        doubletIE0.dump(directory+'/doubletIE0.npy')
        kIE0.dump(directory+'/kIE0.npy')
        rsIE0.dump(directory+'/rsIE0.npy') 
        interceptIE0.dump(directory+'/interceptIE0.npy')

    doubletIE.dump(directory+'/doubletIE.npy')
    dfiredIE.dump(directory+'/dfiredIE.npy')
    kIE.dump(directory+'/kIE.npy')
    rsIE.dump(directory+'/rsIE.npy') 
    interceptIE.dump(directory+'/interceptIE.npy')

if testEI and not singleOnly:
    print "doubletEI"
    #Doublet
    gList = np.hstack((np.array([np.repeat(gList0,ng)]).transpose(),np.tile(-gList0*rI,(1,ng)).transpose())).reshape(ng,ng,2)
    print gList
    spikeTime = np.array([np.array([]),np.array([0.0])])  
    sel = [0, 1]
    pos = np.repeat(pos0,2)
    singletI0 = np.load(singletI0Filename)
    singletE0 = np.load(singletE0Filename)
    singletE = np.load(directory+'/singletE.npy')
    doubletEI = np.empty((nloc,nloc,ng*ng,run_nt))
    doubletEI0 = np.empty((nloc,nloc,ng*ng,run_nt))
    dfiredEI = np.empty((nloc,nloc,ng,ng))
    kEI = np.empty((nloc,nloc,run_nt))
    rsEI = np.empty((nloc,nloc,run_nt))
    interceptEI = np.empty((nloc,nloc,run_nt))
    kEI0 = np.empty((nloc,nloc,run_nt))
    rsEI0 = np.empty((nloc,nloc,run_nt))
    interceptEI0 = np.empty((nloc,nloc,run_nt))
    for i in xrange(nloc):
        for j in xrange(nloc):
            loc = np.array([locE[i], locI[j]])
            for ig in xrange(ng):
                v1 = singletE[i,ig,idt:]
                if plotData:
                    bfigname = 'doublets-E'+str(i)+'-I'+str(j)+'-ig'+str(ig)
                    bfig = pyplot.figure(bfigname,figsize = (8,4))
                    ax1 = bfig.add_subplot(111)         
                    ax1.plot(t[idt:],v1,'r',lw=3)
                    if v0 == -70:
                        bfign0 = 'doublets-E'+str(i)+'-I'+str(j)+'-ig'+str(ig)
                        bfig0 = pyplot.figure(bfign0,figsize = (8,4))
                        ax2 = bfig0.add_subplot(111)         
                        ax2.plot(t[idt:],singletE0[i,ig,idt:],'r',lw=3)
                for jg in xrange(ng):
                    print gList[ig,jg,:]
                    RList = np.zeros((2,2))
                    RList[0,:] = getGH(dt,gList[ig,jg,0])
                    print RList
                    cell, vSL, sL = prepCell(gList[ig,jg,:], loc, pos, 2, -70)
                    vtmp, dfiredEI[i,j,ig,jg] = bproceed(cell, v0, sL, gList[ig,jg,:], RList, vSL, spikeTime, 2, sel, trans, trans + end_t, 0, tstep, '')
                    vtmp = vtmp[ntrans:ntrans+iend] - leakyV[:iend]
                    doubletEI[i,j,ig*ng+jg,idt:] = vtmp
                    addv = v1 + singletI0[j,jg,:iend]
                    vs = vtmp - addv
                    if plotData:
                        ax1.plot(t[idt:],vtmp,'c',lw=2)
                        ax1.plot(t[idt:],singletI0[j,jg,:iend],'b',lw=1)
                        ax1.plot(t[idt:],addv,':g', lw=1)
                        ax1.plot(t[idt:],vs,':k', lw=0.5)
                    if v0 == -70:
                        spikeTime0 = np.array([np.array([0.0]),np.array([dt])])  
                        RList = np.zeros((2,2))
                        vtmp, _ = bproceed(cell, v0, sL, gList[ig,jg,:], RList, vSL, spikeTime0, 2, sel, trans, trans + run_t, 0, tstep, '')
                        vtmp = vtmp[ntrans+idt:ntrans+run_nt] + 70
                        doubletEI0[i,j,ig*ng+jg,idt:] = vtmp
                        addv = singletE0[i,ig,idt:] + singletI0[j,jg,:iend] 
                        vs = vtmp - addv
                        if plotData:
                            ax2.plot(t[idt:],vtmp,'c',lw=2)
                            ax2.plot(t[idt:],singletI0[j,jg,:iend],'b',lw=1)
                            ax2.plot(t[idt:],addv,':g', lw =1)
                            ax2.plot(t[idt:],vs,':k', lw =0.5)

                if savePlot and plotData:
                    pyplot.figure(bfigname)
                    pyplot.savefig(directory+'/'+bfigname+'.png',format='png',bbox_inches='tight',dpi=rdpi);
                    if v0 == -70:
                        pyplot.figure(bfign0)
                        pyplot.savefig(directory+'/'+bfign0+'.png',format='png',bbox_inches='tight',dpi=rdpi);
                    pyplot.close()
            #k
            v1 = np.repeat(singletE[i,:,idt:],ng,axis=0).reshape(ng*ng,iend)
            v2 = np.vstack([singletI0[j,:,:iend],]*ng)
            addv = v1+v2
            coef = v1*v2
            b = doubletEI[i,j,:,idt:] - addv
            arrived = False
            for it in xrange(iend):
                coefmat0 = np.vstack((coef[:,it],np.zeros((ng*ng)))).T
                coefmat1 = np.vstack((coef[:,it],np.ones((ng*ng)))).T
                if np.absolute(coef[:,it]).sum() == 0:
                    print 'inf slope due to delayed arrival of EPSP at soma'
                    kEI[i,j,idt+it] = 0
                    c = [0,0]
                else:
                    if not arrived and np.std(coef[:,it]) > 0:
                        arrived = True
                    kEI[i,j,idt+it] = np.linalg.lstsq(coefmat0,b[:,it])[0][0]
                    c = np.linalg.lstsq(coefmat1,b[:,it])[0]
                interceptEI[i,j,idt+it] = c[1]
                pred0 = kEI[i,j,idt+it]*coef[:,it]
                pred1 = c[0]*coef[:,it]+c[1]
                rsEI[i,j,idt+it] = getRS(pred0,b[:,it])

            if plotData:
                kfigname = 'k-E'+str(i)+'-I'+str(j)
                wfigname = 'midK-E'+str(i)+'-I'+str(j)
                getKfig(kfigname,t[idt:],kEI[i,j,idt:],rsEI[i,j,idt:],interceptEI[i,j,idt:],coef,b,plim,directory,savePlot,'EI')
                getWfig(wfigname,kEI[i,j,idt:],coef,b,directory,savePlot,'EI')
                wfigname = 'worstK-E'+str(i)+'-I'+str(j)
                it = np.argmin(rsEI[i,j,idt:])
                getWfig(wfigname,kEI[i,j,idt:],coef,b,directory,savePlot,'EI',it)
            if v0 == -70:
                v1 = np.repeat(singletE0[i,:,idt:],ng,axis=0).reshape(ng*ng,iend)
                v2 = np.vstack([singletI0[j,:,:iend],]*ng)
                addv = v1+v2
                coef = v1*v2
                b = doubletEI0[i,j,:,idt:] - addv
                arrived = False
                for it in xrange(iend):
                    coefmat0 = np.vstack((coef[:,it],np.zeros((ng*ng)))).T
                    coefmat1 = np.vstack((coef[:,it],np.ones((ng*ng)))).T
                    if np.absolute(coef[:,it]).sum() == 0:
                        print 'inf slope due to delayed arrival of EPSP at soma'
                        kEI0[i,j,idt+it] = 0
                        c = [0,0]
                    else:
                        if not arrived and np.std(coef[:,it]) > 0:
                            arrived = True
                        kEI0[i,j,idt+it] = np.linalg.lstsq(coefmat0,b[:,it])[0][0]
                        c = np.linalg.lstsq(coefmat1,b[:,it])[0]
                    interceptEI0[i,j,idt+it] = c[1]
                    pred0 = kEI0[i,j,idt+it]*coef[:,it]
                    pred1 = c[0]*coef[:,it]+c[1]
                    rsEI0[i,j,idt+it] = getRS(pred0,b[:,it])

                if plotData:
                    kfign0 = 'k0-E'+str(i)+'-I'+str(j)
                    wfign0 = 'midK0-E'+str(i)+'-I'+str(j)
                    getKfig(kfign0,t[idt:],kEI0[i,j,idt:],rsEI0[i,j,idt:],interceptEI0[i,j,idt:],coef,b,plim,directory,savePlot,'EI')
                    getWfig(wfign0,kEI0[i,j,idt:],coef,b,directory,savePlot,'EI')
                    wfign0 = 'worstK0-E'+str(i)+'-I'+str(j)
                    it = np.argmin(rsEI0[i,j,idt:])
                    getWfig(wfign0,kEI0[i,j,idt:],coef,b,directory,savePlot,'EI',it)

    if v0 == -70:
        doubletEI0.dump(directory+'/doubletEI0.npy')
        kEI0.dump(directory+'/kEI0.npy')
        rsEI0.dump(directory+'/rsEI0.npy') 
        interceptEI0.dump(directory+'/interceptEI0.npy')

    doubletEI.dump(directory+'/doubletEI.npy')
    dfiredEI.dump(directory+'/dfiredEI.npy')
    kEI.dump(directory+'/kEI.npy')
    rsEI.dump(directory+'/rsEI.npy') 
    interceptEI.dump(directory+'/interceptEI.npy')
pyplot.show()