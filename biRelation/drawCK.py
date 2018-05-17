import numpy as np
import sys, os, getopt
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot, gridspec
#np.seterr(all='raise')
pyplot.style.use('thin_line')
iv0 = 4
def write_one(data_filename,data,mode='ab'):
    with open(data_filename,mode) as data_file:
        data.tofile(data_file)

def getRs_v(rs,rsE,dtRange,ndt,tstep,nloc,run_nt,k):
    for idt in xrange(ndt):
        iidt = np.round(dtRange[idt]/tstep).astype('int')
        nt = run_nt-iidt
        tmp = rs[:,idt,:,:,iidt:].reshape((nv,nloc*nloc,nt))
        rs = rs + tmp.mean(-1)
        rsE = rsE + tmp.std(-1)
    rs = rs/ndt
    rsE = rsE/ndt
def getdK_v(dk,dkE,dtRange,ndt,tstep,nloc,run_nt,k,k0,typeEI,fdrname,suffix):
    #print typeEI
    for iv in xrange(nv):
        v = vRange[iv]
        dtTmp = np.zeros((ndt,nloc*nloc))
        #print iv, "========"
        for idt in xrange(ndt):
            dt = dtRange[idt]
            iidt = np.round(dt/tstep).astype('int')
            nt = run_nt-iidt
            fdr= fdrname+'-'+str(np.absolute(np.round(v)))+'-'+str(int(dt))+'-'+suffix
            v1 = np.load(fdr+'/singlet'+typeEI[0]+'.npy')
            fdr= fdrname+'-'+str(np.absolute(np.round(v)))+'-0-'+suffix
            v2 = np.load(fdr+'/singlet'+typeEI[1]+'0.npy')
            ng = np.int(v1.shape[1])
            vv = np.empty((nloc*nloc,ng*ng,nt))
            for i in xrange(nloc):
                for j in xrange(nloc):
                    v10 = np.repeat(v1[i,:,iidt:],ng,axis=0).reshape((ng*ng,nt))
                    v20 = np.vstack([v2[j,:,:nt],]*ng)
                    vv[i*nloc+j,:,:] = v10*v20
            index = np.absolute(vv.sum(-2)).argmax(-1)
            tmp = np.take(k[iv,idt,:,:,iidt:].reshape((nloc*nloc,nt)),index)
            tmp0 = np.take(k0[idt,:,:,iidt:].reshape((nloc*nloc,nt)),index)
            dtTmp[idt,:] = (tmp-tmp0)/tmp0*100
            #print np.absolute(dtTmp[idt,:]).max()
            #print "------------"

        dk[iv,:] =  dtTmp.mean(0)
        dkE[iv,:] = dtTmp.std(0)
        #print np.absolute(dk[iv,:]).max()
        #print np.absolute(dkE[iv,:]).max()
        #print "============="
def getRsK(vRange,dtRange,nv,ndt,fdrname,suffix,k,R,typeEI,tstep):
    for iv in xrange(nv):
        for idt in xrange(ndt):
            v = vRange[iv]
            dt = dtRange[idt]
            iidt = np.int(np.round(dt/tstep))
            fdr=fdrname+'-'+str(np.absolute(np.round(v)))+'-'+str(int(dt))+'-'+suffix
            k[iv,idt,:,:,:] = np.load(fdr+'/k'+typeEI+'.npy')
            R[iv,idt,:,:,:] = np.load(fdr+'/rs'+typeEI+'.npy')
    k[np.isnan(k)] = 0
    k[np.isinf(k)] = 0
    R[np.isnan(R)] = 1
    R[np.isinf(R)] = 1
def getRsK0(dtRange,ndt,fdrname,suffix,k0,R0,typeEI,v0,tstep):
    for idt in xrange(ndt):
        dt = dtRange[idt]
        iidt = np.int(np.round(dt/tstep))
        fdr=fdrname+'-'+str(np.absolute(np.round(v0)))+'-'+str(int(dt))+'-'+suffix
        k0[idt,:,:,:] = np.load(fdr+'/k'+typeEI+'0.npy')
        R0[idt,:,:,:] = np.load(fdr+'/rs'+typeEI+'0.npy')
    k0[np.isnan(k0)] = 0
    k0[np.isinf(k0)] = 0
    R0[np.isnan(R0)] = 1
    R0[np.isinf(R0)] = 1
def getArgMax(vRange,dtRange,nv,ndt,fdrname,suffix,rs,rsMax,typeEI,tstep,run_nt):
    nloc = np.int(rs.shape[2])
    for iv in xrange(nv):
        v = vRange[iv]
        for idt in xrange(ndt):
            dt = dtRange[idt]
            iidt = np.int(np.round(dt/tstep))
            nt = run_nt - iidt
            fdr=fdrname+'-'+str(np.absolute(np.round(v)))+'-'+str(int(dt))+'-'+suffix
            v1 = np.load(fdr+'/singlet'+typeEI[0]+'.npy')
            fdr=fdrname+'-'+str(np.absolute(np.round(v)))+'-0-'+suffix
            v2 = np.load(fdr+'/singlet'+typeEI[1]+'0.npy')
            ng = np.int(v1.shape[1])
            vv = np.empty((nloc*nloc,ng*ng,nt))
            for i in xrange(nloc):
                for j in xrange(nloc):
                    v10 = np.repeat(v1[i,:,iidt:],ng,axis=0).reshape(ng*ng,nt)
                    v20 = np.vstack([v2[j,:,:nt],]*ng)
                    vv[i*nloc+j,:,:] = v10*v20
            index0 = np.absolute(vv.sum(-2)).argmax(-1)
            tmp = rs[iv,idt,:,:,:].reshape((nloc*nloc,run_nt))
            index = np.arange(nloc*nloc)*run_nt + (iidt+index0)
            rsMax[iv,idt,:] = np.take(tmp,index)
v0 = -70
idpi = 600
run_t = 240
suffix = '238660'
theme='cK'
tstep = 0.1
nloc = 4
vRange = np.arange(-74,4,-66)
#dtRange = np.array([0,4,8,12,20,30,70,130,210],dtype='float')
dtRange = np.array([0],dtype='float')
nv = vRange.size
run_nt = np.round(run_t/tstep + 1).astype('int')
ndt = dtRange.size
ndtG = min([2,ndt])
ndtD = min([2,ndt]) 
print nv, ndt
kEE = np.empty((nv,ndt,nloc,nloc,run_nt))
rsEE = np.empty((nv,ndt,nloc,nloc,run_nt))
kEI = np.empty((nv,ndt,nloc,nloc,run_nt))
rsEI = np.empty((nv,ndt,nloc,nloc,run_nt))
kIE = np.empty((nv,ndt,nloc,nloc,run_nt))
rsIE = np.empty((nv,ndt,nloc,nloc,run_nt))
kII = np.empty((nv,ndt,nloc,nloc,run_nt))
rsII = np.empty((nv,ndt,nloc,nloc,run_nt))
kEE0 = np.empty((ndt,nloc,nloc,run_nt))
rsEE0 = np.empty((ndt,nloc,nloc,run_nt))
kEI0 = np.empty((ndt,nloc,nloc,run_nt))
rsEI0 = np.empty((ndt,nloc,nloc,run_nt))
kIE0 = np.empty((ndt,nloc,nloc,run_nt))
rsIE0 = np.empty((ndt,nloc,nloc,run_nt))
kII0 = np.empty((ndt,nloc,nloc,run_nt))
rsII0 = np.empty((ndt,nloc,nloc,run_nt))
figname = theme
fdrname = theme
pyplot.figure(figname,figsize = (8,4), dpi=idpi) # inches
getRsK(vRange,dtRange,nv,ndt,fdrname,suffix,kEE,rsEE,'EE',tstep)
getRsK(vRange,dtRange,nv,ndt,fdrname,suffix,kEI,rsEI,'EI',tstep)
getRsK(vRange,dtRange,nv,ndt,fdrname,suffix,kIE,rsIE,'IE',tstep)
getRsK(vRange,dtRange,nv,ndt,fdrname,suffix,kII,rsII,'II',tstep)
getRsK0(dtRange,ndt,fdrname,suffix,kEE0,rsEE0,'EE',v0,tstep)
getRsK0(dtRange,ndt,fdrname,suffix,kEI0,rsEI0,'EI',v0,tstep)
getRsK0(dtRange,ndt,fdrname,suffix,kIE0,rsIE0,'IE',v0,tstep)
getRsK0(dtRange,ndt,fdrname,suffix,kII0,rsII0,'II',v0,tstep)
    
RsKfn = theme+'RsK.bin'
write_one(RsKfn,kEE,'wb')
write_one(RsKfn,rsEE,'wb')
write_one(RsKfn,kEI,'ab')
write_one(RsKfn,rsEI,'ab')
write_one(RsKfn,kIE,'ab')
write_one(RsKfn,rsIE,'ab')
write_one(RsKfn,kII,'ab')
write_one(RsKfn,rsII,'ab')

RsKfn = theme+'RsK0.bin'
write_one(RsKfn,kEE0,'wb')
write_one(RsKfn,rsEE0,'wb')
write_one(RsKfn,kEI0,'ab')
write_one(RsKfn,rsEI0,'ab')
write_one(RsKfn,kIE0,'ab')
write_one(RsKfn,rsIE0,'ab')
write_one(RsKfn,kII0,'ab')
write_one(RsKfn,rsII0,'ab')

gs = gridspec.GridSpec(10, 10)
pyplot.subplot(gs[:4,:5])

dk = np.empty((nv,nloc * nloc))
dkE = np.empty((nv,nloc * nloc))
getdK_v(dk,dkE,dtRange,ndtD,tstep,nloc,run_nt,kEE,kEE0,'EE',fdrname,suffix)
for iloc in xrange(nloc*nloc):
    lee,_,_ = pyplot.errorbar(vRange,dk[:,iloc],dkE[:,iloc],c='r', label='EE')

getdK_v(dk,dkE,dtRange,ndtD,tstep,nloc,run_nt,kEI,kEI0,'EI',fdrname,suffix)
for iloc in xrange(nloc*nloc):
    lei,_,_ = pyplot.errorbar(vRange,dk[:,iloc],dkE[:,iloc],c='m', label='EI')

getdK_v(dk,dkE,dtRange,ndtD,tstep,nloc,run_nt,kIE,kIE0,'IE',fdrname,suffix)
for iloc in xrange(nloc*nloc):
    lie,_,_ = pyplot.errorbar(vRange,dk[:,iloc],dkE[:,iloc],c='c', label='IE')

getdK_v(dk,dkE,dtRange,ndtD,tstep,nloc,run_nt,kII,kII0,'II',fdrname,suffix)
for iloc in xrange(nloc*nloc):
    lii,_,_ = pyplot.errorbar(vRange,dk[:,iloc],dkE[:,iloc],c='b', label='II')
pyplot.legend([lee,lei,lie,lii],['EE','EI','IE','II'])
pyplot.ylabel('% change')
#pyplot.legend([lee,lei],['EE','EI'])
pyplot.axis('tight')

pyplot.subplot(gs[6:,:5])
rsMax=np.empty((nv,ndtG,nloc*nloc))
getArgMax(vRange,dtRange,nv,ndtG,fdrname,suffix,rsEE,rsMax,'EE',tstep,run_nt)
print "EE-----"
if ndtG > 1:
    for idt in xrange(ndtG):
        print rsMax[0,idt,:].mean(-1)
        print rsMax[0,idt,:].std(-1)
rsMax = rsMax.reshape((nv,ndtG*nloc*nloc))
lee,_,_ = pyplot.errorbar(vRange,rsMax.mean(-1),rsMax.std(-1),c='r',label='EE')

rsMax=np.empty((nv,ndtG,nloc*nloc))
getArgMax(vRange,dtRange,nv,ndtG,fdrname,suffix,rsEI,rsMax,'EI',tstep,run_nt)
print "EI-----"
if ndtG > 1:
    for idt in xrange(ndtG):
        print rsMax[0,idt,:].mean(-1)
        print rsMax[0,idt,:].std(-1)
rsMax = rsMax.reshape((nv,ndtG*nloc*nloc))
lei,_,_ = pyplot.errorbar(vRange,rsMax.mean(-1),rsMax.std(-1),c='m',label='EI')

rsMax=np.empty((nv,ndtG,nloc*nloc))
getArgMax(vRange,dtRange,nv,ndtG,fdrname,suffix,rsIE,rsMax,'IE',tstep,run_nt)
print "IE-----"
#for idt in xrange(ndtG):
#    print rsMax[0,idt,:].mean(-1)
#    print rsMax[0,idt,:].std(-1)
rsMax = rsMax.reshape((nv,ndtG*nloc*nloc))
lie,_,_ = pyplot.errorbar(vRange,rsMax.mean(-1),rsMax.std(-1),c='c',label='IE')

rsMax=np.empty((nv,ndtG,nloc*nloc))
getArgMax(vRange,dtRange,nv,ndtG,fdrname,suffix,rsII,rsMax,'II',tstep,run_nt)
print "II-----"
#for idt in xrange(ndtG):
#    print rsMax[0,idt,:].mean(-1)
#    print rsMax[0,idt,:].std(-1)
rsMax = rsMax.reshape((nv,ndtG*nloc*nloc))
lii,_,_ = pyplot.errorbar(vRange,rsMax.mean(-1),rsMax.std(-1),c='b',label='II')
pyplot.legend([lee,lei,lie,lii],['EE','EI','IE','II'])

pyplot.savefig(theme+'.png',format='png',bbox_inches='tight',dpi=idpi)
