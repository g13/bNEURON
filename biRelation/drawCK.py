import numpy as np
import sys, os, getopt
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot, gridspec
#np.seterr(all='raise')
pyplot.style.use('publish_cms')
fmt='eps'
bw=False
def write_one(data_filename,data,mode='ab'):
    with open(data_filename,mode) as data_file:
        data.tofile(data_file)

def getRs_v(rs,rsE,dtRange,ndt,tstep,nloc1,nloc2,nv,run_nt,k):
    for idt in xrange(ndt):
        iidt = np.round(dtRange[idt]/tstep).astype('int')
        nt = run_nt-iidt
        tmp = rs[:,idt,:,:,iidt:].reshape((nv,nloc1*nloc2,nt))
        rs = rs + tmp.mean(-1)
        rsE = rsE + tmp.std(-1)
    rs = rs/ndt
    rsE = rsE/ndt
def getdK_v(dk,tstep,nloc1,nloc2,nv,run_nt,k,k0,typeEI,fdrname,suffix):
    #print typeEI
    ndt = 1
    idt = 0
    dt = 0
    for iv in xrange(nv):
        v = vRange[iv]
        #print iv, "========"
        iidt = np.round(dt/tstep).astype('int')
        nt = run_nt-iidt
        fdr= fdrname+'-'+str(np.absolute(np.round(v)))+'-'+str(int(dt))+'-'+suffix
        v1 = np.load(fdr+'/singlet'+typeEI[0]+'.npy')
        fdr= fdrname+'-'+str(np.absolute(np.round(v)))+'-0-'+suffix
        v2 = np.load(fdr+'/singlet'+typeEI[1]+'0.npy')
        ng = np.int(v1.shape[1])
        vv = np.empty((nloc1*nloc2,ng*ng,nt))
        for i in xrange(nloc1):
            for j in xrange(nloc2):
                v10 = np.repeat(v1[i,:,iidt:],ng,axis=0).reshape((ng*ng,nt))
                v20 = np.vstack([v2[j,:,:nt],]*ng)
                vv[i*nloc2+j,:,:] = v10*v20
        index = np.absolute(vv.sum(-2)).argmax(-1)
        tmp = np.take(k[iv,idt,:,:,iidt:].reshape((nloc1*nloc2,nt)),index)
        tmp0 = np.take(k0[idt,:,:,iidt:].reshape((nloc1*nloc2,nt)),index)
        dk[iv,:] = (tmp-tmp0)/tmp0*100

def getK_v(tstep,nloc1,nloc2,nv,run_nt,k,typeEI,fdrname,suffix):
    #print typeEI, 'getK_v'
    ndt = 1
    idt = 0
    dt = 0
    kpick = np.zeros((nv,nloc1*nloc2))
    for iv in xrange(nv):
        v = vRange[iv]
        iidt = np.round(dt/tstep).astype('int')
        nt = run_nt-iidt
        fdr= fdrname+'-'+str(np.absolute(np.round(v)))+'-'+str(int(dt))+'-'+suffix
        v1 = np.load(fdr+'/singlet'+typeEI[0]+'.npy')
        fdr= fdrname+'-'+str(np.absolute(np.round(v)))+'-0-'+suffix
        v2 = np.load(fdr+'/singlet'+typeEI[1]+'0.npy')
        ng = np.int(v1.shape[1])
        vv = np.empty((nloc1*nloc2,ng*ng,nt))
        for i in xrange(nloc1):
            for j in xrange(nloc2):
                v10 = np.repeat(v1[i,:,iidt:],ng,axis=0).reshape((ng*ng,nt))
                v20 = np.vstack([v2[j,:,:nt],]*ng)
                vv[i*nloc2+j,:,:] = v10*v20
        index = np.absolute(vv).sum(-2).argmax(-1)
        tmp = k[iv,idt,:,:,iidt:].reshape((nloc1*nloc2,nt))
        for i in xrange(nloc1*nloc2):
            kpick[iv,i] = tmp[i,index[i]]
    return kpick

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
    nloc1 = np.int(rs.shape[2])
    nloc2 = np.int(rs.shape[3])
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
            vv = np.empty((nloc1*nloc2,ng*ng,nt))
            for i in xrange(nloc1):
                for j in xrange(nloc2):
                    v10 = np.repeat(v1[i,:,iidt:],ng,axis=0).reshape(ng*ng,nt)
                    v20 = np.vstack([v2[j,:,:nt],]*ng)
                    vv[i*nloc2+j,:,:] = v10*v20
            index0 = np.absolute(vv.sum(-2)).argmax(-1)
            tmp = rs[iv,idt,:,:,:].reshape((nloc1*nloc2,run_nt))
            index = np.arange(nloc1*nloc2)*run_nt + (iidt+index0)
            rsMax[iv,idt,:] = np.take(tmp,index)
v0 = -70
idpi = 900
run_t = 240
suffix = '238660'
theme='cK'
tstep = 0.1
nlocE = 6
nlocI = 3
vRange = np.array([-74,-70,-66,-62])
dtRange = np.array([0,4,8,12,20,30,50,70,130,210],dtype='float')
idtRange = (dtRange/tstep).astype('int')
#dtRange = np.array([0,10,20,30],dtype='float')
nv = vRange.size
run_nt = np.round(run_t/tstep + 1).astype('int')
ndt = dtRange.size-2
ndtRS = 3;
ndtG = min([ndt,ndt])
print nv, ndt
kEE = np.empty((nv,ndt,nlocE,nlocE,run_nt))
rsEE = np.empty((nv,ndt,nlocE,nlocE,run_nt))
interceptEE = np.empty((nv,ndt,nlocE,nlocE,run_nt))
kEI = np.empty((nv,ndt,nlocE,nlocI,run_nt))
rsEI = np.empty((nv,ndt,nlocE,nlocI,run_nt))
interceptEI = np.empty((nv,ndt,nlocE,nlocI,run_nt))
kIE = np.empty((nv,ndt,nlocI,nlocE,run_nt))
rsIE = np.empty((nv,ndt,nlocI,nlocE,run_nt))
interceptIE = np.empty((nv,ndt,nlocI,nlocE,run_nt))
kII = np.empty((nv,ndt,nlocI,nlocI,run_nt))
rsII = np.empty((nv,ndt,nlocI,nlocI,run_nt))
interceptII = np.empty((nv,ndt,nlocI,nlocI,run_nt))
kEE0 = np.empty((ndt,nlocE,nlocE,run_nt))
rsEE0 = np.empty((ndt,nlocE,nlocE,run_nt))
kEI0 = np.empty((ndt,nlocE,nlocI,run_nt))
rsEI0 = np.empty((ndt,nlocE,nlocI,run_nt))
kIE0 = np.empty((ndt,nlocI,nlocE,run_nt))
rsIE0 = np.empty((ndt,nlocI,nlocE,run_nt))
kII0 = np.empty((ndt,nlocI,nlocI,run_nt))
rsII0 = np.empty((ndt,nlocI,nlocI,run_nt))
figname = theme
fdrname = theme
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
write_one(RsKfn,rsEE,'ab')
write_one(RsKfn,kEI,'ab')
write_one(RsKfn,rsEI,'ab')
write_one(RsKfn,kIE,'ab')
write_one(RsKfn,rsIE,'ab')
write_one(RsKfn,kII,'ab')
write_one(RsKfn,rsII,'ab')

RsKfn0 = theme+'RsK0.bin'
write_one(RsKfn0,kEE0,'wb')
write_one(RsKfn0,rsEE0,'ab')
write_one(RsKfn0,kEI0,'ab')
write_one(RsKfn0,rsEI0,'ab')
write_one(RsKfn0,kIE0,'ab')
write_one(RsKfn0,rsIE0,'ab')
write_one(RsKfn0,kII0,'ab')
write_one(RsKfn0,rsII0,'ab')

pyplot.figure(figname,figsize = (25,5)) # inches
tpickK = idtRange[2:-1]
ndtK0 = 0
ndtK1 = 1

gs = gridspec.GridSpec(1, 4)

pyplot.subplot(gs[0])

kEE_all = getK_v(tstep,nlocE,nlocE,nv,run_nt,kEE,'EE',fdrname,suffix)
kIE_all = getK_v(tstep,nlocI,nlocE,nv,run_nt,kIE,'IE',fdrname,suffix)
kEI_all = getK_v(tstep,nlocE,nlocI,nv,run_nt,kEI,'EI',fdrname,suffix)
kII_all = getK_v(tstep,nlocI,nlocI,nv,run_nt,kII,'II',fdrname,suffix)

if bw:
    lee,_,_ = pyplot.errorbar(vRange,kEE_all.mean(-1),kEE_all.std(-1),ls=':',marker='*',c='k',label='EE')
    lei,_,_ = pyplot.errorbar(vRange,kEI_all.mean(-1),kEI_all.std(-1),ls='-',marker='*',c='k',label='EI')
    lie,_,_ = pyplot.errorbar(vRange,kIE_all.mean(-1),kIE_all.std(-1),ls=':',marker='o', mfc='none',c='k',label='IE')
    lii,_,_ = pyplot.errorbar(vRange,kII_all.mean(-1),kII_all.std(-1),ls='-',marker='o', mfc='none',c='k',label='II')
else:
    lee,_,_ = pyplot.errorbar(vRange,kEE_all.mean(-1),kEE_all.std(-1),ls='-',marker='*',c='r',label='EE')
    lei,_,_ = pyplot.errorbar(vRange,kEI_all.mean(-1),kEI_all.std(-1),ls='-',marker='*',c='g',label='EI')
    lie,_,_ = pyplot.errorbar(vRange,kIE_all.mean(-1),kIE_all.std(-1),ls='-',marker='o', mfc='none',c='m',label='IE')
    lii,_,_ = pyplot.errorbar(vRange,kII_all.mean(-1),kII_all.std(-1),ls='-',marker='o', mfc='none',c='b',label='II')
pyplot.legend([lee,lei,lie,lii],['EE','EI','IE','II'])
pyplot.xlabel(r'$V_{0} (mV)$')
pyplot.ylabel(r'$k (mV^{-1})$')
pyplot.xticks(vRange)

pyplot.subplot(gs[1])

dk_EE = np.empty((nv,nlocE * nlocE))
getdK_v(dk_EE,tstep,nlocE,nlocE,nv,run_nt,kEE,kEE0,'EE',fdrname,suffix)
#for iloc in xrange(nlocE*nlocE):
#    lee,_,_ = pyplot.errorbar(vRange,dk_EE[:,iloc],dkE_EE[:,iloc],ls=':',marker='*', label='EE')
if bw:
    lee,_,_ = pyplot.errorbar(vRange,dk_EE.mean(1),dk_EE.std(1),ls=':',marker='*',c='k', label='EE')
else:
    lee,_,_ = pyplot.errorbar(vRange,dk_EE.mean(1),dk_EE.std(1),ls='-',marker='*',c='r', label='EE')

dk_EI = np.empty((nv,nlocE * nlocI))
getdK_v(dk_EI,tstep,nlocE,nlocI,nv,run_nt,kEI,kEI0,'EI',fdrname,suffix)
#for iloc in xrange(nlocE*nlocI):
#    lei,_,_ = pyplot.errorbar(vRange,dk_EI[:,iloc],dkE_EI[:,iloc],ls='-',marker='*', label='EI')
if bw:
    lei,_,_ = pyplot.errorbar(vRange,dk_EI.mean(1),dk_EI.std(1),ls='-',marker='*',c='k', label='EI')
else:
    lei,_,_ = pyplot.errorbar(vRange,dk_EI.mean(1),dk_EI.std(1),ls='-',marker='*',c='g', label='EI')

dk_IE = np.empty((nv,nlocI * nlocE))
getdK_v(dk_IE,tstep,nlocI,nlocE,nv,run_nt,kIE,kIE0,'IE',fdrname,suffix)
#for iloc in xrange(nlocI*nlocE):
#    lie,_,_ = pyplot.errorbar(vRange,dk_IE[:,iloc],dkE_IE[:,iloc],ls=':',marker='o', label='IE')
if bw:
    lie,_,_ = pyplot.errorbar(vRange,dk_IE.mean(1),dk_IE.std(1),ls=':',marker='o',c='k', mfc='none', label='IE')
else:
    lie,_,_ = pyplot.errorbar(vRange,dk_IE.mean(1),dk_IE.std(1),ls='-',marker='o',c='m', mfc='none', label='IE')

dk_II = np.empty((nv,nlocI * nlocI))
getdK_v(dk_II,tstep,nlocI,nlocI,nv,run_nt,kII,kII0,'II',fdrname,suffix)
#for iloc in xrange(nlocI*nlocI):
    #lii,_,_ = pyplot.errorbar(vRange,dk_II[:,iloc],dkE_II[:,iloc],ls='-',marker='o', label='II')
if bw:
    lii,_,_ = pyplot.errorbar(vRange,dk_II.mean(1),dk_II.std(1),ls='-',marker='o',c='k', mfc='none', label='II')
else:
    lii,_,_ = pyplot.errorbar(vRange,dk_II.mean(1),dk_II.std(1),ls='-',marker='o',c='b', mfc='none', label='II')

pyplot.xticks(vRange)
pyplot.legend([lee,lei,lie,lii],['EE','EI','IE','II'])
pyplot.ylabel('% Change')
pyplot.xlabel(r'$V_{0} (mV)$')
#pyplot.legend([lee,lei],['EE','EI'])

#pyplot.subplot(132)
pyplot.subplot(gs[2])

rsMaxEE=np.empty((nv,ndtG,nlocE*nlocE))
getArgMax(vRange,dtRange,nv,ndtG,fdrname,suffix,rsEE,rsMaxEE,'EE',tstep,run_nt)
print "rsEE-----"
rsMax = rsMaxEE.reshape((nv,ndtG*nlocE*nlocE))
if bw:
    lee,_,_ = pyplot.errorbar(vRange,rsMax.mean(-1),rsMax.std(-1),ls=':',marker='*',c='k',label='EE')
else:
    lee,_,_ = pyplot.errorbar(vRange,rsMax.mean(-1),rsMax.std(-1),ls='-',marker='*',c='r',label='EE')

rsMaxEI=np.empty((nv,ndtG,nlocE*nlocI))
getArgMax(vRange,dtRange,nv,ndtG,fdrname,suffix,rsEI,rsMaxEI,'EI',tstep,run_nt)
print "rsEI-----"
rsMax = rsMaxEI.reshape((nv,ndtG*nlocE*nlocI))
if bw:
    lei,_,_ = pyplot.errorbar(vRange,rsMax.mean(-1),rsMax.std(-1),ls='-',marker='*',c='k',label='EI')
else:
    lei,_,_ = pyplot.errorbar(vRange,rsMax.mean(-1),rsMax.std(-1),ls='-',marker='*',c='g',label='EI')

rsMaxIE=np.empty((nv,ndtG,nlocI*nlocE))
getArgMax(vRange,dtRange,nv,ndtG,fdrname,suffix,rsIE,rsMaxIE,'IE',tstep,run_nt)
print "rsIE-----"
rsMax = rsMaxIE.reshape((nv,ndtG*nlocI*nlocE))
if bw:
    lie,_,_ = pyplot.errorbar(vRange,rsMax.mean(-1),rsMax.std(-1),ls=':',marker='o', mfc='none',c='k',label='IE')
else:
    lie,_,_ = pyplot.errorbar(vRange,rsMax.mean(-1),rsMax.std(-1),ls='-',marker='o', mfc='none',c='m',label='IE')

rsMaxII=np.empty((nv,ndtG,nlocI*nlocI))
getArgMax(vRange,dtRange,nv,ndtG,fdrname,suffix,rsII,rsMaxII,'II',tstep,run_nt)
print "rsII-----"
rsMax = rsMaxII.reshape((nv,ndtG*nlocI*nlocI))
if bw:
    lii,_,_ = pyplot.errorbar(vRange,rsMax.mean(-1),rsMax.std(-1),ls='-',marker='o', mfc='none',c='k',label='II')
else:
    lii,_,_ = pyplot.errorbar(vRange,rsMax.mean(-1),rsMax.std(-1),ls='-',marker='o', mfc='none',c='b',label='II')

pyplot.legend([lee,lei,lie,lii],['EE','EI','IE','II'])
pyplot.xlabel(r'$V_{0} (mV)$')
pyplot.ylabel(r'$R^{2}$')
pyplot.xticks(vRange)
pyplot.ylim(0.99,1.01)
pyplot.yticks([0.99,0.995,1.0,1.005,1.01])

#pyplot.subplot(133)
pyplot.subplot(gs[3])

interceptMaxEE=np.empty((nv,ndtG,nlocE*nlocE))
getArgMax(vRange,dtRange,nv,ndtG,fdrname,suffix,interceptEE,interceptMaxEE,'EE',tstep,run_nt)
print "interceptEE-----"
interceptMax = interceptMaxEE.reshape((nv,ndtG*nlocE*nlocE))
if bw:
    lee,_,_ = pyplot.errorbar(vRange,interceptMax.mean(-1),interceptMax.std(-1),ls=':',marker='*',c='k',label='EE')
else:
    lee,_,_ = pyplot.errorbar(vRange,interceptMax.mean(-1),interceptMax.std(-1),ls='-',marker='*',c='r',label='EE')

interceptMaxEI=np.empty((nv,ndtG,nlocE*nlocI))
getArgMax(vRange,dtRange,nv,ndtG,fdrname,suffix,interceptEI,interceptMaxEI,'EI',tstep,run_nt)
print "interceptEI-----"
interceptMax = interceptMaxEI.reshape((nv,ndtG*nlocE*nlocI))
if bw:
    lei,_,_ = pyplot.errorbar(vRange,interceptMax.mean(-1),interceptMax.std(-1),ls='-',marker='*',c='k',label='EI')
else:
    lei,_,_ = pyplot.errorbar(vRange,interceptMax.mean(-1),interceptMax.std(-1),ls='-',marker='*',c='g',label='EI')

interceptMaxIE=np.empty((nv,ndtG,nlocI*nlocE))
getArgMax(vRange,dtRange,nv,ndtG,fdrname,suffix,interceptIE,interceptMaxIE,'IE',tstep,run_nt)
print "interceptIE-----"
interceptMax = interceptMaxIE.reshape((nv,ndtG*nlocI*nlocE))
if bw:
    lie,_,_ = pyplot.errorbar(vRange,interceptMax.mean(-1),interceptMax.std(-1),ls=':',marker='o', mfc='none',c='k',label='IE')
else:
    lie,_,_ = pyplot.errorbar(vRange,interceptMax.mean(-1),interceptMax.std(-1),ls='-',marker='o', mfc='none',c='m',label='IE')

interceptMaxII=np.empty((nv,ndtG,nlocI*nlocI))
getArgMax(vRange,dtRange,nv,ndtG,fdrname,suffix,interceptII,interceptMaxII,'II',tstep,run_nt)
print "interceptII-----"
interceptMax = interceptMaxII.reshape((nv,ndtG*nlocI*nlocI))
if bw:
    lii,_,_ = pyplot.errorbar(vRange,interceptMax.mean(-1),interceptMax.std(-1),ls='-',marker='o', mfc='none',c='k',label='II')
else:
    lii,_,_ = pyplot.errorbar(vRange,interceptMax.mean(-1),interceptMax.std(-1),ls='-',marker='o', mfc='none',c='b',label='II')
pyplot.xticks(vRange)

pyplot.legend([lee,lei,lie,lii],['EE','EI','IE','II'])
pyplot.xlabel(r'$V_{0} (mV)$')
pyplot.ylabel(r'$intercept (mV)$')

pyplot.tight_layout()
pyplot.savefig(theme+'.'+fmt,format=fmt,dpi=idpi)

pyplot.figure(figname+'-t',figsize = (25,5)) # inches

gs = gridspec.GridSpec(1, 4)

pyplot.subplot(gs[0])

#kEE_all = getK_v(tstep,nlocE,nlocE,nv,run_nt,kEE,'EE',fdrname,suffix)
#kIE_all = getK_v(tstep,nlocI,nlocE,nv,run_nt,kIE,'IE',fdrname,suffix)
#kEI_all = getK_v(tstep,nlocE,nlocI,nv,run_nt,kEI,'EI',fdrname,suffix)
#kII_all = getK_v(tstep,nlocI,nlocI,nv,run_nt,kII,'II',fdrname,suffix)

if bw:
    lee,_,_ = pyplot.errorbar(vRange,kEE_all.mean(-1),kEE_all.std(-1),ls=':',marker='*',c='k',label='EE')
    lei,_,_ = pyplot.errorbar(vRange,kEI_all.mean(-1),kEI_all.std(-1),ls='-',marker='*',c='k',label='EI')
    lie,_,_ = pyplot.errorbar(vRange,kIE_all.mean(-1),kIE_all.std(-1),ls=':',marker='o', mfc='none',c='k',label='IE')
    lii,_,_ = pyplot.errorbar(vRange,kII_all.mean(-1),kII_all.std(-1),ls='-',marker='o', mfc='none',c='k',label='II')
else:
    lee,_,_ = pyplot.errorbar(vRange,kEE_all.mean(-1),kEE_all.std(-1),ls='-',marker='*',c='r',label='EE')
    lei,_,_ = pyplot.errorbar(vRange,kEI_all.mean(-1),kEI_all.std(-1),ls='-',marker='*',c='g',label='EI')
    lie,_,_ = pyplot.errorbar(vRange,kIE_all.mean(-1),kIE_all.std(-1),ls='-',marker='o', mfc='none',c='m',label='IE')
    lii,_,_ = pyplot.errorbar(vRange,kII_all.mean(-1),kII_all.std(-1),ls='-',marker='o', mfc='none',c='b',label='II')

pyplot.legend([lee,lei,lie,lii],['EE','EI','IE','II'])
pyplot.xlabel(r'$V_{0} (mV)$')
pyplot.ylabel(r'$k (mV^{-1})$')
pyplot.xticks(vRange)

pyplot.subplot(gs[1])

if bw:
    lee,_,_ = pyplot.errorbar(vRange,dk_EE.mean(-1),dk_EE.std(-1),ls=':',marker='*',c='k',label='EE')
    lei,_,_ = pyplot.errorbar(vRange,dk_EI.mean(-1),dk_EI.std(-1),ls='-',marker='*',c='k',label='EI')
    lie,_,_ = pyplot.errorbar(vRange,dk_IE.mean(-1),dk_IE.std(-1),ls=':',marker='o', mfc='none',c='k',label='IE')
    lii,_,_ = pyplot.errorbar(vRange,dk_II.mean(-1),dk_II.std(-1),ls='-',marker='o', mfc='none',c='k',label='II')
else:
    lee,_,_ = pyplot.errorbar(vRange,dk_EE.mean(-1),dk_EE.std(-1),ls='-',marker='*',c='r',label='EE')
    lei,_,_ = pyplot.errorbar(vRange,dk_EI.mean(-1),dk_EI.std(-1),ls='-',marker='*',c='g',label='EI')
    lie,_,_ = pyplot.errorbar(vRange,dk_IE.mean(-1),dk_IE.std(-1),ls='-',marker='o', mfc='none',c='m',label='IE')
    lii,_,_ = pyplot.errorbar(vRange,dk_II.mean(-1),dk_II.std(-1),ls='-',marker='o', mfc='none',c='b',label='II')

pyplot.legend([lee,lei,lie,lii],['EE','EI','IE','II'])
pyplot.xlabel(r'$V_{0} (mV)$')
pyplot.ylabel(r'$k (mV^{-1})$')
pyplot.xticks(vRange)

tpickRS = idtRange[2:-1]
tpickIntercept = idtRange[2:-1]
ndtIntercept0 = 0 
ndtIntercept1 = 3
ndtRS0 = 0
ndtRS1 = 3

pyplot.subplot(gs[2])

rsEE_all = rsEE[:,ndtRS0:ndtRS1,:,:,tpickRS].mean(1).mean(1).mean(1)
rsEI_all = rsEI[:,ndtRS0:ndtRS1,:,:,tpickRS].mean(1).mean(1).mean(1)
rsIE_all = rsIE[:,ndtRS0:ndtRS1,:,:,tpickRS].mean(1).mean(1).mean(1)
rsII_all = rsII[:,ndtRS0:ndtRS1,:,:,tpickRS].mean(1).mean(1).mean(1)

if bw:
    lee,_,_ = pyplot.errorbar(vRange,rsEE_all.mean(-1),rsEE_all.std(-1),ls=':',marker='*',c='k',label='EE')
    lei,_,_ = pyplot.errorbar(vRange,rsEI_all.mean(-1),rsEI_all.std(-1),ls='-',marker='*',c='k',label='EI')
    lie,_,_ = pyplot.errorbar(vRange,rsIE_all.mean(-1),rsIE_all.std(-1),ls=':',marker='o', mfc='none',c='k',label='IE')
    lii,_,_ = pyplot.errorbar(vRange,rsII_all.mean(-1),rsII_all.std(-1),ls='-',marker='o', mfc='none',c='k',label='II')
else:
    lee,_,_ = pyplot.errorbar(vRange,rsEE_all.mean(-1),rsEE_all.std(-1),ls='-',marker='*',c='r',label='EE')
    lei,_,_ = pyplot.errorbar(vRange,rsEI_all.mean(-1),rsEI_all.std(-1),ls='-',marker='*',c='g',label='EI')
    lie,_,_ = pyplot.errorbar(vRange,rsIE_all.mean(-1),rsIE_all.std(-1),ls='-',marker='o', mfc='none',c='m',label='IE')
    lii,_,_ = pyplot.errorbar(vRange,rsII_all.mean(-1),rsII_all.std(-1),ls='-',marker='o', mfc='none',c='b',label='II')

pyplot.legend([lee,lei,lie,lii],['EE','EI','IE','II'])
pyplot.xlabel(r'$V_{0} (mV)$')
pyplot.ylabel(r'$R^{2}$')
pyplot.xticks(vRange)
pyplot.ylim(0.99,1.01)
pyplot.yticks([0.99,0.995,1.0,1.005,1.01])

pyplot.subplot(gs[3])
interceptEE_all = interceptEE[:,ndtIntercept0:ndtIntercept1,:,:,tpickIntercept].mean(1).mean(1).mean(1)
interceptEI_all = interceptEI[:,ndtIntercept0:ndtIntercept1,:,:,tpickIntercept].mean(1).mean(1).mean(1)
interceptIE_all = interceptIE[:,ndtIntercept0:ndtIntercept1,:,:,tpickIntercept].mean(1).mean(1).mean(1)
interceptII_all = interceptII[:,ndtIntercept0:ndtIntercept1,:,:,tpickIntercept].mean(1).mean(1).mean(1)

if bw:
    lee,_,_ = pyplot.errorbar(vRange,interceptEE_all.mean(-1),interceptEE_all.std(-1),ls=':',marker='*',c='k',label='EE')
    lei,_,_ = pyplot.errorbar(vRange,interceptEI_all.mean(-1),interceptEI_all.std(-1),ls='-',marker='*',c='k',label='EI')
    lie,_,_ = pyplot.errorbar(vRange,interceptIE_all.mean(-1),interceptIE_all.std(-1),ls=':',marker='o', mfc='none',c='k',label='IE')
    lii,_,_ = pyplot.errorbar(vRange,interceptII_all.mean(-1),interceptII_all.std(-1),ls='-',marker='o', mfc='none',c='k',label='II')
else:
    lee,_,_ = pyplot.errorbar(vRange,interceptEE_all.mean(-1),interceptEE_all.std(-1),ls='-',marker='*',c='r',label='EE')
    lei,_,_ = pyplot.errorbar(vRange,interceptEI_all.mean(-1),interceptEI_all.std(-1),ls='-',marker='*',c='g',label='EI')
    lie,_,_ = pyplot.errorbar(vRange,interceptIE_all.mean(-1),interceptIE_all.std(-1),ls='-',marker='o', mfc='none',c='m',label='IE')
    lii,_,_ = pyplot.errorbar(vRange,interceptII_all.mean(-1),interceptII_all.std(-1),ls='-',marker='o', mfc='none',c='b',label='II')


pyplot.legend([lee,lei,lie,lii],['EE','EI','IE','II'])
pyplot.xlabel(r'$V_{0} (mV)$')
pyplot.ylabel(r'$intercept (mV)$')
pyplot.xticks(vRange)

pyplot.tight_layout()
pyplot.savefig(theme+'-t.'+fmt,format=fmt,dpi=idpi)

