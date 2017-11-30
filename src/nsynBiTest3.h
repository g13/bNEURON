#ifndef NB_H
#define NB_H
#include <cmath>
#include <vector>
#include "typedefs.h"
#include "nNeuroSt.h"
#include "nNeuroLib.h"
namespace sB{
    const bool debug = true;
    const bool debug2 = true;
}
template<typename T>
inline void getNear(T *range, size n, double target, double &ratio, size &istart, size &jnext) {
    size i;
    if (target < range[0]) {
        istart = 1; jnext = 0;
        ratio = (target-range[1])/(range[0]-range[1]);
    } else if (target >= range[n-1]) {
        istart = n-2; jnext = n-1;
        ratio = (target-range[istart])/(range[jnext]-range[istart]);
    } else {
        for (i=0;i<n-1;i++) {
            if (target >= range[i] && target < range[i+1] ) {
                istart = i; jnext = i+1;
                ratio = (target-range[istart])/(range[jnext]-range[istart]);
                break;
            }
        }
    }
}

inline void interpPSP(std::vector<double> &v, size vs, double ****PSP, double *vRange, size *idtRange, size iSyn, int nv, size ndt, double vTar, double dtTar, size limit) {
    size iv, jv, idt, jdt;
    size k;
    double rv, rdt, base;
    getNear(idtRange,ndt,dtTar,rdt, idt, jdt);
    getNear(vRange, nv, vTar,rv, iv, jv);
    size iidt = idtRange[idt];
    size jjdt = idtRange[jdt];
    
    size tl = limit - jjdt;
    //cout << " jjdt " << jjdt << " tlimit " << limit << endl;
    //cout << " iv" << iv << " iSyn " << iSyn << endl;
    //cout << " jv" << jv << endl;
    //cout << " idt" << idt << endl;
    //cout << " jdt" << jdt << endl;
    for (k=1;k<tl;k++) {
        //cout << vs + k << ", "  << v[vs+k], " -> ";
        base = PSP[iv][idt][iSyn][iidt+k];
        v[vs+k] += base + rv*(PSP[jv][idt][iSyn][iidt+k]-base) + rdt*(PSP[iv][jdt][iSyn][jjdt+k]-base);
        if (idt>ndt) {
            assert(idt<ndt);
        }
        if (vs+k>v.size()) {
            //cout << "capacity: " << v.size() << endl;
            //cout << vs+k << " > " << v.size() << endl;
            //cout << " length:"  << limit << "-" << jjdt << " = " << tl << endl;
            assert(vs+k<v.size());
            //cout << v[vs+k] << endl;
        }
    }
}

inline void interpPSP0(std::vector<double> &v, size vs, double ****PSP, double *vRange, int nv, double vTar, size tl, int iSyn) {
    size i,j,k;
    double r, base;
    getNear(vRange,nv,vTar,r,i,j);
    
    for (k=1;k<tl;k++) {
        base = PSP[i][0][iSyn][k];
        v[vs+k] += base + r*(PSP[j][0][iSyn][k]-base);
    }
}

inline void interpkV_old(std::vector<double> &v, size vs, double ******kV, double *vRange, size *idtRange, int **nv, size ndt, double vTar, double dtTar0, double dtTar1, size tl, size iSyn, size jSyn) {
    size iv ,jv;
    size idt_,jdt_;
    size idt,jdt;
    size idt0,jdt0;
    size idt1,jdt1;
    size k;
    double rdt,rdt0,rdt1,rv,rdt_;
    double base0,base,dt0,dt,dv0,dv;
    if (sB::debug2) {
        cout << "   dtTar1 " << dtTar1 << endl;
        cout << "   dtTar " << dtTar0 << endl;
    }
    getNear(idtRange,ndt,dtTar0,rdt,idt,jdt);
    if (sB::debug2) {
        cout << "   idt " << idt << ", jdt" << jdt << endl;
    }
    double dtTar = dtTar1 - dtTar0;
    getNear(idtRange,ndt,dtTar,rdt_,idt_,jdt_);
    if (sB::debug2) {
        cout << "   jSyn " << jSyn << " idt_" << idt_  << endl;
    }
    getNear(vRange,nv[jSyn][idt_],vTar,rv,iv,jv);
    if (sB::debug2) {
        cout << "   nv " << nv[jSyn][idt_] << " iv " << iv  << ", jv " << jv << ", rv " << rv << endl;
    }
    
    
    if (dtTar != 0) {
        getNear(idtRange,ndt,idtRange[idt] + dtTar,rdt0,idt0,jdt0);
        getNear(idtRange,ndt,idtRange[jdt] + dtTar,rdt1,idt1,jdt1);
        size iidt0 = idtRange[idt0];
        size jjdt0 = idtRange[jdt0];
        size iidt1 = idtRange[idt1];
        size jjdt1 = idtRange[jdt1];

        for (k=1;k<tl;k++) {
            base0 = kV[iv][idt][iSyn][jSyn][idt0][iidt0+k];
            base = base0 + rdt0*(kV[iv][idt][iSyn][jSyn][jdt0][jjdt0+k] - base0);

            dv0 = kV[jv][idt][iSyn][jSyn][idt0][iidt0+k];
            dv = dv0 + rdt0*(kV[jv][idt][iSyn][jSyn][jdt0][jjdt0+k]-dv0);
            
            dt0 = kV[iv][jdt][iSyn][jSyn][idt1][iidt1+k];
            dt = dt0 + rdt1*(kV[iv][jdt][iSyn][jSyn][jdt1][jjdt1+k]-dt0);
            v[vs+k] += base + rv*(dv-base) + rdt*(dt-base);
        }
    } else {
        size iidt = idtRange[idt];
        size jjdt = idtRange[jdt];
        for (k=1;k<tl;k++) {
            base = kV[iv][idt][iSyn][jSyn][idt][iidt+k];
            //cout << " base " << base;
            v[vs+k] += base + rv  * (kV[jv][idt][iSyn][jSyn][idt][iidt+k]-base)
                            + rdt * (kV[iv][jdt][iSyn][jSyn][jdt][jjdt+k]-base);
        }
    }
}

inline void interpkV(std::vector<double> &v, size vs, double ******kV, double *vRange, size *idtRange, int nv, size ndt, double vTar, double dtTar0, double dtTar1, size tl, size iSyn, size jSyn) {
    size iv ,jv;
    size idt_,jdt_;
    size idt,jdt;
    size idt0,jdt0;
    size idt1,jdt1;
    size k;
    double rdt,rdt0,rdt1,rv,rdt_;
    double base0,base,dt0,dt,dv0,dv;
    getNear(idtRange,ndt,dtTar0,rdt,idt,jdt);
    double dtTar = dtTar1 - dtTar0;
    getNear(vRange,nv,vTar,rv,iv,jv);
    
    if (dtTar != 0) {
        getNear(idtRange,ndt,idtRange[idt] + dtTar,rdt0,idt0,jdt0);
        getNear(idtRange,ndt,idtRange[jdt] + dtTar,rdt1,idt1,jdt1);
        size iidt0 = idtRange[idt0];
        size jjdt0 = idtRange[jdt0];
        size iidt1 = idtRange[idt1];
        size jjdt1 = idtRange[jdt1];

        for (k=1;k<tl;k++) {
            base0 = kV[iv][idt][iSyn][jSyn][idt0][iidt0+k];
            base = base0 + rdt0*(kV[iv][idt][iSyn][jSyn][jdt0][jjdt0+k] - base0);

            dv0 = kV[jv][idt][iSyn][jSyn][idt0][iidt0+k];
            dv = dv0 + rdt0*(kV[jv][idt][iSyn][jSyn][jdt0][jjdt0+k]-dv0);
            
            dt0 = kV[iv][jdt][iSyn][jSyn][idt1][iidt1+k];
            dt = dt0 + rdt1*(kV[iv][jdt][iSyn][jSyn][jdt1][jjdt1+k]-dt0);
            v[vs+k] += base + rv*(dv-base) + rdt*(dt-base);
        }
    } else {
        size iidt = idtRange[idt];
        size jjdt = idtRange[jdt];
        for (k=1;k<tl;k++) {
            base = kV[iv][idt][iSyn][jSyn][idt][iidt+k];
            //cout << " base " << base;
            v[vs+k] += base + rv  * (kV[jv][idt][iSyn][jSyn][idt][iidt+k]-base)
                            + rdt * (kV[iv][jdt][iSyn][jSyn][jdt][jjdt+k]-base);
        }
    }
}

inline void interpPSP_nV(std::vector<double> &v, size vs, double ****PSP, size iSyn, double dtTar, size limit, size iv0) {
    size k;
    size idtTar = static_cast<size>(dtTar);
    for (k=1;k<limit;k++) {
        v[vs+k] += PSP[iv0][0][iSyn][idtTar+k];
    }
}

inline void interpkV0(std::vector<double> &v, size vs, double ****kV0, size *idtRange, size ndt, double dtTar0, double dtTar1, size tl, size iSyn, size jSyn) {
    size idt,jdt;
    size k;
    double rdt;
    double base;
    getNear(idtRange,ndt,dtTar0,rdt,idt,jdt);
    size idt0 = static_cast<size>(dtTar1);
    for (k=1;k<tl;k++) {
        base = kV0[idt][iSyn][jSyn][idtRange[idt]+idt0+k];
        v[vs+k] += base + rdt * (kV0[jdt][iSyn][jSyn][idtRange[jdt]+idt0+k]-base);
        //if (k==1) {
        //    cout << "first value " << base << endl;
        //}
        //if (k==50) {
        //    cout << "5ms value " << base << endl;
        //}
    }
}

inline void interpVinit(std::vector<double> &v, size vs,  double **vLeak, double *vRange, size nv, double vTar, size tl) {
    size i,j,k;
    double r;
    getNear(vRange,nv,vTar,r,i,j);
    for (k=0;k<tl;k++)
        v[vs+k] = vLeak[i][k] + r*(vLeak[j][k]-vLeak[i][k]);
}

unsigned int bilinear_nSyn(std::vector<double> &v, nNL &neuroLib, nNS &neuron, double run_t, double ignore_t, std::vector<double> &tsp, double vCross, double vBack, vector<bool> &ei){
    double vTarget, dtTarget, dtTarget1;
    size tl, vs, ve, vc, i, ii, j, k, ith, ith_old, i_b = 0;
    double tstep = neuroLib.tstep;
    size itref = static_cast<size>(round(neuron.tRef/tstep));
    size nin = neuron.tin.size();
    size nv = neuroLib.nv;
    size ndt = neuroLib.ndt;
    size nt0 = neuroLib.nt;
    size l0 = nt0 - 1;
    size t0 = l0 * tstep;

    size nt = neuroLib.idtRange[ndt-1]+1;
    size l1 = nt - 1;
    size t1 = l1 * tstep;

    size tb = t1 - ignore_t;
    size lb = static_cast<size>(tb/tstep);
    size nb = lb + 1;

    size run_lt = static_cast<size>(run_t/tstep);
    size run_nt = run_lt + 1;

    bool crossed; 
    size ncross = 0;
    
    size limit, it, jt;
    unsigned int spikeCount = 0, spiked;
    
    double vinit = v[0];
    v.assign(run_nt,neuron.vRest);
    if (l0 > run_lt) {
        tl = run_nt;
    } else {
        tl = nt0;
    }
    interpVinit(v,0,neuroLib.vLeak,neuroLib.vRange,nv,vinit,tl);
    if (neuron.tin.size() > 0) {
        vs = static_cast<size>(neuron.tin[0]/tstep);
        for (i=0; i<neuron.tin.size(); i++) {
            if (sB::debug2) {
                cout << "i front " << i << endl;
                cout << "t " << vs*tstep << endl;
            }
            // linear
            if (neuron.tin[i] > run_t) {
                break;
            }
            vTarget = v[vs];

            if (vs+l0 > run_lt) {
                tl = run_nt-vs;
            } else {
                tl = nt0;
            }

            interpPSP0(v, vs, neuroLib.sPSP, neuroLib.vRange, neuroLib.nv, vTarget, tl, neuron.inID[i]);
            if (sB::debug2) {
                cout << " linear ok" << endl;
            }

            if (i<neuron.tin.size()-1) {
                ve = static_cast<size>(neuron.tin[i+1]/tstep);
                if (ve + 1 > vs + nt) {
                    vc = vs + nt;
                } else {
                    vc = ve + 1;
                }
            } else {
                vc = run_nt;
            }
            if (i_b == i && i != 0) cout << "first v after cross " << v[vs] << endl;
            for (k=i; k>i_b; k--) {
                j = k-1; // prevent negative for unsigned int iteration number. 
                dtTarget = neuron.tin[i]-neuron.tin[j];
                if (dtTarget > tb) break;
                dtTarget = dtTarget/tstep; // for interp along idtRange
                if (neuron.tin[j] + tb > run_t)
                    tl = run_nt - vs;
                else tl = static_cast<size>(round(nb - dtTarget));

                if (sB::debug2) {
                    cout << " b " << k << " v " << v[vs] << endl;
                    cout << " tl " << tl << " vs + tl " << vs + tl << endl;
                }
                interpkV(v, vs, neuroLib.kV, neuroLib.vRange, neuroLib.idtRange, neuroLib.nv, ndt, vTarget, dtTarget, dtTarget, tl, neuron.inID[j], neuron.inID[i]);
                if (sB::debug2) {
                    cout << " biv " << v[vs] << endl;
                }
            } 
            if (sB::debug2) {
                cout << " bilinear ok" << endl;
            }

            //check spike
            crossed = false;
            ith = i;
            for (k=vs;k<vc;k++) {
                while(v[k] > vCross) {
                    assert(k<vc);
                    crossed = true;
                    ncross = ncross +1;
                    if (sB::debug) {
                        std::cout << "crossed at " << k*tstep << ", start ith " << ith << " t= " << neuron.tin[ith] << std::endl;
                        std::cout << "v[k0] " << v[k]  << " vCross " << vCross << std::endl;
                    }
                    ith_old = ith; 
                    spiked = 1;
                    tsp.push_back(k*neuroLib.tstep+neuron.tRef/2);
                    vs = k + itref;
                    for (j=0;j<itref;j++) {
                        if (k+j < run_nt) {
                            v[k+j] = vCross;
                        } else {
                            vs = k+j-1;
                            break;
                        }
                    }
                    v[vs] = neuron.vRest;
                    for (ith = ith_old; ith < nin; ith++) {
                        if (neuron.tin[ith] > vs*tstep) {
                            break;                
                        }
                    }
                    ith--;
                    assert(ith>=ith_old);
                    i_b = ith + 1;
                    if (sB::debug) {
                        std::cout << "backed at " << vs*tstep << ", end ith " << ith << " t= " << neuron.tin[ith] << std::endl;
                        std::cout << "v[k1] " << v[vs]  << " vBack " << vBack << std::endl;
                    }
                    if (spiked){
                        spikeCount = spikeCount + 1;
                        neuron.tsp.push_back(tsp.back());
                    } 
                    if (vs + l0 <= run_lt) {
                        tl = nt0;
                    } else {
                        tl = run_nt-vs;
                    }
                    for (j=1; j<tl;j++) {
                        v[vs+j] = neuron.vRest;
                    }
                    vTarget = neuron.vRest;
                    //if (vs < run_lt) {
                    if (0) {
                        for (ii=ith+1; ii>ith+1; ii--) {
                            i = ii-1;
                            it = round(neuron.tin[i]/tstep);
                            dtTarget = vs - it;
                            if (l1 <= dtTarget) {
                                break;
                            }
                            if (neuron.tin[i] + t1 > run_t) {
                                //cout << neuron.tin[i] << " + " << t1 << " > " << run_t << endl;
                                limit = run_nt - it;
                            } else {
                                limit = nt;
                                //cout << neuron.tin[i] << " + " << t1 << " < " << run_t << endl;
                            }
                            interpPSP(v,vs, neuroLib.sPSP, neuroLib.vRange, neuroLib.idtRange, neuron.inID[i], neuroLib.nv, ndt, vTarget, dtTarget, limit);
                            for (k=i; k>0; k--) {
                                j = k-1; // prevent negative for unsigned int iteration number. 
                                jt = neuron.tin[j]/tstep;
                                dtTarget = it - jt;
                                dtTarget1 = vs - jt;
                                if (dtTarget > lb) {
                                    break;
                                }
                                if (neuron.tin[j] + tb > run_t)
                                    tl = run_nt - vs;
                                else tl = static_cast<size>(round(nb - dtTarget1));
                                //cout << "vs + tl " << vs << " + " << tl << " < " << run_nt << endl;
                                interpkV(v, vs, neuroLib.kV, neuroLib.vRange, neuroLib.idtRange, neuroLib.nv, ndt, vTarget, dtTarget, dtTarget1, tl, neuron.inID[j], neuron.inID[i]);
                            } 
                        }
                        vc = vs + l1; 
                        if (ith<neuron.tin.size()-1) {
                            ve = static_cast<size>(neuron.tin[ith+1]/tstep);
                            if (ve + 1 < vc) {
                                vc = ve + 1;
                            }
                        } else {
                            if (vc > run_lt) {
                                vc = run_nt;
                            }
                        }
                        k = vc-1;
                        for (j=vs; j<vc; j++) {
                            if (v[j] > vCross) {
                                k = j;
                                if (sB::debug) {
                                    std::cout << "cross again at " << k << std::endl;
                                }
                                break;
                            }
                        }
                    } else {
                        cout << " no readjust input, only leakage" << endl;
                        if (ith<neuron.tin.size()-1) {
                            ve = static_cast<size>(neuron.tin[ith+1]/tstep);
                        }
                        break;
                    }
                }
                if (crossed) {
                    i = ith;
                    cout << " last input during cross " << i << endl;
                    break;
                }
            }
            vs = ve;
            //cout << "i back" << i << "ve " << vs << endl;
        }
    }
    std::cout << "crossed " << ncross << " times" << std::endl;
    return spikeCount;
}

unsigned int linear_nSyn(std::vector<double> &v, nNL &neuroLib, nNS &neuron, double run_t, double ignore_t, std::vector<double> &tsp, double vCross, double vBack, vector<bool> &ei){
    double vTarget, dtTarget;
    size tl, vs, ve, vc, i, j, k, ith_old, ith;
    double tstep = neuroLib.tstep;
    size itref = static_cast<size>(round(neuron.tRef/tstep));
    size nin = neuron.tin.size();
    size nv = neuroLib.nv;
    size ndt = neuroLib.ndt;
    size nt0 = neuroLib.nt;
    size l0 = nt0 - 1;
    size t0 = l0 * tstep;

    size nt = neuroLib.idtRange[ndt-1]+1;
    size l1 = nt - 1;
    size t1 = l1 * tstep;

    size run_lt = static_cast<size>(run_t/tstep);
    size run_nt = run_lt + 1;

    bool crossed;
    size ncross = 0;
    //vector<bool> inTref(neuron.tin.size(),0);

    size limit, it;
    unsigned int spikeCount = 0, spiked;
    
    double vinit = v[0];
    v.assign(run_nt,neuron.vRest);
    if (l0 > run_lt) {
        tl = run_nt;
    } else {
        tl = nt0;
    }
    interpVinit(v,0,neuroLib.vLeak,neuroLib.vRange,nv,vinit,tl);
    if (neuron.tin.size() > 0) {
        vs = static_cast<size>(neuron.tin[0]/tstep);
        for (i=0; i<neuron.tin.size(); i++) {
            // linear
            if (neuron.tin[i] > run_t) {
                break;
            }
            //if (inTref[i]) continue;
            vTarget = v[vs];

            if (vs+l0 > run_lt) {
                tl = run_nt-vs;
            } else {
                tl = nt0;
            }

            interpPSP0(v, vs, neuroLib.sPSP, neuroLib.vRange, neuroLib.nv, vTarget, tl, neuron.inID[i]);
            
            if (i<neuron.tin.size()-1) {
                ve = static_cast<size>(neuron.tin[i+1]/tstep);
                if (ve + 1 > vs + nt0) {
                    vc = vs + nt0;
                } else {
                    vc = ve + 1;
                }
            } else {
                vc = run_nt;
            }
            
            //deal consecutive nearThreshold until 
            crossed = false;
            ith = i;
            for (k=vs;k<vc;k++) {
                while (v[k] > vCross) {
                    //std::cout << "k " << k << std::endl;
                    crossed = true;
                    ncross = ncross +1;
                    // come in and out in multiples of tstep 
                    if (sB::debug) {
                        std::cout << "crossed at " << k*tstep << ", start ith " << ith << " t= " << neuron.tin[ith] << std::endl;
                        std::cout << "v[k0] " << v[k]  << " vCross " << vCross << std::endl;
                    }
                    ith_old = ith; 
                    spiked = 1;
                    tsp.push_back(k*neuroLib.tstep+neuron.tRef/2);
                    vs = k + itref;
                    for (j=0;j<itref;j++) {
                        if (k+j < run_nt) {
                            v[k+j] = vCross;
                        } else {
                            vs = k+j-1;
                            break;
                        }
                    }
                    v[vs] = neuron.vRest;
                    for (ith = ith_old; ith < nin; ith++) {
                        if (neuron.tin[ith] > vs*tstep) {
                            break;                
                        }
                    }
                    ith--;
                    assert(ith>=ith_old);
                    if (sB::debug) {
                        std::cout << "backed at " << vs*tstep << ", end ith " << ith <<  " t= " << neuron.tin[ith] << std::endl;
                        std::cout << "v[k] " << v[vs]  << " vCross " << vCross << " vBack " << vBack << std::endl;
                    }
                    assert(k*tstep < vs*tstep);
                    if (spiked){
                        spikeCount = spikeCount + 1;
                        neuron.tsp.push_back(tsp.back());
                    } 
                    if (vs + l0 <= run_lt)
                        tl = nt0;
                    else tl = run_nt-vs;
                    for (j=1; j<tl;j++) {
                        v[vs+j] = neuron.vRest;
                    }
                    vTarget = neuron.vRest;
                    //if (vs < run_lt) {
                    if (0) {
                        for (i=ith+1;i>0; i--) {
                            j = i-1;
                            //if (inTref[j]) continue;
                            it = round(neuron.tin[j]/tstep);
                            dtTarget = vs - it;
                            if (l1 <= dtTarget) {
                                break;
                            }

                            if (neuron.tin[j] + t1 > run_t) {
                                limit = run_nt - it;
                            } else {
                                limit = nt;
                            }
                            interpPSP(v,vs, neuroLib.sPSP, neuroLib.vRange, neuroLib.idtRange, neuron.inID[j], neuroLib.nv, ndt, vTarget, dtTarget, limit);
                        }
                        vc = vs + l1; 
                        if (ith<neuron.tin.size()-1) {
                            ve = static_cast<size>(neuron.tin[ith+1]/tstep);
                            if (ve < vc) {
                                vc = ve;
                            }
                        } else {
                            if (vc > run_lt) {
                                vc = run_nt;
                            }
                        }
                        k = vc-1;
                        //cout << " vc " << vc << endl;
                        for (j=vs; j<vc; j++) {
                            if (v[j] > vCross) {
                                k = j;
                                if (sB::debug) {
                                    std::cout << "cross again at " << k << std::endl;
                                }
                                break;
                            }
                        }
                    } else {
                        cout << " no readjust input, only leakage" << endl;
                        if (ith<neuron.tin.size()-1) {
                            ve = static_cast<size>(neuron.tin[ith+1]/tstep);
                        }
                        break;
                    }
                }
                if (crossed) {
                    i = ith;
                    break;
                }
            }
            vs = ve;
        }
    }
    std::cout << "crossed " << ncross << " times" << std::endl;
    return spikeCount;
}

unsigned int bilinear0_nSyn(std::vector<double> &v, nNL &neuroLib, nNS &neuron, double run_t, double ignore_t, std::vector<double> &tsp, double vCross, double vBack, vector<bool> &ei){
    double vTarget, dtTarget, dtTarget1;
    size tl, vs, ve, vc, i, ii, j, k, ith, ith_old, i_b = 0;
    double tstep = neuroLib.tstep;
    size itref = static_cast<size>(round(neuron.tRef/tstep));
    size nin = neuron.tin.size();
    size nv = neuroLib.nv;
    size ndt = neuroLib.ndt;
    size nt0 = neuroLib.nt;
    size l0 = nt0 - 1;
    size t0 = l0 * tstep;

    size nt = neuroLib.idtRange[ndt-1]+1;
    size l1 = nt - 1;
    size t1 = l1 * tstep;

    size tb = t1 - ignore_t;
    size lb = static_cast<size>(tb/tstep);
    size nb = lb + 1;

    size run_lt = static_cast<size>(run_t/tstep);
    size run_nt = run_lt + 1;

    bool crossed; 
    size ncross = 0;
    
    
    size limit, it, jt;
    unsigned int spikeCount = 0, spiked;
    
    double vinit = v[0];
    v.assign(run_nt,neuron.vRest);
    if (l0 > run_lt) {
        tl = run_nt;
    } else {
        tl = nt0;
    }
    vTarget = neuron.vRest;
    if (vinit!=neuron.vRest){
        interpVinit(v,0,neuroLib.vLeak,neuroLib.vRange,nv,vinit,tl);
        cout << " bilinear0 need to start at resting potential for best accuracy" << endl;
    }
    double rv;
    size iv0, jv;
    getNear(neuroLib.vRange, neuroLib.nv, vTarget, rv, iv0, jv);
    cout << " vRest at " << iv0 << "th vRange" << endl;
    if (neuron.tin.size() > 0) {
        vs = static_cast<size>(neuron.tin[0]/tstep);
        for (i=0; i<neuron.tin.size(); i++) {
            //cout << "i front " << i << endl;
            //cout << "t " << vs*tstep << endl;
            // linear
            if (neuron.tin[i] > run_t) {
                break;
            }

            if (vs+l0 > run_lt) {
                tl = run_nt-vs;
            } else {
                tl = nt0;
            }
            interpPSP_nV(v, vs, neuroLib.sPSP, neuron.inID[i], 0, tl, iv0);

            if (i<neuron.tin.size()-1) {
                ve = static_cast<size>(neuron.tin[i+1]/tstep);
                if (ve + 1 > vs + nt) {
                    vc = vs + nt;
                } else {
                    vc = ve + 1;
                }
            } else {
                vc = run_nt;
            }
            if (i_b == i && i != 0) cout << "first v after cross " << v[vs] << endl;
            for (k=i; k>i_b; k--) {
                j = k-1; // prevent negative for unsigned int iteration number
                dtTarget = neuron.tin[i]-neuron.tin[j];
                if (dtTarget > tb) break;
                dtTarget = dtTarget/tstep; // for interp along idtRange
                if (neuron.tin[j] + tb > run_t)
                    tl = run_nt - vs;
                else tl = static_cast<size>(round(nb - dtTarget));
                //cout << j << "th " << neuron.ei[neuron.inID[i]] << ", " << i <<"th " << neuron.ei[neuron.inID[j]] <<  endl;
                interpkV0(v, vs, neuroLib.kV0, neuroLib.idtRange, ndt, dtTarget, 0, tl, neuron.inID[j], neuron.inID[i]);
            } 
            //check spike
            crossed = false;
            ith = i;
            for (k=vs;k<vc;k++) {
                while(v[k] > vCross) {
                    assert(k<vc);
                    crossed = true;
                    ncross = ncross +1;
                    if (sB::debug) {
                        std::cout << "crossed at " << k*tstep << ", start ith " << ith << " t= " << neuron.tin[ith] << std::endl;
                        std::cout << "v[k0] " << v[k]  << " vCross " << vCross << std::endl;
                    }
                    ith_old = ith; 
                    spiked = 1;
                    tsp.push_back(k*neuroLib.tstep+neuron.tRef/2);
                    vs = k + itref;
                    for (j=0;j<itref;j++) {
                        if (k+j < run_nt) {
                            v[k+j] = vCross;
                        } else {
                            vs = k+j-1;
                            break;
                        }
                    }
                    v[vs] = neuron.vRest;
                    for (ith = ith_old; ith < nin; ith++) {
                        if (neuron.tin[ith] > vs*tstep) {
                            break;                
                        }
                    }
                    ith--;
                    assert(ith>=ith_old);
                    i_b = ith + 1;
                    if (sB::debug) {
                        std::cout << "backed at " << vs*tstep << ", end ith " << ith << " t= " << neuron.tin[ith] << std::endl;
                        std::cout << "v[k1] " << v[vs]  << " vBack " << vBack << std::endl;
                    }
                    if (spiked){
                        spikeCount = spikeCount + 1;
                        neuron.tsp.push_back(tsp.back());
                    } 
                    if (vs + l0 <= run_lt)
                        tl = nt0;
                    else tl = run_nt-vs;
                    for (j=1; j<tl;j++) {
                        v[vs+j] = neuron.vRest;
                    }
                    vTarget = neuron.vRest;
                    //if (vs < run_lt) {
                    cout << " no readjust input, only leakage" << endl;
                    if (ith<neuron.tin.size()-1) {
                        ve = static_cast<size>(neuron.tin[ith+1]/tstep);
                    }
                    break;
                }
                if (crossed) {
                    i = ith;
                    cout << " last input during cross " << i << endl;
                    break;
                }
            }
            vs = ve;
            //cout << "i back" << i << "ve " << vs << endl;
        }
    }
    std::cout << "crossed " << ncross << " times" << std::endl;
    return spikeCount;
}
#endif
