#include "linear_bilinear.h"

unsigned int bilinear_nSyn(std::vector<double> &v, nNL &neuroLib, nNS &neuron, double run_t, double ignore_t, std::vector<double> &tsp, double vCross, double vBack, vector<bool> &ei, int afterSpikeBehavior){
    bool debug;
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
            if (lb::debug2) {
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
            if (lb::debug2) {
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
            if (i_b == i && i != 0) cout << "first input after cross " << v[vs] << endl;
            for (k=i; k>i_b; k--) {
                j = k-1; // prevent negative for unsigned int iteration number. 
                dtTarget = neuron.tin[i]-neuron.tin[j];
                if (dtTarget > tb) break;
                dtTarget = dtTarget/tstep; // for interp along idtRange
                if (neuron.tin[j] + tb > run_t) {
                    tl = run_nt - vs;
                } else {
                    tl = static_cast<size>(round(nb - dtTarget));
                }

                if (lb::debug2) {
                    cout << j << "th " << neuron.ei[neuron.inID[i]] << ", " << i <<"th " << neuron.ei[neuron.inID[j]] <<  endl;
                    cout << " b " << k << " v " << v[vs] << endl;
                    cout << " tl " << tl << " vs + tl " << vs + tl << endl;
                }
                interpkV(v, vs, neuroLib.kV, neuroLib.vRange, neuroLib.idtRange, neuroLib.nv, ndt, vTarget, dtTarget, dtTarget, tl, neuron.inID[j], neuron.inID[i]);
                if (lb::debug2) {
                    cout << " biv " << v[vs] << endl;
                }
            } 
            if (lb::debug2) {
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
                    if (lb::debug) {
                        std::cout << "crossed at " << k*tstep << ", start ith " << ith << " t= " << neuron.tin[ith] << std::endl;
                        std::cout << "v[k0] " << v[k]  << " vCross " << vCross << std::endl;
                    }
                    ith_old = ith; 
                    spiked = 1;
                    double tmpTsp = k*neuroLib.tstep+neuron.tRef/2;
                    if (tmpTsp <= run_t) {
                        tsp.push_back(tmpTsp);
                    }
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
                    if (lb::debug) {
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
                    vTarget = v[vs];
                    //if (vs < run_lt) {
                    if (afterSpikeBehavior) {
                        cout << " adding afterSpike sPSP " << endl;
                        for (ii=ith+1; ii>0; ii--) {
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
                            if (afterSpikeBehavior==2) {
                                i_b = 0;
                                cout << " adding afterspike bPSP " << endl;
                                for (k=i; k>0; k--) {
                                    j = k-1; // prevent negative for unsigned int iteration number. 
                                    jt = neuron.tin[j]/tstep;
                                    dtTarget = it - jt;
                                    dtTarget1 = vs - jt;
                                    if (dtTarget > lb) {
                                        break;
                                    } 
                                    if (neuron.tin[j] + tb > run_t) {
                                        tl = run_nt - vs;
                                    } else { 
                                        tl = static_cast<size>(round(nb - dtTarget1));
                                    }
                                    //cout << "vs + tl " << vs << " + " << tl << " < " << run_nt << endl;
                                    interpkV(v, vs, neuroLib.kV, neuroLib.vRange, neuroLib.idtRange, neuroLib.nv, ndt, vTarget, dtTarget, dtTarget1, tl, neuron.inID[j], neuron.inID[i]);
                                } 
                            }
                        }
                        if (ith<neuron.tin.size()-1) {
                            ve = static_cast<size>(neuron.tin[ith+1]/tstep);
                            if (ve + 1 > vs + nt) {
                                vc = vs + nt;
                            } else {
                                vc = ve + 1;
                            }
                        } else {
                            vc = run_nt;
                        }
                        k = vc-1; // initial value in case vs == vc
                        for (j=vs; j<vc; j++) {
                            if (v[j] > vCross) {
                                k = j;
                                if (lb::debug) {
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

unsigned int linear_nSyn(std::vector<double> &v, nNL &neuroLib, nNS &neuron, double run_t, double ignore_t, std::vector<double> &tsp, double vCross, double vBack, vector<bool> &ei, int afterSpikeBehavior){
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
                    assert(k<vc);
                    crossed = true;
                    ncross = ncross +1;
                    if (lb::debug) {
                        std::cout << "crossed at " << k*tstep << ", start ith " << ith << " t= " << neuron.tin[ith] << std::endl;
                        std::cout << "v[k0] " << v[k]  << " vCross " << vCross << std::endl;
                    }
                    ith_old = ith; 
                    spiked = 1;
                    double tmpTsp = k*neuroLib.tstep+neuron.tRef/2;
                    if (tmpTsp <= run_t) {
                        tsp.push_back(tmpTsp);
                    }
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
                    if (lb::debug) {
                        std::cout << "backed at " << vs*tstep << ", end ith " << ith <<  " t= " << neuron.tin[ith] << std::endl;
                        std::cout << "v[k] " << v[vs]  << " vCross " << vCross << " vBack " << vBack << std::endl;
                    }
                    assert(k*tstep < vs*tstep);
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
                    vTarget = v[vs];
                    //if (vs < run_lt) {
                    if (afterSpikeBehavior) {
                        for (i=ith+1;i>0; i--) {
                            j = i-1;
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
                        if (ith<neuron.tin.size()-1) {
                            ve = static_cast<size>(neuron.tin[ith+1]/tstep);
                            if (ve + 1 > vs + nt) {
                                vc = vs + nt;
                            } else {
                                vc = ve + 1;
                            }
                        } else {
                            vc = run_nt;
                        }
                        k = vc-1;
                        //cout << " vc " << vc << endl;
                        for (j=vs; j<vc; j++) {
                            if (v[j] > vCross) {
                                k = j;
                                if (lb::debug) {
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
        }
    }
    std::cout << "crossed " << ncross << " times" << std::endl;
    return spikeCount;
}

unsigned int bilinear0_nSyn(std::vector<double> &v, nNL &neuroLib, nNS &neuron, double run_t, double ignore_t, std::vector<double> &tsp, double vCross, double vBack, vector<bool> &ei, int afterSpikeBehavior, bool kVStyle){
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
            if (lb::debug2) {
                cout << "i front " << i << endl;
                cout << "t " << vs*tstep << endl;
            }
            // linear
            if (neuron.tin[i] > run_t) {
                break;
            }

            if (vs+l0 > run_lt) {
                tl = run_nt-vs;
            } else {
                tl = nt0;
            }
            interpPSP_nV(v, vs, neuroLib.sPSP, neuroLib.idtRange, ndt, neuron.inID[i], 0, tl, iv0);
            if (lb::debug2) {
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
                j = k-1; // prevent negative for unsigned int iteration number
                dtTarget = neuron.tin[i]-neuron.tin[j];
                if (dtTarget > tb) break;
                dtTarget = dtTarget/tstep; // for interp along idtRange
                if (neuron.tin[j] + tb > run_t)
                    tl = run_nt - vs;
                } else  {
                    tl = static_cast<size>(round(nb - dtTarget));
                }

                if (lb::debug2) {
                    cout << j << "th " << neuron.ei[neuron.inID[i]] << ", " << i <<"th " << neuron.ei[neuron.inID[j]] <<  endl;
                    cout << " b " << k << " v " << v[vs] << endl;
                    cout << " tl " << tl << " vs + tl " << vs + tl << endl;
                }
                interpkV0(v, vs, neuroLib.kV0, neuroLib.idtRange, ndt, dtTarget, 0, tl, neuron.inID[j], neuron.inID[i]);
                if (lb::debug2) {
                    cout << " biv0 " << v[vs] << endl;
                }
            } 
            if (lb::debug2) {
                cout << " bilinear0 ok" << endl;
            }

            //check spike
            crossed = false;
            ith = i;
            for (k=vs;k<vc;k++) {
                while(v[k] > vCross) {
                    assert(k<vc);
                    crossed = true;
                    ncross = ncross +1;
                    if (lb::debug) {
                        std::cout << "crossed at " << k*tstep << ", start ith " << ith << " t= " << neuron.tin[ith] << std::endl;
                        std::cout << "v[k0] " << v[k]  << " vCross " << vCross << std::endl;
                    }
                    ith_old = ith; 
                    spiked = 1;
                    double tmpTsp = k*neuroLib.tstep+neuron.tRef/2;
                    if (tmpTsp <= run_t) {
                        tsp.push_back(tmpTsp);
                    }
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
                    if (lb::debug) {
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
                    vTarget = v[vs];
                    //if (vs < run_lt) {
                    if (afterSpikeBehavior) {
                        cout << " adding afterSpike sPSP " << endl;
                        for (ii=ith+1; ii>0; ii--) {
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
                            interpPSP_nV(v, vs, neuroLib.sPSP, neuroLib.idtRange, ndt, neuron.inID[i], dtTarget, limit, iv0, kVStyle);
                            if (afterSpikeBehavior==2) {
                                i_b = 0;
                                cout << " adding afterspike bPSP " << endl;
                                for (k=i; k>0; k--) {
                                    j = k-1; // prevent negative for unsigned int iteration number. 
                                    jt = neuron.tin[j]/tstep;
                                    dtTarget = it - jt;
                                    dtTarget1 = vs - jt;
                                    if (dtTarget > lb) {
                                        break;
                                    } 
                                    if (neuron.tin[j] + tb > run_t) {
                                        tl = run_nt - vs;
                                    } else { 
                                        tl = static_cast<size>(round(nb - dtTarget1));
                                    }
                                    //cout << "vs + tl " << vs << " + " << tl << " < " << run_nt << endl;
                                    if (kVStyle) {
                                        interpkV(v, vs, neuroLib.kV, neuroLib.vRange, neuroLib.idtRange, neuroLib.nv, ndt, neuron.vRest, dtTarget, dtTarget1, tl, neuron.inID[j], neuron.inID[i]);
                                    } else {
                                        interpkV0(v, vs, neuroLib.kV0, neuroLib.idtRange, ndt, dtTarget, dtTarget1, tl, neuron.inID[j], neuron.inID[i]);
                                    }
                                }
                            }
                        }
                        if (ith<neuron.tin.size()-1) {
                            ve = static_cast<size>(neuron.tin[ith+1]/tstep);
                            if (ve + 1 > vs + nt) {
                                vc = vs + nt;
                            } else {
                                vc = ve + 1;
                            }
                        } else {
                            vc = run_nt;
                        }
                        k = vc-1; // initial value in case vs == vc
                        for (j=vs; j<vc; j++) {
                            if (v[j] > vCross) {
                                k = j;
                                if (lb::debug) {
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
