#include "linear_bilinear.h"

unsigned int bilinear_nSyn(Cell &cell, vector<vector<double>> &spikeTrain, vector<double> dendVclamp, double rd, vector<double> &v, nNL &neuroLib, nNS &neuron, double run_t, double ignore_t, vector<double> &tsp, double vCross, double vBack, int afterCrossBehavior, bool spikeShape, bool dtSquare, int itrial, bool sliceDebugPlot){
    double vTarget, dtTarget;
    size tl, vs, ve, vc, i, ii, j, k, ith, ith_old, i_b = 0;
    size nin = neuron.tin.size(), nv = neuroLib.nv, ndt = neuroLib.ndt;
    double tstep = neuroLib.tstep;
    size itref = static_cast<size>(round(neuron.tRef/tstep));
    // upper limits of time lengths:
    // length of PSP0
    size nt0 = neuroLib.nt;
    size l0 = nt0 - 1;
    double t0 = l0 * tstep;

    // length of PSP and kV from interp of dt
    size nt1 = neuroLib.idtRange[ndt-2] + (nt0 - neuroLib.idtRange[ndt-1]);
    size l1 = nt1 - 1;
    double t1 = l1 * tstep;

    // input time difference
    size nt2 = neuroLib.idtRange[ndt-1]+1;
    size l2 = nt2 - 1;
    double t2 = l2 * tstep;

    // ignore some time difference
    double tb = t1 - ignore_t;
    size lb = static_cast<size>(tb/tstep);
    size nb = lb + 1;
    cout << "nt0 " << nt0 << endl;
    cout << "nt1 " << nt1 << endl;
    cout << "nt2 " << nt2 << endl;
    cout << "nb " << nb << endl;

    size run_lt = static_cast<size>(run_t/tstep);
    size run_nt = run_lt + 1;

    bool crossed;
    size ncross = 0;
    unsigned int spikeCount = 0;

    vector<size> s0(neuroLib.nSyn,0);
    vector<size> s1(neuroLib.nSyn,0);

    double tCross_old = 0;
    vector<double> vS(neuron.tin.size(),neuron.vRest);

    double vinit = v[0];
    double vCross_old = vinit;
    v.assign(run_nt, neuron.vRest);
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
            vS[i] = vTarget;

            if (vs+l0 > run_lt) {
                tl = run_nt - vs;
            } else {
                tl = nt0;
            }

            interpPSP0(v, vs, neuroLib.sPSP, neuroLib.vRange, neuroLib.nv, vTarget, tl, neuron.inID[i]);
            if (lb::debug2) {
                cout << " linear ok" << endl;
            }

            if (i<neuron.tin.size()-1) {
                ve = static_cast<size>(neuron.tin[i+1]/tstep);
                if (ve > vs + l2) {
                    vc = vs + nt2;
                } else {
                    vc = ve + 1;
                }
            } else {
                vc = run_nt;
            }
            if (lb::debug2) {
                if (i_b == i && i != 0) cout << "first input after cross " << v[vs] << endl;
            }
            for (k=i; k>i_b; k--) {
                j = k-1; // prevent negative for unsigned int iteration number
                dtTarget = neuron.tin[i]-neuron.tin[j];
                if (dtTarget > tb) {
                    break;
                }
                dtTarget = dtTarget/tstep; // for interp along idtRange
                if (neuron.tin[j] + t1 > run_t) {
                    tl = run_nt - vs;
                } else {
                    tl = static_cast<size>(round(nt1 - dtTarget));
                }

                if (lb::debug2) {
                    cout << j << "th " << neuron.ei[neuron.inID[j]] << ", " << i <<"th " << neuron.ei[neuron.inID[i]] <<  endl;
                    cout << " b " << k << " v " << v[vs] << endl;
                    cout << " tl " << tl << " vs + tl " << vs + tl << endl;
                    cout << " dtTar  " << dtTarget << endl;
                    cout << " run_nt  " << run_nt << endl;
                }
                interpkV(v, vs, neuroLib.kV, neuroLib.vRange, neuroLib.idtRange, neuroLib.nv, ndt, vTarget, dtTarget, dtTarget, tl, neuron.inID[j], neuron.inID[i], dtSquare);
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
                unsigned int spiked;
                while (v[k] > vCross) {
                    assert(k<vc);
                    crossed = true;
                    ncross = ncross +1;
                    if (lb::debug) {
                        cout << "crossed at " << k*tstep << ", start ith " << ith << " t= " << neuron.tin[ith] << endl;
                        cout << "v[vs] " << v[k]  << " vCross " << vCross << endl;
                    }
                    if (spikeShape) {
                        vs = k;
                        clampDendRaw(neuroLib, neuron, t1, vs*tstep, v[k], dendVclamp, ith, rd, tCross_old, vCross_old, vS);
                        string fign;
                        if (sliceDebugPlot) {
                            fign = "bi-" + to_string(itrial) + "-" + to_string(ncross);
                        } else {
                            fign = "";
                        }
                        spiked = neuroAlterB(neuron, neuroLib, v, ith, vs, run_nt, tstep, tsp, v[vs], vBack, cell, spikeTrain, s0, s1, dendVclamp,fign);
                    } else {
                        spiked = 1;
                        tsp.push_back(k*tstep);
                        vs = k + itref;
                        for (j=0;j<=itref;j++) {
                            if (k+j >= run_nt) {
                                vs = run_lt;
                                break;
                            }
                            v[k+j] = neuron.vRest;
                        }
                    }
                    ith_old = ith; 
                    for (ith = ith_old; ith < nin; ith++) {
                        if (neuron.tin[ith] > vs*tstep) {
                            break;                
                        }
                    }
                    ith--;
                    assert(ith>=ith_old);
                    if (lb::debug) {
                        cout << "backed at " << vs*tstep << ", end ith " << ith << " t= " << neuron.tin[ith] << endl;
                        cout << "v[ve] " << v[vs] << " vCross " << vCross << " vBack " << vBack << endl;
                    }
                    tCross_old = vs*tstep;
                    vCross_old = v[vs];
                    if (v[vs] > vBack) {
                        if (vs < run_lt) {
                            cout << " NEURON simulation has some problem" << endl;
                            assert(vs == run_lt);
                        }
                    }
                    if (spiked){
                        spikeCount = spikeCount + spiked;
                        neuron.tsp.push_back(tsp.back());
                    } 
                    if (run_nt-vs == 1) {
                        cout << " time runs out while crossing " << "\n";
                        break;
                    }
                    if (spikeShape) {
                        if (spiked) {
                            if (vs + neuroLib.nvASt-1 <= run_lt) {
                                tl = neuroLib.nvASt;
                            } else {
                                tl = run_nt-vs;
                            }
                            interpVAS(v,vs,neuroLib.vAS,neuroLib.vASrange,neuroLib.nvAS,v[vs],neuron.vRest,tl);
                            cout << "vas ended with " << v[vs+tl-1] << " at " << vs+tl-1 << endl;
                        } else {
                            if (vs + l0 <= run_lt) {
                                tl = nt0;
                            } else {
                                tl = run_nt-vs;
                            }
                            interpVinit(v,vs,neuroLib.vLeak,neuroLib.vRange,nv,v[vs],tl);
                        }
                    } else {
                        if (vs + l0 <= run_lt) {
                            tl = nt0;
                        } else {
                            tl = run_nt-vs;
                        }
                        for (j=1; j<tl;j++) {
                            v[vs+j] = neuron.vRest;
                        }
                    }

                    if (afterCrossBehavior>0) {
                        vTarget = v[vs];
                        for (ii=ith+1; ii>0; ii--) {
                            size i = ii-1;
                            size it = round(neuron.tin[i]/tstep);
                            dtTarget = vs - it;
                            if (lb < dtTarget) {
                                break;
                            }

                            if (run_t < neuron.tin[i] + t1) {
                                tl = run_nt - vs;
                            } else {
                                tl = static_cast<size>(round(nt1 - dtTarget));
                            }
                            interpPSP(v,vs, neuroLib.sPSP, neuroLib.vRange, neuroLib.idtRange, neuron.inID[i], neuroLib.nv, ndt, vTarget, dtTarget, tl);
                            if (afterCrossBehavior==2) {
                                i_b = 0;
                                for (k=i; k>0; k--) {
                                    j = k-1; // prevent negative for unsigned int iteration number. 
                                    size jt = neuron.tin[j]/tstep;
                                    dtTarget = it - jt;
                                    double dtTarget1 = vs - jt;
                                    if (dtTarget > lb || dtTarget1 > lb) {
                                        break;
                                    } 
                                    if (neuron.tin[j] + t1 > run_t) {
                                        tl = run_nt - vs;
                                    } else {
                                        tl = static_cast<size>(round(nt1 - dtTarget1));
                                    }
                                    //cout << "vs + tl " << vs << " + " << tl << " < " << run_nt << endl;
                                    interpkV(v, vs, neuroLib.kV, neuroLib.vRange, neuroLib.idtRange, neuroLib.nv, ndt, vTarget, dtTarget, dtTarget1, tl, neuron.inID[j], neuron.inID[i], dtSquare);
                                } 
                            }
                        }
                        if (lb::debug) {
                            cout << "   finished afterspike sbPSP adjustment" << endl;
                        }
                        if (ith<neuron.tin.size()-1) {
                            ve = static_cast<size>(neuron.tin[ith+1]/tstep);
                            if (ve > vs + l2) {
                                vc = vs + nt2;
                            } else {
                                vc = ve + 1;
                            }
                        } else {
                            vc = run_nt;
                        }
                        // initial value in case vs == vc
                        k = vc-1;
                        if (lb::debug) {
                            cout << " vc " << vc << endl;
                            cout << "v k at vc-1: " << v[k] << endl;
                        }
                        for (j=vs; j<vc; j++) {
                            if (v[j] > vCross) {
                                k = j;
                                if (lb::debug) {
                                    cout << "need to cross again at " << k << endl;
                                }
                                break;
                            }
                        }
                    } else {
                        if (lb::debug) {
                            cout << " no readjust input, only leakage" << endl;
                        }
                        if (ith<neuron.tin.size()-1) {
                            ve = static_cast<size>(neuron.tin[ith+1]/tstep);
                            i_b = ith + 1;
                        }
                        break;
                    }
                }
                if (crossed) {
                    i = ith;
                    if (lb::debug) {
                        cout << " last input during cross " << i << endl;
                    }
                    break;
                }
            }
            vs = ve;
            //cout << "i back" << i << "ve " << vs << endl;
        }
    }
    cout << "crossed " << ncross << " times" << endl;
    return spikeCount;
}

unsigned int linear_nSyn(Cell &cell, vector<vector<double>> &spikeTrain, vector<double> dendVclamp, double rd, vector<double> &v, nNL &neuroLib, nNS &neuron, double run_t, double ignore_t, vector<double> &tsp, double vCross, double vBack, int afterCrossBehavior, bool spikeShape, int itrial, bool sliceDebugPlot){
    double vTarget, dtTarget;
    size tl, vs, ve, vc, i, j, k, ith, ith_old;
    size nin = neuron.tin.size(), nv = neuroLib.nv, ndt = neuroLib.ndt;
    double tstep = neuroLib.tstep;
    size itref = static_cast<size>(round(neuron.tRef/tstep));
    // upper limits of time lengths:
    // length of PSP0
    size nt0 = neuroLib.nt;
    size l0 = nt0 - 1;
    double t0 = l0 * tstep;

    // length of PSP from interp of dt
    size nt1 = neuroLib.idtRange[ndt-2] + (nt0 - neuroLib.idtRange[ndt-1]);
    size l1 = nt1 - 1;
    double t1 = l1 * tstep;

    // input time difference
    size nt2 = neuroLib.idtRange[ndt-1]+1;
    size l2 = nt2 - 1;
    double t2 = l2 * tstep;

    // ignore some time difference
    double tb = t1 - ignore_t;
    size lb = static_cast<size>(tb/tstep);
    size nb = lb + 1;
    cout << "nt0 " << nt0 << endl;
    cout << "nt1 " << nt1 << endl;
    cout << "nt2 " << nt2 << endl;
    cout << "nb " << nb << endl;

    size run_lt = static_cast<size>(run_t/tstep);
    size run_nt = run_lt + 1;

    bool crossed;
    size ncross = 0;
    unsigned int spikeCount = 0;

    vector<size> s0(neuroLib.nSyn,0);
    vector<size> s1(neuroLib.nSyn,0);

    double tCross_old = 0;
    vector<double> vS(neuron.tin.size(),neuron.vRest);

    double vinit = v[0];
    double vCross_old = vinit;
    v.assign(run_nt, neuron.vRest);
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
            vS[i] = vTarget;

            if (vs+l0 > run_lt) {
                tl = run_nt - vs;
            } else {
                tl = nt0;
            }

            interpPSP0(v, vs, neuroLib.sPSP, neuroLib.vRange, neuroLib.nv, vTarget, tl, neuron.inID[i]);
            if (lb::debug2) {
                cout << " linear ok" << endl;
            }

            if (i<neuron.tin.size()-1) {
                ve = static_cast<size>(neuron.tin[i+1]/tstep);
                if (ve > vs + l2) {
                    vc = vs + nt2;
                } else {
                    vc = ve + 1;
                }
            } else {
                vc = run_nt;
            }
            
            crossed = false;
            ith = i;
            for (k=vs;k<vc;k++) {
                unsigned int spiked;
                while (v[k] > vCross) {
                    assert(k<vc);
                    crossed = true;
                    ncross = ncross +1;
                    if (lb::debug) {
                        cout << "crossed at " << k*tstep << ", start ith " << ith << " t= " << neuron.tin[ith] << endl;
                        cout << "v[vs] " << v[k]  << " vCross " << vCross << endl;
                    }
                    if (spikeShape) {
                        vs = k;
                        clampDendRaw(neuroLib, neuron, t1, vs*tstep, v[k], dendVclamp, ith, rd, tCross_old, vCross_old, vS);
                        string fign;
                        if (sliceDebugPlot) {
                            fign = "li-" + to_string(itrial) + "-" + to_string(ncross);
                        } else {
                            fign = "";
                        }
                        spiked = neuroAlterB(neuron, neuroLib, v, ith, vs, run_nt, tstep, tsp, v[vs], vBack, cell, spikeTrain, s0, s1, dendVclamp,fign);
                    } else {
                        spiked = 1;
                        tsp.push_back(k*tstep);
                        vs = k + itref;
                        for (j=0;j<=itref;j++) {
                            if (k+j >= run_nt) {
                                vs = run_lt;
                                break;
                            }
                            v[k+j] = neuron.vRest;
                        }
                    }
                    ith_old = ith; 
                    for (ith = ith_old; ith < nin; ith++) {
                        if (neuron.tin[ith] > vs*tstep) {
                            break;                
                        }
                    }
                    ith--;
                    assert(ith>=ith_old);
                    if (lb::debug) {
                        cout << "backed at " << vs*tstep << ", end ith " << ith << " t= " << neuron.tin[ith] << endl;
                        cout << "v[ve] " << v[vs]  << " vCross " << vCross << " vBack " << vBack << endl;
                    }
                    tCross_old = vs*tstep;
                    vCross_old = v[vs];
                    if (v[vs] > vBack) {
                        if (vs < run_lt) {
                            cout << " NEURON simulation has some problem" << endl;
                            assert(vs == run_lt);
                        }
                    }
                    if (spiked){
                        spikeCount = spikeCount + spiked;
                        neuron.tsp.push_back(tsp.back());
                    } 
                    if (run_nt-vs == 1) {
                        cout << " time runs out while crossing " << "\n";
                        break;
                    }
                    if (spikeShape) {
                        if (spiked) {
                            if (vs + neuroLib.nvASt-1 <= run_lt) {
                                tl = neuroLib.nvASt;
                            } else {
                                tl = run_nt-vs;
                            }
                            interpVAS(v,vs,neuroLib.vAS,neuroLib.vASrange,neuroLib.nvAS,v[vs],neuron.vRest,tl);
                            cout << "vas ended with " << v[vs+tl-1] << " at " << vs+tl-1 << endl;
                        } else {
                            if (vs + l0 <= run_lt) {
                                tl = nt0;
                            } else {
                                tl = run_nt-vs;
                            }
                            interpVinit(v,vs,neuroLib.vLeak,neuroLib.vRange,nv,v[vs],tl);
                        }
                    } else {
                        if (vs + l0 <= run_lt) {
                            tl = nt0;
                        } else {
                            tl = run_nt-vs;
                        }
                        for (j=1; j<tl;j++) {
                            v[vs+j] = neuron.vRest;
                        }
                    }

                    if (afterCrossBehavior>0) {
                        vTarget = v[vs];
                        for (i=ith+1; i>0; i--) {
                            j = i-1;
                            size it = round(neuron.tin[j]/tstep);
                            dtTarget = vs - it;
                            if (lb < dtTarget) {
                                break;
                            }

                            if (run_t < neuron.tin[j] + t1) {
                                tl = run_nt - vs;
                            } else {
                                tl = static_cast<size>(round(nt1 - dtTarget));
                            }
                            interpPSP(v,vs, neuroLib.sPSP, neuroLib.vRange, neuroLib.idtRange, neuron.inID[j], neuroLib.nv, ndt, vTarget, dtTarget, tl);
                        }
                        if (ith<neuron.tin.size()-1) {
                            ve = static_cast<size>(neuron.tin[ith+1]/tstep);
                            if (ve > vs + l2) {
                                vc = vs + nt2;
                            } else {
                                vc = ve + 1;
                            }
                        } else {
                            vc = run_nt;
                        }
                        // initial value in case vs == vc
                        k = vc-1;
                        if (lb::debug) {
                            cout << " vc " << vc << endl;
                            cout << "v k at vc-1: " << v[k] << endl;
                        }
                        for (j=vs; j<vc; j++) {
                            if (v[j] > vCross) {
                                k = j;
                                if (lb::debug) {
                                    cout << "need to cross again at " << k << endl;
                                }
                                break;
                            }
                        }
                    } else {
                        if (lb::debug) {
                            cout << " no readjust input, only leakage" << endl;
                        }
                        if (ith<neuron.tin.size()-1) {
                            ve = static_cast<size>(neuron.tin[ith+1]/tstep);
                        }
                        break;
                    }
                }
                if (crossed) {
                    i = ith;
                    if (lb::debug) {
                        cout << " last input during cross " << i << endl;
                    }
                    break;
                }
            }
            vs = ve;
        }
    }
    cout << "crossed " << ncross << " times" << endl;
    return spikeCount;
}

unsigned int bilinear0_nSyn(Cell &cell, vector<vector<double>> &spikeTrain, vector<double> dendVclamp, double rd, vector<double> &v, nNL &neuroLib, nNS &neuron, double run_t, double ignore_t, vector<double> &tsp, double vCross, double vBack, int afterCrossBehavior, bool spikeShape, bool kVStyle, bool dtSquare, int itrial, bool sliceDebugPlot){
    double vTarget, dtTarget;
    size tl, vs, ve, vc, i, ii, j, k, ith, ith_old, i_b = 0;
    size nin = neuron.tin.size(), nv = neuroLib.nv, ndt = neuroLib.ndt;
    double tstep = neuroLib.tstep;
    size itref = static_cast<size>(round(neuron.tRef/tstep));
    // upper limits of time lengths:
    // length of PSP0
    size nt0 = neuroLib.nt;
    size l0 = nt0 - 1;
    double t0 = l0 * tstep;

    // length of PSP and kV from interp of dt
    size nt1 = neuroLib.idtRange[ndt-2] + (nt0 - neuroLib.idtRange[ndt-1]);
    size l1 = nt1 - 1;
    double t1 = l1 * tstep;

    // input time difference
    size nt2 = neuroLib.idtRange[ndt-1]+1;
    size l2 = nt2 - 1;
    double t2 = l2 * tstep;

    // ignore some time difference
    double tb = t1 - ignore_t;
    size lb = static_cast<size>(tb/tstep);
    size nb = lb + 1;
    cout << "nt0 " << nt0 << endl;
    cout << "nt1 " << nt1 << endl;
    cout << "nt2 " << nt2 << endl;
    cout << "nb " << nb << endl;

    size run_lt = static_cast<size>(run_t/tstep);
    size run_nt = run_lt + 1;

    bool crossed;
    size ncross = 0;

    unsigned int spikeCount = 0;

    vector<size> s0(neuroLib.nSyn,0);
    vector<size> s1(neuroLib.nSyn,0);

    double tCross_old = 0;
    vector<double> vS(neuron.tin.size(),neuron.vRest);

    double vinit = v[0];
    double vCross_old = vinit;
    v.assign(run_nt, neuron.vRest);
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
            vTarget = v[vs];
            vS[i] = vTarget;

            if (vs+l0 > run_lt) {
                tl = run_nt - vs;
            } else {
                tl = nt0;
            }

            interpPSP_nV(v, vs, neuroLib.sPSP, neuroLib.idtRange, ndt, neuron.inID[i], 0, tl, iv0);
            if (lb::debug2) {
                cout << " linear ok" << endl;
            }

            if (i<neuron.tin.size()-1) {
                ve = static_cast<size>(neuron.tin[i+1]/tstep);
                if (ve > vs + l2) {
                    vc = vs + nt2;
                } else {
                    vc = ve + 1;
                }
            } else {
                vc = run_nt;
            }
            if (lb::debug2) {
                if (i_b == i && i != 0) cout << "first input after cross " << v[vs] << endl;
            }
            for (k=i; k>i_b; k--) {
                j = k-1; // prevent negative for unsigned int iteration number
                dtTarget = neuron.tin[i]-neuron.tin[j];
                if (dtTarget > tb) {
                    break;
                }
                dtTarget = dtTarget/tstep; // for interp along idtRange
                if (neuron.tin[j] + t1 > run_t) {
                    tl = run_nt - vs;
                } else {
                    tl = static_cast<size>(round(nt1 - dtTarget));
                }

                if (lb::debug2) {
                    cout << j << "th " << neuron.ei[neuron.inID[j]] << ", " << i <<"th " << neuron.ei[neuron.inID[i]] <<  endl;
                    cout << " b " << k << " v " << v[vs] << endl;
                    cout << " tl " << tl << " vs + tl " << vs + tl << endl;
                    cout << " dtTar  " << dtTarget << endl;
                    cout << " run_nt  " << run_nt << endl;
                }
                interpkV0(v, vs, neuroLib.kV0, neuroLib.idtRange, ndt, dtTarget, tl, neuron.inID[j], neuron.inID[i]);
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
                unsigned int spiked;
                while (v[k] > vCross) {
                    assert(k<vc);
                    crossed = true;
                    ncross = ncross +1;
                    if (lb::debug) {
                        cout << "crossed at " << k*tstep << ", start ith " << ith << " t= " << neuron.tin[ith] << endl;
                        cout << "v[vs] " << v[k]  << " vCross " << vCross << endl;
                    }
                    if (spikeShape) {
                        vs = k;
                        clampDendRaw(neuroLib, neuron, t1, vs*tstep, v[k], dendVclamp, ith, rd, tCross_old, vCross_old, vS);
                        string fign;
                        if (sliceDebugPlot) {
                            fign = "bi0-" + to_string(itrial) + "-" + to_string(ncross);
                        } else {
                            fign = "";
                        }
                        spiked = neuroAlterB(neuron, neuroLib, v, ith, vs, run_nt, tstep, tsp, v[vs], vBack, cell, spikeTrain, s0, s1, dendVclamp,fign);
                    } else {
                        spiked = 1;
                        tsp.push_back(k*tstep);
                        vs = k + itref;
                        for (j=0;j<=itref;j++) {
                            if (k+j >= run_nt) {
                                vs = run_lt;
                                break;
                            }
                            v[k+j] = neuron.vRest;
                        }
                    }
                    ith_old = ith; 
                    for (ith = ith_old; ith < nin; ith++) {
                        if (neuron.tin[ith] > vs*tstep) {
                            break;                
                        }
                    }
                    ith--;
                    assert(ith>=ith_old);
                    if (lb::debug) {
                        cout << "backed at " << vs*tstep << ", end ith " << ith << " t= " << neuron.tin[ith] << endl;
                        cout << "v[ve] " << v[vs] << " vCross " << vCross << " vBack " << vBack << endl;
                    }
                    tCross_old = vs*tstep;
                    vCross_old = v[vs];
                    if (v[vs] > vBack) {
                        if (vs < run_lt) {
                            cout << " NEURON simulation has some problem" << endl;
                            assert(vs == run_lt);
                        }
                    }
                    if (spiked){
                        spikeCount = spikeCount + spiked;
                        neuron.tsp.push_back(tsp.back());
                    } 
                    if (run_nt-vs == 1) {
                        cout << " time runs out while crossing " << "\n";
                        break;
                    }
                    if (spikeShape) {
                        if (spiked) {
                            if (vs + neuroLib.nvASt-1 <= run_lt) {
                                tl = neuroLib.nvASt;
                            } else {
                                tl = run_nt-vs;
                            }
                            interpVAS(v,vs,neuroLib.vAS,neuroLib.vASrange,neuroLib.nvAS,v[vs],neuron.vRest,tl);
                            cout << "vas ended with " << v[vs+tl-1] << " at " << vs+tl-1 << endl;
                        } else {
                            if (vs + l0 <= run_lt) {
                                tl = nt0;
                            } else {
                                tl = run_nt-vs;
                            }
                            interpVinit(v,vs,neuroLib.vLeak,neuroLib.vRange,nv,v[vs],tl);
                        }
                    } else {
                        if (vs + l0 <= run_lt) {
                            tl = nt0;
                        } else {
                            tl = run_nt-vs;
                        }
                        for (j=1; j<tl;j++) {
                            v[vs+j] = neuron.vRest;
                        }
                    }

                    if (afterCrossBehavior>0) {
                        vTarget = v[vs];
                        for (ii=ith+1; ii>0; ii--) {
                            size i = ii-1;
                            size it = round(neuron.tin[i]/tstep);
                            dtTarget = vs - it;
                            if (lb < dtTarget) {
                                break;
                            }

                            if (run_t < neuron.tin[i] + t1) {
                                tl = run_nt - vs;
                            } else {
                                tl = static_cast<size>(round(nt1 - dtTarget));
                            }
                            interpPSP_nV(v, vs, neuroLib.sPSP, neuroLib.idtRange, ndt, neuron.inID[i], dtTarget, tl, iv0);
                            if (afterCrossBehavior==2) {
                                i_b = 0;
                                for (k=i; k>0; k--) {
                                    j = k-1; // prevent negative for unsigned int iteration number. 
                                    size jt = neuron.tin[j]/tstep;
                                    dtTarget = it - jt;
                                    double dtTarget1 = vs - jt;
                                    if (dtTarget > lb || dtTarget1 > lb) {
                                        break;
                                    } 
                                    if (neuron.tin[j] + t1 > run_t) {
                                        tl = run_nt - vs;
                                    } else {
                                        tl = static_cast<size>(round(nt1 - dtTarget1));
                                    }
                                    //cout << "vs + tl " << vs << " + " << tl << " < " << run_nt << endl;
                                    interpkV(v, vs, neuroLib.kV, neuroLib.vRange, neuroLib.idtRange, neuroLib.nv, ndt, neuron.vRest, dtTarget, dtTarget1, tl, neuron.inID[j], neuron.inID[i], dtSquare);
                                }
                            }
                        }
                        if (lb::debug) {
                            cout << "   finished afterspike sbPSP adjustment" << endl;
                        }
                        if (ith<neuron.tin.size()-1) {
                            ve = static_cast<size>(neuron.tin[ith+1]/tstep);
                            if (ve > vs + l2) {
                                vc = vs + nt2;
                            } else {
                                vc = ve + 1;
                            }
                        } else {
                            vc = run_nt;
                        }
                        // initial value in case vs == vc
                        k = vc-1;
                        if (lb::debug) {
                            cout << " vc " << vc << endl;
                            cout << "v k at vc-1: " << v[k] << endl;
                        }
                        for (j=vs; j<vc; j++) {
                            if (v[j] > vCross) {
                                k = j;
                                if (lb::debug) {
                                    cout << "need to cross again at " << k << endl;
                                }
                                break;
                            }
                        }
                    } else {
                        if (lb::debug) {
                            cout << " no readjust input, only leakage" << endl;
                        }
                        if (ith<neuron.tin.size()-1) {
                            ve = static_cast<size>(neuron.tin[ith+1]/tstep);
                            i_b = ith + 1;
                        }
                        break;
                    }
                }
                if (crossed) {
                    i = ith;
                    if (lb::debug) {
                        cout << " last input during cross " << i << endl;
                    }
                    break;
                }
            }
            vs = ve;
            //cout << "i back" << i << "ve " << vs << endl;
        }
    }
    cout << "crossed " << ncross << " times" << endl;
    return spikeCount;
}
