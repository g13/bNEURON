#include "jumpy_linear.h"

bool jl::check_crossing(Input &input, nNL &neuroLib, Cross &cross, nNS &neuron, double tol_tl, double tol_tb, double end_t, size tail_l, size head, jND &jnd, double &t_cross, double vC) {
    // check all tmax after head and find the v > vC
    size i, j;
    double t, dt;
    size nt = static_cast<size>(end_t/neuroLib.tstep);
    double v, v_pass; 
    bool lcrossed = false;
    double tCross = cross.tCross.back();
    for (i=tail_l; i<=head; i++) {
        t = input.tMax[i]+input.t[i];
        if (t <= jnd.t.back() || input.ID[i] >= neuroLib.nE) {
            // ignore inh input and input that maxed out before input.t[head] 
            continue;
        }
        // ignore tmax that go beyond next input
        if (head < neuron.tin.size()-1) {
            if (t >= neuron.tin[head+1]/neuroLib.tstep) {
                continue;
            }
        } else {
            if ( t > nt ) {
                t = nt;
            }
        }
        v = jl::find_v_at_t(input, neuroLib, cross, head, tail_l, t, tCross, tol_tl, tol_tb, neuron.vReset);
        if (v > vC) {
            lcrossed = true;
            break;    
        }
    }
    if (lcrossed) {
        // linear confirm crossed
        // perform linear (and parabola) interp iteration until tolerance reaches
        if (jl::debug) {
            cout << "v " << v << " > vThres, find t and v for cross" << endl;
        }
        if (fabs(v-vC)>neuron.vTol && t-jnd.t.back()>1) {
            t_cross = jl::interp_for_t_cross(v, jnd.v.back(), t, jnd.t.back(), head, tail_l, tCross, tol_tl, tol_tb, neuroLib, cross, input, neuron.vRest, neuron.vTol, vC, v_pass);
            if (jl::debug) {
                cout << "t_left " << jnd.t.back() << " t_cross " << t_cross << ", t_right " << t << endl;
                cout << "v_left " << jnd.v.back() << " v_cross " << v_pass << ", t_right " << v << endl;
                assert(t_cross >= jnd.t.back());
            }
            jnd.t.push_back(t_cross);
            jnd.v.push_back(v_pass);
        } else {
            t_cross = t;
            jnd.t.push_back(t);
            jnd.v.push_back(v);    
        }
        if (jl::debug) {
            cout << "set v: " << jnd.v.back() << endl;
        }
        return true;
    }
    return false;
}

bool jl::update_vinit_of_new_input_check_crossing(Input &input, Cross &cross, nNL &neuroLib, nNS &neuron, size head, size &tail_l, size old_tail_l, double tol_tl, double tol_tb, double end_t, jND &jnd, double &t_cross, double vC, size corrSize) {
    size i, j, idt;
    double dt;
    double tstep = neuroLib.tstep;
    double v_pass;
    double tCross = cross.tCross.back();
    double v;
    if (jl::debug) {
        cout << " adding " << head << " input t " << input.t[head] << " synapse ID " << neuron.inID[head] << endl;
    }
    input.ID.push_back(neuron.inID[head]);
    dt = input.t[head] - tCross;
    if (jl::debug) {
        if (dt < -1e-14 ) {
            cout << dt << " < 0 " << endl;
            cout << head << endl;
            cout << neuron.tin[head] << endl;
            cout << input.t[head] << endl;
            cout << tCross*tstep << endl;
            assert(dt>=0);
        }
    }
    // initialize to 
    if (dt < tol_tl) {
        v = add_vinit_contribution(neuroLib.vLeak, cross.vCross.back(), dt);
    } else {
        v = neuron.vRest;
    }
    for (i=tail_l; i<head; i++) {
        dt = input.t[head] - input.dt[i];
        idt = static_cast<size>(round(dt));
        add_input_i_contribution(i,idt,neuroLib,input,v);
    }
    // update head's v and tmax
    jnd.t.push_back(input.t[head]);
    jnd.v.push_back(v);
    // check if vmax at left bound very unlikely
    if (v > vC) {
        if (jl::debug2) {
            cout << " crossing upon input not very probable at low input rate " << endl;
        }
        jnd.t.pop_back();
        jnd.v.pop_back();
        // perform linear interp iteration until tolerance reaches
        tail_l = old_tail_l;
        t_cross = jl::interp_for_t_cross(v, jnd.v.back(), input.t[head], jnd.t.back(), head-1, tail_l, tCross, tol_tl, tol_tb, neuroLib, cross, input, neuron.vRest, neuron.vTol, vC, v_pass);
        if (jl::debug2) {
            if (t_cross-input.t[head] > 1e-14) {
                assert(t_cross<=input.t[head]);
            }
        }
        jnd.t.push_back(t_cross);
        jnd.v.push_back(v_pass);
        input.t.pop_back();
        input.ID.pop_back();
        return true;
    } else {
        add_new_input_info(neuroLib.vRange, neuroLib.nv, input, cross, neuron.tin[head], neuroLib.tMax, neuroLib.sfireCap, v, corrSize, neuron.inID[head]);
        if (jl::debug2) {
            if (input.t.size() != head+1) {
                cout << "before check " << input.t.size() << " != " << head+1 << endl;
                assert(input.t.size() == head+1);
            }
        }
    }
    // check all tmax after head and find the upper limit of vmax
    return check_crossing(input, neuroLib, cross, neuron, tol_tl, tol_tb, end_t, tail_l, head, jnd, t_cross, vC);
}

void jl::update_info_after_cross(Input &input, nNL &neuroLib, Cross &cross, nNS &neuron, double tCross, double vCross, size i_prior, size &tail, size head, size corrSize, int afterCrossBehavior) {
    size i_, j_;
    double r_;
    double dt, v;
    size i, j;
    size idt;
    // inputs that only matters during cross, can fill junk with those
    size i_cross = i_prior + 1;
    size i_start = i_cross;
    if (jl::debug) {
        if (input.t.size() != i_start) {
            cout << "before update " << input.t.size() << " != " << i_start << endl;
            assert(input.t.size() == i_start);
        }
    }
    if (!afterCrossBehavior) {
        for (i=i_cross; i<=head; i++) {
            input.junk_this();
        }
        // new start after cross
        tail = head+1;
    } else {
        if (i_cross < tail) {
            if (jl::debug) {
                cout << " crossing part lasting longer than PSP sample length, very unlikely" << endl;
                cout << " dead old input " << tail-i_cross << endl;
            }
            for (i=i_cross; i<tail; i++) {
                input.junk_this();
            }
            i_start = tail;
        } else {
            // update for input that come before cross and lingers after cross;
            if (jl::debug) {
                cout << "cross update starting tail " << tail << endl;
                cout << " # lingering old inputs: " << i_cross - tail << endl;
            }
            for (i=tail; i<i_cross; i++) {
                if (jl::debug) {
                    cout << "lingering " << i << " < " << neuron.tin.size() << endl;
                }
                dt = tCross - input.t[i];
                if (jl::debug) {
                    assert(dt>0);
                }
                jb::getNear(neuroLib.idtRange,neuroLib.ndt,
                              dt, input.dTijr[i].r, input.dTijr[i].i, input.dTijr[i].j);
                input.dt[i] = tCross;
                input.cCross[i] = cross.nCross;
                input.Vijr[i] = cross.vCross.back();
                if (input.Vijr[i].j >= neuroLib.sfireCap[input.ID[i]][input.dTijr[i].i] && input.dTijr[i].r > 0.1) {
                    // if sufficiently close to spike case
                    input.tMax[i] = neuroLib.tMax[input.Vijr[i].j][input.dTijr[i].i][input.ID[i]];
                } else {
                    input.tMax[i] = linear_interp_tMax(neuroLib.tMax, cross.vCross.back(), input.dTijr[i],input.ID[i]);
                }
                if (jl::debug) {
                    cout << "   loop ended" << endl;
                }
            }
        }
        // update for new input during the cross that lingers after cross
        if (jl::debug) {
            cout << " # lingering new inputs: " << head-i_start+1 << endl;
            cout << input.t.size() << " == " << input.Vijr.size() << endl;
            if (input.t.size() != i_start) {
                cout << input.t.size() << " != " << i_start << endl;
                assert(input.t.size() == i_start);
            }
        }
        for (i=i_start; i<=head; i++) {
            input.t.push_back(neuron.tin[i]/neuroLib.tstep);
            input.ID.push_back(neuron.inID[i]);
            input.dt.push_back(tCross);
            input.idt.push_back(0); // never used in jl, but useful in jb
            dt = tCross - input.t[i];
            jb::getNear(neuroLib.idtRange, neuroLib.ndt, dt, r_, i_, j_);
            if (jl::debug) {
                assert(input.t.size() == i+1);
                if (dt < -1e-14) {
                    cout << tCross << " - " << input.t[i] << " = " << dt << endl;
                    assert(dt >= 0);
                }
                if (r_-1>1e-14) {
                    cout << tCross  << "-" << input.t[i] << endl;
                    cout << "tail " << tail << ", " << input.t[tail] << endl;
                    cout << i << " head " << head << "total " << neuron.tin.size() << endl;
                    cout << "dt " << dt << " i_ " << i_ << " j_ " << j_ << endl;
                    cout << "r_ " << r_ << endl;
                    assert(r_ <= 1);
                }
            }
            input.dTijr.push_back(IJR(i_,j_,r_));
            input.cCross.push_back(cross.nCross);
            input.Vijr.push_back(cross.vCross.back());
            if (input.Vijr[i].j >= neuroLib.sfireCap[input.ID[i]][input.dTijr[i].i] && input.dTijr[i].r > 0.1) {
                // if sufficiently close to spike case
                input.tMax.push_back(neuroLib.tMax[input.Vijr[i].j][input.dTijr[i].i][input.ID[i]]);
            } else {
                input.tMax.push_back(linear_interp_tMax(neuroLib.tMax, cross.vCross.back(), input.dTijr[i],input.ID[i]));
            }
            input.Tmpijr.push_back(IJR(0,1,0)); 
            input.bir.push_back(BilinearRelationships(corrSize));
        }
    }
    if (jl::debug) {
        input.assert_size();
    }
}

unsigned int nsyn_jLinear(nNS &neuron, nNL &neuroLib, Input &input, jND &jnd, Cross &cross, double end_t, double ignore_t, size corrSize, vector<double> &tsp, double vC, double vB, int afterCrossBehavior, bool spikeShape) {
    size i, j, i_prior_cross;
    double tstep = neuroLib.tstep;
    unsigned int nc_old, nc = 0;
    double t_cross, vBack, tBack, dt;
    double iend = end_t/tstep;
    double tol_tl = neuroLib.nt;
    double tol_tb = round((neuroLib.tol_tb-ignore_t)/tstep);
    cout << "linear corr length " << tol_tl << endl;
    cout << "after cross linear corr length " << tol_tb << endl;
    cout << "total inputs " << neuron.tin.size() << endl;
    size tail_l = 0;
    size old_tail_l;
    bool crossed, spiked;
    double tref = neuron.tRef/tstep;
    size i_, j_;
    double r_;
    size ii;
    jb::getNear(neuroLib.vRange, neuroLib.nv, cross.v[0], r_, i_, j_);
    cross.vCross.push_back(IJR(i_,j_,r_));
    jnd.t.push_back(0);
    jnd.v.push_back(cross.v[0]);
    if (jl::debug) {
        cout << "set v: " << jnd.v.back() << endl;
    }
    for (i=0;i<neuron.tin.size();i++) {
        input.t.push_back(neuron.tin[i]/tstep);
        if (jl::debug) {
            cout << " new input i = " << i << endl;
            assert(input.t.size() == i+1);
        }
        old_tail_l = tail_l;
        dt = input.t[i]-cross.tCross.back();
        if (dt < tol_tb) { 
            move_corr_window_i(input.t, tail_l, input.t[i], tol_tb);
        } else {
            if (dt < tol_tl) {
                move_corr_window_i(input.t, tail_l, input.t[i], dt);
            } else{
                move_corr_window_i(input.t, tail_l, input.t[i], tol_tl);
            }
        }
        crossed = jl::update_vinit_of_new_input_check_crossing(input, cross, neuroLib, neuron, i, tail_l, old_tail_l, tol_tl, tol_tb, end_t, jnd, t_cross, vC, corrSize);
        ii = 0;
        while (crossed) {
            ii++;
            if (jl::debug) {
                cout << " crossing " << endl;
                if (ii > 1) cout << " AGAIN?!!" << endl;
            }
            if (t_cross < neuron.tin[i]/tstep + 1e-14) {
                if (jl::debug) {
                    cout << " unlikely, cross upon or before input" << endl;
                }
                i_prior_cross = i-1;
            } else {
                i_prior_cross = i;
            }
            if (jl::debug) {
                cout << " vCross = " << jnd.v.back();
                cout << " tCross = " << jnd.t.back() << endl;
                cout << " input " << i_prior_cross << ", t" << input.t[i_prior_cross]*tstep << endl;
                if (input.t.size() != i_prior_cross+1) {
                    cout << "before alter " << input.t.size() << " != " << i_prior_cross+1 << endl;
                    assert(input.t.size() == i_prior_cross+1);
                }
            }
            nc_old = nc;
            if (spikeShape) {
                size ith = i_prior_cross;
                spiked = true;
                vBack = cross.v.back();
                i = ith;
                nc = nc + spiked;
            } else {
                tBack = t_cross+tref;
                size tl = tref;
                if (tBack > iend) {
                    tBack = iend;
                    tl = iend - ceil(t_cross);
                }
                for (j=0;j<tl;j++) {
                    cross.v.push_back(vC);
                    cross.t.push_back(ceil(t_cross)+j);
                }
                vBack = neuron.vRest;
                cross.v.push_back(vBack);
                cross.t.push_back(tBack);
                double tmpTsp = t_cross*tstep + neuron.tRef/2;
                if (tmpTsp <= end_t) {
                    tsp.push_back(tmpTsp);
                    nc++;
                }
                spiked = true;
                i = i_prior_cross + 1;
                if (i < neuron.inID.size()) {
                    while (neuron.tin[i]-1e-14 < tBack*tstep+1e-14) {
                        i++;
                        if (i==neuron.inID.size()) {
                            break;
                        }
                    }
                }
                i--;
            }
            jnd.t.push_back(tBack);
            jnd.v.push_back(vBack);
            if (jl::debug) {
                cout << "set v: " << jnd.v.back() << endl;
                cout << " backed at " << tBack*neuroLib.tstep << " with v: " << vBack << endl;
                cout << " input from " << i_prior_cross << " to " << i << endl;
                cout << " input at " << neuron.tin[i] << endl;
                if (i+1 < neuron.tin.size()) {
                    cout << " next input " << neuron.tin[i+1] << endl;
                }
            }
            if (vBack >= vC) {
                break;
            }
            jb::getNear(neuroLib.vRange,neuroLib.nv,
                            vBack, r_, i_, j_);
            cross.vCross.push_back(IJR(i_,j_,r_));
            cross.nCross++;
            cross.iCross.push_back(cross.v.size());
            cross.tCross.push_back(tBack);
            if (jl::debug) {
                cout << "crossed " << endl;
            }

            old_tail_l = tail_l;
            move_corr_window(neuron.tin, tail_l, tBack, tol_tb, tstep);
            jl::update_info_after_cross(input, neuroLib, cross, neuron, tBack, vBack, i_prior_cross, tail_l, i, corrSize, afterCrossBehavior);
            if (jl::debug) {
                cout << " updated" << endl;
            }
            crossed = jl::check_crossing(input, neuroLib, cross, neuron, tol_tl, tol_tb, end_t, tail_l, i, jnd, t_cross, vC);
            if (jl::debug) {
                cout << "------------" << endl;
            }
        }
        if (jl::debug) {
            cout << " dead with " << i << ", v " << jnd.v.back() << " at " << jnd.t.back() << endl;
            input.print_this(i);
            input.assert_size();
        }
    }
    cout << " crossed " << cross.nCross << " times" << endl;
    return nc;
}
