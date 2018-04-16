#include "jumpy_bilinear.h"

double interp_for_t_cross(double v_right, double v_left, double t_right, double t_left, size head, size tail_l, size tail_b, double tCross, double tol_tl, nNL &neuroLib, Cross &cross, Input &input, double v0, double vtol, double vC, double &v, bool debug) {
    size ival = 0;
    double v1, v2;
    double t_cross1, t_cross2, t_cross;
    do {
        if (debug) {
            assert(vC >= v_left);
        }
        if (abs(t_right - t_left)<1e-14) {
            v = v_right;
            return t_right;
        }
        t_cross1 = t_left + (vC - v_left)/(v_right-v_left)*(t_right-t_left);
        //t_cross1 = ceil(t_cross1);
        if (jb::debug2) {
            cout << " find v at " << t_cross1 << endl;
            cout <<"t: " << t_left << ", " << t_right << endl;
            assert(t_cross1 >= t_left);
            assert(t_right > t_left);
        }
        v1 = find_v_at_t(input, neuroLib, cross, head, tail_l, tail_b, t_cross1, tCross, tol_tl, v0, debug);
        if (jb::debug2) {
            cout << v_left <<  ", " << v1 << ", " << v_right << endl;
        }
        ival++;
        if (fabs(v1-vC)<vtol) {
            v = v1;
            t_cross = t_cross1;
            break;
        }
        if (t_cross1-t_left <= 1 || t_right-t_cross1 <=1) {
            if (t_cross1-t_left <= 1) {
                v = v1; 
            } else {
                v = v_right;
            }
            t_cross = t_cross1;
            if (jb::debug2) {
                cout << "temp resol reached " << endl;
            }
            break;
        }
        // solve for a, b of f(t)-v_left = a(t-t_left)^2 + b(t-t_left);
        // solve for t when a(t-t_left)^2 + b(t-t_left) = (vThres - v_left)
        t_cross2 = parabola(t_left,v_left,t_right,v_right,t_cross1, v1, vC);
        //t_cross2 = ceil(t_cross2);
        if (jb::debug2) {
            cout << " v: " << v_left << ", " << v1 << ", " << v_right << endl;
            cout << " t: " << t_left << ", " << t_cross1 << ", " << t_right << endl;
            cout << "t_cross2 " << t_cross2 << endl;
            assert(t_cross2 >= t_left);
            assert(t_cross2 <= t_right);
        }
        // if t interval smaller than sample temp resolution, apply the quadratic interpolation as solution
        if (jb::debug2) {
            cout << " find v at " << t_cross2 << endl;
        }
        v2 = find_v_at_t(input, neuroLib, cross, head, tail_l, tail_b, t_cross2, tCross, tol_tl, v0, debug);
        //if (jb::debug2) {
            ival++;
        //}
        if (fabs(v2-vC)<vtol) {
            v = v2;
            t_cross = t_cross2;
            break;
        }
        if (jb::debug2) {
            cout << " v: " << v_left << ", " << v1 << ", " << v2 << ", " << v_right << endl;
            cout << " t: " << t_left << ", " << t_cross1 << ", " << t_cross2 << ", " << t_right << endl;
        }
        if (v2>v1) {
            getLR(v_left,v1,v2,v_right,t_left,t_cross1,t_cross2,t_right, vC);
        } else {
            getLR(v_left,v2,v1,v_right,t_left,t_cross2,t_cross1,t_right, vC);
        } 
        if (jb::debug2) {
            cout << " v: " << v_left << ", " << v_right << endl;
            cout << " t: " << t_left << ", " << t_right << endl;
        }
    } while (true);
    if (jb::debug2) {
        cout << ival << " evaluations" << endl;
    }
    return t_cross;
}

bool check_crossing(Input &input, nNL &neuroLib, Cross &cross, nNS &neuron, double tol_tl, double tol_tb, double end_t, size tail_l, size tail_b, size head, jND &jnd, double &t_cross, double vC, bool debug) {
    // check all tmax after head and find the upper limit of vmax
    size i, j, it;
    double t, dt;
    size nt = end_t/neuroLib.tstep;
    double v, v_pass; 
    double vmax = jnd.v.back();
    double tmax = jnd.t.back();
    size vi_tail;
    size vi_tail_max;
    double tCross = cross.tCross.back();
    for (i=tail_l; i<=head; i++) {
        t = input.tMax[i]+input.t[i];
        vi_tail = 0;
        if (t <= jnd.t.back() || neuroLib.ei[input.ID[i]]) {
            // ignore inh input and input that uncorred after cross
            continue;
        }
        // ignore tmax that go beyond next input
        if (head < neuron.tin.size()-1) {
            if (t >= neuron.tin[head+1]/neuroLib.tstep) {
                continue;
            }
        } else {
            if ( t > nt ){
                t = nt;
            }
        }
        // initialize with leak from last cross
        dt = t-tCross;
        if (dt < tol_tl) {
            v = add_vinit_contribution(neuroLib.vLeak, cross.vCross.back(), dt);
        } else {
            v = neuron.vRest;
        }
        for (j=tail_l; j<=head; j++) {
            dt = t - input.t[j];
            if (dt > tol_tl) {
                continue;
            } else { 
                // get tail for tmax
                if ( dt < tol_tb ) {
                    if (!vi_tail) {
                        vi_tail = j;
                    }
                }
            }
            dt = t - input.dt[j];
            it = static_cast<size>(round(dt));
            add_input_i_contribution(j, it, neuroLib, input, v);
        }
        if (v > vmax) {
            vmax = v;
            tmax = t; 
            vi_tail_max = vi_tail;
        }
    }
    if (vmax > vC) {
        if (jb::debug2) {
            cout << "linear vmax " << vmax <<  " > vThres" << endl;
        }
        // perform linear interp iteration until tolerance reaches
        if (jb::debug) {
            cout << "l+b vmax > vThres, find t and v for cross" << endl;
        }
        if (fabs(vmax-vC)>neuron.vTol && tmax-jnd.t.back()>1) {
            t_cross = interp_for_t_cross(vmax, jnd.v.back(), tmax, jnd.t.back(), head, tail_l, tail_b, tCross, tol_tl, neuroLib, cross, input, neuron.vRest, neuron.vTol, vC, v_pass, debug);
            if (jb::debug) {
                assert(t_cross > jnd.t.back());
            }
            jnd.t.push_back(t_cross);
            jnd.v.push_back(v_pass);
        } else {
            t_cross = tmax;
            jnd.t.push_back(tmax);    
            jnd.v.push_back(vmax);    
        }
        if (jb::debug) {
            cout << "set v: " << jnd.v.back() << endl;
        }
        return true;
    }
    return false;
}

bool update_vinit_of_new_input_check_crossing(Input &input, Cross &cross, nNL &neuroLib, nNS &neuron, size head, size tail_l, size tail_b, double tol_tl, double tol_tb, double end_t, jND &jnd, double &t_cross, double vC, size corrSize, bool debug) {
    size i, j, ii,jj;
    double dt;
    double tstep = neuroLib.tstep;
    double v_pass;
    size idt;
    double tCross = cross.tCross.back();
    double v;
    //if (jb::debug) {
    //    cout << " adding " << head << " input t " << input.t[head] << " synapse ID " << neuron.inID[head] << endl;
    //}
    input.ID.push_back(neuron.inID[head]);
    dt = input.t[head] - tCross;
    if (jb::debug) {
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
    if (debug) {
        cout << " last v " << jnd.v.back() << endl;
        if (head-1 < tail_l && tail_l != 0) {
            cout << "first v after cross ";
        } else {
            cout << " start with leak";
        }
        cout << v << endl;
    }
    //cout << " tailing input " << tail_l << endl;
    double vtmp = v;
    double dvtmp;
    if (debug) {
        cout << "head: " << head << " at t: " << input.t[head] << endl;
        cout << "tail_l: " << tail_l << " at t: " << input.t[tail_l] << endl;
        cout << "tail_b: " << tail_l << " at t: " << input.t[tail_b] << endl;
        cout << "last cross at t: " << cross.tCross.back() << endl;
    }

    for (ii=head; ii>tail_l; ii--) {
        i = ii-1;
        if (debug) {
            dvtmp = 0;
            cout << "   linear input: " << i << " at t: " << input.t[i] << " <= dt: " << input.dt[i] << endl;
        }
        dt = input.t[head] - input.dt[i];
        input.idt[i] = static_cast<size>(round(dt));
        input.bir[i].ID.push_back(neuron.inID[head]);
        add_input_i_contribution(i,input.idt[i],neuroLib,input,v);
        if (debug) {
            cout << "   contributing " << v-vtmp << " with synapse ID: " << input.ID[i] << endl;
            vtmp = v;
        }
        if (i<tail_b) continue;
        for (jj=head; jj>ii; jj--) {
            j = jj-1;
            if (debug){
                cout << "       +bilinear input: " << j << " at t: " << input.t[j] << " <= dt: " << input.dt[j] << endl;
            }
            add_input_i_j_bilinear_contribution(input, neuroLib, i, j, input.idt[j], v, debug);
            if (debug) {
                dvtmp += v-vtmp;
                cout << "       contributing " << v-vtmp << " with synapse ID " << input.ID[j] << " with dT i: " << input.bir[i].dTijr[j-i-1].i << " j: " << input.bir[i].dTijr[j-i-1].j << " r: " << input.bir[i].dTijr[j-i-1].r << " at V i: " << input.Vijr[j].i << " j: " << input.Vijr[j].j << " r: " << input.Vijr[j].r << endl;
                vtmp = v;
            }
        }
        if (debug) {
            cout << "   bilinear total change: " << dvtmp << endl;
            cout << "   current v: " << v << endl;
            cout << endl;
        }
    }
    // update head's v and tmax
    jnd.t.push_back(input.t[head]);
    jnd.v.push_back(v);
    if (jb::debug) {
        cout << "set v: " << jnd.v.back() << endl;
    }
    // check if vmax at left bound very unlikely
    if (v > vC) {
        if(jb::debug2) {
            cout << " crossing upon input not very probable at low input rate " << endl;
        }
        jnd.t.pop_back();
        jnd.v.pop_back();
        // perform linear interp iteration until tolerance reaches
        if (jb::debug2) {
            cout << "t: " << jnd.t.back()  << " < " << input.t[head] << endl;
        }
        t_cross = interp_for_t_cross(v, jnd.v.back(), input.t[head], jnd.t.back(), head-1, tail_l, tail_b, tCross, tol_tl, neuroLib, cross, input, neuron.vRest, neuron.vTol, vC, v_pass, debug);
        if (jb::debug2) {
            if (t_cross-input.t[head] > 1e-14) {
                assert(t_cross<=input.t[head]);
            }
        }
        jnd.t.push_back(t_cross);
        jnd.v.push_back(v_pass);
        if (jb::debug) {
            cout << "set v: " << jnd.v.back() << endl;
        }
        input.t.pop_back();
        input.ID.pop_back();
        for (ii=head; ii>tail_l; ii--) {
            i = ii-1;
            input.bir[i].ID.pop_back();
        }
        if (jb::debug2) {
            if (input.t.size() != head) {
                cout << "before check " << input.t.size() << " != " << head << endl;
                assert(input.t.size() == head);
            }
        }
        return true;
    } else {
        add_new_input_info(neuroLib.vRange, neuroLib.nv, input, cross, neuron.tin[head], neuroLib.tMax, v, corrSize, neuron.inID[head]);
        for (ii=head; ii>tail_l; ii--) {
            i = ii-1;
            add_relation_between_new_and_ith_input(input, head, i, neuroLib.idtRange, neuroLib.ndt);
        }
        if (jb::debug2) {
            if (input.t.size() != head+1) {
                cout << "before check " << input.t.size() << " != " << head+1 << endl;
                assert(input.t.size() == head+1);
            }
        }
    }
    // check all tmax after head and find the upper limit of vmax
    return check_crossing(input, neuroLib, cross, neuron, tol_tl, tol_tb, end_t, tail_l, tail_b, head, jnd, t_cross, vC, debug);
}

void update_info_after_cross(Input &input, nNL &neuroLib, Cross &cross, nNS &neuron, double tCross, double vCross, size i_prior, size &tail_l, size &tail_b, size head, size corrSize, int afterCrossBehavior) {
    size i_, j_;
    double r_;
    double dt, v;
    size i, j;
    size idt;
    // inputs that only matters during cross, can fill junk with those
    size i_cross = i_prior + 1;
    size i_start = i_cross;
    if (jb::debug) {
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
        tail_l = head+1;
        tail_b = head+1;
    } else {
        if (i_cross < tail_b) {
            if (jb::debug) {
                cout << " crossing part lasting longer than PSP sample length, very unlikely" << endl;
                cout << " dead old input " << tail_l-i_cross << endl;
            }
            for (i=i_cross; i<tail_b; i++) {
                input.junk_this();
            }
            i_start = tail_b;
        } else {
            // update for input that come before cross and lingers after cross;
            if (jb::debug) {
                cout << "cross update starting tail_b " << tail_b << endl;
                cout << " # lingering old inputs: " << i_cross - tail_b << endl;
            }
            for (i=tail_b; i<i_cross; i++) {
                if (jb::debug) {
                    cout << "lingering " << i << " < " << neuron.tin.size() << endl;
                }
                dt = tCross - input.t[i];
                if (jb::debug) {
                    assert(dt>0);
                }
                jb::getNear(neuroLib.idtRange,neuroLib.ndt,
                              dt, input.dTijr[i].r, input.dTijr[i].i, input.dTijr[i].j);
                input.dt[i] = tCross;
                input.cCross[i] = cross.nCross;
                input.Vijr[i] = cross.vCross.back();
                if (input.Vijr[i].j < neuroLib.fireCap[input.ID[i]][input.dTijr[i].i]) {
                    v = neuroLib.vRange[input.Vijr[i].i] + (neuroLib.vRange[input.Vijr[i].j] - neuroLib.vRange[input.Vijr[i].i]) * input.Vijr[i].r;
                    jb::getNear(neuroLib.vRange, neuroLib.fireCap[input.ID[i]][input.dTijr[i].i], v, input.Vijr[i].r, input.Vijr[i].i, input.Vijr[i].j);
                }
                input.tMax[i] = linear_interp_tMax(neuroLib.tMax, cross.vCross.back(), input.dTijr[i],input.ID[i]);
                if (jb::debug) {
                    cout << "   loop ended" << endl;
                }
            }
        }
        // update for new input during the cross that lingers after cross
        if (jb::debug) {
            cout << " # lingering new inputs: " << head-i_start+1 << endl;
            cout << input.t.size() << " == " << input.Vijr.size();
            if (input.t.size() != i_start) {
                cout << input.t.size() << " != " << i_start << endl;
                assert(input.t.size() == i_start);
            }
        }
        for (i=i_start; i<=head; i++) {
            //cout << "i " << i << " < " << neuron.tin.size() << " == " << neuron.inID.size() << endl;
            input.t.push_back(neuron.tin[i]/neuroLib.tstep);
            input.ID.push_back(neuron.inID[i]);
            input.dt.push_back(tCross);
            input.idt.push_back(0);
            dt = tCross - input.t[i];
            jb::getNear(neuroLib.idtRange, neuroLib.ndt, dt, r_, i_, j_);
            if (jb::debug) {
                assert(input.t.size() == i+1);
                if (dt < -1e-14) {
                    cout << tCross << " - " << input.t[i] << " = " << dt << endl;
                    assert(dt >= 0);
                }
                if (r_-1>1e-14) {
                    cout << tCross  << "-" << input.t[i] << endl;
                    cout << "tail_b " << tail_b << ", " << input.t[tail_b] << endl;
                    cout << i << " head " << head << "total " << neuron.tin.size() << endl;
                    cout << "dt " << dt << " i_ " << i_ << " j_ " << j_ << endl;
                    cout << "r_ " << r_ << endl;
                    assert(r_ <= 1);
                }
            }
            input.dTijr.push_back(IJR(i_,j_,r_));
            input.cCross.push_back(cross.nCross);
            input.Vijr.push_back(cross.vCross.back());
            if (input.Vijr[i].j < neuroLib.fireCap[input.ID[i]][input.dTijr[i].i]) {
                v = neuroLib.vRange[input.Vijr[i].i] + (neuroLib.vRange[input.Vijr[i].j] - neuroLib.vRange[input.Vijr[i].i]) * input.Vijr[i].r;
                jb::getNear(neuroLib.vRange, neuroLib.fireCap[input.ID[i]][input.dTijr[i].i], v, input.Vijr[i].r, input.Vijr[i].i, input.Vijr[i].j);
            }
            input.tMax.push_back(linear_interp_tMax(neuroLib.tMax, cross.vCross.back(), input.dTijr[i],input.ID[i]));
            input.Tmpijr.push_back(IJR(0,1,0)); 
            input.bir.push_back(BilinearRelationships(corrSize));
            if (afterCrossBehavior == 2) {
                for (j=tail_b; j<i; j++) {
                    add_relation_between_new_and_ith_input(input, i, j, neuroLib.idtRange, neuroLib.ndt);
                }
            } else {
                // new start for bilinear
                tail_b = head+1;
            }
        }
    }
    if (jb::debug) {
        input.assert_size();
    }
}

unsigned int nsyn_jBilinear(nNS &neuron, nNL &neuroLib, Input &input, jND &jnd, Cross &cross, double end_t, double ignore_t, size corrSize, vector<double> &tsp, double vC, double vB, int afterCrossBehavior, bool spikeShape) {
    size i, j, i_prior_cross;
    double tstep = neuroLib.tstep;
    unsigned int nc_old, nc = 0;
    double t_cross, vBack, tBack;
    double iend = end_t/tstep;
    double tol_tl = neuroLib.nt;
    double tol_tb = round((neuroLib.tol_tb-ignore_t)/tstep);
    cout << "linear corr length " << tol_tl << endl;
    cout << "bilinear corr length " << tol_tb << endl;
    cout << "total inputs " << neuron.tin.size() << endl;
    size tail_b = 0, tail_l = 0;
    size old_tail_b, old_tail_l;
    bool crossed, spiked;
    double tref = neuron.tRef/tstep;
    size i_, j_;
    double r_;
    bool debug = false;
    size ii;
    jb::getNear(neuroLib.vRange, neuroLib.nv, cross.v[0], r_, i_, j_);
    cross.vCross.push_back(IJR(i_,j_,r_));
    jnd.t.push_back(0);
    jnd.v.push_back(cross.v[0]);
    if (debug) {
        cout << "set v: " << jnd.v.back() << endl;
    }
    for (i=0;i<neuron.tin.size();i++) {
        input.t.push_back(neuron.tin[i]/tstep);
        if (jb::debug) {
            cout << " new input i = " << i << endl;
            assert(input.t.size() == i+1);
        }
        if (i>=0) {
            debug = false;
        } else {
            debug = false;
        }
        old_tail_l = tail_l;
        old_tail_b = tail_b;
        move_corr_window(neuron.tin, tail_l, input.t[i], tol_tl, tstep);
        move_corr_window(neuron.tin, tail_b, input.t[i], tol_tb, tstep);
        if (debug) {
            cout << " tail_l " << old_tail_l << " -> " << tail_l << endl;
            cout << " tail_b " << old_tail_b << " -> " << tail_b << endl;
        }
        crossed = update_vinit_of_new_input_check_crossing(input, cross, neuroLib, neuron, i, tail_l, tail_b, tol_tl, tol_tb, end_t, jnd, t_cross, vC, corrSize, debug);
        ii = 0;
        while (crossed) {
            ii++;
            if (t_cross < neuron.tin[i]/tstep + 1e-14) {
                if (jb::debug) {
                    cout << " unlikely, cross upon or before input" << endl;
                }
                tail_l = old_tail_l;
                tail_b = old_tail_b;
                i_prior_cross = i-1;
            } else {
                i_prior_cross = i;
            }
            if (jb::debug) {
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
                spiked = 1;
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
                for (j=0;j<tl; j++) {
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
                spiked = nc - nc_old;
                i = i_prior_cross + 1;
                if (i < neuron.inID.size()) {
                    while (neuron.tin[i]-1e-14 < tBack*tstep + 1e-14) {
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
            if (jb::debug) {
                cout << "set v: " << jnd.v.back() << endl;
                cout << " backed at " << tBack*neuroLib.tstep << " with v: " << vBack << endl;
                cout << " input from " << i_prior_cross << " to " << i;
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
            if (jb::debug) {
                cout << "crossed " << endl;
            }

            old_tail_b = tail_b;
            old_tail_l = tail_l;
            move_corr_window(neuron.tin, tail_b, tBack, tol_tb, tstep);
            tail_l = tail_b;
            if (debug) {
                cout << " tail_l " << old_tail_l << " -> " << tail_l << endl;
                cout << " tail_b " << old_tail_b << " -> " << tail_b << endl;
                old_tail_b = tail_b;
                old_tail_l = tail_l;
            }
            update_info_after_cross(input, neuroLib, cross, neuron, tBack, vBack, i_prior_cross, tail_l, tail_b, i, corrSize, afterCrossBehavior);
            if (jb::debug) {
                cout << " updated" << endl;
            }
            if (debug) {
                cout << " tail_l " << old_tail_l << " -> " << tail_l << endl;
                cout << " tail_b " << old_tail_b << " -> " << tail_b << endl;
            }
            crossed = check_crossing(input, neuroLib, cross, neuron, tol_tl, tol_tb, end_t, tail_l, tail_b, i, jnd, t_cross, vC, debug);
            if (jb::debug) {
                cout << "------------" << endl;
            }
        }
        if (jb::debug) {
            cout << " dead with " << i << ", v " << jnd.v.back() << " at " << jnd.t.back() << endl;
            input.print_this(i);
            input.assert_size();
        }
    }
    cout << " crossed " << cross.nCross << " times" << endl;
    return nc;
}
