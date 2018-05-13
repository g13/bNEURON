#include "jumpy_bilinear.h"

double interp_for_t_cross(double v_right, double v_left, double t_right, double t_left, size head, size tail_l, size tail_b, double tCross, double tol_tl, double tol_tb, nNL &neuroLib, Cross &cross, Input &input, double v0, double vtol, double vC, double &v, bool debug, bool dtSquare) {
    size ival = 0;
    double v1, v2;
    double t_cross1, t_cross2, t_cross;
    do {
        if (jb::debug2) {
            assert(vC >= v_left);
        }
        if (abs(t_right - t_left)<pow(2,-52)) {
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
        v1 = find_v_at_t(input, neuroLib, cross, head, tail_l, tail_b, t_cross1, tCross, tol_tl, tol_tb, v0, debug, dtSquare);
        if (jb::debug2) {
            cout << v_left << ", " << v1 << ", " << v_right << endl;
            ival++;
        }
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
            cout << " t_cross2 " << t_cross2 << endl;
            assert(t_cross2 >= t_left);
            assert(t_cross2 <= t_right);
        }
        // if t interval smaller than sample temp resolution, apply the quadratic interpolation as solution
        if (jb::debug2) {
            cout << " find v at " << t_cross2 << endl;
        }
        v2 = find_v_at_t(input, neuroLib, cross, head, tail_l, tail_b, t_cross2, tCross, tol_tl, tol_tb, v0, debug, dtSquare);
        if (jb::debug2) {
            ival++;
        }
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

bool check_crossing(Input &input, nNL &neuroLib, Cross &cross, nNS &neuron, double tol_tl, double tol_tb, double end_t, size tail_l, size tail_b, size head, jND &jnd, double &t_cross, double vC, bool debug, bool dtSquare) {
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
        v = find_v_at_t(input, neuroLib, cross, head, tail_l, tail_b, t, tCross, tol_tl, tol_tb, neuron.vReset, debug, dtSquare);
        if (v > vC) {
            lcrossed = true;
            break;    
        }
    }
    if (lcrossed) {
        // bilinear confirm crossed
        // perform linear (and parabola) interp iteration until tolerance reaches
        if (jb::debug) {
            cout << "v " << v << " > vThres, find t and v for cross" << endl;
        }
        if (fabs(v-vC)>neuron.vTol && t-jnd.t.back()>1) {
            if (jb::debug) {
                cout << "t_left " << jnd.t.back() << ", t_right " << t << endl;
                cout << "v_left " << jnd.v.back() << ", v_right " << v << endl;
            }
            t_cross = interp_for_t_cross(v, jnd.v.back(), t, jnd.t.back(), head, tail_l, tail_b, tCross, tol_tl, tol_tb, neuroLib, cross, input, neuron.vRest, neuron.vTol, vC, v_pass, debug, dtSquare);
            if (jb::debug) {
                cout << "t_left " << jnd.t.back() << " t_cross " << t_cross << ", t_right " << t << endl;
                cout << "v_left " << jnd.v.back() << " v_cross " << v_pass << ", v_right " << v << endl;
                assert(t_cross >= jnd.t.back());
            }
            jnd.t.push_back(t_cross);
            jnd.v.push_back(v_pass);
        } else {
            t_cross = t;
            jnd.t.push_back(t);
            jnd.v.push_back(v);    
        }
        if (jb::debug) {
            cout << "set v: " << jnd.v.back() << endl;
        }
        return true;
    }
    return false;
}

bool update_vinit_of_new_input_check_crossing(Input &input, Cross &cross, nNL &neuroLib, nNS &neuron, size head, size &tail_l, size &tail_b, size old_tail_l, size old_tail_b, double tol_tl, double tol_tb, double end_t, jND &jnd, double &t_cross, double vC, size corrSize, bool debug, bool dtSquare) {
    size i, j, idt;
    double dt;
    double tstep = neuroLib.tstep;
    double v_pass;
    double tCross = cross.tCross.back();
    double v;
    if (jb::debug) {
        cout << " adding " << head << " input t " << input.t[head] << " synapse ID " << neuron.inID[head] << endl;
    }
    input.ID.push_back(neuron.inID[head]);
    dt = input.t[head] - tCross;
    if (jb::debug) {
        if (dt < -pow(2,-52)) {
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
        cout << "tail_b: " << tail_b << " at t: " << input.t[tail_b] << endl;
        cout << "last cross at t: " << cross.tCross.back() << endl;
    }

    for (i=tail_l; i<head; i++) {
        if (debug) {
            dvtmp = 0;
            cout << "   linear input: " << i << " at t: " << input.t[i] << " <= dt: " << input.dt[i] << endl;
        }
        dt = input.t[head] - input.dt[i];
        idt = static_cast<size>(round(dt));
        if (i>=tail_b) {
            input.bir[i].ID.push_back(neuron.inID[head]);
        }
        add_input_i_contribution(i,idt,neuroLib,input,v);
        if (debug) {
            cout << "   contributing " << v-vtmp << " with synapse ID: " << input.ID[i] << " at V i: " << input.Vijr[i].i << " j: " << input.Vijr[i].j << " r: " << input.Vijr[i].r << endl;
            vtmp = v;
        }
        for (j=tail_b; j<i; j++) {
            if (debug){
                cout << "       +bilinear input: " << j << " at t: " << input.t[j] << " <= dt: " << input.dt[j] << endl;
            }
            add_input_i_j_bilinear_contribution(input, neuroLib, j, i,idt, v, debug);
            if (debug) {
                dvtmp += v-vtmp;
                cout << "       contributing " << v-vtmp << " with synapse ID " << input.ID[j] << " with dT i: " << input.bir[j].dTijr[i-j-1].i << " j: " << input.bir[j].dTijr[i-j-1].j << " r: " << input.bir[j].dTijr[i-j-1].r << endl;
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
        cout << "set v after bilinear: " << jnd.v.back() << endl;
    }
    // check if vmax at left bound very unlikely
    if (v > vC) {
        if (jb::debug2) {
            cout << " crossing upon input not very probable at low input rate " << endl;
        }
        jnd.t.pop_back();
        jnd.v.pop_back();
        // perform linear interp iteration until tolerance reaches
        if (jb::debug2) {
            cout << "t: " << jnd.t.back()  << " < " << input.t[head] << endl;
        }
        tail_l = old_tail_l;
        tail_b = old_tail_b;
        if (jb::debug) {
            cout << "t_left " << jnd.t.back() << ", t_right " << input.t[head] << endl;
            cout << "v_left " << jnd.v.back() << ", v_right " << v << endl;
        }
        t_cross = interp_for_t_cross(v, jnd.v.back(), input.t[head], jnd.t.back(), head-1, tail_l, tail_b, tCross, tol_tl, tol_tb, neuroLib, cross, input, neuron.vRest, neuron.vTol, vC, v_pass, debug, dtSquare);
        if (jb::debug) {
            cout << "t_left " << jnd.t.back() << " t_cross " << t_cross << ", t_right " << input.t[head] << endl;
            cout << "v_left " << jnd.v.back() << " v_cross " << v_pass << ", v_right " << v << endl;
            assert(t_cross >= jnd.t.back());
        }
        if (jb::debug2) {
            int exponent;
            frexp(input.t[head], &exponent);
            if (t_cross > input.t[head] + pow(2,exponent-52)) {
                assert(t_cross<=input.t[head]);
            }
        }
        jnd.t.push_back(t_cross);
        jnd.v.push_back(v_pass);
        if (jb::debug) {
            cout << "set v retrace 1 input: " << jnd.v.back() << endl;
        }
        // pop back the info of input[head]
        input.t.pop_back();
        input.ID.pop_back();
        for (i=tail_b; i<head; i++) {
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
        add_new_input_info(neuroLib.vRange, neuroLib.nv, input, cross, neuron.tin[head], neuroLib.tMax, neuroLib.sfireCap, v, corrSize, neuron.inID[head]);
        for (i=tail_b; i<head; i++) {
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
    return check_crossing(input, neuroLib, cross, neuron, tol_tl, tol_tb, end_t, tail_l, tail_b, head, jnd, t_cross, vC, debug, dtSquare);
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
    if (afterCrossBehavior == 0) {
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
                cout << " dead old input " << tail_b-i_cross << endl;
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
                if (input.Vijr[i].j >= neuroLib.sfireCap[input.dTijr[i].i][input.ID[i]] && input.dTijr[i].r > 0.1) {
                    // if sufficiently close to spike case
                    input.tMax[i] = neuroLib.tMax[input.Vijr[i].j][input.dTijr[i].i][input.ID[i]];
                } else {
                    input.tMax[i] = linear_interp_tMax(neuroLib.tMax, cross.vCross.back(), input.dTijr[i],input.ID[i]);
                }
                if (jb::debug) {
                    cout << "   loop ended" << endl;
                }
            }
        }
        // update for new input during the cross that lingers after cross
        if (jb::debug) {
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
            input.idt.push_back(0);
            dt = tCross - input.t[i];
            jb::getNear(neuroLib.idtRange, neuroLib.ndt, dt, r_, i_, j_);
            if (jb::debug) {
                cout << " new input i = " << i << endl;
                assert(input.t.size() == i+1);
                if (dt < -pow(2,-52)) {
                    cout << tCross << " - " << input.t[i] << " = " << dt << endl;
                    assert(dt >= 0);
                }
                if (r_> 1+pow(2,-52)) {
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
            if (input.Vijr[i].j >= neuroLib.sfireCap[input.dTijr[i].i][input.ID[i]] && input.dTijr[i].r > 0.1) {
                // if sufficiently close to spike case
                input.tMax.push_back(neuroLib.tMax[input.Vijr[i].j][input.dTijr[i].i][input.ID[i]]);
            } else {
                input.tMax.push_back(linear_interp_tMax(neuroLib.tMax, cross.vCross.back(), input.dTijr[i],input.ID[i]));
            }
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

unsigned int nsyn_jBilinear(Cell &cell, vector<vector<double>> &spikeTrain, vector<double> dendVclamp, nNS &neuron, nNL &neuroLib, Input &input, jND &jnd, Cross &cross, double end_t, double ignore_t, size corrSize, vector<double> &tsp, double vC, double vB, int afterCrossBehavior, bool spikeShape, bool dtSquare) {
    size i, j, i_prior_cross;
    double tstep = neuroLib.tstep;
    unsigned int nc_old, nc = 0;
    double t_cross, vBack, tBack, dt;
    double iend = end_t/tstep;
    size ndt = neuroLib.ndt;
    double tol_tl = neuroLib.nt;
    double tol_tb = neuroLib.idtRange[ndt-2] + neuroLib.nt - neuroLib.idtRange[ndt-1] - round(ignore_t/tstep);
    cout << "linear corr length " << tol_tl << endl;
    cout << "bilinear corr length " << tol_tb << endl;
    cout << "total inputs " << neuron.tin.size() << endl;
    size tail_b = 0, tail_l = 0;
    size old_tail_b, old_tail_l;
    bool crossed, spiked;
    double tref = neuron.tRef/tstep;
    vector<size> s0(neuroLib.nSyn,0);
    vector<size> s1(neuroLib.nSyn,0);
    size i_, j_;
    double r_;
    bool debug = false;
    size ii;
    jb::getNear(neuroLib.vRange, neuroLib.nv, cross.v[0], r_, i_, j_);
    cross.vCross.push_back(IJR(i_,j_,r_));
    jnd.t.push_back(0);
    jnd.v.push_back(cross.v[0]);
    if (debug) {
        cout << "set v initial: " << jnd.v.back() << endl;
    }
    for (i=0;i<neuron.tin.size();i++) {
        input.t.push_back(neuron.tin[i]/tstep);
        if (jb::debug) {
            cout << " new input i = " << i << endl;
            assert(input.t.size() == i+1);
        }
        if (i>=7) {
            debug = false;
        } else {
            debug = false;
        }
        old_tail_l = tail_l;
        old_tail_b = tail_b;
        dt = input.t[i]-cross.tCross.back();
        move_corr_window_i(input.t, tail_b, input.t[i], tol_tb);
        if (dt < tol_tb) { 
            move_corr_window_i(input.t, tail_l, input.t[i], tol_tb);
        } else {
            if (dt < tol_tl) {
                move_corr_window_i(input.t, tail_l, input.t[i], dt);
            } else{
                move_corr_window_i(input.t, tail_l, input.t[i], tol_tl);
            }
        }
        if (debug) {
            cout << " tail_l " << old_tail_l << " -> " << tail_l << endl;
            cout << " tail_b " << old_tail_b << " -> " << tail_b << endl;
        }
        crossed = update_vinit_of_new_input_check_crossing(input, cross, neuroLib, neuron, i, tail_l, tail_b, old_tail_l, old_tail_b, tol_tl, tol_tb, end_t, jnd, t_cross, vC, corrSize, debug, dtSquare);
        ii = 0;
        while (crossed) {
            ii++;
            if (jb::debug) {
                cout << " crossing " << endl;
                if (ii > 1) cout << " AGAIN?!!" << endl;
            }
            int exponent;
            frexp(neuron.tin[i]/tstep, &exponent);
            if (t_cross < neuron.tin[i]/tstep + pow(2,exponent-52)) {
                if (jb::debug) {
                    cout << " unlikely, cross upon or before input" << endl;
                }
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
                i = neuroAlter(neuron, neuroLib, cross, i_prior_cross, jnd, end_t, round(t_cross), tBack, vBack, tstep, tsp, vB, nc, cell, spikeTrain, s0, s1, dendVclamp);
                spiked = nc > nc_old;
            } else {
                tBack = t_cross+tref;
                vBack = neuron.vReset;
                size tl = tref;
                if (tBack > iend) {
                    tBack = iend;
                    tl = iend - ceil(t_cross);
                }
                for (j=0;j<tl;j++) {
                    cross.v.push_back(vBack);
                    cross.t.push_back(ceil(t_cross)+j);
                }
                tsp.push_back(t_cross*tstep);
                nc++;
                spiked = true;
                i = i_prior_cross;
                int exponent;
                frexp(tBack*tstep, &exponent);
                while (neuron.tin[i] < tBack*tstep + pow(2,exponent-52)) {
                    i++;
                    if (i==neuron.inID.size()) {
                        break;
                    }
                }
                i--;
            }
            jnd.t.push_back(tBack);
            jnd.v.push_back(vBack);
            if (jb::debug) {
                cout << "set v back: " << jnd.v.back() << endl;
                cout << " backed at " << tBack*tstep << " with v: " << vBack << endl;
                cout << " input from " << i_prior_cross << " to " << i << endl;
                cout << " input at " << neuron.tin[i] << endl;
                if (i+1 < neuron.tin.size()) {
                    cout << " next input " << neuron.tin[i+1] << endl;
                }
            }
            //if (vBack >= vC) {
            //    if (tBack*tstep < end_t) {
            //        cout << " it must crossing at t end, since vBack > vC" << endl;
            //        assert(tBack*tstep >= end_t);
            //    }
            //    break;
            //}
            if (tBack*tstep >= end_t) {
                cout << " ending while crossing " << endl;
                if (tBack*tstep > end_t) {
                    cout << tBack*tstep << " == " << end_t << endl;
                    assert(tBack*tstep == end_t);
                }
                break;
            }
            jb::getNear(neuroLib.vRange,neuroLib.nv,
                            vBack, r_, i_, j_);
            cross.vCross.push_back(IJR(i_,j_,r_));
            cross.nCross++;
            cross.iCross.push_back(cross.v.size());
            cross.tCross.push_back(tBack);
            if (jb::debug) {
                cout << cross.nCross << "th cross return" << endl;
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
                cout << " tail_l " << old_tail_l << " -> " << tail_l << endl;
                cout << " tail_b " << old_tail_b << " -> " << tail_b << endl;
            }
            crossed = check_crossing(input, neuroLib, cross, neuron, tol_tl, tol_tb, end_t, tail_l, tail_b, i, jnd, t_cross, vC, debug, dtSquare);
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
