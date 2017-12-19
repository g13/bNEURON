#include "jumpy_bilinear.h"

double interp_for_t_cross(double v_right, double v_left, double t_right, double t_left, size head, size tail_l, size tail_b, double tCross, double tol_tl, nNL &neuroLib, Cross &cross, Input &input, nNS &neuron, double &v, bool debug) {
    size ival = 0;
    double v1, v2, v0 = neuron.vRest;
    double t_cross1, t_cross2, t_cross;
    do {
        assert(neuron.vThres >= v_left);
        if (abs(t_right - t_left)<1e-10) {
            v = v_right;
            return t_right;
        }
        assert(t_right > t_left);
        t_cross1 = t_left + (neuron.vThres - v_left)/(v_right-v_left)*(t_right-t_left);
        //t_cross1 = ceil(t_cross1);
        assert(t_cross1 >= t_left);
        if (jb::debug2) {
            std::cout << " find v at " << t_cross1 << std::endl;
            std::cout <<"t: " << t_left << ", " << t_right << std::endl;
        }
        v1 = find_v_at_t(input, neuroLib, cross, head, tail_l, tail_b, t_cross1, tCross, tol_tl, v0, debug);
        if (jb::debug2) {
            std::cout << v_left <<  ", " << v1 << ", " << v_right << std::endl;
        }
        ival++;
        if (fabs(v1-neuron.vThres)<neuron.vTol) {
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
                std::cout << "temp resol reached " << std::endl;
            }
            break;
        }
        // solve for a, b of f(t)-v_left = a(t-t_left)^2 + b(t-t_left);
        // solve for t when a(t-t_left)^2 + b(t-t_left) = (vThres - v_left)
        t_cross2 = parabola(t_left,v_left,t_right,v_right,t_cross1, v1, neuron.vThres);
        //t_cross2 = ceil(t_cross2);
        if (jb::debug2) {
            std::cout << " v: " << v_left << ", " << v1 << ", " << v_right << std::endl;
            std::cout << " t: " << t_left << ", " << t_cross1 << ", " << t_right << std::endl;
            std::cout << "t_cross2 " << t_cross2 << std::endl;
            assert(t_cross2 >= t_left);
            assert(t_cross2 <= t_right);
        }
        // if t interval smaller than sample temp resolution, apply the quadratic interpolation as solution
        if (jb::debug2) {
            std::cout << " find v at " << t_cross2 << std::endl;
        }
        v2 = find_v_at_t(input, neuroLib, cross, head, tail_l, tail_b, t_cross2, tCross, tol_tl, v0, debug);
        //if (jb::debug2) {
            ival++;
        //}
        if (fabs(v2-neuron.vThres)<neuron.vTol) {
            v = v2;
            t_cross = t_cross2;
            break;
        }
        if (jb::debug2) {
            std::cout << " v: " << v_left << ", " << v1 << ", " << v2 << ", " << v_right << std::endl;
            std::cout << " t: " << t_left << ", " << t_cross1 << ", " << t_cross2 << ", " << t_right << std::endl;
        }
        if (v2>v1) {
            getLR(v_left,v1,v2,v_right,t_left,t_cross1,t_cross2,t_right, neuron.vThres);
        } else {
            getLR(v_left,v2,v1,v_right,t_left,t_cross2,t_cross1,t_right, neuron.vThres);
        } 
        if (jb::debug2) {
            std::cout << " v: " << v_left << ", " << v_right << std::endl;
            std::cout << " t: " << t_left << ", " << t_right << std::endl;
        }
    } while (true);
    //if (jb::debug2) {
        std::cout << ival << " evaluations" << std::endl;
    //}
    return t_cross;
}

bool check_crossing(Input &input, nNL &neuroLib, Cross &cross, nNS &neuron, double tol_tl, double tol_tb, double end_t, size tail_l, size tail_b, size head, jND &jnd, double &t_cross, bool debug) {
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
        if (t <= jnd.t.back() || neuroLib.ei[input.ID[i]] || input.inTref[i]) {
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
            v = add_vinit_contribution(neuroLib, cross.vCross.back(), dt);
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
            if (!input.inTref[j]) {
                add_input_i_contribution(j, it, neuroLib, input, v);
            }
        }
        if (v > vmax) {
            vmax = v;
            tmax = t; 
            vi_tail_max = vi_tail;
        }
    }
    if (vmax > neuron.vThres) {
        if (jb::debug2) {
            std::cout << "linear vmax " << vmax <<  " > vThres" << std::endl;
        }
        // perform linear interp iteration until tolerance reaches
        if (jb::debug) {
            std::cout << "l+b vmax > vThres, find t and v for cross" << std::endl;
        }
        if (fabs(vmax-neuron.vThres)>neuron.vTol && tmax-jnd.t.back()>1) {
            t_cross = interp_for_t_cross(vmax, jnd.v.back(), tmax, jnd.t.back(), head, tail_l, tail_b, tCross, tol_tl, neuroLib, cross, input, neuron, v_pass, debug);
            assert(t_cross > jnd.t.back());
            jnd.t.push_back(t_cross);
            jnd.v.push_back(v_pass);
        } else {
            t_cross = tmax;
            jnd.t.push_back(tmax);    
            jnd.v.push_back(vmax);    
        }
        return true;
    }
    return false;
}
bool update_vinit_of_new_input_check_crossing(Input &input, Cross &cross, nNL &neuroLib, nNS &neuron, size head, size tail_l, size tail_b, double tol_tl, double tol_tb, double end_t, jND &jnd, double &t_cross, size corrSize, bool debug) {
    long i;
    size j;
    double dt;
    double tstep = neuroLib.tstep;
    double v_pass;
    double ir;
    size idt;
    double tCross = cross.tCross.back();
    double v;
        // if tstep1 isn't 0 then collect v with tstep1 until the time of the new input
    //if (jb::debug) {
    //    cout << " adding " << head << " input t " << input.t[head] << " synapse ID " << neuron.inID[head] << endl;
    //}
    input.ID.push_back(neuron.inID[head]);
    dt = input.t[head] - tCross;
    if (jb::debug) {
        if (dt < -1e-10 ) {
            std::cout << dt << " < 0 " << std::endl;
            std::cout << head << std::endl;
            std::cout << neuron.tin[head] << std::endl;
            std::cout << input.t[head] << std::endl;
            std::cout << tCross*tstep << std::endl;
            assert(dt>=0);
        }
    }
    // initialize to 
    if (dt < tol_tl) {
        v = add_vinit_contribution(neuroLib, cross.vCross.back(), dt);
    } else {
        v = neuron.vRest;
    }
    if (head-1 < tail_l && tail_l != 0) cout << "first v after cross " << v << endl;
    //cout << " tailing input " << tail_l << endl;
    double vtmp = v;
    double dvtmp;
    if (debug) {
        cout << " last v " << jnd.v.back() << endl;
        cout << " start with leak" << v << endl;
    }
    for (i=head-1; i>=tail_l; i--) {
        if (debug) {
            dvtmp = 0;
        }
        dt = input.t[head] - input.dt[i];
        input.bir[i].idt.push_back(static_cast<size>(round(dt)));
        input.bir[i].ID.push_back(neuron.inID[head]);
        if (!input.inTref[i]) {
            add_input_i_contribution(i,input.bir[i].idt.back(),neuroLib,input,v);
            if (debug) {
                cout << " " << neuron.ei[input.ID[i]] << "-" << i << " input contributing " << v-vtmp << ", synapse ID " << input.ID[i] << endl;
                vtmp = v;
            }
            if (i<tail_b) continue;
            for (j=head-1; j>i; j--) {
                ir = head-j-1;
                if (!input.inTref[j]) {
                    add_input_i_j_bilinear_contribution(input, neuroLib, i, j, input.bir[j].idt[ir], v, debug);
                    if (debug) {
                        dvtmp += v-vtmp;
                        cout << "    +" << neuron.ei[input.ID[j]] << "-" << j << " bilinear contributing " << v-vtmp << ", synapse ID " << input.ID[j] << " with dti " << input.bir[i].dTijr[j-i-1].i << ", r" << input.bir[i].dTijr[j-i-1].r << " at vi " << input.Vijr[j].i << ", r" << input.Vijr[j].r << endl;
                        vtmp = v;
                    }
                }
            }
            if (debug) {
                cout << " bilinear total: " << dvtmp << endl;
                cout << " v: " << v << endl;
            }
        }
    }
    // update head's v and tmax
    jnd.t.push_back(input.t[head]);
    jnd.v.push_back(v);
    //if (jb::debug) {
    //    cout << "v_" << head << " = " << jnd.v.back() << endl;
    //}
    // check if vmax at left bound very unlikely
    if (v > neuron.vThres) {
        if(jb::debug2) {
            std::cout << " crossing upon input not very probable at low input rate " << std::endl;
        }
        jnd.t.pop_back();
        jnd.v.pop_back();
        // perform linear interp iteration until tolerance reaches
        cout << "t: " << jnd.t.back()  << " < " << input.t[head] << endl;
        t_cross = interp_for_t_cross(v, jnd.v.back(), input.t[head], jnd.t.back(), head-1, tail_l, tail_b, tCross, tol_tl, neuroLib, cross, input, neuron, v_pass, debug);
        if (t_cross-input.t[head] > 1e-10) {
        
            assert(t_cross<=input.t[head]);
        }
        jnd.t.push_back(t_cross);
        jnd.v.push_back(v_pass);
        input.t.pop_back();
        input.ID.pop_back();
        for (i=head-1; i>=tail_l; i--) {
            input.bir[i].idt.pop_back();
            input.bir[i].ID.pop_back();
        }
        if (input.t.size() != head) {
            cout << "before check " << input.t.size() << " != " << head << endl;
            assert(input.t.size() == head);
        }
        return true;
    } else {
        add_new_input_info(neuroLib, input, cross, neuron.tin[head], neuroLib.tMax, v, corrSize, neuron.inID[head]);
        for (i=head-1; i>=tail_l; i--) {
            add_relation_between_new_and_ith_input(input, head, i, neuroLib);
        }
        if (input.t.size() != head+1) {
            cout << "before check " << input.t.size() << " != " << head+1 << endl;
            assert(input.t.size() == head+1);
        }
    }
    // check all tmax after head and find the upper limit of vmax
    return check_crossing(input, neuroLib, cross, neuron, tol_tl, tol_tb, end_t, tail_l, tail_b, head, jnd, t_cross, debug);
}

void update_info_after_cross(Input &input, nNL &neuroLib, Cross &cross, nNS &neuron, double tCross, double vCross, size i_prior, size tail, size head, size corrSize, std::vector<double>&tsp, bool spiked) {
    size i_, j_;
    double r_;
    double dt, v;
    size i, j;
    size idt;
    // inputs that only matters during cross, can fill junk with those
    size i_cross = i_prior + 1;
    size i_start = i_cross;
    if (input.t.size() != i_start) {
        cout << "before update " << input.t.size() << " != " << i_start << endl;
        assert(input.t.size() == i_start);
    }
    if (i_cross < tail) {
        //if (jb::debug) {
            std::cout << " crossing part lasting longer than PSP sample length, very unlikely" << std::endl;
        //}
        cout << " dead old input " << tail-i_cross << endl;
        for (i=i_cross; i<tail; i++) {
            input.junk_this();
        }
        i_start = tail;
    } else {
        // update for input that come before cross and lingers after cross;
        if (jb::debug) {
          cout << "cross update starting tail " << tail << endl;
          cout << " linger old inputs: " << i_cross - tail << endl;
        }
        for (i=tail; i<i_cross; i++) {
            if (jb::debug) {
                cout << "lingering " << i << " < " << neuron.tin.size() << endl;
            }
            dt = tCross - input.t[i];
            assert(dt>0);
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
            //cout << "   loop ended" << endl;
        }
    }
    // update for new input that lingers after cross
    if (jb::debug) {
        cout << " lingering new inputs: " << head-i_start+1 << endl;
        cout << input.t.size() << " == " << input.Vijr.size();
    }
    if (input.t.size() != i_start) {
        cout << input.t.size() << " != " << i_start << endl;
        assert(input.t.size() == i_start);
    }
    for (i=i_start; i<=head; i++) {
        //std::cout << "i " << i << " < " << neuron.tin.size() << " == " << neuron.inID.size() << std::endl;
        input.t.push_back(neuron.tin[i]/neuroLib.tstep);
        assert(input.t.size() == i+1);
        if (spiked && neuron.tin[i] > tsp.back()-1) {
            input.inTref.push_back(0);
        } else {
            input.inTref.push_back(0);
        }
        input.ID.push_back(neuron.inID[i]);
        if (jb::debug) {
            cout << head << " >= " << i << " > " << i_start << endl;
            cout << i << " == " << input.t.size()-1 << endl;
            cout << input.t[i] << " == " << input.t.back() << endl;
        }
        if (input.t[i] != input.t.back()) {
            std::cout << i << " == " << input.t.size() -1 << std::endl;
            std::cout << input.t[i]*neuroLib.tstep << " == " << neuron.tin[i]/neuroLib.tstep << std::endl;
            std::cout << input.t[i]*neuroLib.tstep << " == " << input.t[i-1]*neuroLib.tstep << std::endl;
            assert(input.t[i] == input.t.back());
        }
        if ( input.t[i] - tCross > 1e-10) {
            cout << "input.t " << i << " = " << input.t[i] << " < " << tCross << endl;
            assert(tCross>=input.t[i]);
        }
        input.dt.push_back(tCross);
        dt = tCross - input.t[i];
        if (dt < -1e-10) {
            cout << tCross << " - " << input.t[i] << " = " << dt << endl;
            assert(dt >= 0);
        }
        jb::getNear(neuroLib.idtRange, neuroLib.ndt,
                       dt, r_, i_, j_);
        if (r_-1>1e-10) {
            cout << tCross  << "-" << input.t[i] << endl;
            cout << "tail " << tail << ", " << input.t[tail] << endl;
            cout << i << " head " << head << "total " << neuron.tin.size() << endl;
            cout << "dt " << dt << " i_ " << i_ << " j_ " << j_ << endl;
            cout << "r_ " << r_ << endl;
            assert(r_ <= 1);
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
        for (j=tail; j<i; j++) {
            dt = input.t[i] - input.dt[j];
            input.bir[j].idt.push_back(static_cast<size>(round(dt)));
            add_relation_between_new_and_ith_input(input, i, j, neuroLib);
        }
        cout << "   bir pushed" << endl;
    }
    input.assert_size();
}

unsigned int nsyn_jBilinear(nNS &neuron, nNL &neuroLib, Input &input, jND &jnd, Cross &cross, double end_t, std::vector<double> &v, size corrSize, std::vector<double> &tsp, double vStop, vector<bool> &ei) {
    size i, j, i_prior_cross;
    double iend = end_t/neuroLib.tstep;
    unsigned int nc_old, nc = 0;
    double t_cross, vBack, tBack;
    double tstep = neuroLib.tstep;
    double tol_tl = round(neuroLib.tol_tl/tstep);
    double tol_tb = round(neuroLib.tol_tb/tstep);
        std::cout << "linear corr length " << tol_tl << std::endl;
        std::cout << "bilinear corr length " << tol_tb << std::endl;
        std::cout << "total inputs " << neuron.tin.size() << std::endl;
    size tail_b = 0;
    size tail_l = 0;
    size old_tail_b;
    size old_tail_l;
    bool crossed, spiked;
    double tref = neuron.tRef/tstep;
    size i_, j_;
    double r_;
    bool debug = false;
    size ii;
    jb::getNear(neuroLib.vRange,neuroLib.nv,
                    cross.v[0], r_, i_, j_);
    cross.vCross.push_back(IJR(i_,j_,r_));
    jnd.t.push_back(0);
    jnd.v.push_back(cross.v[0]);
    for (i=0;i<neuron.tin.size();i++) {
        if (i==2) {
            debug = false;
            //cout << " 39 ith input, just after jump" << endl;
            //cout << endl;
        } else  {
            debug = false;
            //cout << i << " th input" << endl;
        }
        input.assert_size();
        input.t.push_back(neuron.tin[i]/tstep); 
        assert(input.t.size() == i+1);
        old_tail_l = tail_l;
        old_tail_b = tail_b;
        move_corr_window(neuron.tin, tail_l, input.t[i], tol_tl, tstep);
        move_corr_window(neuron.tin, tail_b, input.t[i], tol_tb, tstep);
        crossed = update_vinit_of_new_input_check_crossing(input, cross, neuroLib, neuron, i, tail_l, tail_b, tol_tl, tol_tb, end_t, jnd, t_cross, corrSize, debug);
        ii = 0;
        while (crossed) {
            ii++;
            if (t_cross < neuron.tin[i]/tstep + 1e-10) {
                cout << " cross upon or before input" << endl;
                tail_l = old_tail_l;
                tail_b = old_tail_b;
                //reverse_corr_window(neuron.tin, tail_l, t_cross, tol_tl, tstep);
                i_prior_cross = i-1;
            } else {
                i_prior_cross = i;
            }
            if (jb::debug) {
                std::cout << " vCross = " << jnd.v.back();
                std::cout << " tCross = " << jnd.t.back() << std::endl;
                std::cout << " input " << i_prior_cross << ", t" << input.t[i_prior_cross]*neuroLib.tstep << std::endl;
            }
            if (input.t.size() != i_prior_cross+1) {
                cout << "before alter " << input.t.size() << " != " << i_prior_cross+1 << endl;
                assert(input.t.size() == i_prior_cross+1);
            }
            tBack = t_cross + tref;
            if (tBack>iend) {
                tBack = iend;
            }
            vBack = neuron.vRest;
            double tmpTsp = t_cross*neuroLib.tstep+neuron.tRef/2;
            if (tmpTsp <= end_t) {
                tsp.push_back(tmpTsp);
            }
            nc_old = nc;
            nc++;
            spiked = nc - nc_old;
            i = i_prior_cross + 1;
            if (i < neuron.inID.size()) {
                while (neuron.tin[i]-1e-10 < tBack*tstep+1e-10) {
                    i++;
                    if (i == neuron.inID.size()) break;
                }
            }
            i--;
            jnd.t.push_back(tBack);
            jnd.v.push_back(vBack);
            if (vBack >= neuron.vThres) break;
            if (jb::debug) {
                std::cout << " backed at " << tBack*neuroLib.tstep << std::endl;
                std::cout << " input from " << i_prior_cross << " to " << i << std::endl;
                std::cout << " input at " << neuron.tin[i] << std::endl;
            }
            jb::getNear(neuroLib.vRange,neuroLib.nv,
                            vBack, r_, i_, j_);
            cross.vCross.push_back(IJR(i_,j_,r_));
            cross.nCross++;
            cross.iCross.push_back(cross.v.size());
            cross.tCross.push_back(tBack);
            cout << "crossed " << endl;

            old_tail_b = tail_b;
            old_tail_l = tail_l;
            move_corr_window(neuron.tin, tail_l, tBack, tol_tl, tstep);
            move_corr_window(neuron.tin, tail_b, tBack, tol_tb, tstep);
            update_info_after_cross(input, neuroLib, cross, neuron, tBack, vBack, i_prior_cross, tail_b, i, corrSize, tsp, spiked);
            cout << " updated" << endl;
            tail_b = i + 1;
            tail_l = i + 1;
            crossed = check_crossing(input, neuroLib, cross, neuron, tol_tl, tol_tb, end_t, tail_l, tail_b, i, jnd, t_cross, debug);
            cout << "------------" << endl;
        }
        cout << " dead with " << i << ", v " << jnd.v.back() << " at " << jnd.t.back() << endl;
        input.print_this(i);
    }
    cout << " crossed " << cross.nCross << " times" << endl;
    return nc;
}
