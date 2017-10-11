#ifndef nJL_H 
#define nJL_H
#include <cmath>
#include <vector>
#include "typedefs.h"
#include "nNeuroLib.h"
#include "nNeuroSt2.h"
#include "neuroTest.h"
#include "nsynJtest2.h"
namespace jl{
    bool debug = true;
    bool debug2 = true;
    inline double find_v_at_t(Input &input, nNL &neuroLib, Cross &cross, size head, size tail_l, double t, double tCross, double tol_tl, double v) {
        long i;
        size j, it;
        double dt;
        dt = t - tCross;
        if (dt < tol_tl) {
            v = add_vinit_contribution(neuroLib, cross.vCross.back(), dt);
        }
        for (i=head; i>=tail_l; i--) {
            dt = t - input.dt[i];
            if (dt > tol_tl) break;
            it = static_cast<size>(round(dt));
            if (!input.inTref[i]) {
                add_input_i_contribution(i,it,neuroLib,input,v);
            }
        }
        return v;
    }

    inline double interp_for_t_cross(double v_right, double v_left, double t_right, double t_left, size head, size tail_l, double tCross, double tol_tl, nNL &neuroLib, Cross &cross, Input &input, nNS &neuron, double &v) {
        size ival = 0;
        double v1, v2, v0 = neuron.vRest;
        double t_cross1, t_cross2, t_cross;
        do {
            assert(neuron.vThres >= v_left);
            if (abs(t_right - t_left)<1e-10) {
                v = v_right;
                return t_right;
            }
            if (t_right < t_left) {
                std::cout << t_left << " < " << t_right << std::endl;
                assert(t_right >= t_left);
            }
            t_cross1 = t_left + (neuron.vThres - v_left)/(v_right-v_left)*(t_right-t_left);
            t_cross1 = ceil(t_cross1);
            assert(t_cross1 >= t_left);
            if (debug2) {
                std::cout << " find v at " << t_cross1 << std::endl;
                std::cout <<"t: " << t_left << ", " << t_right << std::endl;
            }
            v1 = find_v_at_t(input, neuroLib, cross, head, tail_l, t_cross1, tCross, tol_tl, v0);
            if (debug2) {
                std::cout << "v1 " << v1 << std::endl;
                ival++;
            }
            if (fabs(v1-neuron.vThres)<neuron.v_tol) {
                v = v1;
                t_cross = t_cross1;
                break;
            }

            if (t_cross1-t_left <= 1 || t_right-t_cross1 <=1) {
                if (t_cross1-t_left <=1 ) {
                    v = v1; 
                } else {
                    v = v_right;
                }
                t_cross = t_cross1;
                if (debug2) {
                    std::cout << "temp resol reached " << std::endl;
                }
                break;
            }
            // solve for a, b of f(t)-v_left = a(t-t_left)^2 + b(t-t_left);
            // solve for t when a(t-t_left)^2 + b(t-t_left) = (vThres - v_left)
            t_cross2 = parabola(t_left,v_left,t_right,v_right,t_cross1, v1, neuron.vThres);
            t_cross2 = ceil(t_cross2);
            if (debug2) {
                std::cout << " v: " << v_left << ", " << v1 << ", " << v_right << std::endl;
                std::cout << " t: " << t_left << ", " << t_cross1 << ", " << t_right << std::endl;
                std::cout << "t_cross2 " << t_cross2 << std::endl;
                assert(t_cross2 >= t_left);
                assert(t_cross2 <= t_right);
            }
            // if t interval smaller than sample temp resolution, apply the quadratic interpolation as solution
            if (debug2) {
                std::cout << " find v at " << t_cross2 << std::endl;
            }
            v2 = find_v_at_t(input, neuroLib, cross, head, tail_l, t_cross2, tCross, tol_tl, v0);
            if (debug2) {
                ival++;
            }
            if (fabs(v2-neuron.vThres)<neuron.v_tol) {
                v = v2;
                t_cross = t_cross2;
                break;
            }
            if (debug2) {
                std::cout << " v: " << v_left << ", " << v1 << ", " << v2 << ", " << v_right << std::endl;
                std::cout << " t: " << t_left << ", " << t_cross1 << ", " << t_cross2 << ", " << t_right << std::endl;
            }
            if (v2>v1) {
                getLR(v_left,v1,v2,v_right,t_left,t_cross1,t_cross2,t_right, neuron.vThres);
            } else {
                getLR(v_left,v2,v1,v_right,t_left,t_cross2,t_cross1,t_right, neuron.vThres);
            } 
            if (debug2) {
                std::cout << " v: " << v_left << ", " << v_right << std::endl;
                std::cout << " t: " << t_left << ", " << t_right << std::endl;
            }
        } while (true);
        if (debug2) {
            std::cout << ival << " evaluations" << std::endl;
        }
        return t_cross;
    }

    bool check_crossing(Input &input, nNL &neuroLib, Cross &cross, nNS &neuron, double tol_tl, double end_t, size tail_l, size head, jND &jnd, double &t_cross) {
        // check all tmax after head and find the upper limit of vmax
        size i,j;
        double t, dt;
        size nt = static_cast<size>(end_t/neuroLib.tstep)+1;
        double v, v_pass; 
        double vmax = jnd.v.back();
        double tmax = jnd.t.back();
        size it;
        double tCross = cross.tCross.back();
        for (i=tail_l; i<=head; i++) {
            t = input.tMax[i]+input.t[i];
            if (t <= jnd.t.back() || !neuroLib.ei[input.ID[i]] || input.inTref[i]) {
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
                dt = t - input.dt[j];
                if (dt > tol_tl) {
                    continue;
                } 
                it = static_cast<size>(round(dt));
                add_input_i_contribution(j, it, neuroLib, input, v);
            }
            if (v > vmax) {
                vmax = v;
                tmax = t; 
                if (vmax > neuron.vThres) {
                    break;    
                }
            }
        }
        if (vmax > neuron.vThres) {
            // perform linear interp iteration until tolerance reaches
            if (debug) {
                std::cout << "vmax " << vmax << " > vThres, find t and v for cross" << std::endl;
            }
            if (fabs(vmax-neuron.vThres)>neuron.v_tol) {
                t_cross = interp_for_t_cross(vmax, jnd.v.back(), tmax, jnd.t.back(), head, tail_l, tCross, tol_tl, neuroLib, cross, input, neuron, v_pass);
                cout << "t_left " << jnd.t.back() << " t_cross " << t_cross << ", t_right " << tmax << endl;
                cout << "v_left " << jnd.v.back() << " v_cross " << v_pass << ", t_right " << vmax << endl;
                assert(t_cross >= jnd.t.back());
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

    inline bool update_vinit_of_new_input_check_crossing(Input &input, Cross &cross, nNL &neuroLib, nNS &neuron, size head, size tail_l, double tol_tl, double end_t, jND &jnd, double &t_cross, size corrSize) {
        long i;
        double dt;
        double tstep = neuroLib.tstep;
        double v_pass;
        double ir;
        size idt;
        double tCross = cross.tCross.back();
        double v;
            // if tstep1 isn't 0 then collect v with tstep1 until the time of the new input
        input.ID.push_back(neuron.inID[head]);
        dt = input.t[head] - tCross;
        if (debug) {
            if (dt < 0 ) {
                std::cout << head << " < " << input.t.size() << std::endl;
                std::cout << input.t[head]*tstep << std::endl;
                std::cout << tCross*tstep << std::endl;
                assert(dt>0);
            }
        }
        // initialize to 
        if (dt < tol_tl) {
            v = add_vinit_contribution(neuroLib, cross.vCross.back(), dt);
        } else {
            v = neuron.vRest;
        }
        for (i=head-1; i>=tail_l; i--) {
            dt = input.t[head] - input.dt[i];
            input.bir[i].idt.push_back(static_cast<size>(round(dt)));
            input.bir[i].ID.push_back(neuron.inID[head]);
            ir = head-i-1;
            //std::cout <<  "input t "<<input.t[i] << std::endl;
            assert(input.bir[i].idt.size() == ir+1);
            if (!input.inTref[i]) {
                add_input_i_contribution(i,input.bir[i].idt[ir],neuroLib,input,v);
            }
        }
        // update head's v and tmax
        jnd.t.push_back(input.t[head]);
        jnd.v.push_back(v);
        // check if vmax at left bound very unlikely
        if (v > neuron.vThres) {
            if(debug2) {
                std::cout << " crossing upon input not very probable at low input rate " << std::endl;
            }
            jnd.t.pop_back();
            jnd.v.pop_back();
            // perform linear interp iteration until tolerance reaches
            t_cross = interp_for_t_cross(v, jnd.v.back(), input.t[head], jnd.t.back(), head-1, tail_l, tCross, tol_tl, neuroLib, cross, input, neuron, v_pass);
            if (t_cross-input.t[head] > 0) {
                assert(t_cross<=input.t[head]);
            }
            jnd.t.push_back(t_cross);
            jnd.v.push_back(v_pass);
            input.t.pop_back();
            input.ID.pop_back();
            for (i=head-1; i>=tail_l; i--) {
                input.bir[i].idt.pop_back();
            }
            return true;
        } else {
            // add new input's info 
            add_new_input_info(neuroLib, input, cross, neuron.tin[head], neuroLib.tMax, v, corrSize, neuron.inID[head]);
            if (input.t.size() != head+1) {
                cout << "before check " << input.t.size() << " != " << head+1 << endl;
                assert(input.t.size() == head+1);
            }
        }
        // check all tmax after head and find the upper limit of vmax
        return check_crossing(input, neuroLib, cross, neuron, tol_tl, end_t, tail_l, head, jnd, t_cross);
    }

    inline void update_info_after_cross(Input &input, nNL &neuroLib, Cross &cross, nNS &neuron, double tCross, double vCross, size i_prior, size tail, size head, size corrSize, std::vector<double>&tsp, bool spiked) {
        size i_, j_;
        double r_;
        double dt, v;
        size i, j;
        size idt;
        // inputs that only matters during cross, can fill junk with those
        size i_cross = i_prior + 1;
        size i_start = i_cross;
        if (i_cross < tail) {
            //if (debug) {
                std::cout << " crossing part lasting longer than PSP sample length, very unlikely" << std::endl;
            //}
            cout << " dead old input " << tail-i_cross << endl;
            for (i=i_cross; i<tail; i++) {
                input.junk_this();
            }
            i_start = tail;
        } else {
            // update for input that come before cross and lingers after cross;
            if (debug) {
              cout << "cross update starting tail " << tail << endl;
              cout << " linger old inputs: " << i_cross - tail << endl;
            }
            for (i=tail; i<i_cross; i++) {
                if (debug) {
                    cout << "lingering " << i << " < " << neuron.tin.size() << endl;
                }
                dt = tCross - input.t[i];
                assert(dt>0);
                jumpy::getNear(neuroLib.idtRange,neuroLib.ndt,
                              dt, input.dTijr[i].r, input.dTijr[i].i, input.dTijr[i].j);
                input.dt[i] = tCross;
                input.cCross[i] = cross.nCross;
                input.Vijr[i] = cross.vCross.back();
                //if (input.Vijr[i].j < neuroLib.fireCap[input.ID[i]][input.dTijr[i].i]) {
                    v = neuroLib.vRange[input.Vijr[i].i] + (neuroLib.vRange[input.Vijr[i].j] - neuroLib.vRange[input.Vijr[i].i]) * input.Vijr[i].r;
                    jumpy::getNear(neuroLib.vRange, neuroLib.fireCap[input.ID[i]][input.dTijr[i].i], v, input.Vijr[i].r, input.Vijr[i].i, input.Vijr[i].j);
                //}
                input.tMax[i] = linear_interp_tMax(neuroLib.tMax, cross.vCross.back(), input.dTijr[i],i);
                //cout << "   loop ended" << endl;
            }
        }
        // update for new input that lingers after cross
        if (debug) {
            cout << " lingering new inputs: " << head-i_start+1 << endl;
            cout << input.t.size() << " == " << input.Vijr.size() << endl;
        }
        if (input.t.size() != i_start) {
            cout << input.t.size() << " != " << i_start << endl;
            assert(input.t.size() == i_start);
        }
        for (i=i_start; i<=head; i++) {
            input.t.push_back(neuron.tin[i]/neuroLib.tstep);
            assert(input.t.size() == i+1);
            if (spiked && neuron.tin[i] > tsp.back()-1) {
                input.inTref.push_back(0);
            } else {
                input.inTref.push_back(0);
            }
            input.ID.push_back(neuron.inID[i]);
            if (debug) {
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
            jumpy::getNear(neuroLib.idtRange, neuroLib.ndt,
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
                jumpy::getNear(neuroLib.vRange, neuroLib.fireCap[input.ID[i]][input.dTijr[i].i], v, input.Vijr[i].r, input.Vijr[i].i, input.Vijr[i].j);
            }
            input.tMax.push_back(linear_interp_tMax(neuroLib.tMax, cross.vCross.back(), input.dTijr[i],i));
            input.Tmpijr.push_back(IJR(0,1,0)); 
            input.bir.push_back(BilinearRelationships(corrSize));
            for (j=tail; j<i; j++) {
                dt = input.t[i] - input.dt[j];
                input.bir[j].idt.push_back(static_cast<size>(round(dt)));
            }
        }
        input.assert_size();
    }
}

unsigned int nsyn_jLinear(nNS &neuron, nNL &neuroLib, Input &input, jND &jnd, Cross &cross, double end_t, std::vector<double> &v, size corrSize, std::vector<double> &tsp, double vStop, Cell &cell, vector<bool> &ei) {
    size i, i_prior_cross;
    double iend = end_t/neuroLib.tstep;
    unsigned int nc_old, nc = 0;
    double t_cross, vBack, tBack;
    double tstep = neuroLib.tstep;
    double tol_tl = round(neuroLib.tol_tl/tstep);
        std::cout << "linear corr length " << tol_tl << std::endl;
        std::cout << "total inputs " << neuron.tin.size() << std::endl;
    size tail_l = 0;
    size old_tail_l;
    bool crossed, spiked;
    double tref = neuron.tref/tstep;
    size i_, j_;
    double r_;
    size ii;
    jumpy::getNear(neuroLib.vRange,neuroLib.nv,
                    cross.v[0], r_, i_, j_);
    cross.vCross.push_back(IJR(i_,j_,r_));
    jnd.t.push_back(0);
    jnd.v.push_back(cross.v[0]);
    for (i=0;i<neuron.tin.size();i++) {
        input.assert_size();
        input.t.push_back(neuron.tin[i]/tstep); 
        assert(input.t.size() == i+1);
        old_tail_l = tail_l;
        move_corr_window(neuron.tin, tail_l, input.t[i], tol_tl,tstep);
        crossed = jl::update_vinit_of_new_input_check_crossing(input, cross, neuroLib, neuron, i, tail_l, tol_tl, end_t, jnd, t_cross, corrSize);
        ii = 0;
        while (crossed) {
            ii++;
            cout << " crossing " << endl;
            if (ii > 1) cout << " AGAIN?!!" << endl;
            if (t_cross <= neuron.tin[i]/tstep) {
                tail_l = old_tail_l;
                //reverse_corr_window(neuron.tin, tail_l, t_cross, tol_tl, tstep);
                i_prior_cross = i-1;
            } else {
                i_prior_cross = i;
                //move_corr_window(neuron.tin, tail_l, t_cross, tol_tl, tstep);
            }
            if (jl::debug) {
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
            tsp.push_back(t_cross*neuroLib.tstep+neuron.tref/2);
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
            if (jl::debug) {
                std::cout << " backed at " << tBack*neuroLib.tstep << " v " << vBack << std::endl;
                std::cout << " input from " << i_prior_cross << " to " << i << std::endl;
                std::cout << " input " << neuron.tin[i] << std::endl;
            }
            jumpy::getNear(neuroLib.vRange,neuroLib.nv,
                            vBack, r_, i_, j_);
            cross.vCross.push_back(IJR(i_,j_,r_));
            cross.nCross++;
            cross.iCross.push_back(cross.v.size());
            cross.tCross.push_back(tBack);
            cout << "crossed " << endl;

            old_tail_l = tail_l;
            move_corr_window(neuron.tin, tail_l, tBack, tol_tl, tstep);
            jl::update_info_after_cross(input, neuroLib, cross, neuron, tBack, vBack, i_prior_cross, tail_l, i, corrSize, tsp, spiked);
            cout << " updated" << endl;
            tail_l = i + 1;
            crossed = jl::check_crossing(input, neuroLib, cross, neuron, tol_tl, end_t, tail_l, i, jnd, t_cross);
            cout << "------------" << endl;
        }
        cout << " dead with " << i << ", v " << jnd.v.back() << " at " << jnd.t.back() << endl;
        input.print_this(i);
    }
    return nc;
}
#endif
