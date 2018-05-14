#ifndef nJL_H 
#define nJL_H
#include <cmath>
#include <vector>
#include "typedefs.h"
#include "jumpy_bilinear.h"
#include "Yale_NEURON_PyAPI.h"
using std::cout;
using std::endl;
using std::vector;

namespace jl {
    const bool debug = true;
    const bool debug2 = true;

    inline double find_v_at_t(Input &input, nNL &neuroLib, Cross &cross, size head, size tail_l, double t, double tCross, double tol_tl, double tol_tb, double v) {
        size i, j, idt;
        double dt;
        dt = t - tCross;
        if (dt < tol_tl) {
            v = add_vinit_contribution(neuroLib.vLeak, cross.vCross.back(), dt);
        }
        // else v = neuron.vReset (presetted)
    
        if (dt < tol_tb) {
            move_corr_window_i(input.t, tail_l, t, tol_tb);
        } else {
            if (dt < tol_tl) {
                move_corr_window_i(input.t, tail_l, t, dt);
            } else {
                move_corr_window_i(input.t, tail_l, t, tol_tl);
            }
        }
    
        for (i=tail_l; i<=head; i++) {
            dt = t - input.dt[i];
            idt = static_cast<size>(round(dt));
            add_input_i_contribution(i,idt,neuroLib,input,v);
        }
        return v;
    }

    inline double interp_for_t_cross(double v_right, double v_left, double t_right, double t_left, size head, size tail_l, double tCross, double tol_tl, double tol_tb, nNL &neuroLib, Cross &cross, Input &input, double v0, double vtol, double vC, double &v) {
        size ival = 0;
        double v1, v2;
        double t_cross1, t_cross2, t_cross;
        do {
            if (debug) {
                assert(vC >= v_left);
            }
            assert(t_right >= t_left);
            if (t_right - t_left <= 1) {
                v = v_right;
                return t_right;
            }
            double dt = (vC - v_left)/(v_right-v_left)*(t_right-t_left);
            if (dt < 1.0) {
                dt = 1.0;
            }
            t_cross1 = t_left + dt;
            //t_cross1 = ceil(t_cross1);
            if (debug2) {
                cout << " find v at " << t_cross1 << endl;
                cout <<"t: " << t_left << ", " << t_right << endl;
                assert(t_cross1 >= t_left);
                assert(t_right > t_left);
            }
            v1 = find_v_at_t(input, neuroLib, cross, head, tail_l, t_cross1, tCross, tol_tl, tol_tb, v0);
            if (debug2) {
                cout << v_left << ", " << v1 << ", " << v_right << endl;
                ival++;
            }
            if (v1-vC<vtol && v1>vC) {
                v = v1;
                t_cross = t_cross1;
                break;
            }
    
            // solve for a, b of f(t)-v_left = a(t-t_left)^2 + b(t-t_left);
            // solve for t when a(t-t_left)^2 + b(t-t_left) = (vThres - v_left)
            t_cross2 = parabola(t_left,v_left,t_right,v_right,t_cross1, v1, vC);
            //t_cross2 = ceil(t_cross2);
            if (debug2) {
                cout << " v: " << v_left << ", " << v1 << ", " << v_right << endl;
                cout << " t: " << t_left << ", " << t_cross1 << ", " << t_right << endl;
                cout << " t_cross2 " << t_cross2 << endl;
                assert(t_cross2 >= t_left);
                assert(t_cross2 <= t_right);
            }
            // if t interval smaller than sample temp resolution, apply the quadratic interpolation as solution
            if (debug2) {
                cout << " find v at " << t_cross2 << endl;
            }
            v2 = find_v_at_t(input, neuroLib, cross, head, tail_l, t_cross2, tCross, tol_tl, tol_tb, v0);
            if (debug2) {
                ival++;
            }
            if (v2-vC<vtol && v2>vC) {
                v = v2;
                t_cross = t_cross2;
                break;
            }
            if (debug2) {
                cout << " v: " << v_left << ", " << v1 << ", " << v2 << ", " << v_right << endl;
                cout << " t: " << t_left << ", " << t_cross1 << ", " << t_cross2 << ", " << t_right << endl;
            }
            if (v2>v1) {
                getLR(v_left,v1,v2,v_right,t_left,t_cross1,t_cross2,t_right, vC);
            } else {
                getLR(v_left,v2,v1,v_right,t_left,t_cross2,t_cross1,t_right, vC);
            } 
            if (debug2) {
                cout << " v: " << v_left << ", " << v_right << endl;
                cout << " t: " << t_left << ", " << t_right << endl;
            }
        } while (true);
        if (debug2) {
            cout << ival << " evaluations" << endl;
        }
        return t_cross;
    }

    bool check_crossing(Input &input, nNL &neuroLib, Cross &cross, nNS &neuron, double tol_tl, double tol_tb, double end_t, size tail_l, size head, jND &jnd, double &t_cross, double vC);

    bool update_vinit_of_new_input_check_crossing(Input &input, Cross &cross, nNL &neuroLib, nNS &neuron, size head, size &tail_l, size old_tail_l, double tol_tl, double tol_tb, double end_t, jND &jnd, double &t_cross, double vC, size corrSize);

    void update_info_after_cross(Input &input, nNL &neuroLib, Cross &cross, nNS &neuron, double tCross, double vCross, size i_prior, size &tail, size head, size corrSize, int afterCrossBehavior);
}

unsigned int nsyn_jLinear(Cell &cell, vector<vector<double>> &spikeTrain, vector<double> dendVclamp, nNS &neuron, nNL &neuroLib, Input &input, jND &jnd, Cross &cross, double end_t, double ignore_t, size corrSize, vector<double> &tsp, double vC, double vB, int afterCrossBehavior, bool spikeShape);

#endif
