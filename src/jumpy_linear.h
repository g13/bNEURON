#ifndef nJL_H 
#define nJL_H
#include <cmath>
#include <vector>
#include "typedefs.h"
#include "jumpy_bilinear.h"
namespace jl{
    const bool debug = false;
    const bool debug2 = false;
    inline double find_v_at_t(Input &input, nNL &neuroLib, Cross &cross, size head, size tail_l, double t, double tCross, double tol_tl, double v) {
        long i;
        size j, it;
        double dt;
        dt = t - tCross;
        if (dt < tol_tl) {
            v = add_vinit_contribution(neuroLib.vLeak, cross.vCross.back(), dt);
        }
        for (i=head; i>=tail_l; i--) {
            dt = t - input.dt[i];
            if (dt > tol_tl) break;
            it = static_cast<size>(round(dt));
            add_input_i_contribution(i,it,neuroLib,input,v);
        }
        return v;
    }

    inline double interp_for_t_cross(double v_right, double v_left, double t_right, double t_left, size head, size tail_l, double tCross, double tol_tl, nNL &neuroLib, Cross &cross, Input &input, nNS &neuron, double &v) {
        size ival = 0;
        double v1, v2, v0 = neuron.vRest;
        double t_cross1, t_cross2, t_cross;
        do {
            assert(neuron.vThres >= v_left);
            if (abs(t_right - t_left)<1e-14) {
                v = v_right;
                return t_right;
            }
            if (t_right < t_left) {
                std::cout << t_left << " < " << t_right << std::endl;
                assert(t_right >= t_left);
            }
            t_cross1 = t_left + (neuron.vThres - v_left)/(v_right-v_left)*(t_right-t_left);
            //t_cross1 = ceil(t_cross1);
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
            if (fabs(v1-neuron.vThres)<neuron.vTol) {
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
            //t_cross2 = ceil(t_cross2);
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
            if (fabs(v2-neuron.vThres)<neuron.vTol) {
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

    bool check_crossing(Input &input, nNL &neuroLib, Cross &cross, nNS &neuron, double tol_tl, double end_t, size tail_l, size head, jND &jnd, double &t_cross);

    bool update_vinit_of_new_input_check_crossing(Input &input, Cross &cross, nNL &neuroLib, nNS &neuron, size head, size tail_l, double tol_tl, double end_t, jND &jnd, double &t_cross, size corrSize);

    void update_info_after_cross(Input &input, nNL &neuroLib, Cross &cross, nNS &neuron, double tCross, double vCross, size i_prior, size tail, size head, size corrSize, std::vector<double>&tsp, bool spiked);
}

unsigned int nsyn_jLinear(nNS &neuron, nNL &neuroLib, Input &input, jND &jnd, Cross &cross, double end_t, std::vector<double> &v, size corrSize, std::vector<double> &tsp, double vStop, vector<bool> &ei, int afterSpikeBehavior, bool spikeShape);

#endif
