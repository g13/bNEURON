#ifndef JB_H 
#define JB_H
#include <cmath>
#include <vector>
#include "typedefs.h"
#include "nNeuroLib.h"
#include "nNeuroSt.h"

namespace jb{
    const bool debug = true;
    const bool debug2 = true;
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
}
inline void move_corr_window(vector<double> &input, size &tail, double t_head, double tol_tb, double tstep) {
    double dt;
    dt = t_head - round(input[tail]/tstep);
    while (dt > tol_tb) {
        tail++;
        dt = t_head - round(input[tail]/tstep);
    }
}
inline void reverse_corr_window(vector<double> &input, size &tail, double t_head, double tol_tb, double tstep) {
    double dt;
    do {
        tail--;
        if (tail<0) {
            break;
        }
        dt = t_head - round(input[tail]/tstep);
    } while (dt <= tol_tb);
    tail++;
}
inline void add_relation_between_new_and_ith_input(Input &input, size inew, size i, nNL &neuroLib) {
    double t = input.t[inew] - input.t[i];
    if (t < 0) {
        std::cout << i << ", " << inew << std::endl;
        std::cout << input.t[i] << ", " << input.t[inew] << std::endl;
    }
    double r_;
    size i_,j_;
    //size idt = static_cast<size>(t);
    jb::getNear(neuroLib.idtRange, neuroLib.ndt, t, r_,i_,j_);
    input.bir[i].dTijr.push_back(IJR(i_,j_,r_));
    //input.bir[i].Tijr.push_back(IJR(idt,idt+1,t-idt));
}
inline double linear_interp_tMax(double ***tMax, IJR &v, IJR &dT, size i) {
    double base = tMax[v.i][dT.i][i];
    return round(base + v.r *  (tMax[v.j][dT.i][i] -base)
                      + dT.r * (tMax[v.i][dT.j][i] -base));
}
inline void add_new_input_info(nNL &neuroLib, Input &input, Cross &cross, double t0, double ***tMax, double v, size corrSize, size ID) {
    double t = input.t.back();
    double r_;
    size i_, j_;
    jb::getNear(neuroLib.vRange, neuroLib.nv, v, r_,i_,j_);
    input.Vijr.push_back(IJR(i_,j_,r_));
    input.dTijr.push_back(IJR(0,1,0));
    input.dt.push_back(t);
    input.Tmpijr.push_back(IJR(0,1,0));
    input.cCross.push_back(0);
    input.tMax.push_back(linear_interp_tMax(tMax, input.Vijr.back(), input.dTijr.back(),ID));
    input.bir.push_back(BilinearRelationships(corrSize));
    input.inTref.push_back(0);
}
inline double add_vinit_contribution(nNL &neuroLib, IJR &v, double dt) { 
    size i = static_cast<size>(round(dt));
    double base = neuroLib.vLeak[v.i][i];
    return base + v.r * (neuroLib.vLeak[v.j][i] -base);
}
inline double linear_interp_leak(double ***PSP, IJR &v, size ID, size idt) {
    double base = PSP[v.i][ID][idt];
    return base + v.r *  (PSP[v.j][ID][idt] -base);
}
inline double linear_interp_PSP(double ****PSP, IJR &v, IJR &dT, size ID, size idt, size *idtRange) {
    size iidt = idt + idtRange[dT.i];
    size jjdt = idt + idtRange[dT.j];
    //cout << " dT.i " << dT.i << endl;
    //cout << " dT.r " << dT.r << endl;
    double base = PSP[v.i][dT.i][ID][iidt];
    return base + v.r *  (PSP[v.j][dT.i][ID][iidt] -base)
                + dT.r * (PSP[v.i][dT.j][ID][jjdt] -base);
}
inline void add_input_i_contribution(size i, size idt, nNL &neuroLib, Input &input, double &v) {
    //cout << "synapse " << input.ID[i] << " contributing " << input.t[i] << endl;
    v = v + linear_interp_PSP(neuroLib.sPSP, 
                              input.Vijr[i], input.dTijr[i], 
                              input.ID[i], idt, neuroLib.idtRange);
}
inline double linear_interp_kV(double ******kV, IJR &v, IJR &ijdT, size dT, size it, size *idtRange, size ndt, size i, size j, bool debug) {
    if (dT == 0) {
        size jt = idtRange[ijdT.j] + it;
        it = idtRange[ijdT.i] + it;
        if (debug) {
            cout << "     dT " << dT << ": " << ijdT.i << ", " << ijdT.j << ", r" << ijdT.r << endl;
            cout << "     it " << it <<", jt" << jt << endl;
        }
        double base = kV[v.i][ijdT.i][i][j][ijdT.i][it];
        double dv = kV[v.j][ijdT.i][i][j][ijdT.i][it] - base;
        double dt = kV[v.i][ijdT.j][i][j][ijdT.j][jt] - base;
        if (debug) {
            cout << "     base " << base << ", dv " << dv << ", dt " << dt << endl;
        }
        return base + v.r * dv + ijdT.r * dt;
    } else {
        size idt0, jdt0, idt1, jdt1;
        double rdt0, rdt1;
        jb::getNear(idtRange,ndt,idtRange[ijdT.i] + dT, rdt0, idt0, jdt0);
        jb::getNear(idtRange,ndt,idtRange[ijdT.j] + dT, rdt1, idt1, jdt1);
        size iidt0 = idtRange[idt0]+it;
        size jjdt0 = idtRange[jdt0]+it;
        size iidt1 = idtRange[idt1]+it;
        size jjdt1 = idtRange[jdt1]+it;
        if (debug) {
            cout << "     dT " << dT << ": " << iidt0 <<", " << jjdt0 << ", r" << rdt0 << "; " << iidt1 << ", " << jjdt1 << ", r" << rdt1 << " it " << it << endl;
        }

        double base0 = kV[v.i][ijdT.i][i][j][idt0][iidt0];
        double base = base0 + rdt0*(kV[v.i][ijdT.i][i][j][jdt0][jjdt0]-base0);

        double dv0 = kV[v.j][ijdT.i][i][j][idt0][iidt0];
        double dv = dv0 + rdt0*(kV[v.j][ijdT.i][i][j][jdt0][jjdt0]-dv0);

        double dt0 = kV[v.i][ijdT.j][i][j][idt1][iidt1];
        double dt = dt0 + rdt1*(kV[v.i][ijdT.j][i][j][jdt1][jjdt1]-dt0);

        if (debug) {
            cout << base0 << ", " << dv0 << ", " << dt0 << endl;
            cout << "inds: "<< v.i << ", "<< ijdT.j<< ", " << i << ", " << j << ", " << idt1 << ", " << iidt1 << endl;
        }
        return base + v.r * (dv-base) + ijdT.r * (dt-base);
    }
}
inline void add_input_i_j_bilinear_contribution(Input &input, nNL &neuroLib, size i, size j, size it, double &v, bool debug) {
    size k = j-i-1;
    if (input.cCross[j] > input.cCross[i]) {
        cout << "i.cCross " << input.cCross[i] << " < j.cCross " << input.cCross[j] << endl;
        cout << i << " t " << input.t[i] << j << " t " << input.t[j] <<  "-> " << input.t[j] - input.t[i] << endl;
        assert(input.cCross[j] <= input.cCross[i]);
    }
    //if (jb::debug) {
    //    cout << "        " << " bilinear " << i << ", " << j << " contributing" << endl;
    //    cout << "        " << input.Vijr[j].i << ", " << input.bir[i].dTijr[k].i << ", " <<
    //                      input.dt[j] << " - " <<input.t[j] << ", " << it << endl;
    //}
    v += linear_interp_kV(neuroLib.kV, input.Vijr[j], input.bir[i].dTijr[k], 
                          input.dt[j]-input.t[j], it, neuroLib.idtRange, neuroLib.ndt, input.ID[i], input.ID[j], debug);
    //if (jb::debug) {
    //    cout << "       " << " v -> " << v << endl;
    //}
}
inline double find_v_at_t(Input &input, nNL &neuroLib, Cross &cross, size head, size tail_l, size tail_b, double t, double tCross, double tol_tl, double v, bool debug) {
    long i;
    size j, it;
    double dt, dt0;
    dt = t - tCross;
    if (dt < tol_tl) {
        v = add_vinit_contribution(neuroLib, cross.vCross.back(), dt);
    }
    for (i=head; i>=tail_l; i--) {
        dt0 = t - input.t[i];
        if (dt0 > tol_tl) break;
        dt = t - input.dt[i];
        it = static_cast<size>(round(dt));
        if (!input.inTref[i]) {
            add_input_i_contribution(i,it,neuroLib,input,v);
            //cout << " i " << i << " cCross" << input.cCross[i] << endl;
            if (i < tail_b) continue;
            for (j=head; j>i; j--) {
                //cout << "   j " << j << " cCross" << input.cCross[j] << endl;
                dt = t - input.dt[j];
                it = static_cast<size>(round(dt));
                if (!input.inTref[j]) {
                    add_input_i_j_bilinear_contribution(input, neuroLib, i, j, it, v, debug);
                }
            }
        }
    }
    return v;
}
inline double parabola(double t_left, double v_left, double t_right, double v_right, double t_mid, double v_mid, double vThres) {
    double t1 = t_mid-t_left;
    double t1_2 = t1*t1;
    double t2 = t_right-t_left;
    double t2_2 = t2*t2; 
    double denorm = t1_2*t2-t2_2*t1;
    assert(denorm != 0.0);
    v_right = v_right-v_left;
    v_mid = v_right-v_left;
    vThres = vThres-v_left;
    double a = (v_mid*t2-v_right*t1)/denorm;
    double b = (v_right*t1_2-v_mid*t2_2)/denorm;
    return (-b+sqrt(b*b+4*a*vThres))/(2*a) + t_left;
}
inline void getLR(double &v_left, double v1, double v2, double &v_right, double &t_left, double t_cross1, double t_cross2, double &t_right, double vTar) {
    if (vTar < v1) {
        v_right = v1;
        t_right = t_cross1;
    } else {
        if (vTar < v2) {
            v_left = v1;
            t_left = t_cross1;
            v_right = v2;
            t_right = t_cross2;
        } else {
            v_left = v2;
            t_left = t_cross2;
        }
    }
}

double interp_for_t_cross(double v_right, double v_left, double t_right, double t_left, size head, size tail_l, size tail_b, double tCross, double tol_tl, nNL &neuroLib, Cross &cross, Input &input, nNS &neuron, double &v, bool debug);

bool check_crossing(Input &input, nNL &neuroLib, Cross &cross, nNS &neuron, double tol_tl, double tol_tb, double end_t, size tail_l, size tail_b, size head, jND &jnd, double &t_cross, bool debug);

bool update_vinit_of_new_input_check_crossing(Input &input, Cross &cross, nNL &neuroLib, nNS &neuron, size head, size tail_l, size tail_b, double tol_tl, double tol_tb, double end_t, jND &jnd, double &t_cross, size corrSize, bool debug);

void update_info_after_cross(Input &input, nNL &neuroLib, Cross &cross, nNS &neuron, double tCross, double vCross, size i_prior, size tail, size head, size corrSize, std::vector<double>&tsp, bool spiked);

unsigned int nsyn_jBilinear(nNS &neuron, nNL &neuroLib, Input &input, jND &jnd, Cross &cross, double end_t, std::vector<double> &v, size corrSize, std::vector<double> &tsp, double vStop, vector<bool> &ei);
#endif
