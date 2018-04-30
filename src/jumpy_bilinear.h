#ifndef JB_H 
#define JB_H
#include <cmath>
#include <vector>
#include "typedefs.h"
#include "nNeuroSt.h"
#include "Yale_NEURON_PyAPI.h"

using std::vector;
using std::cout;
using std::endl;
namespace jb {
    const bool debug = false;
    const bool debug2 = false;
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

inline void move_corr_window(vector<double> &input, size &tail, double t_head, double tol_t, double tstep) {
    double dt;
    double tt_head = t_head*tstep;
    double tol_tt = tol_t*tstep;
    dt = tt_head - input[tail];
    while (dt >= tol_tt && tail < input.size()-1) {
        tail++;
        dt = tt_head - input[tail];
    }
}

inline void move_corr_window_i(vector<double> &input, size &tail, double t_head, double tol_t) {
    double dt;
    dt = t_head - input[tail];
    while (dt >= tol_t && tail < input.size()-1) {
        tail++;
        dt = t_head - input[tail];
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
    } while (dt < tol_tb);
    tail++;
}

inline void add_relation_between_new_and_ith_input(Input &input, size inew, size i, size *idtRange, size ndt) {
    double t = input.t[inew] - input.t[i];
    if (t < 0) {
        cout << i << ", " << inew << endl;
        cout << input.t[i] << ", " << input.t[inew] << endl;
    }
    double r_;
    size i_,j_;
    //size idt = static_cast<size>(t);
    jb::getNear(idtRange, ndt, t, r_,i_,j_);
    input.bir[i].dTijr.push_back(IJR(i_,j_,r_));
    //input.bir[i].Tijr.push_back(IJR(idt,idt+1,t-idt));
}

inline double linear_interp_tMax(size ***tMax, IJR &v, IJR &dT, size i) {
    double base = tMax[v.i][dT.i][i];
    return round(base + v.r *  (tMax[v.j][dT.i][i] -base)
                      + dT.r * (tMax[v.i][dT.j][i] -base));
}

inline void add_new_input_info(double *vRange, size nv, Input &input, Cross &cross, double t0, size ***tMax, int **sfireCap, double v, size corrSize, size ID) {
    double t = input.t.back();
    double r_;
    size i_, j_;
    jb::getNear(vRange, nv, v, r_,i_,j_);
    input.Vijr.push_back(IJR(i_,j_,r_));
    input.dTijr.push_back(IJR(0,1,0));
    if (j_ >= sfireCap[0][ID] && r_ > 0.1) {
        // if sufficiently close to spike case
        input.tMax.push_back(tMax[j_][0][ID]);
    } else {
        input.tMax.push_back(linear_interp_tMax(tMax, input.Vijr.back(), input.dTijr.back(),ID));
    }
    input.dt.push_back(t);
    input.idt.push_back(0);
    input.Tmpijr.push_back(IJR(0,1,0));
    input.cCross.push_back(0);
    input.bir.push_back(BilinearRelationships(corrSize));
}

inline double add_vinit_contribution(double **vLeak, IJR &v, double dt) { 
    size i = static_cast<size>(round(dt));
    double base = vLeak[v.i][i];
    return base + v.r * (vLeak[v.j][i] -base);
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

inline double linear_interp_kV(double ******kV, IJR &v, IJR &ijdT, size dT, size it, size *idtRange, size ndt, size i, size j, bool debug, bool dtSquare) {
    if (dT == 0) {
        if (debug) {
            cout << "           dT" << dT << ", dt " << it << endl;
            cout << "           rdt0 " << ijdT.r << endl;
            cout << "           idt0 " << ijdT.i << ": " << idtRange[ijdT.i] << endl;
            cout << "           jdt0 " << ijdT.j << ": " << idtRange[ijdT.j] << endl;
        }
        size iidt = idtRange[ijdT.i] + it;
        size jjdt = idtRange[ijdT.j] + it;
        double base = kV[v.i][ijdT.i][i][j][ijdT.i][iidt];
        double dv = kV[v.j][ijdT.i][i][j][ijdT.i][iidt] - base;
        double dt = kV[v.i][ijdT.j][i][j][ijdT.j][jjdt] - base;
        if (debug) {
            cout << "           base " << base << ", dv " << dv << ", dt " << dt << endl;
        }
        return base + v.r * dv + ijdT.r * dt;
    } else {
        size idt, jdt;
        double rdt;
        size dT0 = idtRange[ijdT.i] + ijdT.r*(idtRange[ijdT.j]-idtRange[ijdT.i]);
        if (debug) {
            cout << "           dT0 " << dT0 << ", dT " << dT << ", dt " << it << endl;
            cout << "           rdt0 " << ijdT.r << endl;
            cout << "           idt0 " << ijdT.i << ": " << idtRange[ijdT.i] << endl;
            cout << "           jdt0 " << ijdT.j << ": " << idtRange[ijdT.j] << endl;
        }
        size dT1 = dT0 + dT;
        jb::getNear(idtRange,ndt,dT1,rdt,idt,jdt);
        size iidt = idtRange[idt] + it;
        size jjdt = idtRange[jdt] + it;
        if (debug) {
            cout << "           dT1 " << dT <<  ", dt " << it << endl;
            cout << "           rdt1 " << rdt << endl;
            cout << "           idt1 " << idt << ": " << idtRange[idt] <<", iidt1 " << iidt << endl;
            cout << "           jdt1 " << jdt << ": " << idtRange[jdt] <<", jjdt1 " << jjdt << endl;
        }
        if (dtSquare) {

            double base0 = kV[v.i][ijdT.i][i][j][idt][iidt];
            double base = base0 + rdt*(kV[v.i][ijdT.i][i][j][jdt][jjdt] - base0);

            double t0 = kV[v.i][ijdT.j][i][j][idt][iidt];
            double t1 = t0 + rdt*(kV[v.i][ijdT.j][i][j][jdt][jjdt]-t0);
            
            double v0 = kV[v.j][ijdT.i][i][j][idt][iidt];
            double v1 = v0 + rdt*(kV[v.j][ijdT.i][i][j][jdt][jjdt] - v0);

            return base + v.r * (v1-base) + ijdT.r*(t1-base);
        } else {
            double base0 = kV[v.i][ijdT.i][i][j][idt][iidt];
            double dv = v.r*(kV[v.j][ijdT.i][i][j][idt][iidt]-base0);

            double dt0 = ijdT.r*(kV[v.i][ijdT.j][i][j][idt][iidt]-base0);

            double dt1 = rdt*(kV[v.i][ijdT.i][i][j][jdt][jjdt]-base0);
            return base0 + dv + dt0 + dt1;
        }
    }
}

inline double linear_interp_kV_old(double ******kV, IJR &v, IJR &ijdT, size dT, size it, size *idtRange, size ndt, size i, size j, bool debug) {
    if (dT == 0) {
        if (debug) {
            cout << "           dT " << dT << ": i " << ijdT.i << ", j " << ijdT.j << ", r " << ijdT.r << endl;
            cout << "           dt " << it <<": it " << idtRange[ijdT.i] << ", jt " << idtRange[ijdT.j] << endl;
        }
        size jt = idtRange[ijdT.j] + it;
        it = idtRange[ijdT.i] + it;
        double base = kV[v.i][ijdT.i][i][j][ijdT.i][it];
        double dv = kV[v.j][ijdT.i][i][j][ijdT.i][it] - base;
        double dt = kV[v.i][ijdT.j][i][j][ijdT.j][jt] - base;
        if (debug) {
            cout << "           base " << base << ", dv " << dv << ", dt " << dt << endl;
        }
        return base + v.r * dv + ijdT.r * dt;
    } else {
        size idt0, jdt0, idt1, jdt1;
        double rdt0, rdt1;
        if (debug) {
            cout << "           dT " << dT << ": i " << ijdT.i << ", j " << ijdT.j << ", r " << ijdT.r << endl;
            cout << "           dt " << it <<": it " << idtRange[ijdT.i] << ", jt " << idtRange[ijdT.j] << endl;
        }
        jb::getNear(idtRange,ndt,idtRange[ijdT.i] + dT, rdt0, idt0, jdt0);
        jb::getNear(idtRange,ndt,idtRange[ijdT.j] + dT, rdt1, idt1, jdt1);
        size iidt0 = idtRange[idt0]+it;
        size jjdt0 = idtRange[jdt0]+it;
        size iidt1 = idtRange[idt1]+it;
        size jjdt1 = idtRange[jdt1]+it;
        if (debug) {
            cout << "           rdt0 " << rdt0 << " rdt1 " << rdt1 << endl;
            cout << "           idt0 " << idt0 << ": idt " << idtRange[idt0] <<", iidt0 " << iidt0 << endl;
            cout << "           jdt0 " << jdt0 << ": jdt " << idtRange[jdt0] <<", jjdt0 " << jjdt0 << endl;
            cout << "           idt1 " << idt1 << ": idt " << idtRange[idt1] <<", iidt1 " << iidt1 << endl;
            cout << "           jdt1 " << jdt1 << ": jdt " << idtRange[jdt1] <<", jjdt1 " << jjdt1 << endl;
        }

        double base0 = kV[v.i][ijdT.i][i][j][idt0][iidt0];
        double base = base0 + rdt0*(kV[v.i][ijdT.i][i][j][jdt0][jjdt0]-base0);

        double dv0 = kV[v.j][ijdT.i][i][j][idt0][iidt0];
        double dv = dv0 + rdt0*(kV[v.j][ijdT.i][i][j][jdt0][jjdt0]-dv0);

        double dt0 = kV[v.i][ijdT.j][i][j][idt1][iidt1];
        double dt = dt0 + rdt1*(kV[v.i][ijdT.j][i][j][jdt1][jjdt1]-dt0);

        if (debug) {
            cout << "           base0 " << base0 << ", dv0 " << dv0 << ", dt0 " << dt0 << endl;
            cout << "           base " << base << ", dv " << dv << ", dt " << dt << endl;
        }
        return base + v.r * (dv-base) + ijdT.r * (dt-base);
    }
}

inline void add_input_i_contribution(size i, size idt, nNL &neuroLib, Input &input, double &v) {
    //cout << "synapse " << input.ID[i] << " contributing " << input.t[i] << endl;
    v = v + linear_interp_PSP(neuroLib.sPSP, 
                              input.Vijr[i], input.dTijr[i], 
                              input.ID[i], idt, neuroLib.idtRange);
}

inline void add_input_i_j_bilinear_contribution(Input &input, nNL &neuroLib, size i, size j, size it, double &v, bool debug = false, bool dtSquare = true) {
    size k = j-i-1;
    if (debug) {
        if (input.cCross[j] > input.cCross[i]) {
            cout << "i.cCross " << input.cCross[i] << " < j.cCross " << input.cCross[j] << endl;
            cout << i << " t " << input.t[i] << j << " t " << input.t[j] <<  "-> " << input.t[j] - input.t[i] << endl;
            assert(input.cCross[j] <= input.cCross[i]);
        }
    }
    v += linear_interp_kV(neuroLib.kV, input.Vijr[j], input.bir[i].dTijr[k], 
                          input.dt[j]-input.t[j], it, neuroLib.idtRange, neuroLib.ndt, input.ID[i], input.ID[j], debug, dtSquare);
}

inline double find_v_at_t(Input &input, nNL &neuroLib, Cross &cross, size head, size tail_l, size tail_b, double t, double tCross, double tol_tl, double tol_tb, double v, bool debug, bool dtSquare) {
    size i, j, idt;
    double dt;
    dt = t - tCross;
    if (dt < tol_tl) {
        v = add_vinit_contribution(neuroLib.vLeak, cross.vCross.back(), dt);
    }
    // else v = neuron.vReset (presetted)

    move_corr_window_i(input.t, tail_b, t, tol_tb);
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
        for (j=tail_b; j<i; j++) {
            add_input_i_j_bilinear_contribution(input, neuroLib, j, i, idt, v, debug, dtSquare);
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
    if (jb::debug) {
        assert(denorm != 0.0);
    }
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

double interp_for_t_cross(double v_right, double v_left, double t_right, double t_left, size head, size tail_l, size tail_b, double tCross, double tol_tl, double tol_tb, nNL &neuroLib, Cross &cross, Input &input, double v0, double vtol, double vC, double &v, bool debug, bool dtSquare);

bool check_crossing(Input &input, nNL &neuroLib, Cross &cross, nNS &neuron, double tol_tl, double tol_tb, double end_t, size tail_l, size tail_b, size head, jND &jnd, double &t_cross, double vC, bool debug, bool dtSquare);

bool update_vinit_of_new_input_check_crossing(Input &input, Cross &cross, nNL &neuroLib, nNS &neuron, size head, size &tail_l, size &tail_b, size old_tail_l, size old_tail_b, double tol_tl, double tol_tb, double end_t, jND &jnd, double &t_cross, double vC, size corrSize, bool debug, bool dtSquare);

void update_info_after_cross(Input &input, nNL &neuroLib, Cross &cross, nNS &neuron, double tCross, double vCross, size i_prior, size &tail_l, size &tail_b, size head, size corrSize, int afterCrossBehavior);

unsigned int nsyn_jBilinear(Cell &cell, vector<vector<double>> &spikeTrain, vector<double> dendVclamp, nNS &neuron, nNL &neuroLib, Input &input, jND &jnd, Cross &cross, double end_t, double ignore_t, size corrSize, vector<double> &tsp, double vC, double vB, int afterCrossBehavior, bool spikeShape, bool dtSquare);

#endif
