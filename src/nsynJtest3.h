#ifndef NJ_H 
#define NJ_H
#include <cmath>
#include <vector>
#include "typedefs.h"
#include "nNeuroLib.h"
#include "nNeuroSt2.h"
#include "neuroTest.h"

namespace jumpy{
    const bool debug = true;
    const bool debug2 = true;
    bool debug3;
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
    jumpy::getNear(neuroLib.idtRange, neuroLib.ndt, t, r_,i_,j_);
    input.bir[i].dTijr.push_back(IJR(i_,j_,r_));
    //input.bir[i].Tijr.push_back(IJR(idt,idt+1,t-idt));
}
inline double linear_interp_tMax(double ***tMax, struct IJR &v, struct IJR &dT, size i) {
    double base = tMax[v.i][dT.i][i];
    return round(base + v.r *  (tMax[v.j][dT.i][i] -base)
                      + dT.r * (tMax[v.i][dT.j][i] -base));
}
inline void add_new_input_info(nNL &neuroLib, Input &input, Cross &cross, double t0, double ***tMax, double v, size corrSize, size ID) {
    double t = input.t.back();
    double r_;
    size i_, j_;
    jumpy::getNear(neuroLib.vRange, neuroLib.nv, v, r_,i_,j_);
    input.Vijr.push_back(IJR(i_,j_,r_));
    input.dTijr.push_back(IJR(0,1,0));
    input.dt.push_back(t);
    input.Tmpijr.push_back(IJR(0,1,0));
    input.cCross.push_back(0);
    input.tMax.push_back(linear_interp_tMax(tMax, input.Vijr.back(), input.dTijr.back(),ID));
    input.bir.push_back(BilinearRelationships(corrSize));
    input.inTref.push_back(0);
}

inline double add_vinit_contribution(nNL &neuroLib, struct IJR &v, double dt) { 
    size i = static_cast<size>(round(dt));
    double base = neuroLib.vLeak[v.i][i];
    return base + v.r * (neuroLib.vLeak[v.j][i] -base);
}
inline double linear_interp_leak(double ***PSP, struct IJR &v, size ID, size idt) {
    double base = PSP[v.i][ID][idt];
    return base + v.r *  (PSP[v.j][ID][idt] -base);
}
inline double linear_interp_PSP(double ****PSP, struct IJR &v, struct IJR &dT, size ID, size idt, size *idtRange) {
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

inline double linear_interp_kV(double ******kV, struct IJR &v, struct IJR &ijdT, size dT, size it, size *idtRange, size ndt, size i, size j) {
    if (dT == 0) {
        size jt = idtRange[ijdT.j] + it;
        it = idtRange[ijdT.i] + it;
        if (jumpy::debug3) {
            cout << "     dT " << dT << ": " << ijdT.i << ", " << ijdT.j << ", r" << ijdT.r << endl;
            cout << "     it " << it <<", jt" << jt << endl;
        }
        double base = kV[v.i][ijdT.i][i][j][ijdT.i][it];
        double dv = kV[v.j][ijdT.i][i][j][ijdT.i][it] - base;
        double dt = kV[v.i][ijdT.j][i][j][ijdT.j][jt] - base;
        if (jumpy::debug3) {
            cout << "     base " << base << ", dv " << dv << ", dt " << dt << endl;
        }
        return base + v.r * dv + ijdT.r * dt;
    } else {
        size idt0, jdt0, idt1, jdt1;
        double rdt0, rdt1;
        jumpy::getNear(idtRange,ndt,idtRange[ijdT.i] + dT, rdt0, idt0, jdt0);
        jumpy::getNear(idtRange,ndt,idtRange[ijdT.j] + dT, rdt1, idt1, jdt1);
        size iidt0 = idtRange[idt0]+it;
        size jjdt0 = idtRange[jdt0]+it;
        size iidt1 = idtRange[idt1]+it;
        size jjdt1 = idtRange[jdt1]+it;
        if (jumpy::debug3) {
            cout << "     dT " << dT << ": " << iidt0 <<", " << jjdt0 << ", r" << rdt0 << "; " << iidt1 << ", " << jjdt1 << ", r" << rdt1 << " it " << it << endl;
        }

        double base0 = kV[v.i][ijdT.i][i][j][idt0][iidt0];
        double base = base0 + rdt0*(kV[v.i][ijdT.i][i][j][jdt0][jjdt0]-base0);

        double dv0 = kV[v.j][ijdT.i][i][j][idt0][iidt0];
        double dv = dv0 + rdt0*(kV[v.j][ijdT.i][i][j][jdt0][jjdt0]-dv0);

        double dt0 = kV[v.i][ijdT.j][i][j][idt1][iidt1];
        double dt = dt0 + rdt1*(kV[v.i][ijdT.j][i][j][jdt1][jjdt1]-dt0);

        if (jumpy::debug3) {
            cout << base0 << ", " << dv0 << ", " << dt0 << endl;
            cout << "inds: "<< v.i << ", "<< ijdT.j<< ", " << i << ", " << j << ", " << idt1 << ", " << iidt1 << endl;
        }
        return base + v.r * (dv-base) + ijdT.r * (dt-base);
    }
}

inline void add_input_i_j_bilinear_contribution(Input &input, nNL &neuroLib, size i, size j, size it, double &v) {
    size k = j-i-1;
    if (input.cCross[j] > input.cCross[i]) {
        cout << "i.cCross " << input.cCross[i] << " < j.cCross " << input.cCross[j] << endl;
        cout << i << " t " << input.t[i] << j << " t " << input.t[j] <<  "-> " << input.t[j] - input.t[i] << endl;
        assert(input.cCross[j] <= input.cCross[i]);
    }
    //if (jumpy::debug) {
    //    cout << "        " << " bilinear " << i << ", " << j << " contributing" << endl;
    //    cout << "        " << input.Vijr[j].i << ", " << input.bir[i].dTijr[k].i << ", " <<
    //                      input.dt[j] << " - " <<input.t[j] << ", " << it << endl;
    //}
    v += linear_interp_kV(neuroLib.kV, input.Vijr[j], input.bir[i].dTijr[k], 
                          input.dt[j]-input.t[j], it, neuroLib.idtRange, neuroLib.ndt, input.ID[i], input.ID[j]);
    //if (jumpy::debug) {
    //    cout << "       " << " v -> " << v << endl;
    //}
}
inline double find_v_at_t(Input &input, nNL &neuroLib, Cross &cross, size head, size tail_l, size tail_b, double t, double tCross, double tol_tl, double v) {
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
                    add_input_i_j_bilinear_contribution(input, neuroLib, i, j, it, v);
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
inline double interp_for_t_cross(double v_right, double v_left, double t_right, double t_left, size head, size tail_l, size tail_b, double tCross, double tol_tl, nNL &neuroLib, Cross &cross, Input &input, nNS &neuron, double &v) {
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
        t_cross1 = ceil(t_cross1);
        assert(t_cross1 >= t_left);
        if (jumpy::debug2) {
            std::cout << " find v at " << t_cross1 << std::endl;
            std::cout <<"t: " << t_left << ", " << t_right << std::endl;
        }
        v1 = find_v_at_t(input, neuroLib, cross, head, tail_l, tail_b, t_cross1, tCross, tol_tl, v0);
        if (jumpy::debug2) {
            std::cout << v_left <<  ", " << v1 << ", " << v_right << std::endl;
        }
        ival++;
        if (fabs(v1-neuron.vThres)<neuron.v_tol) {
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
            if (jumpy::debug2) {
                std::cout << "temp resol reached " << std::endl;
            }
            break;
        }
        // solve for a, b of f(t)-v_left = a(t-t_left)^2 + b(t-t_left);
        // solve for t when a(t-t_left)^2 + b(t-t_left) = (vThres - v_left)
        t_cross2 = parabola(t_left,v_left,t_right,v_right,t_cross1, v1, neuron.vThres);
        t_cross2 = ceil(t_cross2);
        if (jumpy::debug2) {
            std::cout << " v: " << v_left << ", " << v1 << ", " << v_right << std::endl;
            std::cout << " t: " << t_left << ", " << t_cross1 << ", " << t_right << std::endl;
            std::cout << "t_cross2 " << t_cross2 << std::endl;
            assert(t_cross2 >= t_left);
            assert(t_cross2 <= t_right);
        }
        // if t interval smaller than sample temp resolution, apply the quadratic interpolation as solution
        if (jumpy::debug2) {
            std::cout << " find v at " << t_cross2 << std::endl;
        }
        v2 = find_v_at_t(input, neuroLib, cross, head, tail_l, tail_b, t_cross2, tCross, tol_tl, v0);
        //if (jumpy::debug2) {
            ival++;
        //}
        if (fabs(v2-neuron.vThres)<neuron.v_tol) {
            v = v2;
            t_cross = t_cross2;
            break;
        }
        if (jumpy::debug2) {
            std::cout << " v: " << v_left << ", " << v1 << ", " << v2 << ", " << v_right << std::endl;
            std::cout << " t: " << t_left << ", " << t_cross1 << ", " << t_cross2 << ", " << t_right << std::endl;
        }
        if (v2>v1) {
            getLR(v_left,v1,v2,v_right,t_left,t_cross1,t_cross2,t_right, neuron.vThres);
        } else {
            getLR(v_left,v2,v1,v_right,t_left,t_cross2,t_cross1,t_right, neuron.vThres);
        } 
        if (jumpy::debug2) {
            std::cout << " v: " << v_left << ", " << v_right << std::endl;
            std::cout << " t: " << t_left << ", " << t_right << std::endl;
        }
    } while (true);
    //if (jumpy::debug2) {
        std::cout << ival << " evaluations" << std::endl;
    //}
    return t_cross;
}
bool check_crossing(Input &input, nNL &neuroLib, Cross &cross, nNS &neuron, double tol_tl, double tol_tb, double end_t, size tail_l, size tail_b, size head, jND &jnd, double &t_cross) {
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
        if (jumpy::debug2) {
            std::cout << "linear vmax " << vmax <<  " > vThres" << std::endl;
        }
        // perform linear interp iteration until tolerance reaches
        if (jumpy::debug) {
            std::cout << "l+b vmax > vThres, find t and v for cross" << std::endl;
        }
        if (fabs(vmax-neuron.vThres)>neuron.v_tol && tmax-jnd.t.back()>1) {
            t_cross = interp_for_t_cross(vmax, jnd.v.back(), tmax, jnd.t.back(), head, tail_l, tail_b, tCross, tol_tl, neuroLib, cross, input, neuron, v_pass);
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
inline bool update_vinit_of_new_input_check_crossing(Input &input, Cross &cross, nNL &neuroLib, nNS &neuron, size head, size tail_l, size tail_b, double tol_tl, double tol_tb, double end_t, jND &jnd, double &t_cross, size corrSize) {
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
    //if (jumpy::debug) {
    //    cout << " adding " << head << " input t " << input.t[head] << " synapse ID " << neuron.inID[head] << endl;
    //}
    input.ID.push_back(neuron.inID[head]);
    dt = input.t[head] - tCross;
    if (jumpy::debug) {
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
    if (jumpy::debug3) {
        cout << " last v " << jnd.v.back() << endl;
        cout << " start with leak" << v << endl;
    }
    for (i=head-1; i>=tail_l; i--) {
        if (jumpy::debug3) {
            dvtmp = 0;
        }
        dt = input.t[head] - input.dt[i];
        input.bir[i].idt.push_back(static_cast<size>(round(dt)));
        input.bir[i].ID.push_back(neuron.inID[head]);
        if (!input.inTref[i]) {
            add_input_i_contribution(i,input.bir[i].idt.back(),neuroLib,input,v);
            if (jumpy::debug3) {
                cout << " " << neuron.ei[input.ID[i]] << "-" << i << " input contributing " << v-vtmp << ", synapse ID " << input.ID[i] << endl;
                vtmp = v;
            }
            if (i<tail_b) continue;
            for (j=head-1; j>i; j--) {
                ir = head-j-1;
                if (!input.inTref[j]) {
                    add_input_i_j_bilinear_contribution(input, neuroLib, i, j, input.bir[j].idt[ir], v);
                    if (jumpy::debug3) {
                        dvtmp += v-vtmp;
                        cout << "    +" << neuron.ei[input.ID[j]] << "-" << j << " bilinear contributing " << v-vtmp << ", synapse ID " << input.ID[j] << " with dti " << input.bir[i].dTijr[j-i-1].i << ", r" << input.bir[i].dTijr[j-i-1].r << " at vi " << input.Vijr[j].i << ", r" << input.Vijr[j].r << endl;
                        vtmp = v;
                    }
                }
            }
            if (jumpy::debug3) {
                cout << " bilinear total: " << dvtmp << endl;
                cout << " v: " << v << endl;
            }
        }
    }
    // update head's v and tmax
    jnd.t.push_back(input.t[head]);
    jnd.v.push_back(v);
    //if (jumpy::debug) {
    //    cout << "v_" << head << " = " << jnd.v.back() << endl;
    //}
    // check if vmax at left bound very unlikely
    if (v > neuron.vThres) {
        if(jumpy::debug2) {
            std::cout << " crossing upon input not very probable at low input rate " << std::endl;
        }
        jnd.t.pop_back();
        jnd.v.pop_back();
        // perform linear interp iteration until tolerance reaches
        cout << "t: " << jnd.t.back()  << " < " << input.t[head] << endl;
        t_cross = interp_for_t_cross(v, jnd.v.back(), input.t[head], jnd.t.back(), head-1, tail_l, tail_b, tCross, tol_tl, neuroLib, cross, input, neuron, v_pass);
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
    return check_crossing(input, neuroLib, cross, neuron, tol_tl, tol_tb, end_t, tail_l, tail_b, head, jnd, t_cross);
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
    if (input.t.size() != i_start) {
        cout << "before update " << input.t.size() << " != " << i_start << endl;
        assert(input.t.size() == i_start);
    }
    if (i_cross < tail) {
        //if (jumpy::debug) {
            std::cout << " crossing part lasting longer than PSP sample length, very unlikely" << std::endl;
        //}
        cout << " dead old input " << tail-i_cross << endl;
        for (i=i_cross; i<tail; i++) {
            input.junk_this();
        }
        i_start = tail;
    } else {
        // update for input that come before cross and lingers after cross;
        if (jumpy::debug) {
          cout << "cross update starting tail " << tail << endl;
          cout << " linger old inputs: " << i_cross - tail << endl;
        }
        for (i=tail; i<i_cross; i++) {
            if (jumpy::debug) {
                cout << "lingering " << i << " < " << neuron.tin.size() << endl;
            }
            dt = tCross - input.t[i];
            assert(dt>0);
            jumpy::getNear(neuroLib.idtRange,neuroLib.ndt,
                          dt, input.dTijr[i].r, input.dTijr[i].i, input.dTijr[i].j);
            input.dt[i] = tCross;
            input.cCross[i] = cross.nCross;
            input.Vijr[i] = cross.vCross.back();
            if (input.Vijr[i].j < neuroLib.fireCap[input.ID[i]][input.dTijr[i].i]) {
                v = neuroLib.vRange[input.Vijr[i].i] + (neuroLib.vRange[input.Vijr[i].j] - neuroLib.vRange[input.Vijr[i].i]) * input.Vijr[i].r;
                jumpy::getNear(neuroLib.vRange, neuroLib.fireCap[input.ID[i]][input.dTijr[i].i], v, input.Vijr[i].r, input.Vijr[i].i, input.Vijr[i].j);
            }
            input.tMax[i] = linear_interp_tMax(neuroLib.tMax, cross.vCross.back(), input.dTijr[i],i);
            //cout << "   loop ended" << endl;
        }
    }
    // update for new input that lingers after cross
    if (jumpy::debug) {
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
        if (jumpy::debug) {
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
            add_relation_between_new_and_ith_input(input, i, j, neuroLib);
        }
        cout << "   bir pushed" << endl;
    }
    input.assert_size();
}

unsigned int nsyn_jBilinear(nNS &neuron, nNL &neuroLib, Input &input, jND &jnd, Cross &cross, double end_t, std::vector<double> &v, size corrSize, std::vector<double> &tsp, double vStop, Cell &cell, vector<bool> &ei) {
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
        if (i==2) {
            jumpy::debug3 = false;
            //jumpy::debug3 = false;
            //cout << " 39 ith input, just after jump" << endl;
            //cout << endl;
        } else  {
            jumpy::debug3 = false;
            //cout << i << " th input" << endl;
        }
        input.assert_size();
        input.t.push_back(neuron.tin[i]/tstep); 
        assert(input.t.size() == i+1);
        old_tail_l = tail_l;
        old_tail_b = tail_b;
        move_corr_window(neuron.tin, tail_l, input.t[i], tol_tl, tstep);
        move_corr_window(neuron.tin, tail_b, input.t[i], tol_tb, tstep);
        crossed = update_vinit_of_new_input_check_crossing(input, cross, neuroLib, neuron, i, tail_l, tail_b, tol_tl, tol_tb, end_t, jnd, t_cross, corrSize);
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
            if (jumpy::debug) {
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
            if (jumpy::debug) {
                std::cout << " backed at " << tBack*neuroLib.tstep << std::endl;
                std::cout << " input from " << i_prior_cross << " to " << i << std::endl;
                std::cout << " input at " << neuron.tin[i] << std::endl;
            }
            jumpy::getNear(neuroLib.vRange,neuroLib.nv,
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
            crossed = check_crossing(input, neuroLib, cross, neuron, tol_tl, tol_tb, end_t, tail_l, tail_b, i, jnd, t_cross);
            cout << "------------" << endl;
        }
        cout << " dead with " << i << ", v " << jnd.v.back() << " at " << jnd.t.back() << endl;
        input.print_this(i);
    }
    return nc;
}
#endif
