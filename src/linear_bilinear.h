#ifndef LB_H
#define LB_H
#include <cmath>
#include <vector>
#include "nNeuroSt.h"
#include "Yale_NEURON_PyAPI.h"
using std::cout;
using std::endl;
using std::vector;
namespace lb{
    const bool debug = true;
    const bool debug2 = false;
}

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

inline void interpPSP(vector<double> &v, size vs, double ****PSP, double *vRange, size *idtRange, size iSyn, int nv, size ndt, double vTar, double dtTar, size tl) {
    size iv, jv, idt, jdt;
    size k;
    double rv, rdt, base;
    getNear(idtRange,ndt,dtTar,rdt, idt, jdt);
    getNear(vRange, nv, vTar,rv, iv, jv);
    size iidt = idtRange[idt];
    size jjdt = idtRange[jdt];
    
    //cout << " jjdt " << jjdt << " tl " << tl << endl;
    //cout << " iv" << iv << " iSyn " << iSyn << endl;
    //cout << " jv" << jv << endl;
    //cout << " idt" << idt << endl;
    //cout << " jdt" << jdt << endl;
    for (k=1;k<tl;k++) {
        //cout << vs + k << ", "  << v[vs+k], " -> ";
        base = PSP[iv][idt][iSyn][iidt+k];
        v[vs+k] += base + rv*(PSP[jv][idt][iSyn][iidt+k]-base) + rdt*(PSP[iv][jdt][iSyn][jjdt+k]-base);
        if (idt>ndt) {
            assert(idt<ndt);
        }
        if (vs+k>v.size()) {
            assert(vs+k<v.size());
        }
    }
}

inline void interpPSP0(vector<double> &v, size vs, double ****PSP, double *vRange, int nv, double vTar, size tl, int iSyn) {
    size i,j,k;
    double r, base;
    getNear(vRange,nv,vTar,r,i,j);
    
    for (k=1;k<tl;k++) {
        base = PSP[i][0][iSyn][k];
        v[vs+k] += base + r*(PSP[j][0][iSyn][k]-base);
    }
}

inline void clampDendRaw(nNL &neuroLib, nNS &neuron, double tol_tb, double tCross, double v0, vector<double> &dendVclamp, size head, double rd, double tCross_old, double vCross_old, vector<double>& vS) {
    double dt, vbase, dv;
    size i, ID, idt, iDt;
    double rv, rdt;
    size iv, jv;
    size iidt, jjdt;
    vector<double> dendv(neuron.nSyn,0);
    if (true) {
        for (i = head; i>0; i--) {
            if (tCross - neuron.tin[i] > tol_tb) break;
            if (neuron.tin[i] > tCross_old) {
                dt = tCross - neuron.tin[i];
                getNear(neuroLib.vRange,neuroLib.nv,vS[i],rv,iv,jv); 
                iidt = 0; jjdt = 1; rdt = 0;
            } else {
                dt = tCross - tCross_old;
                iDt = static_cast<size>(round((tCross_old - neuron.tin[i])/neuroLib.tstep));
                getNear(neuroLib.vRange,neuroLib.nv,vCross_old,rv,iv,jv); 
                getNear(neuroLib.idtRange,neuroLib.ndt,iDt,rdt,iidt,jjdt);
            }
            idt = static_cast<size>(round(dt/neuroLib.tstep));
            ID = neuron.inID[i];
            if (lb::debug2) {
                cout << "it " << idt << endl;
                cout << "v " << iv << ", " << jv << ", " << rv << endl;
                cout << "idt " << iidt << ", " << jjdt << ", " << rdt << endl;
            }
            vbase = neuroLib.dendv[iv][iidt][ID][neuroLib.idtRange[iidt]+idt];
            dv = vbase + rv * (neuroLib.dendv[jv][iidt][ID][neuroLib.idtRange[iidt]+idt] - vbase)
                       + rdt * (neuroLib.dendv[iv][jjdt][ID][neuroLib.idtRange[jjdt]+idt] - vbase);
            if (lb::debug2) {
                cout << i << "th input " << ID << " dv = " << dv << endl;
            }
            dendv[ID] += dv;
        }
        for (int i=0; i<neuron.nSyn; i++) {
            if (abs(dendv[i] - 0) > pow(2,-52)) {
                dendVclamp[i] = v0 + dendv[i] * rd;
                if (lb::debug) {
                    cout << "   dend " << i << " will be clamped at " << dendVclamp[i] << endl;
                }
            }
        }
    } else {
        for (int i=0; i<neuron.nSyn; i++) {
            dendVclamp[i] = neuron.vRest + (v0-neuron.vRest) * rd;
            if (lb::debug) {
                cout << "   dend " << i << " will be clamped at " << dendVclamp[i] << endl;
            }
        }
    }
}

inline void interpkV(vector<double> &v, size vs, double ******kV, double *vRange, size *idtRange, int nv, size ndt, double vTar, double dtTar0, double dtTar1, size tl, size iSyn, size jSyn, bool dtSquare) {
    size idt0,jdt0;
    size idt1,jdt1;
    size iv,jv;
    double rdt0,rdt1,rv;
    getNear(vRange,nv,vTar,rv,iv,jv);
    getNear(idtRange,ndt,dtTar0,rdt0,idt0,jdt0);
    
    if (lb::debug2) {
        cout << " idt0 " << idt0 << endl;
        cout << " jdt0 " << jdt0 << endl;
        assert(jdt0 < ndt);
    }
    if (dtTar1 == dtTar0) {
        size iidt0 = idtRange[idt0];
        size jjdt0 = idtRange[jdt0];
        if (lb::debug2) {
            cout << " iidt0 " << iidt0 << endl;
            cout << " jjdt0 " << jjdt0 << endl;
            cout << " jjdt0 + tl " << jjdt0 + tl << endl;
        }
        for (size i=1;i<tl;i++) {
            double base = kV[iv][idt0][iSyn][jSyn][idt0][iidt0+i];
            v[vs+i] += base + rv  * (kV[jv][idt0][iSyn][jSyn][idt0][iidt0+i]-base)
                            + rdt0 * (kV[iv][jdt0][iSyn][jSyn][jdt0][jjdt0+i]-base);
        }
    } else {
        getNear(idtRange,ndt,dtTar1,rdt1,idt1,jdt1);
        size iidt1 = idtRange[idt1];
        size jjdt1 = idtRange[jdt1];
        if (dtSquare) {
            for (size i=1;i<tl;i++) {
                double base0 = kV[iv][idt0][iSyn][jSyn][idt1][iidt1+i];
                double base = base0 + rdt1*(kV[iv][idt0][iSyn][jSyn][jdt1][jjdt1+i] - base0);

                double t0 = kV[iv][jdt0][iSyn][jSyn][idt1][iidt1+i];
                double t1 = t0 + rdt1*(kV[iv][jdt0][iSyn][jSyn][jdt1][jjdt1+i]-t0);
                
                double v0 = kV[jv][idt0][iSyn][jSyn][idt1][iidt1+i];
                double v1 = v0 + rdt1*(kV[jv][idt0][iSyn][jSyn][jdt1][jjdt1+i] - v0);
                v[vs+i] += base + rv*(v1-base) + rdt1*(t1-base);
            }
        } else {
            for (size i=0;i<tl;i++) {
                double base0 = kV[iv][idt0][iSyn][jSyn][idt1][iidt1+i];
                double dv = rv*(kV[jv][idt0][iSyn][jSyn][idt1][iidt1+i]-base0);

                double dt0 = rdt0*(kV[iv][jdt0][iSyn][jSyn][idt1][iidt1+i]-base0);

                double dt1 = rdt1*(kV[iv][idt0][iSyn][jSyn][jdt1][jjdt1+i]-base0);

                v[vs+i] += base0 + dv + dt0 + dt1;
            }
        }
    }
}

inline void interpPSP_nV(vector<double> &v, size vs, double ****PSP, size *idtRange, size ndt, size iSyn, double dtTar, size tl, size iv0) {
    size k;
    if (dtTar!=0) {
        size idt, jdt;
        double rdt;
        double base;
        getNear(idtRange,ndt,dtTar,rdt, idt, jdt);
        size jjdt = idtRange[jdt];
        for (k=1;k<tl;k++) {
            base = PSP[iv0][idt][iSyn][idtRange[idt]+k];
            v[vs+k] += base + rdt*(PSP[iv0][jdt][iSyn][jjdt+k]-base);
        }
    } else {
        for (k=1;k<tl;k++) {
             v[vs+k] += PSP[iv0][0][iSyn][k];
        }
    }
}

inline void interpkV0(vector<double> &v, size vs, double ****kV0, size *idtRange, size ndt, double dtTar, size tl, size iSyn, size jSyn) {
    size idt,jdt;
    size k;
    double rdt;
    double base;
    getNear(idtRange,ndt,dtTar,rdt,idt,jdt);
    for (k=1;k<tl;k++) {
        base = kV0[idt][iSyn][jSyn][idtRange[idt]+k];
        v[vs+k] += base + rdt * (kV0[jdt][iSyn][jSyn][idtRange[jdt]+k]-base);
        if (lb::debug2 && k==1) {
            cout << "first value " << base << endl;
        }
        if (lb::debug2 && k==tl-1) {
            cout << "ending value at " << tl+idtRange[jdt]-1 << ": " << base << endl;
        }
    }
}

inline void interpVinit(vector<double> &v, size vs,  double **vLeak, double *vRange, size nv, double vTar, size tl) {
    size i,j,k;
    double r;
    getNear(vRange,nv,vTar,r,i,j);
    for (k=0;k<tl;k++)
        v[vs+k] = vLeak[i][k] + r*(vLeak[j][k]-vLeak[i][k]);
}

unsigned int bilinear_nSyn(Cell &cell, vector<vector<double>> &spikeTrain, vector<double> dendVclamp, double rd, vector<double> &v, nNL &neuroLib, nNS &neuron, double run_t, double ignore_t, vector<double> &tsp, double vCross, double vBack, int afterCrossBehavior, bool spikeShape, bool dtSquare);

unsigned int linear_nSyn(Cell &cell, vector<vector<double>> &spikeTrain, vector<double> dendVclamp, double rd, vector<double> &v, nNL &neuroLib, nNS &neuron, double run_t, double ignore_t, vector<double> &tsp, double vCross, double vBack, int afterCrossBehavior, bool spikeShape);

unsigned int bilinear0_nSyn(Cell &cell, vector<vector<double>> &spikeTrain, vector<double> dendVclamp, double rd, vector<double> &v, nNL &neuroLib, nNS &neuron, double run_t, double ignore_t, vector<double> &tsp, double vCross, double vBack, int afterCrossBehavior, bool spikeShape, bool kVStyle, bool dtSquare);

#endif
