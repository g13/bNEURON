#ifndef LB_H
#define LB_H
#include <cmath>
#include <vector>
#include "nNeuroSt.h"
#include "Yale_NEURON_PyAPI.h"
using cout;
using endl;
using vector;
namespace lb{
    const bool debug = false;
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

inline void interpPSP(vector<double> &v, size vs, double ****PSP, double *vRange, size *idtRange, size iSyn, int nv, size ndt, double vTar, double dtTar, size limit) {
    size iv, jv, idt, jdt;
    size k;
    double rv, rdt, base;
    getNear(idtRange,ndt,dtTar,rdt, idt, jdt);
    getNear(vRange, nv, vTar,rv, iv, jv);
    size iidt = idtRange[idt];
    size jjdt = idtRange[jdt];
    
    size tl = limit - jjdt;
    //cout << " jjdt " << jjdt << " tlimit " << limit << endl;
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
            //cout << "capacity: " << v.size() << endl;
            //cout << vs+k << " > " << v.size() << endl;
            //cout << " length:"  << limit << "-" << jjdt << " = " << tl << endl;
            assert(vs+k<v.size());
            //cout << v[vs+k] << endl;
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

inline void interpkV(vector<double> &v, size vs, double ******kV, double *vRange, size *idtRange, int nv, size ndt, double vTar, double dtTar0, double dtTar1, size tl, size iSyn, size jSyn) {
    size iv ,jv;
    size idt_,jdt_;
    size idt,jdt;
    size idt0,jdt0;
    size idt1,jdt1;
    size k;
    double rdt,rdt0,rdt1,rv,rdt_;
    double base0,base,dt0,dt,dv0,dv;
    if (iSyn == jSyn) {
        getNear(idtRange+1,ndt-1,dtTar0,rdt,idt,jdt);
        idt++;
        jdt++;
    } else {
        getNear(idtRange,ndt,dtTar0,rdt,idt,jdt);
    }
    double dtTar = dtTar1 - dtTar0;
    getNear(vRange,nv,vTar,rv,iv,jv);
    
    if (dtTar != 0) {
        getNear(idtRange,ndt,idtRange[idt] + dtTar,rdt0,idt0,jdt0);
        getNear(idtRange,ndt,idtRange[jdt] + dtTar,rdt1,idt1,jdt1);
        size iidt0 = idtRange[idt0];
        size jjdt0 = idtRange[jdt0];
        size iidt1 = idtRange[idt1];
        size jjdt1 = idtRange[jdt1];

        for (k=1;k<tl;k++) {
            base0 = kV[iv][idt][iSyn][jSyn][idt0][iidt0+k];
            base = base0 + rdt0*(kV[iv][idt][iSyn][jSyn][jdt0][jjdt0+k] - base0);

            dv0 = kV[jv][idt][iSyn][jSyn][idt0][iidt0+k];
            dv = dv0 + rdt0*(kV[jv][idt][iSyn][jSyn][jdt0][jjdt0+k]-dv0);
            
            dt0 = kV[iv][jdt][iSyn][jSyn][idt1][iidt1+k];
            dt = dt0 + rdt1*(kV[iv][jdt][iSyn][jSyn][jdt1][jjdt1+k]-dt0);
            v[vs+k] += base + rv*(dv-base) + rdt*(dt-base);
        }
    } else {
        size iidt = idtRange[idt];
        size jjdt = idtRange[jdt];
        for (k=1;k<tl;k++) {
            base = kV[iv][idt][iSyn][jSyn][idt][iidt+k];
            //cout << " base " << base;
            v[vs+k] += base + rv  * (kV[jv][idt][iSyn][jSyn][idt][iidt+k]-base)
                            + rdt * (kV[iv][jdt][iSyn][jSyn][jdt][jjdt+k]-base);
        }
    }
}

inline void interpPSP_nV(vector<double> &v, size vs, double ****PSP, size *idtRange, size ndt, size iSyn, double dtTar, size limit, size iv0, bool kVStyle) {
    size k;
    if (dtTar!=0) {
        if (kVStyle) {
            size idt, jdt;
            double rdt;
            double base;
            getNear(idtRange,ndt,dtTar,rdt, idt, jdt);
            size jjdt = idtRange[jdt];
            limit = limit - jjdt;
            for (k=1;k<limit;k++) {
                base = PSP[iv0][idt][iSyn][idtRange[idt]+k];
                v[vs+k] += base + rdt*(PSP[iv0][jdt][iSyn][jjdt+k]-base);
            }
        } else {
            size idtTar = static_cast<size>(round(dtTar));
            limit = limit - idtTar;
            for (k=1;k<limit;k++) {
                 v[vs+k] += PSP[iv0][0][iSyn][idtTar + k];
            }
        }
    } else {
        for (k=1;k<limit;k++) {
             v[vs+k] += PSP[iv0][0][iSyn][k];
        }
    }
}

inline void interpkV0(vector<double> &v, size vs, double ****kV0, size *idtRange, size ndt, double dtTar0, double dtTar1, size tl, size iSyn, size jSyn) {
    size idt,jdt;
    size k;
    double rdt;
    double base;
    if (iSyn == jSyn) {
        getNear(idtRange+1,ndt-1,dtTar0,rdt,idt,jdt);
        idt++;
        jdt++;
    } else {
        getNear(idtRange,ndt,dtTar0,rdt,idt,jdt);
    }
    size idt0 = static_cast<size>(dtTar1);
    for (k=1;k<tl;k++) {
        base = kV0[idt][iSyn][jSyn][idtRange[idt]+idt0+k];
        v[vs+k] += base + rdt * (kV0[jdt][iSyn][jSyn][idtRange[jdt]+idt0+k]-base);
        if (lb::debug2 && k==1) {
            cout << "first value " << base << endl;
        }
        if (lb::debug2 && k==50) {
            cout << "5ms value " << base << endl;
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

unsigned int bilinear_nSyn(Cell &cell, vector<vector<double>> &spikeTrain, vector<double> dendVclamp, vector<double> &v, nNL &neuroLib, nNS &neuron, double run_t, double ignore_t, vector<double> &tsp, double vCross, double vBack, vector<bool> &ei, int afterCrossBehavior, bool spikeShape);

unsigned int linear_nSyn(Cell &cell, vector<vector<double>> &spikeTrain, vector<double> dendVclamp, vector<double> &v, nNL &neuroLib, nNS &neuron, double run_t, double ignore_t, vector<double> &tsp, double vCross, double vBack, vector<bool> &ei, int afterCrossBehavior, bool spikeShape);

unsigned int bilinear0_nSyn(Cell &cell, vector<vector<double>> &spikeTrain, vector<double> dendVclamp, vector<double> &v, nNL &neuroLib, nNS &neuron, double run_t, double ignore_t, vector<double> &tsp, double vCross, double vBack, vector<bool> &ei, int afterCrossBehavior, bool spikeShape, bool kVStyle);

#endif
