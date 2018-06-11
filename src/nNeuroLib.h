#ifndef NNL_H
#define NNL_H
#include <cmath>
#include <cstring>
#include <cassert>
#include "matFunc.h"
#include "typedefs.h"

using std::cout;
using std::endl;
using std::memcpy;

struct nNeuroLib {
    size ndt,nt,nv,nSyn,nvAS,nvASt,nvNS,nvNSt;

    double ****sPSP,    *sPSP_ptr, 
           ******kV,     *kV_ptr,
           ****kV0,     *kV0_ptr,
           *vRange, *dtRange, tstep,
           *vASrange,   **vAS, *vAS_ptr,
           *vNSrange,   **vNS, *vNS_ptr,
           **vLeak, *vLeak_ptr;
    size ***tMax, *tMax_ptr;
    double ***dendVleak, *dendVleak_ptr;
    int *loc;
    int **fireCap,  *fireCap_ptr;
    int **sfireCap,  *sfireCap_ptr;
    int ****bfireCap,  *bfireCap_ptr;
    double ****dendv, *dendv_ptr;
    double *pos;
    double *gList;
    double *dist;
    bool *ei;
    int nE, nI;
    //double vReset, vThres;
    double tol_tl;
    double tol_tb;
    size *idtRange;
    const char *file;

    nNeuroLib(const char *filename);
    void clearLib(); // clear mem
};

typedef struct nNeuroLib nNL;
#endif
