#ifndef NNL_H
#define NNL_H
#include "matFunc.h"
#include <cmath>
using std::cout;
using std::endl;

typedef struct nNeuroLib
{
    size ndt,nt,nv,nSyn;

    double ****sPSP,    *sPSP_ptr, 
           ******kV,     *kV_ptr,
           ****kV0,     *kV0_ptr,
           *vRange, *dtRange, tstep,
           **vLeak, *vLeak_ptr;
    double ***tMax, *tMax_ptr;
    double ***dendVleak, *dendVleak_ptr;
    int *loc;
    int **fireCap,  *fireCap_ptr;
    int **sfireCap,  *sfireCap_ptr;
    int ****bfireCap,  *bfireCap_ptr;
    double ****dendv, *dendv_ptr;
    double *pos;
    double *gList;
    bool *ei;
    //double vReset, vThres;
    double tol_tl;
    double tol_tb;
    size *idtRange;
    const char *file;

    void readLib(const char *filename) {
        MATFile *pmat;
        const char **var;
        mxArray *tmp;
        size arraySize, dimSize[5];
        file = filename;
        int i,n;
        
        openMat(pmat, file);
        
        var = (const char **)matGetDir(pmat, &n);
        if (var == NULL) {
            std::cout << "Error reading content of file: " << file;
            mxFree(var);
            abort();
        } 
        std::cout << "Reading file " << file << "... " << std::endl;
        std::cout << n << " variable(s):" << std::endl;
        for (i=0; i<n; i++) {
            std::cout << var[i];
            if (i<n-1) std::cout << ", ";
        }
        std::cout << std::endl;
        mxFree(var);
    
        readVar(nSyn,"n", pmat, file);
        readVar(tstep,"tstep", pmat, file);

        readArray(tmp,"sPSP", dimSize, arraySize, pmat, file);
        sPSP_ptr = new double[arraySize];
        memcpy((void *) sPSP_ptr, (void *)(mxGetPr(tmp)),arraySize*sizeof(double));
        pointer4d(sPSP,sPSP_ptr,dimSize);
        mxDestroyArray(tmp);

        nt = dimSize[3];
    
        readArray(tmp,"tmax", dimSize, arraySize, pmat, file);
        tMax_ptr = new double[arraySize];
        memcpy((void *) tMax_ptr, (void *)(mxGetPr(tmp)),arraySize*sizeof(double));
        pointer3d(tMax,tMax_ptr,dimSize);
        mxDestroyArray(tmp);

        readArray(tmp,"dendvleak", dimSize, arraySize, pmat, file);
        dendVleak_ptr = new double[arraySize];
        memcpy((void *) dendVleak_ptr, (void *)(mxGetPr(tmp)),arraySize*sizeof(double));
        pointer3d(dendVleak,dendVleak_ptr,dimSize);
        mxDestroyArray(tmp);

        readArray(tmp,"kv", dimSize, arraySize, pmat, file);
        kV_ptr = new double[arraySize];
        memcpy((void *) kV_ptr, (void *)(mxGetPr(tmp)),arraySize*sizeof(double));
        pointer6d(kV,kV_ptr,dimSize);
        mxDestroyArray(tmp);
    
        readArray(tmp,"kv0", dimSize, arraySize, pmat, file);
        kV0_ptr = new double[arraySize];
        memcpy((void *) kV0_ptr, (void *)(mxGetPr(tmp)),arraySize*sizeof(double));
        pointer4d(kV0,kV0_ptr,dimSize);
        mxDestroyArray(tmp);
    
        readArray(tmp,"vleakage", dimSize, arraySize, pmat, file);
        vLeak_ptr = new double[arraySize];
        memcpy((void *) vLeak_ptr, (void *)(mxGetPr(tmp)),arraySize*sizeof(double));
        pointer2d(vLeak,vLeak_ptr,dimSize);
        mxDestroyArray(tmp);

        readArray(tmp,"fireCap", dimSize, arraySize, pmat, file);
        double *dfireCap = new double[arraySize];
        fireCap_ptr = new int[arraySize];
        memcpy((void *) dfireCap, (void *)(mxGetPr(tmp)),arraySize*sizeof(double));
        for (i=0;i<arraySize; i++) {
            fireCap_ptr[i] = static_cast<int>(dfireCap[i]);
        }
        pointer2d(fireCap,fireCap_ptr,dimSize);
        mxDestroyArray(tmp);
        delete []dfireCap;

        readArray(tmp,"sf", dimSize, arraySize, pmat, file);
        double *dsfireCap = new double[arraySize];
        sfireCap_ptr = new int[arraySize];
        memcpy((void *) dsfireCap, (void *)(mxGetPr(tmp)),arraySize*sizeof(double));
        for (i=0;i<arraySize; i++) {
            sfireCap_ptr[i] = static_cast<int>(dsfireCap[i]);
        }
        pointer2d(sfireCap,sfireCap_ptr,dimSize);
        mxDestroyArray(tmp);
        delete []dsfireCap;

        readArray(tmp,"bf", dimSize, arraySize, pmat, file);
        double *dbfireCap = new double[arraySize];
        bfireCap_ptr = new int[arraySize];
        memcpy((void *) dbfireCap, (void *)(mxGetPr(tmp)),arraySize*sizeof(double));
        for (i=0;i<arraySize; i++) {
            bfireCap_ptr[i] = static_cast<int>(dbfireCap[i]);
        }
        pointer4d(bfireCap,bfireCap_ptr,dimSize);
        mxDestroyArray(tmp);
        delete []dbfireCap;

        readArray(tmp,"dendv", dimSize, arraySize, pmat, file);
        dendv_ptr = new double[arraySize];
        memcpy((void *) dendv_ptr, (void *)(mxGetPr(tmp)),arraySize*sizeof(double));
        pointer4d(dendv,dendv_ptr,dimSize);
        mxDestroyArray(tmp);

        readArray(tmp,"vRange", dimSize, arraySize, pmat, file);
        vRange = new double[arraySize];
        memcpy((void *) vRange, (void *)(mxGetPr(tmp)),arraySize*sizeof(double));
        //std::cout << "vRange:" << std::endl; disp1d(vRange,dimSize[0]);
        mxDestroyArray(tmp);
        
        nv = arraySize;
    
        readArray(tmp,"dtRange", dimSize, arraySize, pmat, file);
        dtRange = new double[arraySize];
        memcpy((void *) dtRange, (void *)(mxGetPr(tmp)),arraySize*sizeof(double));
        //std::cout << "dtRange:" << std::endl; disp1d(dtRange,dimSize[0]);
        mxDestroyArray(tmp);
        
        ndt = arraySize;

        readArray(tmp,"loc", dimSize, arraySize, pmat, file);
        double *dloc = new double[arraySize];
        loc = new int[arraySize];
        memcpy((void *) dloc, (void *)(mxGetPr(tmp)),arraySize*sizeof(double));
        //memcpy((void *) loc, (void *)(mxGetPr(tmp)),arraySize*sizeof(int64));
        for (i=0;i<arraySize; i++) {
            loc[i] = static_cast<int>(dloc[i]);
        }
        mxDestroyArray(tmp);
        delete []dloc;
        int nloc = arraySize;
        assert(nloc == nSyn);

        readArray(tmp,"pos", dimSize, arraySize, pmat, file);
        pos = new double[arraySize];
        memcpy((void *) pos, (void *)(mxGetPr(tmp)),arraySize*sizeof(double));
        mxDestroyArray(tmp);

        int npos = arraySize;
        assert(npos == nSyn);

        readArray(tmp,"gList", dimSize, arraySize, pmat, file);
        gList = new double[arraySize];
        memcpy((void *) gList, (void *)(mxGetPr(tmp)),arraySize*sizeof(double));
        mxDestroyArray(tmp);

        int nStrength = arraySize;
        assert(nStrength == nSyn);

        ei = new bool[nSyn];
        for (i=0; i<nSyn; i++) {
            if (gList[i] > 0) {
                ei[i] = true;
            } else {
                ei[i] = false;
            }
        }

        closeMat(pmat, file);
        assert(ndt>=2);
        assert(dtRange[0] == 0.0);
        tol_tb = dtRange[ndt-1];
        tol_tl = dtRange[ndt-1];

        std::cout << "Done" << std::endl;
        std::cout << std::endl;
        std::cout << "nSyn = " << nSyn << std::endl;
        std::cout << "nt = " << nt << std::endl;
        std::cout << "ndt = " << ndt << std::endl;
        std::cout << "nv = " << nv << std::endl;
        std::cout << "tstep = " << tstep << std::endl;

        std::cout << "vRange:" << std::endl; disp1d(vRange,nv);
        std::cout << "dtRange:" << std::endl; disp1d(dtRange,ndt);
        std::cout << "loc:" << std::endl; disp1d(loc,nSyn);
        std::cout << "pos:" << std::endl; disp1d(pos,nSyn);
        std::cout << "gList:" << std::endl; disp1d(gList,nSyn);
        std::cout << "ei:" << std::endl; disp1d(ei,nSyn);
        dimSize[0] = nSyn,
        dimSize[1] = ndt,
        std::cout << "fireCap:" << std::endl; disp2d(fireCap,dimSize);
        dimSize[0] = ndt,
        dimSize[1] = nSyn,
        std::cout << "sfireCap:" << std::endl; disp2d(sfireCap,dimSize);

        idtRange = new size[ndt];
        for (i=0; i<ndt; i++){
            idtRange[i] = static_cast<size>(round(dtRange[i]/tstep));
            assert(idtRange[i]*tstep-dtRange[i] < 1e-12);
        }
        for (i=1; i<ndt; i++) {
            if (nt - idtRange[i] < idtRange[ndt-1] - idtRange[i-1]) {
                std::cout << "check the " << i << "th entry of dtRange" << std::endl;
                assert(nt - idtRange[i] >= idtRange[ndt-1] - idtRange[i-1]);
            }
        }
    }
    // clear mem
    void clearLib() {
        size dimSize[5];

        dimSize[5] = nt;
        dimSize[4] = ndt;
        dimSize[3] = nSyn;
        dimSize[2] = nSyn;
        dimSize[1] = ndt;
        dimSize[0] = nv;
        del6d(kV,dimSize);
        delete []kV_ptr;
        cout << "kV unloaded" << endl;

        dimSize[2] = nSyn;
        dimSize[1] = ndt;
        dimSize[0] = nv;
        del3d(tMax,dimSize);
        delete []tMax_ptr;
        cout << "tMax unloaded" << endl;

        dimSize[3] = nt;
        dimSize[2] = nSyn;
        dimSize[1] = ndt;
        dimSize[0] = nv;
        del4d(sPSP,dimSize);
        delete []sPSP_ptr;
        cout << "sPSP unloaded" << endl;

        dimSize[3] = nt;
        dimSize[2] = nSyn;
        dimSize[1] = ndt;
        dimSize[0] = nv;
        del4d(dendv,dimSize);
        delete []dendv_ptr; 
        cout << "dendv unloaded" << endl;

        dimSize[1] = nt;
        dimSize[0] = nv;
        del2d(vLeak);
        delete []vLeak_ptr;
        cout << "vleak unloaded" << endl;

        dimSize[3] = ndt;
        dimSize[2] = nSyn;
        dimSize[1] = nSyn;
        dimSize[0] = ndt;
        del4d(bfireCap,dimSize);
        delete []bfireCap_ptr;
        cout << "bf unloaded" << endl;

        dimSize[1] = nSyn;
        dimSize[0] = ndt;
        del2d(sfireCap);
        delete []sfireCap_ptr;
        cout << "sf unloaded" << endl;
        
        dimSize[1] = ndt;
        dimSize[0] = nSyn;
        del2d(fireCap);
        delete []fireCap_ptr;
        cout << "fireCap unloaded" << endl;

        dimSize[3] = nt;
        dimSize[2] = nSyn;
        dimSize[1] = nSyn;
        dimSize[0] = ndt;
        del4d(kV0,dimSize);
        delete []kV0_ptr;
        cout << "kV0 unloaded" << endl;

        dimSize[2] = nt;
        dimSize[1] = nSyn;
        dimSize[0] = nv;
        del3d(dendVleak,dimSize);
        delete []dendVleak_ptr;
        cout << "dendVleak unloaded" << endl;

        delete []vRange;
        delete []dtRange;
        delete []idtRange;
        delete []loc;
        delete []pos;
        delete []gList;
        delete []ei;
    }
} nNL;
#endif
