#ifndef YALE_NEURON_PY_API
#define YALE_NEURON_PY_API
#include "nNeuroSt.h"
#include <Python.h>
#include <vector>
#include <iostream>
#include <numpy/arrayobject.h>
#include <cmath>

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
using std::vector;
using std::cout;
using std::endl;

const double tau_er = 0.098814229249;
const double tau_ed = 8.333333333333;
const double tau_ir = 0.98814229249;
const double tau_id = 50.0;
static vector<vector<double>> dummy_dendV;

inline void get_gh(vector<double> &spikeTrain, long j0, long j1, double t, double &g, double &h, double t0, double t1, double f) {
    double nextRel, lastRel = -1e23;
    double etd, etr, c, dt, g0;
    g0 = g;
    for (long j=j0; j<j1; j++) {
        nextRel = spikeTrain[j];
        dt = nextRel - lastRel;
        etd = exp(-dt/t1);
        etr = exp(-dt/t0);
        c = t0/(t1-t0)*(etd-etr);
        g = g*etd + c*h;
        h = h*etr;
        lastRel = nextRel;
        h = h + f/t0;
    }
    dt = t - lastRel;
    if (dt < 0) {
        cout << spikeTrain[j0] << "->" << spikeTrain[j1] << "| t " << t << endl;
        assert(dt > 0);
    }
    etd = exp(-dt/t1);
    etr = exp(-dt/t0);
    c = t0/(t1-t0)*(etd-etr);
    g = g*etd + c*h;
    h = h*etr;
    assert(f>0);
    assert(h>=0);
    assert(g>=0);
}
inline void get_RList(vector<vector<double>> &spikeTrain, vector<long> &j0, vector<long> &j1, double t0, vector<vector<double>> &RList, vector<bool> &ei, double* gList){
    //cout << "ei size " << ei.size() << endl;
    for (long i=0; i<ei.size(); i++){
        //cout << " tail " << j0[i] << " head " << j1[i] << endl;
        //cout << " ei " << ei[i] << " f " <<  gList[i] << endl;
        if (ei[i]>0) {
            get_gh(spikeTrain[i], j0[i], j1[i], t0, RList[i][0], RList[i][1], tau_er, tau_ed, gList[i]);
        } else {
            get_gh(spikeTrain[i], j0[i], j1[i], t0, RList[i][0], RList[i][1], tau_ir, tau_id,-gList[i]);
        }
    //     cout << "[" << RList[i][0] << ", " << RList[i][1] << "]" << endl;
    }
}

struct syn_set{
    int *loc;
    double *pos;
    double *gList;
    size nSyn;
    double v0;
};
typedef struct syn_set SynSet;

struct get_cell{
    PyObject *Py_Cell;
    PyObject *Py_synList;
    PyObject *Py_vecStimList;
    PyObject *pModule;
    PyObject *Py_loc;
    PyObject *Py_pos;
    get_cell(SynSet syn);
};
typedef struct get_cell Cell;

//finialize
void NEURON_cleanup(Cell &cell); 

unsigned int Py_proceed(Cell &cell, double vinit, vector<vector<double>> &RList, vector<long> &s1, vector<vector<double>> &spikeTrain, int n, double trans, double tend, double vBack, double tref, double vThres, long oneGo, vector<double> &v, long &nt, vector<double> &tsp, double t0, double tstep, vector<double> &dendVclamp, bool getDendV = 0, vector<vector<double>> &dendV = dummy_dendV); 

size neuroAlter(nNS &neuron, nNL &neuroLib, Cross &cross, size i_prior_cross, jND &jnd, double end_t, double it, double &tBack, double &vBack, double tstep, std::vector<double> &tsp, double vStop, unsigned int &nc, Cell &cell, vector<bool> &ei, vector<vector<double>> &spikeTrain, vector<long> &s0, vector<long> &s1, vector<double> &dendVclamp);

size neuroAlterB(nNS &neuron, nNL &neuroLib, vector<double> &v, vector<double> &crossv, size &ith, size vs, size run_nt, double tstep, vector<double> &tsp, double vinit, size &ve, double vBack, Cell cell, vector<bool> &ei, vector<vector<double>> &spikeTrain, vector<long> &s0, vector<long> &s1, vector<double> &dendVclamp);

#endif
