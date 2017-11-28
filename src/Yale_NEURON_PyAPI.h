#ifndef YALE_NEURON_PY_API
#define YALE_NEURON_PY_API
#include <Python.h>
#include "nNeuroSt.h"
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

inline void get_gh(vector<double> &spikeTrain, long j0, long j1, double t, double &g, double &h, double t0, double t1, double f);

inline void get_RList(vector<vector<double>> &spikeTrain, vector<long> &j0, vector<long> &j1, double t0, vector<vector<double>> &RList, vector<bool> &ei, double* gList);

struct SynSet{
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
    get_cell(synSet syn);
};
typedef struct get_cell Cell;

//finialize
void NEURON_cleanup(Cell &cell); 

inline unsigned int Py_proceed(Cell &cell, double vinit, vector<vector<double>> &RList, vector<long> &s1, vector<vector<double>> &spikeTrain, int n, double trans, double tend, double vBack, double tref, double vThres, long oneGo, vector<double> &v, long &nt, vector<double> &tsp, double t0, double tstep, vector<double> &dendVclamp, bool getDendV = 0, vector<vector<double>> &dendV = dummy_dendV); 

inline size neuroAlter(nNS &neuron, nNL &neuroLib, Cross &cross, size i_prior_cross, jND &jnd, double end_t, double it, double &tBack, double &vBack, double tstep, std::vector<double> &tsp, double vStop, unsigned int &nc, Cell &cell, vector<bool> &ei, vector<vector<double>> &spikeTrain, vector<long> &s0, vector<long> &s1, vector<double> &dendVclamp);

inline size neuroAlterB(nNS &neuron, nNL &neuroLib, vector<double> &v, vector<double> &crossv, size &ith, size vs, size run_nt, double tstep, vector<double> &tsp, double vinit, size &ve, double vBack, Cell cell, vector<bool> &ei, vector<vector<double>> &spikeTrain, vector<long> &s0, vector<long> &s1, vector<double> &dendVclamp);

#endif
