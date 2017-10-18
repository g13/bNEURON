#ifndef NEUROALTER_H
#define NEUROALTER_H
#include <Python.h>
#include "nNeuroSt2.h"
#include "typedefs.h"
#include <vector>
#include <iostream>
#include <numpy/arrayobject.h>
#include <math.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
using std::vector;
using std::cout;
using std::endl;

const double tau_er = 0.098814229249;
const double tau_ed = 8.333333333333;
const double tau_ir = 0.98814229249;
const double tau_id = 50.0;
static vector<vector<double>> dummy_dendV;

inline void getgh(vector<double> &spikeTrain, long j0, long j1, double t, double &g, double &h, double t0, double t1, double f) {
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
inline void getRList(vector<vector<double>> &spikeTrain, vector<long> &j0, vector<long> &j1, double t0, vector<vector<double>> &RList, vector<bool> &ei, double* gList){
    //cout << "ei size " << ei.size() << endl;
    for (long i=0; i<ei.size(); i++){
        //cout << " tail " << j0[i] << " head " << j1[i] << endl;
        //cout << " ei " << ei[i] << " f " <<  gList[i] << endl;
        if (ei[i]>0) {
            getgh(spikeTrain[i], j0[i], j1[i], t0, RList[i][0], RList[i][1], tau_er, tau_ed, gList[i]);
        } else {
            getgh(spikeTrain[i], j0[i], j1[i], t0, RList[i][0], RList[i][1], tau_ir, tau_id,-gList[i]);
        }
    //     cout << "[" << RList[i][0] << ", " << RList[i][1] << "]" << endl;
    }
}
typedef struct SynSet{
    int *loc;
    double *pos;
    double *gList;
    size nSyn;
    double v0;
} synSet;

typedef struct getCell{
    PyObject *Py_Cell;
    PyObject *Py_synList;
    PyObject *Py_vecStimList;
    PyObject *pModule;
    PyObject *Py_loc;
    PyObject *Py_pos;
    getCell(synSet syn) {
        PyObject *pName, *pFunc;
        PyObject *pArgs, *pValue, *Py_gList;
        //PyObject *Py_loc, *Py_pos;
    
        pName = PyString_FromString("neuroAlter");
        pModule = PyImport_Import(pName);
        Py_DECREF(pName);
        if (pModule != NULL) {
            cout << "neuroAlter loaded" << endl;
            pFunc = PyObject_GetAttrString(pModule,"prepCell");
            cout << " prepCell function access granted" << endl;
            if (pFunc && PyCallable_Check(pFunc)) {
                npy_intp dim = syn.nSyn;
                long n = static_cast<long>(dim);
                Py_gList = PyArray_SimpleNewFromData(1,&dim,NPY_DOUBLE,syn.gList);
                if (!Py_gList) {
                    cout << "gList Wrapper failed" << endl;
                }
                Py_loc = PyArray_SimpleNewFromData(1,&dim,NPY_INT,syn.loc);
                if (!Py_loc) {
                    cout << "loc Wrapper failed" << endl;
                }
                Py_pos = PyArray_SimpleNewFromData(1,&dim,NPY_DOUBLE,syn.pos);
                if (!Py_pos) {
                    cout << "pos Wrapper failed" << endl;
                }
                pArgs = PyTuple_New(5);
                PyTuple_SetItem(pArgs,0,Py_gList);
                Py_INCREF(Py_loc);
                PyTuple_SetItem(pArgs,1,Py_loc);
                Py_INCREF(Py_pos);
                PyTuple_SetItem(pArgs,2,Py_pos);
                PyTuple_SetItem(pArgs,3,PyInt_FromLong(n));
                PyTuple_SetItem(pArgs,4,PyFloat_FromDouble(syn.v0));
                cout << "Args ready" << endl;
                pValue = PyObject_CallObject(pFunc,pArgs);
                cout << "prepCell called" << endl;
                Py_DECREF(pArgs);
                if (pValue!=NULL) {
                    Py_Cell = PyTuple_GetItem(pValue,0);
                    Py_INCREF(Py_Cell);
                    Py_vecStimList = PyTuple_GetItem(pValue,1);
                    Py_INCREF(Py_vecStimList);
                    Py_synList = PyTuple_GetItem(pValue,2);
                    Py_INCREF(Py_synList);
                    Py_DECREF(pValue);
                } else {
                    Py_DECREF(pFunc);
                    cout << "prepCell failed"<< endl;
                }
            }
            Py_XDECREF(pFunc);
        } else {
            cout << " Houston we lost our module" << endl;
        }
    }
} Cell;

//finialize
void NEURON_cleanup(Cell &cell) {
    Py_CLEAR(cell.pModule);
    Py_CLEAR(cell.Py_Cell);
    Py_CLEAR(cell.Py_synList);
    Py_CLEAR(cell.Py_vecStimList);
}
inline unsigned int Py_proceed(Cell &cell, double vinit, vector<vector<double>> &RList, vector<long> &s1, vector<vector<double>> &spikeTrain, int n, double trans, double tend, double vBack, double tref, double vThres, long oneGo, vector<double> &v, long &nt, vector<double> &tsp, double t0, double tstep, vector<double> &dendVclamp, bool getDendV = 0, vector<vector<double>> &dendV = dummy_dendV) {

    PyObject **pRList = new PyObject*[n];
    PyObject *pFunc;
    PyObject *pArgs, *pValue, *Py_dendVclamp;
    PyObject **Py_spikeTrain = new PyObject*[n];
    PyObject *Py_spikeTrainTuple = PyTuple_New(n);
    PyObject *Py_RList = PyTuple_New(n);
    long fired;
    double *v_p, *tsp_p;
    npy_intp dim;
    int countdown = 10;

    if (cell.pModule != NULL) {
        cout << "module intact" << endl;
        pFunc = PyObject_GetAttrString(cell.pModule,"proceed");
        cout << "proceed " << endl;
        if (pFunc!=NULL && PyCallable_Check(pFunc)) {
            if (cell.Py_Cell != NULL) {
                cout << countdown-- << " ";
            } else {
                cout << " no Cell " << endl;
            }
            if (cell.Py_synList != NULL) {
                cout << countdown-- << " ";
            } else {
                cout << " no synapse " << endl;
            }
            if (cell.Py_vecStimList != NULL) {
                cout << countdown-- << " ";
            } else {
                cout << " no vecStim " << endl;
            }
            if (cell.Py_loc != NULL) {
                cout << countdown-- << " ";
            } else {
                cout << " no dend loc " << endl;
            }
            if (cell.Py_pos != NULL) {
                cout << countdown-- << " ";
            } else {
                cout << " no dend pos " << endl;
            }
            cout << "cell prepared" << endl;
            cout << endl;
            cout << "vecTuple = (";
            for (long i=0; i<n; i++){
                dim = spikeTrain[i].size()-s1[i];
                assert(dim>0);
                Py_spikeTrain[i] = PyArray_SimpleNewFromData(1,&dim,NPY_DOUBLE,&spikeTrain[i][s1[i]]);

                cout << "np.array([";
                for (long j=0; j<dim; j++) {
                    cout << spikeTrain[i][s1[i]+j];
                    if (j!=dim-1) cout << ", ";
                }
                cout << "])";
                if (i != n-1) cout << ", ";

                if (!Py_spikeTrain[i]) {
                    cout << "spikeTrain " << i << " wrapper failed" << endl;
                }
                
                PyTuple_SetItem(Py_spikeTrainTuple,i,Py_spikeTrain[i]);
            }
            cout << ")" << endl; 
            cout << "RList = [";
            for (long i=0; i<n; i++){
                dim = RList[i].size();
                pRList[i] = PyArray_SimpleNewFromData(1,&dim,NPY_DOUBLE,RList[i].data());
                if (!pRList[i]) {
                    cout << "RList " << i << " wrapper failed" << endl;
                }
                PyTuple_SetItem(Py_RList,i,pRList[i]);
                cout << "[" << RList[i][0] << ", " << RList[i][1] << "]";
                if (i != n-1) cout << ", ";
            }
            cout << "]" << endl;
            cout << "dendVclamp = np.array([";
            for (long i=0; i<n; i++) {
                cout << dendVclamp[i];
                if (i!=n-1) cout << ", "; 
            }
            cout << "])" << endl;

            dim = dendVclamp.size();     
            Py_dendVclamp = PyArray_SimpleNewFromData(1,&dim,NPY_DOUBLE,dendVclamp.data());
            if (!Py_dendVclamp) {
                cout << "vCap wrapper failed" << endl;
            }
            pArgs = PyTuple_New(21);
            Py_INCREF(cell.Py_Cell);
            PyTuple_SetItem(pArgs, 0,cell.Py_Cell);
            PyTuple_SetItem(pArgs, 1,PyFloat_FromDouble(vinit));
            Py_INCREF(cell.Py_synList);
            PyTuple_SetItem(pArgs, 2,cell.Py_synList);
            PyTuple_SetItem(pArgs, 3,Py_RList);
            Py_INCREF(cell.Py_vecStimList);
            PyTuple_SetItem(pArgs, 4,cell.Py_vecStimList);
            PyTuple_SetItem(pArgs, 5,Py_spikeTrainTuple);
            PyTuple_SetItem(pArgs, 6,PyInt_FromLong(n));
            PyTuple_SetItem(pArgs, 7,PyFloat_FromDouble(trans));
            PyTuple_SetItem(pArgs, 8,PyFloat_FromDouble(tend));
            PyTuple_SetItem(pArgs, 9,PyFloat_FromDouble(vBack));
            PyTuple_SetItem(pArgs,10,PyFloat_FromDouble(tref));
            PyTuple_SetItem(pArgs,11,PyFloat_FromDouble(vThres));
            PyTuple_SetItem(pArgs,12,PyInt_FromLong(oneGo));
            PyTuple_SetItem(pArgs,13,PyFloat_FromDouble(t0));
            PyTuple_SetItem(pArgs,14,PyFloat_FromDouble(tstep));
            Py_INCREF(cell.Py_loc);
            PyTuple_SetItem(pArgs,15,cell.Py_loc);
            Py_INCREF(cell.Py_pos);
            PyTuple_SetItem(pArgs,16,cell.Py_pos);
            PyTuple_SetItem(pArgs,17,Py_dendVclamp);
            PyTuple_SetItem(pArgs,18,Py_False);
            PyTuple_SetItem(pArgs,19,Py_True);
            if (getDendV) {
                PyTuple_SetItem(pArgs,20,Py_True);
            } else {
                PyTuple_SetItem(pArgs,20,Py_False);
            }
            cout << "..." << endl;
            pValue = PyObject_CallObject(pFunc,pArgs);
            Py_DECREF(pArgs);
            cout << " and back to earth"; 
            if (pValue!=NULL) {
                PyObject *Py_v = PyTuple_GetItem(pValue,0);
                v_p = (double *) PyArray_DATA((PyArrayObject *)Py_v);
                nt = PyArray_SIZE(Py_v); 

                PyObject *Py_ntrans = PyTuple_GetItem(pValue,3);
                long gone = PyLong_AsLong(Py_ntrans);

                cout << " safe, crossed " <<  nt-gone <<  " planets, last trans included" << endl;

                for (long i = gone; i<nt; i++) {
                    v.push_back(v_p[i]);
                }
                cout << v_p[static_cast<int>(gone)] << ", " << v_p[static_cast<int>(gone)+1]<<"; " << v.back()  << "==" << v_p[nt-1] << "; " <<  v_p[static_cast<int>(gone)] - v_p[static_cast<int>(gone)+1] << endl;

                PyObject *Py_fired = PyTuple_GetItem(pValue,1);
                fired = PyLong_AsLong(Py_fired);

                PyObject *Py_tsp = PyTuple_GetItem(pValue,2);
                long nspike = PyArray_SIZE(Py_tsp); 
                assert(fired==nspike);
                tsp_p = (double*) PyArray_DATA((PyArrayObject *)Py_tsp);
                for (long i = 0; i<nspike; i++) {
                    tsp.push_back(tsp_p[i]);
                }
                
                if (getDendV) {
                    cout << " getting dend V " << endl;
                    PyObject *Py_dendV = PyTuple_GetItem(pValue,4);
                    v_p = (double *) PyArray_DATA((PyArrayObject *)Py_dendV);
                    for (int i=0; i<n; i++) {
                        //cout << "i " << i << endl;
                        for (long j=gone; j<nt; j++) {
                            //cout << "j " << j << ", " << i*nt+j <<  endl;
                            dendV[i].push_back(v_p[i*nt+j]);
                        }
                    }
                }
                nt = nt - gone;

                Py_DECREF(pValue);
            } else {
                Py_DECREF(pFunc);
                char pause;
                cout << " but, the universe has rejected us" << endl;
                std::cin >> pause;
            }
        }
        Py_XDECREF(pFunc);
    } else {
        cout << " Houston we have a trouble!" << endl;
    }
    //Py_Finalize();
    //cout << "finalized " << endl;
    return static_cast<unsigned int>(fired);
}
inline size neuroAlter(nNS &neuron, nNL &neuroLib, Cross &cross, size i_prior_cross, jND &jnd, double end_t, double it, double &tBack, double &vBack, double tstep, std::vector<double> &tsp, double vStop, unsigned int &nc, Cell &cell, vector<bool> &ei, vector<vector<double>> &spikeTrain, vector<long> &s0, vector<long> &s1, vector<double> &dendVclamp) {
    if (it!=round(it)) {
        cout << "it needs roud off" << it << " != " << round(it) << endl;
        it = round(it);
        assert(it == round(it));
    }
    double t0 = it*tstep;
    int i, lastSize = cross.v.size();
    long appendSize;

    double vinit = jnd.v.back();
    double vThres = neuron.vThres;
    double tref = neuron.tref;
    double trans = neuron.trans;
    vector<vector<double>> RList(neuron.nSyn,vector<double>(2,0));
    long nin = neuron.inID.size();
    
    getRList(spikeTrain, s0, s1, t0, RList, ei, neuroLib.gList);

    npy_intp dim;
    //cross.v.push_back(vinit);
    //cross.t.push_back(it);
    if (cross.t.size() != cross.v.size()) {
        cout << cross.t.size() << "==" << cross.v.size() << endl;
        assert(cross.t.size() == cross.v.size());
    }
    size nt = 0;
    nc = nc + Py_proceed(cell, vinit, RList, s1, spikeTrain, neuron.nSyn, trans, trans + end_t, vStop, tref, vThres, 0, cross.v, nt, tsp, t0, tstep, dendVclamp); 
    assert(nt == cross.v.size() - lastSize);
    tBack = round(it + nt-1);
    assert(int(tBack) == tBack);
    size ith = i_prior_cross + 1;
    while (neuron.tin[ith]-1e-10 < tBack*tstep+1e-10) {
        ith++;
        if ( ith == nin) break;
    }
    ith--;
    //std::cout << "back at i " << ith << std::endl;
    for (i=it; i<=tBack; i++) {
        cross.t.push_back(i);
    }
    if (cross.t.size() != cross.v.size()) {
        cout << cross.t.size() << "==" << cross.v.size() << endl;
        assert(cross.t.size() == cross.v.size());
    }
    cross.tCross.push_back(tBack);
    cross.iCross.push_back(cross.v.size());
    cross.nCross++;
    cout << cross.nCross << "th cross" << endl;
    vBack = cross.v.back(); 
    return ith;
}

inline size neuroAlterB(nNS &neuron, nNL &neuroLib, vector<double> &v, vector<double> &crossv, size &ith, size vs, size run_nt, double tstep, vector<double> &tsp, double vinit, size &ve, double vBack, Cell cell, vector<bool> &ei, vector<vector<double>> &spikeTrain, vector<long> &s0, vector<long> &s1, vector<double> &dendVclamp) {
    //void *Cell, *vecStimList, *synList;
    size i,j,k;
    double t = vs*tstep;
    double vThres = neuron.vThres;
    double tref = neuron.tref;
    double trans = neuron.trans;
    size endCross, lCross, lastCross = crossv.size() -1;
    vector<vector<double>> RList(neuron.nSyn,vector<double>(2,0));
    long nin = neuron.inID.size();

    //cout << " houston, we are now inside the cave" << endl;
    //cout << " will be back if temp goes below: " << vBack << endl;
    for (i=0; i< neuron.nSyn; i++){
        while (spikeTrain[i][s0[i]] < t - neuroLib.tol_tl) {
            s0[i]++; 
        }
        assert(s0[i] < spikeTrain[i].size());
        while (spikeTrain[i][s1[i]] < t) {
            s1[i]++;
        }
        assert(s1[i] < spikeTrain[i].size());
        cout << "s0 " << s0[i] << "s1 " << s1[i] << endl;
    }
    //cout << " ropes ready" << endl;

    getRList(spikeTrain, s0, s1, t, RList, ei, neuroLib.gList);

    //cout << " ropes tightened " << endl;
    double tend = trans + (run_nt-1)*tstep;
    npy_intp dim;
    size nc = Py_proceed(cell, vinit, RList, s1, spikeTrain, neuron.nSyn, trans, tend, vBack, tref, vThres, 0, crossv, lCross, tsp, t, tstep, dendVclamp); 
    assert(crossv.size()-(lastCross+1) == lCross);
    //cout << "lastCross = " << lastCross;
    for (i=0; i<lCross; i++) {
        if (vs+i>=v.size()) {
             lCross = i;
             break;
        }
        //cout << lastCross+i+1 << ", " << crossv.at(lastCross+i+1) << endl;
        v[vs+i] = crossv.at(lastCross+i+1);
    }

    for (ith = ith; ith < nin; ith++) {
        if (neuron.tin[ith] > t+(lCross-1)*tstep) {
            break;                
        }
    }
    ith--;

    ve = vs + lCross - 1;
    assert(neuron.tin[ith] <= t+(lCross-1)*tstep);
    return nc;
}
#endif
