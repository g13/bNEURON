#include "Yale_NEURON_PyAPI.h"


get_cell::get_cell(SynSet syn) {
    PyObject *pName, *pFunc;
    PyObject *pArgs, *pValue, *Py_gList;
    //PyObject *Py_loc, *Py_pos;
    _import_array();

    pName = PyString_FromString("neuroAlter");
    pModule = PyImport_Import(pName);
    Py_DECREF(pName);
    if (pModule != NULL) {
        cout << "neuroAlter loaded" << endl;
        pFunc = PyObject_GetAttrString(pModule,"prepCell");
        if (pFunc && PyCallable_Check(pFunc)) {
            cout << " prepCell function access granted" << endl;
            npy_intp dim = syn.nSyn;
            long n = static_cast<long>(dim);
            //cout << "dim: " << n << endl;
            //cout << syn.gList[0] << ", " << syn.gList[1] << ", " << syn.gList[2] << endl;
            Py_gList = PyArray_SimpleNewFromData(1,&dim,NPY_DOUBLE,syn.gList);
            if (!Py_gList) {
                cout << "gList Wrapper failed" << endl;
            } else {
                //cout << " gList ready" << endl;
            }
            Py_loc = PyArray_SimpleNewFromData(1,&dim,NPY_INT,syn.loc);
            if (!Py_loc) {
                cout << "loc Wrapper failed" << endl;
            } else {
                //cout << " loc ready" << endl;
            }
            Py_pos = PyArray_SimpleNewFromData(1,&dim,NPY_DOUBLE,syn.pos);
            if (!Py_pos) {
                cout << "pos Wrapper failed" << endl;
            } else {
                //cout << " pos ready" << endl;
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

//finialize
void NEURON_cleanup(Cell &cell) {
    Py_CLEAR(cell.pModule);
    Py_CLEAR(cell.Py_Cell);
    Py_CLEAR(cell.Py_synList);
    Py_CLEAR(cell.Py_vecStimList);
}

unsigned int Py_proceed(Cell &cell, double vinit, vector<vector<double>>
&RList, vector<long> &s1, vector<vector<double>> &spikeTrain, int n,
double trans, double tend, double vBack, double tref, double vThres, long
oneGo, vector<double> &v, long &nt, vector<double> &tsp, double t0,
double tstep, vector<double> &dendVclamp, long insert, bool getDendV=0,
vector<vector<double>> &dendV = dummy_dendV) {

    cout << "its here " << vinit << endl;
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
        if (yl::debug) {
            cout << "       module intact" << endl;
        }
        pFunc = PyObject_GetAttrString(cell.pModule,"proceed");
        if (yl::debug) {
            cout << "       proceed" << endl;
        }
        if (pFunc!=NULL && PyCallable_Check(pFunc)) {
            if (yl::debug) {
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
            }
            if (yl::debug) {
                cout << "       vecTuple = (";
            }
            for (long i=0; i<n; i++){
                dim = spikeTrain[i].size()-s1[i];
                if (yl::debug) {
                    assert(dim>0);
                }
                Py_spikeTrain[i] = PyArray_SimpleNewFromData(1,&dim,NPY_DOUBLE,&spikeTrain[i][s1[i]]);

                if (yl::debug) {
                    cout << "       np.array([";
                    for (long j=0; j<dim; j++) {
                        cout << spikeTrain[i][s1[i]+j];
                        if (j!=dim-1) cout << ", ";
                    }
                    cout << "])";
                    if (i != n-1) cout << ", ";
                    if (!Py_spikeTrain[i]) {
                        cout << endl;
                        cout << "-------> spikeTrain " << i << " wrapper failed" << endl;
                    }
                }
                
                PyTuple_SetItem(Py_spikeTrainTuple,i,Py_spikeTrain[i]);
            }
            if (yl::debug) {
                cout << ")" << endl; 
                cout << "       RList = [";
            }
            for (long i=0; i<n; i++){
                dim = RList[i].size();
                pRList[i] = PyArray_SimpleNewFromData(1,&dim,NPY_DOUBLE,RList[i].data());
                if (yl::debug) {
                    if (!pRList[i]) {
                        cout << "-----------> RList " << i << " wrapper failed" << endl;
                    }
                }
                PyTuple_SetItem(Py_RList,i,pRList[i]);
                if (yl::debug) {
                    cout << "[" << RList[i][0] << ", " << RList[i][1] << "]";
                    if (i != n-1) cout << ", ";
                }
            }
            if (yl::debug) {
                cout << "]" << endl;
                cout << "       dendVclamp = np.array([";
                for (long i=0; i<n; i++) {
                    cout << dendVclamp[i];
                    if (i!=n-1) cout << ", "; 
                }
                cout << "])" << endl;
            }

            dim = dendVclamp.size();     
            Py_dendVclamp = PyArray_SimpleNewFromData(1,&dim,NPY_DOUBLE,dendVclamp.data());
            if (yl::debug) {
                if (!Py_dendVclamp) {
                    cout << "-----------> vCap wrapper failed" << endl;
                }
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
            if (yl::debug) {
                cout << "       ..." << endl;
            }
            pValue = PyObject_CallObject(pFunc,pArgs);
            Py_DECREF(pArgs);
            if (yl::debug) {
                cout << " and back to earth"; 
            }
            if (pValue!=NULL) {
                PyObject *Py_v = PyTuple_GetItem(pValue,0);
                v_p = (double *) PyArray_DATA((PyArrayObject *)Py_v);
                nt = PyArray_SIZE(Py_v); 

                PyObject *Py_ntrans = PyTuple_GetItem(pValue,3);
                long gone = PyLong_AsLong(Py_ntrans);

                if (yl::debug) {
                    cout << " safe, crossed " <<  nt-gone <<  " planets, last trans included" << endl;
                }

                if (insert>=0) {
                    for (long i = gone; i<nt; i++) {
                        v[insert+gone-i] = v_p[i];
                    }
                } else {
                    for (long i = gone; i<nt; i++) {
                        v.push_back(v_p[i]);
                    }
                }
                if (yl::debug) {
                    cout << v_p[static_cast<int>(gone)] << ", " << v_p[static_cast<int>(gone)+1]<<"; " << v.back()  << "==" << v_p[nt-1] << "; " <<  v_p[static_cast<int>(gone)] - v_p[static_cast<int>(gone)+1] << endl;
                }

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
                    if (yl::debug) {
                        cout << " getting dend V " << endl;
                    }
                    PyObject *Py_dendV = PyTuple_GetItem(pValue,4);
                    v_p = (double *) PyArray_DATA((PyArrayObject *)Py_dendV);
                    for (int i=0; i<n; i++) {
                        for (long j=gone; j<nt; j++) {
                            dendV[i].push_back(v_p[i*nt+j]);
                        }
                    }
                }
                nt = nt - gone;

                Py_DECREF(pValue);
            } else {
                Py_DECREF(pFunc);
                if (yl::debug) {
                    char pause;
                    cout << " but, the universe has rejected us" << endl;
                    std::cin >> pause;
                }
            }
        }
        Py_XDECREF(pFunc);
    } else {
        if (yl::debug) {
            cout << " Houston we have a trouble!" << endl;
        }
    }
    //Py_Finalize();
    //cout << "finalized " << endl;
    return static_cast<unsigned int>(fired);
}

size neuroAlter(nNS &neuron, nNL &neuroLib, Cross &cross, size i_prior_cross, jND &jnd, double end_t, double it, double &tBack, double &vBack, double tstep, std::vector<double> &tsp, double vStop, unsigned int &nc, Cell &cell, vector<vector<double>> &spikeTrain, vector<long> &s0, vector<long> &s1, vector<double> &dendVclamp) {
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
    double tref = neuron.tRef;
    double trans = neuron.trans;
    vector<vector<double>> RList(neuron.nSyn,vector<double>(2,0));
    long nin = neuron.inID.size();
    
    get_RList(spikeTrain, s0, s1, t0, RList, neuron.ei, neuroLib.gList);

    npy_intp dim;
    //cross.v.push_back(vinit);
    //cross.t.push_back(it);
    if (cross.t.size() != cross.v.size()) {
        cout << cross.t.size() << "==" << cross.v.size() << endl;
        assert(cross.t.size() == cross.v.size());
    }
    size nt = 0;
    nc = nc + Py_proceed(cell, vinit, RList, s1, spikeTrain, neuron.nSyn, trans, trans + end_t, vStop, tref, vThres, 0, cross.v, nt, tsp, t0, tstep, dendVclamp, -1, false); 
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

size neuroAlterB(nNS &neuron, nNL &neuroLib, vector<double> &v, size &ith, size &vs, size run_nt, double tstep, vector<double> &tsp, double vinit, double vBack, Cell &cell, vector<vector<double>> &spikeTrain, vector<long> &s0, vector<long> &s1, vector<double> &dendVclamp) {
    size i,j,k;
    double t = vs*tstep;
    double vThres = neuron.vThres;
    double tref = neuron.tRef;
    double trans = neuron.trans;
    size lCross;
    vector<vector<double>> RList(neuron.nSyn,vector<double>(2,0));
    long nin = neuron.inID.size();

    if (yl::debug) {
        cout << "houston, we are now inside the cave, will retreat if temp goes below: " << vBack << endl;
    }
    for (i=0; i< neuron.nSyn; i++){
        while (spikeTrain[i][s0[i]] < t - neuroLib.tol_tl) {
            s0[i]++; 
        }
        if (s1[i] < s0[i]) {
            s1[i] = s0[i] + 1;
        }
        while (spikeTrain[i][s1[i]] < t) {
            s1[i]++;
        }
        if (yl::debug) {
            cout << "   s0 " << s0[i] << "s1 " << s1[i] << endl;
            assert(s1[i] < spikeTrain[i].size());
        }
    }
    if (yl::debug) {
        cout << "   rope's ends tightened" << endl;
    }

    get_RList(spikeTrain, s0, s1, t, RList, neuron.ei, neuroLib.gList);

    if (yl::debug) {
        cout << "   rope is ready" << endl;
    }
    double tend = trans + (run_nt-1)*tstep;
    npy_intp dim;
    size nc = Py_proceed(cell, vinit, RList, s1, spikeTrain, neuron.nSyn, trans, tend, vBack, tref, vThres, 0, v, lCross, tsp, t, tstep, dendVclamp,vs, false); 

    vs = vs + lCross - 1;
    return nc;
}
