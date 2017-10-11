#ifndef NNS_H
#define NNS_H
#include "nNeuroLib.h"
#include <vector>
#include <random>
#include <fstream>
#include <array>
#include "heap.h"

using std::vector;
using std::cout;
using std::endl;

struct IJR {
    // func(x) = func(i) + r*(func(j)-func(i));
    size i;
    size j;
    double r;
    IJR(){};
    IJR(double i0, double j0, double r0) : i(i0), j(j0), r(r0) { };
    IJR operator+(struct IJR ijr) {
        struct IJR ijr0;
        double r0 = this->r + ijr.r;
        size ir = static_cast<size>(r0);
        ijr0.i = this->i + ijr.i + ir;
        ijr0.r = r0-ir; 
        ijr0.j = ijr0.i+1;
        return ijr0;
    }
    IJR operator+(double r0) {
        struct IJR ijr0;
        double r1 = this->r + r0;
        size ir = static_cast<size>(r1);
        ijr0.i = this->i + ir;
        ijr0.r = r1-ir; 
        ijr0.j = ijr0.i+1;
        return ijr0;
    }
    IJR operator+(size i0) {
        struct IJR ijr0;
        ijr0.i = this->i + i0;
        ijr0.r = this->r; 
        ijr0.j = ijr0.i+1;
        return ijr0;
    }
};
typedef struct jumpyNeuronData {
    std::vector<double> t,v;
    jumpyNeuronData(size rSize) {
        t.reserve(2*rSize);
        v.reserve(2*rSize);
    }
} jND;
typedef struct CrossData {
    // index 0 for vinit
    // counting from 1.
    size nCross;
    std::vector<size> iCross; // next Cross index
    std::vector<struct IJR> vCross;
    std::vector<double> tCross; // t of crossing back
    std::vector<double> t,v,gE,gI,m,n,h,hE,hI;
    CrossData(size nt, double vinit) {
        nCross = 0;
        iCross.reserve(nt/2);
        tCross.reserve(nt/2);
        vCross.reserve(nt/2);
        tCross.push_back(0);
        v.reserve(nt);
        t.reserve(nt);
        iCross.push_back(1);
        v.push_back(vinit);
        t.push_back(0);
    }
} Cross;
typedef struct BilinearRelationships {
    // dt to the second input;
    std::vector<struct IJR> dTijr; 
    std::vector<size> idt; 
    std::vector<unsigned int> ID;
    BilinearRelationships(size corrSize) {
        dTijr.reserve(corrSize);
        idt.reserve(corrSize);
        ID.reserve(corrSize);
    }
} biR;
typedef struct Inputs {
    // 0 for sample tstep;
    // 1 for demanded tstep;
    std::vector<struct IJR> Vijr;
    // dT for cross and K
    std::vector<struct IJR> dTijr,Tmpijr;
    // current index of effecting cross
    std::vector<size> cCross;
    std::vector<double> t,dt,tMax;
    std::vector<biR> bir;
    std::vector<size> ID;
    std::vector<bool> inTref;
    Inputs(size rSize){
        t.reserve(rSize); 
        dt.reserve(rSize); 
        tMax.reserve(rSize); 
        cCross.reserve(rSize);
        dTijr.reserve(rSize); 
        Vijr.reserve(rSize);
        ID.reserve(rSize);
        Tmpijr.reserve(rSize); 
        bir.reserve(rSize);
        inTref.reserve(rSize);
    }
    void junk_this() {
        t.push_back(0); 
        ID.push_back(-1);
        dt.push_back(0);
        tMax.push_back(0); 
        cCross.push_back(0);
        dTijr.push_back(IJR(0,1,0)); 
        Vijr.push_back(IJR(0,1,0));
        Tmpijr.push_back(IJR(0,1,0)); 
        bir.push_back(BilinearRelationships(0));
        inTref.push_back(0);
    }   
    void assert_size() {
        size Size = t.size(); 
        assert(dt.size() == Size); 
        assert(tMax.size() == Size); 
        assert(cCross.size() == Size);
        assert(dTijr.size() == Size); 
        if (Vijr.size() != Size) {
            cout << Vijr.size() << "!=" << Size << endl;
            assert(Vijr.size() == Size);
        }
        assert(Tmpijr.size() == Size); 
        assert(bir.size() == Size);
        assert(inTref.size() == Size);
    }   
    void print_this(int i) {
        cout << "t: " << t[i] << endl;
        cout << "ID: " << ID[i] << endl;
        cout << "V index: " << Vijr[i].i << ", " << Vijr[i].j << ", " << Vijr[i].r << endl;
        cout << "dT index: " << dTijr[i].i << ", " << dTijr[i].j << ", " << dTijr[i].r << endl;
    }
} Input;

typedef struct nNeuroSt{
    double v_tol, run_t, trans;
    double vThres, vReset, vRest, tref, tstep;
    unsigned int nSyn;
    vector<double> tin;
    vector<size> inID; 
    vector<double> tsp;
    size nOut;
	vector<std::minstd_rand> poiGen, ranGen;
    vector<vector<double>> tPoi;
    vector<double> tmp;
    vector<int> heap;
    bool status;
    vector<bool> ei;
    std::uniform_real_distribution<double> uniform0_1 = std::uniform_real_distribution<double>(0.0,1.0);;

    void initialize(unsigned int seed, int nSyn0, vector<double> &t0, double run_t0, vector<double> &rate, vector<double> &max_rate, bool *ei0, bool testing, double tstep0) {
        status = true;
        nOut = 0;
        run_t = run_t0;
        nSyn = nSyn0;
        tstep = tstep0;
        ei.assign(ei0,ei0+nSyn);
            
        poiGen.assign(nSyn,std::minstd_rand());
        ranGen.assign(nSyn,std::minstd_rand());
        tPoi.assign(nSyn,vector<double>());
        //cout << "E or I :";
        for (int i = 0; i < nSyn; i++ ) {
            //cout << ei[i] << ", ";
            ranGen[i].seed(seed+i+1);
            poiGen[i].seed(seed-i-1);
        }
        cout << endl;
        if (!testing) {
            cout << "initialized input heap" << endl;
            double tmpT;
            double dt;
            for (int i=0; i<nSyn; i++) {
                if (rate[i] > 0 ) {
                    tmpT = t0[i];
                    do tmpT = tmpT -log(uniform0_1(poiGen[i]))/max_rate[i];
                    while (uniform0_1(ranGen[i]) > rate[i]/max_rate[i]);
                    //cout << " uniform " <<uniform0_1(poiGen[i]) << endl;
                    //cout << " t0[i] " <<t0[i] << "r[i] " << rate[i] << "rm[i]" << max_rate[i] << endl;
                    //cout << "   " << tmpT << " < " << run_t << endl;   
                    tmpT = static_cast<int>(tmpT/tstep)*tstep;
                } else {
                    tmpT = run_t + 1.0;
                }
                tmp.push_back(tmpT);
                tPoi[i].push_back(tmp.back());
            }
            my::make_heap(tmp.data(), nSyn, heap);
        }
    }
    void getNextInput(vector<double> rate, vector<double> max_rate) {
        int i =  heap[0];
        double dt;
        if (tmp[i] > run_t) {
            status = false;
            return;
        } else {
            //cout << tmp[i] << endl;
            tin.push_back(tmp[i]);
            inID.push_back(i);
        }
        double tmpT = tPoi[i].back();
        do {
            dt = -log(uniform0_1(poiGen[i]))/max_rate[i];
            tmpT = tmpT + dt;
            tmpT = round(tmpT/tstep)*tstep;
        } while (uniform0_1(ranGen[i]) > rate[i]/max_rate[i] || dt < tstep);
        tmp[i] = tmpT;
        tPoi[i].push_back(tmpT);
        my::heapify(tmp.data(),0,tmp.size(),heap);
    }
    void setInputs(vector<vector<double>> &inputs) {
        int index;
        vector<int> stillHas(nSyn,0);
        for (int i=0; i<nSyn; i++) {
            if (i < inputs.size() && inputs[i].size() > stillHas[i]){
                tmp.push_back(inputs[i][stillHas[i]]); 
                stillHas[i] += 1;
            } else {
                tmp.push_back(run_t+1);
            }
            tPoi[i].push_back(tmp[i]);
        }
        my::make_heap(tmp.data(), nSyn, heap);
        cout << "initialized input heap" << endl;
        do {
            index = heap[0];
            if (tmp[index] > run_t) {
                status = false;
                return;
            } else {
                tin.push_back(tmp[index]);
                inID.push_back(index);
            }
            if (inputs[index].size() > stillHas[index]) {
                tmp[index] = inputs[index][stillHas[index]];
                stillHas[index] += 1;
            } else {
                tmp[index] = run_t+1;
            }
            tPoi[index].push_back(tmp[index]);
            my::heapify(tmp.data(),0,tmp.size(),heap);
        } while (status);
    }
    void writeAndUpdateIn(size i, std::ofstream &tIncome_file) {
        // write 
        tIncome_file.write((char*)&i, sizeof(i));
        if (i!=0) {
            tIncome_file.write((char*)&(tin[0]), i * sizeof(tin[0]));
            tIncome_file.write((char*)&(inID[0]), i * sizeof(inID[0]));
        }
        size keeping = tin.size()-i;
        if (keeping) {
            tin.assign(tin.end()-keeping, tin.end());
            inID.assign(inID.end()-keeping, inID.end());
        } else {
            tin.clear();
            inID.clear();
        }
    }
    void writeAndUpdateOut(size i, std::ofstream &raster_file) {
        raster_file.write((char*)&i, sizeof(i));
        if (i!=0) 
            raster_file.write((char*)&(tsp[0]), i*sizeof(tsp[0]));
        nOut = nOut + i; 
        size keeping = tsp.size()-i;
        if (keeping) 
            tsp.assign(tsp.end()-keeping, tsp.end());
        else tsp.clear(); 
    }
    void clear() {
        tin.clear();
        inID.clear();
        status = true;
        tPoi.clear();
        ei.clear();
        tmp.clear();
        poiGen.clear();
        ranGen.clear();
    }
} nNS;
#endif
