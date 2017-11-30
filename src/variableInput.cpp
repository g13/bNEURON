//#include <fenv.h>;
#include "mex.h"
#include <boost/program_options.hpp>
#include <cstring>
//#include "boost_program_options_overload.h"
#include "matFunc.h"
#include "time.h"
#include "nsynJtest3.h"
#include "nsynJltest3.h"
#include "nNeuroLib.h"
#include "nNeuroSt.h"
#include "nsynBiTest3.h"
#include "typedefs.h"
#include "Yale_NEURON_PyAPI.h"
#include "input_args_CA1.h"
#include <fstream>
#include <iostream>

using std::cout;
using std::endl;
using std::vector;
using std::string;
using std::ofstream;
using std::ifstream;
using std::to_string;
using std::ios;
namespace po = boost::program_options;

int main(int argc, char **argv)
{
    //bool win = true;
    //feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
    // prhs[] ~ lib_file, para_file, ith, run_t, ignore_t, vinit
    ofstream tIncome_file, raster_file, data_file, jND_file;
    mxArray *para ;
    double cpu_t_sim, cpu_t_bilinear, cpu_t_linear, cpu_t_bilinear0;
    double cpu_t_jbilinear, cpu_t_jlinear;
    unsigned int i,j,k,nt, nt0;
    clockid_t clk_id = CLOCK_PROCESS_CPUTIME_ID;
    struct timespec tpS, tpE;
    MATFile *matFile;
    vector<double> simV;
    vector<double> biV, biV0;
    vector<double> liV;
    vector<double> jlv, jbv;
    vector<double> tsp_sim, tsp_bi, tsp_jbi, tsp_li, tsp_jli, tsp_bi0;

    //InputArgsCA1 inputArg = input_args_CA1();
    InputArgsCA1 inputArg = input_args_CA1();
    inputArg.read(argc,argv);
    if (inputArg.dtVarLevels) {
        cout << "different dt can only affect YALE NEURON simulation" << endl;
    }

    nNL neuroLib(inputArg.libFile.c_str());
    nNS neuron(inputArg.seed,neuroLib.nSyn,neuroLib.ei,inputArg.trans,inputArg.tRef,inputArg.vTol);

    vector<vector<double>> rate(neuroLib.nSyn,vector<double>());
    for (int i; i<inputArg.inputLevel.size(); i++) {
        rate[i].reserve(neuroLib.nSyn);
        for (int j; j<neuroLib.nSyn; j++) {
            rate[i].push_back(inputArg.inputLevel[i]*inputArg.input[0][i]);
        }
    }
    string prefix = inputArg.theme + "-s" + to_string(static_cast<unsigned int>(inputArg.seed)) + "-at" + to_string(static_cast<unsigned int>(std::time(NULL)));
	raster_file.open(prefix + "Raster.bin", ios::out|ios::binary);
	tIncome_file.open(prefix + "tIn.bin", ios::out|ios::binary);
	data_file.open(prefix + "Data.bin", ios::out|ios::binary);
	jND_file.open(prefix + "jND.bin", ios::out|ios::binary);

    double tstep = neuroLib.tstep;
    cout << "using tstep: " << tstep << " ms" << endl;
    Py_Initialize();
    cout << "Python " << endl;
    _import_array();
    PyRun_SimpleString("import sys");
    PyRun_SimpleString("sys.path.append(\".\")");
    cout << "initialized " << endl;

    SynSet syn = {
        neuroLib.loc, 
        neuroLib.pos, 
        neuroLib.gList, 
        neuroLib.nSyn, 
        inputArg.vRest
    }; 
    Cell cell(syn);

    for (int ii=0; ii<inputArg.inputLevel.size(); ii++) {
        vector<double> t0(neuroLib.nSyn,0);
        double run_t = inputArg.runTime[ii];
        neuron.initialize(inputArg.runTime[ii],tstep,t0,rate[ii]);
        double v0 = neuroLib.vRange[inputArg.vInit];
        while (neuron.status) neuron.getNextInput(rate[ii]);
        cout << " initialized" << endl;
        // get external inputs

        cout << "generating inputs" << endl;
        cout << "done" << endl;
        cout << endl;

        cout << "INPUTS: ";
        vector<vector<double>> spikeTrain(neuroLib.nSyn,vector<double>());
        cout << "(";
        for (i=0; i<neuroLib.nSyn; i++) {
            cout << "[";
            for (j=0; j<neuron.tPoi[i].size(); j++) {
                //spikeTrain[i].push_back(neuron.tPoi[i][j] + trans0);
                spikeTrain[i].push_back(neuron.tPoi[i][j]);
                cout << neuron.tPoi[i][j];
                if (j<neuron.tPoi[i].size()-1) {
                    cout << ", ";
                }
            }
            cout << "]";
            if (i<neuroLib.nSyn-1) {
                cout << ", ";
            } else {
                cout << ")" << endl;
            }
        }
        cout << "total spikes for synapse: ";
        for (i=0; i<neuroLib.nSyn; i++) {
            cout << spikeTrain[i].size(); 
            if (i<neuroLib.nSyn-1) {
                cout << ", ";
            }
        }
        cout << endl;
        assert(inputArg.vBuffer < 0);
        cout << " spikeTrain generated" << endl;

        nt = static_cast<unsigned int>(run_t/tstep)+1;
        nt0 = static_cast<unsigned int>(run_t/inputArg.tstep[ii])+1;
        vector<vector<double>> dendV(neuron.nSyn,vector<double>());
        for (i=0; i<neuron.nSyn; i++) {
            dendV[i].reserve(nt0);
        }
        simV.reserve(nt0);
        //simV.push_back(v0);
        biV.reserve(nt);
        biV.push_back(v0);
        biV0.reserve(nt);
        biV0.push_back(v0);
        liV.reserve(nt);
        liV.push_back(v0);
        clock_gettime(clk_id,&tpS);
        neuron.vReset = inputArg.vRest;
        neuron.vRest = inputArg.vRest;
        double vCrossl = inputArg.vRest + (inputArg.vThres - inputArg.vRest)*inputArg.rLinear;
        double vBackl = inputArg.vRest + (inputArg.vThres - inputArg.vRest)*inputArg.rLinear + inputArg.vBuffer;
        double vCrossb = inputArg.vRest + (inputArg.vThres - inputArg.vRest)*inputArg.rBiLinear;
        double vBackb = inputArg.vRest + (inputArg.vThres - inputArg.vRest)*inputArg.rBiLinear + inputArg.vBuffer;
        cout << " linear -> HH  " << vCrossl << endl;
        cout << " HH -> linear " << vBackl << endl;
        cout << " bilinear -> HH  " << vCrossb << endl;
        cout << " HH -> bilinear " << vBackb << endl;
        size plchldr_size0,plchldr_size1=0;
        double plchldr_double = 0;
        vector<double> plchldr_vec;
        vector<long> s1(neuroLib.nSyn,0);
        vector<vector<double>> RList(neuroLib.nSyn,vector<double>(2,0));
        //for (int ir=0; ir<neuroLib.nSyn; ir++) {
        //    RList[ir][0] = 1.89e-13*(ir*1.0/(neuroLib.nSyn*1.0));
        //    RList[ir][1] = 2.89e-53*(ir*1.0/(neuroLib.nSyn*1.0));
        //}
        unsigned int nc = 0;
        // NEURON
        neuron.vThres = inputArg.vThres;
        cout << " yale NEURON begin" << endl;
        cout << " point of no return unless spike " << inputArg.vThres << endl;
        vector<double> dendVclamp(neuroLib.nSyn,1000); // 1000 default no dend clamp
        nc = Py_proceed(cell, v0, RList, s1,  spikeTrain, neuroLib.nSyn, inputArg.trans0, inputArg.trans0 + run_t, plchldr_double, inputArg.tRef, inputArg.vThres, 1, simV, plchldr_size0, tsp_sim, 0, inputArg.tstep[ii], dendVclamp, 1, dendV);
        //}
        clock_gettime(clk_id,&tpE);
        
        cpu_t_sim = static_cast<double>(tpE.tv_sec-tpS.tv_sec) + static_cast<double>(tpE.tv_nsec - tpS.tv_nsec)/1e9;
        cout << "sim ended, took " << cpu_t_sim << "s" << endl;
        cout << "spikes: " << nc << endl;
        nc = 0;

        cout << endl;

        cout << " bilinear begin" << endl;
        clock_gettime(clk_id,&tpS);

        nc = bilinear_nSyn(biV, neuroLib, neuron, run_t, inputArg.ignoreT, tsp_bi, vCrossb, vBackb, neuron.ei);
        clock_gettime(clk_id,&tpE);
        cpu_t_bilinear = static_cast<double>(tpE.tv_sec-tpS.tv_sec) + static_cast<double>(tpE.tv_nsec - tpS.tv_nsec)/1e9;
        cout << "bilinear est. ended, took " << cpu_t_bilinear << "s" << endl;
        cout << "spikes: " << nc << endl;
        nc = 0;
        
        cout << endl;
        clock_gettime(clk_id,&tpS);
        double totalRate = 0;
        for (i=0;i<rate[ii].size();i++) {
            totalRate += rate[ii][i];
        }
        size rSize = static_cast<size>((totalRate)*run_t);
        size corrSize = static_cast<size>((totalRate)*neuroLib.tol_tb);
        jND jndb(rSize);
        Input inputb(rSize);
        Cross crossb(nt,v0);
        neuron.vThres = vCrossb;

        cout << " jBilinear begin" << endl;
        nc = nsyn_jBilinear(neuron, neuroLib, inputb, jndb, crossb, run_t, jbv, corrSize, tsp_jbi, vBackb, neuron.ei);

        clock_gettime(clk_id,&tpE);
        cpu_t_jbilinear = static_cast<double>(tpE.tv_sec-tpS.tv_sec) + static_cast<double>(tpE.tv_nsec - tpS.tv_nsec)/1e9;
        cout << "jBilinear est. ended, took " << cpu_t_jbilinear << "s" << endl;
        cout << "spikes: " << nc << endl;
        cout << endl;

        cout << " jLinear begin" << endl;
        clock_gettime(clk_id,&tpS);
        rSize = static_cast<size>((totalRate)*run_t);
        corrSize = static_cast<size>((totalRate)*neuroLib.nt);
        jND jndl(rSize);
        Input inputl(rSize);
        Cross crossl(nt,v0);
        neuron.vThres = vCrossl;
        
        nc = nsyn_jLinear(neuron, neuroLib, inputl, jndl, crossl, run_t, jlv, corrSize, tsp_jli, vBackl, neuron.ei);

        clock_gettime(clk_id,&tpE);
        cpu_t_jlinear = static_cast<double>(tpE.tv_sec-tpS.tv_sec) + static_cast<double>(tpE.tv_nsec - tpS.tv_nsec)/1e9;
        cout << "jLinear est. ended, took " << cpu_t_jlinear << "s" << endl;
        cout << "spikes: " << nc << endl;
        
        cout << endl;
        cout << " linear begin: " << endl;
        clock_gettime(clk_id,&tpS);
        nc = 0;

        nc = linear_nSyn(liV, neuroLib, neuron, run_t, inputArg.ignoreT, tsp_li, vCrossl, vBackl, neuron.ei);

        clock_gettime(clk_id,&tpE);
        cpu_t_linear = static_cast<double>(tpE.tv_sec-tpS.tv_sec) + static_cast<double>(tpE.tv_nsec - tpS.tv_nsec)/1e9;
        cout << "linear est. ended, took " << cpu_t_linear << "s" << endl;
        cout << "spikes: " << nc << endl;

        cout << endl;
        cout << " bilinear0 begin: " << endl;
        clock_gettime(clk_id,&tpS);
            nc = bilinear0_nSyn(biV0, neuroLib, neuron, run_t, inputArg.ignoreT, tsp_bi0, vCrossb, vBackb, neuron.ei);
        clock_gettime(clk_id,&tpE);
        cpu_t_bilinear0 = static_cast<double>(tpE.tv_sec-tpS.tv_sec) + static_cast<double>(tpE.tv_nsec - tpS.tv_nsec)/1e9;
        cout << "bilinear0 est. ended, took " << cpu_t_bilinear0 << "s" << endl;
        cout << "spikes: " << nc << endl;
        cout << endl;
        nc = 0;
        vector<double>* output[4] = {&simV,&biV,&liV,&biV0};
        neuron.writeAndUpdateIn(neuron.tin.size(), tIncome_file);
        //neuron.writeAndUpdateOut(neuron.tsp.size(), raster_file);
        data_file.write((char*)output[0]->data(), nt0*sizeof(double));
        for (i=1;i<4;i++) {
            data_file.write((char*)output[i]->data(), nt*sizeof(double));
        }
        for (i=0;i<neuron.nSyn;i++) {
            data_file.write((char*)dendV[i].data(), nt0*sizeof(double));
        }
        size jndSize = jndb.t.size();
        jND_file.write((char*)&(jndSize),sizeof(size));
        jND_file.write((char*)&(jndb.t[0]),jndSize*sizeof(double));
        assert(jndb.v.size() == jndb.t.size());
        jND_file.write((char*)&(jndb.v[0]),jndSize*sizeof(double));
        size ncross = crossb.iCross.size()-1;
        if (crossb.v.size() != crossb.t.size()) {
            cout << crossb.v.size() <<  "!= " << crossb.t.size() << endl;
            assert(crossb.v.size() == crossb.t.size());
        }

        jND_file.write((char*)&(ncross),sizeof(size));
        size tmpSize;
        for (i=0;i<ncross;i++) {
            tmpSize = crossb.iCross[i+1] - crossb.iCross[i];
            jND_file.write((char*)&(tmpSize),sizeof(size));
            jND_file.write((char*)&(crossb.t[crossb.iCross[i]]),tmpSize*sizeof(double));
            jND_file.write((char*)&(crossb.v[crossb.iCross[i]]),tmpSize*sizeof(double));
        }

        jndSize = jndl.t.size();
        jND_file.write((char*)&(jndSize),sizeof(size));
        jND_file.write((char*)&(jndl.t[0]),jndSize*sizeof(double));
        assert(jndl.v.size() == jndl.t.size());
        jND_file.write((char*)&(jndl.v[0]),jndSize*sizeof(double));
        ncross = crossl.iCross.size()-1;
        if (crossl.v.size() != crossl.t.size()) {
            cout << crossl.v.size() <<  "!= " << crossl.t.size() << endl;
            assert(crossl.v.size() == crossl.t.size());
        }
        jND_file.write((char*)&(ncross),sizeof(size));
        tmpSize;
        for ( i=0;i<ncross;i++) {
            tmpSize = crossl.iCross[i+1] - crossl.iCross[i];
            jND_file.write((char*)&(tmpSize),sizeof(size));
            jND_file.write((char*)&(crossl.t[crossl.iCross[i]]),tmpSize*sizeof(double));
            jND_file.write((char*)&(crossl.v[crossl.iCross[i]]),tmpSize*sizeof(double));
        }

        size rasterSize = tsp_sim.size();
        raster_file.write((char*)&rasterSize,sizeof(size));
        raster_file.write((char*)&(tsp_sim[0]),rasterSize*sizeof(double));
        rasterSize = tsp_bi.size();
        raster_file.write((char*)&rasterSize,sizeof(size));
        raster_file.write((char*)&(tsp_bi[0]),rasterSize*sizeof(double));
        rasterSize = tsp_li.size();
        raster_file.write((char*)&rasterSize,sizeof(size));
        raster_file.write((char*)&(tsp_li[0]),rasterSize*sizeof(double));
        rasterSize = tsp_jbi.size();
        raster_file.write((char*)&rasterSize,sizeof(size));
        raster_file.write((char*)&(tsp_jbi[0]),rasterSize*sizeof(double));
        rasterSize = tsp_jli.size();
        raster_file.write((char*)&rasterSize,sizeof(size));
        raster_file.write((char*)&(tsp_jli[0]),rasterSize*sizeof(double));
        rasterSize = tsp_bi0.size();
        raster_file.write((char*)&rasterSize,sizeof(size));
        raster_file.write((char*)&(tsp_bi0[0]),rasterSize*sizeof(double));
        simV.clear();
        biV.clear();
        biV0.clear();
        jlv.clear();
        jbv.clear();
        liV.clear();
        cout << "level " << ii << " finalized " << endl;
        cout << " run time " << inputArg.runTime[ii] << " tstep " << tstep << endl;
        cout << "time table: sim,   bi, jbi,    li, jli,    bi0" << endl;
        cout << "       " << cpu_t_sim << ", ";
        cout << cpu_t_bilinear << ", " << cpu_t_jbilinear << ", ";
        cout << cpu_t_linear << ", " << cpu_t_jlinear << ", ";
        cout << cpu_t_bilinear0;
        cout << endl;
        neuron.clear();
    }
    if (jND_file.is_open())          jND_file.close();
    if (data_file.is_open())        data_file.close();
    if (raster_file.is_open())      raster_file.close();
    if (tIncome_file.is_open())     tIncome_file.close();
    neuroLib.clearLib();
    NEURON_cleanup(cell);
    cout << " cell cleaned " << endl;
    Py_Finalize();
    cout << "size: " <<  sizeof(size) << " bytes" << endl;
    cout << "size_b: "<< sizeof(size_b) << " bytes" << endl;
}
