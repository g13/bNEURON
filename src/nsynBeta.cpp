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
#include "nNeuroSt2.h"
#include "nsynBiTest3.h"
#include "typedefs.h"
#include "neuroTest.h"
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
    ifstream cfg_file;
    string lib_file, para_file;
    mxArray *para ;
    double cpu_t_sim, cpu_t_bilinear, cpu_t_linear, cpu_t_bilinear0;
    double cpu_t_jbilinear, cpu_t_jlinear;
    double rLinear;
    //if (!win) {
        clockid_t clk_id = CLOCK_PROCESS_CPUTIME_ID;
        struct timespec tpS, tpE;
    //}
    unsigned int i,j,k,nt;
    double run_t, ignore_t, tau_ed, tau_er, tau_id, tau_ir;
    double gNa, vNa, gK, vK, gLeak, vLeak, vT, vE, vI, vRest, DeltaT,S;
    double tref, trans, trans0;
    double tmpsum;
    double rd;
    double sampleStep;
    double v_tol, vBuffer;
    unsigned int vinit;
    MATFile *matFile;
    nNL neuroLib;
    nNS neuron;
    po::variables_map vm;
    vector<double> simV;
    vector<double> biV, biV0;
    vector<double> liV;
    vector<double> jlv, jbv;
    vector<double> tsp_sim, tsp_bi, tsp_jbi, tsp_li, tsp_jli, tsp_bi0;
    //vector<double> vargList;
	unsigned int seed0 = static_cast<unsigned int>(std::time(NULL));
    unsigned int seed;
    string config_file;
    string theme;
    bool testing;
    vector<double> test0,test1,test2;

	po::options_description cmd_line_options("cmd line options"),
        config_file_options("config file"),
        Generic("options");
    cmd_line_options.add_options()
        ("config_file,c", po::value<string>(&config_file)->default_value("sR_init.cfg"),"config filename");
	Generic.add_options()
		("lib_file,L",po::value<string>(&lib_file),"neuron library")
		("para_file,p",po::value<string>(&para_file),"parameter file")
        ("seed,s",po::value<unsigned int>(&seed),"seeding")
		("theme,m",po::value<string>(&theme),"parameter file")
		("rate,r", po::value<vector<double>>()->multitoken()->composing(), "poisson rate array")
		("vThres", po::value<double>(), "NEURON vThres")
		("vRest", po::value<double>(), "NEURON vRest")
        ("sampleStep", po::value<double>(&sampleStep), "sampling step for jump data")
		("run_t,t", po::value<double>(&run_t), "sim time")
        ("vinit,v", po::value<unsigned int>(&vinit), "index of intial voltage in vRange")
        ("tref",po::value<double>(&tref),"refractory period")
        ("trans",po::value<double>(&trans),"transient VClamp time")
        ("trans0",po::value<double>(&trans0),"default transient VClamp time")
        ("rLinear",po::value<double>(&rLinear),"linear->HH threshold")
        ("v_tol",po::value<double>(&v_tol),"crossing threshold tolerance")
        ("vBuffer",po::value<double>(&vBuffer),"return threshold buffer")
        ("dendClampRatio",po::value<double>(&rd),"dendrite dV contribution ratio")
        ("test0",po::value<vector<double>>()->multitoken()->composing(),"testSynapse0")
        ("test1",po::value<vector<double>>()->multitoken()->composing(),"testSynapse1")
        ("test2",po::value<vector<double>>()->multitoken()->composing(),"testSynapse2")
        ("testing",po::value<bool>(&testing),"test")
		("ignore_t", po::value<double>(&ignore_t),"ingore time while applying bilinear rules");

	cmd_line_options.add(Generic);
	config_file_options.add(Generic);
	store(po::parse_command_line(argc, argv, cmd_line_options), vm);
	notify(vm);
	cfg_file.open(config_file);
	if (cfg_file) {
		store(po::parse_config_file(cfg_file, config_file_options, true), vm);
		notify(vm);
	} else {
		cout << "cannot open the configuration file: " << config_file << endl;
		return 0;
	}
	cfg_file.close();

    neuroLib.readLib(lib_file.c_str());
    
    vT = vm["vThres"].as<double>();
    vRest = vm["vRest"].as<double>();

    vector<double> r = vm["rate"].as<vector<double>>();
    if (!seed) seed = seed0;
    string prefix = "E"+to_string(static_cast<int>(r[0]))+"-I"+to_string(static_cast<int>(r[1]))+"-corrL"+to_string(static_cast<int>(neuroLib.tol_tb-ignore_t))+"-"+to_string(static_cast<int>(run_t))+"-"+to_string(static_cast<unsigned int>(seed));
    theme = "-" + theme;
	raster_file.open(prefix + "Raster" + theme + ".bin", ios::out|ios::binary);
	tIncome_file.open(prefix + "tIn" + theme + ".bin", ios::out|ios::binary);
	data_file.open(prefix + "Data" + theme + ".bin", ios::out|ios::binary);
	jND_file.open(prefix + "jND" + theme + ".bin", ios::out|ios::binary);

    double tstep = neuroLib.tstep;
    cout << "using tstep: " << tstep << " ms" << endl;
    cout << "E rate " << r[0] << ", I rate " << r[1] << endl;
    Py_Initialize();
    cout << "Python " << endl;
    _import_array();
    PyRun_SimpleString("import sys");
    PyRun_SimpleString("sys.path.append(\".\")");
    cout << "initialized " << endl;
    //Py_Finalize();
    //cout << "finalized " << endl;
    //getCell(neuroLib.loc, neuroLib.pos, neuroLib.gList, neuroLib.nSyn, Py_Cell, Py_synList, Py_vecStimList, pModule, vRest);
    synSet syn = {
        neuroLib.loc, 
        neuroLib.pos, 
        neuroLib.gList, 
        neuroLib.nSyn, 
        vRest
    }; 
    Cell cell(syn);

    vector<double> t0(neuroLib.nSyn,0);
    vector<double> max_rate,rate;
    for (i=0; i<r.size(); i++) {
        r[i] = r[i]/1000;
    }
    vector<double> rm(r.begin(),r.end());
    assert(r.size()==2 && r[0]>=0.0 && r[1] >= 0.0);
    for (i=0; i<neuroLib.nSyn; i++) {
        if (neuroLib.ei[i]) {
            max_rate.push_back(rm[0]);
            rate.push_back(r[0]);
        } else {
            max_rate.push_back(rm[1]);
            rate.push_back(r[1]);
        }
    }
    neuron.initialize(seed,neuroLib.nSyn,t0,run_t,rate,max_rate,neuroLib.ei,testing,tstep);
    neuron.trans = trans;
    neuron.tref = tref;
    neuron.v_tol = v_tol;
    cout << " initialized" << endl;
    // get external inputs

    if (!testing) {
        cout << "generating inputs" << endl;
        while (neuron.status) neuron.getNextInput(rate,max_rate);
        cout << "done" << endl;
    } else {
        test0 = vm["test0"].as<vector<double>>();
        test1 = vm["test1"].as<vector<double>>();
        test2 = vm["test2"].as<vector<double>>();
        cout << "Test Inputs" << endl;
        cout << "0: ";
        for (i=0;i<test0.size();i++) {
            test0[i] = static_cast<int>(test0[i]/tstep)*tstep;
            cout << test0[i] << " ";
        }
        cout << endl;
        cout << "1: ";
        for (i=0;i<test1.size();i++) {
            test1[i] = static_cast<int>(test1[i]/tstep)*tstep;
            cout << test1[i] << " ";
        }
        cout << endl;
        cout << "2: ";
        for (i=0;i<test2.size();i++) {
            test2[i] = static_cast<int>(test2[i]/tstep)*tstep;
            cout << test2[i] << " ";
        }
        cout << endl;
        vector<vector<double>> testInputs;
        testInputs.push_back(test0);
        testInputs.push_back(test1);
        testInputs.push_back(test2);
        neuron.setInputs(testInputs);
    }
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
    assert(vBuffer < 0);
    cout << " spikeTrain generated" << endl;

    nt = static_cast<unsigned int>(run_t/tstep)+1;
    double v0 = neuroLib.vRange[vinit];
    vector<vector<double>> dendV(neuron.nSyn,vector<double>());
    for (i=0; i<neuron.nSyn; i++) {
        dendV[i].reserve(nt);
    }
    simV.reserve(nt);
    //simV.push_back(v0);
    biV.reserve(nt);
    biV.push_back(v0);
    biV0.reserve(nt);
    biV0.push_back(v0);
    liV.reserve(nt);
    liV.push_back(v0);
    clock_gettime(clk_id,&tpS);
    neuron.vReset = vRest;
    neuron.vRest = vRest;
    double rHH = 1.3;
    double vCrossl = vRest + (vT -vRest)*rLinear;
    double vBackl = vRest + (vT -vRest)*rLinear + vBuffer;
    double vCrossb = vCrossl;
    double vBackb = vBackl;
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
    neuron.vThres = vRest + (vT -vRest)*rHH;
    cout << " yale NEURON begin" << endl;
    cout << " point of no return unless spike " << neuron.vThres << endl;
    vector<double> dendVclamp(neuroLib.nSyn,1000); // 1000 default no dend clamp
    //for (int irep=0; irep < 1; irep++) {
    //    cout << "Neuron rep: " << irep << endl;
    //    simV.clear();
    //    simV.reserve(nt);
          nc = Py_proceed(cell, v0, RList, s1,  spikeTrain, neuroLib.nSyn, trans0, trans0+run_t, plchldr_double, tref, neuron.vThres, 1, simV, plchldr_size0, tsp_sim, 0, tstep, dendVclamp, 1, dendV);
    //}
    clock_gettime(clk_id,&tpE);
    
    cpu_t_sim = static_cast<double>(tpE.tv_sec-tpS.tv_sec) + static_cast<double>(tpE.tv_nsec - tpS.tv_nsec)/1e9;
    cout << "sim ended, took " << cpu_t_sim << "s" << endl;
    cout << "spikes: " << nc << endl;
    nc = 0;

    cout << endl;
 
    cout << " bilinear begin" << endl;
    clock_gettime(clk_id,&tpS);
          nc = bilinear_nSyn(biV, neuroLib, neuron, run_t, ignore_t, tsp_bi, vCrossb, vBackb, neuron.ei);
    clock_gettime(clk_id,&tpE);
    cpu_t_bilinear = static_cast<double>(tpE.tv_sec-tpS.tv_sec) + static_cast<double>(tpE.tv_nsec - tpS.tv_nsec)/1e9;
    cout << "bilinear est. ended, took " << cpu_t_bilinear << "s" << endl;
    cout << "spikes: " << nc << endl;
    nc = 0;
    
    cout << endl;
    clock_gettime(clk_id,&tpS);
    double totalRate = 0;
    for (i=0;i<r.size();i++) {
        totalRate += r[i];
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

    nc = linear_nSyn(liV, neuroLib, neuron, run_t, ignore_t, tsp_li, vCrossl, vBackl, neuron.ei);

    clock_gettime(clk_id,&tpE);
    cpu_t_linear = static_cast<double>(tpE.tv_sec-tpS.tv_sec) + static_cast<double>(tpE.tv_nsec - tpS.tv_nsec)/1e9;
    cout << "linear est. ended, took " << cpu_t_linear << "s" << endl;
    cout << "spikes: " << nc << endl;

    cout << endl;
    cout << " bilinear0 begin: " << endl;
    clock_gettime(clk_id,&tpS);
        nc = bilinear0_nSyn(biV0, neuroLib, neuron, run_t, ignore_t, tsp_bi0, vCrossb, vBackb, neuron.ei);
    clock_gettime(clk_id,&tpE);
    cpu_t_bilinear0 = static_cast<double>(tpE.tv_sec-tpS.tv_sec) + static_cast<double>(tpE.tv_nsec - tpS.tv_nsec)/1e9;
    cout << "bilinear0 est. ended, took " << cpu_t_bilinear0 << "s" << endl;
    cout << "spikes: " << nc << endl;
    cout << endl;
    nc = 0;
    if (1) {
    vector<double>* output[4] = {&simV,&biV,&liV,&biV0};
    neuron.writeAndUpdateIn(neuron.tin.size(), tIncome_file);
    //neuron.writeAndUpdateOut(neuron.tsp.size(), raster_file);
    for (i=0;i<4;i++) {
        data_file.write((char*)output[i]->data(), nt*sizeof(double));
    }
    for (i=0;i<neuron.nSyn;i++) {
        data_file.write((char*)dendV[i].data(), nt*sizeof(double));
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
    for (size i=0;i<ncross;i++) {
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
    }
    if (jND_file.is_open())          jND_file.close();
    if (data_file.is_open())        data_file.close();
    if (raster_file.is_open())      raster_file.close();
    if (tIncome_file.is_open())     tIncome_file.close();
    if (1) {
    neuroLib.clearLib();
    NEURON_cleanup(cell);
    cout << " cell cleaned " << endl;
    Py_Finalize();
    }
    cout << "finalized " << endl;
    cout << "size: " <<  sizeof(size) << " bytes" << endl;
    cout << "size_b: "<< sizeof(size_b) << " bytes" << endl;
    cout << "time table: sim,   bi, jbi,    li, jli,    bi0" << endl;
    cout << "       " << cpu_t_sim << ", ";
    cout << cpu_t_bilinear << ", " << cpu_t_jbilinear << ", ";
    cout << cpu_t_linear << ", " << cpu_t_jlinear << ", ";
    cout << cpu_t_bilinear0;
    cout << endl;
}
