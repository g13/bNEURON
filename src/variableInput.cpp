#include "variableInput.h"

namespace po = boost::program_options;

int main(int argc, char **argv)
{
    //bool win = true;
    //feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
    ofstream data_file, tIncome_file, raster_file, jND_file, cpu_file;
    double cpu_t_sim, cpu_t_bilinear, cpu_t_linear, cpu_t_bilinear0;
    double cpu_t_jbilinear, cpu_t_jlinear;
    clockid_t clk_id = CLOCK_PROCESS_CPUTIME_ID;
    struct timespec tpS, tpE;
    vector<double> simV;
    vector<double> biV, biV0;
    vector<double> liV;
    vector<double> tsp_sim, tsp_bi, tsp_jbi, tsp_li, tsp_jli, tsp_bi0;

    InputArgsCA1 inputArg = input_args_CA1();
    inputArg.read(argc,argv);
    if (inputArg.dtVarLevels) {
        cout << "different dt can only affect Yale NEURON simulation" << endl;
    }
    nNL neuroLib(inputArg.libFile.c_str());
    if (neuroLib.nSyn >= inputArg.nInput) {
        neuroLib.nSyn = inputArg.nInput;
    } else {
        cout << " nInput > nSyn, exit" << endl;
        return 0;
    }
    nNS neuron(inputArg.seed,neuroLib.nSyn,neuroLib.ei,inputArg.trans,inputArg.tRef,inputArg.vTol);
    inputArg.reformat_input_table(neuroLib.tstep);

    vector<vector<double>> rate(inputArg.inputLevel.size(),vector<double>());
    cout << "exact rate levels " << endl;
    for (int i=0; i<inputArg.inputLevel.size(); i++) {
        rate[i].reserve(neuroLib.nSyn);
        if (inputArg.input[i].size() != neuroLib.nSyn) {
            cout << inputArg.input[i].size() << "!=" << neuroLib.nSyn << endl;
            assert(inputArg.input[i].size() == neuroLib.nSyn);
        }
        for (int j=0; j<neuroLib.nSyn; j++) {
            rate[i].push_back(inputArg.inputLevel[i]*inputArg.input[i][j]/1000.0);
            cout << rate[i][j] << ", ";
        }
        cout << endl;
    }
    string prefix = inputArg.theme + "-s" + to_string(static_cast<unsigned int>(inputArg.seed));
	data_file.open("Data-" + prefix + ".bin", ios::out|ios::binary);
	raster_file.open("Raster-" + prefix + ".bin", ios::out|ios::binary);
	tIncome_file.open("tIn-" + prefix + ".bin", ios::out|ios::binary);
	jND_file.open("jND-" + prefix + ".bin", ios::out|ios::binary);
	cpu_file.open("cpuTime-" + prefix + ".bin", ios::out|ios::binary);

    double tstep = neuroLib.tstep;
    cout << "using tstep: " << tstep << " ms" << endl;
    Py_Initialize();
    cout << "Python " << endl;
    PyRun_SimpleString("import sys");
    PyRun_SimpleString("sys.path.append(\".\")");

    SynSet syn = {
        neuroLib.loc, 
        neuroLib.pos, 
        neuroLib.gList, 
        neuroLib.nSyn, 
        inputArg.vRest
    }; 
    cout << " NEURON" << endl;
    Cell cell(syn);

    cout << "initialized " << endl;
    for (int ii=0; ii<inputArg.inputLevel.size(); ii++) {
        vector<double> t0(neuroLib.nSyn,0);
        double run_t = inputArg.runTime[ii];
        cout << "generating inputs" << endl;
        neuron.initialize(inputArg.runTime[ii],tstep,t0,rate[ii]);
        if (inputArg.vInit >= neuroLib.nv) {
            cout << "vInit " << inputArg.vInit << " out of range, nv " << neuroLib.nv << endl;
            assert(inputArg.vInit < neuroLib.nv);
        }
        double v0 = neuroLib.vRange[inputArg.vInit];
        cout << "start v: " << v0 << endl;
        while (neuron.status) neuron.getNextInput(rate[ii]);
        cout << "done" << endl;
        cout << endl;
        vector<vector<double>> spikeTrain;
        cout << "INPUTS: ";
        spikeTrain.assign(neuroLib.nSyn,vector<double>());
        cout << "(";
        for (int i=0; i<neuroLib.nSyn; i++) {
            cout << "[";
            for (int j=0; j<neuron.tPoi[i].size(); j++) {
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
        for (int i=0; i<neuroLib.nSyn; i++) {
            cout << spikeTrain[i].size(); 
            if (i<neuroLib.nSyn-1) {
                cout << ", ";
            }
        }
        cout << endl;
        assert(inputArg.vBuffer < 0);
        cout << " spikeTrain generated" << endl;

        unsigned int nt = static_cast<unsigned int>(run_t/tstep)+1;
        unsigned int nt0 = static_cast<unsigned int>(run_t/inputArg.tstep[ii])+1;
        vector<vector<double>> dendV(neuron.nSyn,vector<double>());
        for (int i=0; i<neuron.nSyn; i++) {
            dendV[i].reserve(nt0);
        }
        simV.reserve(nt0);
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
        unsigned int nc = 0;
        // NEURON
        neuron.vThres = inputArg.vThres;
        cout << " yale NEURON begin" << endl;
        cout << " point of no return unless spike " << inputArg.vThres << endl;
        vector<double> dendVclamp(neuroLib.nSyn,1000); // 1000 default no dend clamp
        nc = Py_proceed(cell, v0, RList, s1,  spikeTrain, neuroLib.nSyn, inputArg.trans0, inputArg.trans0 + run_t, plchldr_double, inputArg.tRef, inputArg.vThres, 1, simV, plchldr_size0, tsp_sim, 0, inputArg.tstep[ii], dendVclamp, -1, 1, dendV);
        clock_gettime(clk_id,&tpE);
        
        cpu_t_sim = static_cast<double>(tpE.tv_sec-tpS.tv_sec) + static_cast<double>(tpE.tv_nsec - tpS.tv_nsec)/1e9;
        cout << "sim ended, took " << cpu_t_sim << "s" << endl;
        cout << "spikes: " << nc << endl;
        nc = 0;

        cout << endl;

        cout << " bilinear begin" << endl;
        clock_gettime(clk_id,&tpS);

        nc = bilinear_nSyn(cell, spikeTrain, dendVclamp, biV, neuroLib, neuron, run_t, inputArg.ignoreT, tsp_bi, vCrossb, vBackb, inputArg.afterCrossBehavior, inputArg.spikeShape, inputArg.dtSquare);
        clock_gettime(clk_id,&tpE);
        cpu_t_bilinear = static_cast<double>(tpE.tv_sec-tpS.tv_sec) + static_cast<double>(tpE.tv_nsec - tpS.tv_nsec)/1e9;
        cout << "bilinear est. ended, took " << cpu_t_bilinear << "s" << endl;
        cout << "spikes: " << nc << endl;
        nc = 0;
        
        cout << endl;
        clock_gettime(clk_id,&tpS);
        double totalRate = 0;
        for (int i=0;i<rate[ii].size();i++) {
            totalRate += rate[ii][i];
        }
        size rSize = static_cast<size>((totalRate)*run_t);
        size corrSize = static_cast<size>((totalRate)*neuroLib.tol_tb);
        jND jndb(rSize);
        Input inputb(rSize);
        Cross crossb(nt,v0);
        if (ii>0) {
            jndb.initialize(rSize);
            inputb.initialize(rSize);
            crossb.initialize(nt,v0);
        }

        cout << " jBilinear begin" << endl;
        nc = nsyn_jBilinear(cell, spikeTrain, dendVclamp, neuron, neuroLib, inputb, jndb, crossb, run_t, inputArg.ignoreT, corrSize, tsp_jbi, vCrossb, vBackb, inputArg.afterCrossBehavior, false, inputArg.dtSquare);

        clock_gettime(clk_id,&tpE);
        cpu_t_jbilinear = static_cast<double>(tpE.tv_sec-tpS.tv_sec) + static_cast<double>(tpE.tv_nsec - tpS.tv_nsec)/1e9;
        cout << "jBilinear est. ended, took " << cpu_t_jbilinear << "s" << endl;
        cout << "spikes: " << nc << endl;
        cout << endl;

        clock_gettime(clk_id,&tpS);
        rSize = static_cast<size>((totalRate)*run_t);
        corrSize = static_cast<size>((totalRate)*neuroLib.nt);
        jND jndl(rSize);
        Input inputl(rSize);
        Cross crossl(nt,v0);
        if (ii>0) {
            jndl.initialize(rSize);
            inputl.initialize(rSize);
            crossl.initialize(nt,v0);
        }
        
        cout << " jLinear begin" << endl;
        nc = nsyn_jLinear(cell, spikeTrain, dendVclamp, neuron, neuroLib, inputl, jndl, crossl, run_t,inputArg.ignoreT, corrSize, tsp_jli, vCrossl, vBackl, inputArg.afterCrossBehavior, false); //inputArg.spikeShape

        clock_gettime(clk_id,&tpE);
        cpu_t_jlinear = static_cast<double>(tpE.tv_sec-tpS.tv_sec) + static_cast<double>(tpE.tv_nsec - tpS.tv_nsec)/1e9;
        cout << "jLinear est. ended, took " << cpu_t_jlinear << "s" << endl;
        cout << "spikes: " << nc << endl;
        
        cout << endl;
        cout << " linear begin: " << endl;
        clock_gettime(clk_id,&tpS);
        nc = 0;

        nc = linear_nSyn(cell, spikeTrain, dendVclamp, liV, neuroLib, neuron, run_t, inputArg.ignoreT, tsp_li, vCrossl, vBackl, inputArg.afterCrossBehavior, inputArg.spikeShape);

        clock_gettime(clk_id,&tpE);
        cpu_t_linear = static_cast<double>(tpE.tv_sec-tpS.tv_sec) + static_cast<double>(tpE.tv_nsec - tpS.tv_nsec)/1e9;
        cout << "linear est. ended, took " << cpu_t_linear << "s" << endl;
        cout << "spikes: " << nc << endl;

        cout << endl;
        cout << " bilinear0 begin: " << endl;
        clock_gettime(clk_id,&tpS);
            nc = bilinear0_nSyn(cell, spikeTrain, dendVclamp, biV0, neuroLib, neuron, run_t, inputArg.ignoreT, tsp_bi0, vCrossb, vBackb, inputArg.afterCrossBehavior, inputArg.spikeShape, inputArg.kVStyle, inputArg.dtSquare);
        clock_gettime(clk_id,&tpE);
        cpu_t_bilinear0 = static_cast<double>(tpE.tv_sec-tpS.tv_sec) + static_cast<double>(tpE.tv_nsec - tpS.tv_nsec)/1e9;
        cout << "bilinear0 est. ended, took " << cpu_t_bilinear0 << "s" << endl;
        cout << "spikes: " << nc << endl;
        cout << endl;

        cpu_file.write((char*)&(cpu_t_sim),sizeof(double));
        cpu_file.write((char*)&(cpu_t_bilinear), sizeof(double));
        cpu_file.write((char*)&(cpu_t_linear), sizeof(double));
        cpu_file.write((char*)&(cpu_t_jbilinear), sizeof(double));
        cpu_file.write((char*)&(cpu_t_jlinear), sizeof(double));
        cpu_file.write((char*)&(cpu_t_bilinear0), sizeof(double));

        neuron.writeAndUpdateIn(neuron.tin.size(), tIncome_file);
        vector<double>* output[4] = {&simV,&biV,&liV,&biV0};
        assert(simV.size()==nt0);
        assert(biV.size()==nt);
        assert(dendV[0].size()==nt0);
        write_data_file(data_file, output, dendV, nt0, nt, neuron.nSyn);
        
        vector<size> tmpSize;
        size ncross, jndSize;
        //neuron.writeAndUpdateOut(neuron.tsp.size(), raster_file);
        output[0] = &jndb.t;
        output[1] = &jndb.v;
        assert(jndb.v.size() == jndb.t.size());
        jndSize = jndb.t.size();
        size_data_write(jND_file, output, 2, jndSize, 0);

        tmpSize.reserve(crossb.nCross);
        for (int i=0;i<crossb.nCross;i++) {
            tmpSize.push_back(crossb.iCross[i+1] - crossb.iCross[i]);
            cout << " tmpSize " << i << " : " << tmpSize.back() << endl;
        }
        if (crossb.v.size() != crossb.t.size()) {
            cout << crossb.v.size() <<  "!= " << crossb.t.size() << endl;
            assert(crossb.v.size() == crossb.t.size());
        }
        output[0] = &(crossb.v);
        output[1] = &(crossb.t);
        nsection_write(jND_file, output, 2, crossb.iCross, tmpSize);
        tmpSize.clear();

        output[0] = &jndl.t;
        output[1] = &jndl.v;
        assert(jndl.v.size() == jndl.t.size());
        jndSize = jndl.t.size();
        size_data_write(jND_file, output, 2, jndSize, 0);

        tmpSize.reserve(crossb.nCross);
        for (int i=0;i<crossb.nCross;i++) {
            tmpSize.push_back(crossl.iCross[i+1] - crossl.iCross[i]);
        }
        if (crossl.v.size() != crossl.t.size()) {
            cout << crossl.v.size() <<  "!= " << crossl.t.size() << endl;
            assert(crossl.v.size() == crossl.t.size());
        }
        output[0] = &(crossl.v);
        output[1] = &(crossl.t);
        nsection_write(jND_file, output, 2, crossl.iCross, tmpSize);
        tmpSize.clear();

        size rasterSize = tsp_sim.size();
        output[0] = &tsp_sim;
        size_data_write(raster_file, output, 1, rasterSize, 0);

        output[0] = &tsp_bi;
        rasterSize = tsp_bi.size();
        size_data_write(raster_file, output, 1, rasterSize, 0);

        output[0] = &tsp_li;
        rasterSize = tsp_li.size();
        size_data_write(raster_file, output, 1, rasterSize, 0);

        output[0] = &tsp_jbi;
        rasterSize = tsp_jbi.size();
        size_data_write(raster_file, output, 1, rasterSize, 0);

        output[0] = &tsp_jli;
        rasterSize = tsp_jli.size();
        size_data_write(raster_file, output, 1, rasterSize, 0);

        output[0] = &tsp_bi0;
        rasterSize = tsp_bi0.size();
        size_data_write(raster_file, output, 1, rasterSize, 0);

        simV.clear();
        biV.clear();
        biV0.clear();
        liV.clear();

        tsp_sim.clear();
        tsp_bi.clear();
        tsp_li.clear();
        tsp_jbi.clear();
        tsp_jli.clear();
        tsp_bi0.clear();
        cout << "level " << ii << " finalized " << endl;
        cout << " run time " << inputArg.runTime[ii] << " tstep " << tstep << endl;
        cout << "time table: sim,   bi, jbi,    li, jli,    bi0" << endl;
        cout << "       " << cpu_t_sim << ", ";
        cout << cpu_t_bilinear << ", " << cpu_t_jbilinear << ", ";
        cout << cpu_t_linear << ", " << cpu_t_jlinear << ", ";
        cout << cpu_t_bilinear0;
        cout << endl;
        neuron.clear();
        spikeTrain.clear();
    }
    if (data_file.is_open())        data_file.close();
    if (jND_file.is_open())         jND_file.close();
    if (raster_file.is_open())      raster_file.close();
    if (tIncome_file.is_open())     tIncome_file.close();
    if (cpu_file.is_open())         cpu_file.close();
    neuroLib.clearLib();
    NEURON_cleanup(cell);
    cout << " cell cleaned " << endl;
    Py_Finalize();
    cout << "size: " <<  sizeof(size) << " bytes" << endl;
    cout << "size_b: "<< sizeof(size_b) << " bytes" << endl;
}

void write_data_file(ofstream &data_file, vector<double>** output, vector<vector<double>> &dendV, unsigned int nt0, unsigned int nt, int nSyn) {
    data_file.write((char*)output[0]->data(), nt0*sizeof(double));
    for (int i=1;i<4;i++) {
        data_file.write((char*)output[i]->data(), nt*sizeof(double));
    }
    for (int i=0;i<nSyn;i++) {
        data_file.write((char*)dendV[i].data(), nt0*sizeof(double));
    }
}
