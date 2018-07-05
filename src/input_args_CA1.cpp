#include "input_args_CA1.h"

input_args_CA1::input_args_CA1() {
    po::options_description CA1_options;
    CA1_options.add_options()
		("libFile,l",po::value<string>(&libFile),"neuron library")
        ("vInit,v", po::value<unsigned int>(&vInit), "index of intial voltage in vRange")
        ("vThres", po::value<double>(&vThres), "spiking threshold")
        ("vRest", po::value<double>(&vRest), "resting threshold")
        ("tRef",po::value<double>(&tRef),"refractory period")
        ("trans",po::value<double>(&trans),"transient VClamp time")
        ("dtrans",po::value<double>(&dtrans),"dendClamp time")
        ("dtau",po::value<double>(&dtau),"dendv exp collecting time")
        ("trans0",po::value<double>(&trans0),"default transient VClamp time")
        ("rLinear",po::value<double>(&rLinear),"linear model spiking threshold")
        ("rBiLinear", po::value<double>(&rBiLinear), "bilinear model spiking threshold")
        ("vTol",po::value<double>(&vTol),"crossing threshold tolerance")
        ("vBuffer",po::value<double>(&vBuffer),"return threshold buffer")
        ("dendClampRatio",po::value<double>(&dendClampRatio)->default_value(1.0),"dendrite dV contribution ratio, negative value for dendV, positive value for constants over distance")
        ("afterCrossBehavior",po::value<int>(&afterCrossBehavior)->default_value(0), "0:no linear or bilinear extension after spike. 1:no bilinear extension. 2: all extend.")
        ("kVStyle",po::value<int>(&kVStyle)->default_value(0), "bilinear0 after spike kV style, 0: kV0 style. 1: kV style.")
        ("dtSquare",po::value<bool>(&dtSquare)->default_value(false)," if true, the second dt variable of kV is considered a second tier variable when interpolating")
        ("method,d",po::value<int>(&i)->default_value(-1), "0-5 for sim,bi,jb,jl,li,b0 !!")
        ("getDendV",po::value<bool>(&getDendV)->default_value(false), "record dendV at syn input locations")
        ("tIn",po::value<bool>(&tIn)->default_value(false), "record and output poisson input time")
        ("setV",po::value<bool>(&setV)->default_value(false), "if true cross with uniform membrane potential = vThres (preTran)")
        ("spikeShape",po::value<bool>(&spikeShape)->default_value(true),"if false, crossing is spiking")
        ("sliceDebugPlot",po::value<bool>(&sliceDebugPlot)->default_value(false),"plot every cross")
        ("pas",po::value<bool>(&pas)->default_value(false),"if true, crossing is spiking and reset")
        ("clusterDend", po::value<vector<int>>()->multitoken()->composing(), "sequence of clustering ID")
        ("clusterSizeCDF", po::value<vector<int>>()->multitoken()->composing(), "CDF of cluster size matching with ID sequence")
        ("clusterClampRatio", po::value<vector<double>>()->multitoken()->composing(), "dendClampRatio in clusters")
		("ignoreT", po::value<double>(&ignoreT),"ingore time while applying bilinear rules");
        
    cmdLineOptions.add(CA1_options);
    configFileOptions.add(CA1_options);
}
void input_args_CA1::setbit() {
    mode.set(i);
}
void input_args_CA1::get_rd() {
    vector<int> dendSeq = vm["clusterDend"].as<vector<int>>();
    clusterSizeCDF = vm["clusterSizeCDF"].as<vector<int>>();
    clusterClampRatio = vm["clusterClampRatio"].as<vector<double>>();
    assert(clusterClampRatio.size() == clusterSizeCDF.size());
    assert(dendSeq.size() == nInput);
    clusterDend.assign(clusterSizeCDF.size(),vector<int>());
    int j = 0;
    int i = 0;
    while (j<nInput) {
        clusterDend[i].push_back(dendSeq[j]);
        if (j==clusterSizeCDF[i]-1 && j!=nInput-1) {
            i++;     
            clusterSizeCDF[i] = clusterSizeCDF[i-1] + clusterSizeCDF[i];
        }
        j++;
    }
    assert(clusterSizeCDF.back() == nInput);
    cout << "Clustering Dendrites: {";
    for (i=0; i<clusterClampRatio.size(); i++) {
        cout << "(";
        for (j=0; j<clusterDend[i].size(); j++) {
            cout << clusterDend[i][j];
            if (j<clusterDend[i].size()-1) {
                 cout << ", ";
            } else {
                cout << ")";
            }
        }
        if (i<clusterClampRatio.size()-1) {
            cout << ",";
        } else {
            cout << "}\n";
        }
    }
}
