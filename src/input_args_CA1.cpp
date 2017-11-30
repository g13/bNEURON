#include "input_args_CA1.h"

input_args_CA1::input_args_CA1(){
    configFileOptions.add_options()
		("libFile,l",po::value<string>(&libFile),"neuron library")
		("paraFile,p",po::value<string>(&paraFile),"parameter file")
        ("vInit,v", po::value<unsigned int>(&vInit), "index of intial voltage in vRange")
        ("tRef",po::value<double>(&tRef),"refractory period")
        ("trans",po::value<double>(&trans),"transient VClamp time")
        ("trans0",po::value<double>(&trans0),"default transient VClamp time")
        ("rLinear",po::value<double>(&rLinear),"linear model spiking threshold")
        ("rBiLinear", po::value<double>(&rBiLinear), "bilinear model spiking threshold")
        ("vTol",po::value<double>(&vTol),"crossing threshold tolerance")
        ("vBuffer",po::value<double>(&vBuffer),"return threshold buffer")
        ("dendClampRatio",po::value<double>(&dendClampRatio),"dendrite dV contribution ratio")
		("ignoreT", po::value<double>(&ignoreT),"ingore time while applying bilinear rules");
        
    cmdLineOptions.add(configFileOptions);
}
