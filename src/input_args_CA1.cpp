#include "input_args_CA1.h"

input_args_CA1::input_args_CA1() {
    po::options_description CA1_options;
    CA1_options.add_options()
		("libFile,l",po::value<string>(&(this->libFile)),"neuron library")
        ("vInit,v", po::value<unsigned int>(&(this->vInit)), "index of intial voltage in vRange")
        ("tRef",po::value<double>(&(this->tRef)),"refractory period")
        ("trans",po::value<double>(&(this->trans)),"transient VClamp time")
        ("trans0",po::value<double>(&(this->trans0)),"default transient VClamp time")
        ("rLinear",po::value<double>(&(this->rLinear)),"linear model spiking threshold")
        ("rBiLinear", po::value<double>(&(this->rBiLinear)), "bilinear model spiking threshold")
        ("vTol",po::value<double>(&(this->vTol)),"crossing threshold tolerance")
        ("vBuffer",po::value<double>(&(this->vBuffer)),"return threshold buffer")
        ("dendClampRatio",po::value<double>(&(this->dendClampRatio)),"dendrite dV contribution ratio")
		("ignoreT", po::value<double>(&(this->ignoreT)),"ingore time while applying bilinear rules");
        
    this->cmdLineOptions.add(CA1_options);
    this->configFileOptions.add(CA1_options);
}
