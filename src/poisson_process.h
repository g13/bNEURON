// generate poisson process
//
#include<cmath>
#include<random> 
#include<vector>
#include<cstring>
using std::string;
using std::vector;

double get_next_spike_time(double lastSpikeTime, double currentPoissonRate, double maxPoissonRate);

void predetermine_poisson_spikes(vector<double> &spikeTime, vector<double> &inputRate, double maxPoissonRate, string outputFn="", bool saveInputToFile = false);
