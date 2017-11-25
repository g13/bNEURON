#include "poisson_process.h"

double get_next_spike_time(vector<std::minstd_rand> &poiGen, double lastSpikeTime, double currentPoissonRate, double maxPoissonRate) {
    double spikeTime = lastSpikeTime;
    do spikeTime = spikeTime - log
    poiGen
}

void predetermine_poisson_spikes(vector<double> &spikeTime, vector<double> &inputRate, double maxPoissonRate, string outputFn="", bool saveInputToFile = false);
