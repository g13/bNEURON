#ifndef INPUT_ARGS_CA1_H
#define INPUT_ARGS_CA1
#include "input_args.h"
#include <bitset>

struct input_args_CA1 : input_args {
    string libFile;
    double vThres;
    double vRest;
    double vTol;
    double rLinear;
    double vBuffer;
    double dendClampRatio;
    double ignoreT;
    double tRef;
    double trans;
    double dtrans;
    double dtau;
    double trans0;
    double rBiLinear;
    int kVStyle;
    int afterCrossBehavior;
    bool spikeShape;
    bool pas;
    bool sliceDebugPlot;
    int i;
    bool dtSquare;
    bool getDendV;
    bool tIn;
    unsigned int vInit;
    std::bitset<7> mode;
    vector<double> clusterClampRatio;
    vector<int> clusterSizeCDF;
    vector<vector<int>> clusterDend;
    input_args_CA1();
    void setbit();
    using input_args::read;
    void get_rd(); 
};

typedef struct input_args_CA1 InputArgsCA1;
#endif
