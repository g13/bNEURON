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
    unsigned int vInit;
    std::bitset<6> mode;
    input_args_CA1();
    void setbit();
    using input_args::read;
};

typedef struct input_args_CA1 InputArgsCA1;
#endif
