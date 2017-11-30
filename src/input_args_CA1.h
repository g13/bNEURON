#ifndef INPUT_ARGS_CA1_H
#define INPUT_ARGS_CA1
#include "input_args.h"

struct input_args_CA1 : input_args {
    string libFile;
    string paraFile; 
    double vThres;
    double vRest;
    double vTol;
    double rLinear;
    double vBuffer;
    double dendClampRatio;
    double ignoreT;
    double tRef;
    double trans;
    double trans0;
    double rBiLinear;
    unsigned int vInit;
    input_args_CA1();
};

typedef struct input_args_CA1 InputArgsCA1;
#endif
