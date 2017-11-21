#include <boost/program_options.hpp>
#include <cstring>
#include <fstream>
#include <iostream>
#include <vector>
// my naming convention, functions and methods "_", variables "caMel", types "CaMel"
using std::cout;
using std::endl;
using std::string;
using std::vector;
namespace po = boost::program_options;

struct input_args {
	po::options_description initialGenerics, cmdLineOptions, configFileOptions;
    po::variables_map vm;
    string theme;
    unsigned int seed;
    unsigned long nInput;
    vector<double> tstep, runTime, inputLevel;
    vector<vector<double>> input;
    string inputMode;
    bool irregInputLevels, dtVarLevels, tVarLevels;
    bool tVar, pVar, extendVar;
    string configFn, inputLevelsFn, inputFn;
    vector<vector<unsigned long>> postID;
    vector<vector<unsigned long>> preID;
    vector<vector<unsigned int>> dendID;
    input_args();
    int read(int argc, char **argv);
};
typedef struct input_args InputArgs;
//struct input_args_B : input_args {
//
//}

int read_input_table(string tableFn, vector<unsigned long> columns, vector<vector<double>> &vars);
