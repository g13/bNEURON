#include "input_args.h"
#include <cstdlib>
#include <cmath>
using std::ifstream;
using std::streampos;
using std::ios;

input_args::input_args() {
    irregInputLevels = false;
    dtVarLevels = false;
    tVarLevels = false;
    tVar = false;
    pVar = false;
    extendVar = false;
    initialGenerics.add_options()
        ("configFn,c", po::value<string>(&configFn)->default_value("config"),"filename of configuration");
    configFileOptions.add_options()
        ("theme,m", po::value<string>(&theme), "theme of simulation")
        ("irregInputLevels", "irregular input levels")
        ("tVarLevels", "variable input level sim time")
        ("dtVarLevels", " use different tstep for each input level")
        ("nInput", po::value<unsigned long>(&nInput)->default_value(1), " number of input sites")
        ("inputLevelsFn", po::value<string>(&inputLevelsFn),"filename of irregular input levels (need to set irregInputLevels, dtVarLevels or tVarLevels flag)")
        ("inputLinspace,i", po::value<double>()->multitoken()->composing(), "evenly spaced input level, start, total steps, end")
		("runTime,t", po::value<double>(), "simulation time for each input level")
        ("dt", po::value<double>(), "dt time step for each run")
        ("tVar", " use temporally variable input")
        ("pVar", " use positionally variable input")
        ("extendVar", " extend pVar or tVar to other input levels")
        ("inputFn", po::value<string>(&inputFn),"filename of variable input table (need to set tVar or pVar flag)");
    cmdLineOptions.add(configFileOptions);
    cmdLineOptions.add(initialGenerics);
}
int input_args::read(int argc, char **argv) {
    ifstream cfgFile;
    vector<vector<double>> vars;
    vector<unsigned long> column;
    vector<double> inputLinspace;
    po::store(po::parse_command_line(argc, argv, cmdLineOptions), vm);
    po::notify(vm);
    cfgFile.open(configFn);
    if (cfgFile) {
        store(po::parse_config_file(cfgFile, configFileOptions, true), vm);
    } else {
        cout << "cannot open configuration file: " << configFn << endl; 
        return 0;
    }
    // input levels
    if (vm.count("irregInputLevels")) {
        irregInputLevels = true;
        vars.push_back(inputLevel);
        if (vm.count("inputLinspace")) {
            cout << " using irregInputLevels flag, ignoring inputLinspace" << endl;
        }
    } else {
        if (!vm.count("inputLinspace")) {
            cout << "must provide a inputLinspace {begin,nlevels(,end)}" << endl;
            return 0;
        } else {
            inputLinspace = vm["inputLinspace"].as<vector<double>>();
            inputLevel[0] = inputLinspace[0];
            if (inputLinspace[1] > 1.5) {
                // inputLinspace[1] should be a "double" integer
                double dLevel = (inputLinspace[2] - inputLinspace[0])/(inputLinspace[1]-1);
                for (int i=1; i<inputLinspace[1]; i++) {
                    inputLevel.push_back(inputLinspace[0] + dLevel*i);
                }
                assert(abs(inputLevel.back()-inputLinspace[2]) < 1e-10);
            }
        }
    }
    if (vm.count("tVarLevel")) {
        vars.push_back(runTime);
        tVarLevels = true;
        if (vm.count("runTime")) {
            cout << " using tVarLevel flag, ignoring runTime" << endl;
        }
    } else {
        if (!vm.count("runTime")) {
            cout << " must provide a single runTime or an array of runTime" << endl;
            return 0;
        } else {
            runTime.assign(inputLevel.size(), vm["runTime"].as<double>());
            if (runTime[0] <= 0.0 ) {
                cout << " runTime must be positive" << endl;
                return 0;
            }
        }
    }
    if (vm.count("dtVar")) {
        vars.push_back(tstep);
        dtVarLevels = true;
        if (vm.count("dt")) {
            cout << " using dtVar flag, ignoring single dt" << endl;
        }
    } else {
        if (!vm.count("dt")) {
            cout << "must provide a single dt or an array of dt" << endl;
            return 0;
        } else {
            tstep.assign(inputLevel.size(), vm["dt"].as<double>());
            if (tstep[0] <= 0.0) {
                cout << " dt must be positive" << endl;
                return 0;
            }
        }
    }
    if (irregInputLevels || tVarLevels || dtVarLevels) {
        if (!read_input_table(inputLevelsFn,column,vars)) {
            return 0;
        } else {
            cout << "# input levels: " << inputLevel.size() << endl;
            if (tVarLevels) {
                assert(runTime.size()==inputLevel.size());
            }
            if (dtVarLevels) {
                assert(tstep.size()==inputLevel.size());
            }
        }
        vars.clear();
    }
    // per input level
    if (vm.count("pVar")) {
        pVar = true;
    }
    if (vm.count("tVar")) {
        tVar = true;
    }
    if (vm.count("extendVar")) {
        extendVar = true;
    }
    if (tVar && pVar && (tVarLevels || dtVarLevels || extendVar)) {
        cout << " (tVar + pVar) is not compatible to tVarLevels, dtVarLevels, extendVar flag" << endl;
        cout << " run additional instances of simulation instead" << endl;
        return 0;
    }
    if (inputFn.empty()){
        cout << "no file for tVar available" << endl;
    } else {
        column.clear();
        if (pVar && tVar && !tVarLevels && !dtVarLevels && !extendVar) {
            inputMode = "pt";
            column.push_back(round(runTime[0]/tstep[0]));
            if (!read_input_table(inputFn, column, input)) {
                return 0;
            }
            assert(input.size() == inputLevel.size());
        }
        if (!pVar && tVar && (tVarLevels || dtVarLevels)) {
            if (tVarLevels && !dtVarLevels) {
                inputMode = "t-T";
            }
            if (!tVarLevels && dtVarLevels) {
                inputMode = "t-dt";
            }
            for (int i=0; i<inputLevel.size(); i++) {
                column.push_back(round(runTime[i]/tstep[i]));
            }
            input.assign(inputLevel.size(),vector<double>());
            if (!read_input_table(inputFn, column, input)) {
                return 0;
            }
        }
        if (pVar && !tVar && (tVarLevels || dtVarLevels)) {
            if (extendVar) {
                if (tVarLevels && !dtVarLevels) {
                    inputMode = "P-T";
                }
                if (!tVarLevels && dtVarLevels) {
                    inputMode = "P-dt";
                }
                input.assign(inputLevel.size(), vector<double>());
            } else {
                if (tVarLevels && !dtVarLevels) {
                    inputMode = "p-T";
                }
                if (!tVarLevels && dtVarLevels) {
                    inputMode = "p-dt";
                }
                column.push_back(nInput);
            }
            if (!read_input_table(inputFn, column, input)) {
                return 0;
            }
            if (inputMode == "P-T" || inputMode == "P-dt") {
                for (int i=0; i<inputLevel.size(); i++) {
                    assert(input[i].size() == nInput);
                }
            }
            if (inputMode == "p-T" || inputMode == "p-dt") {
                assert(input.size() == 1);
                assert(input[0].size() == nInput);
            }
        }
        if ((!pVar && tVar || !tVar && pVar) && !tVarLevels && !dtVarLevels) {
            column.clear();
            if (!extendVar) {
                if (pVar) {
                    inputMode = "p";
                    column.push_back(nInput);
                }
                if (tVar) {
                    inputMode = "t";
                    column.push_back(round(runTime[0]/tstep[0]));
                }
            } else {
                if (pVar) {
                    inputMode = "P";
                    input.assign(inputLevel.size(),vector<double>());
                }
                if (tVar) {
                    inputMode = "T";
                    input.assign(inputLevel.size(),vector<double>());
                }
            }
            if (!read_input_table(inputFn, column, input)) {
                return 0;
            }
            if (inputMode == "t") {
                assert(input.size() == 1);
            }
            if (inputMode == "T") {
                assert(input.size() == inputLevel.size());
                for (int i=0; i<input.size(); i++) {
                    assert(input[i].size() == input[0].size());
                }
                assert(input[0].size() == round(runTime[0]/tstep[0]));
                assert(input[0].size() == round(runTime.back()/tstep.back()));
            }
            if (inputMode == "p") {
                assert(input.size() == 1);
            }
            if (inputMode == "P") {
                assert(input.size() == inputLevel.size());
                for (int i=0; i<input.size(); i++) {
                    assert(input[i].size() == nInput);
                }
            }
        }
    }
    return 1;
}
int read_input_table(string tableFn, vector<unsigned long> column, vector<vector<double>> &arr) {
    // arr are arrays needing to be fed.
    ifstream tableFile;
    tableFile.open(tableFn, ios::binary);
    streampos fBegin, fEnd;
    unsigned long length;
    if (tableFile) {
        if (column.size() == 1 && arr.size() == 0) {
            // uniform column table not knowing rows
            for (int i=0; i<arr.size(); i++) {
                arr.push_back(vector<double>(column[0],0));
                tableFile.read((char*)(&arr[i][0]),sizeof(double)*column[i]);
                cout << i << ": " << arr[i][0];
                for (int j=1; j<column[i]; j++) {
                    cout << ", " << arr[i][j]; 
                }
                cout << endl;
            }
        } else {
            // rows known
            if (column.size() == 0) {
                // known # rows, get uniform columns
                tableFile.seekg(0,ios::end);
                fEnd = tableFile.tellg();
                tableFile.seekg(0,ios::beg);
                fBegin = tableFile.tellg();
                length = (fEnd - fBegin)/(sizeof(double)*arr.size());
                column.assign(arr.size(),length);
            } else {
                if (column.size() == 1 && arr.size() > 1) {
                    column.assign(arr.size(),column[0]);
                }
            } // else non-uniform columns
            assert( column.size() == arr.size() );
            for (int i=0; i<arr.size(); i++) {
                arr[i].assign(column[i],0);
                tableFile.read((char*)(&arr[i][0]),sizeof(double)*column[i]);
                cout << i << ": " << arr[i][0];
                for (int j=1; j<column[i]; j++) {
                    cout << ", " << arr[i][j]; 
                }
                cout << endl;
            }
        }
    }
    tableFile.close();
    return 1;
}
