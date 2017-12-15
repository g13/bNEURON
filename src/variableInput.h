#include "Yale_NEURON_PyAPI.h"
#include "input_args_CA1.h"
#include "jumpy_linear.h"
#include "linear_bilinear.h"
#include <ctime>
#include <fstream>
#include <iostream>
#include <boost/program_options.hpp>
#include <cstring>

using std::cout;
using std::endl;
using std::vector;
using std::string;
using std::ofstream;
using std::ifstream;
using std::to_string;
using std::ios;

template<typename T>
void size_data_write(ofstream &file, vector<T>** data, int nrow, size nSize, size istart){
    file.write((char*)&(nSize), sizeof(size));
    if (nSize) {
        for (int i=0; i<nrow; i++) {
            //cout << "size: " << nSize << " vec size " << data[i]->size() << " start at " << istart << endl;
            file.write((char*)&(data[i]->at(istart)), nSize*sizeof(T));
        }
    } else {
        cout << "empty array" << endl;
    }
}
template<typename T>
void nsection_write(ofstream &file, vector<T>** data, int nrow, vector<size> &istart, vector<size> &nSize) {
    int nsection = nSize.size();
    file.write((char*)&(nsection), sizeof(int));
    if  (nsection) {
        for (int i=0; i<nsection; i++) {
            size_data_write(file, data, nrow, nSize[i], istart[i]);
        }
    } else {
        cout << "no section recorded" << endl;
    }
}

void write_data_file(ofstream &data_file, vector<double>** output, vector<vector<double>> &dendV, unsigned int nt0, unsigned int nt, int nSyn);

int main(int argc, char **argv);
