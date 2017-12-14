#include "Yale_NEURON_PyAPI.h"
#include "input_args_CA1.h"
#include "jumpy_linear.h"
#include "linear_bilinear.h"
#include <ctime>
#include <fstream>
#include <iostream>
#include <boost/program_options.hpp>
#include <cstring>

template<T>
void size_data_write(ofstream file, vector<T> *data, int nrow, size nSize, size istart){
    file.write((char*)&(nSize), nSize, sizeof(size));
    for (int i=0; i<nrow; i++) {
        file.write((char*)&(data[i]->operator[istart]), nSize*sizeof(T));
    }
}
template<T>
void nsection_write(ofstream file, vector<T>* data, int nrow, vector<size> istart, vector<size> nSize) {
    int nsection = start.size();
    file.write((char*)&(nrow), sizeof(int));
    for (int i=0; i<nsection; i++) {
        size_data_write(file, data, nrow, nSize[i], start[i]);
    }
}
write_data_file(ofstream data_file, vector<double> *output, vector<double> dendV, unsigned int nt0, unsigned int nt, int nSyn);
int main(int argc, char **argv);
