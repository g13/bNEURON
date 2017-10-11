#ifndef MATFUNC_H
#define MATFUNC_H
#include <iostream>
#include "mat.h"
#include "matrix.h"
#include "assert.h"
#include "typedefs.h"

template<typename T>
void del2d(T **ptr) {
    delete []ptr;
}

template<typename T, typename U>
void del3d(T ***ptr, U *dimSize) {
    U i;
    for (i=0; i<dimSize[0]; i++) delete []ptr[i];
    delete []ptr;
}

template<typename T, typename U>
void del4d(T ****ptr, U *dimSize) {
    U i,j;
    for (i=0; i<dimSize[0]; i++) {
        for (j=0; j<dimSize[1]; j++)
            delete []ptr[i][j];
        delete []ptr[i];
    }
    delete []ptr;
}

template<typename T, typename U>
void del5d(T *****ptr, U *dimSize) {
    U i,j,k;
    for (i=0; i<dimSize[0]; i++) {
        for (j=0; j<dimSize[1]; j++){
            for (k=0; k<dimSize[2]; k++)  
                delete []ptr[i][j][k];
            delete []ptr[i][j];
        }
        delete []ptr[i];
    }
    delete []ptr;
}

template<typename T, typename U>
void del6d(T ******ptr, U *dimSize) {
    U i,j,k,l;
    for (i=0; i<dimSize[0]; i++) {
        for (j=0; j<dimSize[1]; j++) {
            for (k=0; k<dimSize[2]; k++) {
                for (l=0; l<dimSize[3]; l++)
                    delete []ptr[i][j][k][l];
                delete []ptr[i][j][k];
            }
            delete []ptr[i][j];
        }
        delete []ptr[i];
    }
    delete []ptr;
}

template<typename T, typename U>
void disp1d(T v, U dimSize) {
    U i;
    std::cout << "{";
    for (i=0; i<dimSize; i++) {
        std::cout << v[i];
        if (i<dimSize-1) std::cout << ", ";
    } 
    std::cout << "}" << std::endl;
    std::cout << std::endl;
}
template<typename T, typename U>
void disp2d(T v, U *dimSize) {
    U i,j;
    std::cout << "{" << std::endl;
    for (i=0; i<dimSize[0]; i++) {
        std::cout << "{";
        for (j=0; j<dimSize[1]; j++) {
            std::cout << v[i][j];
            if (j<dimSize[1]-1) std::cout << ", ";
        }
        if (i==dimSize[0]-1) std::cout << "}}" << std::endl;
        else std::cout << "}," << std::endl;
    } 
}
template<typename T, typename U>
void disp3d(T v, U *dimSize) {
    U i,j,k;
    std::cout << "{" << std::endl;
    for (i=0; i<dimSize[0]; i++) {
        std::cout << "{";
        for (j=0; j<dimSize[1]; j++) {
            std::cout << "{";
            for (k=0; k<dimSize[2]; k++) {
                std::cout << v[i][j][k];
                if (k<dimSize[2]-1) std::cout << ", ";
            }
            if (j==dimSize[1]-1) std::cout << "}";
            else std::cout << "},";
        }
        if (i==dimSize[0]-1) std::cout << "}}" << std::endl;
        else std::cout << "}," << std::endl;
    } 
}
template<typename T, typename U>
void disp4d(T v, U *dimSize) {
    U i,j,k,l;
    std::cout << "{" << std::endl;
    for (i=0; i<dimSize[0]; i++) {
        std::cout << "{";
        for (j=0; j<dimSize[1]; j++) {
            std::cout << "{";
            for (k=0; k<dimSize[2]; k++) {
                std::cout << "{";
                for (l=0; l<dimSize[3]; l++) {
                    std::cout << v[i][j][k][l];
                    if (l<dimSize[3]-1) std::cout << ", ";
                }
                if (k==dimSize[2]-1) std::cout << "}";
                else std::cout << "},";
            }
            if (j==dimSize[1]-1) std::cout << "}";
            else std::cout << "},";
        } 
        if (i==dimSize[0]-1) std::cout << "}}" << std::endl;
        else std::cout << "}," << std::endl;
    }
}

template<typename T> 
void pointer2d(T** &ptr, T* source, size *dimSize) {
    size i;
    ptr = new T*[dimSize[0]];
    for (i=0; i<dimSize[0]; i++) {
        ptr[i] = source + dimSize[1]*i;
    }
}
template<typename T> 
void pointer3d(T*** &ptr, T* source, size *dimSize) {
    size i,j;
    ptr = new T**[dimSize[0]];
    for (i=0; i<dimSize[0]; i++) {
        ptr[i] = new T*[dimSize[1]];
        for (j=0; j<dimSize[1]; j++) {
            ptr[i][j] = source + dimSize[1]*dimSize[2]*i \
                                          + dimSize[2]*j;
        }
    }
}
template<typename T> 
void pointer4d(T**** &ptr, T* source, size *dimSize) {
    size i,j,k;
    ptr = new T***[dimSize[0]];
    for (i=0; i<dimSize[0]; i++) {
        ptr[i] = new T**[dimSize[1]];
        for (j=0; j<dimSize[1]; j++) {
            ptr[i][j] = new T*[dimSize[2]];
            for (k=0; k<dimSize[2]; k++) {
                ptr[i][j][k] = source + dimSize[1]*dimSize[2]*dimSize[3]*i \
                                                 + dimSize[2]*dimSize[3]*j \
                                                            + dimSize[3]*k;
            }
        }
    }
}

template<typename T> 
void pointer5d(T***** &ptr, T* source, size *dimSize) {
    size i,j,k,l;
    ptr = new T****[dimSize[0]];
    for (i=0; i<dimSize[0]; i++) {
        ptr[i] = new T***[dimSize[1]];
        for (j=0; j<dimSize[1]; j++) {
            ptr[i][j] = new T**[dimSize[2]];
            for (k=0; k<dimSize[2]; k++) {
                ptr[i][j][k] = new T*[dimSize[3]];
                for (l=0; l<dimSize[3]; l++) {
                    ptr[i][j][k][l] = source + dimSize[1]*dimSize[2]*dimSize[3]*dimSize[4]*i \
                                                     + dimSize[2]*dimSize[3]*dimSize[4]*j \
                                                                + dimSize[3]*dimSize[4]*k \
                                                                           + dimSize[4]*l;
                }
            }
        }
    }
}

template<typename T> 
void pointer6d(T****** &ptr, T* source, size *dimSize) {
    size i,j,k,l,m;
    ptr = new T*****[dimSize[0]];
    for (i=0; i<dimSize[0]; i++) {
        ptr[i] = new T****[dimSize[1]];
        for (j=0; j<dimSize[1]; j++) {
            ptr[i][j] = new T***[dimSize[2]];
            for (k=0; k<dimSize[2]; k++) {
                ptr[i][j][k] = new T**[dimSize[3]];
                for (l=0; l<dimSize[3]; l++) {
                    ptr[i][j][k][l] = new T*[dimSize[4]];
                    for (m=0; m<dimSize[4]; m++) {
                        ptr[i][j][k][l][m] = source + dimSize[1]*dimSize[2]*dimSize[3]*dimSize[4]*dimSize[5]*i \
                                                                +dimSize[2]*dimSize[3]*dimSize[4]*dimSize[5]*j \
                                                                           +dimSize[3]*dimSize[4]*dimSize[5]*k \
                                                                                      +dimSize[4]*dimSize[5]*l \
                                                                                                 +dimSize[5]*m;
                    }
                }
            }
        }
    }
}

void openMat(MATFile* &pmat, const char* file) {
    pmat = matOpen(file, "r");
    if (pmat == NULL) {
        std::cout << "Error opening file: " << file << std::endl;
        abort();
    } 
}
void closeMat(MATFile* pmat, const char* file) {
    if (matClose(pmat) != 0) {
        std::cout << "Error closing file: " << file << std::endl;
        abort();
    }
}
void getArrayDims(mxArray* &tmp, size* dimSize, size &arraySize, const char* name)
{
    unsigned int i,ndim;
    ndim = mxGetNumberOfDimensions(tmp);
    std::cout << "Dimension of array " << name << ": ";
    for (i=0; i<ndim; i++){
        dimSize[i] = mxGetDimensions(tmp)[ndim-i-1];
        std::cout << dimSize[i];
        if (i<ndim-1) std::cout << ", ";
        else std::cout << std::endl;
    }
    arraySize = mxGetNumberOfElements(tmp);
}
void readArray(mxArray* &tmp, const char* name, size* dimSize, size &arraySize, MATFile* &pmat, const char* file)
{
    tmp = matGetVariable(pmat, name);
    if (tmp == NULL) {
        std::cout << "Error reading array " << name << " in " << file << std::endl;
        mxDestroyArray(tmp);
        closeMat(pmat, file);
        abort();
    }
    getArrayDims(tmp, dimSize, arraySize, name);
}
template<typename T>
void readVar(T &var, const char *name, MATFile* pmat, const char* file)
{
    mxArray* tmp;
    tmp = matGetVariable(pmat, name);
    if (tmp == NULL) {
        std::cout << "Error reading variable " << name << " in " << file << std::endl;
        mxDestroyArray(tmp);
        closeMat(pmat, file);
        abort();
    }
    assert(mxGetNumberOfElements(tmp)==1);
    var = static_cast<T>(mxGetScalar(tmp));
    std::cout << "Variable " << name << " = " << var << std::endl;
    mxDestroyArray(tmp);
}

template<typename T>
void getFieldVar(T &var, mxArray *matStruct, unsigned int index, const char* fieldname)
{
    mxArray *tmp;
    tmp = mxGetField(matStruct,index,fieldname);
    var = static_cast<T>(mxGetScalar(tmp));
    std::cout << fieldname << ": " << var << std::endl;
    //mxDestroyArray(tmp);
}
template<typename T>
void getIthElementOfFieldArray(T &var, mxArray* para, unsigned int ith,  const char* fieldname) {
    mxArray* tmp;
    tmp = mxGetField(para,0,fieldname);
    var = static_cast<T>(*(mxGetPr(tmp)+ith));
    std::cout << fieldname << ": " << var << std::endl;
    //mxDestroyArray(tmp);
}
#endif
