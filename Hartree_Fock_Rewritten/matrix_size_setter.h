#ifndef MATRIX_SIZE_SETTER_H
#define MATRIX_SIZE_SETTER_H

#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

class matrix_size_setter
{
public:
    matrix_size_setter(vec zZz, string sys_string, int n_N);
    vec Z;
    string System_String;
    int Set_Matrix_Size();
    int n_Nuclei;
    double value, tempvalue;

};

#endif // MATRIX_SIZE_SETTER_H
