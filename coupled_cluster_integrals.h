#ifndef COUPLED_CLUSTER_INTEGRALS_H
#define COUPLED_CLUSTER_INTEGRALS_H

#include <iostream>
#define ARMA_NO_DEBUG
#include "../../home/ole/Desktop/include/armadillo"

using namespace std;
using namespace arma;

class coupled_cluster_integrals
{
public:
    coupled_cluster_integrals(int mat_siz, mat ccc, int n_E);
    int Matrix_Size;

    vec Indexed_integrals;

    void Mount_Indexed_Integrals(vec intsrr);

    void Return_Integrals(vec Integrals, int Version);
    void Prepear_Integrals();
    void Prepear_Integralsv2();
    int Return_Integral_Index(int a, int b, int i, int j);

    double EqualFunc(int a, int b);
    vector<cube> Integrals;
    mat Return_FS(vec eigval);
    int n_Electrons;
    mat c;
    field<mat> temp_cube33;
    field<mat> temp_cube3;
    vector <cube> temp4;
    field<mat> Rear_cube, Rear_cube2, Rear_cube3, Rear_cube4, Rear_cube6, Rear_cube7, Rear_cube8, Rear_cube9, Rear_cube10;
    field<mat> Return_Rearranged_Integrals2();
    field<mat> Return_Rearranged_Integrals3();
    field<mat> Return_Rearranged_Integrals4();
    field<mat> Return_Rearranged_Integrals5();
    field<mat> Return_Rearranged_Integrlas6();
    field<mat> Return_Rearranged_Integrals7();
    field<mat> Return_Rearranged_Integrals8();
    field<mat> Return_Rearranged_Integrals9();
    field<mat> Return_Rearranged_Integrals10();
    field<mat> Return_Integrals_Part_2();
    field<mat> Return_V2_integrals();

    void Delete_Everything();
};

#endif // COUPLED_CLUSTER_INTEGRALS_H
