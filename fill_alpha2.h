#ifndef FILL_ALPHA2_H
#define FILL_ALPHA2_H

#include <iostream>
#include "../../home/ole/Desktop/include/armadillo"

using namespace std;
using namespace arma;

class fill_alpha2
{
public:
    fill_alpha2(int n_N, vec zZz, string b_s, int mat_size, int max_num_b);

    int n_Nuclei;
    vec n_Basis;
    vec n_Orbitals;
    vec Z;
    mat alpha;
    mat c;
    mat Potenser;
    vec Basis_Per_Orbital;
    int Matrix_Size;
    int max_numb_bas;

    int tot_Basis_Functions;
    string Basis_Set;

    mat Fyll_Opp_Alpha();
    mat Fyll_Opp_c();
    vec Fyll_Opp_Nr_Basis_Functions();
    vec Fyll_Opp_Antall_Orbitaler();
    mat Fyll_Opp_Potenser();
};

#endif // FILL_ALPHA2_H
