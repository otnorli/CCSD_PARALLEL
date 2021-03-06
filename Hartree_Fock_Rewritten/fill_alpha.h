#ifndef FILL_ALPHA_H
#define FILL_ALPHA_H

#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

class Fill_Alpha
{
public:
    Fill_Alpha(int n_N, vec zZz, string b_s, int mat_size);
    int n_Nuclei;
    vec n_Basis;
    vec n_Orbitals;
    vec Z;
    mat alpha;
    mat c;
    mat Potenser;
    vec Basis_Per_Orbital;
    int Matrix_Size;

    int tot_Basis_Functions;
    string Basis_Set;

    mat Fyll_Opp_Alpha();
    mat Fyll_Opp_c();
    vec Fyll_Opp_Nr_Basis_Functions();
    vec Fyll_Opp_Antall_Orbitaler();
    mat Fyll_Opp_Potenser();

};

#endif // FILL_ALPHA_H
