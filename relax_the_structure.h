#ifndef RELAX_THE_STRUCTURE_H
#define RELAX_THE_STRUCTURE_H

#include <iostream>
#include "../../home/ole/Desktop/include/armadillo"
#include "coupled_cluster.h"

using namespace std;
using namespace arma;

class Relax_The_Structure
{
public:
    Relax_The_Structure(int n_N, int n_elec, mat rr, string B_s, vec zz, string meth);
    int n_Nuclei,n_Electrons, Matrix_Size;
    vec Z, n_Basis, Number_Of_Orbitals,eigenvalues_O, eigenvalues_F;
    string Basis_Set;
    mat alpha, c, R, EK, O, D, V, F, eigenvektors_O, eigenvektors_F, P, C, Potenser;
    string Method;

    mat Forces_R; // "Relaxation forces"
    mat return_R;

    mat Relaxation_HF(double toler);
    mat Relaxation_CCSD(double toler);
};

#endif // RELAX_THE_STRUCTURE_H
