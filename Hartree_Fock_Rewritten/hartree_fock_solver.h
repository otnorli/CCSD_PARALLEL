#ifndef HARTREE_FOCK_SOLVER_H
#define HARTREE_FOCK_SOLVER_H

#include <iostream>
#include <armadillo>
#include "matrix_size_setter.h"
#include "fill_alpha.h"
#include "hartree_integrals.h"

using namespace std;
using namespace arma;

class Hartree_Fock_Solver
{
public:
    // Input parametre
    Hartree_Fock_Solver(int n_N, vec zz, mat rr, string B_s);
    int n_Nuclei;
    vec Z;
    mat R;
    string Basis_Set;

    // Globale parametre
    vec n_Basis;
    int n_Electrons;
    mat alpha;
    mat c;
    vec Number_Of_Orbitals;
    int Matrix_Size;
    mat EK, O, D, V, F, eigenvektors_O, eigenvektors_F, P;
    vec C;
    vec eigenvalues_O, eigenvalues_F;
    mat Potenser;


    // Funksjoner
    double get_Energy(int n_S);
    vec Normalize(vec Vektorn, mat Matrisen);


};

#endif // HARTREE_FOCK_SOLVER_H
