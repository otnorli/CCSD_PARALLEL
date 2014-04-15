#ifndef INITIALIZER_H
#define INITIALIZER_H

#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include <iomanip>
#include "../../home/ole/Desktop/include/armadillo"
#include "hartree_fock_solver.h"
#include "hartree_integrals.h"
#include "ccsd_memory_optimized.h"

using namespace std;
using namespace arma;

class Initializer
{
public:
    Initializer(int n_N, int n_E, string B_S, string met, bool r_pos, mat rr, vec zz, double con_crit, int ran, int siz, bool us_ang);
    double Go(bool printie);

    int n_Nuclei, n_Electrons;
    double convergance_criteria;
    string Basis_Set, Method;
    mat R;
    vec Z;

    bool use_angstrom;
    bool Relax_Pos;

    // MPI
    int rank, size;

};

#endif // INITIALIZER_H
