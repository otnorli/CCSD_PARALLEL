#ifndef INITIALIZER_H
#define INITIALIZER_H

#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include <iomanip>
#include "../../home/ole/Desktop/include/armadillo"
#include "hartree_fock_solver.h"
#include "coupled_cluster.h"
#include "relax_the_structure.h"
#include "hartree_integrals.h"
#include "coupled_cluster_v2.h"
#include "ccsd_even.h"
#include "ccsd_memory_optimized.h"
#include "ccsd_v2_optimized.h"

using namespace std;
using namespace arma;

class Initializer
{
public:
    Initializer(int n_N, int n_E, string B_S, string met, bool r_pos, mat rr, vec zz, double con_crit, string r_uns);
    double Go();

    int n_Nuclei, n_Electrons;
    double convergance_criteria;
    string Basis_Set, Method;
    mat R;
    vec Z;
    string R_Units;

    bool Relax_Pos;

};

#endif // INITIALIZER_H
