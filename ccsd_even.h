#ifndef CCSD_EVEN_H
#define CCSD_EVEN_H

#include <iostream>
#include "../../home/ole/Desktop/include/armadillo"
#include "hartree_fock_solver.h"
#include "coupled_cluster_integrals.h"

using namespace std;
using namespace arma;

class ccsd_even
{
public:
    ccsd_even(int n_N, vec zz, mat rr, string B_s, int n_Elec);

    // Input variables:
    int n_Nuclei;
    vec Z;
    mat R;
    string Basis_Set;
    int Matrix_Size;
    int n_Electrons;
    int unocc_orb;
    int iter;
    int Speed_Elec, Speed_Occ;

    // Ting vi kommer til å trenge som globale variabler i denne klassen for beregningene:
    vector<cube> denom_abij; //  integ;
    mat t1, fs, denom_ai, t1_new;
    field<mat> T_2_new;
    field<mat> integ;

    // Hypercompact form
    field<mat> W_4, t2, t2_new, tau3, tau4, W_1, W_2, W_3, DEN_ABIJ;
    mat D2, D3, T_1, T_1_new, D1;
    mat FS_AI, FS_IJ, FS_AB, DEN_AI;
    mat W_22;

    // Split up integral for much less storage, even if meny more variables, storage reduced by >60%
    field<mat> integ2, integ3, integ4, integ5, integ6, integ7, integ8, integ9, integ10;

    // Energier:
    double E_old, E_new;

    // Generelle CCSD funksjoner
    double CCSD(double toler, bool print_stuff);
    int EqualFunc(int a, int b);
    double Calc_Energy();

    // Intermediates
    void Fill_W1();
    void Fill_W2();
    void Fill_W3();
    void Fill_W4();
    void Fill_F1();
    void Fill_F2();
    void Fill_F3();
    void Fill_tau();

    // Amplitudes
    void Fill_t1_new();
    void Fill_t2_new();
};

#endif // CCSD_EVEN_H
