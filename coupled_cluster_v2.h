#ifndef COUPLED_CLUSTER_V2_H
#define COUPLED_CLUSTER_V2_H

#include <iostream>
#include "../../home/ole/Desktop/include/armadillo"
#include "hartree_fock_solver.h"
#include "coupled_cluster_integrals.h"

using namespace std;
using namespace arma;


class coupled_cluster_v2
{
public:
    coupled_cluster_v2(int n_N, vec zz, mat rr, string B_s, int n_Elec);
    
    // Input variables:
    int n_Nuclei;
    vec Z;
    mat R;
    string Basis_Set;
    int Matrix_Size;
    int n_Electrons;
    int unocc_orb;
    int iter;

    // Ting vi kommer til å trenge som globale variabler i denne klassen for beregningene:
    vector<cube> denom_abij, MO; //  integ;
    mat t1, fs, denom_ai, t1_new;
    mat D1, D2, D3;
    field<mat> W_4, W_3, W_2, W_1;
    field<mat> t2;
    field<mat> T_2, T_2_new;
    field<mat> tau, tau2, integ, integ2;
    
    // Energier:
    double E_old, E_new;

    // Generelle CCSD funksjoner
    double CCSD(double toler, bool print_stuff);
    int EqualFunc(int a, int b);
    double Calc_Energy();

    // Smart t2 manouvering, not so smart jet
    void Init_t2();
    vec t2_b_values;
    vec t2_j_values;

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

#endif // COUPLED_CLUSTER_V2_H
