#ifndef COUPLED_CLUSTER_H
#define COUPLED_CLUSTER_H

#include <iostream>
#include "../../home/ole/Desktop/include/armadillo"
#include "hartree_fock_solver.h"
#include "coupled_cluster_integrals.h"

using namespace std;
using namespace arma;

class coupled_cluster
{
public:
    coupled_cluster(int n_N, vec zz, mat rr, string B_s, int n_Elec);
    double E_old, E_new;

    // De ulike coupled cluster metodene implementert
    double CCSD(double toler, bool print_stuff);
    double CCSDt(double toler);

    // Generelle intput ting
    int n_Nuclei;
    vec Z;
    mat R;
    string Basis_Set;
    int Matrix_Size;
    int n_Electrons;
    int unocc_orb;

    // Ting vi kommer til å trenge som globale variabler i denne klassen for beregningene:
    vector<cube> t2, denom_abij, t2_new, tau, integ, G2; //  integ;
    mat t1, fs, denom_ai, t1_new;

    mat F1, F2, F3;
    mat D1, D2, D3;
    vector<cube> G1, G3, G4;
    field<mat> W_4, W_1;
    field <mat> T_2, T_2_new;

    // CCSD(T) ting
   //cube<cube> t3;

    // Funksjoner
    double tau1(int a, int b, int i, int j);
    double tau2(int a, int b, int i, int j);
    double EqualFunc(int a, int b);
    void Fill_G1();
    void Fill_F1();
    void Fill_G3();
    void Fill_F3();
    void Fill_G2();
    void Fill_F2();
    void Fill_G4();
    void Fill_t1_new();
    void Fill_t2_new();
    double Calc_Energy();
    void Fill_tau();

    // DIIS
    void DIIS(int t_degree);
    void Initialize_DIIS();
    int iter, number_elements_DIIS;
    vector<mat> Error_Stored1;
    vector<mat> Answer_Stored1;
    vector<mat> Error_Stored2;
    vector<mat> Answer_Stored2;
    mat DIIS_B1, delta_p1;
    vec DIIS_Z1, DIIS_c1;
    mat DIIS_B2, delta_p2, DIIS_t2;
    vec DIIS_Z2, DIIS_c2;
    bool DIIS_does;


    // Speedup variabler
    mat speed;
};

#endif // COUPLED_CLUSTER_H
