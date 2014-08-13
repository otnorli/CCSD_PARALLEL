#ifndef CCSDT_H
#define CCSDT_H

#include <iostream>
#define ARMA_NO_DEBUG
#include "../../home/ole/Desktop/include/armadillo"
#include "mpi.h"
#include "hartree_fock_solver.h"
#include <time.h>

using namespace std;
using namespace arma;

class CCSDT
{
public:
    CCSDT(int n_N, vec zz, mat rr, string B_s, int n_Elec, int ran, int siz, Hartree_Fock_Solver *Hartfock, bool frez);

    // Input
    int Method_Nr;
    int n_Nuclei;
    vec Z;
    mat R;
    string Basis_Set;
    int Matrix_Size;
    int n_Electrons;
    int unocc_orb;
    int iter;
    int Speed_Elec, Speed_Occ;
    Hartree_Fock_Solver *HartFock;
    bool freeze_core;
    int rank, size;
    double****** t3;
    double****** t3_new;
    mat t1, t1_new;
    field<mat> t2, t2_new;

    mat den_ai;
    field<mat> MOs;

    int CCSDT_n;

    // Gogo
    double Calc_Energy(double toler, int numbah);

    double Update_Energy();
    void Fill_T1();
    void Fill_T2();
    void Fill_T3();

    int EqualFunc(int a, int b);

    void AOtoMO();
    mat c;
    mat fs;

    // Intermediates
    void Fill_X6();

    void Fill_X14();
    void Fill_X15();

    void Fill_X7();
    void Fill_X8();
    void Fill_X4();
    void Fill_X1();
    void Fill_X2();
    void Fill_X3();
    void Fill_Tau();


    void Fill_X12();
    void Fill_X13();

    void Fill_W1();
    void Fill_W2();
    void Fill_W3();
    void Fill_W4();
    void Fill_F1();
    void Fill_F2();
    void Fill_F3();
    void Fill_tau3();

    field<mat> X1, X2, X3, X4, X6, X12, X13, X14, X15;
    mat X7, X8;

    field<mat> tau3, W1, W2, W3, W4;
    mat D1, D2, D3;


};

#endif // CCSDT_H
