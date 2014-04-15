#ifndef CCSD_MEMORY_OPTIMIZED_H
#define CCSD_MEMORY_OPTIMIZED_H

#include <iostream>
#define ARMA_NO_DEBUG
#include "../../home/ole/Desktop/include/armadillo"
#include "hartree_fock_solver.h"
#include <time.h>
#include "mpi.h"

using namespace std;
using namespace arma;

class CCSD_Memory_optimized
{
public:
    CCSD_Memory_optimized(int n_N, vec zz, mat rr, string B_s, int n_Elec, int ran, int siz, Hartree_Fock_Solver *Hartfock);

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
    mat fs;
    Hartree_Fock_Solver *HartFock;

    // Hypercompact form
    field<mat> W_1, W_2, W_4, tau3;
    mat D2, D3, D1;
    mat FS_AI, FS_IJ, FS_AB, DEN_AI;
    mat tau1;
    mat Zero_Matrix;

    // Extreme hypercompact form
    field<mat> W_3;

    // Split up integral for much less storage, even if meny more variables, storage reduced by >60%
    field<mat> integ2, integ3, integ4, integ5, integ6, integ7, integ8, integ9, integ10;

    // Energier:
    double E_old, E_new;

    // Generelle CCSD funksjoner
    double CCSD(double toler, bool print_stuff);
    int EqualFunc(int a, int b);
    double Calc_Energy();

    // Intermediates
    void Fill_W1_and_W3();
    void Fill_W2();
    void Fill_W4();
    void Fill_F1();
    void Fill_F2();
    void Fill_F3();
    void Fill_tau();
    void Fill_2D_tau(int a, int b);

    // Amplitudes
    void Fill_t1_new();
    void Fill_t2_new();
    void Map_T_new();
    mat T_1, T_1_new;
    field<mat>  t2;

    // MOs
    field<mat> MO3, MO4, MO6_1, MO6_2, MO8, MO9, MO10, MO2, MOLeftovers;
    cube compact_mo,compact_mo2;
    void Prepear_AOs();
    long Return_Integral_Index(int a, int b, int i, int j);
    mat c;
    //vec Integrals;
    mat Fill_FS(vec eigval);

    mat integ2_2D, integ3_2D, integ4_2D, integ5_2D, integ6_2D, integ7_2D, integ8_2D, integ9_2D, integ10_2D;
    void Fill_integ2_2D(int a, int i);
    void Fill_integ3_2D(int a, int i);
    void Fill_integ4_2D(int a, int i);
    void Fill_integ5_2D(int a, int i);
    void Fill_integ6_2D(int a, int i);
    void Fill_integ7_2D(int a, int i);
    void Fill_integ8_2D(int a, int i);
    void Fill_integ9_2D(int a, int i);
    void Fill_integ10_2D(int a, int i);
    void Fill_integ11_2D(int a, int i);

    void Fill_integ2_2D_even_even(int a, int i);
    void Fill_integ2_2D_even_odd(int a, int i);
    void Fill_integ2_2D_odd_even(int a, int i);
    void Fill_integ2_2D_odd_odd(int a, int i);

    // MPI
    int rank, size;
    void Map_Part1_For_MPI();
    void Distribute_Part1();
    double** Part1_MPI;
    void Map_T2_For_MPI();
    double** T2_MPI; // Contains T2 new
    vec WORK_EACH_NODE;
    vec WORK_EACH_NODE_Part1;
    vec WORK_EACH_NODE_Part2;
    void Fill_W4_MPI();
};

#endif // CCSD_MEMORY_OPTIMIZED_H
