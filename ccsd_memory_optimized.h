#ifndef CCSD_MEMORY_OPTIMIZED_H
#define CCSD_MEMORY_OPTIMIZED_H

#include <iostream>
#define ARMA_NO_DEBUG
#include "../../home/ole/Desktop/include/armadillo"
#include "hartree_fock_solver.h"
#include <time.h>
//#include "mpi.h"
#include "ccsd_non_iterative_part.h"

using namespace std;
using namespace arma;

class CCSD_Memory_optimized
{
public:
    CCSD_Memory_optimized(int n_N, vec zz, mat rr, string B_s, int n_Elec, int ran, int siz, Hartree_Fock_Solver *Hartfock, bool frez);

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
    field<mat> t2;

    // MOs
    field<mat> MO3, MO4, MO6_1, MO6_2, MO8, MO9, MO10, MO2, MOLeftovers;
    //cube compact_mo2;
    void Prepear_AOs(int nr_freeze);
    long Return_Integral_Index(int a, int b, int i, int j);
    mat c;
    //vec Integrals;
    mat Fill_FS(vec eigval, int nr_freeze);

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

    // Memory distribution arrays
    field<mat> W_5;
    void Fill_W5();

    // MPI
    int rank, size;
    void Distribute_Part1();
    double* MY_OWN_MPI;
    double* SHARED_INFO_MPI;
    vec WORK_EACH_NODE;
    vec WORK_EACH_NODE_Part1;
    vec WORK_EACH_NODE_Part2;
    void Fill_W4_MPI();
    mat MO_Grid_A_B(int frez);
    mat MO_Grid_I_J(int frez);
    void Where_To_Start_On_What();
    mat Start_Pos;
    //MPI_Comm* MPI_WORLD_GROUPS_1;
    //MPI_Comm* MPI_WORLD_GROUPS_2;

    MPI_Datatype* mpi_types_array;

    int** Global_Displacement_1;
    int** Global_Displacement_2;
    int** Global_Worksize_1;
    int** Global_Worksize_2;

    int** Global_Displacement_1_1;
    int** Global_Worksize_1_1;
    int** Global_Displacement_2_1;
    int** Global_Worksize_2_1;

    double**** Index_Swapping_W_4;
    int** What_I_Goes_Where;
    int Nr_Parallel_Operations;
    int Double_Size;
    double* Copy_Matrix;
    int security_size;
    int Local_Displacement1;

    int* Displacement_Each_Node_T2_Parallel;
    int* Work_Each_Node_T2_Parallel;

    int* Displacement_Each_Node_part1_Parallel;
    int* Work_Each_Node_part1_Parallel;

    // Frozen Core
    bool freeze_core;

    // Distribution of work stuff
    int Calc_sum_a_n(int a);
    mat Where_To_Start_Part2;
    int jump;

    int Freeze_How_Meny_Orbitals();
};

#endif // CCSD_MEMORY_OPTIMIZED_H
