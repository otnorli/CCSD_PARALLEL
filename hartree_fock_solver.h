#ifndef HARTREE_FOCK_SOLVER_H
#define HARTREE_FOCK_SOLVER_H

#include <iostream>
#define ARMA_NO_DEBUG
#include "../../home/ole/Desktop/include/armadillo"
#include "matrix_size_setter.h"
#include "fill_alpha.h"
#include "hartree_integrals.h"
#include "mpi.h"

using namespace std;
using namespace arma;

class Hartree_Fock_Solver
{
public:
    // Input parametre
    Hartree_Fock_Solver(int n_N, vec zz, mat rr, string B_s, int n_elec, bool pstf, int ran, int siz);
    int n_Nuclei, steppp;
    vec Z;
    mat R;
    string Basis_Set;
    bool print_stuff;

    // MPI
    int rank, size;
    double** Fock_MPI;
    void Map_MPI();
    vec WORK_PER_NODE;
    double Calc_Integrals_On_The_Fly(int orb1, int orb3, int orb2, int orb4);
    int Calc_Which_Atom_We_Are_Dealing_With(int orb1);
    mat E_index;

    // Globale parametre
    vec n_Basis;
    int n_Electrons;
    mat alpha;
    mat c;
    vec Number_Of_Orbitals;
    int Matrix_Size;
    mat EK, O, D, V, F, eigenvektors_O, eigenvektors_F, P,C;
    vec eigenvalues_O, eigenvalues_F;
    mat Potenser;
    vector<cube> return_Q;
    field<mat> field_Q;
    Hartree_Integrals HartInt;

    mat Return_Field_Q(int i, int j);

    mat ReturnC();
    mat ReturnAlpha();
    int ReturnMatrixSize();
    mat ReturnPotenser();
    vec ReturnNumberOfOrbitals();
    vector<cube> ReturnQ();
    vec Return_Indexed_Q();
    mat return_eigval_F();
    mat returnC();
    mat Energy_Fock_Matrix;

    int Get_Integral_Index(int a, int b, int i, int j);
    vec Stored_Indexed_Q;

    void Normalize_small_c();


    double Calc_Density();
    double Calc_Energy();
    double Make_Fock_Matrix();

    double Single_E_Energy;
    double Two_E_Energy;
    double Nuclei_Rep_Energy;

    double get_Energy(double toler);
    vec Normalize(vec Vektorn, mat Matrisen);
    void Set_New_R(mat rr);

    // DIIS
    void DIIS();
    void Initialize_DIIS();
    mat DIIS_B, delta_p;
    vec DIIS_Z, DIIS_c;
    int number_elements_DIIS;
    vector<mat> Stored_Error;
    vector<mat> Stored_F;

    void Delete_Everything();

    // Spin unrestricted function
    void Unrestricted_Fock_Matrix();
    void Unrestricted_P_Matrix();
    void Unrestricted_Energy();
    mat F_up, F_down;
    mat P_up, P_down;



};

#endif // HARTREE_FOCK_SOLVER_H
