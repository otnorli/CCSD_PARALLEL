#ifndef HARTREE_INTEGRALS_H
#define HARTREE_INTEGRALS_H

#define ARMA_NO_DEBUG
#include "../../home/ole/Desktop/include/armadillo"
#include <iostream>
#include <iostream>
#include <iomanip>

using namespace std;
using namespace arma;

class Hartree_Integrals
{
public:
    Hartree_Integrals();
    int n_Nuclei, Matrix_Size, n_Electrons;
    vec n_Basis, n_Orbitals, Z;
    mat alpha, c, R, Potenser, Ek, O;
    void Set_Input(mat alp, int n_tot_bf, mat cc, vec n_O, mat pot, vec n_B, int n_N, mat pos, vec zZz, int n_E);

    // Overlap
    mat Overlap_Matrix();
    double Overlap_Integral_Single(int ind1, int ind2, int bas1, int bas2, int E_counter);
    mat E_index;
    mat get_E_index();
    vector<cube> E_ij;
    void Fill_E_ij();

    // Kinetic energy
    mat Kinetic_Energy();
    double Kinetic_Energy_Single(int ind1, int ind2, int bas1, int bas2, int E_counter);

    // Nuclei-Electron interaction
    mat Nuclei_Electron_Interaction();
    double Nuclei_Electron_Interaction_Single(int ind1, int ind2, int a1, int a2, int bas1, int bas2, int E_counter);
    double Nuclei_Electron_Interaction_Single_1d(double p, int t, int u, int v, rowvec R1, rowvec R2, double n);
    double Boys(double x, double n);
    double A_function(double n, double input);
    double Gamma_Function(double n, double input);
    double Boys_Start;
    int Boys_N;
    vec F_Boys;
    int Double_Factorial(int TALL);
    void Set_Boys_Start(int N);
    int Factorial(int N);

    // Boys stuff
    mat F_tabulated;

    // Nuclei-Nuclei repulsion
    double Nuclei_Nuclei_Interaction();

    // Electron-Electron interaction
    double Electron_Electron_Interaction_Single(int ind1, int ind2, int ind3, int ind4, int a1, int a2, int a3, int a4, int bas1, int bas2, int bas3, int bas4, int E_counter1, int E_counter2);
    vector<cube> R_ijk;
    void Set_R_ijk(double p, int t, int u, int v, rowvec R1, rowvec R2);

    // Cleanup
    void Delete_Everything();
    void Fill_F_Tabulated();
};

#endif // HARTREE_INTEGRALS_H
