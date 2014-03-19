#include <iostream>
#include "../../home/ole/Desktop/include/armadillo"
#include <initializer.h>
#include <time.h>
#include "../../home/ole/Desktop/OpenBLAS/cblas.h"
#include "../../home/ole/Desktop/OpenBLAS/common.h"

using namespace arma;
using namespace std;

/*
  This is a program applying the Hartree-Fock method to determine an approximation to the wavefunction and the ground state energy.
  This program is also capable of calculating CCSD energies.
  This program also has DIIS implemented for speedup and this is always used (unless you change the code somewhat).

  Assumptions:
  - Approximate the wavefunction by a single slater determinant for fermions
  - Born-Opphenheimer
  - Mean field approximation
  - Relativistic effects neglected
  - The variational solution is assumed to be a linaer combination of a finite number of basis functions, which is usually chosen to be orthogonal
  - The combination of gaussian which is used is supposed to in combination resembel a Slater orbital.

  Programming language:
  - C++
  - C
  - Library package Armadillo

  Comment language:
  - Mostly english
  - Some native norwegian

  Inner workings:
  - Filling up alpha automaticly, or manually if you want, add it manually in fill_alpha.cpp. Just do the same as is done
  - Implemented meny basis sets, should be reasonably easy to understand how it works from the if tests
  - See hartree_fock_solver.cpp for the algorithm and hartree_integrals.cpp for the mathematics
  - coupled_cluster.cpp for algorithm and coupled_cluster_integrals.cpp for some rearranging for speedups

  Brukermanual:
  - Dersom du har en krisesituasjon gå inn i klassen fill_alpha og oppdater max_size_bas til et større tall, skal automatisere i en fremtidig versjon
  - Vennligst send eventuelle klagebrev til nærmeste søppelkasse,  adresse Dont Careveien 45
  - Ta Ctrl + F og let etter (<-- !) for å finne viktige kommentarer underveis i koden

  Input to give:
  - Method
  - Basis set
  - Nuclei charges and positions in a.u.
  - Relaxation on or off

  Current bugs:
  - Boys_Start N undetermined
  - DIIS buggy :O
  - Some of the D, F, ... orbitals in fill_alpha can be skipped because of symmetry?

  Program written by Ole Norli
*/

int main()
{
    // old    ,  new   ,   hf                                                    old-hf     new-hf      (new-hf)/(old-hf)
    //   0.99s ,   0.92s ,   0.81s , 19 iterasions ,  26 orbitals , 4-31G         0.18       0.11           0,61
    //   4.39s ,   4.00s ,   3.46s , 18 iterations ,  38 orbitals , dzvp          0,93       0,54           0,58
    //  16.82s ,  15.54s ,  13.47s , 18 iterations ,  50 orbitals , cc-pVDZ       3,35       2.06           0,61
    //  15.01s ,  12.05s ,   7.58s , 20 iterations ,  62 orbitals , 6-311ss       7,43       4.47           0,60
    //  76.69s ,  58.38s ,  30.33s , 18 iterations ,  98 orbitals , 6-311-2d2p   46,36      28.05           0,60
    // 367.33s , 295.75s , 196.15s , 18 iterations , 130 orbitals , cc-pVTZ     171,18      99.60           0,58

    // Definerer noen variabler, ikke bry deg om at de er veldig rare:
    int n_Nuclei=5000, n_Electrons;bool Test_Ongoing, Relax_Pos;
    string Basis_Set, Method;double Energy, convergance_criteria;
    vec Z = zeros(n_Nuclei);mat R = zeros(n_Nuclei, 3);
    clock_t start = clock(); string R_Units;
    clock_t slutt;

    // Litt input her sånn generelt
    convergance_criteria = pow(10.0, -8.0); // Criteria at which calculations broken
    Test_Ongoing = false;

    // Velg basis set

    //Basis_Set = "STO-3G";
    Basis_Set = "cc-pVTZ";
    //Basis_Set = "4-31G";
    //Basis_Set = "3-21G";
    //Basis_Set = "cc-pVDZ";
    //Basis_Set = "Thjiisen";
    //Basis_Set = "6-311ss";
    //Basis_Set = "6-311-2d2p";
    //Basis_Set = "dzvp";

    // Velg metode
    Method = "CCSD";
    //Method = "CCSDt";
    //Method = "HF";

    // Velg ekstras
    Relax_Pos = false; // This option dont work
    //R_Units = "Angstrom"; // If you give units in Angstrom activate this. If not then dont

    // Input posisjon og ladninger:

    //Z(0) = 10;

    Z(0) = 1;
    Z(1) = 8;
    Z(2) = 1;

    R(0,0) = 0;
    R(0,1) = 1.079252144093028;
    R(0,2) = 1.474611055780858;
    R(2,0) = 0;
    R(2,1) = 1.079252144093028;
    R(2,2) = -1.474611055780858;
    R(1,0) = 0;
    R(1,1) = 0;
    R(1,2) = 0;

    if (Test_Ongoing == false)
    {
        // Lets go :D
        Initializer Init(n_Nuclei, n_Electrons, Basis_Set, Method, Relax_Pos, R, Z, convergance_criteria, R_Units);
        Energy = Init.Go();
    }

    if (Test_Ongoing == true)
    {
        // Lets not go.. :(
        cout << "Test activated" << endl;

        mat A = randu(4,6);
        cout << A << endl;
        mat B = zeros(4,3);

        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 6; j++)
            {
                B(i/2, j/2) = A(i,j);
                j++;
            }
            i++;
        }

        for (int i = 1; i < 4; i++)
        {
            for (int j = 0; j < 6; j++)
            {
                B(i/2+2, j/2) = A(i,j);
            }
        }

        cout << B << endl;






















    }

    // Finalize
    slutt = clock();
    cout << "Tid brukt: " << (double)(slutt - start)/CLOCKS_PER_SEC << "s" << endl;
    return 0;
}
