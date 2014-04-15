#include "initializer.h"

Initializer::Initializer(int n_N, int n_E, string B_S, string met, bool r_pos, mat rr, vec zz, double con_crit, int ran, int siz, bool us_ang)
{
    Basis_Set = B_S;
    n_Nuclei = n_N;
    n_Electrons = n_E;
    Relax_Pos = r_pos;
    Method = met;
    R = rr;
    Z = zz;
    convergance_criteria = con_crit;
    rank = ran;
    size = siz;
    use_angstrom = us_ang;
}

double Initializer::Go(bool printie)
{
    double Energy;
    int temp = 0;
    for (int i = 0; i < n_Nuclei; i++)
    {
        if (Z(i) > 0)
        {
            temp += 1;
        }
    }

    if (rank != 0)
    {
        printie = false;
    }

    n_Nuclei = temp;
    vec Z_temp = Z;
    mat R_temp = R;
    Z = zeros(n_Nuclei);
    R = zeros(n_Nuclei,3);

    n_Electrons = 0;
    for (int i = 0; i < n_Nuclei; i++)
    {
        Z(i) = Z_temp(i);
        R.row(i) = R_temp.row(i);
        n_Electrons += Z(i);
    }

    if (use_angstrom == true)
    {
        R = R * 1.889725989;
    }

    if (printie==true)
    {
        cout << "Method used: " << Method << endl;
    }

    if (Relax_Pos)
    {
        if (printie == true)
        {
            cout << "Relaxation activated" << endl;
        }
    }

    else
    {
        if (printie == true)
        {
            cout << "Relaxation not activated" << endl;
        }
    }

    if (printie == true)
    {
        cout << "Basis set used: " << Basis_Set << endl;
        cout << "Number of electrons: " << n_Electrons << endl;
        cout << "R = " << endl << R << endl << "Z = " << endl << Z << endl;
        cout << "Starting calculations..." << endl;
    }

    // Lets go

    Hartree_Fock_Solver HartFock(n_Nuclei, Z, R, Basis_Set, n_Electrons, printie, rank, size);
        if (Method == "HF")
        {
            // Kaller på Hartree Fock

            Energy = HartFock.get_Energy(convergance_criteria);
            if (printie == true)
            {
                cout << "Energy = " << Energy << " med " << convergance_criteria << " som convergens kriterie i Hartree Fock Metoden" << endl;
            }
        }

        if (Method == "CCSD")
        {
            // Kaller på Coupled Cluster, CCSD
            CCSD_Memory_optimized CC(n_Nuclei, Z, R, Basis_Set, n_Electrons, rank, size, &HartFock);
            Energy = CC.CCSD(convergance_criteria, printie);
            if (printie == true)
            {
                cout << "Energy = " << Energy << " med CCSD metoden" << endl;
            }
        }

        return Energy;

}
