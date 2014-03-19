#include "initializer.h"

Initializer::Initializer(int n_N, int n_E, string B_S, string met, bool r_pos, mat rr, vec zz, double con_crit, string r_uns)
{
    Basis_Set = B_S;
    n_Nuclei = n_N;
    n_Electrons = n_E;
    Relax_Pos = r_pos;
    Method = met;
    R = rr;
    Z = zz;
    convergance_criteria = con_crit;
}

double Initializer::Go()
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


    cout << "Method used: " << Method << endl;
    if (Relax_Pos)
    {
        cout << "Relaxation activated" << endl;
    }
    else
    {
        cout << "Relaxation not activated" << endl;
    }

    cout << "Basis set used: " << Basis_Set << endl;
    cout << "Number of electrons: " << n_Electrons << endl;
    cout << "R = " << endl << R << endl << "Z = " << endl << Z << endl;
    cout << "Starting calculations..." << endl;

    // Lets go

        if (Relax_Pos)
        {
            if (Method == "HF")
            {
                Relax_The_Structure RTS(n_Nuclei, n_Electrons, R, Basis_Set, Z, Method);
                R = RTS.Relaxation_HF(convergance_criteria);
                cout << "After relaxation, new R = " << endl << R << endl;
            }

            if (Method == "CCSD")
            {
                Relax_The_Structure RTS(n_Nuclei, n_Electrons, R, Basis_Set, Z, Method);
                R = RTS.Relaxation_CCSD(convergance_criteria);
                cout << "After relaxation, new R = " << endl << R << endl;
            }
        }

        if (Method == "HF")
        {
            // Kaller på Hartree Fock
            Hartree_Fock_Solver HartFock(n_Nuclei, Z, R, Basis_Set, n_Electrons, true);
            Energy = HartFock.get_Energy(convergance_criteria);
            cout << "Energy = " << Energy << " med " << convergance_criteria << " som convergens kriterie i Hartree Fock Metoden" << endl;
        }

        if (Method == "CCSD")
        {
            // Kaller på Coupled Cluster, CCSD
            //ccsd_v2_optimized CC(n_Nuclei, Z, R, Basis_Set, n_Electrons);
            CCSD_Memory_optimized CC(n_Nuclei, Z, R, Basis_Set, n_Electrons);
            Energy = CC.CCSD(convergance_criteria, true);
            cout << "Energy = " << Energy << " med CCSD metoden" << endl;
        }

        return Energy;

}
