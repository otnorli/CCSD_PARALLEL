#include "initializer.h"

Initializer::Initializer(int n_N, int n_E, string B_S, string met, bool r_pos, mat rr, vec zz, double con_crit, int ran, int siz, bool us_ang, bool frez)
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
    frozen_core = frez;
}

double Initializer::Go(bool printie)
{
    double Energy;
    int half_elec;
    double temp_en;
    int temp = 0;

    // No need to give input how many atoms are present.
    // We use this code to ensure number of nuclei is the correct number here
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

    // If angstrom is used, scale it to atomic units
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

    // Check which CCSDT-n method is in use
    int method_nr;

    if (Method == "CCSDT")
    {
        method_nr = 6;
    }

    if (Method == "CCSD1")
    {
        method_nr = -1;
        Method = "CCSDT";
    }

    if (Method == "CCSDT-1a")
    {
        method_nr = 1;
        Method = "CCSDT";
    }

    if (Method == "CCSDT-1b")
    {
        method_nr = 2;
        Method = "CCSDT";
    }

    if (Method == "CCSDT-2")
    {
        Method = "CCSDT";
        method_nr = 3;
    }

    if (Method == "CCSDT-3")
    {
        method_nr = 4;
        Method = "CCSDT";
    }

    if (Method == "CCSDT-4")
    {
        method_nr = 5;
        Method = "CCSDT";
    }

    if (Method == "CCSDT--Q")
    {
        cout << "Not implemented" << endl;
        //method_nr = 7;
        Method = "CCSDT";
    }

    if (Method == "UCCSDT")
    {
        // Not working jet
        Method = "CCSDT";
        method_nr = 7;
    }


    // Lets go. Call correct function and start calculations

    Hartree_Fock_Solver HartFock(n_Nuclei, Z, R, Basis_Set, n_Electrons, printie, rank, size, frozen_core);
    if (Method == "HF")
    {
        // Kaller på Hartree Fock

        Energy = HartFock.get_Energy(convergance_criteria, 0);
        if (printie == true)
        {
            cout << "Energy = " << Energy << " med " << convergance_criteria << " som convergens kriterie i Hartree Fock Metoden" << endl;
        }
    }

    if (Method == "UHF")
    {
        half_elec = (int)n_Electrons/2;

        if (n_Electrons % 2 == 1)
        {
            HartFock.Set_UP_DOWN(half_elec+1, half_elec);
            Energy = HartFock.get_Energy(convergance_criteria, 1);
            if (printie == true)
            {
                cout << "Energy = " << Energy << " med " << convergance_criteria << " som convergens kriterie med UHF metoden" << endl;
            }
        }

        else
        {
            HartFock.Set_UP_DOWN(half_elec, half_elec);
            temp_en = HartFock.get_Energy(convergance_criteria, 1);
            HartFock.Set_UP_DOWN(half_elec + 1, half_elec - 1);
            Energy = HartFock.get_Energy(convergance_criteria, 1);
            if (printie == true)
            {
                cout << "Energies calculated with UHF" << endl;
                cout << "Energy = " << Energy << " med tripplet configuration. Energy = " << temp_en << " med singlet configuration." << endl;
            }
        }
    }

    if (Method == "CCSD")
    {
        // Kaller på Coupled Cluster, CCSD
        CCSD_Memory_optimized CC(n_Nuclei, Z, R, Basis_Set, n_Electrons, rank, size, &HartFock, frozen_core);
        Energy = CC.CCSD(convergance_criteria, printie);
        if (printie == true)
        {
            cout << "Energy = " << Energy << " med CCSD metoden" << endl;
        }
    }

    if (Method == "CCSDT")
    {
        // Call CCSDT
        CCSDT CC(n_Nuclei, Z, R, Basis_Set, n_Electrons, rank, size, &HartFock, frozen_core);
        Energy = CC.Calc_Energy(convergance_criteria, method_nr);
        if (printie == true)
        {
            cout << "Energy = " << Energy << " med CCSDT methoden" << endl;
        }
    }

    return Energy;

}
