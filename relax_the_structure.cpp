#include "relax_the_structure.h"

Relax_The_Structure::Relax_The_Structure(int n_N, int n_elec, mat rr, string B_s, vec zz, string meth)
{
    // This class should put the system in a local minimum in terms of positions for the nuclei

    n_Nuclei = n_N;
    Z = zz;
    R = rr;
    Basis_Set = B_s;
    n_Electrons = n_elec;
    Method = meth;
}

mat Relax_The_Structure::Relaxation_HF(double toler)
{
    double force_tolerance = 0.1; // NOT SAME AS toler VARIABLE, THIS IS NOT REALISTICALY 10e-10!!!!
    double max_force = 1000000;
    double dr,dt;
    double e_new, e_old;
    double iter = 0.0;
    int max_iter = 5;
    double constant_term = 0.01;

    // Initialise system
    Hartree_Fock_Solver Object(n_Nuclei, Z, R, Basis_Set, n_Electrons, false);
    Forces_R = zeros(R.size()/3, 3);
    dr = 0.1;
    dt = 0.01;

    while ((max_force > force_tolerance) && (iter < max_iter))
    {
        iter += 1.0;
        Object.Set_New_R(R);
        e_old = Object.get_Energy(toler);

        for (int i = 0; i < n_Nuclei; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                R(i,j) += dr;
                Object.Set_New_R(R);
                e_new = Object.get_Energy(toler);
                R(i,j) -= 2*dr;
                Object.Set_New_R(R);
                e_new += Object.get_Energy(toler);
                R(i,j) += dr;
                Forces_R(i,j) = constant_term*(e_new - 2*e_old)/((iter)*dr*dr);
            }
        }

        cout << sqrt(R(0)*R(0) + R(1)*R(1) + R(2)* R(2)) << " " << e_old << endl;

        return_R = R + Forces_R*dt;
        max_force = max(max(Forces_R))*dt;

        R = return_R;
    }

    return return_R;
}

mat Relax_The_Structure::Relaxation_CCSD(double toler)
{
    Forces_R = zeros(R.size()/2, 3);
    double force_tolerance = toler;
    double max_force = 1000000;
    double dr,dt;
    double e_new, e_old;
    double iter = 0.0;
    int max_iter = 5;
    double constant_term = 0.01;

    // Do some calculations here

    return Forces_R;
}
