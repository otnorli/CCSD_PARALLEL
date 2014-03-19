#include "coupled_cluster.h"

coupled_cluster::coupled_cluster(int n_N, vec zz, mat rr, string B_s, int n_Elec)
{
    n_Nuclei = n_N;
    Z = zz;
    R = rr;
    Basis_Set = B_s;
    n_Electrons = n_Elec;
}

double coupled_cluster::CCSD(double toler, bool print_stuff=true)
{
    return 0;
}
    /*
    bool continue_ccsd = true;
    double convergance_check;
    int Itermax = 100;
    iter = 0;
    E_old = 0;
    E_new = 0;

    // Initializion
    double E_HF;
    Hartree_Fock_Solver HartFock(n_Nuclei, Z, R, Basis_Set, n_Electrons, true);
    E_HF = HartFock.get_Energy(toler); // Calc hartree fock energy
    Matrix_Size = 2*HartFock.ReturnMatrixSize(); // This will be twice as big as in Hatree Fock, because of spin up and down
    unocc_orb = Matrix_Size - n_Electrons; // Number of unocupied orbitals
    Initialize_DIIS();
    t1 = zeros(Matrix_Size, Matrix_Size); // Amplitudes
    fs = zeros(Matrix_Size, Matrix_Size);
    denom_ai = zeros(Matrix_Size, Matrix_Size);
    t1_new = zeros(Matrix_Size, Matrix_Size);
    F1 = zeros(Matrix_Size, Matrix_Size);
    F2 = zeros(Matrix_Size, Matrix_Size);
    F3 = zeros(Matrix_Size, Matrix_Size);
    D1 = zeros(unocc_orb, Matrix_Size);
    D2 = zeros(Matrix_Size, Matrix_Size);
    D3 = zeros(Matrix_Size, Matrix_Size);
    W_4.set_size(unocc_orb, n_Electrons);
    T_2.set_size(unocc_orb, n_Electrons);
    T_2_new.set_size(unocc_orb, n_Electrons);
    mat temp_matrix = zeros(unocc_orb, n_Electrons);

    // Initialize four dimentional arrays
    for (int i = 0; i < Matrix_Size; i++)
    {
        cube temp_mat = zeros(Matrix_Size, Matrix_Size, Matrix_Size);
        //t2.push_back(temp_mat);
        //t2_new.push_back(temp_mat);
        integ.push_back(temp_mat);
        denom_abij.push_back(temp_mat);
        //G1.push_back(temp_mat);
        //G3.push_back(temp_mat);
        //G2.push_back(temp_mat);
        //G4.push_back(temp_mat);
        tau.push_back(temp_mat);
    }

    for (int i = 0; i < n_Electrons; i++)
    {
        cube temp_mat = zeros(n_Electrons, n_Electrons, n_Electrons);
        G2.push_back(temp_mat);
    }

    for (int i = 0; i < unocc_orb; i++)
    {
        cube temp_mat = zeros(n_Electrons, n_Electrons, n_Electrons);
        G4.push_back(temp_mat);
        G1.push_back(temp_mat);
        temp_mat = zeros(n_Electrons, unocc_orb, n_Electrons);
        G3.push_back(temp_mat);
        for (int j = 0; j < n_Electrons; j++)
        {
            W_4(i,j) = temp_matrix;
            T_2(i,j) = temp_matrix;
            T_2_new(i,j) = temp_matrix;
        }
    }

    coupled_cluster_integrals ccint(Matrix_Size/2, HartFock.ReturnC(), n_Electrons);
    //integ = ccint.Return_Integrals(HartFock.ReturnQ()); // 2 elektron integralene her
    fs = ccint.Return_FS(HartFock.return_eigval_F()); // alt annet her

    // Coupled Cluster is an iterative process, meaning we make an initial guess of t1 and t2 and calculate the energy
    // Then update t1 and t2 and check for convergance in the energy

    // Find denominator matrix, this should not change during our calculations
    for (int a = n_Electrons; a < Matrix_Size; a++)
    {
        for (int i = 0; i < n_Electrons; i++)
        {
            denom_ai(a,i) = fs(i,i) - fs(a,a);
            for (int b = n_Electrons; b < Matrix_Size; b++)
            {
                for (int j = 0; j < n_Electrons; j++)
                {
                    denom_abij.at(a)(b,i,j) = fs(i,i) + fs(j,j) - fs(a,a) - fs(b,b);
                }
            }
        }
    }

    // The initial guess of t1 and t2 is made here, t1 is guessed to be 0

    //for (int a = n_Electrons; a < Matrix_Size; a++)
    for (int a = 0; a < unocc_orb; a++)
    {
        //for (int b = n_Electrons; b < Matrix_Size; b++)
        for (int b = 0; b < unocc_orb; b++)
        {
            for (int i = 0; i < n_Electrons; i++)
            {
                for (int j = 0; j < n_Electrons; j++)
                {
                    //t2.at(a)(b,i,j) = integ.at(i)(j,a,b) / denom_abij.at(a)(b,i,j);
                    T_2.at(a,i)(b,j) = integ.at(i)(j,a,b) / denom_abij.at(a)(b,i,j);
                }
            }
        }
    }

    // Here starts the iteration, most calculations are taken into functions, even tho they are only called at one place
    // This is to get a more compact algo which is in my opinion more easy to read

    E_new = Calc_Energy(); // Starting energy, with t1=0 guess

    if (print_stuff == true)
    {
        cout << "Energi: " << E_new << " Steg: " << iter << endl;
    }

    while (continue_ccsd == true)
    {
        E_old = E_new;

        // Update intermediates
        Fill_tau();

        Fill_G4();
        Fill_G1();
        Fill_G2();
        Fill_G3();

        Fill_F1();
        Fill_F2();
        Fill_F3();

        // Find new amplitudes
        Fill_t1_new();
        Fill_t2_new();

        // Do DIIS
        if (DIIS_does == true)
        {
            //DIIS(1);
            //DIIS(2);
        }

        // Update amplitudes
        t1 = t1_new;
        //t2 = t2_new;
        T_2 = T_2_new;

        // Find energy
        E_new = Calc_Energy();
        iter += 1;

        // Output energy
        if (print_stuff == true)
        {
            cout << "Energi: " << E_new << " Steg: " << iter << endl;
        }

        // Check convergance
        convergance_check = sqrt((E_new - E_old) * (E_new - E_old));
        if (convergance_check < toler || iter >= Itermax)
        {
            continue_ccsd = false;
        }
    }

    return E_new+E_HF;
}

double coupled_cluster::CCSDt(double toler)
{
    cout << "Method not implemented << endl;" << endl;
    return 0;
}

double coupled_cluster::tau1(int a, int b, int i, int j)
{
    double value;
    //value = t2.at(a)(b,i,j) + 0.5 * (t1(a,i) * t1(b,j) - t1(b,i)*t1(a,j));
    value = T_2.at(a,i).at(b,j) + 0.5 * (t1(a,i) * t1(b,j) - t1(b,i)*t1(a,j));
    return value;
}

double coupled_cluster::tau2(int a, int b, int i, int j)
{
    double value;
    //value = t2.at(a)(b,i,j) + t1(a,i)*t1(b,j) - t1(a,j)*t1(b,i);
    value = T_2.at(a,i).at(b,j) + t1(a,i)*t1(b,j) - t1(a,j)*t1(b,i);
    return value;
}

double coupled_cluster::EqualFunc(int a, int b)
{
    if (a == b)
    {
        return 1;
    }

    else
    {
        return 0;
    }
}

void coupled_cluster::Fill_G1()
{
    // [W_2]

    /*
    G1 = integ;
    for (int a = n_Electrons; a < Matrix_Size; a++)
    {
        for (int b = n_Electrons; b < Matrix_Size; b++)
        {
            for (int g = n_Electrons; g < Matrix_Size; g++)
            {
                for (int h = n_Electrons; h < Matrix_Size; h++)
                {
                    for (int i = 0; i < n_Electrons; i++)
                    {
                        G1.at(a)(b,g,h) -= t1.at(b,i) * integ.at(a)(i,g,h);
                        G1.at(a)(b,g,h) += t1.at(a,i) * integ.at(b)(i,g,h);

                        for (int j = 0; j < n_Electrons; j++)
                        {
                            G1.at(a)(b,g,h) += 0.25*tau2(a,b,i,j)*integ.at(i)(j,g,h);
                        }
                    }
                }
            }
        }
    }

    double temp;
    int temp_a;

    //for (int a = n_Electrons; a < Matrix_Size; a++)
    for (int a = 0; a < unocc_orb; a++)
    {
        temp_a = a + n_Electrons;
        for (int i = 0; i < n_Electrons; i++)
        {
            for (int j = 0; j < n_Electrons; j++)
            {
                for (int m = 0; m < n_Electrons; m++)
                {
                    G1.at(a)(i,j,m) = integ.at(temp_a)(m,i,j);
                    for (int e = n_Electrons; e < Matrix_Size; e++)
                    {
                        G1.at(a)(i,j,m) += integ.at(temp_a)(m,i,e) * t1(e,j);
                        G1.at(a)(i,j,m) -= integ.at(temp_a)(m,j,e) * t1(e,i);
                        temp = 0;
                        for (int f = n_Electrons; f < Matrix_Size; f++)
                        {
                            temp += integ.at(temp_a)(m,e,f) * tau.at(e)(f,i,j);
                        }
                        G1.at(a)(i,j,m) += 0.5*temp;
                    }
                }
            }
        }
    }
}

void coupled_cluster::Fill_F1()
{
    for (int a = n_Electrons; a < Matrix_Size; a++)
    {
        for (int b = n_Electrons; b < Matrix_Size; b++)
        {
            F1.at(a,b) = (1 - EqualFunc(a,b))*fs.at(a,b);

            for (int i = 0; i < n_Electrons; i++)
            {
                F1.at(a,b) -= 0.5*fs.at(i,b)*t1.at(a,i);

                for (int g = n_Electrons; g < Matrix_Size; g++)
                {
                    F1.at(a,b) += integ.at(i)(a,g,b) * t1.at(g,i);

                    for (int j = 0; j < n_Electrons; j++)
                    {
                        F1.at(a,b) -= 0.5*tau1(a,g,i,j) * integ.at(i)(j,b,g);
                    }
                }
            }
        }
    }

    int temp_a;
    for (int a = 0; a < unocc_orb; a++)
    {
        temp_a = a + n_Electrons;
        for (int m = 0; m < n_Electrons; m++)
        {
            D1(a,m) = fs(temp_a,m);
            for (int f = n_Electrons; f < Matrix_Size; f++)
            {
                for (int n = 0; n < n_Electrons; n++)
                {
                    D1(a,m) += integ.at(m)(n,temp_a,f) * t1(f,n);
                }
            }
        }
    }
}

void coupled_cluster::Fill_G3()
{
    // [W_4]

    /*
    G3 = integ;

    for (int i = 0; i < n_Electrons; i++)
    {
        for (int j = 0; j < n_Electrons; j++)
        {
            for (int g = n_Electrons; g<Matrix_Size; g++)
            {
                G3.at(i)(span(), span(), span(j)) += t1(g,j)*integ.at(i)(span(), span(), span(g));
            }
        }
    }

    for (int i = 0; i < n_Electrons; i++)
    {
        for (int j = 0; j < n_Electrons; j++)
        {
            for (int a = n_Electrons; a < Matrix_Size; a++)
            {
                for (int b = n_Electrons; b < Matrix_Size; b++)
                {
                    for (int k = 0; k < n_Electrons; k++)
                    {
                        G3.at(i)(a,b,j) -= t1.at(a,k) * integ.at(i)(k,b,j);

                        for (int g = n_Electrons; g < Matrix_Size; g++)
                        {
                            G3.at(i)(a,b,j) -= integ.at(i)(k,b,g) * (0.5*t2.at(g)(a,j,k) + t1.at(g,j) * t1.at(a,k));
                        }
                    }
                }
            }
        }
    }

    double temp;
    int temp_a, temp_e;
    //for (int a = n_Electrons; a < Matrix_Size; a++)
    for (int a = 0; a < unocc_orb; a++)
    {
        temp_a = a + n_Electrons;
        //for (int e = n_Electrons; e < Matrix_Size; e++)
        for (int e = 0; e < unocc_orb; e++)
        {
            temp_e = e + n_Electrons;
            for (int i = 0; i < n_Electrons; i++)
            {
                for (int m = 0; m < n_Electrons; m++)
                {
                    // <kb||cj> t_ik^ac
                    //G3.at(a)(e,i,m) = integ.at(temp_a)(m,i,temp_e);
                    //G3.at(a)(i,e,m) = integ.at(temp_a)(m, i, temp_e);
                    W_4.at(a,i).at(e,m) = integ.at(temp_a)(m,i,temp_e);

                    for (int n = 0; n < n_Electrons; n++)
                    {
                        //G3.at(a)(i,e,m) -= G4.at(e)(n,m,i) * t1(temp_a,n);
                        W_4.at(a,i).at(e,m) -= G4.at(e)(n,m,i) * t1(temp_a, n);
                    }

                    for (int f = n_Electrons; f < Matrix_Size; f++)
                    {
                        //G3.at(a)(i,e,m) += integ.at(m)(temp_a,temp_e,f) * t1(f,i);
                        W_4.at(a,i).at(e,m) += integ.at(m)(temp_a,temp_e,f) * t1(f,i);
                        temp = 0;
                        for (int n = 0; n < n_Electrons; n++)
                        {
                            //temp += integ.at(m)(n,temp_e,f) * t2.at(temp_a)(f,i,n);
                            temp += integ.at(m)(n,temp_e,f) * T_2.at(a, i).at(f-n_Electrons, n);
                        }
                        //G3.at(a)(i,e,m) += 0.5*temp;
                        W_4.at(a,i).at(e,m) += 0.5*temp;
                    }
                }
            }
        }
    }
}

void coupled_cluster::Fill_F3()
{
    F3 = fs;
    for (int i = 0; i < n_Electrons; i++)
    {
        for (int a = n_Electrons; a < Matrix_Size; a++)
        {
            for (int j = 0; j < n_Electrons; j++)
            {
                for (int b = n_Electrons; b < Matrix_Size; b++)
                {
                    F3.at(i,j) += t1.at(b,j)*integ.at(i)(j,a,b);
                }
            }
        }
    }

    double temp;
    for (int a = n_Electrons; a < Matrix_Size; a++)
    {
        for (int e = n_Electrons; e < Matrix_Size; e++)
        {
            D3(a,e) = (1 - EqualFunc(a,e)) * fs(a,e);

            for (int m = 0; m < n_Electrons; m++)
            {
                D3(a,e) -= D1(e-n_Electrons,m) * t1(a,m);
                for (int f = n_Electrons; f < Matrix_Size; f++)
                {
                    D3(a,e) += integ.at(a)(m,e,f) * t1(f,m);
                    temp = 0;
                    for (int n = 0; n < n_Electrons; n++)
                    {
                        //temp += integ.at(n)(m,e,f) * t2.at(a)(f,n,m);
                        temp += integ.at(n)(m,e,f) * T_2.at(a-n_Electrons, n).at(f-n_Electrons,m);
                    }
                    D3(a,e) -= 0.5*temp;
                }
            }
        }
    }
}

void coupled_cluster::Fill_G2()
{
    // [W_1]

    /*
    G2 = integ;
    for (int i = 0; i < n_Electrons; i++)
    {
        for (int j = 0; j < n_Electrons; j++)
        {
            for (int k = 0; k < n_Electrons; k++)
            {
                for (int l = 0; l < n_Electrons; l++)
                {
                    for (int a = n_Electrons; a < Matrix_Size; a++)
                    {
                        G2.at(i)(j,k,l) += t1.at(a,l) * integ.at(i)(j,k,a);
                        G2.at(i)(j,k,l) -= t1.at(a,k) * integ.at(i)(j,l,a);

                        for (int b = n_Electrons; b < Matrix_Size; b++)
                        {
                            G2.at(i)(j,k,l) += 0.25*tau2(a,b,k,l) * integ.at(i)(j,a,b);
                        }
                    }
                }
            }
        }
    }

    double temp;
    for (int i = 0; i < n_Electrons; i++)
    {
        for (int j = 0; j < n_Electrons; j++)
        {
            for (int k = 0; k < n_Electrons; k++)
            {
                for (int l = 0; l < n_Electrons; l++)
                {
                    G2.at(i)(j,k,l) = integ.at(k)(l,i,j);
                    for (int e = n_Electrons; e < Matrix_Size; e++)
                    {
                        G2.at(i)(j,k,l) += integ.at(k)(l,i,e) * t1(e,j);
                        G2.at(i)(j,k,l) -= integ.at(k)(l,j,e) * t1(e,i);
                        temp = 0;
                        for (int f = n_Electrons; f < Matrix_Size; f++)
                        {
                            temp += integ.at(k)(l,e,f) * tau.at(e)(f,i,j);
                        }
                        G2.at(i)(j,k,l) += 0.5*temp;
                    }
                }
            }
        }
    }
}

void coupled_cluster::Fill_F2()
{
    for (int i = 0; i < n_Electrons; i++)
    {
        for (int j = 0; j < n_Electrons; j++)
        {
            F2.at(i,j) = (1 - EqualFunc(i,j)) * fs.at(i,j);

            for (int a = n_Electrons; a < Matrix_Size; a++)
            {
                F2.at(i,j) += 0.5*t1.at(a,j) * fs.at(i,a);

                for (int k=0; k < n_Electrons; k++)
                {
                    F2.at(i,j) += t1.at(a,k) * integ.at(i)(k,j,a);

                    for (int g=n_Electrons; g < Matrix_Size; g++)
                    {
                        F2.at(i,j) += 0.5 * integ.at(i)(k,a,g) * tau1(a,g,j,k);
                    }
                }
            }
        }
    }

    double temp;
    for (int m = 0; m < n_Electrons; m++)
    {
        for (int i = 0; i < n_Electrons; i++)
        {
            D2(m,i) = (1 - EqualFunc(m,i)) * fs(m,i);
            for (int e = n_Electrons; e < Matrix_Size; e++)
            {
                D2(m,i) += D1(e - n_Electrons,m) * t1(e,i);
                for (int n = 0; n < n_Electrons; n++)
                {
                    D2(m,i) += integ.at(m)(n,i,e) * t1(e,n);
                    temp = 0;
                    for (int f = n_Electrons; f < Matrix_Size; f++)
                    {
                        //temp += integ.at(m)(n,f,e) * t2.at(f)(e,i,n);
                        temp += integ.at(m)(n,f,e) * T_2.at(f-n_Electrons,i).at(e-n_Electrons,n);
                    }
                    D2(m,i) += 0.5*temp;
                }
            }
        }
    }
}

void coupled_cluster::Fill_G4()
{
    // [W_3]

    for (int m = 0; m < n_Electrons; m++)
    {
        for (int n = 0; n < n_Electrons; n++)
        {
            for (int i = 0; i < n_Electrons; i++)
            {
                //for (int e = n_Electrons; e < Matrix_Size; e++)
                for (int e = 0; e < unocc_orb; e++)
                {
                    G4.at(e)(n,m,i) = integ.at(m)(n,e+n_Electrons,i);
                    for (int f = n_Electrons; f < Matrix_Size; f++)
                    {
                        G4.at(e)(n,m,i) += integ.at(m)(n,e+n_Electrons,f)* t1(f,i);
                    }
                }
            }
        }
    }
}

void coupled_cluster::Fill_t1_new()
{

    t1_new = fs;
    for (int a = n_Electrons; a < Matrix_Size; a++)
    {
        for (int i = 0; i < n_Electrons; i++)
        {
            for (int b = n_Electrons; b < Matrix_Size; b++)
            {
                t1_new.at(a,i) += t1.at(b,i)*F1.at(a,b);
            }

            for (int j = 0; j < n_Electrons; j++)
            {
                t1_new.at(a,i) -= t1.at(a,j) * F2.at(j,i);

                for (int g = n_Electrons; g < Matrix_Size; g++)
                {
                    //t1_new.at(a,i) += F3.at(j,g) * t2.at(a)(g,i,j);
                    t1_new.at(a,i) += F3.at(j,g) * T_2.at(a-n_Electrons,i)(g-n_Electrons,j);
                    t1_new.at(a,i) -= t1.at(g,j)*integ.at(j)(a,i,g);

                    for (int h = n_Electrons; h < Matrix_Size; h++)
                    {
                        //t1_new.at(a,i) -= 0.5*t2.at(g)(h,i,j) * integ.at(j)(a,g,h);
                        t1_new.at(a,i) -= 0.5*T_2.at(g-n_Electrons,i)(h-n_Electrons,j) * integ.at(j)(a,g,h);
                    }

                    for (int k = 0; k < n_Electrons; k++)
                    {
                        //t1_new.at(a,i) -= 0.5*t2.at(a)(g,j,k) * integ.at(k)(j,g,i);
                        t1_new.at(a,i) -= 0.5*T_2.at(a-n_Electrons,j)(g-n_Electrons,k) * integ.at(k)(j,g,i);
                    }
                }
            }
            t1_new.at(a,i) = t1_new.at(a,i)/denom_ai.at(a,i);
        }
    }


    // -0.0501234

/*
    for (int a = n_Electrons; a < Matrix_Size; a++)
    {
        for (int i = 0; i < n_Electrons; i++)
        {
            t1_new(a,i) = fs(a,i);

            for (int c = n_Electrons; c < Matrix_Size; c++)
            {
                t1_new(a,i) += (1 - EqualFunc(a,c)) * fs(a,c) * t1(c,i);
                for (int k = 0; k < n_Electrons; k++)
                {
                    t1_new(a,i) -= fs(c,k) * t1(c,i) * t1(a,k);
                    t1_new(a,i) += D1(c-n_Electrons,k) * t2.at(a)(c,i,k);
                    for (int d = n_Electrons; d < Matrix_Size; d++)
                    {
                        t1_new(a,i) += 0.5*integ.at(k)(a,c,d) * t2.at(c)(d,k,i);
                    }
                }
            }

            for (int k = 0; k < n_Electrons; k++)
            {
                t1_new(a,i) -= t1(a,k) * F2(i,k);
                for (int c = n_Electrons; c < Matrix_Size; c++)
                {
                    for (int l = 0; l < n_Electrons; l++)
                    {
                        t1_new(a,i) -= 0.5 * t2.at(c)(a,k,l) * G4.at(c-n_Electrons)(k,l,i);
                    }

                    t1_new(a,i) += integ.at(k)(a,c,i) * t1(c,k);
                    for (int d = n_Electrons; d < Matrix_Size; d++)
                    {
                        t1_new(a,i) += t1(c,k) * t1(d,i) * integ.at(k)(a,c,d);
                    }
                }
            }

            t1_new(a,i) /= denom_ai(a,i);
        }
    }

}

void coupled_cluster::Fill_t2_new()
{
    /*
    for (int a = n_Electrons; a < Matrix_Size; a++)
    {
        for (int b = n_Electrons; b < Matrix_Size; b++)
        {
            for (int i = 0; i < n_Electrons; i++)
            {
                for (int j = 0; j < n_Electrons; j++)
                {
                    t2_new.at(a)(b,i,j) = integ.at(i)(j,a,b);

                    for (int g = n_Electrons; g < Matrix_Size; g++)
                    {
                        t2_new.at(a)(b,i,j) += t2.at(a)(g,i,j) * F1.at(b,g);
                        t2_new.at(a)(b,i,j) -= t2.at(b)(g,i,j) * F1.at(a,g);

                        t2_new.at(a)(b,i,j) += t1.at(g,i) * integ.at(a)(b,g,j);
                        t2_new.at(a)(b,i,j) -= t1.at(g,j) * integ.at(a)(b,g,i);

                        for (int h = n_Electrons; h < Matrix_Size; h++)
                        {
                            t2_new.at(a)(b,i,j) += 0.5*tau2(g,h,i,j) * G1.at(a)(b,g,h);
                        }

                        for (int k = 0; k < n_Electrons; k++)
                        {
                            t2_new.at(a)(b,i,j) += 0.5 * t2.at(a)(g,i,j) * F3.at(k,g) * (t1.at(a,k) - t1.at(b,k));
                            t2_new.at(a)(b,i,j) += 0.5 * F3.at(k,g) * t2.at(a)(b,i,k) * (t1.at(g,i) - t1.at(g,j));

                            t2_new.at(a)(b,i,j) += t2.at(a)(g,i,k) * G3.at(k)(b,g,j);
                            t2_new.at(a)(b,i,j) += t2.at(b)(g,j,k) * G3.at(k)(a,g,i);
                            t2_new.at(a)(b,i,j) -= t2.at(a)(g,j,k) * G3.at(k)(b,g,i);
                            t2_new.at(a)(b,i,j) -= t2.at(b)(g,i,k) * G3.at(k)(a,g,j);

                            t2_new.at(a)(b,i,j) -= t1.at(g,i) * t1.at(a,k) * integ.at(k)(b,g,j);
                            t2_new.at(a)(b,i,j) -= t1.at(g,j) * t1.at(b,k) * integ.at(k)(a,g,i);
                            t2_new.at(a)(b,i,j) += t1.at(g,j) * t1.at(a,k) * integ.at(k)(b,g,i);
                            t2_new.at(a)(b,i,j) -= t1.at(g,i) * t1.at(b,k) * integ.at(k)(a,g,j);
                        }
                    }

                    for (int k = 0; k < n_Electrons; k++)
                    {
                        t2_new.at(a)(b,i,j) += F2.at(k,i) * t2.at(a)(b,j,k);
                        t2_new.at(a)(b,i,j) -= F2.at(k,j) * t2.at(a)(b,i,k);
                        t2_new.at(a)(b,i,j) += t1.at(b,k) * integ.at(k)(a,i,j);
                        t2_new.at(a)(b,i,j) -= t1.at(a,k) * integ.at(k)(b,i,j);

                        for (int l = 0; l < n_Electrons; l++)
                        {
                            t2_new.at(a)(b,i,j) += 0.5*tau2(a,b,k,l) * G2.at(k)(l,i,j);
                        }
                    }
                    t2_new.at(a)(b,i,j) = t2_new.at(a)(b,i,j) / denom_abij.at(a)(b,i,j);
                }
            }
        }
    }

    double temp;
    int temp_a;

    //for (int a = n_Electrons; a < Matrix_Size; a++)
    for (int a = 0; a < unocc_orb; a++)
    {
        // General case of random a,b,i,j
        //for (int b = n_Electrons; b < Matrix_Size; b++)
        for (int b = 0; b < unocc_orb; b++)
        //for (int b = a+1; b < Matrix_Size; b++)
        {
            for (int i = 0; i < n_Electrons; i++)
            {
                for (int j = 0; j < n_Electrons; j++)
                //for (int j = i+1; j < n_Electrons; j++)
                {
                    // I_ab^ij
                    T_2_new.at(a,i)(b,j) = integ.at(a+n_Electrons)(b,i,j);

                    // P(ab) P(ij) [W_4] t_jk^bc
                    for (int e = n_Electrons; e < Matrix_Size; e++)
                    {
                        for (int m = 0; m < n_Electrons; m++)
                        {
                            T_2_new.at(a,i)(b,j) += T_2.at(e,m)(b,j) * W_4.at(a,i).at(e-n_Electrons,m);
                            T_2_new.at(a,i)(b,j) -= T_2.at(e,m)(b,i) * W_4.at(a,j).at(e-n_Electrons,m);
                            T_2_new.at(a,i)(b,j) -= T_2.at(e,m)(a,j) * W_4.at(b,i).at(e-n_Electrons,m);
                            T_2_new.at(a,i)(b,j) += T_2.at(e,m)(a,i) * W_4.at(b,j).at(e-n_Electrons,m);

                        }
                    }

                    // - P(ab) [W_2] t_k^b
                    for (int m = 0; m < n_Electrons; m++)
                    {
                        T_2_new.at(a,i)(b,j) -= G1.at(a)(i,j,m) * t1(b+n_Electrons,m);
                        T_2_new.at(a,i)(b,j) += G1.at(b)(i,j,m) * t1(a+n_Electrons,m);
                    }

                    // P(ij) I_ab^cj t_i^c + 1/2 I_ab^cd tau_ij^cd
                    for (int e = n_Electrons; e < Matrix_Size; e++)
                    {
                        T_2_new.at(a,i)(b,j) += integ.at(a+n_Electrons)(b+n_Electrons,i,e) * t1(e,j);
                        T_2_new.at(a,i)(b,j) -= integ.at(a+n_Electrons)(b+n_Electrons,j,e) * t1(e,i);

                        temp = 0;
                        for (int f = n_Electrons; f < Matrix_Size; f++)
                        {
                            temp += integ.at(a+n_Electrons)(b+n_Electrons,e,f) * tau.at(e)(f,i,j);
                        }
                        T_2_new.at(a,i)(b,j) += 0.5*temp;
                    }

                    // 1/2 [W_1] tau_kl^ab
                    temp = 0;
                    for (int m = 0; m < n_Electrons; m++)
                    {
                        for (int n = 0; n < n_Electrons; n++)
                        {
                            temp += G2.at(i)(j,m,n) * tau.at(a+n_Electrons)(b+n_Electrons,m,n);
                        }
                    }
                    T_2_new.at(a,i)(b,j) += 0.5*temp;

                    // P(ij) t_jk^ab [F_2]
                    for (int m = 0; m < n_Electrons; m++)
                    {
                        T_2_new.at(a,i)(b,j) += T_2.at(a,j)(b,m) * D2(i,m);
                        T_2_new.at(a,i)(b,j) -= T_2.at(a,i)(b,m) * D2(j,m);
                    }

                    // P(ab) t_ij^bc [F_3]_c^a, akta her. i følge utledninga er det motsatt fortegn...
                    for (int e = n_Electrons; e < Matrix_Size; e++)
                    {
                        T_2_new.at(a,i)(b,j) -= T_2.at(b,i)(e,j) * D3(a,e);
                        T_2_new.at(a,i)(b,j) += T_2.at(a,i)(e,j) * D3(b,e);
                    }

                    T_2_new.at(a,i)(b,j) = T_2_new.at(a,i)(b,j) / denom_abij.at(a+n_Electrons)(b+n_Electrons,i,j);

                    T_2_new.at(b,i)(a,j) = -T_2_new.at(a,i)(b,j);
                    T_2_new.at(b,j)(a,i) = T_2_new.at(a,i)(b,j);
                    T_2_new.at(a,j)(b,i) = T_2_new.at(b,i)(a,j);
                }
            }
        }
    }

/*

        // Special case where a = b, some of the terms will be 0
        for (int i = 0; i < n_Electrons; i++)
        {
            for (int j = i+1; i < n_Electrons; i++)
            {
                // I_ab^ij
                T_2_new.at(a,i)(a,j) = integ.at(a)(a,i,j);

                // P(ij) I_ab^cj t_i^c + 1/2 I_ab^cd tau_ij^cd
                for (int e = n_Electrons; e < Matrix_Size; e++)
                {
                    T_2_new.at(a,i)(a,j) += integ.at(a)(a,i,e) * t1(e,j);
                    T_2_new.at(a,i)(a,j) -= integ.at(a)(a,j,e) * t1(e,i);

                    temp = 0;
                    for (int f = n_Electrons; f < Matrix_Size; f++)
                    {
                        temp += integ.at(a)(a,e,f) * tau.at(e)(f,i,j);
                    }
                    T_2_new.at(a,i)(a,j) += 0.5*temp;
                }

                // 1/2 [W_1] tau_kl^ab
                temp = 0;
                for (int m = 0; m < n_Electrons; m++)
                {
                    for (int n = 0; n < n_Electrons; n++)
                    {
                        temp += G2.at(i)(j,m,n) * tau.at(a)(a,m,n);
                    }
                }
                T_2_new.at(a,i)(a,j) += 0.5*temp;

                // P(ij) t_jk^ab [F_2]
                for (int m = 0; m < n_Electrons; m++)
                {
                    T_2_new.at(a,i)(a,j) += T_2.at(a,j)(a,m) * D2(i,m);
                    T_2_new.at(a,i)(a,j) -= T_2.at(a,i)(a,m) * D2(j,m);
                }

                T_2_new.at(a,i)(a,j) = T_2_new.at(a,i)(a,j) / denom_abij.at(a)(a,i,j);

                T_2_new.at(a,j)(a,i) = -T_2_new.at(a,i)(a,j);
            }
        }

        // Special case where i = j, some of the terms will be 0
        for (int b = a+1; b < Matrix_Size; b++)
        {
            for (int i = 0; i < n_Electrons; i++)
            {
                // I_ab^ij
                T_2_new.at(a,i)(b,i) = integ.at(a)(b,i,i);

                // - P(ab) [W_2] t_k^b
                for (int m = 0; m < n_Electrons; m++)
                {
                    T_2_new.at(a,i)(b,i) -= G1.at(a-n_Electrons)(i,i,m) * t1(b,m);
                    T_2_new.at(a,i)(b,i) += G1.at(b-n_Electrons)(i,i,m) * t1(a,m);
                }

                // 1/2 I_ab^cd tau_ij^cd
                for (int e = n_Electrons; e < Matrix_Size; e++)
                {
                    temp = 0;
                    for (int f = n_Electrons; f < Matrix_Size; f++)
                    {
                        temp += integ.at(a)(b,e,f) * tau.at(e)(f,i,i);
                    }
                    T_2_new.at(a,i)(b,i) += 0.5*temp;
                }

                // 1/2 [W_1] tau_kl^ab
                temp = 0;
                for (int m = 0; m < n_Electrons; m++)
                {
                    for (int n = 0; n < n_Electrons; n++)
                    {
                        temp += G2.at(i)(i,m,n) * tau.at(a)(b,m,n);
                    }
                }
                T_2_new.at(a,i)(b,i) += 0.5*temp;

                // P(ab) t_ij^bc [F_3]_c^a, akta her. i følge utledninga er det motsatt fortegn...
                for (int e = n_Electrons; e < Matrix_Size; e++)
                {
                    T_2_new.at(a,i)(b,i) -= T_2.at(b,i)(e,i) * D3(a,e);
                    T_2_new.at(a,i)(b,i) += T_2.at(a,i)(e,i) * D3(b,e);
                }

                T_2_new.at(a,i)(b,i) = T_2_new.at(a,i)(b,i) / denom_abij.at(a)(b,i,i);

                T_2_new.at(b,i)(a,i) = -T_2_new.at(a,i)(b,i);
            }
        }

        // Special case where a = b AND i = j, meny of the terms will be 0
        for (int i = 0; i < n_Electrons; i++)
        {
            // I_ab^ij
            T_2_new.at(a,i)(a,i) = integ.at(i)(i,a,a);

            // 1/2 I_ab^cd tau_ij^cd
            for (int e = n_Electrons; e < Matrix_Size; e++)
            {
                temp = 0;
                for (int f = n_Electrons; f < Matrix_Size; f++)
                {
                    temp += integ.at(e)(f,a,a) * tau.at(e)(f,i,i);
                }
                T_2_new.at(a,i)(a,i) += 0.5*temp;
            }

            // 1/2 [W_1] tau_kl^ab
            temp = 0;
            for (int m = 0; m < n_Electrons; m++)
            {
                for (int n = 0; n < n_Electrons; n++)
                {
                    temp += G2.at(i)(i,m,n) * tau.at(a)(a,m,n);
                }
            }
            T_2_new.at(a,i)(a,i) += 0.5*temp;

            T_2_new.at(a,i)(a,i) = T_2_new.at(a,i)(a,i) / denom_abij.at(a)(a,i,i);
        }
    }


    /*
    for (int a = n_Electrons; a < Matrix_Size; a++)
    {
        // General case of random a,b,i,j
        //for (int b = n_Electrons; b < Matrix_Size; b++)
        for (int b = a+1; b < Matrix_Size; b++)
        {
            for (int i = 0; i < n_Electrons; i++)
            {
                //for (int j = 0; j < n_Electrons; j++)
                for (int j = i+1; j < n_Electrons; j++)
                {
                    // I_ab^ij
                    t2_new.at(a)(b,i,j) = integ.at(a)(b,i,j);

                    // P(ab) P(ij) [W_4] t_jk^bc
                    for (int e = n_Electrons; e < Matrix_Size; e++)
                    {
                        for (int m = 0; m < n_Electrons; m++)
                        {
                            //t2_new.at(a)(b,i,j) += t2.at(e)(b,m,j) * G3.at(a-n_Electrons)(i,e-n_Electrons,m);
                            //t2_new.at(a)(b,i,j) -= t2.at(e)(b,m,i) * G3.at(a-n_Electrons)(j,e-n_Electrons,m);
                            //t2_new.at(a)(b,i,j) -= t2.at(e)(a,m,j) * G3.at(b-n_Electrons)(i,e-n_Electrons,m);
                            //t2_new.at(a)(b,i,j) += t2.at(e)(a,m,i) * G3.at(b-n_Electrons)(j,e-n_Electrons,m);

                            t2_new.at(a)(b,i,j) += t2.at(e)(b,m,j) * W_4.at(a-n_Electrons,i).at(e-n_Electrons,m);
                            t2_new.at(a)(b,i,j) -= t2.at(e)(b,m,i) * W_4.at(a-n_Electrons,j).at(e-n_Electrons,m);
                            t2_new.at(a)(b,i,j) -= t2.at(e)(a,m,j) * W_4.at(b-n_Electrons,i).at(e-n_Electrons,m);
                            t2_new.at(a)(b,i,j) += t2.at(e)(a,m,i) * W_4.at(b-n_Electrons,j).at(e-n_Electrons,m);

                        }
                    }

                    // - P(ab) [W_2] t_k^b
                    for (int m = 0; m < n_Electrons; m++)
                    {
                        t2_new.at(a)(b,i,j) -= G1.at(a-n_Electrons)(i,j,m) * t1(b,m);
                        t2_new.at(a)(b,i,j) += G1.at(b-n_Electrons)(i,j,m) * t1(a,m);
                    }

                    // P(ij) I_ab^cj t_i^c + 1/2 I_ab^cd tau_ij^cd
                    for (int e = n_Electrons; e < Matrix_Size; e++)
                    {
                        t2_new.at(a)(b,i,j) += integ.at(a)(b,i,e) * t1(e,j);
                        t2_new.at(a)(b,i,j) -= integ.at(a)(b,j,e) * t1(e,i);

                        temp = 0;
                        for (int f = n_Electrons; f < Matrix_Size; f++)
                        {
                            temp += integ.at(a)(b,e,f) * tau.at(e)(f,i,j);
                        }
                        t2_new.at(a)(b,i,j) += 0.5*temp;
                    }

                    // 1/2 [W_1] tau_kl^ab
                    temp = 0;
                    for (int m = 0; m < n_Electrons; m++)
                    {
                        for (int n = 0; n < n_Electrons; n++)
                        {
                            temp += G2.at(i)(j,m,n) * tau.at(a)(b,m,n);
                        }
                    }
                    t2_new.at(a)(b,i,j) += 0.5*temp;

                    // P(ij) t_jk^ab [F_2]
                    for (int m = 0; m < n_Electrons; m++)
                    {
                        t2_new.at(a)(b,i,j) += t2.at(a)(b,j,m) * D2(i,m);
                        t2_new.at(a)(b,i,j) -= t2.at(a)(b,i,m) * D2(j,m);
                    }

                    // P(ab) t_ij^bc [F_3]_c^a, akta her. i følge utledninga er det motsatt fortegn...
                    for (int e = n_Electrons; e < Matrix_Size; e++)
                    {
                        t2_new.at(a)(b,i,j) -= t2.at(b)(e,i,j) * D3(a,e);
                        t2_new.at(a)(b,i,j) += t2.at(a)(e,i,j) * D3(b,e);
                    }

                    t2_new.at(a)(b,i,j) = t2_new.at(a)(b,i,j) / denom_abij.at(a)(b,i,j);

                    t2_new.at(b)(a,i,j) = -t2_new.at(a)(b,i,j);
                    t2_new.at(b)(a,j,i) = t2_new.at(a)(b,i,j);
                    t2_new.at(a)(b,j,i) = t2_new.at(b)(a,i,j);
                }
            }
        }

        // Special case where a = b, some of the terms will be 0
        for (int i = 0; i < n_Electrons; i++)
        {
            for (int j = i+1; i < n_Electrons; i++)
            {
                // I_ab^ij
                t2_new.at(a)(a,i,j) = integ.at(i)(j,a,a);

                // P(ij) I_ab^cj t_i^c + 1/2 I_ab^cd tau_ij^cd
                for (int e = n_Electrons; e < Matrix_Size; e++)
                {
                    t2_new.at(a)(a,i,j) += integ.at(a)(a,e,j) * t1(e,i);
                    t2_new.at(a)(a,i,j) -= integ.at(a)(a,e,i) * t1(e,j);

                    temp = 0;
                    for (int f = n_Electrons; f < Matrix_Size; f++)
                    {
                        temp += integ.at(e)(f,a,a) * tau.at(e)(f,i,j);
                    }
                    t2_new.at(a)(a,i,j) += 0.5*temp;
                }

                // 1/2 [W_1] tau_kl^ab
                temp = 0;
                for (int m = 0; m < n_Electrons; m++)
                {
                    for (int n = 0; n < n_Electrons; n++)
                    {
                        temp += G2.at(i)(j,m,n) * tau.at(a)(a,m,n);
                    }
                }
                t2_new.at(a)(a,i,j) += 0.5*temp;

                // P(ij) t_jk^ab [F_2]
                for (int m = 0; m < n_Electrons; m++)
                {
                    t2_new.at(a)(a,i,j) += t2.at(a)(a,j,m) * D2(i,m);
                    t2_new.at(a)(a,i,j) -= t2.at(a)(a,i,m) * D2(j,m);
                }

                t2_new.at(a)(a,i,j) = t2_new.at(a)(a,i,j) / denom_abij.at(a)(a,i,j);

                t2_new.at(a)(a,j,i) = -t2_new.at(a)(a,i,j);
            }
        }

        // Special case where i = j, some of the terms will be 0
        for (int b = a+1; b < Matrix_Size; b++)
        {
            for (int i = 0; i < n_Electrons; i++)
            {
                // I_ab^ij
                t2_new.at(a)(b,i,i) = integ.at(i)(i,a,b);

                // - P(ab) [W_2] t_k^b
                for (int m = 0; m < n_Electrons; m++)
                {
                    t2_new.at(a)(b,i,i) -= G1.at(a-n_Electrons)(i,i,m) * t1(b,m);
                    t2_new.at(a)(b,i,i) += G1.at(b-n_Electrons)(i,i,m) * t1(a,m);
                }

                // 1/2 I_ab^cd tau_ij^cd
                for (int e = n_Electrons; e < Matrix_Size; e++)
                {
                    temp = 0;
                    for (int f = n_Electrons; f < Matrix_Size; f++)
                    {
                        temp += integ.at(e)(f,a,b) * tau.at(e)(f,i,i);
                    }
                    t2_new.at(a)(b,i,i) += 0.5*temp;
                }

                // 1/2 [W_1] tau_kl^ab
                temp = 0;
                for (int m = 0; m < n_Electrons; m++)
                {
                    for (int n = 0; n < n_Electrons; n++)
                    {
                        temp += G2.at(i)(i,m,n) * tau.at(a)(b,m,n);
                    }
                }
                t2_new.at(a)(b,i,i) += 0.5*temp;

                // P(ab) t_ij^bc [F_3]_c^a, akta her. i følge utledninga er det motsatt fortegn...
                for (int e = n_Electrons; e < Matrix_Size; e++)
                {
                    t2_new.at(a)(b,i,i) -= t2.at(b)(e,i,i) * D3(a,e);
                    t2_new.at(a)(b,i,i) += t2.at(a)(e,i,i) * D3(b,e);
                }

                t2_new.at(a)(b,i,i) = t2_new.at(a)(b,i,i) / denom_abij.at(a)(b,i,i);

                t2_new.at(b)(a,i,i) = -t2_new.at(a)(b,i,i);
            }
        }

        // Special case where a = b AND i = j, meny of the terms will be 0
        for (int i = 0; i < n_Electrons; i++)
        {
            // I_ab^ij
            t2_new.at(a)(a,i,i) = integ.at(i)(i,a,a);

            // 1/2 I_ab^cd tau_ij^cd
            for (int e = n_Electrons; e < Matrix_Size; e++)
            {
                temp = 0;
                for (int f = n_Electrons; f < Matrix_Size; f++)
                {
                    temp += integ.at(e)(f,a,a) * tau.at(e)(f,i,i);
                }
                t2_new.at(a)(a,i,i) += 0.5*temp;
            }

            // 1/2 [W_1] tau_kl^ab
            temp = 0;
            for (int m = 0; m < n_Electrons; m++)
            {
                for (int n = 0; n < n_Electrons; n++)
                {
                    temp += G2.at(i)(i,m,n) * tau.at(a)(a,m,n);
                }
            }
            t2_new.at(a)(a,i,i) += 0.5*temp;

            t2_new.at(a)(a,i,i) = t2_new.at(a)(a,i,i) / denom_abij.at(a)(a,i,i);
        }
    }
}

double coupled_cluster::Calc_Energy()
{
    double E1=0;
    for (int a = n_Electrons; a < Matrix_Size; a++)
    {
        for (int i = 0; i < n_Electrons; i++)
        {
            E1 += fs(i,a) * t1(a,i);

            for (int b = n_Electrons; b < Matrix_Size; b++)
            {
                for (int j = 0; j < n_Electrons; j++)
                {
                    E1 += 0.25*integ.at(i)(j,a,b) * T_2_new.at(a,i)(b,j);
                    E1 += 0.5*integ.at(i)(j,a,b) * t1(a,i) * t1(b,j);
                }
            }
        }
    }
    return E1;
}

void coupled_cluster::Fill_tau()
{
    for (int a = n_Electrons; a < Matrix_Size; a++)
    {
        for (int b = n_Electrons; b < Matrix_Size; b++)
        {
            for (int i = 0; i < n_Electrons; i++)
            {
                for (int j = 0; j < n_Electrons; j++)
                {
                    tau.at(a)(b,i,j) = T_2_new.at(a,i)(b,j) + t1(a,i)*t1(b,j) - t1(a,j)*t1(b,i);
                }
            }
        }
    }
}

void coupled_cluster::DIIS(int t_degree)
{    
    // Method to speed up CCSD, the iteration part, T1 and T2 amplitudes
    // Algorithm based on: http://vergil.chemistry.gatech.edu/notes/diis/node2.html
    // Same exact algorithm for both methods. For now they are seperated by if tests but should be possible to combine

    // DIIS algorithm can cause errors, we put some try and except up in here to fix that

    double temp_trace;

    if (t_degree == 1)
    {
            if (iter > number_elements_DIIS)
            {
                // Store new values here
                Answer_Stored1.push_back(t1_new);
                Error_Stored1.push_back(t1_new - Answer_Stored1.at(number_elements_DIIS-1));

                // Push old values out here
                for (int i = 0; i < number_elements_DIIS; i++)
                {
                    Error_Stored1.at(i) = Error_Stored1.at(i+1);
                    Answer_Stored1.at(i) = Answer_Stored1.at(i+1);
                }

                // Delete values we dont need anymore here
                Error_Stored1.resize(number_elements_DIIS);
                Answer_Stored1.resize(number_elements_DIIS);

                // Make B matrix, which is <i|j>
                mat mat1, mat2;
                for (int i = 0; i < number_elements_DIIS; i++)
                {
                    for (int j = 0; j < number_elements_DIIS; j++)
                    {
                        mat1 = Error_Stored1.at(i);
                        mat2 = Error_Stored1.at(j);
                        temp_trace = trace(mat1.t() * mat2);
                        DIIS_B1(i,j) = temp_trace;
                    }
                }

                // Solve this DIIS_B * DIIS_C = DIIS_Z matrix problem with easy armadillo solve() function
                DIIS_c1 = solve(DIIS_B1, DIIS_Z1);

                // Construct t1 based on a pre defined past number of terms, one term outside loop to make sure no += problems with un-defined sized matrix
                t1_new = DIIS_c1.at(0) * Answer_Stored1.at(0);
                for (int i = 1; i < number_elements_DIIS; i++)
                {
                    t1_new += DIIS_c1.at(i) * Answer_Stored1.at(i);
                }
            }

            else if (iter == 0)
            {
                // Store first term, first guess for t1 is 0 so error is different here
                Answer_Stored1.push_back(t1_new);
                Error_Stored1.push_back(t1_new);
            }

            else if (iter == number_elements_DIIS)
            {
                // Since first error term is a little weird, we do this one extra term for good measure
                if (iter > number_elements_DIIS)
                {
                    // Store new values here
                    Answer_Stored1.push_back(t1_new);
                    Error_Stored1.push_back(t1_new - Answer_Stored1.at(number_elements_DIIS-1));

                    // Push old values out here
                    for (int i = 0; i < number_elements_DIIS; i++)
                    {
                        Error_Stored1.at(i) = Error_Stored1.at(i+1);
                        Answer_Stored1.at(i) = Answer_Stored1.at(i+1);
                    }

                    // Delete values we dont need anymore here
                    Error_Stored1.resize(number_elements_DIIS);
                    Answer_Stored1.resize(number_elements_DIIS);
                }
            }

            else
            {
                // Store new values here
                Answer_Stored1.push_back(t1_new);
                Error_Stored1.push_back(t1_new - Answer_Stored1.at(iter-1));
            }


    }

    else if (t_degree == 2)
    {
        // Amplitudes for T2, first rescale 4D array to 2D so we can use the exact same algorithm
        for (int i = 0; i < Matrix_Size; i++)
        {
            for (int j = 0; j < Matrix_Size; j++)
            {
                for (int k = 0; k < Matrix_Size; k++)
                {
                    for (int l = 0; l < Matrix_Size; l++)
                    {
                        DIIS_t2(i*Matrix_Size + j, k*Matrix_Size + l) = t2_new.at(i)(j,k,l);
                    }
                }
            }
        }

        // Only start after the first iteration, since P is needed and not calculated before our secound round in this function
        if (iter > number_elements_DIIS)
        {
            // Calculate error, this term is 0 when we have converged
            Error_Stored2.push_back(DIIS_t2 - Answer_Stored2.at(number_elements_DIIS-1));

            // Add the new matrixes
            Answer_Stored2.push_back(DIIS_t2);

            // This for loop ensures we only store the pre defined number of elements (this easy algo is actually predefined +1)
            for (int i = 0; i < number_elements_DIIS; i++)
            {
                Error_Stored2.at(i) = Error_Stored2.at(i+1);
                Answer_Stored2.at(i) = Answer_Stored2.at(i+1);
            }

            // Delete terms we dont need
            Error_Stored2.resize(number_elements_DIIS);
            Answer_Stored2.resize(number_elements_DIIS);

            // Make B matrix, which is <i|j>
            mat mat1, mat2;
            for (int i = 0; i < number_elements_DIIS; i++)
            {
                for (int j = 0; j < number_elements_DIIS; j++)
                {
                    mat1 = Error_Stored2.at(i);
                    mat2 = Error_Stored2.at(j);
                    temp_trace = trace(mat1.t() * mat2);
                    DIIS_B1(i,j) = temp_trace;
                }
            }

            // Solve this DIIS_B * DIIS_C = DIIS_Z matrix problem with easy armadillo solve() function
            DIIS_c1 = solve(DIIS_B1, DIIS_Z1);

            // Construct t1 based on a pre defined past number of terms, one term outside loop to make sure no += problems with un-defined sized matrix
            DIIS_t2 = DIIS_c1.at(0) * Answer_Stored2.at(0);
            for (int i = 1; i < number_elements_DIIS; i++)
            {
                DIIS_t2 += DIIS_c1.at(i) * Answer_Stored2.at(i);
            }

            // Rescale back to 4D array instead of 2D which was used for calculations
            for (int i = 0; i < Matrix_Size; i++)
            {
                for (int j = 0; j < Matrix_Size; j++)
                {
                    for (int k = 0; k < Matrix_Size; k++)
                    {
                        for (int l = 0; l < Matrix_Size; l++)
                        {
                            t2_new.at(i)(j,k,l) = DIIS_t2(i*Matrix_Size + j, k*Matrix_Size + l);
                        }
                    }
                }
            }
        }

        else if (iter == 0)
        {
            // Store first term, first guess for t1 is 0 so error is different here
            Answer_Stored2.push_back(DIIS_t2);
            Error_Stored2.push_back(DIIS_t2);
        }

        else if (iter == number_elements_DIIS)
        {
            // Calculate error, this term is 0 when we have converged
            Error_Stored2.push_back(DIIS_t2 - Answer_Stored2.at(number_elements_DIIS-1));

            // Add the new matrixes
            Answer_Stored2.push_back(DIIS_t2);

            // This for loop ensures we only store the pre defined number of elements (this easy algo is actually predefined +1)
            for (int i = 0; i < number_elements_DIIS; i++)
            {
                Error_Stored2.at(i) = Error_Stored2.at(i+1);
                Answer_Stored2.at(i) = Answer_Stored2.at(i+1);
            }

            // Delete terms we dont need
            Error_Stored2.resize(number_elements_DIIS);
            Answer_Stored2.resize(number_elements_DIIS);
        }

        else
        {
            // Store new values here
            Answer_Stored2.push_back(DIIS_t2);
            Error_Stored2.push_back(DIIS_t2 - Answer_Stored2.at(iter-1));
        }
    }
}

void coupled_cluster::Initialize_DIIS()
{
    // Initialize DIIS method, define some stuff here for speedup
    // Please see function DIIS() for info on the method
    number_elements_DIIS = 3; // This term is pre defined, but i choose 3 and not to take it as input (<-- !)
    DIIS_B1 = zeros(number_elements_DIIS+1, number_elements_DIIS+1);
    DIIS_Z1 = zeros(number_elements_DIIS+1);
    DIIS_Z1(number_elements_DIIS) = -1.0;
    DIIS_does = true;
    DIIS_t2 = zeros(Matrix_Size*Matrix_Size, Matrix_Size*Matrix_Size);
    for (int i = 0; i < number_elements_DIIS; i++)
    {
        DIIS_B1(i,number_elements_DIIS) = -1.0;
        DIIS_B1(number_elements_DIIS,i) = -1.0;
    }
}
*/
