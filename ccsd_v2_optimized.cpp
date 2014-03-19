#include "ccsd_v2_optimized.h"

ccsd_v2_optimized::ccsd_v2_optimized(int n_N, vec zz, mat rr, string B_s, int n_Elec)
{
    n_Nuclei = n_N;
    Z = zz;
    R = rr;
    Basis_Set = B_s;
    n_Electrons = n_Elec;
}

double ccsd_v2_optimized::CCSD(double toler, bool print_stuff)
{
    // General starting stuff
    bool continue_ccsd = true;
    double convergance_check;
    int Itermax = 1000;
    iter = 1;
    E_old = 0;
    E_new = 0;
    int B=0,J; // dummie variables

    // Initializion
    double E_HF;
    Hartree_Fock_Solver HartFock(n_Nuclei, Z, R, Basis_Set, n_Electrons, true);
    E_HF = HartFock.get_Energy(toler); // Calc hartree fock energy
    Matrix_Size = 2*HartFock.ReturnMatrixSize(); // This will be twice as big as in Hatree Fock, because of spin up and down
    unocc_orb = Matrix_Size - n_Electrons; // Number of unocupied orbitals
    mat temp_matrix;
    Speed_Elec = (int) n_Electrons/2;
    Speed_Occ = (int) unocc_orb/2;

    // Transform to MO basis inside the object ccint
    coupled_cluster_integrals ccint(Matrix_Size/2, HartFock.ReturnC(), n_Electrons);
    ccint.Mount_Indexed_Integrals(HartFock.Return_Indexed_Q());
    ccint.Prepear_Integralsv2(); // Most of the 2 elekctron integrals contained here
    fs = ccint.Return_FS(HartFock.return_eigval_F()); // eigenvalues of fock matrix
    HartFock.Delete_Everything();

    // We now pull out the parts of the MO basis we need, we split up the this 4D array which is a
    // n^4 array, n^4 = (n_u + n_o)^4, where n_u = unoccupied orbitals in HF basis and n_o is
    // occupied orbitals in HF basis. To the right of each term is the nr of elements stored inside it
    integ2 = ccint.Return_Rearranged_Integrals2();    // n_o^2 n_u^2   |
    integ3 = ccint.Return_Rearranged_Integrals3();    // 1/2 n_u^4     |
    integ4 = ccint.Return_Rearranged_Integrals4();    // n_o^2 n_u^2   |
    integ5 = ccint.Return_Rearranged_Integrals5();    // n_u^3 n_o     |
    integ6 = ccint.Return_Rearranged_Integrlas6();    // n_o^3 n_u     |
    integ7 = ccint.Return_Rearranged_Integrals7();    // n_u^2 n_o^2   |
    integ8 = ccint.Return_Rearranged_Integrals8();    // 1/2 n_o^4     |
    integ9 = ccint.Return_Rearranged_Integrals9();    // 1/2 n_u^3 n_o |
    integ10 = ccint.Return_Rearranged_Integrals10();  // n_u^2 n_o^2   |

    // Delete all variables inside that ccint object,
    // hopefully this is helpful before we create our variables here
    ccint.Delete_Everything();

    // Define variables used in the CCSD method:

    // Optimized 4D arrays
    tau4.set_size(unocc_orb, unocc_orb); // n_o^2 n_u^2
    t2.set_size(unocc_orb, n_Electrons); // n_o^2 n_u^2
    t2_new.set_size(unocc_orb, n_Electrons); // n_o^2 n_u^2
    W_3.set_size(n_Electrons, n_Electrons); // n_o^3 n_u
    W_4.set_size(unocc_orb, n_Electrons); // n_o^2 n_u^2
    DEN_ABIJ.set_size(unocc_orb, n_Electrons); // n_o^2 n_u^2

    // Optimized 2D arrays
    W_1 = zeros(n_Electrons, n_Electrons); // This is 4D array but can be stored as 2D if updated inside T2 amplitude loops
    W_2 = zeros(unocc_orb, n_Electrons); // This is 4D array but can be stored as 2D if updated inside T2 amplitude loops
    D3 = zeros(unocc_orb, unocc_orb);
    D2 = zeros(n_Electrons, n_Electrons);
    T_1 = zeros(unocc_orb, n_Electrons);
    T_1_new = zeros(unocc_orb, n_Electrons);
    D1 = zeros(unocc_orb, n_Electrons);
    FS_AI = zeros(unocc_orb, n_Electrons);
    FS_AB = zeros(unocc_orb, unocc_orb);
    FS_IJ = zeros(n_Electrons, n_Electrons);
    DEN_AI = zeros(unocc_orb, n_Electrons);
    tau1 = zeros(unocc_orb, unocc_orb); // This is a 2D mapping of tau

    // Initialize our 4D arrays
    for (int i = 0; i < n_Electrons; i++)
    {
        for (int j = 0; j < n_Electrons; j++)
        {
            temp_matrix = zeros(unocc_orb, n_Electrons);
            W_3(i,j) = temp_matrix;
        }
    }

    for (int i = 0; i < unocc_orb; i++)
    {
        for (int j = 0; j < n_Electrons; j++)
        {
            temp_matrix = zeros(unocc_orb, n_Electrons);
            t2(i,j) = temp_matrix;
            t2_new(i,j) = temp_matrix;
            W_4(i,j) = temp_matrix;
            DEN_ABIJ(i,j) = temp_matrix;
        }

        for (int j = i+1; j < unocc_orb; j++)
        {
            temp_matrix = zeros(n_Electrons, n_Electrons);
            tau4(i,j) = temp_matrix;
        }
    }

    for (int a = n_Electrons; a < Matrix_Size; a++)
    {
        for (int i = 0; i < n_Electrons; i++)
        {
            for (int b = n_Electrons; b < Matrix_Size; b++)
            {
                for (int j = 0; j < n_Electrons; j++)
                {
                    if ((a+b+i+j)%2 == 0)
                    {
                        DEN_ABIJ(a-n_Electrons, i)((b-n_Electrons), j) = 1/(fs(i,i) + fs(j,j) - fs(a,a) - fs(b,b));
                    }
                }
            }
        }
    }

    for (int a = 0; a < unocc_orb; a++)
    {
        for (int i = 0; i < n_Electrons; i++)
        {
            DEN_AI(a,i) = 1/(fs(i,i) - fs(a+n_Electrons,a+n_Electrons));
            i++;
        }
        a++;
    }

    for (int a = 0+1; a < unocc_orb; a++)
    {
        for (int i = 1; i < n_Electrons; i++)
        {
            DEN_AI(a,i) = 1/(fs(i,i) - fs(a+n_Electrons,a+n_Electrons));
            i++;
        }
        a++;
    }

    for (int i = 0; i < n_Electrons; i++)
    {
        for (int j = 0; j < n_Electrons; j++)
        {
            FS_IJ(i, j) = (1 - EqualFunc(i,j)) * fs(i,j);
            j++;
        }
        i++;
    }

    for (int i = 1; i < n_Electrons; i++)
    {
        for (int j = 1; j < n_Electrons; j++)
        {
            FS_IJ(i, j) = (1 - EqualFunc(i,j)) * fs(i,j);
            j++;
        }
        i++;
    }

    for (int a = 0; a < unocc_orb; a++)
    {
        for (int i = 0; i < n_Electrons; i++)
        {
            FS_AI(a,i) = fs(a+n_Electrons,i);
            i++;
        }

        for (int b = 0; b < unocc_orb; b++)
        {
            FS_AB(a,b) = (1 - EqualFunc(a,b)) * fs(a+n_Electrons,b+n_Electrons);
            b++;
        }
        a++;
    }

    for (int a = 1; a < unocc_orb; a++)
    {
        for (int i = 1; i < n_Electrons; i++)
        {
            FS_AI(a, i) = fs(a+n_Electrons,i);
            i++;
        }

        for (int b = 1; b < unocc_orb; b++)
        {
            FS_AB(a, b) = (1 - EqualFunc(a,b)) * fs(a+n_Electrons,b+n_Electrons);
            b++;
        }
        a++;
    }

    fs.clear(); // Dont need fs anymore

    for (int a = 0; a < unocc_orb; a++)
    {
        for (int b = 0; b < unocc_orb; b++)
        {
            for (int i = 0; i < n_Electrons; i++)
            {
                if ((a-b+i)%2 == 0)
                {
                    for (int j = 0; j < n_Electrons; j++)
                    {
                        t2(a,i)(b,j) = integ4.at(i,j)(a,b) * DEN_ABIJ(a,i)(b,j);
                        j++;
                    }
                }

                else
                {
                    for (int j = 1; j < n_Electrons; j++)
                    {
                        t2(a,i)(b,j) = integ4.at(i,j)(a,b) * DEN_ABIJ(a,i)(b,j);
                        j++;
                    }
                }
            }
            b++;
        }

        for (int b = 1; b < unocc_orb; b++)
        {
            for (int i = 0; i < n_Electrons; i++)
            {
                if ((a-b+i)%2 == 0)
                {
                    for (int j = 0; j < n_Electrons; j++)
                    {
                        t2(a,i)(b,j) = integ4.at(i,j)(a,b) * DEN_ABIJ(a,i)(b,j);
                        j++;
                    }
                }

                else
                {
                    J = 0;
                    for (int j = 1; j < n_Electrons; j++)
                    {
                        t2(a,i)(b,j) = integ4.at(i,j)(a,b) * DEN_ABIJ(a,i)(b,j);
                        j++;
                    }
                }
            }
            b++;
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
        Fill_W3();
        Fill_W4();
        Fill_F1();
        Fill_F2();
        Fill_F3();

        // Find new amplitudes
        Fill_t1_new();
        Fill_t2_new();

        // Update amplitudes
        T_1 = T_1_new;
        t2 = t2_new;

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

    cout << "CCSD correction [au] = " << E_new << endl;

    //cout << t2(1,1) << endl;
    return E_new+E_HF;
}

int ccsd_v2_optimized::EqualFunc(int a, int b)
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

double ccsd_v2_optimized::Calc_Energy()
{
    double E1 = 0;
    for (int i = 0; i < n_Electrons; i++)
    {
        for (int j = i+1; j < n_Electrons; j++)
        {
            Fill_2D_tau(i,j);
            E1 += accu(tau1 % integ4(i,j));
            E1 -= accu(tau1 % integ4(j,i));
        }
    }
    E1 *= 0.25;
    E1 += accu(FS_AI % T_1);
    return E1;
}

void ccsd_v2_optimized::Fill_W1(int i, int j)
{
    // Matrix symmetric with W_1(i,j) = W_1(j,i), the terms on the "diagonal" of i==j will be equal to 0 (<-- !)
    // We actually dont even need to store anything in W_1(j,i) since it is accessed symmetricly later on also
    // This halves our storage requirements. Also compress matrix. This reduce storage by 75% + not stored diagonal.

            if ((i+j)%2 == 1)
            {
                for (int k = 0; k < n_Electrons; k++)
                {
                    for (int l = 1; l < n_Electrons; l++)
                    {
                        W_1(k,l) = integ8.at(i,j)(k,l);
                        W_1(k,l) += accu(integ6.at(k,l)(span(),j) % T_1(span(), i));
                        W_1(k,l) -= accu(integ6.at(k,l)(span(),i) % T_1(span(), j));
                        W_1(k,l) += 0.5*accu(integ4.at(k,l) % tau1);
                        l++;
                    }
                    k++;
                }

                for (int k = 1; k < n_Electrons; k++)
                {
                    for (int l = 0; l < n_Electrons; l++)
                    {
                        W_1(k,l) = integ8.at(i,j)(k,l);
                        W_1(k,l) += accu(integ6.at(k,l)(span(),j) % T_1(span(), i));
                        W_1(k,l) -= accu(integ6.at(k,l)(span(),i) % T_1(span(), j));
                        W_1(k,l) += 0.5*accu(integ4.at(k,l) % tau1);
                        l++;
                    }
                    k++;
                }
            }

            else
            {
                for (int k = 0; k < n_Electrons; k++)
                {
                    for (int l = 0; l < n_Electrons; l++)
                    {
                        W_1(k,l) = integ8.at(i,j)(k, l);
                        W_1(k,l) += accu(integ6.at(k,l)(span(),j) % T_1(span(), i));
                        W_1(k,l) -= accu(integ6.at(k,l)(span(),i) % T_1(span(), j));
                        W_1(k,l) += 0.5*accu(integ4.at(k,l) % tau1);
                        l++;
                    }
                    k++;
                }

                for (int k = 1; k < n_Electrons; k++)
                {
                    for (int l = 1; l < n_Electrons; l++)
                    {
                        W_1(k,l) = integ8.at(i,j)(k, l);
                        W_1(k,l) += accu(integ6.at(k,l)(span(),j) % T_1(span(), i));
                        W_1(k,l) -= accu(integ6.at(k,l)(span(),i) % T_1(span(), j));
                        W_1(k,l) += 0.5*accu(integ4.at(k,l) % tau1);
                        l++;
                    }
                    k++;
                }
            }
}

void ccsd_v2_optimized::Fill_W2(int i, int j)
{
    W_2 = integ6.at(i,j);

    if ((i+j)%2==0)
    {
        for (int a = 0; a < unocc_orb; a++)
        {
            for (int m = 0; m < n_Electrons; m++)
            {
                if (i%2 != 1 && j%2 != 1)
                {
                    W_2(a,m) += accu(integ7.at(a,m)(span(), j) % T_1(span(),i));
                    W_2(a,m) -= accu(integ7.at(a,m)(span(), i) % T_1(span(),j));
                    W_2(a,m) += 0.5*accu(integ5.at(a,m) % tau1);
                }
                m++;
            }
            a++;
        }

        for (int a = 1; a < unocc_orb; a++)
        {
            for (int m = 1; m < n_Electrons; m++)
            {
                //  disse kommer til å ha en slags symmetri
                W_2(a,m) += accu(integ7.at(a,m)(span(), j) % T_1(span(),i));
                W_2(a,m) -= accu(integ7.at(a,m)(span(), i) % T_1(span(),j));
                W_2(a,m) += 0.5*accu(integ5.at(a,m) % tau1);
                m++;
            }
            a++;
        }
    }

    else
    {
        for (int a = 1; a < unocc_orb; a++)
        {
            for (int m = 0; m < n_Electrons; m++)
            {
                W_2(a,m) += accu(integ7.at(a,m)(span(), j) % T_1(span(),i));
                W_2(a,m) -= accu(integ7.at(a,m)(span(), i) % T_1(span(),j));
                W_2(a,m) += 0.5*accu(integ5.at(a,m) % tau1);
                m++;
            }
            a++;
        }

        for (int a = 0; a < unocc_orb; a++)
        {
            for (int m = 1; m < n_Electrons; m++)
            {
                W_2(a,m) += accu(integ7.at(a,m)(span(), j) % T_1(span(),i));
                W_2(a,m) -= accu(integ7.at(a,m)(span(), i) % T_1(span(),j));
                W_2(a,m) += 0.5*accu(integ5.at(a,m) % tau1);
                m++;
            }
            a++;
        }
    }
}

void ccsd_v2_optimized::Fill_W3()
{
    for (int m = 0; m < n_Electrons; m++)
    {
        for (int i = 0; i < n_Electrons; i++)
        {
            if ((m+i)%2 == 0)
            {
                for (int e = 0; e < unocc_orb; e++)
                {
                    for (int n = 0; n < n_Electrons; n++)
                    {
                        W_3(i,m)(e,n) = integ6.at(m,n)(e, i);
                        W_3(i,m)(e,n) -= accu(integ4.at(m,n)(span(), e) % T_1(span(), i));
                        n++;
                    }
                    e++;
                }

                for (int e = 1; e < unocc_orb; e++)
                {
                    for (int n = 1; n < n_Electrons; n++)
                    {
                        W_3(i,m)(e,n) = integ6.at(m,n)(e, i);
                        W_3(i,m)(e,n) -= accu(integ4.at(m,n)(span(), e) % T_1(span(), i));
                        n++;
                    }
                    e++;
                }
            }

            else
            {
                if (i%2 != 0 && m%2 != 1)
                {
                    for (int e = 0; e < unocc_orb; e++)
                    {
                        for (int n = 1; n < n_Electrons; n++)
                        {
                            W_3(i,m)(e,n) = integ6.at(m,n)(e, i);
                            W_3(i,m)(e,n) -= accu(integ4.at(m,n)(span(), e) % T_1(span(), i));
                            n++;
                        }
                        e++;
                    }
                }

                if (i%2 != 1 && m%2 != 0)
                {
                    for (int e = 1; e < unocc_orb; e++)
                    {
                        for (int n = 0; n < n_Electrons; n++)
                        {
                            W_3(i,m)(e,n) = integ6.at(m,n)(e, i);
                            W_3(i,m)(e,n) -= accu(integ4.at(m,n)(span(), e) % T_1(span(), i));
                            n++;
                        }
                        e++;
                    }
                }
            }
        }
    }
}

void ccsd_v2_optimized::Fill_W4()
{
    for (int a = 0; a < unocc_orb; a++)
    {
        for (int i = 0; i < n_Electrons; i++)
        {
            for (int e = 0; e < unocc_orb; e++)
            {
                if ((a-e+i)%2 == 0)
                {
                    for (int m = 0; m < n_Electrons; m++)
                    {
                        W_4.at(a,i).at(e,m) = -integ7.at(a,m)(e,i);
                        W_4.at(a,i).at(e,m) -= accu(W_3.at(i,m)(e,span()) % T_1.row(a));
                        W_4.at(a,i).at(e,m) += accu(integ5.at(a,m)(span(),e) % T_1(span(),i));
                        W_4.at(a,i).at(e,m) += 0.5*accu(integ2.at(e,m) % t2.at(a,i));
                        m++;
                    }
                }

                else
                {
                    for (int m = 1; m < n_Electrons; m++)
                    {
                        W_4.at(a,i).at(e,m) = -integ7.at(a,m)(e,i);

                        if (i%2 != 0)
                        {
                            W_4.at(a,i).at(e,m) -= accu(W_3.at(i,m)(e,span()) % T_1.row(a));
                        }

                        if (a%2 != 1)
                        {
                            W_4.at(a,i).at(e,m) += accu(integ5.at(a,m)(span(),e) % T_1(span(),i));
                        }
                        W_4.at(a,i).at(e,m) += 0.5*accu(integ2.at(e,m) % t2.at(a,i));
                        m++;
                    }
                }
                e++;
            }
        }

        for (int e = 1; e < unocc_orb; e++)
        {
            for (int i = 0; i < n_Electrons; i++)
            {
                if ((a-e+i)%2 == 0)
                {
                    for (int m = 0; m < n_Electrons; m++)
                    {
                        if (a%2 != 0 && i%2 != 1)
                        {
                            W_4.at(a,i).at(e,m) = -integ7.at(a,m)(e,i);
                            W_4.at(a,i).at(e,m) -= accu(W_3.at(i,m)(e,span()) % T_1.row(a));
                            W_4.at(a,i).at(e,m) += accu(integ5.at(a,m)(span(),e) % T_1(span(),i));
                            W_4.at(a,i).at(e,m) += 0.5*accu(integ2.at(e,m) % t2.at(a,i));
                        }
                        m++;
                    }
                }

                else
                {
                    for (int m = 1; m < n_Electrons; m++)
                    {
                        W_4.at(a,i).at(e,m) = -integ7.at(a,m)(e,i);
                        W_4.at(a,i).at(e,m) -= accu(W_3.at(i,m)(e,span()) % T_1.row(a));
                        W_4.at(a,i).at(e,m) += accu(integ5.at(a,m)(span(),e) % T_1(span(),i));
                        W_4.at(a,i).at(e,m) += 0.5*accu(integ2.at(e,m) % t2.at(a,i));
                        m++;
                    }
                }
            }
            e++;
        }
    }
}

void ccsd_v2_optimized::Fill_F1()
{
    D1 = FS_AI;

    for (int m = 0; m < n_Electrons; m++)
    {
        for (int a = 0; a < unocc_orb; a++)
        {
            D1.at(a,m) += accu(integ2.at(a,m) % T_1);
            a++;
        }
        m++;
    }

    for (int m = 1; m < n_Electrons; m++)
    {
        for (int a = 1; a < unocc_orb; a++)
        {
            D1.at(a,m) += accu(integ2.at(a,m) % T_1);
            a++;
        }
        m++;
    }
}

void ccsd_v2_optimized::Fill_F2()
{
    D2 = FS_IJ;

    for (int m = 0; m < n_Electrons; m++)
    {
        for (int i = 0; i < n_Electrons; i++)
        {
            D2.at(m,i) += accu(D1(span(), m) % T_1(span(), i));

            for (int e = 0; e < unocc_orb; e++)
            {
                D2.at(m,i) += 0.5*accu(integ2.at(e,m) % t2.at(e,i));
            }

            for (int k = 0; k < n_Electrons; k++)
            {
                D2.at(m,i) -= accu(integ6.at(m,k)(span(),i) % T_1(span(),k));
            }
            i++;
        }
        m++;
    }

    for (int m = 1; m < n_Electrons; m++)
    {
        for (int i = 1; i < n_Electrons; i++)
        {
            D2.at(m,i) += accu(D1(span(), m) % T_1(span(), i));
            for (int e = 0; e < unocc_orb; e++)
            {
                D2.at(m,i) += 0.5*accu(integ2.at(e,m) % t2.at(e,i));
            }
            for (int k = 0; k < n_Electrons; k++)
            {
                D2.at(m,i) -= accu(integ6.at(m,k)(span(), i) % T_1(span(), k));
            }
            i++;
        }
        m++;
    }
}

void ccsd_v2_optimized::Fill_F3()
{
    D3 = FS_AB;
    for (int a = 0; a < unocc_orb; a++)
    {
        for (int e = 0; e < unocc_orb; e++)
        {
            D3.at(a,e) -= accu(D1.row(e) % T_1.row(a));

            for (int m = 1; m < n_Electrons; m++)
            {
                D3.at(a,e) += accu(integ5.at(a,m)(e, span()) % T_1(span(), m).t());
                D3.at(a,e) -= 0.5*accu(integ2.at(e,m) % t2.at(a,m));
                m++;
            }

            for (int m = 0; m < n_Electrons; m++)
            {
                D3.at(a,e) += accu(integ5.at(a,m)(e, span()) % T_1(span(), m).t());
                D3.at(a,e) -= 0.5*accu(integ2.at(e,m) % t2.at(a,m));
                m++;
            }
            e++;
        }
        a++;
    }

    for (int a = 1; a < unocc_orb; a++)
    {
        for (int e = 1; e < unocc_orb; e++)
        {
            D3.at(a,e) -= accu(D1.row(e) % T_1.row(a));

            for (int m = 1; m < n_Electrons; m++)
            {
                D3.at(a,e) += accu(integ5.at(a,m)(e, span()) % T_1(span(), m).t());
                D3.at(a,e) -= 0.5*accu(integ2.at(e,m) % t2.at(a,m));
                m++;
            }

            for (int m = 0; m < n_Electrons; m++)
            {
                D3.at(a,e) += accu(integ5.at(a,m)(e, span()) % T_1(span(), m).t());
                D3.at(a,e) -= 0.5*accu(integ2.at(e,m) % t2.at(a,m));
                m++;
            }
            e++;
        }
        a++;
    }
}

void ccsd_v2_optimized::Fill_tau()
{
    for (int a = 0; a < unocc_orb; a++)
    {
        for (int b = a+1; b < unocc_orb; b++)
        {
            if ((a+b)%2 == 1)
            {
                for (int i = 0; i < n_Electrons; i++)
                {
                    for (int j = 1; j < n_Electrons; j++)
                    {
                        tau4(a,b)(i,j) = t2.at(a,i)(b,j);
                        if ((a+i)%2 == 0 && (b+j)%2 == 0)
                        {
                            tau4(a,b)(i,j) += T_1(a, i) * T_1(b, j);
                        }

                        if ((a+j)%2 ==0 && (b+i)%2 == 0)
                        {
                            tau4(a,b)(i,j) -= T_1(a, j) * T_1(b, i);
                        }
                        j++;
                    }
                    i++;
                }

                for (int i = 1; i < n_Electrons; i++)
                {
                    for (int j = 0; j < n_Electrons; j++)
                    {
                        tau4(a,b)(i,j) = t2.at(a,i)(b,j);
                        if ((a+i)%2 == 0 && (b+j)%2 == 0)
                        {
                            tau4(a,b)(i,j) += T_1(a, i) * T_1(b, j);
                        }

                        if ((a+j)%2 ==0 && (b+i)%2 == 0)
                        {
                            tau4(a,b)(i,j) -= T_1(a, j) * T_1(b, i);
                        }
                        j++;
                    }
                    i++;
                }
            }

            else
            {
                for (int i = 0; i < n_Electrons; i++)
                {
                    for (int j = 0; j < n_Electrons; j++)
                    {
                        tau4(a,b)(i, j) = t2.at(a,i)(b,j);
                        if ((a+i)%2 == 0 && (b+j)%2 == 0)
                        {
                            tau4(a,b)(i,j) += T_1(a, i) * T_1(b, j);
                        }

                        if ((a+j)%2 ==0 && (b+i)%2 == 0)
                        {
                            tau4(a,b)(i,j) -= T_1(a, j) * T_1(b, i);
                        }
                        j++;
                    }
                    i++;
                }

                for (int i = 1; i < n_Electrons; i++)
                {
                    for (int j = 1; j < n_Electrons; j++)
                    {
                        tau4(a,b)(i,j) = t2.at(a,i)(b,j);

                        if ((a+i)%2 == 0 && (b+j)%2 == 0)
                        {
                            tau4(a,b)(i,j) += T_1(a, i) * T_1(b, j);
                        }

                        if ((a+j)%2 ==0 && (b+i)%2 == 0)
                        {
                            tau4(a,b)(i,j) -= T_1(a, j) * T_1(b, i);
                        }
                        j++;
                    }
                    i++;
                }
            }
        }
    }
}

void ccsd_v2_optimized::Fill_2D_tau(int i, int j)
{
    if ((i+j)%2 == 1)
    {
        for (int a = 0; a < unocc_orb; a++)
        {
            for (int b = 1; b < unocc_orb; b++)
            {
                tau1(a,b) = t2.at(a,i)(b,j);

                if ((a+i)%2 == 0 && (b+j)%2 == 0)
                {
                    tau1(a,b) += T_1(a, i) * T_1(b, j);
                }

                if ((a+j)%2 ==0 && (b+i)%2 == 0)
                {
                    tau1(a,b) -= T_1(a, j) * T_1(b, i);
                }
                b++;
            }
            a++;
        }

        for (int a = 1; a < unocc_orb; a++)
        {
            for (int b = 0; b < unocc_orb; b++)
            {
                tau1(a,b) = t2.at(a,i)(b,j);

                if ((a+i)%2 == 0 && (b+j)%2 == 0)
                {
                    tau1(a,b) += T_1(a, i) * T_1(b, j);
                }

                if ((a+j)%2 ==0 && (b+i)%2 == 0)
                {
                    tau1(a,b) -= T_1(a, j) * T_1(b, i);
                }
                b++;
            }
            a++;
        }
    }

    else
    {
        for (int a = 0; a < unocc_orb; a++)
        {
            for (int b = 0; b < unocc_orb; b++)
            {
                tau1(a,b) = t2.at(a,i)(b,j);

                if ((a+i)%2 == 0 && (b+j)%2 == 0)
                {
                    tau1(a,b) += T_1(a, i) * T_1(b, j);
                }

                if ((a+j)%2 ==0 && (b+i)%2 == 0)
                {
                    tau1(a,b) -= T_1(a, j) * T_1(b, i);
                }
                b++;
            }
            a++;
        }

        for (int a = 1; a < unocc_orb; a++)
        {
            for (int b = 1; b < unocc_orb; b++)
            {
                tau1(a,b) = t2.at(a,i)(b,j);

                if ((a+i)%2 == 0 && (b+j)%2 == 0)
                {
                    tau1(a,b) += T_1(a, i) * T_1(b, j);
                }

                if ((a+j)%2 ==0 && (b+i)%2 == 0)
                {
                    tau1(a,b) -= T_1(a, j) * T_1(b, i);
                }

                b++;
            }
            a++;
        }
    }
}

void ccsd_v2_optimized::Fill_t1_new()
{
    double temp;

    // f_ai
    T_1_new = FS_AI;

    // Both even numbers != 0
    for (int a = 0; a < unocc_orb; a++)
    {
        for (int i = 0; i < n_Electrons; i++)
        {
            // t_ik^ac [F_1]_c^k
            T_1_new(a,i) += accu(D1 % t2(a,i));

            // - t_k^a [F_2]_i^k
            T_1_new(a,i) -= accu(D2(span(), i) % T_1.row(a).t());

            // f_ac t_i^c
            T_1_new(a,i) += accu(FS_AB(a, span()) % T_1(span(), i).t());

            // I_ka^ci t_k^c
            T_1_new(a,i) -= accu(integ10(a,i) % T_1);

            // 1/2 I_ka^cd tau_ij^cd
            temp = 0;
            for (int k = 0; k < i; k++)
            {
                Fill_2D_tau(k,i);
                temp -= accu(integ5.at(a,k) % tau1);
            }
            T_1_new(a,i) += 0.5*temp;

            // - 1/2 I_kl^ci t_kl^ca - 1/2 I_kl^cd t_kl^ca t_i^d
            temp = 0;
            for (int k = 0; k < n_Electrons; k++)
            {
                temp += accu(W_3(i,k) % t2(a,k));
            }
            T_1_new(a,i) += 0.5*temp;

            i++;
        }
        a++;
    }

    // Both odd numbers also != 0
    for (int a = 1; a < unocc_orb; a++)
    {
        for (int i = 1; i < n_Electrons; i++)
        {
            // t_ik^ac [F_1]_c^k
            T_1_new(a,i) += accu(D1 % t2(a,i));

            // - t_k^a [F_2]_i^k
            T_1_new(a,i) -= accu(D2(span(), i) % T_1.row(a).t());

            // f_ac t_i^c
            T_1_new(a,i) += accu(FS_AB(a, span()) % T_1(span(), i).t());

            // I_ka^ci t_k^c
            T_1_new(a,i) -= accu(integ10(a,i) % T_1);

            // 1/2 I_ka^cd tau_ij^cd
            temp = 0;
            for (int k = 0; k < i; k++)
            {
                Fill_2D_tau(k,i);
                temp -= accu(integ5.at(a,k) % tau1);
            }
            T_1_new(a,i) += 0.5*temp;

            // - 1/2 I_kl^ci t_kl^ca - 1/2 I_kl^cd t_kl^ca t_i^d
            temp = 0;
            for (int k = 0; k < n_Electrons; k++)
            {
                temp += accu(W_3(i,k) % t2(a,k));
            }
            T_1_new(a,i) += 0.5*temp;

            i++;
        }
        a++;
    }
}

void ccsd_v2_optimized::Fill_t2_new()
{
    Fill_tau();

    for (int i = 0; i < n_Electrons; i++)
    {
        for (int j = i+1; j < n_Electrons; j++)
        {
            Fill_2D_tau(i,j);
            Fill_W1(i,j);
            Fill_W2(i,j);
            for (int a = 0; a < unocc_orb; a++)
            {
                if ((a+i)%2 == 0)
                {
                    T_1_new(a, i) += 0.5*accu(integ5.at(a,j) % tau1); // One part of T1 calculated here because tau1 is stored 2D
                }
                for (int b = a+1; b < unocc_orb; b++)
                {
                    if ((a+i+b+j)%2 == 0)
                    {
                        // I_ab^ij
                        t2_new(a,i)(b, j) = integ4.at(i,j)(a,b);

                        // P(ab) P(ij) [W_4] t_jk^bc, ARMADILLO
                        t2_new(a,i)(b, j) += accu(t2.at(b,j) % W_4.at(a,i));
                        t2_new(a,i)(b, j) -= accu(t2.at(b,i) % W_4.at(a,j));
                        t2_new(a,i)(b, j) -= accu(t2.at(a,j) % W_4.at(b,i));
                        t2_new(a,i)(b, j) += accu(t2.at(a,i) % W_4.at(b,j));

                        // - P(ab) [W_2] t_k^b, ARMADILLO
                        t2_new(a,i)(b, j) -= accu(W_2(a, span()) % T_1.row(b));
                        t2_new(a,i)(b, j) += accu(W_2(b, span()) % T_1.row(a));

                        // 1/2 I_ab^cd tau_ij^cd, ARMADILLO
                        t2_new(a,i)(b, j) += 0.5*accu(integ3.at(a,b) % tau1);

                        // P(ij) I_ab^cj t_i^c, ARMADILLO
                        t2_new(a,i)(b, j) -= accu(integ9.at(a,b)(span(), i) % T_1(span(), j));
                        t2_new(a,i)(b, j) += accu(integ9.at(a,b)(span(), j) % T_1(span(), i));

                        // 1/2 [W_1] tau_kl^ab, ARMADILLO
                        t2_new(a,i)(b, j) += 0.5*accu(W_1 % tau4.at(a,b));

                        // P(ij) t_jk^ab [F_2] + P(ab) t_ij^bc [F_3]_c^a, ARMADILLO
                        t2_new(a,i)(b, j) += accu(t2.at(a,j)(b, span()) % D2.row(i));
                        t2_new(a,i)(b, j) -= accu(t2.at(a,i)(b, span()) % D2.row(j));
                        t2_new(a,i)(b, j) -= accu(t2.at(b,i)(span(), j) % D3.row(a).t());
                        t2_new(a,i)(b, j) += accu(t2.at(a,i)(span(), j) % D3.row(b).t());

                        // Symmetries
                        t2_new(b,i)(a, j) = -t2_new(a,i)(b, j);
                        t2_new(a,j)(b, i) = t2_new(b,i)(a, j);
                        t2_new(b,j)(a, i) = t2_new(a,i)(b, j);
                    }
                    b++;
                }
            }
            j++;
        }
    }

    for (int i = 0; i < n_Electrons; i++)
    {
        for (int j = i+2; j < n_Electrons; j++)
        {
            Fill_2D_tau(i,j);
            Fill_W1(i,j);
            Fill_W2(i,j);
            E_new += accu(tau1 % integ4(i,j));
            E_new -= accu(tau1 % integ4(j,i));
            for (int a = 0; a < unocc_orb; a++)
            {
                if ((a+i)%2 == 0)
                {
                    T_1_new(a, i) += 0.5*accu(integ5.at(a,j) % tau1); // One part of T1 calculated here because tau1 is stored 2D
                }

                for (int b = a+2; b < unocc_orb; b++)
                {
                    if (a%2 != 0 || b%2 != 0 || i%2 != 1 || j%2 != 1)
                    {
                    if ((a+i+b+j)%2 == 0)
                    {
                        // I_ab^ij
                        t2_new(a,i)(b, j) = integ4.at(i,j)(a,b);

                        // P(ab) P(ij) [W_4] t_jk^bc, ARMADILLO
                        t2_new(a,i)(b, j) += accu(t2.at(b,j) % W_4.at(a,i));
                        t2_new(a,i)(b, j) -= accu(t2.at(b,i) % W_4.at(a,j));
                        t2_new(a,i)(b, j) -= accu(t2.at(a,j) % W_4.at(b,i));
                        t2_new(a,i)(b, j) += accu(t2.at(a,i) % W_4.at(b,j));

                        // - P(ab) [] t_k^b, ARMADILLO
                        t2_new(a,i)(b, j) -= accu(W_2(a, span()) % T_1.row(b));
                        t2_new(a,i)(b, j) += accu(W_2(b, span()) % T_1.row(a));

                        // 1/2 I_ab^cd tau_ij^cd, ARMADILLO
                        t2_new(a,i)(b, j) += 0.5*accu(integ3.at(a,b) % tau1);

                        // P(ij) I_ab^cj t_i^c, ARMADILLO
                        t2_new(a,i)(b, j) -= accu(integ9.at(a,b)(span(), i) % T_1(span(), j));
                        t2_new(a,i)(b, j) += accu(integ9.at(a,b)(span(), j) % T_1(span(), i));

                        // 1/2 [W_1] tau_kl^ab, ARMADILLO
                        t2_new(a,i)(b, j) += 0.5*accu(W_1 % tau4.at(a,b));

                        // P(ij) t_jk^ab [F_2] + P(ab) t_ij^bc [F_3]_c^a, ARMADILLO
                        t2_new(a,i)(b, j) -= accu(t2.at(b,i)(span(), j) % D3.row(a).t());
                        t2_new(a,i)(b, j) += accu(t2.at(a,i)(span(), j) % D3.row(b).t());
                        t2_new(a,i)(b, j) += accu(t2.at(a,j)(b, span()) % D2.row(i));
                        t2_new(a,i)(b, j) -= accu(t2.at(a,i)(b, span()) % D2.row(j));

                        // Symmetries
                        t2_new(b,i)(a, j) = -t2_new(a,i)(b, j);
                        t2_new(a,j)(b, i) = t2_new(b,i)(a, j);
                        t2_new(b,j)(a, i) = t2_new(a,i)(b, j);
                    }
                    }

                    b++;
                }
            }
            j++;
        }
    }

    // Get denominator matrix
    T_1_new = T_1_new % DEN_AI;
    for (int a = 0; a < unocc_orb; a++)
    {
        for (int i = 0; i < n_Electrons; i++)
        {
            t2_new(a,i) = t2_new(a,i) % DEN_ABIJ(a,i);
        }
    }
}








