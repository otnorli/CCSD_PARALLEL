#include "ccsd_even.h"

ccsd_even::ccsd_even(int n_N, vec zz, mat rr, string B_s, int n_Elec)
{
    n_Nuclei = n_N;
    Z = zz;
    R = rr;
    Basis_Set = B_s;
    n_Electrons = n_Elec;
}

double ccsd_even::CCSD(double toler, bool print_stuff)
{
    bool continue_ccsd = true;
    double convergance_check;
    int Itermax = 1000;
    iter = 0;
    E_old = 0;
    E_new = 0;
    int B=0,J; // dummie variables

    // Initializion
    double E_HF;
    Hartree_Fock_Solver HartFock(n_Nuclei, Z, R, Basis_Set, n_Electrons, true);
    E_HF = HartFock.get_Energy(toler); // Calc hartree fock energy
    Matrix_Size = 2*HartFock.ReturnMatrixSize(); // This will be twice as big as in Hatree Fock, because of spin up and down
    unocc_orb = Matrix_Size - n_Electrons; // Number of unocupied orbitals
    mat temp_matrix; cube temp_mat;
    Speed_Elec = (int) n_Electrons/2;
    Speed_Occ = (int) unocc_orb/2;

    // Defining arrays to store stuff, these are the ones that fill up all the memory
    denom_ai = zeros(unocc_orb, n_Electrons);

    // Transform to MO basis
    coupled_cluster_integrals ccint(Matrix_Size/2, HartFock.ReturnC(), n_Electrons);
    //ccint.Return_Integrals(HartFock.ReturnQ(), 2);    // Most of the 2 elekctron integrals contained here
    fs = ccint.Return_FS(HartFock.return_eigval_F()); // eigenvalues of fock matrix
    integ2 = ccint.Return_Rearranged_Integrals2();    // Splitting up 2 elektron integrals for easier matrix multiplication use later on
    integ3 = ccint.Return_Rearranged_Integrals3();    // Each of these contain for example
    integ4 = ccint.Return_Rearranged_Integrals4();    // <ab||ij> where a = unoccupied
    integ5 = ccint.Return_Rearranged_Integrals5();    // and i = unoccupied in HF basis
    integ6 = ccint.Return_Rearranged_Integrlas6();    // Each of these are approx something like
    integ7 = ccint.Return_Rearranged_Integrals7();    // 1/2 n_o^2 * n_u^2 big, where n_o = nr of occupied
    integ8 = ccint.Return_Rearranged_Integrals8();    // and n_u is nr of unoccupied
    integ9 = ccint.Return_Rearranged_Integrals9();    // Together all these 9 arrays are much much
    integ10 = ccint.Return_Rearranged_Integrals10();  // smaller than n^4 = (n_o + n_u)^4 which is very big
    ccint.Delete_Everything();
    // Delete all variables inside that ccint object,
    // hopefully this is helpful before we create our variables here

    // Optimized versions
    W_22 = zeros(n_Electrons, Speed_Elec);
    W_1.set_size(n_Electrons, n_Electrons);
    W_3.set_size(n_Electrons, n_Electrons);
    W_2.set_size(n_Electrons, n_Electrons);
    t2.set_size(unocc_orb, n_Electrons);
    t2_new.set_size(unocc_orb, n_Electrons);
    W_4.set_size(unocc_orb, n_Electrons);
    D3 = zeros(unocc_orb, Speed_Occ);
    D2 = zeros(n_Electrons, Speed_Elec);
    T_1 = zeros(unocc_orb, Speed_Elec);
    T_1_new = zeros(unocc_orb, Speed_Elec);
    D1 = zeros(unocc_orb, Speed_Elec);
    tau3.set_size(n_Electrons, n_Electrons);
    tau4.set_size(unocc_orb, unocc_orb);
    FS_AI = zeros(unocc_orb, Speed_Elec);
    FS_AB = zeros(unocc_orb, Speed_Occ);
    FS_IJ = zeros(n_Electrons, Speed_Elec);
    DEN_AI = zeros(unocc_orb, Speed_Elec);
    DEN_ABIJ.set_size(unocc_orb, n_Electrons);

    // Initialize our arrays
    for (int i = 0; i < n_Electrons; i++)
    {
        for (int j = 0; j < n_Electrons; j++)
        {
            temp_matrix = zeros(unocc_orb, Speed_Occ);
            tau3(i,j) = temp_matrix;

            temp_matrix = zeros(unocc_orb, Speed_Elec);
            W_3(i,j) = temp_matrix;
        }
    }

    for (int i = 0; i < n_Electrons; i++)
    {
        for (int j = i+1; j < n_Electrons; j++)
        {
            temp_matrix = zeros(unocc_orb, Speed_Elec);
            W_2(i,j) = temp_matrix;

            temp_matrix = zeros(n_Electrons, Speed_Elec);
            W_1(i,j) = temp_matrix;
        }
    }

    for (int i = 0; i < unocc_orb; i++)
    {
        temp_mat = zeros(unocc_orb, n_Electrons, n_Electrons);
        denom_abij.push_back(temp_mat);

        for (int j = 0; j < n_Electrons; j++)
        {
            temp_matrix = zeros(unocc_orb, Speed_Elec);
            t2(i,j) = temp_matrix;
            t2_new(i,j) = temp_matrix;
            W_4(i,j) = temp_matrix;
            DEN_ABIJ(i,j) = temp_matrix;
        }

        for (int j = i+1; j < unocc_orb; j++)
        {
            temp_matrix = zeros(n_Electrons, Speed_Elec);
            tau4(i,j) = temp_matrix;
        }
    }

    // Coupled Cluster is an iterative process, meaning we make an initial guess of t1 and t2 and calculate the energy
    // Then update t1 and t2 and check for convergance in the energy

    // Find denominator matrix, this should not change during our calculations

    for (int a = n_Electrons; a < Matrix_Size; a++)
    {
        for (int i = 0; i < n_Electrons; i++)
        {
            denom_ai(a-n_Electrons,i) = fs(i,i) - fs(a,a);
            for (int b = n_Electrons; b < Matrix_Size; b++)
            {
                for (int j = 0; j < n_Electrons; j++)
                {
                    denom_abij.at(a-n_Electrons)(b-n_Electrons,i,j) = fs(i,i) + fs(j,j) - fs(a,a) - fs(b,b);

                    if ((a+b+i+j)%2 == 0)
                    {
                        DEN_ABIJ(a-n_Electrons, i)((b-n_Electrons)/2 + (b-n_Electrons)%2 * Speed_Occ, j/2) = 1/(fs(i,i) + fs(j,j) - fs(a,a) - fs(b,b));
                    }
                }
            }
        }
    }

    for (int a = 0; a < unocc_orb; a++)
    {
        for (int i = 0; i < n_Electrons; i++)
        {
            DEN_AI(a/2,i/2) = 1/(fs(i,i) - fs(a+n_Electrons,a+n_Electrons));
            i++;
        }
        a++;
    }

    for (int a = 0+1; a < unocc_orb; a++)
    {
        for (int i = 1; i < n_Electrons; i++)
        {
            DEN_AI(a/2+Speed_Occ,i/2) = 1/(fs(i,i) - fs(a+n_Electrons,a+n_Electrons));
            i++;
        }
        a++;
    }


    // Set up (1 - delta_pq) * fs, these calculations are the same throughout our entire calculation and
    // does not need to be recalculated. Since it is 2 dimensional arrays we can store them
    // Storing in this way reduce the number of elements stored by a factor of n_Electrons * Unoccupied_Orbitals
    // And then later we reduce this by a factor of 1/2. These terms are as far as i know always 0 in restricted HF basis CCSD with closed shells :O
    // The reason is that eigenvalues of fock matrix (fs) is on diagonal matrix form and we have 1 - delta()
    // And also accomodates matrix matrix multiplication operations without any mapping

    for (int i = 0; i < n_Electrons; i++)
    {
        for (int j = 0; j < n_Electrons; j++)
        {
            FS_IJ(i/2, j/2) = (1 - EqualFunc(i,j)) * fs(i,j);
            j++;
        }
        i++;
    }

    for (int i = 1; i < n_Electrons; i++)
    {
        for (int j = 1; j < n_Electrons; j++)
        {
            FS_IJ(i/2+Speed_Elec, j/2) = (1 - EqualFunc(i,j)) * fs(i,j);
        }
    }

    for (int a = 0; a < unocc_orb; a++)
    {
        for (int i = 0; i < n_Electrons; i++)
        {
            FS_AI(a/2,i/2) = fs(a+n_Electrons,i);
            i++;
        }

        for (int b = 0; b < unocc_orb; b++)
        {
            FS_AB(a/2,b/2) = (1 - EqualFunc(a,b)) * fs(a+n_Electrons,b+n_Electrons);
            b++;
        }
        a++;
    }

    for (int a = 1; a < unocc_orb; a++)
    {
        for (int i = 0; i < n_Electrons; i++)
        {
            FS_AI(a/2+Speed_Occ, i/2) = fs(a+n_Electrons,i);
            i++;
        }

        for (int b = 0; b < unocc_orb; b++)
        {
            FS_AB(a/2+Speed_Occ, b/2) = (1 - EqualFunc(a,b)) * fs(a+n_Electrons,b+n_Electrons);
        }
        a++;
    }

    fs.clear(); // Dont need fs anymore

    // The initial guess of t1 and t2 is made here, t1 is guessed to be 0
    // Hypercompact form of T2
    for (int a = 0; a < unocc_orb; a++)
    {
        B = 0;
        for (int b = 0; b < unocc_orb; b++)
        {
            for (int i = 0; i < n_Electrons; i++)
            {
                if ((a-b+i)%2 == 0)
                {
                    J = 0;
                    for (int j = 0; j < n_Electrons; j++)
                    {
                        t2(a,i)(B,J) = integ4.at(i,j)(a/2+a%2*Speed_Occ,b/2) / denom_abij.at(a)(b,i,j);
                        j++;
                        J++;
                    }
                }

                else
                {
                    J = 0;
                    for (int j = 1; j < n_Electrons; j++)
                    {
                        t2(a,i)(B,J) = integ4.at(i,j)(a/2+a%2*Speed_Occ,b/2) / denom_abij.at(a)(b,i,j);
                        j++;
                        J++;
                    }
                }
            }
            b++;
            B++;
        }

        B = 0;
        for (int b = 1; b < unocc_orb; b++)
        {
            for (int i = 0; i < n_Electrons; i++)
            {
                if ((a-b+i)%2 == 0)
                {
                    J = 0;
                    for (int j = 0; j < n_Electrons; j++)
                    {
                        t2(a,i)(B+Speed_Occ,J) = integ4.at(i,j)(a/2+a%2*Speed_Occ,b/2) / denom_abij.at(a)(b,i,j);
                        j++;
                        J++;
                    }
                }

                else
                {
                    J = 0;
                    for (int j = 1; j < n_Electrons; j++)
                    {
                        t2(a,i)(B+Speed_Occ,J) = integ4.at(i,j)(a/2+a%2*Speed_Occ,b/2) / denom_abij.at(a)(b,i,j);
                        j++;
                        J++;
                    }
                }
            }
            b++;
            B++;
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
        Fill_W1();
        Fill_W2();
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



    cout << "Energy fra CCSD: " << E_new << endl;
    return E_new+E_HF;
}

int ccsd_even::EqualFunc(int a, int b)
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

// E done
double ccsd_even::Calc_Energy()
{
    double E1 = 0;
    for (int i = 0; i < n_Electrons; i++)
    {
        for (int j = 0; j < n_Electrons; j++)
        {
            E1 += accu(tau3(i,j) % integ4(i,j));
        }
    }
    E1 *= 0.25;
    E1 += accu(FS_AI % T_1);// + 0.25*temp;
    return E1;
}

// W_1 done
void ccsd_even::Fill_W1()
{
    int sped, sped_j;
    // Matrix symmetric with W_1(i,j) = W_1(j,i), the terms on the "diagonal" of i==j will be equal to 0 (<-- !)
    // We actually dont even need to store anything in W_1(j,i) since it is accessed symmetricly later on also
    // This halves our storage requirements. Also compress matrix. This reduce storage by 75% + not stored diagonal.

    for (int i = 0; i < n_Electrons; i++)
    {
        sped = i%2;
        for (int j = i+1; j < n_Electrons; j++)
        {
            sped_j = j%2;
            if ((i+j)%2 == 1)
            {
                for (int k = 0; k < n_Electrons; k++)
                {
                    for (int l = 1; l < n_Electrons; l++)
                    {
                        W_1(i,j)(k/2,l/2) = integ8.at(i,j)(k/2+k%2*Speed_Elec, l/2);
                        W_1(i,j)(k/2,l/2) += accu(integ6.at(k,l)(span(sped * Speed_Occ, sped*Speed_Occ - 1 + Speed_Occ),j/2) % T_1(span(sped*Speed_Occ, sped*Speed_Occ -1 + Speed_Occ), i/2));
                        W_1(i,j)(k/2,l/2) -= accu(integ6.at(k,l)(span(sped_j*Speed_Occ, sped_j*Speed_Occ-1+Speed_Occ), i/2) % T_1(span(sped_j*Speed_Occ, sped_j*Speed_Occ-1+Speed_Occ), j/2));
                        W_1(i,j)(k/2,l/2) += 0.5*accu(integ4.at(k,l) % tau3.at(i,j));
                        l++;
                    }
                    k++;
                }

                for (int k = 1; k < n_Electrons; k++)
                {
                    for (int l = 0; l < n_Electrons; l++)
                    {
                        W_1(i,j)(k/2+Speed_Elec,l/2) = integ8.at(i,j)(k/2+k%2*Speed_Elec, l/2);
                        W_1(i,j)(k/2+Speed_Elec,l/2) += accu(integ6.at(k,l)(span(sped * Speed_Occ, sped*Speed_Occ - 1 + Speed_Occ),j/2) % T_1(span(sped*Speed_Occ, sped*Speed_Occ -1 + Speed_Occ), i/2));
                        W_1(i,j)(k/2+Speed_Elec,l/2) -= accu(integ6.at(k,l)(span(sped_j*Speed_Occ, sped_j*Speed_Occ-1+Speed_Occ), i/2) % T_1(span(sped_j*Speed_Occ, sped_j*Speed_Occ-1+Speed_Occ), j/2));
                        W_1(i,j)(k/2+Speed_Elec,l/2) += 0.5*accu(integ4.at(k,l) % tau3.at(i,j));
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
                        W_1(i,j)(k/2,l/2) = integ8.at(i,j)(k/2+k%2*Speed_Elec, l/2);
                        W_1(i,j)(k/2,l/2) += accu(integ6.at(k,l)(span(sped * Speed_Occ, sped*Speed_Occ - 1 + Speed_Occ),j/2) % T_1(span(sped*Speed_Occ, sped*Speed_Occ -1 + Speed_Occ), i/2));
                        W_1(i,j)(k/2,l/2) -= accu(integ6.at(k,l)(span(sped_j*Speed_Occ, sped_j*Speed_Occ-1+Speed_Occ), i/2) % T_1(span(sped_j*Speed_Occ, sped_j*Speed_Occ-1+Speed_Occ), j/2));
                        W_1(i,j)(k/2,l/2) += 0.5*accu(integ4.at(k,l) % tau3.at(i,j));
                        l++;
                    }
                    k++;
                }

                for (int k = 1; k < n_Electrons; k++)
                {
                    for (int l = 1; l < n_Electrons; l++)
                    {
                        W_1(i,j)(k/2+Speed_Elec,l/2) = integ8.at(i,j)(k/2+k%2*Speed_Elec, l/2);
                        W_1(i,j)(k/2+Speed_Elec,l/2) += accu(integ6.at(k,l)(span(sped * Speed_Occ, sped*Speed_Occ - 1 + Speed_Occ),j/2) % T_1(span(sped*Speed_Occ, sped*Speed_Occ -1 + Speed_Occ), i/2));
                        W_1(i,j)(k/2+Speed_Elec,l/2) -= accu(integ6.at(k,l)(span(sped_j*Speed_Occ, sped_j*Speed_Occ-1+Speed_Occ), i/2) % T_1(span(sped_j*Speed_Occ, sped_j*Speed_Occ-1+Speed_Occ), j/2));
                        W_1(i,j)(k/2+Speed_Elec,l/2) += 0.5*accu(integ4.at(k,l) % tau3.at(i,j));
                        l++;
                    }
                    k++;
                }
            }
        }
    }
}

// W_2 done
void ccsd_even::Fill_W2()
{
    // Matrix symmetric with W_2(i,j) = -W_2(j,i) and always 0 on the diagonal where i == j (<-- !)
    // Since symmetry is used later on also W_2(j,i) can be whatever, we do not need to access this ever
    // This halves our storage requirements, also compress matrix = reduce storage 75% + nothing on diagonal

    for (int i = 0; i < n_Electrons; i++)
    {
        for (int j = i+1; j < n_Electrons; j++)
        {
            W_2.at(i,j) = integ6.at(i,j);

            if ((i+j)%2==0)
            {
                for (int a = 0; a < unocc_orb; a++)
                {
                    for (int m = 0; m < n_Electrons; m++)
                    {
                        if (i%2 != 1 && j%2 != 1)
                        {
                            W_2(i,j)(a/2+a%2*Speed_Occ,m/2) += accu(integ7.at(a,m)(span(i%2*Speed_Occ, i%2*Speed_Occ+Speed_Occ-1), j/2) % T_1(span(i%2*Speed_Occ, i%2*Speed_Occ+Speed_Occ-1),i/2));
                            W_2(i,j)(a/2+a%2*Speed_Occ,m/2) -= accu(integ7.at(a,m)(span(j%2*Speed_Occ, j%2*Speed_Occ+Speed_Occ-1), i/2) % T_1(span(j%2*Speed_Occ, j%2*Speed_Occ+Speed_Occ-1),j/2));
                            W_2(i,j)(a/2+a%2*Speed_Occ,m/2) += 0.5*accu(integ5.at(a,m) % tau3.at(i,j));
                        }
                        m++;
                    }
                    a++;
                }

                for (int a = 1; a < unocc_orb; a++)
                {
                    for (int m = 1; m < n_Electrons; m++)
                    {
                        W_2(i,j)(a/2+a%2*Speed_Occ,m/2) += accu(integ7.at(a,m)(span(i%2*Speed_Occ, i%2*Speed_Occ+Speed_Occ-1), j/2) % T_1(span(i%2*Speed_Occ, i%2*Speed_Occ+Speed_Occ-1),i/2));
                        W_2(i,j)(a/2+a%2*Speed_Occ,m/2) -= accu(integ7.at(a,m)(span(j%2*Speed_Occ, j%2*Speed_Occ+Speed_Occ-1), i/2) % T_1(span(j%2*Speed_Occ, j%2*Speed_Occ+Speed_Occ-1),j/2));
                        W_2(i,j)(a/2+a%2*Speed_Occ,m/2) += 0.5*accu(integ5.at(a,m) % tau3.at(i,j));
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
                        W_2(i,j)(a/2+a%2*Speed_Occ,m/2) += accu(integ7.at(a,m)(span(i%2*Speed_Occ, i%2*Speed_Occ+Speed_Occ-1), j/2) % T_1(span(i%2*Speed_Occ, i%2*Speed_Occ+Speed_Occ-1),i/2));
                        W_2(i,j)(a/2+a%2*Speed_Occ,m/2) -= accu(integ7.at(a,m)(span(j%2*Speed_Occ, j%2*Speed_Occ+Speed_Occ-1), i/2) % T_1(span(j%2*Speed_Occ, j%2*Speed_Occ+Speed_Occ-1),j/2));
                        W_2(i,j)(a/2+a%2*Speed_Occ,m/2) += 0.5*accu(integ5.at(a,m) % tau3.at(i,j));
                        m++;
                    }
                    a++;
                }

                for (int a = 0; a < unocc_orb; a++)
                {
                    for (int m = 1; m < n_Electrons; m++)
                    {
                        W_2(i,j)(a/2+a%2*Speed_Occ,m/2) += accu(integ7.at(a,m)(span(i%2*Speed_Occ, i%2*Speed_Occ+Speed_Occ-1), j/2) % T_1(span(i%2*Speed_Occ, i%2*Speed_Occ+Speed_Occ-1),i/2));
                        W_2(i,j)(a/2+a%2*Speed_Occ,m/2) -= accu(integ7.at(a,m)(span(j%2*Speed_Occ, j%2*Speed_Occ+Speed_Occ-1), i/2) % T_1(span(j%2*Speed_Occ, j%2*Speed_Occ+Speed_Occ-1),j/2));
                        W_2(i,j)(a/2+a%2*Speed_Occ,m/2) += 0.5*accu(integ5.at(a,m) % tau3.at(i,j));
                        m++;
                    }
                    a++;
                }
            }
        }
    }
}

// W_3 done
void ccsd_even::Fill_W3()
{
    // Intersting function here. m and n is symmetric, which we will utilize
    // at the same time as we compress the matrixes to half size as usual
    int E,N;

    for (int m = 0; m < n_Electrons; m++)
    {
        for (int i = 0; i < n_Electrons; i++)
        {
            if ((m+i)%2 == 0)
            {
                E = 0;
                for (int e = 0; e < unocc_orb; e++)
                {
                    N = 0;
                    for (int n = 0; n < n_Electrons; n++)
                    {
                        W_3(i,m)(E,N) = integ6.at(m,n)(e/2+e%2*Speed_Occ, i/2);
                        W_3(i,m)(E,N) -= accu(integ4.at(m,n)(span(i%2*Speed_Occ, i%2*Speed_Occ + Speed_Occ-1), e/2) % T_1(span(i%2*Speed_Occ, i%2*Speed_Occ+Speed_Occ-1), i/2));
                        N++;
                        n++;
                    }
                    E++;
                    e++;
                }

                E = 0;
                for (int e = 1; e < unocc_orb; e++)
                {
                    N = 0;
                    for (int n = 1; n < n_Electrons; n++)
                    {
                        W_3(i,m)(E+Speed_Occ,N) = integ6.at(m,n)(e/2+e%2*Speed_Occ, i/2);
                        W_3(i,m)(E+Speed_Occ,N) -= accu(integ4.at(m,n)(span(i%2*Speed_Occ, i%2*Speed_Occ + Speed_Occ-1), e/2) % T_1(span(i%2*Speed_Occ, i%2*Speed_Occ+Speed_Occ-1), i/2));
                        N++;
                        n++;
                    }
                    E++;
                    e++;
                }
            }

            else
            {
                E = 0;
                for (int e = 0; e < unocc_orb; e++)
                {
                    N = 0;
                    for (int n = 1; n < n_Electrons; n++)
                    {
                        W_3(i,m)(E,N) = integ6.at(m,n)(e/2+e%2*Speed_Occ, i/2);
                        W_3(i,m)(E,N) -= accu(integ4.at(m,n)(span(i%2*Speed_Occ, i%2*Speed_Occ + Speed_Occ-1), e/2) % T_1(span(i%2*Speed_Occ, i%2*Speed_Occ+Speed_Occ-1), i/2));
                        N++;
                        n++;
                    }
                    E++;
                    e++;
                }

                E = 0;
                for (int e = 1; e < unocc_orb; e++)
                {
                    N = 0;
                    for (int n = 0; n < n_Electrons; n++)
                    {
                        W_3(i,m)(E+Speed_Occ,N) = integ6.at(m,n)(e/2+e%2*Speed_Occ, i/2);
                        W_3(i,m)(E+Speed_Occ,N) -= accu(integ4.at(m,n)(span(i%2*Speed_Occ, i%2*Speed_Occ + Speed_Occ-1), e/2) % T_1(span(i%2*Speed_Occ, i%2*Speed_Occ+Speed_Occ-1), i/2));
                        N++;
                        n++;
                    }
                    E++;
                    e++;
                }
            }
        }
    }
}

// W_4 done
void ccsd_even::Fill_W4()
{
    int M, E;

    // We only need to calculate the terms of W_4 of which t2 amplitudes is not equal to 0 later on
    // since W_4 only appear one place, in which it is multiplied by t2 on a accumulative basis

    for (int a = 0; a < unocc_orb; a++)
    {
        for (int i = 0; i < n_Electrons; i++)
        {
            E = 0;
            for (int e = 0; e < unocc_orb; e++)
            {
                if ((a-e+i)%2 == 0)
                {
                    M = 0;
                    for (int m = 0; m < n_Electrons; m++)
                    {
                        W_4.at(a,i).at(E,M) = -integ7.at(a,m)(e/2,i/2);
                        W_4.at(a,i).at(E,M) -= accu(W_3.at(i,m)(e/2+e%2*Speed_Occ,span()) % T_1.row(a/2+a%2*Speed_Occ));
                        W_4.at(a,i).at(E,M) += accu(integ5.at(a,m)(span(i%2*Speed_Occ, i%2*Speed_Occ+Speed_Occ-1),e/2) % T_1(span(i%2*Speed_Occ, i%2*Speed_Occ+Speed_Occ-1),i/2));
                        W_4.at(a,i).at(E,M) += 0.5*accu(integ2.at(e,m) % t2.at(a,i));
                        m++;
                        M++;
                    }
                }

                else
                {

                    M = 0;
                    for (int m = 1; m < n_Electrons; m++)
                    {
                        W_4.at(a,i).at(E,M) = -integ7.at(a,m)(e/2,i/2);
                        W_4.at(a,i).at(E,M) -= accu(W_3.at(i,m)(e/2+e%2*Speed_Occ,span()) % T_1.row(a/2+a%2*Speed_Occ));
                        W_4.at(a,i).at(E,M) += accu(integ5.at(a,m)(span(i%2*Speed_Occ, i%2*Speed_Occ+Speed_Occ-1),e/2) % T_1(span(i%2*Speed_Occ, i%2*Speed_Occ+Speed_Occ-1),i/2));
                        W_4.at(a,i).at(E,M) += 0.5*accu(integ2.at(e,m) % t2.at(a,i));
                        m++;
                        M++;
                    }
                }
                e++;
                E++;
            }
        }

        E = 0;
        for (int e = 1; e < unocc_orb; e++)
        {
            for (int i = 0; i < n_Electrons; i++)
            {
                if ((a-e+i)%2 == 0)
                {
                    M = 0;
                    for (int m = 0; m < n_Electrons; m++)
                    {
                        W_4.at(a,i).at(E+Speed_Occ,M) = -integ7.at(a,m)(e/2+Speed_Occ,i/2);
                        W_4.at(a,i).at(E+Speed_Occ,M) -= accu(W_3.at(i,m)(e/2+e%2*Speed_Occ,span()) % T_1.row(a/2+a%2*Speed_Occ));
                        W_4.at(a,i).at(E+Speed_Occ,M) += accu(integ5.at(a,m)(span(i%2*Speed_Occ, i%2*Speed_Occ+Speed_Occ-1),e/2) % T_1(span(i%2*Speed_Occ, i%2*Speed_Occ+Speed_Occ-1),i/2));
                        W_4.at(a,i).at(E+Speed_Occ,M) += 0.5*accu(integ2.at(e,m) % t2.at(a,i));
                        m++;
                        M++;
                    }
                }

                else
                {
                    M = 0;
                    for (int m = 1; m < n_Electrons; m++)
                    {
                        W_4.at(a,i).at(E+Speed_Occ,M) = -integ7.at(a,m)(e/2+Speed_Occ,i/2);
                        W_4.at(a,i).at(E+Speed_Occ,M) -= accu(W_3.at(i,m)(e/2+e%2*Speed_Occ,span()) % T_1.row(a/2+a%2*Speed_Occ));
                        W_4.at(a,i).at(E+Speed_Occ,M) += accu(integ5.at(a,m)(span(i%2*Speed_Occ, i%2*Speed_Occ+Speed_Occ-1),e/2) % T_1(span(i%2*Speed_Occ, i%2*Speed_Occ+Speed_Occ-1),i/2));
                        W_4.at(a,i).at(E+Speed_Occ,M) += 0.5*accu(integ2.at(e,m) % t2.at(a,i));
                        m++;
                        M++;
                    }
                }
            }
            e++;
            E++;
        }
    }
}

// F1 done
void ccsd_even::Fill_F1()
{
    D1 = FS_AI;

    int A,M;
    M = 0;
    for (int m = 0; m < n_Electrons; m++)
    {
        A=0;
        for (int a = 0; a < unocc_orb; a++)
        {
            D1.at(A,M) += accu(integ2.at(a,m) % T_1);
            a++;
            A++;
        }
        m++;
        M++;
    }

    M=0;
    for (int m = 1; m < n_Electrons; m++)
    {
        A=0;
        for (int a = 1; a < unocc_orb; a++)
        {
            D1.at(A+Speed_Occ,M) += accu(integ2.at(a,m) % T_1);
            a++;
            A++;
        }
        m++;
        M++;
    }
}

// F2 done
void ccsd_even::Fill_F2()
{
    int I,M;
    int sped;
    M = 0;
    D2 = FS_IJ;

    for (int m = 0; m < n_Electrons; m++)
    {
        I = 0;
        for (int i = 0; i < n_Electrons; i++)
        {
            D2.at(M,I) += accu(D1(span(0, Speed_Occ-1), m/2) % T_1(span(0, Speed_Occ-1), i/2));

            for (int e = 0; e < unocc_orb; e++)
            {
                D2.at(M,I) += 0.5*accu(integ2.at(e,m) % t2.at(e,i));
            }

            for (int k = 0; k < n_Electrons; k++)
            {
                sped = k%2;
                D2.at(M,I) -= accu(integ6.at(m,k)(span(sped*Speed_Occ, sped * Speed_Occ + Speed_Occ-1), i/2) % T_1(span(sped*Speed_Occ, sped*Speed_Occ+Speed_Occ-1), k/2));
            }
            i++;
            I++;
        }
        m++;
        M++;
    }

    M = 0;
    for (int m = 1; m < n_Electrons; m++)
    {
        I = 0;
        for (int i = 1; i < n_Electrons; i++)
        {
            D2.at(M+Speed_Elec,I) += accu(D1(span(Speed_Occ, unocc_orb-1), m/2) % T_1(span(Speed_Occ, unocc_orb-1), i/2));
            for (int e = 0; e < unocc_orb; e++)
            {
                D2.at(M+Speed_Elec,I) += 0.5*accu(integ2.at(e,m) % t2.at(e,i));
            }
            for (int k = 0; k < n_Electrons; k++)
            {
                sped = k%2;
                D2.at(M,I) -= accu(integ6.at(m,k)(span(sped*Speed_Occ, sped * Speed_Occ + Speed_Occ-1), i/2) % T_1(span(sped*Speed_Occ, sped*Speed_Occ+Speed_Occ-1), k/2));
            }
            i++;
            I++;
        }
        m++;
        M++;
    }
}

// F3 done
void ccsd_even::Fill_F3()
{
    // Sleng på en if (m - e % 2 == 0) her, kan skippe noen kalkulasjoner når du tensor contracter ut leddene

    D3 = FS_AB;
    int E;
    for (int a = 0; a < unocc_orb; a++)
    {
        E = 0;
        for (int e = 0; e < unocc_orb; e++)
        {
            D3.at(a,E) -= accu(D1.row(e/2) % T_1.row(a/2));

            for (int m = 1; m < n_Electrons; m++)
            {
                D3.at(a,E) += accu(integ5.at(a,m)(e/2, span()) % T_1(span(Speed_Occ, unocc_orb-1), m/2).t());
                D3.at(a,E) -= 0.5*accu(integ2.at(e,m)(span(Speed_Occ, unocc_orb-1), span()) % t2.at(a,m)(span(Speed_Occ, unocc_orb-1), span()));
                m++;
            }

            for (int m = 0; m < n_Electrons; m++)
            {
                D3.at(a,E) += accu(integ5.at(a,m)(e/2, span()) % T_1(span(0, Speed_Occ-1), m/2).t());
                D3.at(a,E) -= 0.5*accu(integ2.at(e,m) % t2.at(a,m));
                m++;
            }
            e++;
            E++;
        }
        a++;
    }

    for (int a = 1; a < unocc_orb; a++)
    {
        E = 0;
        for (int e = 1; e < unocc_orb; e++)
        {
            D3.at(a,E) -= accu(D1.row(e/2+Speed_Occ) % T_1.row(a/2+Speed_Occ));

            for (int m = 1; m < n_Electrons; m++)
            {
                D3.at(a,E) += accu(integ5.at(a,m)(e/2+Speed_Occ, span()) % T_1(span(Speed_Occ, unocc_orb-1), m/2).t());
                D3.at(a,E) -= 0.5*accu(integ2.at(e,m) % t2.at(a,m));
                m++;
            }

            for (int m = 0; m < n_Electrons; m++)
            {
                D3.at(a,E) += accu(integ5.at(a,m)(e/2+Speed_Occ, span()) % T_1(span(0, Speed_Occ-1), m/2).t());
                D3.at(a,E) -= 0.5*accu(integ2.at(e,m)(span(0,Speed_Occ-1), span()) % t2.at(a,m)(span(0,Speed_Occ-1), span()));
                m++;
            }
            e++;
            E++;
        }
        a++;
    }
}

// tau is on compact form! However must vectorize multiplications somehow
void ccsd_even::Fill_tau()
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
                        tau4(a,b)(i/2,j/2) = t2.at(a,i)(b/2+b%2*Speed_Occ,j/2);
                        if ((a+i)%2 == 0 && (b+j)%2 == 0)
                        {
                            tau4(a,b)(i/2,j/2) += T_1(a/2+a%2*Speed_Occ, i/2) * T_1(b/2+b%2*Speed_Occ, j/2);
                        }

                        if ((a+j)%2 ==0 && (b+i)%2 == 0)
                        {
                            tau4(a,b)(i/2,j/2) -= T_1(a/2+a%2*Speed_Occ, j/2) * T_1(b/2+b%2*Speed_Occ, i/2);
                        }
                        j++;
                    }
                    i++;
                }

                for (int i = 1; i < n_Electrons; i++)
                {
                    for (int j = 0; j < n_Electrons; j++)
                    {
                        tau4(a,b)(i/2+Speed_Elec,j/2) = t2.at(a,i)(b/2+b%2*Speed_Occ,j/2);
                        if ((a+i)%2 == 0 && (b+j)%2 == 0)
                        {
                            tau4(a,b)(i/2+Speed_Elec,j/2) += T_1(a/2+a%2*Speed_Occ, i/2) * T_1(b/2+b%2*Speed_Occ, j/2);
                        }

                        if ((a+j)%2 ==0 && (b+i)%2 == 0)
                        {
                            tau4(a,b)(i/2+Speed_Elec,j/2) -= T_1(a/2+a%2*Speed_Occ, j/2) * T_1(b/2+b%2*Speed_Occ, i/2);
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
                        tau4(a,b)(i/2, j/2) = t2.at(a,i)(b/2+b%2*Speed_Occ,j/2);
                        if ((a+i)%2 == 0 && (b+j)%2 == 0)
                        {
                            tau4(a,b)(i/2,j/2) += T_1(a/2+a%2*Speed_Occ, i/2) * T_1(b/2+b%2*Speed_Occ, j/2);
                        }

                        if ((a+j)%2 ==0 && (b+i)%2 == 0)
                        {
                            tau4(a,b)(i/2,j/2) -= T_1(a/2+a%2*Speed_Occ, j/2) * T_1(b/2+b%2*Speed_Occ, i/2);
                        }
                        j++;
                    }
                    i++;
                }

                for (int i = 1; i < n_Electrons; i++)
                {
                    for (int j = 1; j < n_Electrons; j++)
                    {
                        tau4(a,b)(i/2+Speed_Elec,j/2) = t2.at(a,i)(b/2+b%2*Speed_Occ,j/2);

                        if ((a+i)%2 == 0 && (b+j)%2 == 0)
                        {
                            tau4(a,b)(i/2+Speed_Elec,j/2) += T_1(a/2+a%2*Speed_Occ, i/2) * T_1(b/2+b%2*Speed_Occ, j/2);
                        }

                        if ((a+j)%2 ==0 && (b+i)%2 == 0)
                        {
                            tau4(a,b)(i/2+Speed_Elec,j/2) -= T_1(a/2+a%2*Speed_Occ, j/2) * T_1(b/2+b%2*Speed_Occ, i/2);
                        }
                        j++;
                    }
                    i++;
                }
            }
        }
    }

    for (int i = 0; i < n_Electrons; i++)
    {
        for (int j = 0; j < n_Electrons; j++)
        {
            if ((i+j)%2 == 1)
            {
                for (int a = 0; a < unocc_orb; a++)
                {
                    for (int b = 1; b < unocc_orb; b++)
                    {
                        tau3(i,j)(a/2,b/2) = t2.at(a,i)(b/2+b%2*Speed_Occ,j/2);

                        if ((a+i)%2 == 0 && (b+j)%2 == 0)
                        {
                            tau3(i,j)(a/2,b/2) += T_1(a/2+a%2*Speed_Occ, i/2) * T_1(b/2+b%2*Speed_Occ, j/2);
                        }

                        if ((a+j)%2 ==0 && (b+i)%2 == 0)
                        {
                            tau3(i,j)(a/2,b/2) -= T_1(a/2+a%2*Speed_Occ, j/2) * T_1(b/2+b%2*Speed_Occ, i/2);
                        }
                        b++;
                    }
                    a++;
                }

                for (int a = 1; a < unocc_orb; a++)
                {
                    for (int b = 0; b < unocc_orb; b++)
                    {
                        tau3(i,j)(a/2+Speed_Occ,b/2) = t2.at(a,i)(b/2+b%2*Speed_Occ,j/2);

                        if ((a+i)%2 == 0 && (b+j)%2 == 0)
                        {
                            tau3(i,j)(a/2+Speed_Occ,b/2) += T_1(a/2+a%2*Speed_Occ, i/2) * T_1(b/2+b%2*Speed_Occ, j/2);
                        }

                        if ((a+j)%2 ==0 && (b+i)%2 == 0)
                        {
                            tau3(i,j)(a/2+Speed_Occ,b/2) -= T_1(a/2+a%2*Speed_Occ, j/2) * T_1(b/2+b%2*Speed_Occ, i/2);
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
                        tau3(i,j)(a/2,b/2) = t2.at(a,i)(b/2+b%2*Speed_Occ,j/2);

                        if ((a+i)%2 == 0 && (b+j)%2 == 0)
                        {
                            tau3(i,j)(a/2,b/2) += T_1(a/2+a%2*Speed_Occ, i/2) * T_1(b/2+b%2*Speed_Occ, j/2);
                        }

                        if ((a+j)%2 ==0 && (b+i)%2 == 0)
                        {
                            tau3(i,j)(a/2,b/2) -= T_1(a/2+a%2*Speed_Occ, j/2) * T_1(b/2+b%2*Speed_Occ, i/2);
                        }
                        b++;
                    }
                    a++;
                }

                for (int a = 1; a < unocc_orb; a++)
                {
                    for (int b = 1; b < unocc_orb; b++)
                    {
                        tau3(i,j)(a/2+Speed_Occ,b/2) = t2.at(a,i)(b/2+b%2*Speed_Occ,j/2);

                        if ((a+i)%2 == 0 && (b+j)%2 == 0)
                        {
                            tau3(i,j)(a/2+Speed_Occ,b/2) += T_1(a/2+a%2*Speed_Occ, i/2) * T_1(b/2+b%2*Speed_Occ, j/2);
                        }

                        if ((a+j)%2 ==0 && (b+i)%2 == 0)
                        {
                            tau3(i,j)(a/2+Speed_Occ,b/2) -= T_1(a/2+a%2*Speed_Occ, j/2) * T_1(b/2+b%2*Speed_Occ, i/2);
                        }

                        b++;
                    }
                    a++;
                }
            }
        }
       tau3(i,i) = zeros(unocc_orb, Speed_Occ);
    }
}

// T1 done
void ccsd_even::Fill_t1_new()
{
    // Optimal version of T1 amplitudes:
    int A=0, I;
    double temp;

    // f_ai
    T_1_new = FS_AI;

    // Both even numbers != 0
    for (int a = 0; a < unocc_orb; a++)
    {
        I = 0;
        for (int i = 0; i < n_Electrons; i++)
        {
            // t_ik^ac [F_1]_c^k
            T_1_new(A,I) += accu(D1 % t2(a,i));

            // - t_k^a [F_2]_i^k
            T_1_new(A,I) -= accu(D2(span(i%2*Speed_Elec, i%2*Speed_Elec+Speed_Elec-1), i/2) % T_1.row(a/2+a%2*Speed_Occ).t());

            // f_ac t_i^c
            T_1_new(A,I) += accu(FS_AB(A, span()) % T_1(span(0, Speed_Occ-1), I).t());

            // I_ka^ci t_k^c
            T_1_new(A,I) -= accu(integ10(a,i) % T_1);

            // 1/2 I_ka^cd tau_ij^cd
            temp = 0;
            for (int k = 0; k < n_Electrons; k++)
            {
                temp += accu(integ5.at(a,k) % tau3(i,k));
            }
            T_1_new(A,I) += 0.5*temp;

            // - 1/2 I_kl^ci t_kl^ca - 1/2 I_kl^cd t_kl^ca t_i^d
            temp = 0;
            for (int k = 0; k < n_Electrons; k++)
            {
                temp += accu(W_3(i,k) % t2(a,k));
            }
            T_1_new(A,I) += 0.5*temp;

            i++;
            I++;
        }
        a++;
        A++;
    }

    // Both odd numbers also != 0
    for (int a = 1; a < unocc_orb; a++)
    {
        I = 0;
        for (int i = 1; i < n_Electrons; i++)
        {
            // t_ik^ac [F_1]_c^k
            T_1_new(A,I) += accu(D1 % t2(a,i));

            // - t_k^a [F_2]_i^k
            T_1_new(A,I) -= accu(D2(span(i%2*Speed_Elec, i%2*Speed_Elec+Speed_Elec-1), i/2) % T_1.row(a/2+a%2*Speed_Occ).t());

            // f_ac t_i^c
            T_1_new(A,I) += accu(FS_AB(A, span()) % T_1(span(0, Speed_Occ-1), I).t());

            // I_ka^ci t_k^c
            T_1_new(A,I) -= accu(integ10(a,i) % T_1);

            // 1/2 I_ka^cd tau_ij^cd
            temp = 0;
            for (int k = 0; k < n_Electrons; k++)
            {
                temp += accu(integ5.at(a,k) % tau3(i,k));
            }
            T_1_new(A,I) += 0.5*temp;

            // - 1/2 I_kl^ci t_kl^ca - 1/2 I_kl^cd t_kl^ca t_i^d
            temp = 0;
            for (int k = 0; k < n_Electrons; k++)
            {
                temp += accu(W_3(i,k) % t2(a,k));
            }
            T_1_new(A,I) += 0.5*temp;

            i++;
            I++;
        }
        a++;
        A++;
    }
    T_1_new = T_1_new % DEN_AI;
}

// T2 almost done, slightly less efficient than with T_2 kept as is
void ccsd_even::Fill_t2_new()
{
    // -0.0501273

    // T2 amplitudes calculated here
    int J;

    // One symmetry that holds and can be proven computationally (by printing the numbers) is that
    // if T2(a,b,i,j) should be not equal to 0 then the following must hold:
    // if (a - b + i) is an even number then j must also be an even number
    // if (a - b + i) is an odd number then j must also be an odd number
    // Everything else is 0. This is implemented here. (a-b+i) = 0 is counted as an even number,
    // also if (a-b+i) is a negative odd or even number then j must be an odd or even number, but positive

    // Also j > i is the only thing worth calculating since symmetry holds. The odd,even thing holds threw symmetry also

    // This function is implemented slightly inefficient
    // However the speedup will all be made up for when program writes new program by itself



    for (int a = 0; a < unocc_orb; a++)
    {
        // a != b and i != j, we have meny symmetries t2(a,b) = -t2(b,a), t2(i,j) = -t2(j,i), and the combination



        // Here a - b is always an odd number:
        for (int b = a+1; b < unocc_orb; b++)
        {
            // Then i + j must be an odd number also:
            for (int i = 0; i < n_Electrons; i++)
            {
                J = (int) (i+1)/2;
                for (int j = i+1; j < n_Electrons; j++)
                {
                    if ((a+i+b+j)%2 == 0)
                    {
                        // I_ab^ij
                        t2_new(a,i)(b/2+b%2*Speed_Occ, j/2) = integ4.at(i,j)(a/2+a%2*Speed_Occ,b/2);

                        // P(ab) P(ij) [W_4] t_jk^bc, ARMADILLO
                        t2_new(a,i)(b/2+b%2*Speed_Occ, j/2) += accu(t2.at(b,j) % W_4.at(a,i));
                        t2_new(a,i)(b/2+b%2*Speed_Occ, j/2) -= accu(t2.at(b,i) % W_4.at(a,j));
                        t2_new(a,i)(b/2+b%2*Speed_Occ, j/2) -= accu(t2.at(a,j) % W_4.at(b,i));
                        t2_new(a,i)(b/2+b%2*Speed_Occ, j/2) += accu(t2.at(a,i) % W_4.at(b,j));

                        // - P(ab) [W_2] t_k^b, ARMADILLO
                        t2_new(a,i)(b/2+b%2*Speed_Occ, j/2) -= accu(W_2.at(i,j)(a/2+a%2*Speed_Occ, span()) % T_1.row(b/2+b%2*Speed_Occ));
                        t2_new(a,i)(b/2+b%2*Speed_Occ, j/2) += accu(W_2.at(i,j)(b/2+b%2*Speed_Occ, span()) % T_1.row(a/2+a%2*Speed_Occ));

                        // 1/2 I_ab^cd tau_ij^cd, ARMADILLO
                        t2_new(a,i)(b/2+b%2*Speed_Occ, j/2) += 0.5*accu(integ3.at(a,b) % tau3.at(i,j));

                        // P(ij) I_ab^cj t_i^c, ARMADILLO
                        t2_new(a,i)(b/2+b%2*Speed_Occ, j/2) -= accu(integ9.at(a,b)(span(j%2*Speed_Occ, j%2*Speed_Occ+Speed_Occ-1), i/2) % T_1(span(j%2*Speed_Occ, j%2*Speed_Occ+Speed_Occ-1), j/2));
                        t2_new(a,i)(b/2+b%2*Speed_Occ, j/2) += accu(integ9.at(a,b)(span(i%2*Speed_Occ, i%2*Speed_Occ+Speed_Occ-1), j/2) % T_1(span(i%2*Speed_Occ, i%2*Speed_Occ+Speed_Occ-1), i/2));

                        // 1/2 [W_1] tau_kl^ab, ARMADILLO
                        t2_new(a,i)(b/2+b%2*Speed_Occ, j/2) += 0.5*accu(W_1.at(i,j) % tau4.at(a,b));

                        // P(ij) t_jk^ab [F_2] + P(ab) t_ij^bc [F_3]_c^a, ARMADILLO
                        if (a % 2 == 0)
                        {
                            if (i % 2 == 0)
                            {
                                t2_new(a,i)(b/2+b%2*Speed_Occ, j/2) += accu(t2.at(a,j)((int)b/2+Speed_Occ, span()) % D2.row(i/2));
                                t2_new(a,i)(b/2+b%2*Speed_Occ, j/2) -= accu(t2.at(a,i)((int)b/2+Speed_Occ, span()) % D2.row((int)(j/2+Speed_Elec)));
                            }

                            else
                            {
                                t2_new(a,i)(b/2+b%2*Speed_Occ, j/2) += accu(t2.at(a,j)((int)b/2+Speed_Occ, span()) % D2.row((int)i/2+Speed_Elec));
                                t2_new(a,i)(b/2+b%2*Speed_Occ, j/2) -= accu(t2.at(a,i)((int)b/2+Speed_Occ, span()) % D2.row((int)j/2));
                            }

                            t2_new(a,i)(b/2+b%2*Speed_Occ, j/2) -= accu(t2.at(b,i)(span(0, Speed_Occ-1), J) % D3.row(a).t());
                            t2_new(a,i)(b/2+b%2*Speed_Occ, j/2) += accu(t2.at(a,i)(span(Speed_Occ, unocc_orb-1), J) % D3.row(b).t());
                        }

                        else
                        {
                            if (i % 2 == 0)
                            {
                                t2_new(a,i)(b/2+b%2*Speed_Occ, j/2) += accu(t2.at(a,j)(b/2,span()) % D2.row(i/2));
                                t2_new(a,i)(b/2+b%2*Speed_Occ, j/2) -= accu(t2.at(a,i)((int)b/2, span()) % D2.row(j/2+Speed_Elec));
                            }

                            else
                            {
                                t2_new(a,i)(b/2+b%2*Speed_Occ, j/2) += accu(t2.at(a,j)((int)b/2, span()) % D2.row(i/2+Speed_Elec));
                                t2_new(a,i)(b/2+b%2*Speed_Occ, j/2) -= accu(t2.at(a,i)((int)b/2, span()) % D2.row(j/2));
                            }

                            t2_new(a,i)(b/2+b%2*Speed_Occ, j/2) -= accu(t2.at(b,i)(span(Speed_Occ, unocc_orb-1), J) % D3.row(a).t());
                            t2_new(a,i)(b/2+b%2*Speed_Occ, j/2) += accu(t2.at(a,i)(span(0, Speed_Occ-1), J) % D3.row(b).t());
                        }

                        // Symmetries
                        t2_new(b,i)(a/2+a%2*Speed_Occ, j/2) = -t2_new(a,i)(b/2+b%2*Speed_Occ, j/2);
                        t2_new(a,j)(b/2+b%2*Speed_Occ, i/2) = t2_new(b,i)(a/2+a%2*Speed_Occ, j/2);
                        t2_new(b,j)(a/2+a%2*Speed_Occ, i/2) = t2_new(a,i)(b/2+b%2*Speed_Occ, j/2);
                    }

                    j++;
                    J++;
                }
            }
            b++;
        }

        // Here a - b is always an even number
        for (int b = a+2; b < unocc_orb; b++)
        {
            // Then i-j must be an even number also
            for (int i = 0; i < n_Electrons; i++)
            {
                J = (int) (i/2);
                J++;
                for (int j = i+2; j < n_Electrons; j++)
                {
                    if (a%2 != 0 || b%2 != 0 || i%2 != 1 || j%2 != 1)
                    {
                    if ((a+i+b+j)%2 == 0)
                    {
                        // I_ab^ij
                        t2_new(a,i)(b/2+b%2*Speed_Occ, j/2) = integ4.at(i,j)(a/2+a%2*Speed_Occ,b/2);

                        // P(ab) P(ij) [W_4] t_jk^bc, ARMADILLO
                        t2_new(a,i)(b/2+b%2*Speed_Occ, j/2) += accu(t2.at(b,j) % W_4.at(a,i));
                        t2_new(a,i)(b/2+b%2*Speed_Occ, j/2) -= accu(t2.at(b,i) % W_4.at(a,j));
                        t2_new(a,i)(b/2+b%2*Speed_Occ, j/2) -= accu(t2.at(a,j) % W_4.at(b,i));
                        t2_new(a,i)(b/2+b%2*Speed_Occ, j/2) += accu(t2.at(a,i) % W_4.at(b,j));

                        // - P(ab) [] t_k^b, ARMADILLO
                        t2_new(a,i)(b/2+b%2*Speed_Occ, j/2) -= accu(W_2.at(i,j)(a/2+a%2*Speed_Occ, span()) % T_1.row(b/2+b%2*Speed_Occ));
                        t2_new(a,i)(b/2+b%2*Speed_Occ, j/2) += accu(W_2.at(i,j)(b/2+b%2*Speed_Occ, span()) % T_1.row(a/2+a%2*Speed_Occ));

                        // 1/2 I_ab^cd tau_ij^cd, ARMADILLO
                        t2_new(a,i)(b/2+b%2*Speed_Occ, j/2) += 0.5*accu(integ3.at(a,b) % tau3.at(i,j));

                        // P(ij) I_ab^cj t_i^c, ARMADILLO
                        t2_new(a,i)(b/2+b%2*Speed_Occ, j/2) -= accu(integ9.at(a,b)(span(j%2*Speed_Occ, j%2*Speed_Occ+Speed_Occ-1), i/2) % T_1(span(j%2*Speed_Occ, j%2*Speed_Occ+Speed_Occ-1), j/2));
                        t2_new(a,i)(b/2+b%2*Speed_Occ, j/2) += accu(integ9.at(a,b)(span(i%2*Speed_Occ, i%2*Speed_Occ+Speed_Occ-1), j/2) % T_1(span(i%2*Speed_Occ, i%2*Speed_Occ+Speed_Occ-1), i/2));

                        // 1/2 [W_1] tau_kl^ab, ARMADILLO
                        t2_new(a,i)(b/2+b%2*Speed_Occ, j/2) += 0.5*accu(W_1.at(i,j) % tau4.at(a,b));

                        // P(ij) t_jk^ab [F_2] + P(ab) t_ij^bc [F_3]_c^a, ARMADILLO
                        if (a % 2 == 0)
                        {
                            t2_new(a,i)(b/2+b%2*Speed_Occ, j/2) -= accu(t2.at(b,i)(span(0, Speed_Occ-1), J) % D3.row(a).t());
                            t2_new(a,i)(b/2+b%2*Speed_Occ, j/2) += accu(t2.at(a,i)(span(0, Speed_Occ-1), J) % D3.row(b).t());

                            if (i % 2 == 0)
                            {
                                t2_new(a,i)(b/2+b%2*Speed_Occ, j/2) += accu(t2.at(a,j)((int)b/2, span()) % D2.row(i/2));
                                t2_new(a,i)(b/2+b%2*Speed_Occ, j/2) -= accu(t2.at(a,i)((int)b/2, span()) % D2.row(j/2+Speed_Elec));
                            }

                            else
                            {
                                t2_new(a,i)(b/2+b%2*Speed_Occ, j/2) += accu(t2.at(a,j)((int)b/2, span()) % D2.row((int)i/2+Speed_Elec));
                                t2_new(a,i)(b/2+b%2*Speed_Occ, j/2) -= accu(t2.at(a,i)((int)b/2, span()) % D2.row(j/2));
                            }
                        }

                        else
                        {
                            t2_new(a,i)(b/2+b%2*Speed_Occ, j/2) -= accu(t2.at(b,i)(span(j%2*Speed_Occ, j%2*Speed_Occ+Speed_Occ-1), J) % D3.row(a).t());
                            t2_new(a,i)(b/2+b%2*Speed_Occ, j/2) += accu(t2.at(a,i)(span(j%2*Speed_Occ, j%2*Speed_Occ+Speed_Occ-1), J) % D3.row(b).t());


                            if (i % 2 == 0)
                            {
                                t2_new(a,i)(b/2+b%2*Speed_Occ, j/2) += accu(t2.at(a,j)(b/2+Speed_Occ,span()) % D2.row(i/2));
                                t2_new(a,i)(b/2+b%2*Speed_Occ, j/2) -= accu(t2.at(a,i)((int)b/2+Speed_Occ, span()) % D2.row(j/2+Speed_Elec));
                            }

                            else
                            {
                                t2_new(a,i)(b/2+b%2*Speed_Occ, j/2) += accu(t2.at(a,j)((int)b/2+Speed_Occ, span()) % D2.row(i/2+Speed_Elec));
                                t2_new(a,i)(b/2+b%2*Speed_Occ, j/2) -= accu(t2.at(a,i)((int)b/2+Speed_Occ, span()) % D2.row(j/2));
                            }

                        }

                        // Symmetries
                        t2_new(b,i)(a/2+a%2*Speed_Occ, j/2) = -t2_new(a,i)(b/2+b%2*Speed_Occ, j/2);
                        t2_new(a,j)(b/2+b%2*Speed_Occ, i/2) = t2_new(b,i)(a/2+a%2*Speed_Occ, j/2);
                        t2_new(b,j)(a/2+a%2*Speed_Occ, i/2) = t2_new(a,i)(b/2+b%2*Speed_Occ, j/2);
                    }
                    }

                    j++;
                    J++;
                }
            }
            b++;
        }
    }

    // Get denominator matrix

    for (int a = 0; a < unocc_orb; a++)
    {
        for (int i = 0; i < n_Electrons; i++)
        {
            t2_new(a,i) = t2_new(a,i) % DEN_ABIJ(a,i);
        }
    }
}
