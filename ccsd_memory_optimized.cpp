#include "ccsd_memory_optimized.h"

CCSD_Memory_optimized::CCSD_Memory_optimized(int n_N, vec zz, mat rr, string B_s, int n_Elec)
{
    n_Nuclei = n_N;
    Z = zz;
    R = rr;
    Basis_Set = B_s;
    n_Electrons = n_Elec;
}

double CCSD_Memory_optimized::CCSD(double toler, bool print_stuff)
{
    /*
     *
     * *Optimeringstanker:
     * - W3 er symmetrisk på en måte som er inni if testene i kommentar i koden. Finn ut hvordan skippe beregningene du ikke trenger
     * - W4 symmetrisk?
     * - W3(0,1) er 0 i øvre halvdel, trenger ikke lagres
     * - W4(0,1) er 0 i nedre halvdel, trenger ikke lagres
     *
     *  - Sjekk om ALLE kalkulasjoner = 0 for bestemte indekser... :O og fjern dette minnet
     *
     *  - integ4, integ6 kan lagres for i > j, [integ(i,j)]
     *  - integ10 kan lagres for (a+i)%2=0... ? :O
     *
     * - W1 = integ8 ikke innefor for loopen
     *
     * */

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
    ccint.Prepear_Integrals(); // Most of the 2 elekctron integrals contained here
    fs = ccint.Return_FS(HartFock.return_eigval_F()); // eigenvalues of fock matrix
    HartFock.Delete_Everything();

    // We now pull out the parts of the MO basis we need, we split up the this 4D array which is a
    // n^4 array, n^4 = (n_u + n_o)^4, where n_u = unoccupied orbitals in HF basis and n_o is
    // occupied orbitals in HF basis. To the right of each term is the nr of elements stored inside it
    integ2 = ccint.Return_Rearranged_Integrals2();    // 1/2 n_o^2 n_u^2 |
    integ3 = ccint.Return_Rearranged_Integrals3();    // 1/4 n_u^4       |
    integ4 = ccint.Return_Rearranged_Integrals4();    // 1/2 n_o^2 n_u^2 |
    integ5 = ccint.Return_Rearranged_Integrals5();    // 1/2 n_u^3 n_o   |
    integ6 = ccint.Return_Rearranged_Integrlas6();    // 1/2 n_o^3 n_u   |
    integ7 = ccint.Return_Rearranged_Integrals7();    // 1/2 n_u^2 n_o^2 |
    integ8 = ccint.Return_Rearranged_Integrals8();    // 1/4 n_o^4       |
    integ9 = ccint.Return_Rearranged_Integrals9();    // 1/4 n_u^3 n_o   |
    integ10 = ccint.Return_Rearranged_Integrals10();  // 1/2 n_u^2 n_o^2 |

    /*
     * Storing in this manner is
     * 1/4 n_u^4 + 1/4 n_o^4 + n_u^3 n_o + 1/2 n_o^3 n_u + 2 n_o^2 n_u^2
     *
     * Which is a decrease in storage of >50 % compared to storing everything, exact percentage depends on system
     *
     * combined with all intermediates, we store:
     *
     *   1/4 n_u^4
     * + 1/4 n_o^4
     * + 3/4 n_u^3 n_o
     * +     n_o^3 n_u
     * + 4/2 n_o^2 n_u^2
     * + 5/2 n_o^2 n_u^2
     *
     * + bunch of 2D arrays which is approx 0 in comparison
     *
     * Current bottleneck in transfer from coupled_cluster_integrals to here, since we store 2 (or 3) times the nr of integrals,
     * and the intermediates right now are much smaller than the integrals. Possible non-storage of integrals next up,
     * with storing the AO integrals (scales 1/8 * (1/2 n_o + 1/2 n_u)^4) would remove this bottleneck, but recalculation of
     * integX's would be required in some way
     *
     * Further improvements in storage of intermediates uneccasary untill bottleneck removed.
     *
     */

    // Delete all variables inside that ccint object,
    // hopefully this is helpful before we create our variables here
    ccint.Delete_Everything();

    // Define variables used in the CCSD method:

    // Optimized 4D arrays
    tau4.set_size(unocc_orb, unocc_orb); // 1/2 n_o^2 n_u^2
    t2.set_size(unocc_orb, n_Electrons); // 1/2 n_o^2 n_u^2
    t2_new.set_size(unocc_orb, n_Electrons); // 1/2 n_o^2 n_u^2
    W_3.set_size(n_Electrons, n_Electrons); // 1/2 n_o^3 n_u
    W_4.set_size(unocc_orb, n_Electrons); // 1/2 n_o^2 n_u^2
    DEN_ABIJ.set_size(unocc_orb, n_Electrons); // 1/2 n_o^2 n_u^2

    // Optimized 2D arrays
    W_1 = zeros(n_Electrons, Speed_Elec); // This is 4D array but can be stored as 2D if updated inside T2 amplitude loops
    W_2 = zeros(unocc_orb, Speed_Elec); // This is 4D array but can be stored as 2D if updated inside T2 amplitude loops
    D3 = zeros(unocc_orb, Speed_Occ);
    D2 = zeros(n_Electrons, Speed_Elec);
    T_1 = zeros(unocc_orb, Speed_Elec);
    T_1_new = zeros(unocc_orb, Speed_Elec);
    D1 = zeros(unocc_orb, Speed_Elec);
    FS_AI = zeros(unocc_orb, Speed_Elec);
    FS_AB = zeros(unocc_orb, Speed_Occ);
    FS_IJ = zeros(n_Electrons, Speed_Elec);
    DEN_AI = zeros(unocc_orb, Speed_Elec);
    tau1 = zeros(unocc_orb, Speed_Occ); // This is a 2D mapping of tau :-O

    // Initialize our 4D arrays
    for (int i = 0; i < n_Electrons; i++)
    {
        for (int j = 0; j < n_Electrons; j++)
        {
            temp_matrix = zeros(unocc_orb, Speed_Elec);
            W_3(i,j) = temp_matrix;
        }
    }

    for (int i = 0; i < unocc_orb; i++)
    {
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

    // Find denominator matrix, this should not change during our calculations, stored as denom^(-1) for later matrix multiplication use

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


    // Code generator up till here, also need to add AO2MO functions (<-- !)

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
            j++;
        }
        i++;
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
        for (int i = 1; i < n_Electrons; i++)
        {
            FS_AI(a/2+Speed_Occ, i/2) = fs(a+n_Electrons,i);
            i++;
        }

        for (int b = 1; b < unocc_orb; b++)
        {
            FS_AB(a/2+Speed_Occ, b/2) = (1 - EqualFunc(a,b)) * fs(a+n_Electrons,b+n_Electrons);
            b++;
        }
        a++;
    }

    fs.clear(); // Dont need fs anymore
    // Continued storage of 2D array fs makes available recalculation of DEN_ABIJ each iteration and removes required storage of this 4D array

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
                        t2(a,i)(B,J) = integ4.at(i,j)(a/2+a%2*Speed_Occ,b/2) * DEN_ABIJ(a,i)(B,J);
                        j++;
                        J++;
                    }
                }

                else
                {
                    J = 0;
                    for (int j = 1; j < n_Electrons; j++)
                    {
                        t2(a,i)(B,J) = integ4.at(i,j)(a/2+a%2*Speed_Occ,b/2) * DEN_ABIJ(a,i)(B,J);
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
                        t2(a,i)(B+Speed_Occ,J) = integ4.at(i,j)(a/2+a%2*Speed_Occ,b/2) * DEN_ABIJ(a,i)(B+Speed_Occ,J);
                        j++;
                        J++;
                    }
                }

                else
                {
                    J = 0;
                    for (int j = 1; j < n_Electrons; j++)
                    {
                        t2(a,i)(B+Speed_Occ,J) = integ4.at(i,j)(a/2+a%2*Speed_Occ,b/2) * DEN_ABIJ(a,i)(B+Speed_Occ,J);
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
        //E_new = 0;

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

        // Find energy, calc energy function has been put into Fill_t2_new() for more efficient storage of variables for now
        E_new = Calc_Energy();
        //E_new *= 0.25;
        //E_new += accu(FS_AI % T_1);
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


    //cout << W_3(0,1) << endl; // Nedre halvdel 0
    //cout << W_4(0,1) << endl; // Øvre halvdel 0
    //cout << integ5(1,1) << endl; // Øvre halvdel 0
    //cout << integ7(1,1) << endl; // Øvre halvdel 0
    //cout << integ6(1,1) << endl; // Øvre halvdel 0


    cout << "CCSD correction [au] = " << E_new << endl;
    return E_new+E_HF;
}

int CCSD_Memory_optimized::EqualFunc(int a, int b)
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

double CCSD_Memory_optimized::Calc_Energy()
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
    E1 += accu(FS_AI % T_1); // This function only called once, when T1 is 0 in the optimized version to avoid recalculation of tau1
    return E1;
}

void CCSD_Memory_optimized::Fill_W1(int i, int j)
{
    int sped, sped_j;
    // Matrix symmetric with W_1(i,j) = W_1(j,i), the terms on the "diagonal" of i==j will be equal to 0 (<-- !)
    // We actually dont even need to store anything in W_1(j,i) since it is accessed symmetricly later on also
    // This halves our storage requirements. Also compress matrix. This reduce storage by 75% + not stored diagonal.

    sped = i%2;
    sped_j = j%2;
    W_1 = integ8(i,j);

            if ((i+j)%2 == 1)
            {

                for (int k = 0; k < n_Electrons; k++)
                {
                    for (int l = 1; l < n_Electrons; l++)
                    {
                        W_1(k/2,l/2) += accu(integ6.at(k,l)(span(sped * Speed_Occ, sped*Speed_Occ - 1 + Speed_Occ),j/2) % T_1(span(sped*Speed_Occ, sped*Speed_Occ -1 + Speed_Occ), i/2));
                        W_1(k/2,l/2) -= accu(integ6.at(k,l)(span(sped_j*Speed_Occ, sped_j*Speed_Occ-1+Speed_Occ), i/2) % T_1(span(sped_j*Speed_Occ, sped_j*Speed_Occ-1+Speed_Occ), j/2));
                        W_1(k/2,l/2) += 0.5*accu(integ4.at(k,l) % tau1);
                        l++;
                    }
                    k++;
                }

                for (int k = 1; k < n_Electrons; k++)
                {
                    for (int l = 0; l < n_Electrons; l++)
                    {
                        W_1(k/2+Speed_Elec,l/2) += accu(integ6.at(k,l)(span(sped * Speed_Occ, sped*Speed_Occ - 1 + Speed_Occ),j/2) % T_1(span(sped*Speed_Occ, sped*Speed_Occ -1 + Speed_Occ), i/2));
                        W_1(k/2+Speed_Elec,l/2) -= accu(integ6.at(k,l)(span(sped_j*Speed_Occ, sped_j*Speed_Occ-1+Speed_Occ), i/2) % T_1(span(sped_j*Speed_Occ, sped_j*Speed_Occ-1+Speed_Occ), j/2));
                        W_1(k/2+Speed_Elec,l/2) += 0.5*accu(integ4.at(k,l) % tau1);
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
                        W_1(k/2,l/2) += accu(integ6.at(k,l)(span(sped * Speed_Occ, sped*Speed_Occ - 1 + Speed_Occ),j/2) % T_1(span(sped*Speed_Occ, sped*Speed_Occ -1 + Speed_Occ), i/2));
                        W_1(k/2,l/2) -= accu(integ6.at(k,l)(span(sped_j*Speed_Occ, sped_j*Speed_Occ-1+Speed_Occ), i/2) % T_1(span(sped_j*Speed_Occ, sped_j*Speed_Occ-1+Speed_Occ), j/2));
                        W_1(k/2,l/2) += 0.5*accu(integ4.at(k,l) % tau1);
                        l++;
                    }
                    k++;
                }

                for (int k = 1; k < n_Electrons; k++)
                {
                    for (int l = 1; l < n_Electrons; l++)
                    {
                        W_1(k/2+Speed_Elec,l/2) += accu(integ6.at(k,l)(span(sped * Speed_Occ, sped*Speed_Occ - 1 + Speed_Occ),j/2) % T_1(span(sped*Speed_Occ, sped*Speed_Occ -1 + Speed_Occ), i/2));
                        W_1(k/2+Speed_Elec,l/2) -= accu(integ6.at(k,l)(span(sped_j*Speed_Occ, sped_j*Speed_Occ-1+Speed_Occ), i/2) % T_1(span(sped_j*Speed_Occ, sped_j*Speed_Occ-1+Speed_Occ), j/2));
                        W_1(k/2+Speed_Elec,l/2) += 0.5*accu(integ4.at(k,l) % tau1);

                        //if (accu(integ6.at(k,l)(span(sped * Speed_Occ, sped*Speed_Occ - 1 + Speed_Occ),j/2) % T_1(span(sped*Speed_Occ, sped*Speed_Occ -1 + Speed_Occ), i/2)) == 0)
                        //{
                            //cout << "sup?" << endl;
                        //}

                        l++;
                    }
                    k++;
                }
            }
}

void CCSD_Memory_optimized::Fill_W2(int i, int j)
{
    // Matrix symmetric with W_2(i,j) = -W_2(j,i) and always 0 on the diagonal where i == j (<-- !)
    // Since symmetry is used later on also W_2(j,i) can be whatever, we do not need to access this ever
    // This halves our storage requirements, also compress matrix = reduce storage 75% + nothing on diagonal

            W_2 = integ6.at(i,j);

            if ((i+j)%2==0)
            {
                for (int a = 0; a < unocc_orb; a++)
                {
                    for (int m = 0; m < n_Electrons; m++)
                    {
                        if (i%2 != 1 && j%2 != 1)
                        {
                            W_2(a/2+a%2*Speed_Occ,m/2) += accu(integ7.at(a,m)(span(i%2*Speed_Occ, i%2*Speed_Occ+Speed_Occ-1), j/2) % T_1(span(i%2*Speed_Occ, i%2*Speed_Occ+Speed_Occ-1),i/2));
                            W_2(a/2+a%2*Speed_Occ,m/2) -= accu(integ7.at(a,m)(span(j%2*Speed_Occ, j%2*Speed_Occ+Speed_Occ-1), i/2) % T_1(span(j%2*Speed_Occ, j%2*Speed_Occ+Speed_Occ-1),j/2));
                            W_2(a/2+a%2*Speed_Occ,m/2) += 0.5*accu(integ5.at(a,m) % tau1);
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
                        W_2(a/2+a%2*Speed_Occ,m/2) += accu(integ7.at(a,m)(span(i%2*Speed_Occ, i%2*Speed_Occ+Speed_Occ-1), j/2) % T_1(span(i%2*Speed_Occ, i%2*Speed_Occ+Speed_Occ-1),i/2));
                        W_2(a/2+a%2*Speed_Occ,m/2) -= accu(integ7.at(a,m)(span(j%2*Speed_Occ, j%2*Speed_Occ+Speed_Occ-1), i/2) % T_1(span(j%2*Speed_Occ, j%2*Speed_Occ+Speed_Occ-1),j/2));
                        W_2(a/2+a%2*Speed_Occ,m/2) += 0.5*accu(integ5.at(a,m) % tau1);
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
                        W_2(a/2+a%2*Speed_Occ,m/2) += accu(integ7.at(a,m)(span(i%2*Speed_Occ, i%2*Speed_Occ+Speed_Occ-1), j/2) % T_1(span(i%2*Speed_Occ, i%2*Speed_Occ+Speed_Occ-1),i/2));
                        W_2(a/2+a%2*Speed_Occ,m/2) -= accu(integ7.at(a,m)(span(j%2*Speed_Occ, j%2*Speed_Occ+Speed_Occ-1), i/2) % T_1(span(j%2*Speed_Occ, j%2*Speed_Occ+Speed_Occ-1),j/2));
                        W_2(a/2+a%2*Speed_Occ,m/2) += 0.5*accu(integ5.at(a,m) % tau1);
                        m++;
                    }
                    a++;
                }

                for (int a = 0; a < unocc_orb; a++)
                {
                    for (int m = 1; m < n_Electrons; m++)
                    {
                        W_2(a/2+a%2*Speed_Occ,m/2) += accu(integ7.at(a,m)(span(i%2*Speed_Occ, i%2*Speed_Occ+Speed_Occ-1), j/2) % T_1(span(i%2*Speed_Occ, i%2*Speed_Occ+Speed_Occ-1),i/2));
                        W_2(a/2+a%2*Speed_Occ,m/2) -= accu(integ7.at(a,m)(span(j%2*Speed_Occ, j%2*Speed_Occ+Speed_Occ-1), i/2) % T_1(span(j%2*Speed_Occ, j%2*Speed_Occ+Speed_Occ-1),j/2));
                        W_2(a/2+a%2*Speed_Occ,m/2) += 0.5*accu(integ5.at(a,m) % tau1);
                        m++;
                    }
                    a++;
                }
            }
}

void CCSD_Memory_optimized::Fill_W3()
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
                        //if (n != m)
                        //{
                        //    W_3(i,n)(E,m/2) = -W_3(i,m)(E,N);
                        //}
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
                        //if (n != m)
                        //{
                        //    W_3(i,n)(E+Speed_Occ,m/2) = -W_3(i,m)(E+Speed_Occ,N);
                        //}
                        N++;
                        n++;
                    }
                    E++;
                    e++;
                }
            }

            else
            {
                if (i%2 != 0 && m%2 != 1)
                {
                    E = 0;
                    for (int e = 0; e < unocc_orb; e++)
                    {
                        N = 0;
                        for (int n = 1; n < n_Electrons; n++)
                        {
                            W_3(i,m)(E,N) = integ6.at(m,n)(e/2+e%2*Speed_Occ, i/2);
                            W_3(i,m)(E,N) -= accu(integ4.at(m,n)(span(i%2*Speed_Occ, i%2*Speed_Occ + Speed_Occ-1), e/2) % T_1(span(i%2*Speed_Occ, i%2*Speed_Occ+Speed_Occ-1), i/2));
                            //if (n != m)
                            //{
                            //    W_3(i,n)(E,m/2) = -W_3(i,m)(E,N);
                            //}
                            N++;
                            n++;
                        }
                        E++;
                        e++;
                    }
                }

                if (i%2 != 1 && m%2 != 0)
                {
                    E = 0;
                    for (int e = 1; e < unocc_orb; e++)
                    {
                        N = 0;
                        for (int n = 0; n < n_Electrons; n++)
                        {
                            W_3(i,m)(E+Speed_Occ,N) = integ6.at(m,n)(e/2+e%2*Speed_Occ, i/2);
                            W_3(i,m)(E+Speed_Occ,N) -= accu(integ4.at(m,n)(span(i%2*Speed_Occ, i%2*Speed_Occ + Speed_Occ-1), e/2) % T_1(span(i%2*Speed_Occ, i%2*Speed_Occ+Speed_Occ-1), i/2));
                            //if (n != m)
                            //{
                            //    W_3(i,n)(E+Speed_Occ,m/2) = -W_3(i,m)(E+Speed_Occ,N);
                            //}
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
}

void CCSD_Memory_optimized::Fill_W4()
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

                        if (i%2 != 0)
                        {
                            W_4.at(a,i).at(E,M) -= accu(W_3.at(i,m)(e/2+e%2*Speed_Occ,span()) % T_1.row(a/2+a%2*Speed_Occ));
                        }

                        if (a%2 != 1)
                        {
                            W_4.at(a,i).at(E,M) += accu(integ5.at(a,m)(span(i%2*Speed_Occ, i%2*Speed_Occ+Speed_Occ-1),e/2) % T_1(span(i%2*Speed_Occ, i%2*Speed_Occ+Speed_Occ-1),i/2));
                        }
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
                        if (a%2 != 0 && i%2 != 1)
                        {
                            W_4.at(a,i).at(E+Speed_Occ,M) = -integ7.at(a,m)(e/2+Speed_Occ,i/2);
                            W_4.at(a,i).at(E+Speed_Occ,M) -= accu(W_3.at(i,m)(e/2+e%2*Speed_Occ,span()) % T_1.row(a/2+a%2*Speed_Occ));
                            W_4.at(a,i).at(E+Speed_Occ,M) += accu(integ5.at(a,m)(span(i%2*Speed_Occ, i%2*Speed_Occ+Speed_Occ-1),e/2) % T_1(span(i%2*Speed_Occ, i%2*Speed_Occ+Speed_Occ-1),i/2));
                            W_4.at(a,i).at(E+Speed_Occ,M) += 0.5*accu(integ2.at(e,m) % t2.at(a,i));
                        }
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

void CCSD_Memory_optimized::Fill_F1()
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

void CCSD_Memory_optimized::Fill_F2()
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

void CCSD_Memory_optimized::Fill_F3()
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

void CCSD_Memory_optimized::Fill_tau()
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
}

void CCSD_Memory_optimized::Fill_2D_tau(int i, int j)
{
    if ((i+j)%2 == 1)
    {
        for (int a = 0; a < unocc_orb; a++)
        {
            for (int b = 1; b < unocc_orb; b++)
            {
                //tau1(a/2, b/2) = tau4(a,b)(i/2+i%2*Speed_Elec,j/2);
                tau1(a/2,b/2) = t2.at(a,i)(b/2+b%2*Speed_Occ,j/2);

                if ((a+i)%2 == 0 && (b+j)%2 == 0)
                {
                    tau1(a/2,b/2) += T_1(a/2+a%2*Speed_Occ, i/2) * T_1(b/2+b%2*Speed_Occ, j/2);
                }

                if ((a+j)%2 ==0 && (b+i)%2 == 0)
                {
                    tau1(a/2,b/2) -= T_1(a/2+a%2*Speed_Occ, j/2) * T_1(b/2+b%2*Speed_Occ, i/2);
                }
                b++;
            }
            a++;
        }

        for (int a = 1; a < unocc_orb; a++)
        {
            for (int b = 0; b < unocc_orb; b++)
            {
                tau1(a/2+Speed_Occ,b/2) = t2.at(a,i)(b/2+b%2*Speed_Occ,j/2);

                if ((a+i)%2 == 0 && (b+j)%2 == 0)
                {
                    tau1(a/2+Speed_Occ,b/2) += T_1(a/2+a%2*Speed_Occ, i/2) * T_1(b/2+b%2*Speed_Occ, j/2);
                }

                if ((a+j)%2 ==0 && (b+i)%2 == 0)
                {
                    tau1(a/2+Speed_Occ,b/2) -= T_1(a/2+a%2*Speed_Occ, j/2) * T_1(b/2+b%2*Speed_Occ, i/2);
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
                tau1(a/2,b/2) = t2.at(a,i)(b/2+b%2*Speed_Occ,j/2);

                if ((a+i)%2 == 0 && (b+j)%2 == 0)
                {
                    tau1(a/2,b/2) += T_1(a/2+a%2*Speed_Occ, i/2) * T_1(b/2+b%2*Speed_Occ, j/2);
                }

                if ((a+j)%2 ==0 && (b+i)%2 == 0)
                {
                    tau1(a/2,b/2) -= T_1(a/2+a%2*Speed_Occ, j/2) * T_1(b/2+b%2*Speed_Occ, i/2);
                }
                b++;
            }
            a++;
        }

        for (int a = 1; a < unocc_orb; a++)
        {
            for (int b = 1; b < unocc_orb; b++)
            {
                tau1(a/2+Speed_Occ,b/2) = t2.at(a,i)(b/2+b%2*Speed_Occ,j/2);

                if ((a+i)%2 == 0 && (b+j)%2 == 0)
                {
                    tau1(a/2+Speed_Occ,b/2) += T_1(a/2+a%2*Speed_Occ, i/2) * T_1(b/2+b%2*Speed_Occ, j/2);
                }

                if ((a+j)%2 ==0 && (b+i)%2 == 0)
                {
                    tau1(a/2+Speed_Occ,b/2) -= T_1(a/2+a%2*Speed_Occ, j/2) * T_1(b/2+b%2*Speed_Occ, i/2);
                }

                b++;
            }
            a++;
        }
    }
}

void CCSD_Memory_optimized::Fill_t1_new()
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
            for (int k = 0; k < i; k++)
            {
                Fill_2D_tau(k,i);
                temp -= accu(integ5.at(a,k) % tau1); // Double calculation here? Rest is in Fill_t2()
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
            for (int k = 0; k < i; k++)
            {
                Fill_2D_tau(k,i); // Double calculation here? Rest is in Fill_t2()
                temp -= accu(integ5.at(a,k) % tau1);
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
}

void CCSD_Memory_optimized::Fill_t2_new()
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

    Fill_tau();

    for (int i = 0; i < n_Electrons; i++)
    {
        J = (int) (i+1)/2;
        for (int j = i+1; j < n_Electrons; j++)
        {
            Fill_2D_tau(i,j);
            Fill_W1(i,j);
            Fill_W2(i,j);
            //E_new += accu(tau1 % integ4(i,j));
            //E_new -= accu(tau1 % integ4(j,i));
            for (int a = 0; a < unocc_orb; a++)
            {
                if ((a+i)%2 == 0)
                {
                    T_1_new(a/2+a%2*Speed_Occ, i/2) += 0.5*accu(integ5.at(a,j) % tau1); // One part of T1 calculated here because tau1 is stored 2D
                }
                for (int b = a+1; b < unocc_orb; b++)
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
                        t2_new(a,i)(b/2+b%2*Speed_Occ, j/2) -= accu(W_2(a/2+a%2*Speed_Occ, span()) % T_1.row(b/2+b%2*Speed_Occ));
                        t2_new(a,i)(b/2+b%2*Speed_Occ, j/2) += accu(W_2(b/2+b%2*Speed_Occ, span()) % T_1.row(a/2+a%2*Speed_Occ));

                        // 1/2 I_ab^cd tau_ij^cd, ARMADILLO
                        t2_new(a,i)(b/2+b%2*Speed_Occ, j/2) += 0.5*accu(integ3.at(a,b) % tau1);

                        // P(ij) I_ab^cj t_i^c, ARMADILLO
                        t2_new(a,i)(b/2+b%2*Speed_Occ, j/2) -= accu(integ9.at(a,b)(span(j%2*Speed_Occ, j%2*Speed_Occ+Speed_Occ-1), i/2) % T_1(span(j%2*Speed_Occ, j%2*Speed_Occ+Speed_Occ-1), j/2));
                        t2_new(a,i)(b/2+b%2*Speed_Occ, j/2) += accu(integ9.at(a,b)(span(i%2*Speed_Occ, i%2*Speed_Occ+Speed_Occ-1), j/2) % T_1(span(i%2*Speed_Occ, i%2*Speed_Occ+Speed_Occ-1), i/2));

                        // 1/2 [W_1] tau_kl^ab, ARMADILLO
                        t2_new(a,i)(b/2+b%2*Speed_Occ, j/2) += 0.5*accu(W_1 % tau4.at(a,b));

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


                    b++;
                }
            }
            j++;
            J++;
        }
    }

    for (int i = 0; i < n_Electrons; i++)
    {
        J = (int) (i/2);
        J++;
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
                    T_1_new(a/2+a%2*Speed_Occ, i/2) += 0.5*accu(integ5.at(a,j) % tau1); // One part of T1 calculated here because tau1 is stored 2D
                }

                for (int b = a+2; b < unocc_orb; b++)
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
                        t2_new(a,i)(b/2+b%2*Speed_Occ, j/2) -= accu(W_2(a/2+a%2*Speed_Occ, span()) % T_1.row(b/2+b%2*Speed_Occ));
                        t2_new(a,i)(b/2+b%2*Speed_Occ, j/2) += accu(W_2(b/2+b%2*Speed_Occ, span()) % T_1.row(a/2+a%2*Speed_Occ));

                        // 1/2 I_ab^cd tau_ij^cd, ARMADILLO
                        t2_new(a,i)(b/2+b%2*Speed_Occ, j/2) += 0.5*accu(integ3.at(a,b) % tau1);

                        // P(ij) I_ab^cj t_i^c, ARMADILLO
                        t2_new(a,i)(b/2+b%2*Speed_Occ, j/2) -= accu(integ9.at(a,b)(span(j%2*Speed_Occ, j%2*Speed_Occ+Speed_Occ-1), i/2) % T_1(span(j%2*Speed_Occ, j%2*Speed_Occ+Speed_Occ-1), j/2));
                        t2_new(a,i)(b/2+b%2*Speed_Occ, j/2) += accu(integ9.at(a,b)(span(i%2*Speed_Occ, i%2*Speed_Occ+Speed_Occ-1), j/2) % T_1(span(i%2*Speed_Occ, i%2*Speed_Occ+Speed_Occ-1), i/2));

                        // 1/2 [W_1] tau_kl^ab, ARMADILLO
                        t2_new(a,i)(b/2+b%2*Speed_Occ, j/2) += 0.5*accu(W_1 % tau4.at(a,b));

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

                    b++;
                }
            }
            j++;
            J++;
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


