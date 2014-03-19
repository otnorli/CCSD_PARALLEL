#include "coupled_cluster_v2.h"

coupled_cluster_v2::coupled_cluster_v2(int n_N, vec zz, mat rr, string B_s, int n_Elec)
{    
    n_Nuclei = n_N;
    Z = zz;
    R = rr;
    Basis_Set = B_s;
    n_Electrons = n_Elec;
}

double coupled_cluster_v2::CCSD(double toler, bool print_stuff)
{
    bool continue_ccsd = true;
    double convergance_check;
    int Itermax = 1000;
    iter = 0;
    E_old = 0;
    E_new = 0;

    // Initializion
    double E_HF;
    Hartree_Fock_Solver HartFock(n_Nuclei, Z, R, Basis_Set, n_Electrons, true);
    E_HF = HartFock.get_Energy(toler); // Calc hartree fock energy
    Matrix_Size = 2*HartFock.ReturnMatrixSize(); // This will be twice as big as in Hatree Fock, because of spin up and down
    unocc_orb = Matrix_Size - n_Electrons; // Number of unocupied orbitals
    mat temp_matrix; cube temp_mat;

    // Defining arrays to store stuff, these are the ones that fill up all the memory
    fs = zeros(Matrix_Size, Matrix_Size);
    t1 = zeros(unocc_orb, n_Electrons);
    t1_new = zeros(unocc_orb, n_Electrons);
    T_2.set_size(unocc_orb, n_Electrons);
    T_2_new.set_size(unocc_orb, n_Electrons);
    tau.set_size(n_Electrons, n_Electrons);
    tau2.set_size(unocc_orb, unocc_orb);
    D1 = zeros(unocc_orb, n_Electrons);
    D2 = zeros(n_Electrons, n_Electrons);
    D3 = zeros(unocc_orb, unocc_orb);
    W_1.set_size(n_Electrons, n_Electrons);
    W_2.set_size(unocc_orb, n_Electrons);
    W_3.set_size(unocc_orb, n_Electrons);
    W_4.set_size(unocc_orb, n_Electrons);
    denom_ai = zeros(unocc_orb, n_Electrons);

    t2.set_size(unocc_orb, n_Electrons);

    // Initialize our arrays
    for (int i = 0; i < n_Electrons; i++)
    {
        for (int j = 0; j < n_Electrons; j++)
        {
            temp_matrix = zeros(n_Electrons, n_Electrons);
            W_1(i,j) = temp_matrix;

            temp_matrix = zeros(unocc_orb, unocc_orb);
            tau(i,j) = temp_matrix;
        }
    }

    for (int i = 0; i < unocc_orb; i++)
    {
        temp_mat = zeros(unocc_orb, n_Electrons, n_Electrons);
        denom_abij.push_back(temp_mat);

        for (int j = 0; j < n_Electrons; j++)
        {
            temp_matrix = zeros(unocc_orb, n_Electrons);
            W_4(i,j) = temp_matrix;
            T_2(i,j) = temp_matrix;
            T_2_new(i,j) = temp_matrix;

            temp_matrix = zeros(n_Electrons, n_Electrons);
            W_2(i,j) = temp_matrix;
            W_3(i,j) = temp_matrix;
        }

        for (int j = 0; j < unocc_orb; j++)
        {
            temp_matrix = zeros(n_Electrons, n_Electrons);
            tau2(i,j) = temp_matrix;
        }
    }

    // Transform to MO basis
    coupled_cluster_integrals ccint(Matrix_Size/2, HartFock.ReturnC(), n_Electrons);
    //integ = ccint.Return_Integrals(HartFock.ReturnQ()); // 2 elektron integralene her
    integ2 = ccint.Return_Integrals_Part_2(); // Splitting up 2 elektron integrals for easier matrix multiplication use later on
    fs = ccint.Return_FS(HartFock.return_eigval_F()); // alt annet her

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
                }
            }
        }
    }

    // The initial guess of t1 and t2 is made here, t1 is guessed to be 0
    for (int a = 0; a < unocc_orb; a++)
    {
        //for (int b = n_Electrons; b < Matrix_Size; b++)
        for (int b = 0; b < unocc_orb; b++)
        {
            for (int i = 0; i < n_Electrons; i++)
            {
                for (int j = 0; j < n_Electrons; j++)
                {
                    T_2.at(a,i)(b,j) = integ.at(i,j)(a+n_Electrons,b+n_Electrons) / denom_abij.at(a)(b,i,j);
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
        t1 = t1_new;
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

int coupled_cluster_v2::EqualFunc(int a, int b)
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

double coupled_cluster_v2::Calc_Energy()
{
    double E1=0;
    int temp_a;
    for (int a = 0; a < unocc_orb; a++)
    {
        temp_a = a + n_Electrons;
        for (int i = 0; i < n_Electrons; i++)
        {
            E1 += fs(i,temp_a) * t1(a,i);
            E1 += 0.25*accu(integ2.at(a,i) % T_2_new.at(a,i));

            for (int b = 0; b < unocc_orb; b++)
            {
                E1 += 0.5*t1.at(a,i) * accu(integ2.at(a,i).row(b) % t1.row(b));
            }
        }
    }
    return E1;
}

void coupled_cluster_v2::Fill_W1()
{
    // Matrix symmetric with W_1(i,j) = W_1(j,i), also SOME of the terms on the "diagonal" of i==j will be equal to 0 (<-- !)
    // We actually dont even need to store anything in W_1(j,i) since it is accessed symmetricly later on also

    for (int i = 0; i < n_Electrons; i++)
    {
        for (int l = 0; l < n_Electrons; l++)
        {
            for (int k = 0; k < n_Electrons; k++)
            {
                for (int j = i+1; j < n_Electrons; j++)
                {
                    W_1.at(i,j).at(k,l) = integ.at(k,l).at(i,j);
                    W_1.at(i,j).at(k,l) += accu(integ.at(k,l)(span(n_Electrons, Matrix_Size-1),j) % t1.col(i));
                    W_1.at(i,j).at(k,l) -= accu(integ.at(k,l)(span(n_Electrons, Matrix_Size-1),i) % t1.col(j));
                    W_1.at(i,j).at(k,l) += 0.5*accu((integ.at(k,l).submat(span(n_Electrons, Matrix_Size-1),span(n_Electrons, Matrix_Size-1))) % tau.at(i,j));
                }

                // "Diagonal" i==j
                W_1.at(i,i).at(k,l) = integ.at(k,l).at(i,i);
                W_1.at(i,i).at(k,l) += 0.5*accu((integ.at(k,l).submat(span(n_Electrons, Matrix_Size-1),span(n_Electrons, Matrix_Size-1))) % tau.at(i,i));
            }
        }
    }
}

void coupled_cluster_v2::Fill_W2()
{
    // Matrix symmetric with W_2(i,j) = -W_2(j,i) and always 0 on the diagonal where i == j (<-- !)
    // Since symmetry is used later on also W_2(j,i) can be whatever, we do not need to access this ever

    int temp_a;
    for (int a = 0; a < unocc_orb; a++)
    {
        temp_a = a + n_Electrons;
        for (int m = 0; m < n_Electrons; m++)
        {
            for (int i = 0; i < n_Electrons; i++)
            {
                for (int j = i+1; j < n_Electrons; j++)
                {
                    W_2.at(a,i).at(j,m) = integ.at(temp_a,m).at(i,j);
                    W_2.at(a,i).at(j,m) += accu(integ.at(temp_a,m)(span(n_Electrons, Matrix_Size-1),j) % t1.col(i));
                    W_2.at(a,i).at(j,m) -= accu(integ.at(temp_a,m)(span(n_Electrons, Matrix_Size-1),i) % t1.col(j));
                    W_2.at(a,i).at(j,m) += 0.5*accu((integ.at(temp_a,m).submat(span(n_Electrons, Matrix_Size-1), span(n_Electrons, Matrix_Size-1))) % tau.at(i,j));
                }
            }
        }
    }
}

void coupled_cluster_v2::Fill_W3()
{
    // Matrix symmetric with W_3(m,n) = -W_3(n,m) however not 0 at diagonal (<-- !)

    for (int e = 0; e < unocc_orb; e++)
    {
        for (int i = 0; i < n_Electrons; i++)
        {
            for (int n = 0; n < n_Electrons; n++)
            {
                for (int m = n+1; m < n_Electrons; m++)
                {
                    W_3.at(e,i).at(m,n) = integ.at(m,n).at(e+n_Electrons,i);
                    W_3.at(e,i).at(m,n) -= accu(integ.at(m,n)(span(n_Electrons, Matrix_Size-1), e+n_Electrons) % t1.col(i));
                    W_3.at(e,i).at(n,m) = -W_3.at(e,i).at(m,n);
                }
                W_3.at(e,i).at(n,n) = integ.at(n,n).at(e+n_Electrons,i);
                W_3.at(e,i).at(n,n) -= accu(integ.at(n,n)(span(n_Electrons, Matrix_Size-1), e+n_Electrons) % t1.col(i));
            }
        }
    }
}

void coupled_cluster_v2::Fill_W4()
{
    int temp_a, temp_e;
    for (int a = 0; a < unocc_orb; a++)
    {
        temp_a = a + n_Electrons;
        for (int e = 0; e < unocc_orb; e++)
        {
            temp_e = e + n_Electrons;
            for (int i = 0; i < n_Electrons; i++)
            {
                for (int m = 0; m < n_Electrons; m++)
                {
                    W_4.at(a,i).at(e,m) = integ.at(temp_a,m).at(i,temp_e);
                    W_4.at(a,i).at(e,m) -= accu(W_3.at(e,i)(m, span()) % t1.row(a));
                    W_4.at(a,i).at(e,m) += accu(integ.at(temp_a,m)(span(n_Electrons, Matrix_Size-1), temp_e) % t1.col(i));
                    W_4.at(a,i).at(e,m) += 0.5*accu(integ2.at(e,m) % T_2.at(a,i));
                }
            }
        }
    }
}

void coupled_cluster_v2::Fill_F1()
{
    int temp_a;
    for (int a = 0; a < unocc_orb; a++)
    {
        temp_a = a + n_Electrons;
        for (int m = 0; m < n_Electrons; m++)
        {
            D1.at(a,m) = fs.at(temp_a,m);
            D1.at(a,m) += accu(integ2.at(a, m) % t1);
        }
    }
}

void coupled_cluster_v2::Fill_F2()
{
    for (int m = 0; m < n_Electrons; m++)
    {
        for (int i = 0; i < n_Electrons; i++)
        {
            D2.at(m,i) = (1 - EqualFunc(m,i)) * fs(m,i);
            D2.at(m,i) += accu(D1.col(m) % t1.col(i));
            for (int e = 0; e < unocc_orb; e++)
            {
                D2.at(m,i) -= accu(integ.at(e+n_Electrons,i)(m, span(0,n_Electrons-1)) % t1.row(e));
                D2.at(m,i) += 0.5 * accu(integ2.at(e,m) % T_2.at(e,i));
            }
        }
    }
}

void coupled_cluster_v2::Fill_F3()
{
    int temp_a, temp_e;
    for (int a = 0; a < unocc_orb; a++)
    {
        temp_a = a + n_Electrons;
        for (int e = 0; e < unocc_orb; e++)
        {
            temp_e = e + n_Electrons;
            D3(a,e) = (1 - EqualFunc(a,e)) * fs(temp_a, temp_e);
            D3.at(a,e) -= accu(D1.row(e) % t1.row(a));
            for (int m = 0; m < n_Electrons; m++)
            {
                D3.at(a,e) += accu(integ.at(temp_a,m)(temp_e, span(n_Electrons, Matrix_Size-1)) % t1.col(m).t());
                D3.at(a,e) -= 0.5*accu(integ2.at(e,m) % T_2.at(a,m));
            }
        }
    }
}

void coupled_cluster_v2::Fill_tau()
{
    for (int a = 0; a < unocc_orb; a++)
    {
        for (int b = 0; b < unocc_orb; b++)
        {
            for (int i = 0; i < n_Electrons; i++)
            {
                for (int j = 0; j < n_Electrons; j++)
                {
                    tau.at(i,j)(a,b) = T_2_new.at(a,i)(b,j) + t1(a,i)*t1(b,j) - t1(a,j)*t1(b,i);
                    tau2.at(a,b)(i,j) = tau.at(i,j)(a,b);
                }
            }
        }
    }
}

void coupled_cluster_v2::Fill_t1_new()
{
    // T1 is NOT optimized jet

    int temp_a;
    for (int a = 0; a < unocc_orb; a++)
    {
        temp_a = a + n_Electrons;
        for (int i = 0; i < n_Electrons; i++)
        {
            // f_ai
            t1_new(a,i) = fs(temp_a,i);

            // t_ik^ac [F_1]_c^k
            for (int k = 0; k < n_Electrons; k++)
            {
                for (int e = n_Electrons; e < Matrix_Size; e++)
                {
                    t1_new(a,i) += D1(e-n_Electrons,k) * T_2.at(a,i)(e-n_Electrons,k);
                }
            }

            // - t_k^a [F_2]_i^k
            for (int k = 0; k < n_Electrons; k++)
            {
                t1_new(a,i) -= D2(k,i) * t1(a,k);
            }

            // I_ka^ci t_k^c
            for (int k = 0; k < n_Electrons; k++)
            {
                for (int e = n_Electrons; e < Matrix_Size; e++)
                {
                    t1_new(a,i) += integ.at(e,i)(k, temp_a) * t1(e-n_Electrons,k);
                }
            }

            // 1/2 I_ka^cd t_ki^cd + I_ka^cd t_k^c t_i^d
            for (int k = 0; k < n_Electrons; k++)
            {
                for (int e = n_Electrons; e < Matrix_Size; e++)
                {
                    for (int d = n_Electrons; d < Matrix_Size; d++)
                    {
                        t1_new(a,i) -= integ.at(temp_a,k)(e,d) * 0.5 * T_2.at(e-n_Electrons,k)(d-n_Electrons,i);
                        t1_new(a,i) -= integ.at(temp_a,k)(e,d) * t1(e-n_Electrons,k) * t1(d-n_Electrons,i);
                    }
                }
            }

            // - 1/2 I_kl^ci t_kl^ca - 1/2 I_kl^cd t_kl^ca t_i^d
            for (int k = 0; k < n_Electrons; k++)
            {
                for (int e = n_Electrons; e < Matrix_Size; e++)
                {
                    for (int l = 0; l < n_Electrons; l++)
                    {
                        t1_new(a,i) -= 0.5 * integ.at(k,l)(e,i) * T_2.at(e-n_Electrons,k)(a,l);
                        for (int d = n_Electrons; d < Matrix_Size; d++)
                        {
                            t1_new(a,i) -= 0.5 * integ.at(k,l)(e,d) * T_2.at(e-n_Electrons, k)(a,l) * t1(d-n_Electrons,i);
                        }
                    }
                }
            }

        // f_ac t_i^c
        for (int e = n_Electrons; e < Matrix_Size; e++)
        {
            t1_new(a,i) += (1 - EqualFunc(temp_a, e)) * fs(temp_a, e) * t1(e-n_Electrons,i);
        }

        t1_new(a,i) = t1_new(a,i)/denom_ai(a,i);
        }
    }
}

void coupled_cluster_v2::Fill_t2_new()
{
    // T2 amplitudes calculated here
    int temp_a, temp_b;

    // One symmetry that holds and can be proven computationally (by printing the numbers) is that
    // if T2(a,b,i,j) should be not equal to 0 then the following must hold:
    // if (a - b + i) is an even number then j must also be an even number
    // if (a - b + i) is an odd number then j must also be an odd number
    // Everything else is 0. This is implemented here. (a-b+i) = 0 is counted as an even number,
    // also if (a-b+i) is a negative odd or even number then j must be an odd or even number, but positive

    // Also j > i is the only thing worth calculating since symmetry holds. The odd,even thing holds threw symmetry also

    // This function currently calculates a few to meny numbers when a - b and i - j are EVEN numbers

    for (int a = 0; a < unocc_orb; a++)
    {
        temp_a = a + n_Electrons;

        // a != b and i != j, we have meny symmetries t2(a,b) = -t2(b,a), t2(i,j) = -t2(j,i), and the combination

        // Here a - b is always an odd number:
        for (int b = a+1; b < unocc_orb; b++)
        {
            temp_b = b + n_Electrons;

            // Then i + j must be an odd number also:
            for (int i = 0; i < n_Electrons; i++)
            {
                for (int j = i+1; j < n_Electrons; j++)
                {
                    // I_ab^ij
                    T_2_new.at(a,i).at(b,j) = integ.at(i,j)(temp_a,temp_b);

                    // P(ab) P(ij) [W_4] t_jk^bc, ARMADILLO
                    T_2_new.at(a,i).at(b,j) += accu(T_2.at(b,j) % W_4.at(a,i));
                    T_2_new.at(a,i).at(b,j) -= accu(T_2.at(b,i) % W_4.at(a,j));
                    T_2_new.at(a,i).at(b,j) -= accu(T_2.at(a,j) % W_4.at(b,i));
                    T_2_new.at(a,i).at(b,j) += accu(T_2.at(a,i) % W_4.at(b,j));

                    // - P(ab) [W_2] t_k^b, ARMADILLO
                    T_2_new.at(a,i).at(b,j) -= accu(W_2.at(a,i)(j, span(0, n_Electrons-1)) % t1.row(b));
                    T_2_new.at(a,i).at(b,j) += accu(W_2.at(b,i)(j, span(0, n_Electrons-1)) % t1.row(a));

                    // 1/2 I_ab^cd tau_ij^cd, ARMADILLO
                    T_2_new.at(a,i).at(b,j) += 0.5 * accu(integ.at(temp_a, temp_b).submat(span(n_Electrons, Matrix_Size-1), span(n_Electrons, Matrix_Size-1)) % tau.at(i,j));

                    // P(ij) I_ab^cj t_i^c, ARMADILLO
                    T_2_new.at(a,i).at(b,j) -= accu(integ.at(temp_a, temp_b)(span(n_Electrons, Matrix_Size-1), i) % t1.col(j));
                    T_2_new.at(a,i).at(b,j) += accu(integ.at(temp_a, temp_b)(span(n_Electrons, Matrix_Size-1), j) % t1.col(i));


                    // 1/2 [W_1] tau_kl^ab, ARMADILLO
                    T_2_new.at(a,i).at(b,j) += 0.5 * accu(W_1.at(i,j) % tau2.at(a,b));

                    // P(ij) t_jk^ab [F_2], ARMADILLO
                    T_2_new.at(a,i).at(b,j) += accu(T_2.at(a,j)(b, span(0,n_Electrons-1)) % D2.row(i));
                    T_2_new.at(a,i).at(b,j) -= accu(T_2.at(a,i)(b, span(0,n_Electrons-1)) % D2.row(j));

                    // P(ab) t_ij^bc [F_3]_c^a, ARMADILLO
                    T_2_new.at(a,i).at(b,j) -= accu(T_2.at(b,i)(span(), j) % D3.row(a).t());
                    T_2_new.at(a,i).at(b,j) += accu(T_2.at(a,i)(span(), j) % D3.row(b).t());

                    // Divided by Denominator
                    T_2_new.at(a,i).at(b,j) = T_2_new.at(a,i).at(b,j) / denom_abij.at(a)(b,i,j);

                    // Symmetries
                    T_2_new.at(b,i).at(a,j) = -T_2_new.at(a,i).at(b,j);
                    T_2_new.at(a,j).at(b,i) = T_2_new.at(b,i).at(a,j);
                    T_2_new.at(b,j).at(a,i) = T_2_new.at(a,i).at(b,j);
                    j++;
                }
            }
            b++;
        }

        // Here a - b is always an even number
        for (int b = a+2; b < unocc_orb; b++)
        {
            temp_b = b + n_Electrons;

            // Then i-j must be an even number also
            for (int i = 0; i < n_Electrons; i++)
            {
                for (int j = i+2; j < n_Electrons; j++)
                {
                    // I_ab^ij
                    T_2_new.at(a,i).at(b,j) = integ.at(i,j)(temp_a,temp_b);

                    // P(ab) P(ij) [W_4] t_jk^bc, ARMADILLO
                    T_2_new.at(a,i).at(b,j) += accu(T_2.at(b,j) % W_4.at(a,i));
                    T_2_new.at(a,i).at(b,j) -= accu(T_2.at(b,i) % W_4.at(a,j));
                    T_2_new.at(a,i).at(b,j) -= accu(T_2.at(a,j) % W_4.at(b,i));
                    T_2_new.at(a,i).at(b,j) += accu(T_2.at(a,i) % W_4.at(b,j));

                    // - P(ab) [W_2] t_k^b, ARMADILLO
                    T_2_new.at(a,i).at(b,j) -= accu(W_2.at(a,i)(j, span(0, n_Electrons-1)) % t1.row(b));
                    T_2_new.at(a,i).at(b,j) += accu(W_2.at(b,i)(j, span(0, n_Electrons-1)) % t1.row(a));

                    // 1/2 I_ab^cd tau_ij^cd, ARMADILLO
                    T_2_new.at(a,i).at(b,j) += 0.5 * accu(integ.at(temp_a, temp_b).submat(span(n_Electrons, Matrix_Size-1), span(n_Electrons, Matrix_Size-1)) % tau.at(i,j));

                    // P(ij) I_ab^cj t_i^c, ARMADILLO
                    T_2_new.at(a,i).at(b,j) -= accu(integ.at(temp_a, temp_b)(span(n_Electrons, Matrix_Size-1), i) % t1.col(j));
                    T_2_new.at(a,i).at(b,j) += accu(integ.at(temp_a, temp_b)(span(n_Electrons, Matrix_Size-1), j) % t1.col(i));

                    // 1/2 [W_1] tau_kl^ab, ARMADILLO
                    T_2_new.at(a,i).at(b,j) += 0.5 * accu(W_1.at(i,j) % tau2.at(a,b));

                    // P(ij) t_jk^ab [F_2], ARMADILLO
                    T_2_new.at(a,i).at(b,j) += accu(T_2.at(a,j)(b, span(0,n_Electrons-1)) % D2.row(i));
                    T_2_new.at(a,i).at(b,j) -= accu(T_2.at(a,i)(b, span(0,n_Electrons-1)) % D2.row(j));

                    // P(ab) t_ij^bc [F_3]_c^a, ARMADILLO
                    T_2_new.at(a,i).at(b,j) -= accu(T_2.at(b,i)(span(), j) % D3.row(a).t());
                    T_2_new.at(a,i).at(b,j) += accu(T_2.at(a,i)(span(), j) % D3.row(b).t());

                    // Divided by Denominator
                    T_2_new.at(a,i).at(b,j) = T_2_new.at(a,i).at(b,j) / denom_abij.at(a)(b,i,j);

                    // Symmetries
                    T_2_new.at(b,i).at(a,j) = -T_2_new.at(a,i).at(b,j);
                    T_2_new.at(a,j).at(b,i) = T_2_new.at(b,i).at(a,j);
                    T_2_new.at(b,j).at(a,i) = T_2_new.at(a,i).at(b,j);

                    j++;
                }
            }
            b++;
        }

        // Symmetric case where a == b, this means a - b = 0 which is an even number
        for (int i = 0; i < n_Electrons; i++)
        {
            for (int j = i+2; j < n_Electrons; j++)
            {
                // I_ab^ij + 1/2 I_ab^cd tau_ij^cd + P(ij) I_ab^cj t_i^c, ARMADILLO
                T_2_new.at(a,i).at(a,j) = integ.at(i,j)(temp_a,temp_a);
                T_2_new.at(a,i).at(a,j) += 0.5 * accu(integ.at(temp_a, temp_a).submat(span(n_Electrons, Matrix_Size-1), span(n_Electrons, Matrix_Size-1)) % tau.at(i,j));
                T_2_new.at(a,i).at(a,j) -= accu(integ.at(temp_a, temp_a)(span(n_Electrons, Matrix_Size-1), i) % t1.col(j));
                T_2_new.at(a,i).at(a,j) += accu(integ.at(temp_a, temp_a)(span(n_Electrons, Matrix_Size-1), j) % t1.col(i));

                // 1/2 [W_1] tau_kl^ab, ARMADILLO
                T_2_new.at(a,i).at(a,j) += 0.5 * accu(W_1.at(i,j) % tau2.at(a,a));

                // P(ij) t_jk^ab [F_2], ARMADILLO
                T_2_new.at(a,i).at(a,j) += accu(T_2.at(a,j)(a, span(0,n_Electrons-1)) % D2.row(i));
                T_2_new.at(a,i).at(a,j) -= accu(T_2.at(a,i)(a, span(0,n_Electrons-1)) % D2.row(j));

                // Denominator
                T_2_new.at(a,i)(a,j) = T_2_new.at(a,i)(a,j) / denom_abij.at(a)(a,i,j);

                // Symmetry
                T_2_new.at(a,j)(a,i) = -T_2_new.at(a,i)(a,j);

                j++;
            }
        }

        // Symmetric case where i == j, which is an EVEN number
        for (int b = a+2; b < unocc_orb; b++)
        {
            temp_b = b + n_Electrons;
            for (int i = 0; i < n_Electrons; i++)
            {
                // I_ab^ij  + 1/2 I_ab^cd tau_ij^cd
                T_2_new.at(a,i).at(b,i) = integ.at(i,i)(temp_a,temp_b);
                T_2_new.at(a,i).at(b,i) += 0.5 * accu(integ.at(temp_a, temp_b).submat(span(n_Electrons, Matrix_Size-1), span(n_Electrons, Matrix_Size-1)) % tau.at(i,i));

                // 1/2 [W_1] tau_kl^ab, ARMADILLO
                T_2_new.at(a,i).at(b,i) += 0.5 * accu(W_1.at(i,i) % tau2.at(a,b));

                // P(ab) t_ij^bc [F_3]_c^a, ARMADILLO
                T_2_new.at(a,i).at(b,i) -= accu(T_2.at(b,i)(span(), i) % D3.row(a).t());
                T_2_new.at(a,i).at(b,i) += accu(T_2.at(a,i)(span(), i) % D3.row(b).t());

                // Denominator
                T_2_new.at(a,i).at(b,i) = T_2_new.at(a,i).at(b,i) / denom_abij.at(a)(b,i,i);

                // Symmetry
                T_2_new.at(b,i).at(a,i) = -T_2_new.at(a,i).at(b,i);
            }
            b++;
        }

        // Symmetric case where a == b AND i == j, here both are EVEN numbers automaticly
        for (int i = 0; i < n_Electrons; i++)
        {
            // I_ab^ij
            T_2_new.at(a,i).at(a,i) = integ.at(i,i)(temp_a,temp_a);
            T_2_new.at(a,i).at(a,i) += 0.5 * accu(integ.at(temp_a, temp_a).submat(span(n_Electrons, Matrix_Size-1), span(n_Electrons, Matrix_Size-1)) % tau.at(i,i));

            // 1/2 [W_1] tau_kl^ab, ARMADILLO
            T_2_new.at(a,i).at(a,i) += 0.5 * accu(W_1.at(i,i) % tau2.at(a,a));

            // Denominator
            T_2_new.at(a,i)(a,i) = T_2_new.at(a,i)(a,i) / denom_abij.at(a)(a,i,i);
        }
    }
}








