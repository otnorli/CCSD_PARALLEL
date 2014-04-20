#include "cc_general.h"

CC_General::CC_General(int n_N, vec zz, mat rr, string B_s, int n_Elec)
{
    n_Nuclei = n_N;
    Z = zz;
    R = rr;
    Basis_Set = B_s;
    n_Electrons = n_Elec;
}

double CC_General::CCSD(double toler)
{
    return 1;
}

double CC_General::CCSDT(double toler)
{
    return 1;
}

double CC_General::CCSDTQ(double toler)
{
    // General starting stuff
    bool continue_ccsd = true;
    double convergance_check;
    int Itermax = 1000;
    iter = 1;
    E_old = 0;
    E_new = 0;

    // Initializion
    double E_HF;
    Hartree_Fock_Solver HartFock(n_Nuclei, Z, R, Basis_Set, n_Electrons, true);
    E_HF = HartFock.get_Energy(toler);
    Matrix_Size = 2*HartFock.ReturnMatrixSize();

    // Transform to MO basis
    c = HartFock.ReturnC();
    construct_fs(HartFock.return_eigval_F());
    construct_integrals(HartFock.Return_Indexed_Q());
    HartFock.Delete_Everything();

    // Define amplitudes
    t1 = zeros(Matrix_Size, Matrix_Size);
    t1_new = zeros(Matrix_Size, Matrix_Size);
    t2.set_size(Matrix_Size, Matrix_Size);
    t2_new.set_size(Matrix_Size, Matrix_Size);
    for (int i = 0; i < Matrix_Size; i++)
    {
        for (int j = 0; j < Matrix_Size; j++)
        {
            t2(i,j) = t1;
            t2_new(i,j) = t1;
        }
    }

    t3 = new double*****[Matrix_Size];
    t3_new = new double*****[Matrix_Size];

    t4 = new double*******[Matrix_Size];
    t4_new = new double*******[Matrix_Size];

    for (int i = n_Electrons; i < Matrix_Size; i++) // Leave un-needed parts empty
    {
        t3[i] = new double ****[Matrix_Size];
        t3_new[i] = new double ****[Matrix_Size];
        t4[i] = new double ******[Matrix_Size];
        t4_new[i] = new double ******[Matrix_Size];

        for (int j = n_Electrons; j < Matrix_Size; j++)
        {
            t3[i][j] = new double ***[Matrix_Size];
            t3_new[i][j] = new double ***[Matrix_Size];
            t4[i][j] = new double *****[Matrix_Size];
            t4_new[i][j] = new double *****[Matrix_Size];

            for (int k = n_Electrons; k < Matrix_Size; k++)
            {
                t3[i][j][k] = new double **[Matrix_Size];
                t3_new[i][j][k] = new double **[Matrix_Size];
                t4[i][j][k] = new double ****[Matrix_Size];
                t4_new[i][j][k] = new double ****[Matrix_Size];

                for (int l = 0; l < n_Electrons; l++) // Leave un-needed parts empty
                {
                    t3[i][j][k][l] = new double *[n_Electrons];
                    t3_new[i][j][k][l] = new double *[n_Electrons];

                    for (int m = 0; m < n_Electrons; m++)
                    {
                        t3[i][j][k][l][m] = new double [n_Electrons];
                        t3_new[i][j][k][l][m] = new double [n_Electrons];
                    }
                }


                for (int l = n_Electrons; l < Matrix_Size; l++) // Leave un-needed parts empty
                {
                    t4[i][j][k][l] = new double ***[n_Electrons];
                    t4_new[i][j][k][l] = new double ***[n_Electrons];

                    for (int m = 0; m < n_Electrons; m++)
                    {

                        t4[i][j][k][l][m] = new double **[n_Electrons];
                        t4_new[i][j][k][l][m] = new double **[n_Electrons];

                        for (int a = 0; a < n_Electrons; a++)
                        {
                            t4[i][j][k][l][m][a] = new double *[n_Electrons];
                            t4_new[i][j][k][l][m][a] = new double *[n_Electrons];
                            for (int b = 0; b < n_Electrons; b++)
                            {
                                t4[i][j][k][l][m][a][b] = new double [n_Electrons];
                                t4_new[i][j][k][l][m][a][b] = new double [n_Electrons];
                            }
                        }
                    }
                }
            }
        }
    }

    // Start calculations
    E_new = Calc_Energy(); // Starting energy
    cout << "Energi: " << E_new << " Steg: " << iter << endl;

    while (continue_ccsd == true)
    {
        // Update old energy
        E_old = E_new;

        // Calc new amplitudes
        Fill_t1_new();
        Fill_t2_new();
        Fill_t3_new();
        Fill_t4_new();

        // Update amplitudes
        t1 = t1_new;
        t2 = t2_new;
        t3 = t3_new;
        t4 = t4_new;

        // Calc energy
        E_new = Calc_Energy();
        iter += 1;

        // Output energy
        cout << "Energi: " << E_new << " Steg: " << iter << endl;

        // Check convergance
        convergance_check = sqrt((E_new - E_old) * (E_new - E_old));
        if (convergance_check < toler || iter >= Itermax)
        {
            continue_ccsd = false;
        }
    }

    // Iterations finished, cout result and return
    cout << "CCSDTQ correction [au] = " << E_new << endl;
    return E_new+E_HF;
}

double CC_General::Calc_Energy()
{
    // Calculate the energy here, currently calculations are "simplified" to always return 1
    return 1;
}

void CC_General::Fill_t1_new()
{
    // T1 amplitudes
    for (int a = n_Electrons; a < Matrix_Size; a++)
    {
        for (int i = 0; i < n_Electrons; i++)
        {
            // + fs_ia
            t1_new(a,i) = fs(i,a);

            // + f_ab t_bi
            for (int b = n_Electrons; b < Matrix_Size; b++)
            {
                t1_new(a,i) += EqualFunc(a,b) * fs(a,b);
            }

            // - f_ji t_aj
            for (int j = 0; j < n_Electrons; j++)
            {
                t1_new(a,i) -= EqualFunc(i,j) * fs(i,j);
            }

            // + <ja||ib> t_bj
            for (int b = n_Electrons; b < Matrix_Size; b++)
            {
                for (int j = 0; j < n_Electrons; j++)
                {
                    t1_new(a,i) -= MOs(j,a)(i,b) * t1(b,j);
                }
            }

            // + f_jb t_abij
            for (int b = n_Electrons; b < Matrix_Size; b++)
            {
                for (int j = 0; j < n_Electrons; j++)
                {
                    t1_new(a,i) += fs(j,b) * t2(a,b)(i,j);
                }
            }

            // + 1/2 <aj||bc t_bcij
            for (int b = n_Electrons; b < Matrix_Size; b++)
            {
                for (int c = n_Electrons; c < Matrix_Size; c++)
                {
                    for (int j = 0; j < n_Electrons; j++)
                    {
                        t1_new(a,i) += 0.5 * MOs(a,j)(b,c) * t2(b,c)(i,j);
                    }
                }
            }

            // - 1/2 <jk||ib> t_abjk
            for (int j = 0; j < n_Electrons; j++)
            {
                for (int k = 0; k < n_Electrons; k++)
                {
                    for (int b = n_Electrons; b < Matrix_Size; b++)
                    {
                        t1_new(a,i) -= 0.5 * MOs(j,k)(i,b) * t2(a,b)(j,k);
                    }
                }
            }

            // + 1/4 <jk||bc> t_abcijk
            for (int b = n_Electrons; b < Matrix_Size; b++)
            {
                for (int c = n_Electrons; c < Matrix_Size; c++)
                {
                    for (int j = 0; j < n_Electrons; j++)
                    {
                        for (int k = 0; k < n_Electrons; k++)
                        {
                            t1_new(a,i) += 0.25 * MOs(j,k)(b,c) * t3[a][b][c][i][j][k];
                        }
                    }
                }
            }

            // - f_jb t_ib t_ja
            for (int j = 0; j < n_Electrons; j++)
            {
                for (int b = n_Electrons; b < Matrix_Size; b++)
                {
                    t1_new(a,i) -= fs(j,b) * t1(b,i) * t1(a,j);
                }
            }

            // + <aj||bc> t_bi t_cj
            for (int b = n_Electrons; b < Matrix_Size; b++)
            {
                for (int j = 0; j < n_Electrons; j++)
                {
                    for (int c = n_Electrons; c < Matrix_Size; c++)
                    {
                        t1_new(a,i) += MOs(a,j)(b,c) * t1(b,i) * t1(c,j);
                    }
                }
            }

            // - <jk||ib> t_aj t_bk
            for (int b = n_Electrons; b < Matrix_Size; b++)
            {
                for (int j = 0; j < n_Electrons; j++)
                {
                    for (int k = 0; k < n_Electrons; k++)
                    {
                        t1_new(a,i) -= MOs(j,k)(i,b) * t1(a,j) * t1(b,k);
                    }
                }
            }

            // + <jk||bc> t_ck t_abij
            for (int b = n_Electrons; b < Matrix_Size; b++)
            {
                for (int j = 0; j < n_Electrons; j++)
                {
                    for (int k = 0; k < n_Electrons; k++)
                    {
                        for (int c = n_Electrons; c < Matrix_Size; c++)
                        {
                            t1_new(a,i) += MOs(j,k)(b,c) * t1(c,k) * t2(a,b)(i,j);
                        }
                    }
                }
            }

            // - 1/2 <jk||bc> t_bi t_acjk
            for (int b = n_Electrons; b < Matrix_Size; b++)
            {
                for (int j = 0; j < n_Electrons; j++)
                {
                    for (int k = 0; k < n_Electrons; k++)
                    {
                        for (int c = n_Electrons; c < Matrix_Size; c++)
                        {
                            t1_new(a,i) -= 0.5*MOs(j,k)(b,c) * t1(b,i) * t2(a,c)(j,k);
                        }
                    }
                }
            }

            // - 1/2 <jk||bc> t_ak t_bcji
            for (int b = n_Electrons; b < Matrix_Size; b++)
            {
                for (int j = 0; j < n_Electrons; j++)
                {
                    for (int k = 0; k < n_Electrons; k++)
                    {
                        for (int c = n_Electrons; c < Matrix_Size; c++)
                        {
                            t1_new(a,i) -= 0.5*MOs(j,k)(b,c) * t1(a,k) * t2(b,c)(j,i);
                        }
                    }
                }
            }

            // - <jk||bc> t_bj t_ci t_ak
            for (int b = n_Electrons; b < Matrix_Size; b++)
            {
                for (int j = 0; j < n_Electrons; j++)
                {
                    for (int k = 0; k < n_Electrons; k++)
                    {
                        for (int c = n_Electrons; c < Matrix_Size; c++)
                        {
                            t1_new(a,i) -= MOs(j,k)(b,c) * t1(b,j) * t1(c,i) * t1(a,k);
                        }
                    }
                }
            }
            t1_new(a,i) = t1_new(a,i) / (fs(i,i) - fs(a,a));
        }
    }
}

void CC_General::Fill_t2_new()
{
    // T2 amplitudes
    for (int a = n_Electrons; a < Matrix_Size; a++)
    {
        for (int b = a+1; b < Matrix_Size; b++)
        {
            for (int i = 0; i < n_Electrons; i++)
            {
                for (int j = i+1; j< n_Electrons; j++)
                {
                    // + <ab||ij>
                    t2_new(a,b)(i,j) = MOs(a,b)(i,j);

                    // + P(ij) <cj||ab> t_ci
                    for (int c = n_Electrons; c < Matrix_Size; c++)
                    {
                        t2_new(a,b)(i,j) += MOs(c,j)(a,b) * t1(c,i);
                        t2_new(a,b)(i,j) -= MOs(c,i)(a,b) * t1(c,j);
                    }

                    // - P(ab) <kb||ij> t_ak
                    for (int k = 0; k < n_Electrons; k++)
                    {
                        t2_new(a,b)(i,j) -= MOs(k,b)(i,j) * t1(a,k);
                        t2_new(a,b)(i,j) += MOs(k,a)(i,j) * t1(b,k);
                    }

                    // + P(ab) f_bc t_acij
                    for (int c = n_Electrons; c < Matrix_Size; c++)
                    {
                        t2_new(a,b)(i,j) += EqualFunc(b,c) * fs(b,c) * t2(a,c)(i,j);
                        t2_new(a,b)(i,j) -= EqualFunc(a,c) * fs(a,c) * t2(b,c)(i,j);
                    }

                    // - P(ij) f_kj t_abik
                    for (int k = 0; k < n_Electrons; k++)
                    {
                        t2_new(a,b)(i,j) -= EqualFunc(k,j) * fs(k,j) * t2(a,b)(i,k);
                        t2_new(a,b)(i,j) += EqualFunc(k,i) * fs(k,i) * t2(a,b)(j,k);
                    }

                    // - P(ij) P(ab) <kb||jc> t_acik
                    for (int k = 0; k < n_Electrons; k++)
                    {
                        for (int c = n_Electrons; c < Matrix_Size; c++)
                        {
                            t2_new(a,b)(i,j) -= MOs(k,b)(j,c) * t2(a,c)(i,k);
                            t2_new(a,b)(i,j) += MOs(k,a)(j,c) * t2(b,c)(i,k);
                            t2_new(a,b)(i,j) += MOs(k,b)(i,c) * t2(a,c)(j,k);
                            t2_new(a,b)(i,j) -= MOs(k,a)(i,c) * t2(b,c)(j,k);
                        }
                    }

                    // + 1/2 <kl||ij> t_abkl
                    for (int k = 0; k < n_Electrons; k++)
                    {
                        for (int l = 0; l < n_Electrons; l++)
                        {
                            t2_new(a,b)(i,j) += MOs(k,l)(i,j) * t2(a,b)(k,l);
                        }
                    }

                    // + 1/2 <ab||cd> t_cdij
                    for (int c = n_Electrons; c < Matrix_Size; c++)
                    {
                        for (int d = n_Electrons; d < Matrix_Size; d++)
                        {
                            t2_new(a,b)(i,j) += 0.5 * MOs(a,b)(c,d) * t2(c,d)(i,j);
                        }
                    }

                    // + f_kc t_abcijk
                    for (int c = n_Electrons; c < Matrix_Size; c++)
                    {
                        for (int k = 0; k < n_Electrons; k++)
                        {
                            t2_new(a,b)(i,j) += fs(k,c) * t3[a][b][c][i][j][k];
                        }
                    }

                    // + 1/2 P(ab) <bk||cd> t_acdijk
                    for (int c = n_Electrons; c < Matrix_Size; c++)
                    {
                        for (int d = n_Electrons; d < Matrix_Size; d++)
                        {
                            for (int k = 0; k < n_Electrons; k++)
                            {
                                t2_new(a,b)(i,j) += 0.5 * MOs(b,k)(c,d) * t3[a][c][d][i][j][k];
                                t2_new(a,b)(i,j) -= 0.5 * MOs(a,k)(c,d) * t3[b][c][d][i][j][k];
                            }
                        }
                    }

                    // - 1/2 P(ij) <kl||jc> t_abcikl
                    for (int c = n_Electrons; c < Matrix_Size; c++)
                    {
                        for (int l = 0; l < n_Electrons; l++)
                        {
                            for (int k = 0; k < n_Electrons; k++)
                            {
                                t2_new(a,b)(i,j) -= 0.5 * MOs(k,l)(j,c) * t3[a][b][c][i][k][l];
                                t2_new(a,b)(i,j) += 0.5 * MOs(k,l)(i,c) * t3[a][b][c][j][k][l];
                            }
                        }
                    }

                    // + 1/4 <kl||cd> t_abcdijkl
                    for (int c = n_Electrons; c < Matrix_Size; c++)
                    {
                        for (int l = 0; l < n_Electrons; l++)
                        {
                            for (int k = 0; k < n_Electrons; k++)
                            {
                                for (int d = n_Electrons; d < Matrix_Size; d++)
                                {
                                    t2_new(a,b)(i,j) += 0.25 * MOs(k,l)(c,d) * t4[a][b][c][d][i][j][k][l];
                                }
                            }
                        }
                    }

                    // - P(ab) P(ij) <ak||cj> t_ci t_bk
                    for (int c = n_Electrons; c < Matrix_Size; c++)
                    {
                        for (int k = 0; k < n_Electrons; k++)
                        {
                            t2_new(a,b)(i,j) -= MOs(a,k)(c,j) * t1(c,i) * t1(b,k);
                            t2_new(a,b)(i,j) += MOs(b,k)(c,j) * t1(c,i) * t1(a,k);
                            t2_new(a,b)(i,j) += MOs(b,k)(c,i) * t1(c,j) * t1(a,k);
                            t2_new(a,b)(i,j) -= MOs(a,k)(c,i) * t1(c,j) * t1(b,k);
                        }
                    }

                    // + <kl||ij> t_ak t_bl
                    for (int k = 0; k < n_Electrons; k++)
                    {
                        for (int l = 0; l < n_Electrons; l++)
                        {
                            t2_new(a,b)(i,j) += MOs(k,l)(i,j) * t1(a,k) * t1(b,l);
                        }
                    }

                    // + <ab||cd> t_ci t_dj
                    for (int c = n_Electrons; c < Matrix_Size; c++)
                    {
                        for (int d = n_Electrons; d < Matrix_Size; d++)
                        {
                            t2_new(a,b)(i,j) += MOs(a,b)(c,d) * t1(c,i) * t1(d,j);
                        }
                    }

                    // - P(ij) f_kc t_abkj t_ci
                    for (int c = n_Electrons; c < Matrix_Size; c++)
                    {
                        for (int k = 0; k < n_Electrons; k++)
                        {
                            t2_new(a,b)(i,j) -= fs(k,c) * t2(a,b)(k,j) * t1(c,i);
                            t2_new(a,b)(i,j) += fs(k,c) * t2(a,b)(k,i) * t1(c,j);
                        }
                    }

                    // - P(ab) f_kc t_acij t_bk
                    for (int c = n_Electrons; c < Matrix_Size; c++)
                    {
                        for (int k = 0; k < n_Electrons; k++)
                        {
                            t2_new(a,b)(i,j) -= fs(k,c) * t2(a,c)(i,j) * t1(b,k);
                            t2_new(a,b)(i,j) += fs(k,c) * t2(b,c)(i,j) * t1(a,k);
                        }
                    }

                    // + P(ab) P(ij) <kb||cd> t_acik t_dj
                    for (int c = n_Electrons; c < Matrix_Size; c++)
                    {
                        for (int k = 0; k < n_Electrons; k++)
                        {
                            for (int d = n_Electrons; d < Matrix_Size; d++)
                            {
                                t2_new(a,b)(i,j) += MOs(k,b)(c,d) * t2(a,c)(i,k) * t1(d,j);
                                t2_new(a,b)(i,j) -= MOs(k,a)(c,d) * t2(b,c)(i,k) * t1(d,j);
                                t2_new(a,b)(i,j) -= MOs(k,b)(c,d) * t2(a,c)(j,k) * t1(d,i);
                                t2_new(a,b)(i,j) += MOs(k,a)(c,d) * t2(b,c)(j,k) * t1(d,i);
                            }
                        }
                    }

                    // - P(ab) P(ij) <kl||cj t_acik t_bl
                    for (int c = n_Electrons; c < Matrix_Size; c++)
                    {
                        for (int k = 0; k < n_Electrons; k++)
                        {
                            for (int l = 0; l < n_Electrons; l++)
                            {
                                t2_new(a,b)(i,j) -= MOs(k,l)(c,j) * t2(a,c)(i,k) * t1(b,l);
                                t2_new(a,b)(i,j) += MOs(k,l)(c,i) * t2(a,c)(j,k) * t1(b,l);
                                t2_new(a,b)(i,j) += MOs(k,l)(c,j) * t2(b,c)(i,k) * t1(a,l);
                                t2_new(a,b)(i,j) -= MOs(k,l)(c,i) * t2(b,c)(j,k) * t1(a,l);
                            }
                        }
                    }

                    // - P(ij) <kl||jc> t_abik t_cl
                    for (int c = n_Electrons; c < Matrix_Size; c++)
                    {
                        for (int k = 0; k < n_Electrons; k++)
                        {
                            for (int l = 0; l < n_Electrons; l++)
                            {
                                t2_new(a,b)(i,j) -= MOs(k,l)(j,c) * t2(a,b)(i,k) * t1(c,l);
                                t2_new(a,b)(i,j) += MOs(k,l)(i,c) * t2(a,b)(j,k) * t1(c,l);
                            }
                        }
                    }

                    // + P(ab) <bk||cd> t_acij t_dk
                    for (int c = n_Electrons; c < Matrix_Size; c++)
                    {
                        for (int d = n_Electrons; d < Matrix_Size; d++)
                        {
                            for (int k = 0; k < n_Electrons; k++)
                            {
                                t2_new(a,b)(i,j) += MOs(b,k)(c,d) * t2(a,c)(i,j) * t1(d,k);
                                t2_new(a,b)(i,j) -= MOs(a,b)(c,d) * t2(b,c)(i,j) * t1(d,k);
                            }
                        }
                    }

                    // + 1/2 P(ij) <kl||cj> t_abkl t_ci
                    for (int c = n_Electrons; c < Matrix_Size; c++)
                    {
                        for (int k = 0; k < n_Electrons; k++)
                        {
                            for (int l = 0; l < n_Electrons; l++)
                            {
                                t2_new(a,b)(i,j) += MOs(k,l)(c,j) * t2(a,b)(k,l) * t1(c,i);
                                t2_new(a,b)(i,j) -= MOs(k,l)(c,i) * t2(a,b)(k,l) * t1(c,j);
                            }
                        }
                    }

                    // - 1/2 P(ab) <ak||cd> t_cdij t_bk
                    for (int c = n_Electrons; c < Matrix_Size; c++)
                    {
                        for (int d = n_Electrons; d < Matrix_Size; d++)
                        {
                            for (int k = 0; k < n_Electrons; k++)
                            {
                                t2_new(a,b)(i,j) -= MOs(a,k)(c,d) * t2(c,d)(i,j) * t1(b,k);
                                t2_new(a,b)(i,j) += MOs(b,k)(c,d) * t2(c,d)(i,j) * t1(a,k);
                            }
                        }
                    }

                    // + <kl||cd> t_abcijk t_dl
                    for (int c = n_Electrons; c < Matrix_Size; c++)
                    {
                        for (int d = n_Electrons; d < Matrix_Size; d++)
                        {
                            for (int k = 0; k < n_Electrons; k++)
                            {
                                for (int l = 0; l < n_Electrons; l++)
                                {
                                    t2_new(a,b)(i,j) += MOs(k,l)(c,d) * t3[a][b][c][i][j][k] * t1(d,l);
                                }
                            }
                        }
                    }

                    // - 1/2 P(ab) <kl||cd> t_acdikj t_bl
                    for (int c = n_Electrons; c < Matrix_Size; c++)
                    {
                        for (int d = n_Electrons; d < Matrix_Size; d++)
                        {
                            for (int k = 0; k < n_Electrons; k++)
                            {
                                for (int l = 0; l < n_Electrons; l++)
                                {
                                    t2_new(a,b)(i,j) -= 0.5 * MOs(k,l)(c,d) * t3[a][c][d][i][k][j] * t1(b,l);
                                    t2_new(a,b)(i,j) += 0.5 * MOs(k,l)(c,d) * t3[b][c][d][i][k][j] * t1(a,l);
                                }
                            }
                        }
                    }

                    // - 1/2 P(ij) <kl||cd> t_adbklj t_ci
                    for (int c = n_Electrons; c < Matrix_Size; c++)
                    {
                        for (int d = n_Electrons; d < Matrix_Size; d++)
                        {
                            for (int k = 0; k < n_Electrons; k++)
                            {
                                for (int l = 0; l < n_Electrons; l++)
                                {
                                    t2_new(a,b)(i,j) -= 0.5 * MOs(k,l)(c,d) * t3[a][d][b][k][l][j] * t1(c,i);
                                    t2_new(a,b)(i,j) += 0.5 * MOs(k,l)(c,d) * t3[a][d][b][k][l][i] * t1(c,j);
                                }
                            }
                        }
                    }

                    // - 1/2 P(ij) <kl||cd> t_cdki t_ablj
                    for (int c = n_Electrons; c < Matrix_Size; c++)
                    {
                        for (int d = n_Electrons; d < Matrix_Size; d++)
                        {
                            for (int k = 0; k < n_Electrons; k++)
                            {
                                for (int l = 0; l < n_Electrons; l++)
                                {
                                    t2_new(a,b)(i,j) -= 0.5 * MOs(k,l)(c,d) * t2(c,d)(k,i) * t2(a,b)(l,j);
                                    t2_new(a,b)(i,j) += 0.5 * MOs(k,l)(c,d) * t2(c,d)(k,j) * t2(a,b)(l,i);
                                }
                            }
                        }
                    }

                    // - 1/2 P(ab) <kl||cd> t_acij t_bdkl
                    for (int c = n_Electrons; c < Matrix_Size; c++)
                    {
                        for (int d = n_Electrons; d < Matrix_Size; d++)
                        {
                            for (int k = 0; k < n_Electrons; k++)
                            {
                                for (int l = 0; l < n_Electrons; l++)
                                {
                                    t2_new(a,b)(i,j) -= 0.5 * MOs(k,l)(c,d) * t2(a,c)(i,j) * t2(b,d)(k,l);
                                    t2_new(a,b)(i,j) += 0.5 * MOs(k,l)(c,d) * t2(b,c)(i,j) * t2(a,d)(k,l);
                                }
                            }
                        }
                    }

                    // + P(ij) <kl||cd> t_acik t_dblj
                    for (int c = n_Electrons; c < Matrix_Size; c++)
                    {
                        for (int d = n_Electrons; d < Matrix_Size; d++)
                        {
                            for (int k = 0; k < n_Electrons; k++)
                            {
                                for (int l = 0; l < n_Electrons; l++)
                                {
                                    t2_new(a,b)(i,j) += MOs(k,l)(c,d) * t2(a,c)(i,k) * t2(d,b)(l,j);
                                    t2_new(a,b)(i,j) -= MOs(k,l)(c,d) * t2(a,c)(j,k) * t2(d,b)(l,i);
                                }
                            }
                        }
                    }

                    // + 1/4 <kl||cd> t_cdij t_abkl
                    for (int c = n_Electrons; c < Matrix_Size; c++)
                    {
                        for (int d = n_Electrons; d < Matrix_Size; d++)
                        {
                            for (int k = 0; k < n_Electrons; k++)
                            {
                                for (int l = 0; l < n_Electrons; l++)
                                {
                                    t2_new(a,b)(i,j) += 0.5 * MOs(k,l)(c,d) * t2(c,d)(i,j) * t2(a,b)(k,l);
                                }
                            }
                        }
                    }

                    // + P(ij) <kl||cj> t_ci t_ak t_bl
                    for (int c = n_Electrons; c < Matrix_Size; c++)
                    {
                        for (int k = 0; k < n_Electrons; k++)
                        {
                            for (int l = 0; l < n_Electrons; l++)
                            {
                                t2_new(a,b)(i,j) += MOs(k,l)(c,j) * t1(c,i) * t1(a,k) * t1(b,l);
                                t2_new(a,b)(i,j) -= MOs(k,l)(c,i) * t1(c,j) * t1(a,k) * t1(b,l);
                            }
                        }
                    }

                    // - P(ab) <ak||cd> t_ci t_dj t_bk
                    for (int c = n_Electrons; c < Matrix_Size; c++)
                    {
                        for (int k = 0; k < n_Electrons; k++)
                        {
                            for (int d = n_Electrons; d < Matrix_Size; d++)
                            {
                                t2_new(a,b)(i,j) -= MOs(a,k)(c,d) * t1(c,i) * t1(d,j) * t1(b,k);
                                t2_new(a,b)(i,j) += MOs(b,k)(c,d) * t1(c,i) * t1(d,j) * t1(a,k);
                            }
                        }
                    }

                    // - P(ab) P(ij) <kl||cd> t_dblj t_ci t_ak
                    for (int c = n_Electrons; c < Matrix_Size; c++)
                    {
                        for (int d = n_Electrons; d < Matrix_Size; d++)
                        {
                            for (int k = 0; k < n_Electrons; k++)
                            {
                                for (int l = 0; l < n_Electrons; l++)
                                {
                                    t2_new(a,b)(i,j) -= MOs(k,l)(c,d) * t2(d,b)(l,j) * t1(c,i) * t1(a,k);
                                    t2_new(a,b)(i,j) += MOs(k,l)(c,d) * t2(d,a)(l,j) * t1(c,i) * t1(b,k);
                                    t2_new(a,b)(i,j) += MOs(k,l)(c,d) * t2(d,b)(l,i) * t1(c,j) * t1(a,k);
                                    t2_new(a,b)(i,j) -= MOs(k,l)(c,d) * t2(d,a)(l,i) * t1(c,j) * t1(b,k);
                                }
                            }
                        }
                    }

                    // - P(ij) <kl||cd> t_ablj t_ck t_di
                    for (int c = n_Electrons; c < Matrix_Size; c++)
                    {
                        for (int d = n_Electrons; d < Matrix_Size; d++)
                        {
                            for (int k = 0; k < n_Electrons; k++)
                            {
                                for (int l = 0; l < n_Electrons; l++)
                                {
                                    t2_new(a,b)(i,j) -= MOs(k,l)(c,d) * t2(a,b)(l,j) * t1(c,k) * t1(d,i);
                                    t2_new(a,b)(i,j) += MOs(k,l)(c,d) * t2(a,b)(l,i) * t1(c,k) * t1(d,j);
                                }
                            }
                        }
                    }

                    // - P(ab) <kl||cd> t_acij t_bk t_dl
                    for (int c = n_Electrons; c < Matrix_Size; c++)
                    {
                        for (int d = n_Electrons; d < Matrix_Size; d++)
                        {
                            for (int k = 0; k < n_Electrons; k++)
                            {
                                for (int l = 0; l < n_Electrons; l++)
                                {
                                    t2_new(a,b)(i,j) -= MOs(k,l)(c,d) * t2(a,c)(i,j) * t1(b,k) * t1(d,l);
                                    t2_new(a,b)(i,j) += MOs(k,l)(c,d) * t2(b,c)(i,j) * t1(a,k) * t1(d,l);
                                }
                            }
                        }
                    }

                    // + <kl||cd> t_abkl t_ci t_dj
                    for (int c = n_Electrons; c < Matrix_Size; c++)
                    {
                        for (int d = n_Electrons; d < Matrix_Size; d++)
                        {
                            for (int k = 0; k < n_Electrons; k++)
                            {
                                for (int l = 0; l < n_Electrons; l++)
                                {
                                    t2_new(a,b)(i,j) += MOs(k,l)(c,d) * t2(a,b)(k,l) * t1(c,i) * t1(d,j);
                                }
                            }
                        }
                    }

                    // + <kl||cd> t_cdij t_ak t_bl
                    for (int c = n_Electrons; c < Matrix_Size; c++)
                    {
                        for (int d = n_Electrons; d < Matrix_Size; d++)
                        {
                            for (int k = 0; k < n_Electrons; k++)
                            {
                                for (int l = 0; l < n_Electrons; l++)
                                {
                                    t2_new(a,b)(i,j) += MOs(k,l)(c,d) * t2(c,d)(i,j) * t1(a,k) * t1(b,l);
                                }
                            }
                        }
                    }

                    // + <kl||cd> t_ci t_ak t_dj t_bl
                    for (int c = n_Electrons; c < Matrix_Size; c++)
                    {
                        for (int d = n_Electrons; d < Matrix_Size; d++)
                        {
                            for (int k = 0; k < n_Electrons; k++)
                            {
                                for (int l = 0; l < n_Electrons; l++)
                                {
                                    t2_new(a,b)(i,j) += MOs(k,l)(c,d) * t1(c,i) * t1(a,k) * t1(d,j) * t1(b,l);
                                }
                            }
                        }
                    }

                    // Denominator
                    t2_new(a,b)(i,j) = t2_new(a,b)(i,j) / (fs(i,i) + fs(j,j) - fs(a,a) - fs(b,b));

                    // Symmetries
                    t2_new(b,a)(i,j) = -t2_new(a,b)(i,j);
                    t2_new(a,b)(j,i) = -t2_new(a,b)(i,j);
                    t2_new(b,a)(j,i) =  t2_new(a,b)(i,j);

                }
            }
        }
    }
}

void CC_General::Fill_t3_new()
{
    // T3 amplitudes
    for (int a = n_Electrons; a < Matrix_Size; a++)
    {
        for (int b = a+1; b < Matrix_Size; b++)
        {
            for (int c = b+1; c < Matrix_Size; c++)
            {
                for (int i = 0; i < n_Electrons; i++)
                {
                    for (int j = i+1; j < n_Electrons; j++)
                    {
                        for (int k = j+1; k < n_Electrons; k++)
                        {
                            // Reset
                            t3_new[a][b][c][i][j][k] = 0;

                            // 1: + P(k/ij) P(a/bc) <bc||dk> t_adij
                            for (int d = n_Electrons; d < Matrix_Size; d++)
                            {
                               t3_new[a][b][c][i][j][k] += MOs(b,c)(d,k) * t2(a,d)(i,j);
                               t3_new[a][b][c][i][j][k] -= MOs(b,c)(d,i) * t2(a,d)(k,j);
                               t3_new[a][b][c][i][j][k] -= MOs(b,c)(d,j) * t2(a,d)(i,k);
                            }

                            // 2: - P(i/jk) P(ab/c) <lc||jk> t_abil
                            for (int l = 0; l < n_Electrons; l++)
                            {
                                t3_new[a][b][c][i][j][k] -= MOs(l,c)(j,k) * t2(a,b)(i,l);
                                t3_new[a][b][c][i][j][k] += MOs(l,c)(i,k) * t2(a,b)(j,l);
                                t3_new[a][b][c][i][j][k] += MOs(l,c)(j,i) * t2(a,b)(k,l);

                                t3_new[a][b][c][i][j][k] += MOs(l,a)(j,k) * t2(c,b)(i,l);
                                t3_new[a][b][c][i][j][k] += MOs(l,b)(j,k) * t2(a,c)(i,l);

                                t3_new[a][b][c][i][j][k] -= MOs(l,a)(i,k) * t2(c,b)(j,l);
                                t3_new[a][b][c][i][j][k] -= MOs(l,b)(i,k) * t2(a,c)(j,l);

                                t3_new[a][b][c][i][j][k] -= MOs(l,a)(j,i) * t2(c,b)(k,l);
                                t3_new[a][b][c][i][j][k] -= MOs(l,b)(j,i) * t2(a,c)(k,l);
                            }

                            // 3: - P(ij/k) f_lk t_abcijl
                            for (int l = 0; l < n_Electrons; l++)
                            {
                                t3_new[a][b][c][i][j][k] -= EqualFunc(l,k) * fs(l,k) * t3[a][b][c][i][j][l];
                                t3_new[a][b][c][i][j][k] += EqualFunc(l,i) * fs(l,i) * t3[a][b][c][k][j][l];
                                t3_new[a][b][c][i][j][k] += EqualFunc(j,l) * fs(l,j) * t3[a][b][c][i][k][l];
                            }

                            // 4: + P(ab/c) f_cd t_abdijk
                            for (int d = n_Electrons; d < Matrix_Size; d++)
                            {
                                t3_new[a][b][c][i][j][k] += EqualFunc(c,d) * fs(c,d) * t3[a][b][d][i][j][k];
                                t3_new[a][b][c][i][j][k] -= EqualFunc(a,d) * fs(a,d) * t3[c][b][d][i][j][k];
                                t3_new[a][b][c][i][j][k] -= EqualFunc(b,d) * fs(b,d) * t3[a][c][d][i][j][k];
                            }

                            // 5: - P(ab/c) P(ij/k) <lc||kd> t_abdijl
                            for (int d = n_Electrons; d < Matrix_Size; d++)
                            {
                                for (int l = 0; l < n_Electrons; l++)
                                {
                                    t3_new[a][b][c][i][j][k] -= MOs(l,c)(k,d) * t3[a][b][d][i][j][l];
                                    t3_new[a][b][c][i][j][k] += MOs(l,a)(k,d) * t3[c][b][d][i][j][l];
                                    t3_new[a][b][c][i][j][k] += MOs(l,b)(k,d) * t3[a][c][d][i][j][l];

                                    t3_new[a][b][c][i][j][k] += MOs(l,c)(i,d) * t3[a][b][d][k][j][l];
                                    t3_new[a][b][c][i][j][k] += MOs(l,c)(j,d) * t3[a][b][d][i][k][l];

                                    t3_new[a][b][c][i][j][k] -= MOs(l,a)(i,d) * t3[c][b][d][k][j][l];
                                    t3_new[a][b][c][i][j][k] -= MOs(l,a)(j,d) * t3[c][b][d][i][k][l];

                                    t3_new[a][b][c][i][j][k] -= MOs(l,b)(i,d) * t3[a][c][d][k][j][l];
                                    t3_new[a][b][c][i][j][k] -= MOs(l,b)(j,d) * t3[a][c][d][i][k][l];
                                }
                            }

                            // 6: + 1/2 P(i/jk) <lm||jk> t_abcilm
                            for (int l = 0; l < n_Electrons; l++)
                            {
                                for (int m = 0; m < n_Electrons; m++)
                                {
                                    t3_new[a][b][c][i][j][k] += 0.5 * MOs(l,m)(j,k) * t3[a][b][c][i][l][m];
                                    t3_new[a][b][c][i][j][k] -= 0.5 * MOs(l,m)(i,k) * t3[a][b][c][j][l][m];
                                    t3_new[a][b][c][i][j][k] -= 0.5 * MOs(l,m)(j,i) * t3[a][b][c][k][l][m];
                                }
                            }

                            // 7: + 1/2 P(a/bc) <bc||de> t_adeijk
                            for (int d = n_Electrons; d < Matrix_Size; d++)
                            {
                                for (int e = n_Electrons; e < Matrix_Size; e++)
                                {
                                    t3_new[a][b][c][i][j][k] += 0.5 * MOs(b,c)(d,e) * t3[a][d][e][i][j][k];
                                    t3_new[a][b][c][i][j][k] -= 0.5 * MOs(a,c)(d,e) * t3[b][d][e][i][j][k];
                                    t3_new[a][b][c][i][j][k] -= 0.5 * MOs(b,a)(d,e) * t3[c][d][e][i][j][k];
                                }
                            }

                            // 8: + f_ld t_abcdijkl
                            for (int d = n_Electrons; d < Matrix_Size; d++)
                            {
                                for (int l = 0; l < n_Electrons; l++)
                                {
                                    t3_new[a][b][c][i][j][k] += fs(l,d) * t4[a][b][c][d][i][j][k][l];
                                }
                            }

                            // 9: + 1/2 P(ab/c) <cl||de> t_abdeijkl
                            for (int d = n_Electrons; d < Matrix_Size; d++)
                            {
                                for (int l = 0; l < n_Electrons; l++)
                                {
                                    for (int e = n_Electrons; e < Matrix_Size; e++)
                                    {
                                        t3_new[a][b][c][i][j][k] += 0.5 * MOs(c,l)(d,e) * t4[a][b][d][e][i][j][k][l];
                                        t3_new[a][b][c][i][j][k] -= 0.5 * MOs(b,l)(d,e) * t4[a][c][d][e][i][j][k][l];
                                        t3_new[a][b][c][i][j][k] -= 0.5 * MOs(a,l)(d,e) * t4[c][b][d][e][i][j][k][l];
                                    }
                                }
                            }

                            // 10: - 1/2 <lm||kd> t_abcdijlm
                            for (int d = n_Electrons; d < Matrix_Size; d++)
                            {
                                for (int l = 0; l < n_Electrons; l++)
                                {
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        t3_new[a][b][c][i][j][k] -= 0.5 * MOs(l,m)(k,d) * t4[a][b][c][d][i][j][l][m];
                                    }
                                }
                            }

                            // 11: - P(abc) P(ij/k) <bl||dk> t_adij t_cl
                            for (int d = n_Electrons; d < Matrix_Size; d++)
                            {
                                for (int l = 0; l < n_Electrons; l++)
                                {
                                    t3_new[a][b][c][i][j][k] -= MOs(b,l)(d,k) * t2(a,d)(i,j) * t1(c,l);
                                    t3_new[a][b][c][i][j][k] += MOs(b,l)(d,i) * t2(a,d)(k,j) * t1(c,l);
                                    t3_new[a][b][c][i][j][k] += MOs(b,l)(d,j) * t2(a,d)(i,k) * t1(c,l);

                                    t3_new[a][b][c][i][j][k] += MOs(a,l)(d,k) * t2(b,d)(i,j) * t1(c,l);
                                    t3_new[a][b][c][i][j][k] += MOs(b,l)(d,k) * t2(c,d)(i,j) * t1(a,l);
                                    t3_new[a][b][c][i][j][k] += MOs(c,l)(d,k) * t2(a,d)(i,j) * t1(b,l);
                                    t3_new[a][b][c][i][j][k] -= MOs(a,l)(d,k) * t2(c,d)(i,j) * t1(b,l);
                                    t3_new[a][b][c][i][j][k] -= MOs(c,l)(d,k) * t2(b,d)(i,j) * t1(a,l);

                                    t3_new[a][b][c][i][j][k] -= MOs(a,l)(d,i) * t2(b,d)(k,j) * t1(c,l);
                                    t3_new[a][b][c][i][j][k] -= MOs(b,l)(d,i) * t2(c,d)(k,j) * t1(a,l);
                                    t3_new[a][b][c][i][j][k] -= MOs(c,l)(d,i) * t2(a,d)(k,j) * t1(b,l);
                                    t3_new[a][b][c][i][j][k] += MOs(a,l)(d,i) * t2(c,d)(k,j) * t1(b,l);
                                    t3_new[a][b][c][i][j][k] += MOs(c,l)(d,i) * t2(b,d)(k,j) * t1(b,l);

                                    t3_new[a][b][c][i][j][k] -= MOs(a,l)(d,j) * t2(b,d)(i,k) * t1(c,l);
                                    t3_new[a][b][c][i][j][k] -= MOs(b,l)(d,j) * t2(c,d)(i,k) * t1(a,l);
                                    t3_new[a][b][c][i][j][k] -= MOs(c,l)(d,j) * t2(a,d)(i,k) * t1(b,l);
                                    t3_new[a][b][c][i][j][k] += MOs(a,l)(d,j) * t2(c,d)(i,k) * t1(b,l);
                                    t3_new[a][b][c][i][j][k] += MOs(c,l)(d,j) * t2(b,d)(i,k) * t1(a,l);
                                }
                            }

                            // 12: + P(a/bc) P(ij/k) <bc||de> t_adij t_ek
                            for (int d = n_Electrons; d < Matrix_Size; d++)
                            {
                                for (int e = n_Electrons; e < Matrix_Size; e++)
                                {
                                    t3_new[a][b][c][i][j][k] += MOs(b,c)(d,e) * t2(a,d)(i,j) * t1(e,k);
                                    t3_new[a][b][c][i][j][k] -= MOs(a,c)(d,e) * t2(b,d)(i,j) * t1(e,k);
                                    t3_new[a][b][c][i][j][k] -= MOs(b,a)(d,e) * t2(c,d)(i,j) * t1(e,k);

                                    t3_new[a][b][c][i][j][k] -= MOs(b,c)(d,e) * t2(a,d)(k,j) * t1(e,i);
                                    t3_new[a][b][c][i][j][k] -= MOs(b,c)(d,e) * t2(a,d)(i,k) * t1(e,j);

                                    t3_new[a][b][c][i][j][k] += MOs(a,c)(d,e) * t2(b,d)(k,j) * t1(e,i);
                                    t3_new[a][b][c][i][j][k] += MOs(a,c)(d,e) * t2(b,d)(i,k) * t1(e,j);

                                    t3_new[a][b][c][i][j][k] += MOs(b,a)(d,e) * t2(c,d)(k,j) * t1(e,i);
                                    t3_new[a][b][c][i][j][k] += MOs(b,a)(d,e) * t2(c,d)(i,k) * t1(e,j);
                                }
                            }

                            // 13: + P(ab/c) P(ijk) <lm||jk> t_abil t_cm
                            for (int l = 0; l < n_Electrons; l++)
                            {
                                for (int m = 0; m < n_Electrons; m++)
                                {
                                    t3_new[a][b][c][i][j][k] += MOs(l,m)(j,k) * t2(a,b)(i,l) * t1(c,m);
                                    t3_new[a][b][c][i][j][k] -= MOs(l,m)(j,k) * t2(c,b)(i,l) * t1(a,m);
                                    t3_new[a][b][c][i][j][k] -= MOs(l,m)(j,k) * t2(a,c)(i,l) * t1(b,m);

                                    t3_new[a][b][c][i][j][k] -= MOs(l,m)(i,k) * t2(a,b)(j,l) * t1(c,m);
                                    t3_new[a][b][c][i][j][k] -= MOs(l,m)(j,i) * t2(a,b)(k,l) * t1(c,m);
                                    t3_new[a][b][c][i][j][k] -= MOs(l,m)(k,j) * t2(a,b)(i,l) * t1(c,m);
                                    t3_new[a][b][c][i][j][k] += MOs(l,m)(k,i) * t2(a,b)(j,l) * t1(c,m);
                                    t3_new[a][b][c][i][j][k] += MOs(l,m)(i,j) * t2(a,b)(k,l) * t1(c,m);

                                    t3_new[a][b][c][i][j][k] += MOs(l,m)(i,k) * t2(c,b)(j,l) * t1(a,m);
                                    t3_new[a][b][c][i][j][k] += MOs(l,m)(j,i) * t2(c,b)(k,l) * t1(a,m);
                                    t3_new[a][b][c][i][j][k] += MOs(l,m)(k,j) * t2(c,b)(i,l) * t1(a,m);
                                    t3_new[a][b][c][i][j][k] -= MOs(l,m)(k,i) * t2(c,b)(j,l) * t1(a,m);
                                    t3_new[a][b][c][i][j][k] -= MOs(l,m)(i,j) * t2(c,b)(k,l) * t1(a,m);

                                    t3_new[a][b][c][i][j][k] += MOs(l,m)(i,k) * t2(a,c)(j,l) * t1(b,m);
                                    t3_new[a][b][c][i][j][k] += MOs(l,m)(j,i) * t2(a,c)(k,l) * t1(b,m);
                                    t3_new[a][b][c][i][j][k] += MOs(l,m)(k,j) * t2(a,c)(i,l) * t1(b,m);
                                    t3_new[a][b][c][i][j][k] -= MOs(l,m)(k,i) * t2(a,c)(j,l) * t1(b,m);
                                    t3_new[a][b][c][i][j][k] -= MOs(l,m)(i,j) * t2(a,c)(k,l) * t1(b,m);
                                }
                            }

                            // 14: - P(ab/c) P(ijk) <lc||jd> t_abil t_dk
                            for (int d = n_Electrons; d < Matrix_Size; d++)
                            {
                                for (int l = 0; l < n_Electrons; l++)
                                {
                                    t3_new[a][b][c][i][j][k] -= MOs(l,c)(j,d) * t2(a,b)(i,l) * t1(d,k);
                                    t3_new[a][b][c][i][j][k] += MOs(l,a)(j,d) * t2(c,b)(i,l) * t1(d,k);
                                    t3_new[a][b][c][i][j][k] += MOs(l,b)(j,d) * t2(a,c)(i,l) * t1(d,k);

                                    t3_new[a][b][c][i][j][k] += MOs(l,c)(i,d) * t2(a,b)(j,l) * t1(d,k);
                                    t3_new[a][b][c][i][j][k] += MOs(l,c)(k,d) * t2(a,b)(i,l) * t1(d,j);
                                    t3_new[a][b][c][i][j][k] += MOs(l,c)(j,d) * t2(a,b)(k,l) * t1(d,i);
                                    t3_new[a][b][c][i][j][k] -= MOs(l,c)(k,d) * t2(a,b)(j,l) * t1(d,i);
                                    t3_new[a][b][c][i][j][k] -= MOs(l,c)(i,d) * t2(a,b)(k,l) * t1(d,j);

                                    t3_new[a][b][c][i][j][k] -= MOs(l,a)(i,d) * t2(c,b)(j,l) * t1(d,k);
                                    t3_new[a][b][c][i][j][k] -= MOs(l,a)(k,d) * t2(c,b)(i,l) * t1(d,j);
                                    t3_new[a][b][c][i][j][k] -= MOs(l,a)(j,d) * t2(c,b)(k,l) * t1(d,i);
                                    t3_new[a][b][c][i][j][k] += MOs(l,a)(k,d) * t2(c,b)(j,l) * t1(d,i);
                                    t3_new[a][b][c][i][j][k] += MOs(l,a)(i,d) * t2(c,b)(k,l) * t1(d,j);

                                    t3_new[a][b][c][i][j][k] -= MOs(l,b)(i,d) * t2(a,c)(j,l) * t1(d,k);
                                    t3_new[a][b][c][i][j][k] -= MOs(l,b)(k,d) * t2(a,c)(i,l) * t1(d,j);
                                    t3_new[a][b][c][i][j][k] -= MOs(l,b)(j,d) * t2(a,c)(k,l) * t1(d,i);
                                    t3_new[a][b][c][i][j][k] += MOs(l,b)(k,d) * t2(a,c)(j,l) * t1(d,i);
                                    t3_new[a][b][c][i][j][k] += MOs(l,b)(i,d) * t2(a,c)(k,l) * t1(d,j);
                                }
                            }

                            // 15: - P(ab/c) f_ld t_abdijk t_cl
                            for (int d = n_Electrons; d < Matrix_Size; d++)
                            {
                                for (int l = 0; l < n_Electrons; l++)
                                {
                                    t3_new[a][b][c][i][j][k] -= fs(l,d) * t3[a][b][d][i][j][k] * t1(c,l);
                                    t3_new[a][b][c][i][j][k] += fs(l,d) * t3[c][b][d][i][j][k] * t1(a,l);
                                    t3_new[a][b][c][i][j][k] += fs(l,d) * t3[a][c][d][i][j][k] * t1(b,l);
                                }
                            }

                            // 16: - P(i/jk) f_ld t_abcljk t_di
                            for (int d = n_Electrons; d < Matrix_Size; d++)
                            {
                                for (int l = 0; l < n_Electrons; l++)
                                {
                                    t3_new[a][b][c][i][j][k] -= fs(l,d) * t3[a][b][c][l][j][k] * t1(d,i);
                                    t3_new[a][b][c][i][j][k] += fs(l,d) * t3[a][b][c][l][i][k] * t1(d,j);
                                    t3_new[a][b][c][i][j][k] += fs(l,d) * t3[a][b][c][l][j][i] * t1(d,k);
                                }
                            }

                            // 17: + P(ab/c) <cl||de> t_abdijk t_el
                            for (int d = n_Electrons; d < Matrix_Size; d++)
                            {
                                for (int l = 0; l < n_Electrons; l++)
                                {
                                    for (int e = n_Electrons; e < Matrix_Size; e++)
                                    {
                                        t3_new[a][b][c][i][j][k] += MOs(c,l)(d,e) * t3[a][b][d][i][j][k] * t1(e,l);
                                        t3_new[a][b][c][i][j][k] -= MOs(a,l)(d,e) * t3[c][b][d][i][j][k] * t1(e,l);
                                        t3_new[a][b][c][i][j][k] -= MOs(b,l)(d,e) * t3[a][c][d][i][j][k] * t1(e,l);
                                    }
                                }
                            }

                            // 18: - P(ij/k) <lm||kd> t_abcijl t_dm
                            for (int d = n_Electrons; d < Matrix_Size; d++)
                            {
                                for (int l = 0; l < n_Electrons; l++)
                                {
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        t3_new[a][b][c][i][j][k] -= MOs(l,m)(k,d) * t3[a][b][c][i][j][l] * t1(d,m);
                                        t3_new[a][b][c][i][j][k] += MOs(l,m)(i,d) * t3[a][b][c][k][j][l] * t1(d,m);
                                        t3_new[a][b][c][i][j][k] += MOs(l,m)(j,d) * t3[a][b][c][i][k][l] * t1(d,m);
                                    }
                                }
                            }

                            // 19: - P(ab/c) P(ij/k) <lm||dk> t_abdijl t_cm
                            for (int d = n_Electrons; d < Matrix_Size; d++)
                            {
                                for (int l = 0; l < n_Electrons; l++)
                                {
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        t3_new[a][b][c][i][j][k] -= MOs(l,m)(d,k) * t3[a][b][d][i][j][l] * t1(c,m);
                                        t3_new[a][b][c][i][j][k] += MOs(l,m)(d,k) * t3[c][b][d][i][j][l] * t1(a,m);
                                        t3_new[a][b][c][i][j][k] += MOs(l,m)(d,k) * t3[a][c][d][i][j][l] * t1(b,m);

                                        t3_new[a][b][c][i][j][k] += MOs(l,m)(d,i) * t3[a][b][d][k][j][l] * t1(c,m);
                                        t3_new[a][b][c][i][j][k] += MOs(l,m)(d,j) * t3[a][b][d][i][k][l] * t1(c,m);

                                        t3_new[a][b][c][i][j][k] -= MOs(l,m)(d,i) * t3[c][b][d][k][j][l] * t1(a,m);
                                        t3_new[a][b][c][i][j][k] -= MOs(l,m)(d,j) * t3[c][b][d][i][k][l] * t1(a,m);

                                        t3_new[a][b][c][i][j][k] -= MOs(l,m)(d,i) * t3[a][c][d][k][j][l] * t1(b,m);
                                        t3_new[a][b][c][i][j][k] -= MOs(l,m)(d,j) * t3[a][c][d][i][k][l] * t1(b,m);
                                    }
                                }
                            }

                            // 20: + P(ab/c) P(ij/k) <lc||de> t_abdijl t_ek
                            for (int d = n_Electrons; d < Matrix_Size; d++)
                            {
                                for (int l = 0; l < n_Electrons; l++)
                                {
                                    for (int e = n_Electrons; e < Matrix_Size; e++)
                                    {
                                        t3_new[a][b][c][i][j][k] += MOs(l,c)(d,e) * t3[a][b][d][i][j][l] * t1(e,k);
                                        t3_new[a][b][c][i][j][k] -= MOs(l,a)(d,e) * t3[c][b][d][i][j][l] * t1(e,k);
                                        t3_new[a][b][c][i][j][k] -= MOs(l,b)(d,e) * t3[a][c][d][i][j][l] * t1(e,k);

                                        t3_new[a][b][c][i][j][k] -= MOs(l,c)(d,e) * t3[a][b][d][k][j][l] * t1(e,i);
                                        t3_new[a][b][c][i][j][k] -= MOs(l,c)(d,e) * t3[a][b][d][i][k][l] * t1(e,j);

                                        t3_new[a][b][c][i][j][k] += MOs(l,c)(d,e) * t3[c][b][d][k][j][l] * t1(e,i);
                                        t3_new[a][b][c][i][j][k] += MOs(l,c)(d,e) * t3[c][b][d][i][k][l] * t1(e,j);

                                        t3_new[a][b][c][i][j][k] += MOs(l,c)(d,e) * t3[a][c][d][k][j][l] * t1(e,i);
                                        t3_new[a][b][c][i][j][k] += MOs(l,c)(d,e) * t3[a][c][d][i][k][l] * t1(e,j);
                                    }
                                }
                            }

                            // 21: - 1/2 P(abc) <bl||de> t_adeijk t_cl
                            for (int d = n_Electrons; d < Matrix_Size; d++)
                            {
                                for (int l = 0; l < n_Electrons; l++)
                                {
                                    for (int e = n_Electrons; e < Matrix_Size; e++)
                                    {
                                        t3_new[a][b][c][i][j][k] -= 0.5 * MOs(b,l)(d,e) * t3[a][d][e][i][j][k] * t1(c,l);
                                        t3_new[a][b][c][i][j][k] += 0.5 * MOs(a,l)(d,e) * t3[b][d][e][i][j][k] * t1(c,l);
                                        t3_new[a][b][c][i][j][k] += 0.5 * MOs(b,l)(d,e) * t3[c][d][e][i][j][k] * t1(a,l);
                                        t3_new[a][b][c][i][j][k] += 0.5 * MOs(c,l)(d,e) * t3[a][d][e][i][j][k] * t1(b,l);
                                        t3_new[a][b][c][i][j][k] -= 0.5 * MOs(a,l)(d,e) * t3[c][d][e][i][j][k] * t1(b,l);
                                        t3_new[a][b][c][i][j][k] -= 0.5 * MOs(c,l)(d,e) * t3[b][d][e][i][j][k] * t1(a,l);
                                    }
                                }
                            }

                            // 22: + 1/2 P(ijk) <lm||dj> t_abclmk t_di
                            for (int d = n_Electrons; d < Matrix_Size; d++)
                            {
                                for (int l = 0; l < n_Electrons; l++)
                                {
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        t3_new[a][b][c][i][j][k] += 0.5 * MOs(l,m)(d,j) * t3[a][b][c][l][m][k] * t1(d,i);
                                        t3_new[a][b][c][i][j][k] -= 0.5 * MOs(l,m)(d,i) * t3[a][b][c][l][m][k] * t1(d,j);
                                        t3_new[a][b][c][i][j][k] -= 0.5 * MOs(l,m)(d,k) * t3[a][b][c][l][m][j] * t1(d,i);
                                        t3_new[a][b][c][i][j][k] -= 0.5 * MOs(l,m)(d,j) * t3[a][b][c][l][m][i] * t1(d,j);
                                        t3_new[a][b][c][i][j][k] += 0.5 * MOs(l,m)(d,i) * t3[a][b][c][l][m][j] * t1(d,k);
                                        t3_new[a][b][c][i][j][k] += 0.5 * MOs(l,m)(d,k) * t3[a][b][c][l][m][i] * t1(d,j);
                                    }
                                }
                            }

                            // 23: + <lm||de> t_abcdijkl t_em
                            for (int d = n_Electrons; d < Matrix_Size; d++)
                            {
                                for (int l = 0; l < n_Electrons; l++)
                                {
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        for (int e = n_Electrons; e < Matrix_Size; e++)
                                        {
                                            t3_new[a][b][c][i][j][k] += MOs(l,m)(d,e) * t4[a][b][c][d][i][j][k][l];
                                        }
                                    }
                                }
                            }

                            // 24: - 1/2 P(ab/c) <lm||de> t_abdeijlk t_cm
                            for (int d = n_Electrons; d < Matrix_Size; d++)
                            {
                                for (int l = 0; l < n_Electrons; l++)
                                {
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        for (int e = n_Electrons; e < Matrix_Size; e++)
                                        {
                                            t3_new[a][b][c][i][j][k] -= 0.5 * MOs(l,m)(d,e) * t4[a][b][d][e][i][j][l][k] * t1(c,m);
                                            t3_new[a][b][c][i][j][k] += 0.5 * MOs(l,m)(d,e) * t4[c][b][d][e][i][j][l][k] * t1(a,m);
                                            t3_new[a][b][c][i][j][k] += 0.5 * MOs(l,m)(d,e) * t4[a][c][d][e][i][j][l][k] * t1(b,m);
                                        }
                                    }
                                }
                            }

                            // 25: - 1/2 P(i/jk) <lm||de> t_aebclmjk t_di
                            for (int d = n_Electrons; d < Matrix_Size; d++)
                            {
                                for (int l = 0; l < n_Electrons; l++)
                                {
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        for (int e = n_Electrons; e < Matrix_Size; e++)
                                        {
                                            t3_new[a][b][c][i][j][k] -= 0.5 * MOs(l,m)(d,e) * t4[a][e][b][c][l][m][j][k] * t1(d,i);
                                            t3_new[a][b][c][i][j][k] += 0.5 * MOs(l,m)(d,e) * t4[a][e][b][c][l][m][i][k] * t1(d,j);
                                            t3_new[a][b][c][i][j][k] += 0.5 * MOs(l,m)(d,e) * t4[a][e][b][c][l][m][j][i] * t1(d,k);
                                        }
                                    }
                                }
                            }

                            // 26: - P(a/bc) P(ij/k) f_ld t_adij t_bclk
                            for (int d = n_Electrons; d < Matrix_Size; d++)
                            {
                                for (int l = 0; l < n_Electrons; l++)
                                {
                                    t3_new[a][b][c][i][j][k] -= fs(l,d) * t2(a,d)(i,j) * t2(b,c)(l,k);
                                    t3_new[a][b][c][i][j][k] += fs(l,d) * t2(b,d)(i,j) * t2(a,c)(l,k);
                                    t3_new[a][b][c][i][j][k] += fs(l,d) * t2(c,d)(i,j) * t2(b,a)(l,k);

                                    t3_new[a][b][c][i][j][k] += fs(l,d) * t2(a,d)(k,j) * t2(b,c)(l,i);
                                    t3_new[a][b][c][i][j][k] += fs(l,d) * t2(a,d)(i,k) * t2(b,c)(l,j);

                                    t3_new[a][b][c][i][j][k] -= fs(l,d) * t2(b,d)(k,j) * t2(a,c)(l,i);
                                    t3_new[a][b][c][i][j][k] -= fs(l,d) * t2(b,d)(i,k) * t2(a,c)(l,j);

                                    t3_new[a][b][c][i][j][k] -= fs(l,d) * t2(c,d)(k,j) * t2(b,a)(l,i);
                                    t3_new[a][b][c][i][j][k] -= fs(l,d) * t2(c,d)(i,k) * t2(b,a)(l,j);
                                }
                            }

                            // 27: - P(ab/c) P(ijk) <lm||jd> t_abil t_dcmk
                            for (int d = n_Electrons; d < Matrix_Size; d++)
                            {
                                for (int l = 0; l < n_Electrons; l++)
                                {
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        t3_new[a][b][c][i][j][k] -= MOs(l,m)(j,d) * t2(a,b)(i,l) * t2(d,c)(m,k);
                                        t3_new[a][b][c][i][j][k] += MOs(l,m)(k,d) * t2(a,b)(i,l) * t2(d,c)(m,j);
                                        t3_new[a][b][c][i][j][k] += MOs(l,m)(i,d) * t2(a,b)(j,l) * t2(d,c)(m,k);
                                        t3_new[a][b][c][i][j][k] += MOs(l,m)(j,d) * t2(a,b)(k,l) * t2(d,c)(m,i);
                                        t3_new[a][b][c][i][j][k] -= MOs(l,m)(k,d) * t2(a,b)(j,l) * t2(d,c)(m,i);
                                        t3_new[a][b][c][i][j][k] -= MOs(l,m)(i,d) * t2(a,b)(k,l) * t2(d,c)(m,j);

                                        t3_new[a][b][c][i][j][k] += MOs(l,m)(j,d) * t2(c,b)(i,l) * t2(d,a)(m,k);
                                        t3_new[a][b][c][i][j][k] -= MOs(l,m)(k,d) * t2(c,b)(i,l) * t2(d,a)(m,j);
                                        t3_new[a][b][c][i][j][k] -= MOs(l,m)(i,d) * t2(c,b)(j,l) * t2(d,a)(m,k);
                                        t3_new[a][b][c][i][j][k] -= MOs(l,m)(j,d) * t2(c,b)(k,l) * t2(d,a)(m,i);
                                        t3_new[a][b][c][i][j][k] += MOs(l,m)(k,d) * t2(c,b)(j,l) * t2(d,a)(m,i);
                                        t3_new[a][b][c][i][j][k] += MOs(l,m)(i,d) * t2(c,b)(k,l) * t2(d,a)(m,j);

                                        t3_new[a][b][c][i][j][k] += MOs(l,m)(j,d) * t2(a,c)(i,l) * t2(d,b)(m,k);
                                        t3_new[a][b][c][i][j][k] -= MOs(l,m)(k,d) * t2(a,c)(i,l) * t2(d,b)(m,j);
                                        t3_new[a][b][c][i][j][k] -= MOs(l,m)(i,d) * t2(a,c)(j,l) * t2(d,b)(m,k);
                                        t3_new[a][b][c][i][j][k] -= MOs(l,m)(j,d) * t2(a,c)(k,l) * t2(d,b)(m,i);
                                        t3_new[a][b][c][i][j][k] += MOs(l,m)(k,d) * t2(a,c)(j,l) * t2(d,b)(m,i);
                                        t3_new[a][b][c][i][j][k] += MOs(l,m)(i,d) * t2(a,c)(k,l) * t2(d,b)(m,j);
                                    }
                                }
                            }

                            // 28: + P(abc) P(ij/k) <bl||de> t_adij t_eclk
                            for (int d = n_Electrons; d < Matrix_Size; d++)
                            {
                                for (int l = 0; l < n_Electrons; l++)
                                {
                                    for (int e = n_Electrons; e < Matrix_Size; e++)
                                    {
                                        t3_new[a][b][c][i][j][k] += MOs(b,l)(d,e) * t2(a,d)(i,j) * t2(e,c)(l,k);
                                        t3_new[a][b][c][i][j][k] -= MOs(a,l)(d,e) * t2(b,d)(i,j) * t2(e,c)(l,k);
                                        t3_new[a][b][c][i][j][k] -= MOs(c,l)(d,e) * t2(a,d)(i,j) * t2(e,b)(l,k);
                                        t3_new[a][b][c][i][j][k] -= MOs(b,l)(d,e) * t2(c,d)(i,j) * t2(e,a)(l,k);
                                        t3_new[a][b][c][i][j][k] += MOs(c,l)(d,e) * t2(b,d)(i,j) * t2(e,a)(l,k);
                                        t3_new[a][b][c][i][j][k] += MOs(a,l)(d,e) * t2(c,d)(i,j) * t2(e,b)(l,k);

                                        t3_new[a][b][c][i][j][k] -= MOs(b,l)(d,e) * t2(a,d)(k,j) * t2(e,c)(l,i);
                                        t3_new[a][b][c][i][j][k] += MOs(a,l)(d,e) * t2(b,d)(k,j) * t2(e,c)(l,i);
                                        t3_new[a][b][c][i][j][k] += MOs(c,l)(d,e) * t2(a,d)(k,j) * t2(e,b)(l,i);
                                        t3_new[a][b][c][i][j][k] += MOs(b,l)(d,e) * t2(c,d)(k,j) * t2(e,a)(l,i);
                                        t3_new[a][b][c][i][j][k] -= MOs(c,l)(d,e) * t2(b,d)(k,j) * t2(e,a)(l,i);
                                        t3_new[a][b][c][i][j][k] -= MOs(a,l)(d,e) * t2(c,d)(k,j) * t2(e,b)(l,i);

                                        t3_new[a][b][c][i][j][k] -= MOs(b,l)(d,e) * t2(a,d)(i,k) * t2(e,c)(l,j);
                                        t3_new[a][b][c][i][j][k] += MOs(a,l)(d,e) * t2(b,d)(i,k) * t2(e,c)(l,j);
                                        t3_new[a][b][c][i][j][k] += MOs(c,l)(d,e) * t2(a,d)(i,k) * t2(e,b)(l,j);
                                        t3_new[a][b][c][i][j][k] += MOs(b,l)(d,e) * t2(c,d)(i,k) * t2(e,a)(l,j);
                                        t3_new[a][b][c][i][j][k] -= MOs(c,l)(d,e) * t2(b,d)(i,k) * t2(e,a)(l,j);
                                        t3_new[a][b][c][i][j][k] -= MOs(a,l)(d,e) * t2(c,d)(i,k) * t2(e,b)(l,j);
                                    }
                                }
                            }

                            // 29: - 1/2 P(a/bc) P(ij/k) <al||de> t_deij t_bclk
                            for (int d = n_Electrons; d < Matrix_Size; d++)
                            {
                                for (int l = 0; l < n_Electrons; l++)
                                {
                                    for (int e = n_Electrons; e < Matrix_Size; e++)
                                    {
                                        t3_new[a][b][c][i][j][k] -= 0.5 * MOs(a,l)(d,e) * t2(d,e)(i,j) * t2(b,c)(l,k);
                                        t3_new[a][b][c][i][j][k] += 0.5 * MOs(b,l)(d,e) * t2(d,e)(i,j) * t2(a,c)(l,k);
                                        t3_new[a][b][c][i][j][k] += 0.5 * MOs(c,l)(d,e) * t2(d,e)(i,j) * t2(b,c)(l,k);

                                        t3_new[a][b][c][i][j][k] += 0.5 * MOs(a,l)(d,e) * t2(d,e)(k,j) * t2(b,c)(l,i);
                                        t3_new[a][b][c][i][j][k] += 0.5 * MOs(a,l)(d,e) * t2(d,e)(i,k) * t2(b,c)(l,j);

                                        t3_new[a][b][c][i][j][k] -= 0.5 * MOs(b,l)(d,e) * t2(d,e)(k,j) * t2(a,c)(l,i);
                                        t3_new[a][b][c][i][j][k] -= 0.5 * MOs(b,l)(d,e) * t2(d,e)(i,k) * t2(a,c)(l,j);

                                        t3_new[a][b][c][i][j][k] -= 0.5 * MOs(c,l)(d,e) * t2(d,e)(k,j) * t2(b,c)(l,i);
                                        t3_new[a][b][c][i][j][k] -= 0.5 * MOs(c,l)(d,e) * t2(d,e)(i,k) * t2(b,c)(l,j);
                                    }
                                }
                            }

                            // 30: + 1/2 P(a/bc) P(ij/k) <lm||dk> t_adij t_bclm
                            for (int d = n_Electrons; d < Matrix_Size; d++)
                            {
                                for (int l = 0; l < n_Electrons; l++)
                                {
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        t3_new[a][b][c][i][j][k] += 0.5 * MOs(l,m)(d,k) * t2(a,d)(i,j) * t2(b,c)(l,m);
                                        t3_new[a][b][c][i][j][k] -= 0.5 * MOs(l,m)(d,k) * t2(b,d)(i,j) * t2(a,c)(l,m);
                                        t3_new[a][b][c][i][j][k] -= 0.5 * MOs(l,m)(d,k) * t2(c,d)(i,j) * t2(b,a)(l,m);

                                        t3_new[a][b][c][i][j][k] -= 0.5 * MOs(l,m)(d,i) * t2(a,d)(k,j) * t2(b,c)(l,m);
                                        t3_new[a][b][c][i][j][k] -= 0.5 * MOs(l,m)(d,j) * t2(a,d)(i,k) * t2(b,c)(l,m);

                                        t3_new[a][b][c][i][j][k] += 0.5 * MOs(l,m)(d,i) * t2(b,d)(k,j) * t2(a,c)(l,m);
                                        t3_new[a][b][c][i][j][k] += 0.5 * MOs(l,m)(d,j) * t2(b,d)(i,k) * t2(a,c)(l,m);

                                        t3_new[a][b][c][i][j][k] += 0.5 * MOs(l,m)(d,i) * t2(c,d)(k,j) * t2(b,a)(l,m);
                                        t3_new[a][b][c][i][j][k] += 0.5 * MOs(l,m)(d,j) * t2(c,d)(i,k) * t2(b,a)(l,m);
                                    }
                                }
                            }

                            // 31: + 1/4 P(a/bc) <lm||de> t_adeijk t_bclm
                            for (int d = n_Electrons; d < Matrix_Size; d++)
                            {
                                for (int l = 0; l < n_Electrons; l++)
                                {
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        for (int e = n_Electrons; e < Matrix_Size; e++)
                                        {
                                            t3_new[a][b][c][i][j][k] += 0.25 * MOs(l,m)(d,e) * t3[a][d][e][i][j][k] * t2(b,c)(l,m);
                                            t3_new[a][b][c][i][j][k] -= 0.25 * MOs(l,m)(d,e) * t3[b][d][e][i][j][k] * t2(a,c)(l,m);
                                            t3_new[a][b][c][i][j][k] -= 0.25 * MOs(l,m)(d,e) * t3[c][d][e][i][j][k] * t2(b,a)(l,m);
                                        }
                                    }
                                }
                            }

                            // 32: + 1/4 P(ij/k) <lm||de> t_abclmk t_deij
                            for (int d = n_Electrons; d < Matrix_Size; d++)
                            {
                                for (int l = 0; l < n_Electrons; l++)
                                {
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        for (int e = n_Electrons; e < Matrix_Size; e++)
                                        {
                                            t3_new[a][b][c][i][j][k] += 0.25 * MOs(l,m)(d,e) * t2(d,e)(i,j) * t3[a][b][c][l][m][k];
                                            t3_new[a][b][c][i][j][k] -= 0.25 * MOs(l,m)(d,e) * t2(d,e)(k,j) * t3[a][b][c][l][m][i];
                                            t3_new[a][b][c][i][j][k] -= 0.25 * MOs(l,m)(d,e) * t2(d,e)(i,k) * t3[a][b][c][l][m][j];
                                        }
                                    }
                                }
                            }

                            // 33: - 1/2 P(a/bc) P(ij/k) <lm||de> t_adeilj t_bcmk
                            for (int d = n_Electrons; d < Matrix_Size; d++)
                            {
                                for (int l = 0; l < n_Electrons; l++)
                                {
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        for (int e = n_Electrons; e < Matrix_Size; e++)
                                        {
                                            t3_new[a][b][c][i][j][k] -= 0.5 * MOs(l,m)(d,e) * t3[a][d][e][i][l][j] * t2(b,c)(m,k);
                                            t3_new[a][b][c][i][j][k] += 0.5 * MOs(l,m)(d,e) * t3[b][d][e][i][l][j] * t2(a,c)(m,k);
                                            t3_new[a][b][c][i][j][k] += 0.5 * MOs(l,m)(d,e) * t3[c][d][e][i][l][j] * t2(b,a)(m,k);

                                            t3_new[a][b][c][i][j][k] += 0.5 * MOs(l,m)(d,e) * t3[a][d][e][k][l][j] * t2(b,c)(m,i);
                                            t3_new[a][b][c][i][j][k] += 0.5 * MOs(l,m)(d,e) * t3[a][d][e][i][l][k] * t2(b,c)(m,j);

                                            t3_new[a][b][c][i][j][k] -= 0.5 * MOs(l,m)(d,e) * t3[b][d][e][k][l][j] * t2(a,c)(m,i);
                                            t3_new[a][b][c][i][j][k] -= 0.5 * MOs(l,m)(d,e) * t3[b][d][e][i][l][k] * t2(a,c)(m,j);

                                            t3_new[a][b][c][i][j][k] -= 0.5 * MOs(l,m)(d,e) * t3[c][d][e][k][l][j] * t2(b,a)(m,i);
                                            t3_new[a][b][c][i][j][k] -= 0.5 * MOs(l,m)(d,e) * t3[c][d][e][i][l][k] * t2(b,a)(m,j);
                                        }
                                    }
                                }
                            }

                            // 34: - 1/2 P(a/bc) P(ij/k) <lm||de> t_beclmk t_adij
                            for (int d = n_Electrons; d < Matrix_Size; d++)
                            {
                                for (int l = 0; l < n_Electrons; l++)
                                {
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        for (int e = n_Electrons; e < Matrix_Size; e++)
                                        {
                                            t3_new[a][b][c][i][j][k] -= 0.5 * MOs(l,m)(d,e) * t3[b][e][c][l][m][k] * t2(a,d)(i,j);
                                            t3_new[a][b][c][i][j][k] += 0.5 * MOs(l,m)(d,e) * t3[a][e][c][l][m][k] * t2(b,d)(i,j);
                                            t3_new[a][b][c][i][j][k] += 0.5 * MOs(l,m)(d,e) * t3[b][e][a][l][m][k] * t2(c,d)(i,j);

                                            t3_new[a][b][c][i][j][k] += 0.5 * MOs(l,m)(d,e) * t3[b][e][c][l][m][i] * t2(a,d)(k,j);
                                            t3_new[a][b][c][i][j][k] += 0.5 * MOs(l,m)(d,e) * t3[b][e][c][l][m][j] * t2(a,d)(i,k);

                                            t3_new[a][b][c][i][j][k] -= 0.5 * MOs(l,m)(d,e) * t3[a][e][c][l][m][i] * t2(b,d)(k,j);
                                            t3_new[a][b][c][i][j][k] -= 0.5 * MOs(l,m)(d,e) * t3[a][e][c][l][m][j] * t2(b,d)(i,k);

                                            t3_new[a][b][c][i][j][k] -= 0.5 * MOs(l,m)(d,e) * t3[b][e][a][l][m][i] * t2(b,d)(k,j);
                                            t3_new[a][b][c][i][j][k] -= 0.5 * MOs(l,m)(d,e) * t3[b][e][a][l][m][j] * t2(b,d)(i,k);
                                        }
                                    }
                                }
                            }

                            // 35: - 1/2 P(ab/c) <lm||de> t_abdijk t_celm
                            for (int d = n_Electrons; d < Matrix_Size; d++)
                            {
                                for (int l = 0; l < n_Electrons; l++)
                                {
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        for (int e = n_Electrons; e < Matrix_Size; e++)
                                        {
                                            t3_new[a][b][c][i][j][k] -= 0.5 * MOs(l,m)(d,e) * t3[a][b][d][i][j][k] * t2(c,e)(l,m);
                                            t3_new[a][b][c][i][j][k] += 0.5 * MOs(l,m)(d,e) * t3[c][b][d][i][j][k] * t2(a,e)(l,m);
                                            t3_new[a][b][c][i][j][k] += 0.5 * MOs(l,m)(d,e) * t3[a][c][d][i][j][k] * t2(b,e)(l,m);
                                        }
                                    }
                                }
                            }

                            // 36: - 1/2 P(i/jk) <lm||de> t_abcmjk t_deli
                            for (int d = n_Electrons; d < Matrix_Size; d++)
                            {
                                for (int l = 0; l < n_Electrons; l++)
                                {
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        for (int e = n_Electrons; e < Matrix_Size; e++)
                                        {
                                            t3_new[a][b][c][i][j][k] -= 0.5 * MOs(l,m)(d,e) * t3[a][b][c][m][j][k] * t2(d,e)(l,i);
                                            t3_new[a][b][c][i][j][k] += 0.5 * MOs(l,m)(d,e) * t3[a][b][c][m][i][k] * t2(d,e)(l,j);
                                            t3_new[a][b][c][i][j][k] += 0.5 * MOs(l,m)(d,e) * t3[a][b][c][m][j][i] * t2(d,e)(l,k);
                                        }
                                    }
                                }
                            }

                            // 37: + P(ab/c) P(ij/k) <lm||de> t_abdijl t_ecmk
                            for (int d = n_Electrons; d < Matrix_Size; d++)
                            {
                                for (int l = 0; l < n_Electrons; l++)
                                {
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        for (int e = n_Electrons; e < Matrix_Size; e++)
                                        {
                                            t3_new[a][b][c][i][j][k] += MOs(l,m)(d,e) * t3[a][b][d][i][j][l] * t2(e,c)(m,k);
                                            t3_new[a][b][c][i][j][k] -= MOs(l,m)(d,e) * t3[c][b][d][i][j][l] * t2(e,a)(m,k);
                                            t3_new[a][b][c][i][j][k] -= MOs(l,m)(d,e) * t3[a][c][d][i][j][l] * t2(e,b)(m,k);

                                            t3_new[a][b][c][i][j][k] -= MOs(l,m)(d,e) * t3[a][b][d][k][j][l] * t2(e,c)(m,i);
                                            t3_new[a][b][c][i][j][k] -= MOs(l,m)(d,e) * t3[a][b][d][i][k][l] * t2(e,c)(m,j);

                                            t3_new[a][b][c][i][j][k] += MOs(l,m)(d,e) * t3[c][b][d][k][j][l] * t2(e,a)(m,i);
                                            t3_new[a][b][c][i][j][k] += MOs(l,m)(d,e) * t3[c][b][d][i][k][l] * t2(e,a)(m,j);

                                            t3_new[a][b][c][i][j][k] += MOs(l,m)(d,e) * t3[a][c][d][k][j][l] * t2(e,b)(m,i);
                                            t3_new[a][b][c][i][j][k] += MOs(l,m)(d,e) * t3[a][c][d][i][k][l] * t2(e,b)(m,j);
                                        }
                                    }
                                }
                            }

                            // 38: + P(ab/c) P(ijk) <lm||jd> t_abil t_dk t_cm BUGGY
                            for (int d = n_Electrons; d < Matrix_Size; d++)
                            {
                                for (int l = 0; l < n_Electrons; l++)
                                {
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        t3_new[a][b][c][i][j][k] += MOs(l,m)(j,d) * t2(a,b)(i,l) * t1(d,k) * t1(c,m);
                                        t3_new[a][b][c][i][j][k] -= MOs(l,m)(j,d) * t2(a,b)(k,l) * t1(d,i) * t1(c,m);
                                        t3_new[a][b][c][i][j][k] -= MOs(l,m)(i,d) * t2(a,b)(j,l) * t1(d,k) * t1(c,m);
                                        t3_new[a][b][c][i][j][k] -= MOs(l,m)(k,d) * t2(a,b)(i,l) * t1(d,j) * t1(c,m);
                                        t3_new[a][b][c][i][j][k] += MOs(l,m)(i,d) * t2(a,b)(k,l) * t1(d,j) * t1(c,m);
                                        t3_new[a][b][c][i][j][k] += MOs(l,m)(k,d) * t2(a,b)(j,l) * t1(d,i) * t1(c,m);

                                        t3_new[a][b][c][i][j][k] -= MOs(l,m)(j,d) * t2(c,b)(i,l) * t1(d,k) * t1(a,m);
                                        t3_new[a][b][c][i][j][k] += MOs(l,m)(j,d) * t2(c,b)(k,l) * t1(d,i) * t1(a,m);
                                        t3_new[a][b][c][i][j][k] += MOs(l,m)(i,d) * t2(c,b)(j,l) * t1(d,k) * t1(a,m);
                                        t3_new[a][b][c][i][j][k] += MOs(l,m)(k,d) * t2(c,b)(i,l) * t1(d,j) * t1(a,m);
                                        t3_new[a][b][c][i][j][k] -= MOs(l,m)(i,d) * t2(c,b)(k,l) * t1(d,j) * t1(a,m);
                                        t3_new[a][b][c][i][j][k] -= MOs(l,m)(k,d) * t2(c,b)(j,l) * t1(d,i) * t1(a,m);

                                        t3_new[a][b][c][i][j][k] -= MOs(l,m)(j,d) * t2(a,c)(i,l) * t1(d,k) * t1(b,m);
                                        t3_new[a][b][c][i][j][k] += MOs(l,m)(j,d) * t2(a,c)(k,l) * t1(d,i) * t1(b,m);
                                        t3_new[a][b][c][i][j][k] += MOs(l,m)(i,d) * t2(a,c)(j,l) * t1(d,k) * t1(b,m);
                                        t3_new[a][b][c][i][j][k] += MOs(l,m)(k,d) * t2(a,c)(i,l) * t1(d,j) * t1(b,m);
                                        t3_new[a][b][c][i][j][k] -= MOs(l,m)(i,d) * t2(a,c)(k,l) * t1(d,j) * t1(b,m);
                                        t3_new[a][b][c][i][j][k] -= MOs(l,m)(k,d) * t2(a,c)(j,l) * t1(d,i) * t1(b,m);

                                    }
                                }
                            }

                            // 39: - P(abc) P(ij/k) <bl||de> t_adij t_ek t_cl
                            for (int d = n_Electrons; d < Matrix_Size; d++)
                            {
                                for (int l = 0; l < n_Electrons; l++)
                                {
                                    for (int e = n_Electrons; e < Matrix_Size; e++)
                                    {
                                        t3_new[a][b][c][i][j][k] -= MOs(b,l)(d,e) * t2(a,d)(i,j) * t1(e,k) * t1(c,l);
                                        t3_new[a][b][c][i][j][k] += MOs(a,l)(d,e) * t2(b,d)(i,j) * t1(e,k) * t1(c,l);
                                        t3_new[a][b][c][i][j][k] += MOs(c,l)(d,e) * t2(a,d)(i,j) * t1(e,k) * t1(b,l);
                                        t3_new[a][b][c][i][j][k] += MOs(b,l)(d,e) * t2(c,d)(i,j) * t1(e,k) * t1(a,l);
                                        t3_new[a][b][c][i][j][k] -= MOs(c,l)(d,e) * t2(b,d)(i,j) * t1(e,k) * t1(a,l);
                                        t3_new[a][b][c][i][j][k] -= MOs(a,l)(d,e) * t2(c,d)(i,j) * t1(e,k) * t1(b,l);

                                        t3_new[a][b][c][i][j][k] += MOs(b,l)(d,e) * t2(a,d)(k,j) * t1(e,i) * t1(c,l);
                                        t3_new[a][b][c][i][j][k] -= MOs(a,l)(d,e) * t2(b,d)(k,j) * t1(e,i) * t1(c,l);
                                        t3_new[a][b][c][i][j][k] -= MOs(c,l)(d,e) * t2(a,d)(k,j) * t1(e,i) * t1(b,l);
                                        t3_new[a][b][c][i][j][k] -= MOs(b,l)(d,e) * t2(c,d)(k,j) * t1(e,i) * t1(a,l);
                                        t3_new[a][b][c][i][j][k] += MOs(c,l)(d,e) * t2(b,d)(k,j) * t1(e,i) * t1(a,l);
                                        t3_new[a][b][c][i][j][k] += MOs(a,l)(d,e) * t2(c,d)(k,j) * t1(e,i) * t1(b,l);

                                        t3_new[a][b][c][i][j][k] += MOs(b,l)(d,e) * t2(a,d)(i,k) * t1(e,j) * t1(c,l);
                                        t3_new[a][b][c][i][j][k] -= MOs(a,l)(d,e) * t2(b,d)(i,k) * t1(e,j) * t1(c,l);
                                        t3_new[a][b][c][i][j][k] -= MOs(c,l)(d,e) * t2(a,d)(i,k) * t1(e,j) * t1(b,l);
                                        t3_new[a][b][c][i][j][k] -= MOs(b,l)(d,e) * t2(c,d)(i,k) * t1(e,j) * t1(a,l);
                                        t3_new[a][b][c][i][j][k] += MOs(c,l)(d,e) * t2(b,d)(i,k) * t1(e,j) * t1(a,l);
                                        t3_new[a][b][c][i][j][k] += MOs(a,l)(d,e) * t2(c,d)(i,k) * t1(e,j) * t1(b,l);
                                    }
                                }
                            }

                            // 40: + P(a/bc) P(ij/k) <lm||dk> t_adij t_bl t_cm
                            for (int d = n_Electrons; d < Matrix_Size; d++)
                            {
                                for (int l = 0; l < n_Electrons; l++)
                                {
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        t3_new[a][b][c][i][j][k] += MOs(l,m)(d,k) * t2(a,d)(i,j) * t1(b,l) * t1(c,m);
                                        t3_new[a][b][c][i][j][k] -= MOs(l,m)(d,k) * t2(b,d)(i,j) * t1(a,l) * t1(c,m);
                                        t3_new[a][b][c][i][j][k] -= MOs(l,m)(d,k) * t2(c,d)(i,j) * t1(b,l) * t1(a,m);

                                        t3_new[a][b][c][i][j][k] -= MOs(l,m)(d,i) * t2(a,d)(k,j) * t1(b,l) * t1(c,m);
                                        t3_new[a][b][c][i][j][k] -= MOs(l,m)(d,j) * t2(a,d)(i,k) * t1(b,l) * t1(c,m);

                                        t3_new[a][b][c][i][j][k] += MOs(l,m)(d,i) * t2(b,d)(k,j) * t1(a,l) * t1(c,m);
                                        t3_new[a][b][c][i][j][k] += MOs(l,m)(d,j) * t2(b,d)(i,k) * t1(a,l) * t1(c,m);

                                        t3_new[a][b][c][i][j][k] += MOs(l,m)(d,i) * t2(c,d)(k,j) * t1(b,l) * t1(a,m);
                                        t3_new[a][b][c][i][j][k] += MOs(l,m)(d,j) * t2(c,d)(i,k) * t1(b,l) * t1(a,m);
                                    }
                                }
                            }

                            // 41: - P(a/bc) P(ij/k) <al||de> t_bclk t_di t_ej
                            for (int d = n_Electrons; d < Matrix_Size; d++)
                            {
                                for (int l = 0; l < n_Electrons; l++)
                                {
                                    for (int e = n_Electrons; e < Matrix_Size; e++)
                                    {
                                        t3_new[a][b][c][i][j][k] -= MOs(a,l)(d,e) * t2(b,c)(l,k) * t1(d,i) * t1(e,j);
                                        t3_new[a][b][c][i][j][k] += MOs(b,l)(d,e) * t2(a,c)(l,k) * t1(d,i) * t1(e,j);
                                        t3_new[a][b][c][i][j][k] += MOs(c,l)(d,e) * t2(b,a)(l,k) * t1(d,i) * t1(e,j);

                                        t3_new[a][b][c][i][j][k] += MOs(a,l)(d,e) * t2(b,c)(l,i) * t1(d,k) * t1(e,j);
                                        t3_new[a][b][c][i][j][k] += MOs(a,l)(d,e) * t2(b,c)(l,j) * t1(d,i) * t1(e,k);

                                        t3_new[a][b][c][i][j][k] -= MOs(b,l)(d,e) * t2(a,c)(l,i) * t1(d,k) * t1(e,j);
                                        t3_new[a][b][c][i][j][k] -= MOs(b,l)(d,e) * t2(a,c)(l,j) * t1(d,i) * t1(e,k);

                                        t3_new[a][b][c][i][j][k] -= MOs(c,l)(d,e) * t2(b,a)(l,i) * t1(d,k) * t1(e,j);
                                        t3_new[a][b][c][i][j][k] -= MOs(c,l)(d,e) * t2(b,a)(l,j) * t1(d,i) * t1(d,j);
                                    }
                                }
                            }

                            // 42: + 1/2 P(a/bc) <lm||de> t_adeijk t_bl t_cm
                            for (int d = n_Electrons; d < Matrix_Size; d++)
                            {
                                for (int l = 0; l < n_Electrons; l++)
                                {
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        for (int e = n_Electrons; e < Matrix_Size; e++)
                                        {
                                            t3_new[a][b][c][i][j][k] += 0.5 * MOs(l,m)(d,e) * t3[a][d][e][i][j][k] * t1(b,l) * t1(c,m);
                                            t3_new[a][b][c][i][j][k] -= 0.5 * MOs(l,m)(d,e) * t3[b][d][e][i][j][k] * t1(a,l) * t1(c,m);
                                            t3_new[a][b][c][i][j][k] -= 0.5 * MOs(l,m)(d,e) * t3[c][d][e][i][j][k] * t1(b,l) * t1(a,m);
                                        }
                                    }
                                }
                            }

                            // 43: + 1/2 P(ij/k) <lm||de> t_abclmk t_di t_ej
                            for (int d = n_Electrons; d < Matrix_Size; d++)
                            {
                                for (int l = 0; l < n_Electrons; l++)
                                {
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        for (int e = n_Electrons; e < Matrix_Size; e++)
                                        {
                                            t3_new[a][b][c][i][j][k] += 0.5 * MOs(l,m)(d,e) * t3[a][b][c][l][m][k] * t1(d,i) * t1(e,j);
                                            t3_new[a][b][c][i][j][k] -= 0.5 * MOs(l,m)(d,e) * t3[a][b][c][l][m][i] * t1(d,k) * t1(e,j);
                                            t3_new[a][b][c][i][j][k] -= 0.5 * MOs(l,m)(d,e) * t3[a][b][c][l][m][j] * t1(d,i) * t1(e,k);
                                        }
                                    }
                                }
                            }

                            // 44: - P(ab/c) <lm||de> t_abdijk t_cl t_em
                            for (int d = n_Electrons; d < Matrix_Size; d++)
                            {
                                for (int l = 0; l < n_Electrons; l++)
                                {
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        for (int e = n_Electrons; e < Matrix_Size; e++)
                                        {
                                            t3_new[a][b][c][i][j][k] -= MOs(l,m)(d,e) * t3[a][b][d][i][j][k] * t1(c,l) * t1(e,m);
                                            t3_new[a][b][c][i][j][k] += MOs(l,m)(d,e) * t3[c][b][d][i][j][k] * t1(a,l) * t1(e,m);
                                            t3_new[a][b][c][i][j][k] += MOs(l,m)(d,e) * t3[a][c][d][i][j][k] * t1(b,l) * t1(e,m);
                                        }
                                    }
                                }
                            }

                            // 45: - P(i/jk) <lm||de> t_abcmjk t_dl t_ei
                            for (int d = n_Electrons; d < Matrix_Size; d++)
                            {
                                for (int l = 0; l < n_Electrons; l++)
                                {
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        for (int e = n_Electrons; e < Matrix_Size; e++)
                                        {
                                            t3_new[a][b][c][i][j][k] -= MOs(l,m)(d,e) * t3[a][b][c][m][j][k] * t1(e,i);
                                            t3_new[a][b][c][i][j][k] += MOs(l,m)(d,e) * t3[a][b][c][m][i][k] * t1(e,j);
                                            t3_new[a][b][c][i][j][k] += MOs(l,m)(d,e) * t3[a][b][c][m][j][i] * t1(e,k);
                                        }
                                    }
                                }
                            }

                            // 46: - P(ab/c) P(ij/k) <lm||de> t_abdijl t_ek t_cm
                            for (int d = n_Electrons; d < Matrix_Size; d++)
                            {
                                for (int l = 0; l < n_Electrons; l++)
                                {
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        for (int e = n_Electrons; e < Matrix_Size; e++)
                                        {
                                            t3_new[a][b][c][i][j][k] -= MOs(l,m)(d,e) * t3[a][b][d][i][j][l] * t1(e,k) * t1(c,m);
                                            t3_new[a][b][c][i][j][k] += MOs(l,m)(d,e) * t3[c][b][d][i][j][l] * t1(e,k) * t1(a,m);
                                            t3_new[a][b][c][i][j][k] += MOs(l,m)(d,e) * t3[a][c][d][i][j][l] * t1(e,k) * t1(b,m);

                                            t3_new[a][b][c][i][j][k] += MOs(l,m)(d,e) * t3[a][b][d][k][j][l] * t1(e,i) * t1(c,m);
                                            t3_new[a][b][c][i][j][k] += MOs(l,m)(d,e) * t3[a][b][d][i][k][l] * t1(e,j) * t1(c,m);

                                            t3_new[a][b][c][i][j][k] -= MOs(l,m)(d,e) * t3[c][b][d][k][j][l] * t1(e,i) * t1(a,m);
                                            t3_new[a][b][c][i][j][k] -= MOs(l,m)(d,e) * t3[c][b][d][i][k][l] * t1(e,j) * t1(a,m);

                                            t3_new[a][b][c][i][j][k] -= MOs(l,m)(d,e) * t3[a][c][d][k][j][l] * t1(e,i) * t1(b,m);
                                            t3_new[a][b][c][i][j][k] -= MOs(l,m)(d,e) * t3[a][c][d][i][k][l] * t1(e,j) * t1(b,m);
                                        }
                                    }
                                }
                            }

                            // 47: + 1/2 P(a/bc) P(ij/k) <lm||de> t_adij t_bclm t_ek
                            for (int d = n_Electrons; d < Matrix_Size; d++)
                            {
                                for (int l = 0; l < n_Electrons; l++)
                                {
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        for (int e = n_Electrons; e < Matrix_Size; e++)
                                        {
                                            t3_new[a][b][c][i][j][k] += 0.5 * MOs(l,m)(d,e) * t2(a,d)(i,j) * t2(b,c)(l,m) * t1(e,k);
                                            t3_new[a][b][c][i][j][k] -= 0.5 * MOs(l,m)(d,e) * t2(b,d)(i,j) * t2(a,c)(l,m) * t1(e,k);
                                            t3_new[a][b][c][i][j][k] -= 0.5 * MOs(l,m)(d,e) * t2(c,d)(i,j) * t2(b,a)(l,m) * t1(e,k);

                                            t3_new[a][b][c][i][j][k] -= 0.5 * MOs(l,m)(d,e) * t2(a,d)(k,j) * t2(b,c)(l,m) * t1(e,i);
                                            t3_new[a][b][c][i][j][k] -= 0.5 * MOs(l,m)(d,e) * t2(a,d)(i,k) * t2(b,c)(l,m) * t1(e,j);

                                            t3_new[a][b][c][i][j][k] += 0.5 * MOs(l,m)(d,e) * t2(b,d)(k,j) * t2(a,c)(l,m) * t1(e,i);
                                            t3_new[a][b][c][i][j][k] += 0.5 * MOs(l,m)(d,e) * t2(b,d)(i,k) * t2(a,c)(l,m) * t1(e,j);

                                            t3_new[a][b][c][i][j][k] += 0.5 * MOs(l,m)(d,e) * t2(c,d)(k,j) * t2(b,a)(l,m) * t1(e,i);
                                            t3_new[a][b][c][i][j][k] += 0.5 * MOs(l,m)(d,e) * t2(c,d)(i,k) * t2(b,a)(l,m) * t1(e,j);
                                        }
                                    }
                                }
                            }

                            // 48: + 1/2 P(a/bc) P(ij/k) <lm||de> t_deij t_bcmk t_al
                            for (int d = n_Electrons; d < Matrix_Size; d++)
                            {
                                for (int l = 0; l < n_Electrons; l++)
                                {
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        for (int e = n_Electrons; e < Matrix_Size; e++)
                                        {
                                            t3_new[a][b][c][i][j][k] += 0.5 * MOs(l,m)(d,e) * t2(d,e)(i,j) * t2(b,c)(m,k) * t1(a,l);
                                            t3_new[a][b][c][i][j][k] -= 0.5 * MOs(l,m)(d,e) * t2(d,e)(i,j) * t2(a,c)(m,k) * t1(b,l);
                                            t3_new[a][b][c][i][j][k] -= 0.5 * MOs(l,m)(d,e) * t2(d,e)(i,j) * t2(b,a)(m,k) * t1(c,l);

                                            t3_new[a][b][c][i][j][k] -= 0.5 * MOs(l,m)(d,e) * t2(d,e)(k,j) * t2(b,c)(m,i) * t1(a,l);
                                            t3_new[a][b][c][i][j][k] -= 0.5 * MOs(l,m)(d,e) * t2(d,e)(i,k) * t2(b,c)(m,j) * t1(a,l);

                                            t3_new[a][b][c][i][j][k] += 0.5 * MOs(l,m)(d,e) * t2(d,e)(k,j) * t2(a,c)(m,i) * t1(b,l);
                                            t3_new[a][b][c][i][j][k] += 0.5 * MOs(l,m)(d,e) * t2(d,e)(i,k) * t2(a,c)(m,j) * t1(b,l);

                                            t3_new[a][b][c][i][j][k] += 0.5 * MOs(l,m)(d,e) * t2(d,e)(k,j) * t2(b,a)(m,i) * t1(c,l);
                                            t3_new[a][b][c][i][j][k] += 0.5 * MOs(l,m)(d,e) * t2(d,e)(i,k) * t2(b,a)(m,j) * t1(c,l);
                                        }
                                    }
                                }
                            }

                            // 49: - P(abc) P(ij/k) <lm||de> t_adij t_ecmk t_bl
                            for (int d = n_Electrons; d < Matrix_Size; d++)
                            {
                                for (int l = 0; l < n_Electrons; l++)
                                {
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        for (int e = n_Electrons; e < Matrix_Size; e++)
                                        {
                                            t3_new[a][b][c][i][j][k] -= MOs(l,m)(d,e) * t2(a,d)(i,j) * t2(e,c)(m,k) * t1(b,l);
                                            t3_new[a][b][c][i][j][k] += MOs(l,m)(d,e) * t2(b,d)(i,j) * t2(e,c)(m,k) * t1(a,l);
                                            t3_new[a][b][c][i][j][k] += MOs(l,m)(d,e) * t2(c,d)(i,j) * t2(e,a)(m,k) * t1(b,l);
                                            t3_new[a][b][c][i][j][k] += MOs(l,m)(d,e) * t2(a,d)(i,j) * t2(e,b)(m,k) * t1(c,l);
                                            t3_new[a][b][c][i][j][k] -= MOs(l,m)(d,e) * t2(c,d)(i,j) * t2(e,b)(m,k) * t1(a,l);
                                            t3_new[a][b][c][i][j][k] -= MOs(l,m)(d,e) * t2(b,d)(i,j) * t2(e,a)(m,k) * t1(c,l);

                                            t3_new[a][b][c][i][j][k] += MOs(l,m)(d,e) * t2(a,d)(k,j) * t2(e,c)(m,i) * t1(b,l);
                                            t3_new[a][b][c][i][j][k] -= MOs(l,m)(d,e) * t2(b,d)(k,j) * t2(e,c)(m,i) * t1(a,l);
                                            t3_new[a][b][c][i][j][k] -= MOs(l,m)(d,e) * t2(c,d)(k,j) * t2(e,a)(m,i) * t1(b,l);
                                            t3_new[a][b][c][i][j][k] -= MOs(l,m)(d,e) * t2(a,d)(k,j) * t2(e,b)(m,i) * t1(c,l);
                                            t3_new[a][b][c][i][j][k] += MOs(l,m)(d,e) * t2(c,d)(k,j) * t2(e,b)(m,i) * t1(a,l);
                                            t3_new[a][b][c][i][j][k] += MOs(l,m)(d,e) * t2(b,d)(k,j) * t2(e,a)(m,i) * t1(c,l);

                                            t3_new[a][b][c][i][j][k] += MOs(l,m)(d,e) * t2(a,d)(i,k) * t2(e,c)(m,j) * t1(b,l);
                                            t3_new[a][b][c][i][j][k] -= MOs(l,m)(d,e) * t2(b,d)(i,k) * t2(e,c)(m,j) * t1(a,l);
                                            t3_new[a][b][c][i][j][k] -= MOs(l,m)(d,e) * t2(c,d)(i,k) * t2(e,a)(m,j) * t1(b,l);
                                            t3_new[a][b][c][i][j][k] -= MOs(l,m)(d,e) * t2(a,d)(i,k) * t2(e,b)(m,j) * t1(c,l);
                                            t3_new[a][b][c][i][j][k] += MOs(l,m)(d,e) * t2(c,d)(i,k) * t2(e,b)(m,j) * t1(a,l);
                                            t3_new[a][b][c][i][j][k] += MOs(l,m)(d,e) * t2(b,d)(i,k) * t2(e,a)(m,j) * t1(c,l);
                                        }
                                    }
                                }
                            }

                            // 50: - P(a/bc) P(ijk) <lm||de> t_adil t_bcmk t_ej
                            for (int d = n_Electrons; d < Matrix_Size; d++)
                            {
                                for (int l = 0; l < n_Electrons; l++)
                                {
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        for (int e = n_Electrons; e < Matrix_Size; e++)
                                        {
                                            t3_new[a][b][c][i][j][k] -= MOs(l,m)(d,e) * t2(a,d)(i,l) * t2(b,c)(m,k) * t1(e,j);
                                            t3_new[a][b][c][i][j][k] += MOs(l,m)(d,e) * t2(a,d)(k,l) * t2(b,c)(m,i) * t1(e,j);
                                            t3_new[a][b][c][i][j][k] += MOs(l,m)(d,e) * t2(a,d)(j,l) * t2(b,c)(m,k) * t1(e,i);
                                            t3_new[a][b][c][i][j][k] += MOs(l,m)(d,e) * t2(a,d)(i,l) * t2(b,c)(m,j) * t1(e,k);
                                            t3_new[a][b][c][i][j][k] -= MOs(l,m)(d,e) * t2(a,d)(k,l) * t2(b,c)(m,j) * t1(e,i);
                                            t3_new[a][b][c][i][j][k] -= MOs(l,m)(d,e) * t2(a,d)(j,l) * t2(b,c)(m,i) * t1(e,k);

                                            t3_new[a][b][c][i][j][k] += MOs(l,m)(d,e) * t2(b,d)(i,l) * t2(a,c)(m,k) * t1(e,j);
                                            t3_new[a][b][c][i][j][k] -= MOs(l,m)(d,e) * t2(b,d)(k,l) * t2(a,c)(m,i) * t1(e,j);
                                            t3_new[a][b][c][i][j][k] -= MOs(l,m)(d,e) * t2(b,d)(j,l) * t2(a,c)(m,k) * t1(e,i);
                                            t3_new[a][b][c][i][j][k] -= MOs(l,m)(d,e) * t2(b,d)(i,l) * t2(a,c)(m,j) * t1(e,k);
                                            t3_new[a][b][c][i][j][k] += MOs(l,m)(d,e) * t2(b,d)(k,l) * t2(a,c)(m,j) * t1(e,i);
                                            t3_new[a][b][c][i][j][k] += MOs(l,m)(d,e) * t2(b,d)(j,l) * t2(a,c)(m,i) * t1(e,k);

                                            t3_new[a][b][c][i][j][k] += MOs(l,m)(d,e) * t2(c,d)(i,l) * t2(b,a)(m,k) * t1(e,j);
                                            t3_new[a][b][c][i][j][k] -= MOs(l,m)(d,e) * t2(c,d)(k,l) * t2(b,a)(m,i) * t1(e,j);
                                            t3_new[a][b][c][i][j][k] -= MOs(l,m)(d,e) * t2(c,d)(j,l) * t2(b,a)(m,k) * t1(e,i);
                                            t3_new[a][b][c][i][j][k] -= MOs(l,m)(d,e) * t2(c,d)(i,l) * t2(b,a)(m,j) * t1(e,k);
                                            t3_new[a][b][c][i][j][k] += MOs(l,m)(d,e) * t2(c,d)(k,l) * t2(b,a)(m,j) * t1(e,i);
                                            t3_new[a][b][c][i][j][k] += MOs(l,m)(d,e) * t2(c,d)(j,l) * t2(b,a)(m,i) * t1(e,k);
                                        }
                                    }
                                }
                            }

                            // 51: - P(a/bc) P(ij/k) <lm||de> t_aeij t_bcmk t_dl
                            for (int d = n_Electrons; d < Matrix_Size; d++)
                            {
                                for (int l = 0; l < n_Electrons; l++)
                                {
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        for (int e = n_Electrons; e < Matrix_Size; e++)
                                        {
                                            t3_new[a][b][c][i][j][k] -= MOs(l,m)(d,e) * t2(a,e)(i,j) * t2(b,c)(m,k) * t1(d,l);
                                            t3_new[a][b][c][i][j][k] += MOs(l,m)(d,e) * t2(b,e)(i,j) * t2(a,c)(m,k) * t1(d,l);
                                            t3_new[a][b][c][i][j][k] += MOs(l,m)(d,e) * t2(c,e)(i,j) * t2(b,a)(m,k) * t1(d,l);

                                            t3_new[a][b][c][i][j][k] += MOs(l,m)(d,e) * t2(a,e)(k,j) * t2(b,c)(m,i) * t1(d,l);
                                            t3_new[a][b][c][i][j][k] += MOs(l,m)(d,e) * t2(a,e)(i,k) * t2(b,c)(m,j) * t1(d,l);

                                            t3_new[a][b][c][i][j][k] -= MOs(l,m)(d,e) * t2(b,e)(k,j) * t2(a,c)(m,i) * t1(d,l);
                                            t3_new[a][b][c][i][j][k] -= MOs(l,m)(d,e) * t2(b,e)(i,k) * t2(a,c)(m,j) * t1(d,l);

                                            t3_new[a][b][c][i][j][k] -= MOs(l,m)(d,e) * t2(c,e)(k,j) * t2(b,a)(m,i) * t1(d,l);
                                            t3_new[a][b][c][i][j][k] -= MOs(l,m)(d,e) * t2(c,e)(i,k) * t2(b,a)(m,j) * t1(d,l);
                                        }
                                    }
                                }
                            }

                            // 52: + P(a/bc) P(ij/k) <lm||de> t_bcmk t_di t_al t_ej
                            for (int d = n_Electrons; d < Matrix_Size; d++)
                            {
                                for (int l = 0; l < n_Electrons; l++)
                                {
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        for (int e = n_Electrons; e < Matrix_Size; e++)
                                        {
                                            t3_new[a][b][c][i][j][k] += MOs(l,m)(d,e) * t2(b,c)(m,k) * t1(d,i) * t1(a,l) * t1(e,j);
                                            t3_new[a][b][c][i][j][k] -= MOs(l,m)(d,e) * t2(a,c)(m,k) * t1(d,i) * t1(b,l) * t1(e,j);
                                            t3_new[a][b][c][i][j][k] -= MOs(l,m)(d,e) * t2(b,a)(m,k) * t1(d,i) * t1(c,l) * t1(e,j);

                                            t3_new[a][b][c][i][j][k] -= MOs(l,m)(d,e) * t2(b,c)(m,i) * t1(d,k) * t1(a,l) * t1(e,j);
                                            t3_new[a][b][c][i][j][k] -= MOs(l,m)(d,e) * t2(b,c)(m,j) * t1(d,i) * t1(a,l) * t1(e,k);

                                            t3_new[a][b][c][i][j][k] += MOs(l,m)(d,e) * t2(a,c)(m,i) * t1(d,k) * t1(b,l) * t1(e,j);
                                            t3_new[a][b][c][i][j][k] += MOs(l,m)(d,e) * t2(a,c)(m,j) * t1(d,i) * t1(b,l) * t1(e,k);

                                            t3_new[a][b][c][i][j][k] += MOs(l,m)(d,e) * t2(b,a)(m,i) * t1(d,k) * t1(c,l) * t1(e,j);
                                            t3_new[a][b][c][i][j][k] += MOs(l,m)(d,e) * t2(b,a)(m,j) * t1(d,i) * t1(c,l) * t1(e,k);
                                        }
                                    }
                                }
                            }

                            // 53: + P(a/bc) P(ij/k) <lm||de> t_adij t_bl t_ek t_cm
                            for (int d = n_Electrons; d < Matrix_Size; d++)
                            {
                                for (int l = 0; l < n_Electrons; l++)
                                {
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        for (int e = n_Electrons; e < Matrix_Size; e++)
                                        {
                                            t3_new[a][b][c][i][j][k] += MOs(l,m)(d,e) * t2(a,d)(i,j) * t1(b,l) * t1(e,k) * t1(c,m);
                                            t3_new[a][b][c][i][j][k] -= MOs(l,m)(d,e) * t2(b,d)(i,j) * t1(a,l) * t1(e,k) * t1(c,m);
                                            t3_new[a][b][c][i][j][k] -= MOs(l,m)(d,e) * t2(c,d)(i,j) * t1(b,l) * t1(e,k) * t1(a,m);

                                            t3_new[a][b][c][i][j][k] -= MOs(l,m)(d,e) * t2(a,d)(k,j) * t1(b,l) * t1(e,i) * t1(c,m);
                                            t3_new[a][b][c][i][j][k] -= MOs(l,m)(d,e) * t2(a,d)(i,k) * t1(b,l) * t1(e,j) * t1(c,m);

                                            t3_new[a][b][c][i][j][k] += MOs(l,m)(d,e) * t2(b,d)(k,j) * t1(a,l) * t1(e,i) * t1(c,m);
                                            t3_new[a][b][c][i][j][k] += MOs(l,m)(d,e) * t2(b,d)(i,k) * t1(a,l) * t1(e,j) * t1(c,m);

                                            t3_new[a][b][c][i][j][k] += MOs(l,m)(d,e) * t2(c,d)(k,j) * t1(b,l) * t1(e,i) * t1(a,m);
                                            t3_new[a][b][c][i][j][k] += MOs(l,m)(d,e) * t2(c,d)(i,k) * t1(b,l) * t1(e,j) * t1(a,m);
                                        }
                                    }
                                }
                            }

                            t3_new[a][b][c][i][j][k] = t3_new[a][b][c][i][j][k] / (fs(i,i) + fs(j,j) + fs(k,k) - fs(a,a) - fs(b,b) - fs(c,c));

                            // Symmetries
                            t3_new[b][a][c][i][j][k] = -t3_new[a][b][c][i][j][k];
                            t3_new[a][b][c][j][i][k] = -t3_new[a][b][c][i][j][k];
                            t3_new[b][a][c][j][i][k] = t3_new[a][b][c][i][j][k];
                            t3_new[c][b][a][i][j][k] = -t3_new[a][b][c][i][j][k];
                            t3_new[c][b][a][j][i][k] = t3_new[a][b][c][i][j][k];
                            t3_new[c][a][b][i][j][k] = t3_new[a][b][c][i][j][k];
                            t3_new[c][a][b][j][i][k] = -t3_new[a][b][c][i][j][k];
                            t3_new[a][c][b][i][j][k] = -t3_new[a][b][c][i][j][k];
                            t3_new[a][c][b][j][i][k] = t3_new[a][b][c][i][j][k];

                            // k symmetries here
                            t3_new[a][b][c][i][k][j] = -t3_new[a][b][c][i][j][k];
                            t3_new[a][b][c][k][j][i] = -t3_new[a][b][c][i][j][k];
                            t3_new[b][a][c][k][j][i] = -t3_new[b][a][c][i][j][k];
                            t3_new[b][a][c][i][k][j] = -t3_new[b][a][c][i][j][k];
                            t3_new[a][b][c][k][i][j] = -t3_new[a][b][c][j][i][k];
                            t3_new[a][b][c][j][k][i] = -t3_new[a][b][c][j][i][k];
                            t3_new[b][a][c][k][i][j] = -t3_new[b][a][c][j][i][k];
                            t3_new[b][a][c][j][k][i] = -t3_new[b][a][c][j][i][k];
                            t3_new[c][b][a][k][j][i] = -t3_new[c][b][a][i][j][k];
                            t3_new[c][b][a][i][k][j] = -t3_new[c][b][a][i][j][k];
                            t3_new[c][b][a][k][i][j] = -t3_new[c][b][a][j][i][k];
                            t3_new[c][b][a][j][k][i] = -t3_new[c][b][a][j][i][k];
                            t3_new[c][a][b][k][j][i] = -t3_new[c][a][b][i][j][k];
                            t3_new[c][a][b][i][k][j] = -t3_new[c][a][b][i][j][k];
                            t3_new[c][a][b][k][i][j] = -t3_new[c][a][b][j][i][k];
                            t3_new[c][a][b][j][k][i] = -t3_new[c][a][b][j][i][k];
                            t3_new[a][c][b][k][j][i] = -t3_new[a][c][b][i][j][k];
                            t3_new[a][c][b][i][k][j] = -t3_new[a][c][b][i][j][k];
                            t3_new[a][c][b][k][i][j] = -t3_new[a][c][b][j][i][k];
                            t3_new[a][c][b][j][k][i] = -t3_new[a][c][b][j][i][k];
                        }
                    }
                }
            }
        }
    }


}

void CC_General::Fill_t4_new()
{
    // Unfinished T4 amplitudes

    for (int a = n_Electrons; a < Matrix_Size; a++)
    {
        for (int b = a+1; b < Matrix_Size; b++)
        {
            for (int c = n_Electrons; c < Matrix_Size; c++) // b+1
            {
                for (int d = n_Electrons; d < Matrix_Size; d++) // c+1
                {
                    for (int i = 0; i < n_Electrons; i++)
                    {
                        for (int j = 0; j < n_Electrons; j++) // i+1
                        {
                            for (int k = 0; k < n_Electrons; k++) // j+1
                            {
                                for (int l = 0; l < n_Electrons; l++) // k+1
                                {
                                    // Reset, listed are the number of the equation in referance and pulled out the part needed for programmer
                                    t4_new[a][b][c][d][i][j][k][l] = 0;
                                    
                                    // 1 + P(l/ijk) P(ab/cd) <cd||el> t_abeijk
                                    for (int e = n_Electrons; e < Matrix_Size; e++)
                                    {
                                        t4_new[a][b][c][d][i][j][k][l] = 0;
                                    }
                                    
                                    // 2 - P(ij/kl) P(abc/d) <md||kl> t_abcijm
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        t4_new[a][b][c][d][i][j][k][l] = 0;
                                    }
                                    
                                    // 3 + P(abc/d) f_de t_abceijkl
                                    for (int e = n_Electrons; e < Matrix_Size; e++)
                                    {
                                        t4_new[a][b][c][d][i][j][k][l] += EqualFunc(d,e) * fs(d,e) * t4[a][b][c][e][i][j][k][l];
                                        t4_new[a][b][c][d][i][j][k][l] -= EqualFunc(a,e) * fs(a,e) * t4[d][b][c][e][i][j][k][l];
                                        t4_new[a][b][c][d][i][j][k][l] -= EqualFunc(b,e) * fs(b,e) * t4[a][d][c][e][i][j][k][l];
                                        t4_new[a][b][c][d][i][j][k][l] -= EqualFunc(c,e) * fs(c,e) * t4[a][b][d][e][i][j][k][l];
                                    }
                                    
                                    // 4 - P(ijk/l) f_ml t_abcdijkm
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        t4_new[a][b][c][d][i][j][k][l] += EqualFunc(m,l) * fs(m,l) * t4[a][b][c][d][i][j][k][m];
                                        t4_new[a][b][c][d][i][j][k][l] -= EqualFunc(m,i) * fs(m,i) * t4[a][b][c][d][l][j][k][m];
                                        t4_new[a][b][c][d][i][j][k][l] -= EqualFunc(m,j) * fs(m,j) * t4[a][b][c][d][i][l][k][m];
                                        t4_new[a][b][c][d][i][j][k][l] -= EqualFunc(m,k) * fs(m,k) * t4[a][b][c][d][i][j][l][m];
                                    }
                                    
                                    // 5 - P(ijk/l) P(abc/d) <md||le> t_abceijkm
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        for (int e = n_Electrons; e < Matrix_Size; e++)
                                        {
                                            t4_new[a][b][c][d][i][j][k][l] -= MOs(m,d)(l,e) * t4[a][b][c][e][i][j][k][m];
                                            t4_new[a][b][c][d][i][j][k][l] += MOs(m,d)(i,e) * t4[a][b][c][e][l][j][k][m];
                                            t4_new[a][b][c][d][i][j][k][l] += MOs(m,d)(j,e) * t4[a][b][c][e][i][l][k][m];
                                            t4_new[a][b][c][d][i][j][k][l] += MOs(m,d)(k,e) * t4[a][b][c][e][i][j][l][m];

                                            t4_new[a][b][c][d][i][j][k][l] += MOs(m,a)(l,e) * t4[d][b][c][e][i][j][k][m];
                                            t4_new[a][b][c][d][i][j][k][l] -= MOs(m,a)(i,e) * t4[d][b][c][e][l][j][k][m];
                                            t4_new[a][b][c][d][i][j][k][l] -= MOs(m,a)(j,e) * t4[d][b][c][e][i][l][k][m];
                                            t4_new[a][b][c][d][i][j][k][l] -= MOs(m,a)(k,e) * t4[d][b][c][e][i][j][l][m];

                                            t4_new[a][b][c][d][i][j][k][l] += MOs(m,b)(l,e) * t4[a][d][c][e][i][j][k][m];
                                            t4_new[a][b][c][d][i][j][k][l] -= MOs(m,b)(i,e) * t4[a][d][c][e][l][j][k][m];
                                            t4_new[a][b][c][d][i][j][k][l] -= MOs(m,b)(j,e) * t4[a][d][c][e][i][l][k][m];
                                            t4_new[a][b][c][d][i][j][k][l] -= MOs(m,b)(k,e) * t4[a][d][c][e][i][j][l][m];

                                            t4_new[a][b][c][d][i][j][k][l] += MOs(m,c)(l,e) * t4[a][b][d][e][i][j][k][m];
                                            t4_new[a][b][c][d][i][j][k][l] -= MOs(m,c)(i,e) * t4[a][b][d][e][l][j][k][m];
                                            t4_new[a][b][c][d][i][j][k][l] -= MOs(m,c)(j,e) * t4[a][b][d][e][i][l][k][m];
                                            t4_new[a][b][c][d][i][j][k][l] -= MOs(m,c)(k,e) * t4[a][b][d][e][i][j][l][m];
                                        }
                                    }
                                    
                                    // 6 + 1/2 P(ab/cd) <cd||ef> t_abefijkl
                                    for (int e = n_Electrons; e < Matrix_Size; e++)
                                    {
                                        for (int f = n_Electrons; f < Matrix_Size; f++)
                                        {
                                            t4_new[a][b][c][d][i][j][k][l] = 0;
                                        }

                                    }
                                    
                                    // 7 + 1/2 P(ij/kl) <mn||kl> t_abcdijmn
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        for (int n = 0; n < n_Electrons; n++)
                                        {
                                            t4_new[a][b][c][d][i][j][k][l] = 0;
                                        }

                                    }
                                    
                                    // 8 - P(ab/c/d) P(ijk/l) <cm||el> t_abeijk t_dm
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        for (int e = n_Electrons; e < Matrix_Size; e++)
                                        {
                                            t4_new[a][b][c][d][i][j][k][l] = 0;
                                        }
                                    }
                                    
                                    // 9 - P(abc/d) P(ij/k/l) <md||ke> t_abcijm t_el
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        for (int e = n_Electrons; e < Matrix_Size; e++)
                                        {
                                            t4_new[a][b][c][d][i][j][k][l] = 0;
                                        }
                                    }
                                    
                                    // 10 + P(ab/cd) P(ijk/l) <cd||ef> t_abeijk t_fl
                                    for (int e = n_Electrons; e < Matrix_Size; e++)
                                    {
                                        for (int f = n_Electrons; f < Matrix_Size; f++)
                                        {
                                            t4_new[a][b][c][d][i][j][k][l] = 0;
                                        }

                                    }
                                    
                                    // 11 + P(abc/d) P(ij/kl) <mn||kl> t_abcijm t_dn
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        for (int n = 0; n < n_Electrons; n++)
                                        {
                                            t4_new[a][b][c][d][i][j][k][l] = 0;
                                        }

                                    }
                                    
                                    // 12 - P(abc/d) f_me t_abceijkl t_dm
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        for (int e = n_Electrons; e < Matrix_Size; e++)
                                        {
                                            t4_new[a][b][c][d][i][j][k][l] -= fs(m,e) * t4[a][b][c][e][i][j][k][l] * t1(d,m);
                                            t4_new[a][b][c][d][i][j][k][l] += fs(m,e) * t4[d][b][c][e][i][j][k][l] * t1(a,m);
                                            t4_new[a][b][c][d][i][j][k][l] += fs(m,e) * t4[a][d][c][e][i][j][k][l] * t1(b,m);
                                            t4_new[a][b][c][d][i][j][k][l] += fs(m,e) * t4[a][b][d][e][i][j][k][l] * t1(c,m);
                                        }
                                    }
                                    
                                    // 13 - P(i/jkl) f_me t_abcdmjkl t_ei
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        for (int e = n_Electrons; e < Matrix_Size; e++)
                                        {
                                            t4_new[a][b][c][d][i][j][k][l] -= fs(m,e) * t4[a][b][c][d][m][j][k][l] * t1(e,i);
                                            t4_new[a][b][c][d][i][j][k][l] += fs(m,e) * t4[a][b][c][d][m][i][k][l] * t1(e,j);
                                            t4_new[a][b][c][d][i][j][k][l] += fs(m,e) * t4[a][b][c][d][m][j][i][l] * t1(e,k);
                                            t4_new[a][b][c][d][i][j][k][l] += fs(m,e) * t4[a][b][c][d][m][j][k][i] * t1(e,l);
                                        }
                                    }
                                    
                                    // 14 + P(abc/d) <md||ef> t_abceijkl t_fm
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        for (int e = n_Electrons; e < Matrix_Size; e++)
                                        {
                                            for (int f = n_Electrons; f < Matrix_Size; f++)
                                            {
                                                t4_new[a][b][c][d][i][j][k][l] += MOs(m,d)(e,f) * t4[a][b][c][e][i][j][k][l] * t1(f,m);
                                                t4_new[a][b][c][d][i][j][k][l] -= MOs(m,a)(e,f) * t4[d][b][c][e][i][j][k][l] * t1(f,m);
                                                t4_new[a][b][c][d][i][j][k][l] -= MOs(m,b)(e,f) * t4[a][d][c][e][i][j][k][l] * t1(f,m);
                                                t4_new[a][b][c][d][i][j][k][l] -= MOs(m,c)(e,f) * t4[a][b][d][e][i][j][k][l] * t1(f,m);
                                            }
                                        }
                                    }
                                    
                                    // 15 - 1/2 P(ab/c/d) <cm||ef> t_abefijkl t_dm
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        for (int e = n_Electrons; e < Matrix_Size; e++)
                                        {
                                            for (int f = n_Electrons; f < Matrix_Size; f++)
                                            {
                                                t4_new[a][b][c][d][i][j][k][l] = 0;
                                            }
                                        }
                                    }
                                    
                                    // 16 - P(abc/d) P(ijk/l) <mn||el> t_abceijkm t_dn
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        for (int e = n_Electrons; e < Matrix_Size; e++)
                                        {
                                            for (int n = 0; n < n_Electrons; n++)
                                            {
                                                t4_new[a][b][c][d][i][j][k][l] = 0;
                                            }
                                        }
                                    }
                                    
                                    // 17 - P(ijk/l) <mn||le> t_abcdijkm t_en
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        for (int e = n_Electrons; e < Matrix_Size; e++)
                                        {
                                            for (int n = 0; n < n_Electrons; n++)
                                            {
                                                t4_new[a][b][c][d][i][j][k][l] = 0;
                                            }
                                        }
                                    }
                                    
                                    // 18 + 1/2 P(i/j/kl) <mn||ej> t_abcdmnkl t_ei
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        for (int e = n_Electrons; e < Matrix_Size; e++)
                                        {
                                            for (int n = 0; n < n_Electrons; n++)
                                            {
                                                t4_new[a][b][c][d][i][j][k][l] = 0;
                                            }
                                        }
                                    }
                                    
                                    // 19 + P(abc/d) P(ijk/l) <md||ef> t_abceijkm t_fl
                                    for (int l = 0; l < n_Electrons; l++)
                                    {
                                        for (int e = n_Electrons; e < Matrix_Size; e++)
                                        {
                                            for (int f = n_Electrons; f < Matrix_Size; f++)
                                            {
                                                t4_new[a][b][c][d][i][j][k][l] = 0;
                                            }
                                        }
                                    }
                                    
                                    // 20 - P(a/b/cd) P(ij/k/l) <bm||ek> t_aeij t_cdml
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        for (int e = n_Electrons; e < Matrix_Size; e++)
                                        {
                                            t4_new[a][b][c][d][i][j][k][l] = 0;
                                        }
                                    }
                                    
                                    // 21 + P(bc/ad) P(ij/kl) <bc||ef> t_aeij t_fdkl
                                    for (int e = n_Electrons; e < Matrix_Size; e++)
                                    {
                                        for (int f = n_Electrons; f < Matrix_Size; f++)
                                        {
                                             t4_new[a][b][c][d][i][j][k][l] = 0;
                                        }
                                    }

                                    
                                    // 22 + P(ab/cd) P(il/jk) <mn||jk> t_abim t_cdnl
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        for (int n = 0; n < n_Electrons; n++)
                                        {
                                            t4_new[a][b][c][d][i][j][k][l] = 0;
                                        }

                                    }
                                    
                                    // 23 - P(ab/cd) P(ijk/l) f_me t_abeijk t_cdml
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        for (int e = n_Electrons; e < Matrix_Size; e++)
                                        {
                                            t4_new[a][b][c][d][i][j][k][l] = 0;
                                        }
                                    }
                                    
                                    // 24 - P(a/bcd) P(ij/kl) f_me t_bcdmkl t_aeij
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        for (int e = n_Electrons; e < Matrix_Size; e++)
                                        {
                                            t4_new[a][b][c][d][i][j][k][l] = 0;
                                        }
                                    }
                                    
                                    // 25 + P(ab/c/d) P(ijk/l) <cm||ef> t_abeijk t_fdml
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        for (int e = n_Electrons; e < Matrix_Size; e++)
                                        {
                                            for (int f = n_Electrons; f < Matrix_Size; f++)
                                            {
                                                t4_new[a][b][c][d][i][j][k][l] = 0;
                                            }
                                        }
                                    }
                                    
                                    // 26 + P(ab/c/d) P(ij/kl) <mc||ef> t_abeijm t_fdkl
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        for (int e = n_Electrons; e < Matrix_Size; e++)
                                        {
                                            for (int f = n_Electrons; f < Matrix_Size; f++)
                                            {
                                                t4_new[a][b][c][d][i][j][k][l] = 0;
                                            }
                                        }
                                    }
                                    
                                    // 27 - 1/2 P(a/bcd) P(ij/kl) <am||ef> t_bcdmkl t_efij
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        for (int e = n_Electrons; e < Matrix_Size; e++)
                                        {
                                            for (int f = n_Electrons; f < Matrix_Size; f++)
                                            {
                                                t4_new[a][b][c][d][i][j][k][l] = 0;
                                            }
                                        }
                                    }
                                    
                                    // 28 - 1/2 P(a/b/cd) P(ijk/l) <bm||ef> t_aefijk t_cdml
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        for (int e = n_Electrons; e < Matrix_Size; e++)
                                        {
                                            for (int f = n_Electrons; f < Matrix_Size; f++)
                                            {
                                                t4_new[a][b][c][d][i][j][k][l] = 0;
                                            }
                                        }
                                    }
                                    
                                    // 29 - P(abc/d) P(ij/k/l) <mn||ke> t_abcijm t_ednl
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        for (int e = n_Electrons; e < Matrix_Size; e++)
                                        {
                                            for (int n = 0; n < n_Electrons; n++)
                                            {
                                                t4_new[a][b][c][d][i][j][k][l] = 0;
                                            }
                                        }
                                    }
                                    
                                    // 30 - P(ab/cd) P(ij/k/l) <mn||ek> t_abeijm t_cdnl
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        for (int e = n_Electrons; e < Matrix_Size; e++)
                                        {
                                            for (int n = 0; n < n_Electrons; n++)
                                            {
                                                t4_new[a][b][c][d][i][j][k][l] = 0;
                                            }
                                        }
                                    }
                                    
                                    // 31 + 1/2 P(ab/cd) P(ijk/l) <mn||el> t_abeijk t_cdmn
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        for (int e = n_Electrons; e < Matrix_Size; e++)
                                        {
                                            for (int n = 0; n < n_Electrons; n++)
                                            {
                                                t4_new[a][b][c][d][i][j][k][l] = 0;
                                            }
                                        }
                                    }
                                    
                                    // 32 + 1/2 P(a/bcd) P(ij/k/l) <mn||ek> t_bcdmnl t_aeij
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        for (int e = n_Electrons; e < Matrix_Size; e++)
                                        {
                                            for (int n = 0; n < n_Electrons; n++)
                                            {
                                                t4_new[a][b][c][d][i][j][k][l] = 0;
                                            }
                                        }
                                    }
                                    
                                    // 33 + P(abc/d) P(ijk/l) <mn||ef> t_abceijkm t_fdnl
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        for (int e = n_Electrons; e < Matrix_Size; e++)
                                        {
                                            for (int n = 0; n < n_Electrons; n++)
                                            {
                                                for (int f = n_Electrons; f < Matrix_Size; f++)
                                                {
                                                    t4_new[a][b][c][d][i][j][k][l] = 0;
                                                }
                                            }
                                        }
                                    }
                                    
                                    // 34 + 1/4 P(ab/cd) <mn||ef> t_abefijkl t_cdmn
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        for (int e = n_Electrons; e < Matrix_Size; e++)
                                        {
                                            for (int n = 0; n < n_Electrons; n++)
                                            {
                                                for (int f = n_Electrons; f < Matrix_Size; f++)
                                                {
                                                    t4_new[a][b][c][d][i][j][k][l] = 0;
                                                }
                                            }
                                        }
                                    }
                                    
                                    // 35 + 1/4 P(ij/kl) <mn||ef> t_abcdmnkl t_efij
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        for (int e = n_Electrons; e < Matrix_Size; e++)
                                        {
                                            for (int n = 0; n < n_Electrons; n++)
                                            {
                                                for (int f = n_Electrons; f < Matrix_Size; f++)
                                                {
                                                    t4_new[a][b][c][d][i][j][k][l] = 0;
                                                }
                                            }
                                        }
                                    }
                                    
                                    // 36 - 1/2 P(ab/cd) P(ijk/l) <mn||ef> t_abefijmk t_cdnl
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        for (int e = n_Electrons; e < Matrix_Size; e++)
                                        {
                                            for (int n = 0; n < n_Electrons; n++)
                                            {
                                                for (int f = n_Electrons; f < Matrix_Size; f++)
                                                {
                                                    t4_new[a][b][c][d][i][j][k][l] = 0;
                                                }
                                            }
                                        }
                                    }
                                    
                                    // 37 - 1/2 P(a/bcd) P(ij/kl) <mn||ef> t_bfedmnkl t_aeij
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        for (int e = n_Electrons; e < Matrix_Size; e++)
                                        {
                                            for (int n = 0; n < n_Electrons; n++)
                                            {
                                                for (int f = n_Electrons; f < Matrix_Size; f++)
                                                {
                                                    t4_new[a][b][c][d][i][j][k][l] = 0;
                                                }
                                            }
                                        }
                                    }
                                    
                                    // 38 - 1/2 P(abc/d) <mn||ef> t_abceijkl t_dfmn
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        for (int e = n_Electrons; e < Matrix_Size; e++)
                                        {
                                            for (int n = 0; n < n_Electrons; n++)
                                            {
                                                for (int f = n_Electrons; f < Matrix_Size; f++)
                                                {
                                                    t4_new[a][b][c][d][i][j][k][l] -= 0.5 * MOs(m,n)(e,f) * t4[a][b][c][e][i][j][k][l] * t2(d,f)(m,n);
                                                    t4_new[a][b][c][d][i][j][k][l] += 0.5 * MOs(m,n)(e,f) * t4[d][b][c][e][i][j][k][l] * t2(a,f)(m,n);
                                                    t4_new[a][b][c][d][i][j][k][l] += 0.5 * MOs(m,n)(e,f) * t4[a][d][c][e][i][j][k][l] * t2(b,f)(m,n);
                                                    t4_new[a][b][c][d][i][j][k][l] += 0.5 * MOs(m,n)(e,f) * t4[a][b][d][e][i][j][k][l] * t2(c,f)(m,n);
                                                }
                                            }
                                        }
                                    }
                                    
                                    // 39 - 1/2 P(i/jkl) <mn||ef> t_abcdnjkl t_efmi
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        for (int e = n_Electrons; e < Matrix_Size; e++)
                                        {
                                            for (int n = 0; n < n_Electrons; n++)
                                            {
                                                for (int f = n_Electrons; f < Matrix_Size; f++)
                                                {
                                                    t4_new[a][b][c][d][i][j][k][l] -= 0.5 * MOs(m,n)(e,f) * t4[a][b][c][d][n][j][k][l] * t2(e,f)(m,i);
                                                    t4_new[a][b][c][d][i][j][k][l] += 0.5 * MOs(m,n)(e,f) * t4[a][b][c][d][n][i][k][l] * t2(e,f)(m,j);
                                                    t4_new[a][b][c][d][i][j][k][l] += 0.5 * MOs(m,n)(e,f) * t4[a][b][c][d][n][j][i][l] * t2(e,f)(m,k);
                                                    t4_new[a][b][c][d][i][j][k][l] += 0.5 * MOs(m,n)(e,f) * t4[a][b][c][d][n][j][k][i] * t2(e,f)(m,l);
                                                }
                                            }
                                        }
                                    }
                                    
                                    // 40 + P(a/cd) P(ij/kl) <mn||ef> t_abeijm t_fcdnkl
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        for (int e = n_Electrons; e < Matrix_Size; e++)
                                        {
                                            for (int n = 0; n < n_Electrons; n++)
                                            {
                                                for (int f = n_Electrons; f < Matrix_Size; f++)
                                                {
                                                    t4_new[a][b][c][d][i][j][k][l] = 0;
                                                }
                                            }
                                        }
                                    }
                                    
                                    // 41 + 1/4 P(a/bcd) P(ijk/l) <mn||ef> t_aefijk t_bcdmnl
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        for (int e = n_Electrons; e < Matrix_Size; e++)
                                        {
                                            for (int n = 0; n < n_Electrons; n++)
                                            {
                                                for (int f = n_Electrons; f < Matrix_Size; f++)
                                                {
                                                    t4_new[a][b][c][d][i][j][k][l] = 0;
                                                }
                                            }
                                        }
                                    }
                                    
                                    // 42 - 1/2 P(a/bcd) P(ij/kl) <mn||ef> t_aefimj t_bcdnkl
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        for (int e = n_Electrons; e < Matrix_Size; e++)
                                        {
                                            for (int n = 0; n < n_Electrons; n++)
                                            {
                                                for (int f = n_Electrons; f < Matrix_Size; f++)
                                                {
                                                    t4_new[a][b][c][d][i][j][k][l] = 0;
                                                }
                                            }
                                        }
                                    }
                                    
                                    // 43 - 1/2 P(ab/cd) P(ijk/l) <mn||ef> t_abeijk t_cfdmnl
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        for (int e = n_Electrons; e < Matrix_Size; e++)
                                        {
                                            for (int n = 0; n < n_Electrons; n++)
                                            {
                                                for (int f = n_Electrons; f < Matrix_Size; f++)
                                                {
                                                    t4_new[a][b][c][d][i][j][k][l] = 0;
                                                }
                                            }
                                        }
                                    }
                                    
                                    // 44 - P(a/bcd) P(ij/kl) <am||ef> t_bcdmkl t_ei t_fj
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        for (int e = n_Electrons; e < Matrix_Size; e++)
                                        {
                                            for (int f = n_Electrons; f < Matrix_Size; f++)
                                            {
                                                t4_new[a][b][c][d][i][j][k][l] = 0;
                                            }
                                        }
                                    }
                                    
                                    // 45 - P(ab/c/d) P(ijk/l) <cm||ef> t_abeijk t_fl t_dm
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        for (int e = n_Electrons; e < Matrix_Size; e++)
                                        {
                                            for (int f = n_Electrons; f < Matrix_Size; f++)
                                            {
                                                t4_new[a][b][c][d][i][j][k][l] = 0;
                                            }
                                        }
                                    }
                                    
                                    // 46 + P(ab/cd) P(ijk/l) <mn||el t_abeijk t_cm t_dn
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        for (int e = n_Electrons; e < Matrix_Size; e++)
                                        {
                                            for (int n = 0; n < n_Electrons; n++)
                                            {
                                                t4_new[a][b][c][d][i][j][k][l] = 0;
                                            }
                                        }
                                    }
                                    
                                    // 47 + P(abc/d) P(ij/k/l) <mn||ke> t_abcijm t_el t_dn
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        for (int e = n_Electrons; e < Matrix_Size; e++)
                                        {
                                            for (int n = 0; n < n_Electrons; n++)
                                            {
                                                t4_new[a][b][c][d][i][j][k][l] = 0;
                                            }
                                        }
                                    }
                                    
                                    // 48 - P(a/bcd) P(i/jkl) <mn||ef> t_fbcdnjkl t_ei t_am
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        for (int e = n_Electrons; e < Matrix_Size; e++)
                                        {
                                            for (int n = 0; n < n_Electrons; n++)
                                            {
                                                for (int f = n_Electrons; f < Matrix_Size; f++)
                                                {
                                                    t4_new[a][b][c][d][i][j][k][l] = 0;
                                                }
                                            }
                                        }
                                    }
                                    
                                    // 49 - P(abc/d) <mn||ef> t_abceijkl t_dm t_fn
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        for (int e = n_Electrons; e < Matrix_Size; e++)
                                        {
                                            for (int n = 0; n < n_Electrons; n++)
                                            {
                                                for (int f = n_Electrons; f < Matrix_Size; f++)
                                                {
                                                    t4_new[a][b][c][d][i][j][k][l] -= MOs(m,n)(e,f) * t4[a][b][c][e][i][j][k][l] * t1(d,m) * t1(f,n);
                                                    t4_new[a][b][c][d][i][j][k][l] += MOs(m,n)(e,f) * t4[d][b][c][e][i][j][k][l] * t1(a,m) * t1(f,n);
                                                    t4_new[a][b][c][d][i][j][k][l] += MOs(m,n)(e,f) * t4[a][d][c][e][i][j][k][l] * t1(b,m) * t1(f,n);
                                                    t4_new[a][b][c][d][i][j][k][l] += MOs(m,n)(e,f) * t4[a][b][d][e][i][j][k][l] * t1(c,m) * t1(f,n);
                                                }
                                            }
                                        }
                                    }
                                    
                                    // 50 - P(i/jkl) <mn||ef> t_abcdnjkl t_em t_fi
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        for (int e = n_Electrons; e < Matrix_Size; e++)
                                        {
                                            for (int n = 0; n < n_Electrons; n++)
                                            {
                                                for (int f = n_Electrons; f < Matrix_Size; f++)
                                                {
                                                    t4_new[a][b][c][d][i][j][k][l] -= MOs(m,n)(e,f) * t4[a][b][c][d][n][j][k][l] * t1(e,m) * t1(f,i);
                                                    t4_new[a][b][c][d][i][j][k][l] += MOs(m,n)(e,f) * t4[a][b][c][d][n][i][k][l] * t1(e,m) * t1(f,j);
                                                    t4_new[a][b][c][d][i][j][k][l] += MOs(m,n)(e,f) * t4[a][b][c][d][n][j][i][l] * t1(e,m) * t1(f,k);
                                                    t4_new[a][b][c][d][i][j][k][l] += MOs(m,n)(e,f) * t4[a][b][c][d][n][j][k][i] * t1(e,m) * t1(f,l);
                                                }
                                            }
                                        }
                                    }
                                    
                                    // 51 + 1/2 P(ij/kl) <mn||ef> t_abcdmnkl t_ei t_fj
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        for (int e = n_Electrons; e < Matrix_Size; e++)
                                        {
                                            for (int n = 0; n < n_Electrons; n++)
                                            {
                                                for (int f = n_Electrons; f < Matrix_Size; f++)
                                                {
                                                    t4_new[a][b][c][d][i][j][k][l] = 0;
                                                }
                                            }
                                        }
                                    }
                                    
                                    // 52 + 1/2 P(ab/cd) <mn||ef> t_abefijkl t_cm t_dn
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        for (int e = n_Electrons; e < Matrix_Size; e++)
                                        {
                                            for (int n = 0; n < n_Electrons; n++)
                                            {
                                                for (int f = n_Electrons; f < Matrix_Size; f++)
                                                {
                                                    t4_new[a][b][c][d][i][j][k][l] = 0;
                                                }
                                            }
                                        }
                                    }
                                    
                                    // 53 - P(ac/b/d) P(ij/kl) <bm||ef> t_aeij t_cfkl t_dm
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        for (int e = n_Electrons; e < Matrix_Size; e++)
                                        {
                                            for (int f = n_Electrons; f < Matrix_Size; f++)
                                            {
                                                t4_new[a][b][c][d][i][j][k][l] = 0;
                                            }
                                        }
                                    }
                                    
                                    // 54 + P(ab/c/d) P(i/j/kl) <mn||je> t_abim t_cekl t_dn
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        for (int e = n_Electrons; e < Matrix_Size; e++)
                                        {
                                            for (int n = 0; n < n_Electrons; n++)
                                            {
                                                t4_new[a][b][c][d][i][j][k][l] = 0;
                                            }
                                        }
                                    }
                                    
                                    // 55 + P(ab/cd) P(i/k/jl) <mn||ek> t_abmj t_cdnl t_ei
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        for (int e = n_Electrons; e < Matrix_Size; e++)
                                        {
                                            for (int n = 0; n < n_Electrons; n++)
                                            {
                                                t4_new[a][b][c][d][i][j][k][l] = 0;
                                            }
                                        }
                                    }
                                    
                                    // 56 - P(ab/c/d) P(i/j/kl) <mc||ef> t_abmj t_fdkl t_ei
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        for (int e = n_Electrons; e < Matrix_Size; e++)
                                        {
                                            for (int f = n_Electrons; f < Matrix_Size; f++)
                                            {
                                                t4_new[a][b][c][d][i][j][k][l] = 0;
                                            }
                                        }
                                    }
                                    
                                    // 57 - P(ab/cd) P(ijk/l) <mn||ef> t_abfijk t_cdnl t_em
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        for (int e = n_Electrons; e < Matrix_Size; e++)
                                        {
                                            for (int n = 0; n < n_Electrons; n++)
                                            {
                                                for (int f = n_Electrons; f < Matrix_Size; f++)
                                                {
                                                    t4_new[a][b][c][d][i][j][k][l] = 0;
                                                }
                                            }
                                        }
                                    }
                                    
                                    // 58 - P(a/bcd) P(ij/kl) <mn||ef> t_bcdnkl t_afij t_em
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        for (int e = n_Electrons; e < Matrix_Size; e++)
                                        {
                                            for (int n = 0; n < n_Electrons; n++)
                                            {
                                                for (int f = n_Electrons; f < Matrix_Size; f++)
                                                {
                                                    t4_new[a][b][c][d][i][j][k][l] = 0;
                                                }
                                            }
                                        }
                                    }
                                    
                                    // 59 - P(ab/c/d) P(ijk/l) <mn||ef> t_abeijk t_fdnl t_cm
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        for (int e = n_Electrons; e < Matrix_Size; e++)
                                        {
                                            for (int n = 0; n < n_Electrons; n++)
                                            {
                                                for (int f = n_Electrons; f < Matrix_Size; f++)
                                                {
                                                    t4_new[a][b][c][d][i][j][k][l] = 0;
                                                }
                                            }
                                        }
                                    }
                                    
                                    // 60 - P(a/bcd) P(i/j/kl) <mn||ef> t_bcdnkl t_aeim t_fj
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        for (int e = n_Electrons; e < Matrix_Size; e++)
                                        {
                                            for (int n = 0; n < n_Electrons; n++)
                                            {
                                                for (int f = n_Electrons; f < Matrix_Size; f++)
                                                {
                                                    t4_new[a][b][c][d][i][j][k][l] = 0;
                                                }
                                            }
                                        }
                                    }
                                    
                                    // 61 + 1/2 P(a/b/cd) P(ijk/l) <mn||ef> t_aefijk t_cdnl t_bm
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        for (int e = n_Electrons; e < Matrix_Size; e++)
                                        {
                                            for (int n = 0; n < n_Electrons; n++)
                                            {
                                                for (int f = n_Electrons; f < Matrix_Size; f++)
                                                {
                                                    t4_new[a][b][c][d][i][j][k][l] = 0;
                                                }
                                            }
                                        }
                                    }
                                    
                                    // 62 + 1/2 P(ab/cd) P(ijk/l) <mn||ef> t_abeijk t_cdmn t_fl
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        for (int e = n_Electrons; e < Matrix_Size; e++)
                                        {
                                            for (int n = 0; n < n_Electrons; n++)
                                            {
                                                for (int f = n_Electrons; f < Matrix_Size; f++)
                                                {
                                                    t4_new[a][b][c][d][i][j][k][l] = 0;
                                                }
                                            }
                                        }
                                    }
                                    
                                    // 63 + 1/2 P(a/bcd) P(ij/kl) <mn||ef> t_bcdnkl t_efij t_am
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        for (int e = n_Electrons; e < Matrix_Size; e++)
                                        {
                                            for (int n = 0; n < n_Electrons; n++)
                                            {
                                                for (int f = n_Electrons; f < Matrix_Size; f++)
                                                {
                                                    t4_new[a][b][c][d][i][j][k][l] = 0;
                                                }
                                            }
                                        }
                                    }
                                    
                                    // 64 - P(a/b/cd) P(ij/kl) <mn||ef> t_fcdnkl t_aeij t_bm
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        for (int e = n_Electrons; e < Matrix_Size; e++)
                                        {
                                            for (int n = 0; n < n_Electrons; n++)
                                            {
                                                for (int f = n_Electrons; f < Matrix_Size; f++)
                                                {
                                                    t4_new[a][b][c][d][i][j][k][l] = 0;
                                                }
                                            }
                                        }
                                    }
                                    
                                    // 65 - P(ab/cd) P(i/j/kl) <mn||ef> t_fcdnkl t_abmj t_ei
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        for (int e = n_Electrons; e < Matrix_Size; e++)
                                        {
                                            for (int n = 0; n < n_Electrons; n++)
                                            {
                                                for (int f = n_Electrons; f < Matrix_Size; f++)
                                                {
                                                    t4_new[a][b][c][d][i][j][k][l] = 0;
                                                }
                                            }
                                        }
                                    }
                                    
                                    // 66 + 1/2 P(a/bcd) P(ij/k/l) <mn||ef> t_bcdmnl t_aeij t_fk
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        for (int e = n_Electrons; e < Matrix_Size; e++)
                                        {
                                            for (int n = 0; n < n_Electrons; n++)
                                            {
                                                for (int f = n_Electrons; f < Matrix_Size; f++)
                                                {
                                                    t4_new[a][b][c][d][i][j][k][l] = 0;
                                                }
                                            }
                                        }
                                    }
                                    
                                    // 67 - P(a/b/cd) P(i/jk/l) <mn||ef> t_aeim t_bfjk t_cdnl
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        for (int e = n_Electrons; e < Matrix_Size; e++)
                                        {
                                            for (int n = 0; n < n_Electrons; n++)
                                            {
                                                for (int f = n_Electrons; f < Matrix_Size; f++)
                                                {
                                                    t4_new[a][b][c][d][i][j][k][l] = 0;
                                                }
                                            }
                                        }
                                    }
                                    
                                    // 68 + 1/2 P(ab/cd) P(ij/kl) <mn||ef> t_aeij t_bfkl t_cdmn
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        for (int e = n_Electrons; e < Matrix_Size; e++)
                                        {
                                            for (int n = 0; n < n_Electrons; n++)
                                            {
                                                for (int f = n_Electrons; f < Matrix_Size; f++)
                                                {
                                                    t4_new[a][b][c][d][i][j][k][l] = 0;
                                                }
                                            }
                                        }
                                    }
                                    
                                    // 69 + 1/2 P(ab/cd) P(ij/kl) <mn||ef> t_efij t_abmk t_cdnl
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        for (int e = n_Electrons; e < Matrix_Size; e++)
                                        {
                                            for (int n = 0; n < n_Electrons; n++)
                                            {
                                                for (int f = n_Electrons; f < Matrix_Size; f++)
                                                {
                                                    t4_new[a][b][c][d][i][j][k][l] = 0;
                                                }
                                            }
                                        }
                                    }
                                    
                                    // 70 + P(a/bcd) P(ij/kl) <mn||ef> t_bcdnkl t_ei t_am t_fj
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        for (int e = n_Electrons; e < Matrix_Size; e++)
                                        {
                                            for (int n = 0; n < n_Electrons; n++)
                                            {
                                                for (int f = n_Electrons; f < Matrix_Size; f++)
                                                {
                                                    t4_new[a][b][c][d][i][j][k][l] = 0;
                                                }
                                            }
                                        }
                                    }
                                    
                                    // 71 + P(ab/cd) P(ijk/l) <mn||ef> t_abeijk t_cm t_fl t_dn
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        for (int e = n_Electrons; e < Matrix_Size; e++)
                                        {
                                            for (int n = 0; n < n_Electrons; n++)
                                            {
                                                for (int f = n_Electrons; f < Matrix_Size; f++)
                                                {
                                                    t4_new[a][b][c][d][i][j][k][l] = 0;
                                                }
                                            }
                                        }
                                    }
                                    
                                    // 72 + P(a/b/cd) P(i/jk/l) <mn||ef> t_bfjk t_cdnl t_ei t_am
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        for (int e = n_Electrons; e < Matrix_Size; e++)
                                        {
                                            for (int n = 0; n < n_Electrons; n++)
                                            {
                                                for (int f = n_Electrons; f < Matrix_Size; f++)
                                                {
                                                    t4_new[a][b][c][d][i][j][k][l] = 0;
                                                }
                                            }
                                        }
                                    }
                                    
                                    // 73 + P(ab/cd) P(ik/jl) <mn||ef> t_abmj t_cdnl t_ei t_fk
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        for (int e = n_Electrons; e < Matrix_Size; e++)
                                        {
                                            for (int n = 0; n < n_Electrons; n++)
                                            {
                                                for (int f = n_Electrons; f < Matrix_Size; f++)
                                                {
                                                    t4_new[a][b][c][d][i][j][k][l] = 0;
                                                }
                                            }
                                        }
                                    }
                                    
                                    // 74 + P(ac/bd) P(ij/kl) <mn||ef> t_aeij t_cfkl t_bm t_dn
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        for (int e = n_Electrons; e < Matrix_Size; e++)
                                        {
                                            for (int n = 0; n < n_Electrons; n++)
                                            {
                                                for (int f = n_Electrons; f < Matrix_Size; f++)
                                                {
                                                    t4_new[a][b][c][d][i][j][k][l] = 0;
                                                }
                                            }
                                        }
                                    }

                                    // Denominator
                                    t4_new[a][b][c][d][i][j][k][l] = t4_new[a][b][c][d][i][j][k][l] / (fs(i,i) + fs(j,j) + fs(k,k) + fs(l,l) - fs(a,a) - fs(b,b) - fs(c,c) - fs(d,d));

                                    // Symmetries (very important here), long list comming
                                    t4_new[b][a][c][d][i][j][k][l] = t4_new[a][b][c][d][i][j][k][l];
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void CC_General::construct_fs(vec f_vec)
{
    fs = zeros(Matrix_Size, Matrix_Size);
    for (int i = 0; i < Matrix_Size; i++)
    {
        fs(i,i) = f_vec(i/2);
    }
}

void CC_General::construct_integrals(vec integrals)
{
    int matsize = Matrix_Size/2;
    vector<cube> temp;
    vector <cube> temp2;
    vector <cube> temp3;
    vector <cube> temp4;
    cube temp_cube = zeros(matsize, matsize, matsize);
    mat temp_matrix = zeros(Matrix_Size, Matrix_Size);
    MOs.set_size(Matrix_Size, Matrix_Size);

    for (int i = 0; i < matsize; i++)
    {
        temp.push_back(temp_cube);
        temp2.push_back(temp_cube);
        temp3.push_back(temp_cube);
        temp4.push_back(temp_cube);
    }

    for (int i = 0; i < Matrix_Size; i++)
    {
        for (int j = 0; j < Matrix_Size; j++)
        {
            MOs(i,j) = temp_matrix;
        }
    }


    for (int a = 0; a < matsize; a++)
    {
        for (int i = 0; i < matsize; i++)
        {
            for (int j = 0; j < matsize; j++)
            {
                for (int k = 0; k < matsize; k++)
                {
                    for (int l = 0; l < matsize; l++)
                    {
                        temp.at(a)(j,k,l) += integrals(Return_Integral_index(i,j,k,l)) * c(i,a);
                    }
                }
            }
        }
    }

    for (int a = 0; a < matsize; a++)
    {
        for (int b = 0; b < matsize; b++)
        {
            for (int i = 0; i < matsize; i++)
            {
                for (int j = 0; j < matsize; j++)
                {
                    for (int k = 0; k < matsize; k++)
                    {
                        temp2.at(a)(b,j,k) += c(i,b) * temp.at(a)(i,j,k);
                    }
                }
            }
        }
    }
    temp.clear();

    for (int a = 0; a < matsize; a++)
    {
        for (int b = 0; b < matsize; b++)
        {
            for (int g = 0; g < matsize; g++)
            {
                for (int i = 0; i < matsize; i++)
                {
                    for (int j = 0; j < matsize; j++)
                    {
                        temp3.at(a)(b,g,j) += c(i,g) * temp2.at(a)(b,i,j);
                    }
                }
            }
        }
    }

    temp2.clear();

    for (int a = 0; a < matsize; a++)
    {
        for (int b = 0; b < matsize; b++)
        {
            for (int g = 0; g < matsize; g++)
            {
                for (int h = 0; h < matsize; h++)
                {
                    for (int i = 0; i < matsize; i++)
                    {
                        temp4.at(a)(b,g,h) += c(i,h) * temp3.at(a)(b,g,i);
                    }
                }
            }
        }
    }

    temp3.clear();

    double val1,val2;

    for (int a = 0; a < Matrix_Size; a++)
    {
        for (int b = 0; b < Matrix_Size; b++)
        {
            for (int i = 0; i < Matrix_Size; i++)
            {
                for (int j = 0; j < Matrix_Size; j++)
                {
                    val1 = EqualFunc(a%2, b%2) * EqualFunc(i%2, j%2) * temp4.at(a/2)(b/2,i/2,j/2);
                    val2 = EqualFunc(a%2, j%2) * EqualFunc(i%2, b%2) * temp4.at(a/2)(j/2,i/2,b/2);
                    MOs(a,b)(i,j) = val1-val2;
                }
            }
        }
    }
}

int CC_General::Return_Integral_index(int a, int b, int i, int j)
{
    long int ab, ij;

    if (a > b)
    {
        ab = (a*(a+1))/2 + b;
    }

    else
    {
        ab = (b*(b+1))/2 + a;
    }

    if (i > j)
    {
        ij = (i*(i+1))/2 + j;
    }

    else
    {
        ij = (j*(j+1))/2 + i;
    }

    if (ab > ij)
    {
        return ((ab*(ab+1)/2) + ij);
    }

    else
    {
        return ((ij*(ij+1)/2) + ab);
    }
}

int CC_General::EqualFunc(int a, int b)
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
