#include "ccsdt.h"

CCSDT::CCSDT(int n_N, vec zz, mat rr, string B_s, int n_Elec, int ran, int siz, Hartree_Fock_Solver *Hartfock, bool frez)
{
    n_Nuclei = n_N;
    Z = zz;
    R = rr;
    Basis_Set = B_s;
    n_Electrons = n_Elec;
    rank = ran;
    size = siz;
    HartFock = Hartfock;
    freeze_core = frez;
}

/*

  Answers to get on benchmark R-CCSDT:

  R = R
- 0.147 594 (2)

  R = 1.5R
- 0.209 519 (470)

            CCSDT program

*/


double CCSDT::Calc_Energy(double toler, int numbah)
{
    // Initiate system
    Method_Nr = numbah;
    double Energy;
    bool continue_ccsdt = true;
    iter = 1;
    int IterMax = 1000;
    double E_old, E_new, E_HF;
    double temp_en;
    int half_elec = (int) n_Electrons/2;

    // Do some Hartree Fock for referance
    if (Method_Nr == 7)
    {
        HartFock->Set_UP_DOWN(half_elec+1, half_elec);
        E_HF = HartFock->get_Energy(toler, 1);
    }

    else
    {
        E_HF = HartFock->get_Energy(toler, 0);
    }

    Matrix_Size = 2*HartFock->ReturnMatrixSize();
    unocc_orb = Matrix_Size - n_Electrons;

    // Put in some MOs and Fock eigenvalues
    AOtoMO();

    // Allocate t3
    t3 = (double******)malloc(Matrix_Size * sizeof(double*****));
    t3_new = (double******)malloc(Matrix_Size * sizeof(double*****));
    for (int i = n_Electrons; i < Matrix_Size; i++)
    {
        t3[i] = (double*****)malloc(Matrix_Size * sizeof(double****));
        t3_new[i] = (double*****)malloc(Matrix_Size * sizeof(double****));
        for (int j = n_Electrons; j < Matrix_Size; j++)
        {
            t3[i][j] = (double****)malloc(Matrix_Size*sizeof(double***));
            t3_new[i][j] = (double****)malloc(Matrix_Size*sizeof(double***));

            for (int k = n_Electrons; k < Matrix_Size; k++)
            {
                t3[i][j][k] = (double***)malloc(n_Electrons*sizeof(double**));
                t3_new[i][j][k] = (double***)malloc(n_Electrons*sizeof(double**));

                for (int l = 0; l < n_Electrons; l++)
                {
                    t3[i][j][k][l] = (double**)malloc(n_Electrons*sizeof(double*));
                    t3_new[i][j][k][l] = (double**)malloc(n_Electrons*sizeof(double*));

                    for (int n = 0; n < n_Electrons; n++)
                    {
                        t3[i][j][k][l][n] = (double*)malloc(n_Electrons*sizeof(double));
                        t3_new[i][j][k][l][n] = (double*)malloc(n_Electrons*sizeof(double));
                        for (int zz = 0; zz < n_Electrons; zz++)
                        {
                            t3[i][j][k][l][n][zz] = 0;
                            t3_new[i][j][k][l][n][zz] = 0;
                        }
                    }
                }
            }
        }
    }

    // Allocate t1
    t1 = zeros(Matrix_Size, Matrix_Size);
    t1_new = zeros(Matrix_Size, Matrix_Size);

    // Allocate t2
    t2.set_size(Matrix_Size, Matrix_Size);
    t2_new.set_size(Matrix_Size, Matrix_Size);
    for (int i = n_Electrons; i < Matrix_Size; i++)
    {
        for (int j = n_Electrons; j < Matrix_Size; j++)
        {
            t2(i,j) = zeros(n_Electrons, n_Electrons);
            t2_new(i,j) = zeros(n_Electrons, n_Electrons);
        }
    }

    // Allocate intermediates
    X1.set_size(Matrix_Size, Matrix_Size);
    X2.set_size(Matrix_Size, n_Electrons);
    X3.set_size(n_Electrons, n_Electrons);
    X4.set_size(Matrix_Size, Matrix_Size);
    X6.set_size(Matrix_Size, n_Electrons);
    X7 = zeros(n_Electrons, n_Electrons);
    X8 = zeros(Matrix_Size, Matrix_Size);
    X12.set_size(Matrix_Size, Matrix_Size);
    X13.set_size(Matrix_Size, n_Electrons);
    X14.set_size(Matrix_Size, Matrix_Size);
    X15.set_size(Matrix_Size, n_Electrons);

    tau3.set_size(Matrix_Size, Matrix_Size);
    den_ai = zeros(Matrix_Size, Matrix_Size);
    W1.set_size(n_Electrons, n_Electrons);
    W2.set_size(Matrix_Size, n_Electrons);
    W3.set_size(Matrix_Size, n_Electrons);
    W4.set_size(Matrix_Size, n_Electrons);
    D1 = zeros(Matrix_Size, Matrix_Size);
    D2 = zeros(Matrix_Size, Matrix_Size);
    D3 = zeros(Matrix_Size, Matrix_Size);

    // Allocate 4D intermediates
    for (int i = 0; i < Matrix_Size; i++)
    {
        for (int j = 0; j < Matrix_Size; j++)
        {
            if ((i >= n_Electrons) && (j >= n_Electrons))
            {
                X1(i,j) = zeros(n_Electrons, Matrix_Size);
                X12(i,j) = zeros(n_Electrons, Matrix_Size);
                X14(i,j) = zeros(n_Electrons, Matrix_Size);
                X4(i,j) = zeros(Matrix_Size, Matrix_Size);

                tau3(i,j) = zeros(n_Electrons, n_Electrons);
            }

            if ((i >= n_Electrons) && (j < n_Electrons))
            {
                X2(i,j) = zeros(n_Electrons, n_Electrons);
                X13(i,j) = zeros(n_Electrons, n_Electrons);
                X15(i,j) = zeros(n_Electrons, n_Electrons);
                X6(i,j) = zeros(n_Electrons, Matrix_Size);

                W2(i,j) = zeros(n_Electrons, n_Electrons);
                W3(i,j) = zeros(n_Electrons, n_Electrons);
                W4(i,j) = zeros(n_Electrons, Matrix_Size);
            }

            if ((j < n_Electrons) && (i < n_Electrons))
            {
                X3(i,j) = zeros(n_Electrons, n_Electrons);
                W1(i,j) = zeros(n_Electrons, n_Electrons);
            }
        }
    }

    // Initial Guess
    for (int i = 0; i < n_Electrons; i++)
    {
        for (int a = n_Electrons; a < Matrix_Size; a++)
        {
            den_ai(a,i) = fs(i,i) - fs(a,a);
        }
    }

    for (int a = n_Electrons; a < Matrix_Size; a++)
    {
        for (int b = n_Electrons; b < Matrix_Size; b++)
        {
            for (int i = 0; i < n_Electrons; i++)
            {
                for (int j = 0; j < n_Electrons; j++)
                {
                    t2(a,b)(i,j) = MOs(a,b)(i,j) / (den_ai(a,i) + den_ai(b,j));
                }
            }
        }
    }

    // Pre iteration
    Fill_tau3();
    E_new = Update_Energy();
    cout << "Starting CCSDT iterations: " << endl;
    cout << "Energy = " << E_new << " iter: " << iter << endl;

    time_t start, slutt;
    // Start iterating
    double check;
    while (continue_ccsdt == true)
    {
        start = clock();
        // Update old energy
        E_old = E_new;
        iter += 1;

        // Update intermediates for t1 and t2
        Fill_F1();
        Fill_F2();
        Fill_F3();
        Fill_W1();
        Fill_W2();
        Fill_W3();
        Fill_W4();

        // CCSDT intermediates

        // 1a & 1b
        if (Method_Nr >= 0)
        {
            Fill_X1();
            Fill_X2();
        }

        if (Method_Nr > 2)
        {
            // 2
            Fill_X12();
            Fill_X13();
        }

        if (Method_Nr > 3)
        {
            // 3
            Fill_X14();
            Fill_X15();
        }

        if (Method_Nr > 4)
        {
            // 4
            Fill_X6();
            Fill_X3();
            Fill_X4();
        }

        if (Method_Nr > 5)
        {
            // Full
            Fill_X7();
            Fill_X8();
        }

        // Calculate new amplitudes
        Fill_T1();
        Fill_T2();

        if (Method_Nr >= 0)
        {
            Fill_T3();
        }

        // Update amplitudes
        t1 = t1_new;
        t2 = t2_new;
        // t3 = t3_new updated inside Fill_T3()
        Fill_tau3();

        // Find new energy
        E_new = Update_Energy();
        check = sqrt((E_new - E_old) * (E_new - E_old));
        slutt = clock();
        cout << "Energy = " << E_new << " iter: " << iter << " " << " time = " << (double)(slutt-start)/CLOCKS_PER_SEC << endl;

        // Check for convergance, or max iter
        if ((check < toler) || (iter == IterMax))
        {
            continue_ccsdt = false;
        }
    }

    // Put together energies
    Energy = E_HF + E_new;

    // Return CCSDT energy
    return Energy;
}

// Check
double CCSDT::Update_Energy()
{
    double ener = 0;

    for (int a = n_Electrons; a < Matrix_Size; a++)
    {
        for (int i = 0; i < n_Electrons; i++)
        {
            ener += fs(a,i) * t1(a,i);

            for (int b = n_Electrons; b < Matrix_Size; b++)
            {
                for (int j = 0; j < n_Electrons; j++)
                {
                    ener += 0.25 * MOs(i,j)(a,b) * tau3(a,b)(i,j);
                }
            }
        }
    }

    return ener;
}

// CCSD part Check
void CCSDT::Fill_T1()
{
    // Update t1_new med nye amplituder for CCSDT
    // Fjerne noen av leddene her pga her inkluderes CCSDTQ ledd mulig ikke alle disse er med i CCSDT

    if (Method_Nr == 7)
    {

        for (int a = n_Electrons; a < Matrix_Size; a++)
        {
            for (int i = 0; i < n_Electrons; i++)
            {
                    t1_new(a,i) = fs(a,i);

                    for (int k = 0; k < n_Electrons; k++)
                    {
                        t1_new(a,i) -= t1(a,k) * D2(k,i);
                    }

                    for (int c = n_Electrons; c < Matrix_Size; c++)
                    {
                        t1_new(a,i) += (1 - EqualFunc(a,c)) * fs(a,c);

                        for (int k = 0; k < n_Electrons; k++)
                        {
                            t1_new(a,i) += MOs(c,i)(k,a) * t1(c,k);
                            t1_new(a,i) += t2(a,c)(i,k) * D1(c,k);

                            for (int l = 0; l < n_Electrons; l++)
                            {
                                t1_new(a,i) -= 0.5 * t2(c,a)(k,l) * W3(c,i)(k,l);
                            }

                            for (int d = n_Electrons; d < Matrix_Size; d++)
                            {
                                t1_new(a,i) += 0.5 * MOs(k,a)(c,d) * tau3(c,d)(k,i);
                            }
                        }
                    }

                    // CCSDT-1a
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


                    t1_new(a,i) = t1_new(a,i) / den_ai(a,i);

            }
        }

    }

    else
    {

    for (int a = n_Electrons; a < Matrix_Size; a++)
    {
        for (int i = 0; i < n_Electrons; i++)
        {
            if ((a+i)%2 == 0)
            {
                t1_new(a,i) = fs(a,i);

                for (int k = 0; k < n_Electrons; k++)
                {
                    t1_new(a,i) -= t1(a,k) * D2(k,i);
                }

                for (int c = n_Electrons; c < Matrix_Size; c++)
                {
                    t1_new(a,i) += (1 - EqualFunc(a,c)) * fs(a,c);

                    for (int k = 0; k < n_Electrons; k++)
                    {
                        t1_new(a,i) += MOs(c,i)(k,a) * t1(c,k);
                        t1_new(a,i) += t2(a,c)(i,k) * D1(c,k);

                        for (int l = 0; l < n_Electrons; l++)
                        {
                            t1_new(a,i) -= 0.5 * t2(c,a)(k,l) * W3(c,i)(k,l);
                        }

                        for (int d = n_Electrons; d < Matrix_Size; d++)
                        {
                            t1_new(a,i) += 0.5 * MOs(k,a)(c,d) * tau3(c,d)(k,i);
                        }
                    }
                }

                if (Method_Nr >= 0)
                {
                // CCSDT-1a
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
                }


                t1_new(a,i) = t1_new(a,i) / den_ai(a,i);
            }
        }
    }

    }
}

// CCSD part Check
void CCSDT::Fill_T2()
{
    // Update t2_new med nye amplituder for CCSDT
    // Fjerne noen av disse leddene pga her inkluderes alt fra CCSDTQ unntatt T4 amplitude leddene

    if (Method_Nr == 7)
    {
        for (int a = n_Electrons; a < Matrix_Size; a++)
        {
            for (int b = a+1; b < Matrix_Size; b++)
            {
                for (int i = 0; i < n_Electrons; i++)
                {
                    for (int j = i+1; j< n_Electrons; j++)
                    {
                            t2_new(a,b)(i,j) = MOs(i,j)(a,b);

                            for (int k = 0; k < n_Electrons; k++)
                            {
                                for (int l = 0; l < n_Electrons; l++)
                                {
                                    t2_new(a,b)(i,j) += 0.5 * tau3(a,b)(k,l) * W1(i,j)(k,l);
                                }

                                t2_new(a,b)(i,j) -= t1(b,k) * W2(a,k)(i,j);
                                t2_new(a,b)(i,j) += t1(a,k) * W2(b,k)(i,j);

                                t2_new(a,b)(i,j) += t2(a,b)(j,k) * D2(k,i);
                                t2_new(a,b)(i,j) -= t2(a,b)(i,k) * D2(k,j);

                                for (int c = n_Electrons; c < Matrix_Size; c++)
                                {
                                    t2_new(a,b)(i,j) += t2(b,c)(j,k) * W4(a,k)(i,c);
                                    t2_new(a,b)(i,j) -= t2(a,c)(j,k) * W4(b,k)(i,c);
                                    t2_new(a,b)(i,j) -= t2(b,c)(i,k) * W4(a,k)(j,c);
                                    t2_new(a,b)(i,j) += t2(a,c)(i,k) * W4(b,k)(j,c);
                                }
                            }

                            for (int c = n_Electrons; c < Matrix_Size; c++)
                            {
                                t2_new(a,b)(i,j) -= t2(b,c)(i,j) * D3(a,c);
                                t2_new(a,b)(i,j) += t2(a,c)(i,j) * D3(b,c);

                                t2_new(a,b)(i,j) += MOs(a,b)(c,j) * t1(c,i);
                                t2_new(a,b)(i,j) -= MOs(a,b)(c,i) * t1(c,j);

                                for (int d = n_Electrons; d < Matrix_Size; d++)
                                {
                                    t2_new(a,b)(i,j) += 0.5 * MOs(a,b)(c,d) * tau3(c,d)(i,j);
                                }
                            }


                            /////////////////////
                            // CCSDT new terms //
                            /////////////////////

                            // Term in CCSDT-1a

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

                            // Term in CCSDT-1a

                            // - 1/2 P(ij) <kl||jc> t_abcikl
                            for (int c = n_Electrons; c < Matrix_Size; c++)
                            {
                                for (int m = 0; m < n_Electrons; m++)
                                {
                                    for (int k = 0; k < n_Electrons; k++)
                                    {
                                        t2_new(a,b)(i,j) -= 0.5* MOs(k,m)(j,c) * t3[a][b][c][i][k][m];
                                        t2_new(a,b)(i,j) += 0.5* MOs(k,m)(i,c) * t3[a][b][c][j][k][m];
                                    }
                                }
                            }

                            if (Method_Nr > 1)
                            {
                                // Term in CCSDT-1b

                                // - 1/2 P(ij) <kl||cd> t_adbklj t_ci
                                // - 1/2 P(ab) <kl||cd> t_acdikj t_bl
                                // + <kl||cd> t_abcijk t_dl
                                for (int c = n_Electrons; c < Matrix_Size; c++)
                                {
                                    for (int d = n_Electrons; d < Matrix_Size; d++)
                                    {
                                        for (int k = 0; k < n_Electrons; k++)
                                        {
                                            for (int l = 0; l < n_Electrons; l++)
                                            {
                                                t2_new(a,b)(i,j) += 0.5 * MOs(k,l)(c,d) *
                                                        (
                                                             -t3[a][c][d][i][k][j] * t1(b,l)
                                                             + t3[b][c][d][i][k][j] * t1(a,l)

                                                             - t3[a][d][b][k][l][j] * t1(c,i)
                                                             + t3[a][d][b][k][l][i] * t1(c,j)

                                                             + 2*t3[a][b][c][i][j][k] * t1(d,l)
                                                         );
                                            }
                                        }
                                    }
                                }

                                // This is ZERO in case where HF is basis
                                // + f_kc t_abcijk
                                for (int c = n_Electrons; c < Matrix_Size; c++)
                                {
                                    for (int k = 0; k < n_Electrons; k++)
                                    {
                                        t2_new(a,b)(i,j) += fs(k,c) * t3[a][b][c][i][j][k];
                                    }
                                }
                            }

                            // Denominator
                            t2_new(a,b)(i,j) = t2_new(a,b)(i,j) / (den_ai(a,i) + den_ai(b,j));


                    }
                }
            }
        }
    }

    else
    {


        for (int a = n_Electrons; a < Matrix_Size; a++)
        {
            for (int b = a+1; b < Matrix_Size; b++)
            {
                for (int i = 0; i < n_Electrons; i++)
                {
                    for (int j = i+1; j< n_Electrons; j++)
                    {
                        if ((a+b+i+j)%2 == 0)
                        {
                            t2_new(a,b)(i,j) = MOs(i,j)(a,b);

                            for (int k = 0; k < n_Electrons; k++)
                            {
                                for (int l = 0; l < n_Electrons; l++)
                                {
                                    t2_new(a,b)(i,j) += 0.5 * tau3(a,b)(k,l) * W1(i,j)(k,l);
                                }

                                t2_new(a,b)(i,j) -= t1(b,k) * W2(a,k)(i,j);
                                t2_new(a,b)(i,j) += t1(a,k) * W2(b,k)(i,j);

                                t2_new(a,b)(i,j) += t2(a,b)(j,k) * D2(k,i);
                                t2_new(a,b)(i,j) -= t2(a,b)(i,k) * D2(k,j);

                                for (int c = n_Electrons; c < Matrix_Size; c++)
                                {
                                    t2_new(a,b)(i,j) += t2(b,c)(j,k) * W4(a,k)(i,c);
                                    t2_new(a,b)(i,j) -= t2(a,c)(j,k) * W4(b,k)(i,c);
                                    t2_new(a,b)(i,j) -= t2(b,c)(i,k) * W4(a,k)(j,c);
                                    t2_new(a,b)(i,j) += t2(a,c)(i,k) * W4(b,k)(j,c);
                                }
                            }

                            for (int c = n_Electrons; c < Matrix_Size; c++)
                            {
                                t2_new(a,b)(i,j) -= t2(b,c)(i,j) * D3(a,c);
                                t2_new(a,b)(i,j) += t2(a,c)(i,j) * D3(b,c);

                                t2_new(a,b)(i,j) += MOs(a,b)(c,j) * t1(c,i);
                                t2_new(a,b)(i,j) -= MOs(a,b)(c,i) * t1(c,j);

                                for (int d = n_Electrons; d < Matrix_Size; d++)
                                {
                                    t2_new(a,b)(i,j) += 0.5 * MOs(a,b)(c,d) * tau3(c,d)(i,j);
                                }
                            }


                            if (Method_Nr >= 0)
                            {
                            /////////////////////
                            // CCSDT new terms //
                            /////////////////////

                            // Term in CCSDT-1a

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

                            // Term in CCSDT-1a

                            // - 1/2 P(ij) <kl||jc> t_abcikl
                            for (int c = n_Electrons; c < Matrix_Size; c++)
                            {
                                for (int m = 0; m < n_Electrons; m++)
                                {
                                    for (int k = 0; k < n_Electrons; k++)
                                    {
                                        t2_new(a,b)(i,j) -= 0.5* MOs(k,m)(j,c) * t3[a][b][c][i][k][m];
                                        t2_new(a,b)(i,j) += 0.5* MOs(k,m)(i,c) * t3[a][b][c][j][k][m];
                                    }
                                }
                            }

                            if (Method_Nr > 1)
                            {
                                // Term in CCSDT-1b

                                // - 1/2 P(ij) <kl||cd> t_adbklj t_ci
                                // - 1/2 P(ab) <kl||cd> t_acdikj t_bl
                                // + <kl||cd> t_abcijk t_dl
                                for (int c = n_Electrons; c < Matrix_Size; c++)
                                {
                                    for (int d = n_Electrons; d < Matrix_Size; d++)
                                    {
                                        for (int k = 0; k < n_Electrons; k++)
                                        {
                                            for (int l = 0; l < n_Electrons; l++)
                                            {
                                                t2_new(a,b)(i,j) += 0.5 * MOs(k,l)(c,d) *
                                                        (
                                                             -t3[a][c][d][i][k][j] * t1(b,l)
                                                             + t3[b][c][d][i][k][j] * t1(a,l)

                                                             - t3[a][d][b][k][l][j] * t1(c,i)
                                                             + t3[a][d][b][k][l][i] * t1(c,j)

                                                             + 2*t3[a][b][c][i][j][k] * t1(d,l)
                                                         );
                                            }
                                        }
                                    }
                                }

                                // This is ZERO in case where HF is basis
                                // + f_kc t_abcijk
                                for (int c = n_Electrons; c < Matrix_Size; c++)
                                {
                                    for (int k = 0; k < n_Electrons; k++)
                                    {
                                        t2_new(a,b)(i,j) += fs(k,c) * t3[a][b][c][i][j][k];
                                    }
                                }
                            }

                            }
                            // Denominator
                            t2_new(a,b)(i,j) = t2_new(a,b)(i,j) / (den_ai(a,i) + den_ai(b,j));

                            // Symmetries
                            t2_new(b,a)(i,j) = -t2_new(a,b)(i,j);
                            t2_new(a,b)(j,i) = -t2_new(a,b)(i,j);
                            t2_new(b,a)(j,i) = t2_new(a,b)(i,j);
                        }
                    }
                }
            }
        }

    }
}

// Correct!
void CCSDT::Fill_T3()
{
    // Oppdater t3_new med nye amplituder for CCSDT
    double temp;

    if (Method_Nr == 7)
    {
        for (int a = n_Electrons; a < Matrix_Size; a++)
        {
            for (int b = n_Electrons; b < Matrix_Size; b++)
            {
                for (int c = n_Electrons; c < Matrix_Size; c++)
                {
                    for (int i = 0; i < n_Electrons; i++)
                    {
                        for (int j = 0; j < n_Electrons; j++)
                        {
                            for (int k = 0; k < n_Electrons; k++)
                            {
                                    t3_new[a][b][c][i][j][k] = 0;

                                    // CCSDT-1/2/3
                                    for (int e = n_Electrons; e < Matrix_Size; e++)
                                    {
                                        // P(a/bc) P(k/ij)
                                        t3_new[a][b][c][i][j][k] += X1(b,c)(k,e) * t2(a,e)(i,j);
                                        t3_new[a][b][c][i][j][k] -= X1(a,c)(k,e) * t2(b,e)(i,j);
                                        t3_new[a][b][c][i][j][k] -= X1(b,a)(k,e) * t2(c,e)(i,j);

                                        t3_new[a][b][c][i][j][k] -= X1(b,c)(i,e) * t2(a,e)(k,j);
                                        t3_new[a][b][c][i][j][k] += X1(a,c)(i,e) * t2(b,e)(k,j);
                                        t3_new[a][b][c][i][j][k] += X1(b,a)(i,e) * t2(c,e)(k,j);

                                        t3_new[a][b][c][i][j][k] -= X1(b,c)(j,e) * t2(a,e)(i,k);
                                        t3_new[a][b][c][i][j][k] += X1(a,c)(j,e) * t2(b,e)(i,k);
                                        t3_new[a][b][c][i][j][k] += X1(b,a)(j,e) * t2(c,e)(i,k);

                                        if (Method_Nr > 3)
                                        {
                                            // P(c/ab) P(i/jk)
                                            t3_new[a][b][c][i][j][k] += X14(a,b)(i,e) * t2(e,c)(j,k);
                                            t3_new[a][b][c][i][j][k] -= X14(c,b)(i,e) * t2(e,a)(j,k);
                                            t3_new[a][b][c][i][j][k] -= X14(a,c)(i,e) * t2(e,b)(j,k);

                                            t3_new[a][b][c][i][j][k] -= X14(a,b)(j,e) * t2(e,c)(i,k);
                                            t3_new[a][b][c][i][j][k] += X14(c,b)(j,e) * t2(e,a)(i,k);
                                            t3_new[a][b][c][i][j][k] += X14(a,c)(j,e) * t2(e,b)(i,k);

                                            t3_new[a][b][c][i][j][k] -= X14(a,b)(k,e) * t2(e,c)(j,i);
                                            t3_new[a][b][c][i][j][k] += X14(c,b)(k,e) * t2(e,a)(j,i);
                                            t3_new[a][b][c][i][j][k] += X14(a,c)(k,e) * t2(e,b)(j,i);
                                        }
                                    }

                                    // CCSDT-1/2/3
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        // P(c/ab) P(i/jk)
                                        t3_new[a][b][c][i][j][k] -= X2(c,m)(j,k) * t2(a,b)(i,m);
                                        t3_new[a][b][c][i][j][k] += X2(a,m)(j,k) * t2(c,b)(i,m);
                                        t3_new[a][b][c][i][j][k] += X2(b,m)(j,k) * t2(a,c)(i,m);

                                        t3_new[a][b][c][i][j][k] += X2(c,m)(i,k) * t2(a,b)(j,m);
                                        t3_new[a][b][c][i][j][k] -= X2(a,m)(i,k) * t2(c,b)(j,m);
                                        t3_new[a][b][c][i][j][k] -= X2(b,m)(i,k) * t2(a,c)(j,m);

                                        t3_new[a][b][c][i][j][k] += X2(c,m)(j,i) * t2(a,b)(k,m);
                                        t3_new[a][b][c][i][j][k] -= X2(a,m)(j,i) * t2(c,b)(k,m);
                                        t3_new[a][b][c][i][j][k] -= X2(b,m)(j,i) * t2(a,c)(k,m);

                                        if (Method_Nr > 3)
                                        {
                                            // P(a/bc) P(k/ij)
                                            t3_new[a][b][c][i][j][k] += X15(a,m)(i,j) * t2(b,c)(m,k);
                                            t3_new[a][b][c][i][j][k] -= X15(b,m)(i,j) * t2(a,c)(m,k);
                                            t3_new[a][b][c][i][j][k] -= X15(c,m)(i,j) * t2(b,a)(m,k);

                                            t3_new[a][b][c][i][j][k] -= X15(a,m)(k,j) * t2(b,c)(m,i);
                                            t3_new[a][b][c][i][j][k] += X15(b,m)(k,j) * t2(a,c)(m,i);
                                            t3_new[a][b][c][i][j][k] += X15(c,m)(k,j) * t2(b,a)(m,i);

                                            t3_new[a][b][c][i][j][k] -= X15(a,m)(i,k) * t2(b,c)(m,j);
                                            t3_new[a][b][c][i][j][k] += X15(b,m)(i,k) * t2(a,c)(m,j);
                                            t3_new[a][b][c][i][j][k] += X15(c,m)(i,k) * t2(b,a)(m,j);
                                        }
                                    }

                                    if (Method_Nr > 2)
                                    {
                                        // CCSDT-3
                                        for (int e = n_Electrons; e < Matrix_Size; e++)
                                        {
                                            // P(i/jk) P(abc)
                                            t3_new[a][b][c][i][j][k] += X12(a,b)(i,e) * t2(e,c)(j,k);
                                            t3_new[a][b][c][i][j][k] -= X12(b,a)(i,e) * t2(e,c)(j,k);
                                            t3_new[a][b][c][i][j][k] -= X12(c,b)(i,e) * t2(e,a)(j,k);
                                            t3_new[a][b][c][i][j][k] -= X12(a,c)(i,e) * t2(e,b)(j,k);
                                            t3_new[a][b][c][i][j][k] += X12(b,c)(i,e) * t2(e,a)(j,k);
                                            t3_new[a][b][c][i][j][k] += X12(c,a)(i,e) * t2(e,b)(j,k);

                                            t3_new[a][b][c][i][j][k] -= X12(a,b)(j,e) * t2(e,c)(i,k);
                                            t3_new[a][b][c][i][j][k] += X12(b,a)(j,e) * t2(e,c)(i,k);
                                            t3_new[a][b][c][i][j][k] += X12(c,b)(j,e) * t2(e,a)(i,k);
                                            t3_new[a][b][c][i][j][k] += X12(a,c)(j,e) * t2(e,b)(i,k);
                                            t3_new[a][b][c][i][j][k] -= X12(b,c)(j,e) * t2(e,a)(i,k);
                                            t3_new[a][b][c][i][j][k] -= X12(c,a)(j,e) * t2(e,b)(i,k);

                                            t3_new[a][b][c][i][j][k] -= X12(a,b)(k,e) * t2(e,c)(j,i);
                                            t3_new[a][b][c][i][j][k] += X12(b,a)(k,e) * t2(e,c)(j,i);
                                            t3_new[a][b][c][i][j][k] += X12(c,b)(k,e) * t2(e,a)(j,i);
                                            t3_new[a][b][c][i][j][k] += X12(a,c)(k,e) * t2(e,b)(j,i);
                                            t3_new[a][b][c][i][j][k] -= X12(b,c)(k,e) * t2(e,a)(j,i);
                                            t3_new[a][b][c][i][j][k] -= X12(c,a)(k,e) * t2(e,b)(j,i);
                                        }

                                        // CCSDT-3
                                        for (int m = 0; m < n_Electrons; m++)
                                        {
                                            t3_new[a][b][c][i][j][k] += X13(a,m)(i,j) * t2(b,c)(m,k);
                                            t3_new[a][b][c][i][j][k] -= X13(a,m)(j,i) * t2(b,c)(m,k);
                                            t3_new[a][b][c][i][j][k] -= X13(a,m)(k,j) * t2(b,c)(m,i);
                                            t3_new[a][b][c][i][j][k] -= X13(a,m)(i,k) * t2(b,c)(m,j);
                                            t3_new[a][b][c][i][j][k] += X13(a,m)(j,k) * t2(b,c)(m,i);
                                            t3_new[a][b][c][i][j][k] += X13(a,m)(k,i) * t2(b,c)(m,j);

                                            t3_new[a][b][c][i][j][k] -= X13(b,m)(i,j) * t2(a,c)(m,k);
                                            t3_new[a][b][c][i][j][k] += X13(b,m)(j,i) * t2(a,c)(m,k);
                                            t3_new[a][b][c][i][j][k] += X13(b,m)(k,j) * t2(a,c)(m,i);
                                            t3_new[a][b][c][i][j][k] += X13(b,m)(i,k) * t2(a,c)(m,j);
                                            t3_new[a][b][c][i][j][k] -= X13(b,m)(j,k) * t2(a,c)(m,i);
                                            t3_new[a][b][c][i][j][k] -= X13(b,m)(k,i) * t2(a,c)(m,j);

                                            t3_new[a][b][c][i][j][k] -= X13(c,m)(i,j) * t2(b,a)(m,k);
                                            t3_new[a][b][c][i][j][k] += X13(c,m)(j,i) * t2(b,a)(m,k);
                                            t3_new[a][b][c][i][j][k] += X13(c,m)(k,j) * t2(b,a)(m,i);
                                            t3_new[a][b][c][i][j][k] += X13(c,m)(i,k) * t2(b,a)(m,j);
                                            t3_new[a][b][c][i][j][k] -= X13(c,m)(j,k) * t2(b,a)(m,i);
                                            t3_new[a][b][c][i][j][k] -= X13(c,m)(k,i) * t2(b,a)(m,j);
                                        }
                                    }

                                    // CCSDT-4 stuff, Correct Up Till Here!

                                    if (Method_Nr > 4)
                                    {
                                        // T3
                                        for (int e = n_Electrons; e < Matrix_Size; e++)
                                        {
                                            for (int d = n_Electrons; d < Matrix_Size; d++)
                                            {
                                                // P(c/ab) X4(a,b)(d,e)
                                                t3_new[a][b][c][i][j][k] += X4(a,b)(d,e) * t3[d][e][c][i][j][k];
                                                t3_new[a][b][c][i][j][k] -= X4(c,b)(d,e) * t3[d][e][a][i][j][k];
                                                t3_new[a][b][c][i][j][k] -= X4(a,c)(d,e) * t3[d][e][b][i][j][k];
                                            }
                                        }

                                        // T3
                                        for (int m = 0; m < n_Electrons; m++)
                                        {
                                            for (int l = 0; l < n_Electrons; l++)
                                            {
                                                // P(k/ij) X(i,j)(l,m)
                                                t3_new[a][b][c][i][j][k] += X3(i,j)(l,m) * t3[a][b][c][l][m][k];
                                                t3_new[a][b][c][i][j][k] -= X3(k,j)(l,m) * t3[a][b][c][l][m][i];
                                                t3_new[a][b][c][i][j][k] -= X3(i,k)(l,m) * t3[a][b][c][l][m][j];
                                            }
                                        }

                                        // T3
                                        for (int l = 0; l < n_Electrons; l++)
                                        {
                                            for (int d = n_Electrons; d < Matrix_Size; d++)
                                            {

                                                // P(i/jk) P(a/bc)
                                                t3_new[a][b][c][i][j][k] += X6(a,l)(i,d) * t3[d][b][c][l][j][k];
                                                t3_new[a][b][c][i][j][k] -= X6(b,l)(i,d) * t3[d][a][c][l][j][k];
                                                t3_new[a][b][c][i][j][k] -= X6(c,l)(i,d) * t3[d][b][a][l][j][k];

                                                t3_new[a][b][c][i][j][k] -= X6(a,l)(j,d) * t3[d][b][c][l][i][k];
                                                t3_new[a][b][c][i][j][k] += X6(b,l)(j,d) * t3[d][a][c][l][i][k];
                                                t3_new[a][b][c][i][j][k] += X6(c,l)(j,d) * t3[d][b][a][l][i][k];

                                                t3_new[a][b][c][i][j][k] -= X6(a,l)(k,d) * t3[d][b][c][l][j][i];
                                                t3_new[a][b][c][i][j][k] += X6(b,l)(k,d) * t3[d][a][c][l][j][i];
                                                t3_new[a][b][c][i][j][k] += X6(c,l)(k,d) * t3[d][b][a][l][j][i];

                                            }
                                        }
                                    }

                                    if (Method_Nr > 5)
                                    {
                                        // T3 - 5
                                        for (int d = n_Electrons; d < Matrix_Size; d++)
                                        {
                                            // P(a/bc)
                                            t3_new[a][b][c][i][j][k] += X8(a,d) * t3[d][b][c][i][j][k];
                                            t3_new[a][b][c][i][j][k] -= X8(b,d) * t3[d][a][c][i][j][k];
                                            t3_new[a][b][c][i][j][k] -= X8(c,d) * t3[d][b][a][i][j][k];
                                        }

                                        // T3 - 5
                                        for (int m = 0; m < n_Electrons; m++)
                                        {
                                            // P(i/jk)
                                            t3_new[a][b][c][i][j][k] += X7(i,m) * t3[a][b][c][m][j][k];
                                            t3_new[a][b][c][i][j][k] -= X7(j,m) * t3[a][b][c][m][i][k];
                                            t3_new[a][b][c][i][j][k] -= X7(k,m) * t3[a][b][c][m][j][i];
                                        }
                                    }

                                    /*
                                    // 0 when HF basis
                                    for (int d = n_Electrons; d < Matrix_Size; d++)
                                    {
                                        t3_new[a][b][c][i][j][k] += (1-EqualFunc(c,d)) * fs(c,d) * t3[a][b][d][i][j][k];
                                        t3_new[a][b][c][i][j][k] -= (1-EqualFunc(b,d)) * fs(b,d) * t3[a][c][d][i][j][k];
                                        t3_new[a][b][c][i][j][k] -= (1-EqualFunc(a,d)) * fs(a,d) * t3[c][b][d][i][j][k];
                                    }

                                    // 0 when HF basis
                                    for (int l = 0; l < n_Electrons; l++)
                                    {
                                        t3_new[a][b][c][i][j][k] -= (1 - EqualFunc(l, k)) * fs(l,k) * t3[a][b][c][i][j][l];
                                        t3_new[a][b][c][i][j][k] += (1 - EqualFunc(l, i)) * fs(l,i) * t3[a][b][c][k][j][l];
                                        t3_new[a][b][c][i][j][k] += (1 - EqualFunc(l, j)) * fs(l,j) * t3[a][b][c][i][k][l];

                                        for (int d = n_Electrons; d < Matrix_Size; d++)
                                        {
                                            t3_new[a][b][c][i][j][k] -= fs(l,d) * t1(d,i) * t3[a][b][c][l][j][k];
                                            t3_new[a][b][c][i][j][k] += fs(l,d) * t1(d,j) * t3[a][b][c][l][i][k];
                                            t3_new[a][b][c][i][j][k] += fs(l,d) * t1(d,k) * t3[a][b][c][l][j][i];
                                        }
                                    }
                                    */


                                    if (Method_Nr > 6)
                                    {
                                        // Add "CCSDT--Q" metoden her :-D
                                    }


                                    // Denominator
                                    temp = t3_new[a][b][c][i][j][k] / (den_ai(a,i) + den_ai(b,j) + den_ai(c,k));

                                    // Symmetries
                                    t3_new[a][b][c][i][j][k] = temp;
                            }
                        }
                    }
                }
            }
        }
    }

    else
    {

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
                            if ((i+j+k+a+b+c)%2 == 0)
                            {
                                t3_new[a][b][c][i][j][k] = 0;

                                // CCSDT-1/2/3
                                for (int e = n_Electrons; e < Matrix_Size; e++)
                                {
                                    // P(a/bc) P(k/ij)
                                    t3_new[a][b][c][i][j][k] += X1(b,c)(k,e) * t2(a,e)(i,j);
                                    t3_new[a][b][c][i][j][k] -= X1(a,c)(k,e) * t2(b,e)(i,j);
                                    t3_new[a][b][c][i][j][k] -= X1(b,a)(k,e) * t2(c,e)(i,j);

                                    t3_new[a][b][c][i][j][k] -= X1(b,c)(i,e) * t2(a,e)(k,j);
                                    t3_new[a][b][c][i][j][k] += X1(a,c)(i,e) * t2(b,e)(k,j);
                                    t3_new[a][b][c][i][j][k] += X1(b,a)(i,e) * t2(c,e)(k,j);

                                    t3_new[a][b][c][i][j][k] -= X1(b,c)(j,e) * t2(a,e)(i,k);
                                    t3_new[a][b][c][i][j][k] += X1(a,c)(j,e) * t2(b,e)(i,k);
                                    t3_new[a][b][c][i][j][k] += X1(b,a)(j,e) * t2(c,e)(i,k);

                                    if (Method_Nr > 3)
                                    {
                                        // P(c/ab) P(i/jk)
                                        t3_new[a][b][c][i][j][k] += X14(a,b)(i,e) * t2(e,c)(j,k);
                                        t3_new[a][b][c][i][j][k] -= X14(c,b)(i,e) * t2(e,a)(j,k);
                                        t3_new[a][b][c][i][j][k] -= X14(a,c)(i,e) * t2(e,b)(j,k);

                                        t3_new[a][b][c][i][j][k] -= X14(a,b)(j,e) * t2(e,c)(i,k);
                                        t3_new[a][b][c][i][j][k] += X14(c,b)(j,e) * t2(e,a)(i,k);
                                        t3_new[a][b][c][i][j][k] += X14(a,c)(j,e) * t2(e,b)(i,k);

                                        t3_new[a][b][c][i][j][k] -= X14(a,b)(k,e) * t2(e,c)(j,i);
                                        t3_new[a][b][c][i][j][k] += X14(c,b)(k,e) * t2(e,a)(j,i);
                                        t3_new[a][b][c][i][j][k] += X14(a,c)(k,e) * t2(e,b)(j,i);
                                    }
                                }

                                // CCSDT-1/2/3
                                for (int m = 0; m < n_Electrons; m++)
                                {
                                    // P(c/ab) P(i/jk)
                                    t3_new[a][b][c][i][j][k] -= X2(c,m)(j,k) * t2(a,b)(i,m);
                                    t3_new[a][b][c][i][j][k] += X2(a,m)(j,k) * t2(c,b)(i,m);
                                    t3_new[a][b][c][i][j][k] += X2(b,m)(j,k) * t2(a,c)(i,m);

                                    t3_new[a][b][c][i][j][k] += X2(c,m)(i,k) * t2(a,b)(j,m);
                                    t3_new[a][b][c][i][j][k] -= X2(a,m)(i,k) * t2(c,b)(j,m);
                                    t3_new[a][b][c][i][j][k] -= X2(b,m)(i,k) * t2(a,c)(j,m);

                                    t3_new[a][b][c][i][j][k] += X2(c,m)(j,i) * t2(a,b)(k,m);
                                    t3_new[a][b][c][i][j][k] -= X2(a,m)(j,i) * t2(c,b)(k,m);
                                    t3_new[a][b][c][i][j][k] -= X2(b,m)(j,i) * t2(a,c)(k,m);

                                    if (Method_Nr > 3)
                                    {
                                        // P(a/bc) P(k/ij)
                                        t3_new[a][b][c][i][j][k] += X15(a,m)(i,j) * t2(b,c)(m,k);
                                        t3_new[a][b][c][i][j][k] -= X15(b,m)(i,j) * t2(a,c)(m,k);
                                        t3_new[a][b][c][i][j][k] -= X15(c,m)(i,j) * t2(b,a)(m,k);

                                        t3_new[a][b][c][i][j][k] -= X15(a,m)(k,j) * t2(b,c)(m,i);
                                        t3_new[a][b][c][i][j][k] += X15(b,m)(k,j) * t2(a,c)(m,i);
                                        t3_new[a][b][c][i][j][k] += X15(c,m)(k,j) * t2(b,a)(m,i);

                                        t3_new[a][b][c][i][j][k] -= X15(a,m)(i,k) * t2(b,c)(m,j);
                                        t3_new[a][b][c][i][j][k] += X15(b,m)(i,k) * t2(a,c)(m,j);
                                        t3_new[a][b][c][i][j][k] += X15(c,m)(i,k) * t2(b,a)(m,j);
                                    }
                                }

                                if (Method_Nr > 2)
                                {
                                    // CCSDT-3
                                    for (int e = n_Electrons; e < Matrix_Size; e++)
                                    {
                                        // P(i/jk) P(abc)
                                        t3_new[a][b][c][i][j][k] += X12(a,b)(i,e) * t2(e,c)(j,k);
                                        t3_new[a][b][c][i][j][k] -= X12(b,a)(i,e) * t2(e,c)(j,k);
                                        t3_new[a][b][c][i][j][k] -= X12(c,b)(i,e) * t2(e,a)(j,k);
                                        t3_new[a][b][c][i][j][k] -= X12(a,c)(i,e) * t2(e,b)(j,k);
                                        t3_new[a][b][c][i][j][k] += X12(b,c)(i,e) * t2(e,a)(j,k);
                                        t3_new[a][b][c][i][j][k] += X12(c,a)(i,e) * t2(e,b)(j,k);

                                        t3_new[a][b][c][i][j][k] -= X12(a,b)(j,e) * t2(e,c)(i,k);
                                        t3_new[a][b][c][i][j][k] += X12(b,a)(j,e) * t2(e,c)(i,k);
                                        t3_new[a][b][c][i][j][k] += X12(c,b)(j,e) * t2(e,a)(i,k);
                                        t3_new[a][b][c][i][j][k] += X12(a,c)(j,e) * t2(e,b)(i,k);
                                        t3_new[a][b][c][i][j][k] -= X12(b,c)(j,e) * t2(e,a)(i,k);
                                        t3_new[a][b][c][i][j][k] -= X12(c,a)(j,e) * t2(e,b)(i,k);

                                        t3_new[a][b][c][i][j][k] -= X12(a,b)(k,e) * t2(e,c)(j,i);
                                        t3_new[a][b][c][i][j][k] += X12(b,a)(k,e) * t2(e,c)(j,i);
                                        t3_new[a][b][c][i][j][k] += X12(c,b)(k,e) * t2(e,a)(j,i);
                                        t3_new[a][b][c][i][j][k] += X12(a,c)(k,e) * t2(e,b)(j,i);
                                        t3_new[a][b][c][i][j][k] -= X12(b,c)(k,e) * t2(e,a)(j,i);
                                        t3_new[a][b][c][i][j][k] -= X12(c,a)(k,e) * t2(e,b)(j,i);
                                    }

                                    // CCSDT-3
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        t3_new[a][b][c][i][j][k] += X13(a,m)(i,j) * t2(b,c)(m,k);
                                        t3_new[a][b][c][i][j][k] -= X13(a,m)(j,i) * t2(b,c)(m,k);
                                        t3_new[a][b][c][i][j][k] -= X13(a,m)(k,j) * t2(b,c)(m,i);
                                        t3_new[a][b][c][i][j][k] -= X13(a,m)(i,k) * t2(b,c)(m,j);
                                        t3_new[a][b][c][i][j][k] += X13(a,m)(j,k) * t2(b,c)(m,i);
                                        t3_new[a][b][c][i][j][k] += X13(a,m)(k,i) * t2(b,c)(m,j);

                                        t3_new[a][b][c][i][j][k] -= X13(b,m)(i,j) * t2(a,c)(m,k);
                                        t3_new[a][b][c][i][j][k] += X13(b,m)(j,i) * t2(a,c)(m,k);
                                        t3_new[a][b][c][i][j][k] += X13(b,m)(k,j) * t2(a,c)(m,i);
                                        t3_new[a][b][c][i][j][k] += X13(b,m)(i,k) * t2(a,c)(m,j);
                                        t3_new[a][b][c][i][j][k] -= X13(b,m)(j,k) * t2(a,c)(m,i);
                                        t3_new[a][b][c][i][j][k] -= X13(b,m)(k,i) * t2(a,c)(m,j);

                                        t3_new[a][b][c][i][j][k] -= X13(c,m)(i,j) * t2(b,a)(m,k);
                                        t3_new[a][b][c][i][j][k] += X13(c,m)(j,i) * t2(b,a)(m,k);
                                        t3_new[a][b][c][i][j][k] += X13(c,m)(k,j) * t2(b,a)(m,i);
                                        t3_new[a][b][c][i][j][k] += X13(c,m)(i,k) * t2(b,a)(m,j);
                                        t3_new[a][b][c][i][j][k] -= X13(c,m)(j,k) * t2(b,a)(m,i);
                                        t3_new[a][b][c][i][j][k] -= X13(c,m)(k,i) * t2(b,a)(m,j);
                                    }
                                }

                                // CCSDT-4 stuff, Correct Up Till Here!

                                if (Method_Nr > 4)
                                {
                                    // T3
                                    for (int e = n_Electrons; e < Matrix_Size; e++)
                                    {
                                        for (int d = n_Electrons; d < Matrix_Size; d++)
                                        {
                                            // P(c/ab) X4(a,b)(d,e)
                                            t3_new[a][b][c][i][j][k] += X4(a,b)(d,e) * t3[d][e][c][i][j][k];
                                            t3_new[a][b][c][i][j][k] -= X4(c,b)(d,e) * t3[d][e][a][i][j][k];
                                            t3_new[a][b][c][i][j][k] -= X4(a,c)(d,e) * t3[d][e][b][i][j][k];
                                        }
                                    }

                                    // T3
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        for (int l = 0; l < n_Electrons; l++)
                                        {
                                            // P(k/ij) X(i,j)(l,m)
                                            t3_new[a][b][c][i][j][k] += X3(i,j)(l,m) * t3[a][b][c][l][m][k];
                                            t3_new[a][b][c][i][j][k] -= X3(k,j)(l,m) * t3[a][b][c][l][m][i];
                                            t3_new[a][b][c][i][j][k] -= X3(i,k)(l,m) * t3[a][b][c][l][m][j];
                                        }
                                    }

                                    // T3
                                    for (int l = 0; l < n_Electrons; l++)
                                    {
                                        for (int d = n_Electrons; d < Matrix_Size; d++)
                                        {

                                            // P(i/jk) P(a/bc)
                                            t3_new[a][b][c][i][j][k] += X6(a,l)(i,d) * t3[d][b][c][l][j][k];
                                            t3_new[a][b][c][i][j][k] -= X6(b,l)(i,d) * t3[d][a][c][l][j][k];
                                            t3_new[a][b][c][i][j][k] -= X6(c,l)(i,d) * t3[d][b][a][l][j][k];

                                            t3_new[a][b][c][i][j][k] -= X6(a,l)(j,d) * t3[d][b][c][l][i][k];
                                            t3_new[a][b][c][i][j][k] += X6(b,l)(j,d) * t3[d][a][c][l][i][k];
                                            t3_new[a][b][c][i][j][k] += X6(c,l)(j,d) * t3[d][b][a][l][i][k];

                                            t3_new[a][b][c][i][j][k] -= X6(a,l)(k,d) * t3[d][b][c][l][j][i];
                                            t3_new[a][b][c][i][j][k] += X6(b,l)(k,d) * t3[d][a][c][l][j][i];
                                            t3_new[a][b][c][i][j][k] += X6(c,l)(k,d) * t3[d][b][a][l][j][i];

                                        }
                                    }
                                }

                                if (Method_Nr > 5)
                                {
                                    // T3 - 5
                                    for (int d = n_Electrons; d < Matrix_Size; d++)
                                    {
                                        // P(a/bc)
                                        t3_new[a][b][c][i][j][k] += X8(a,d) * t3[d][b][c][i][j][k];
                                        t3_new[a][b][c][i][j][k] -= X8(b,d) * t3[d][a][c][i][j][k];
                                        t3_new[a][b][c][i][j][k] -= X8(c,d) * t3[d][b][a][i][j][k];
                                    }

                                    // T3 - 5
                                    for (int m = 0; m < n_Electrons; m++)
                                    {
                                        // P(i/jk)
                                        t3_new[a][b][c][i][j][k] += X7(i,m) * t3[a][b][c][m][j][k];
                                        t3_new[a][b][c][i][j][k] -= X7(j,m) * t3[a][b][c][m][i][k];
                                        t3_new[a][b][c][i][j][k] -= X7(k,m) * t3[a][b][c][m][j][i];
                                    }
                                }

                                /*
                                // 0 when HF basis
                                for (int d = n_Electrons; d < Matrix_Size; d++)
                                {
                                    t3_new[a][b][c][i][j][k] += (1-EqualFunc(c,d)) * fs(c,d) * t3[a][b][d][i][j][k];
                                    t3_new[a][b][c][i][j][k] -= (1-EqualFunc(b,d)) * fs(b,d) * t3[a][c][d][i][j][k];
                                    t3_new[a][b][c][i][j][k] -= (1-EqualFunc(a,d)) * fs(a,d) * t3[c][b][d][i][j][k];
                                }

                                // 0 when HF basis
                                for (int l = 0; l < n_Electrons; l++)
                                {
                                    t3_new[a][b][c][i][j][k] -= (1 - EqualFunc(l, k)) * fs(l,k) * t3[a][b][c][i][j][l];
                                    t3_new[a][b][c][i][j][k] += (1 - EqualFunc(l, i)) * fs(l,i) * t3[a][b][c][k][j][l];
                                    t3_new[a][b][c][i][j][k] += (1 - EqualFunc(l, j)) * fs(l,j) * t3[a][b][c][i][k][l];

                                    for (int d = n_Electrons; d < Matrix_Size; d++)
                                    {
                                        t3_new[a][b][c][i][j][k] -= fs(l,d) * t1(d,i) * t3[a][b][c][l][j][k];
                                        t3_new[a][b][c][i][j][k] += fs(l,d) * t1(d,j) * t3[a][b][c][l][i][k];
                                        t3_new[a][b][c][i][j][k] += fs(l,d) * t1(d,k) * t3[a][b][c][l][j][i];
                                    }
                                }
                                */


                                if (Method_Nr > 6)
                                {
                                    // Add "CCSDT--Q" metoden her :-D
                                }


                                // Denominator
                                temp = t3_new[a][b][c][i][j][k] / (den_ai(a,i) + den_ai(b,j) + den_ai(c,k));

                                // Symmetries
                                t3_new[a][b][c][i][j][k] = temp;
                                t3_new[c][a][b][i][j][k] = temp;
                                t3_new[b][c][a][i][j][k] = temp;
                                t3_new[b][a][c][i][j][k] = -temp;
                                t3_new[c][b][a][i][j][k] = -temp;
                                t3_new[a][c][b][i][j][k] = -temp;

                                t3_new[a][b][c][j][k][i] = temp;
                                t3_new[c][a][b][j][k][i] = temp;
                                t3_new[b][c][a][j][k][i] = temp;
                                t3_new[b][a][c][j][k][i] = -temp;
                                t3_new[c][b][a][j][k][i] = -temp;
                                t3_new[a][c][b][j][k][i] = -temp;

                                t3_new[a][b][c][k][i][j] = temp;
                                t3_new[c][a][b][k][i][j] = temp;
                                t3_new[b][c][a][k][i][j] = temp;
                                t3_new[b][a][c][k][i][j] = -temp;
                                t3_new[c][b][a][k][i][j] = -temp;
                                t3_new[a][c][b][k][i][j] = -temp;

                                t3_new[a][b][c][j][i][k] = -temp;
                                t3_new[c][a][b][j][i][k] = -temp;
                                t3_new[b][c][a][j][i][k] = -temp;
                                t3_new[b][a][c][j][i][k] = temp;
                                t3_new[c][b][a][j][i][k] = temp;
                                t3_new[a][c][b][j][i][k] = temp;

                                t3_new[a][b][c][k][j][i] = -temp;
                                t3_new[c][a][b][k][j][i] = -temp;
                                t3_new[b][c][a][k][j][i] = -temp;
                                t3_new[b][a][c][k][j][i] = temp;
                                t3_new[c][b][a][k][j][i] = temp;
                                t3_new[a][c][b][k][j][i] = temp;

                                t3_new[a][b][c][i][k][j] = -temp;
                                t3_new[c][a][b][i][k][j] = -temp;
                                t3_new[b][c][a][i][k][j] = -temp;
                                t3_new[b][a][c][i][k][j] = temp;
                                t3_new[c][b][a][i][k][j] = temp;
                                t3_new[a][c][b][i][k][j] = temp;
                            }

                            else
                            {
                                t3_new[a][b][c][i][j][k] = 0;
                            }
                        }
                    }
                }
            }
        }
    }

    }

    // Update t3
    for (int a = n_Electrons; a < Matrix_Size; a++)
    {
        for (int b = n_Electrons; b < Matrix_Size; b++)
        {
            for (int c = n_Electrons; c < Matrix_Size; c++)
            {
                for (int i = 0; i < n_Electrons; i++)
                {
                    for (int j = 0; j < n_Electrons; j++)
                    {
                        for (int k = 0; k < n_Electrons; k++)
                        {
                            t3[a][b][c][i][j][k] = t3_new[a][b][c][i][j][k];
                        }
                    }
                }
            }
        }
    }
}

// Correct!
int CCSDT::EqualFunc(int a, int b)
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

// Correct!
void CCSDT::AOtoMO()
{
    int matsize = Matrix_Size / 2;

    MOs.set_size(Matrix_Size, Matrix_Size);

    for (int i = 0; i < Matrix_Size; i++)
    {
        for (int j = 0; j < Matrix_Size; j++)
        {
            MOs(i,j) = zeros(Matrix_Size, Matrix_Size);
        }
    }

    mat recieve_matrix = zeros(matsize, matsize);
    mat temp = zeros(matsize, matsize);
    vector<mat> compact_mo5;
    for (int i = 0; i < matsize; i++)
    {
        compact_mo5.push_back(temp);
    }

    if (Method_Nr == 7)
    {
        c = HartFock->return_C_down();
    }

    else
    {
        c = HartFock->ReturnC();
    }

    field<mat> temp_mo;
    temp_mo.set_size(matsize, matsize);

    field<mat> temp2_mo, temp3_mo, temp4_mo;
    temp2_mo.set_size(matsize, matsize);
    temp3_mo.set_size(matsize, matsize);
    temp4_mo.set_size(matsize, matsize);
    for (int i = 0; i < matsize; i++)
    {
        for (int j = 0; j < matsize; j++)
        {
            temp_mo(i,j) = zeros(matsize, matsize);

            temp2_mo(i,j) = zeros(matsize, matsize);
            temp3_mo(i,j) = zeros(matsize, matsize);
            temp4_mo(i,j) = zeros(matsize, matsize);
        }
    }


    if (Method_Nr != 7)
    {
    for (int a = 0; a < matsize; a++)
    {
        for (int i = 0; i < matsize; i++)
        {
            for (int j = 0; j < matsize; j++)
            {
                // Get AOs
                if (i >= j)
                {
                    recieve_matrix = HartFock->Return_Field_Q(i,j);
                }

                else
                {
                    recieve_matrix = HartFock->Return_Field_Q(j,i);
                }

                for (int k = 0; k < matsize; k++)
                {
                    for (int l = 0; l < matsize; l++)
                    {
                        temp2_mo(a,j)(k,l) += c(i,a) * recieve_matrix(k,l);
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
                        temp3_mo(a,b)(j,k) += c(i,b) * temp2_mo(a,i)(j,k);
                    }
                }
            }
        }
    }

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
                        temp4_mo(a,b)(g,j) += c(i,g) * temp3_mo(a,b)(i,j);
                    }
                }
            }
        }
    }

    for (int a = 0; a < matsize; a++)
    {
        for (int b = 0; b < matsize; b++)
        {
            for (int g = 0; g < matsize; g++)
            {
                for (int d = 0; d < matsize; d++)
                {
                    for (int i = 0; i < matsize; i++)
                    {
                        temp_mo(a,b)(g,d) += c(i,d) * temp4_mo(a,b)(g,i);
                    }
                }
            }
        }
    }

    double val1, val2;
    for (int a = 0; a < Matrix_Size; a++)
    {
        for (int b = 0; b < Matrix_Size; b++)
        {
            for (int i = 0; i < Matrix_Size; i++)
            {
                for (int j = 0; j < Matrix_Size; j++)
                {
                        val1 = EqualFunc(a%2,b%2) * EqualFunc(i%2,j%2) * temp_mo(a/2,b/2)(i/2,j/2);
                        val2 = EqualFunc(a%2,j%2) * EqualFunc(i%2,b%2) * temp_mo(a/2,j/2)(i/2,b/2);
                        MOs(a,i)(b,j) = val1 - val2;

                }
            }
        }
    }
    }

    vec eigval;
    fs = zeros(Matrix_Size, Matrix_Size);

    if (Method_Nr == 7)
    {
        // Unrestricted CCSDT in progress
        // Implementation not finished jet....

        eigval = HartFock->eigenvalues_F_up;
        for (int i = 0; i < Matrix_Size; i++)
        {
            fs(i,i) = eigval(i/2);
            i++;
        }

        eigval = HartFock->eigenvalues_F_down;
        for (int i = 1; i < Matrix_Size; i++)
        {
            fs(i,i) = eigval(i/2);
            i++;
        }
    }



    else
    {
        eigval = HartFock->eigenvalues_F;
        for (int i = 0; i < Matrix_Size; i++)
        {
            fs(i,i) = eigval(i/2);
        }
    }

    if (Method_Nr == 7)
    {
        // Try Unrestricted Spin simple solution

        double val1,val2;
        c = HartFock->return_C_up();

        for (int a = 0; a < matsize; a++)
        {
            for (int i = 0; i < matsize; i++)
            {
                for (int j = 0; j < matsize; j++)
                {
                    // Get AOs
                    if (i >= j)
                    {
                        recieve_matrix = HartFock->Return_Field_Q(i,j);
                    }

                    else
                    {
                        recieve_matrix = HartFock->Return_Field_Q(j,i);
                    }

                    for (int k = 0; k < matsize; k++)
                    {
                        for (int l = 0; l < matsize; l++)
                        {
                            temp2_mo(a,j)(k,l) += c(i,a) * recieve_matrix(k,l);
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
                            temp3_mo(a,b)(j,k) += c(i,b) * temp2_mo(a,i)(j,k);
                        }
                    }
                }
            }
        }

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
                            temp4_mo(a,b)(g,j) += c(i,g) * temp3_mo(a,b)(i,j);
                        }
                    }
                }
            }
        }

        for (int a = 0; a < matsize; a++)
        {
            for (int b = 0; b < matsize; b++)
            {
                for (int g = 0; g < matsize; g++)
                {
                    for (int d = 0; d < matsize; d++)
                    {
                        for (int i = 0; i < matsize; i++)
                        {
                            temp_mo(a,b)(g,d) += c(i,d) * temp4_mo(a,b)(g,i);
                        }
                    }
                }
            }
        }

        for (int a = 0; a < Matrix_Size; a++)
        {
            for (int b = 0; b < Matrix_Size; b++)
            {
                for (int i = 0; i < Matrix_Size; i++)
                {
                    for (int j = 0; j < Matrix_Size; j++)
                    {
                        val1 = EqualFunc(a%2,(b)%2) * EqualFunc(i%2,(j)%2) * temp_mo(a/2,b/2)(i/2,j/2);
                        val2 = EqualFunc(a%2,(j)%2) * EqualFunc(i%2,(b)%2) * temp_mo(a/2,j/2)(i/2,b/2);

                        if ((val1-val2) != 0)
                        {
                            MOs(a,i)(b,j) = val1 - val2;
                        }
                        j++;
                    }
                    i++;
                }
                b++;
            }
            a++;
        }









        c = HartFock->return_C_down();

        for (int a = 0; a < matsize; a++)
        {
            for (int i = 0; i < matsize; i++)
            {
                for (int j = 0; j < matsize; j++)
                {
                    // Get AOs
                    if (i >= j)
                    {
                        recieve_matrix = HartFock->Return_Field_Q(i,j);
                    }

                    else
                    {
                        recieve_matrix = HartFock->Return_Field_Q(j,i);
                    }

                    for (int k = 0; k < matsize; k++)
                    {
                        for (int l = 0; l < matsize; l++)
                        {
                            temp2_mo(a,j)(k,l) += c(i,a) * recieve_matrix(k,l);
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
                            temp3_mo(a,b)(j,k) += c(i,b) * temp2_mo(a,i)(j,k);
                        }
                    }
                }
            }
        }

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
                            temp4_mo(a,b)(g,j) += c(i,g) * temp3_mo(a,b)(i,j);
                        }
                    }
                }
            }
        }

        for (int a = 0; a < matsize; a++)
        {
            for (int b = 0; b < matsize; b++)
            {
                for (int g = 0; g < matsize; g++)
                {
                    for (int d = 0; d < matsize; d++)
                    {
                        for (int i = 0; i < matsize; i++)
                        {
                            temp_mo(a,b)(g,d) += c(i,d) * temp4_mo(a,b)(g,i);
                        }
                    }
                }
            }
        }

        for (int a = 1; a < Matrix_Size; a++)
        {
            for (int b = 1; b < Matrix_Size; b++)
            {
                for (int i = 1; i < Matrix_Size; i++)
                {
                    for (int j = 1; j < Matrix_Size; j++)
                    {
                        val1 = EqualFunc(a%2,(b)%2) * EqualFunc(i%2,(j)%2) * temp_mo(a/2,b/2)(i/2,j/2);
                        val2 = EqualFunc(a%2,(j)%2) * EqualFunc(i%2,(b)%2) * temp_mo(a/2,j/2)(i/2,b/2);
                        MOs(a,i)(b,j) = val1 - val2;
                        j++;
                    }
                    i++;
                }
                b++;
            }
            a++;
        }




        c = HartFock->return_C_up();

        for (int a = 0; a < matsize; a++)
        {
            for (int i = 0; i < matsize; i++)
            {
                for (int j = 0; j < matsize; j++)
                {
                    // Get AOs
                    if (i >= j)
                    {
                        recieve_matrix = HartFock->Return_Field_Q(i,j);
                    }

                    else
                    {
                        recieve_matrix = HartFock->Return_Field_Q(j,i);
                    }

                    for (int k = 0; k < matsize; k++)
                    {
                        for (int l = 0; l < matsize; l++)
                        {
                            temp2_mo(a,j)(k,l) += c(i,a) * recieve_matrix(k,l);
                        }
                    }
                }
            }
        }

        c = HartFock->return_C_down();

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
                            temp3_mo(a,b)(j,k) += c(i,b) * temp2_mo(a,i)(j,k);
                        }
                    }
                }
            }
        }

        c = HartFock->return_C_up();

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
                            temp4_mo(a,b)(g,j) += c(i,g) * temp3_mo(a,b)(i,j);
                        }
                    }
                }
            }
        }

        c = HartFock->return_C_down();

        for (int a = 0; a < matsize; a++)
        {
            for (int b = 0; b < matsize; b++)
            {
                for (int g = 0; g < matsize; g++)
                {
                    for (int d = 0; d < matsize; d++)
                    {
                        for (int i = 0; i < matsize; i++)
                        {
                            temp_mo(a,b)(g,d) += c(i,d) * temp4_mo(a,b)(g,i);
                        }
                    }
                }
            }
        }

        for (int a = 0; a < Matrix_Size; a++)
        {
            for (int b = 1; b < Matrix_Size; b++)
            {
                for (int i = 0; i < Matrix_Size; i++)
                {
                    for (int j = 1; j < Matrix_Size; j++)
                    {
                            val1 = EqualFunc(a%2,(b)%2) * EqualFunc(i%2,j%2) * temp_mo(a/2,b/2)(i/2,j/2);
                            val2 = EqualFunc(a%2,j%2) * EqualFunc(i%2,(b)%2) * temp_mo(a/2,j/2)(i/2,b/2);
                            MOs(a,i)(b,j) = val1 - val2;
                            j++;

                    }
                    i++;
                }
                b++;
            }
            a++;
        }

        c = HartFock->return_C_down();

        for (int a = 0; a < matsize; a++)
        {
            for (int i = 0; i < matsize; i++)
            {
                for (int j = 0; j < matsize; j++)
                {
                    // Get AOs
                    if (i >= j)
                    {
                        recieve_matrix = HartFock->Return_Field_Q(i,j);
                    }

                    else
                    {
                        recieve_matrix = HartFock->Return_Field_Q(j,i);
                    }

                    for (int k = 0; k < matsize; k++)
                    {
                        for (int l = 0; l < matsize; l++)
                        {
                            temp2_mo(a,j)(k,l) += c(i,a) * recieve_matrix(k,l);
                        }
                    }
                }
            }
        }

        c = HartFock->return_C_up();

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
                            temp3_mo(a,b)(j,k) += c(i,b) * temp2_mo(a,i)(j,k);
                        }
                    }
                }
            }
        }

        c = HartFock->return_C_down();

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
                            temp4_mo(a,b)(g,j) += c(i,g) * temp3_mo(a,b)(i,j);
                        }
                    }
                }
            }
        }

        c = HartFock->return_C_up();

        for (int a = 0; a < matsize; a++)
        {
            for (int b = 0; b < matsize; b++)
            {
                for (int g = 0; g < matsize; g++)
                {
                    for (int d = 0; d < matsize; d++)
                    {
                        for (int i = 0; i < matsize; i++)
                        {
                            temp_mo(a,b)(g,d) += c(i,d) * temp4_mo(a,b)(g,i);
                        }
                    }
                }
            }
        }

        for (int a = 1; a < Matrix_Size; a++)
        {
            for (int b = 0; b < Matrix_Size; b++)
            {
                for (int i = 1; i < Matrix_Size; i++)
                {
                    for (int j = 0; j < Matrix_Size; j++)
                    {
                            val1 = EqualFunc(a%2,(b)%2) * EqualFunc(i%2,j%2) * temp_mo(a/2,b/2)(i/2,j/2);
                            val2 = EqualFunc(a%2,j%2) * EqualFunc(i%2,(b)%2) * temp_mo(a/2,j/2)(i/2,b/2);
                            MOs(a,i)(b,j) = val1 - val2;
                            j++;

                    }
                    i++;
                }
                b++;
            }
            a++;
        }
    }
}

// X1, Correct!
void CCSDT::Fill_X1()
{
    for (int c = n_Electrons; c < Matrix_Size; c++)
    {
        for (int b = n_Electrons; b < Matrix_Size; b++)
        {
            for (int k = 0; k < n_Electrons; k++)
            {
                for (int d = n_Electrons; d < Matrix_Size; d++)
                {
                    X1(b,c)(k,d) = MOs(b,c)(d,k);

                    if (Method_Nr > 2)
                    {
                        for (int l = 0; l < n_Electrons; l++)
                        {
                            // ADD THIS
                            //X1(b,c)(k,d) -= fs(l,d) * t2(b,c)(l,k);

                            for (int m = 0; m < n_Electrons; m++)
                            {
                                X1(b,c)(k,d) += 0.5 * MOs(l,m)(d,k) * t2(b,c)(l,m);

                                if (Method_Nr > 5)
                                {
                                    for (int e = n_Electrons; e < Matrix_Size; e++)
                                    {
                                        X1(b,c)(k,d) -= 0.5 * MOs(l,m)(d,e) * t3[b][e][c][l][m][k];
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

// X2, Correct!
void CCSDT::Fill_X2()
{
    for (int c = n_Electrons; c < Matrix_Size; c++)
    {
        for (int l = 0; l < n_Electrons; l++)
        {
            for (int k = 0; k < n_Electrons; k++)
            {
                for (int j = 0; j < n_Electrons; j++)
                {
                    X2(c,l)(j,k) = MOs(l,c)(j,k);

                    if (Method_Nr > 2)
                    {
                        for (int e = n_Electrons; e < Matrix_Size; e++)
                        {
                            for (int d = n_Electrons; d < Matrix_Size; d++)
                            {
                                X2(c,l)(j,k) += 0.5 * t2(d,e)(j,k) * MOs(l,c)(d,e);
                            }
                        }

                        if (Method_Nr > 5)
                        {
                            for (int m = 0; m < n_Electrons; m++)
                            {
                                for (int d = n_Electrons; d < Matrix_Size; d++)
                                {
                                    for (int e = n_Electrons; e < Matrix_Size; e++)
                                    {
                                        X2(c,l)(j,k) += 0.5 * MOs(l,m)(d,e) * t3[d][e][c][j][m][k];
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

//
void CCSDT::Fill_X14()
{
    for (int i = 0; i < n_Electrons; i++)
    {
        for (int a = n_Electrons; a < Matrix_Size; a++)
        {
            for (int b = n_Electrons; b < Matrix_Size; b++)
            {
                for (int d = n_Electrons; d < Matrix_Size; d++)
                {
                    X14(a,b)(i,d) = 0;

                    for (int e = n_Electrons; e < Matrix_Size; e++)
                    {
                        X14(a,b)(i,d) += MOs(a,b)(e,d) * t1(e,i);

                        for (int l = 0; l < n_Electrons; l++)
                        {
                            for (int m = 0; m < n_Electrons; m++)
                            {
                                X14(a,b)(i,d) += MOs(l,m)(e,d) * t1(e,i) * t1(a,l) * t1(b,m);
                                X14(a,b)(i,d) -= MOs(l,m)(e,d) * t1(e,l) * t2(a,b)(i,m);
                                X14(a,b)(i,d) += 0.5*MOs(l,m)(e,d) * t1(e,i) * t2(a,b)(l,m);
                            }
                        }
                    }

                    for (int l = 0; l < n_Electrons; l++)
                    {
                        for (int m = 0; m < n_Electrons; m++)
                        {
                            X14(a,b)(i,d) += MOs(l,m)(i,d) * t1(a,l) * t1(b,m);
                        }
                    }
                }
            }
        }
    }
}

//
void CCSDT::Fill_X15()
{
    for (int a = n_Electrons; a < Matrix_Size; a++)
    {
        for (int l = 0; l < n_Electrons; l++)
        {
            for (int i =0; i < n_Electrons; i++)
            {
                for (int j = 0; j < n_Electrons; j++)
                {
                    X15(a,l)(i,j) = 0;

                    for (int d = n_Electrons; d < Matrix_Size; d++)
                    {
                        //X15(a,l)(i,j) -= fs(m,d) * t2(a,d)(i,j); // 0 in HF

                        for (int e = n_Electrons; e < Matrix_Size; e++)
                        {
                            X15(a,l)(i,j) -= MOs(a,l)(d,e) * t1(d,i) * t1(e,j);
                        }
                    }

                    // P(a/bc) P(k/ij)

                    for (int m = 0; m < n_Electrons; m++)
                    {
                        X15(a,l)(i,j) += MOs(m,l)(i,j) * t1(a,m);

                        for (int d = n_Electrons; d <Matrix_Size; d++)
                        {
                            for (int e = n_Electrons; e < Matrix_Size; e++)
                            {
                                X15(a,l)(i,j) += 0.5 * MOs(m,l)(d,e) * t1(a,m) * tau3(d,e)(i,j);
                               // X15(a,l)(i,j) += MOs(m,l)(d,e) * t1(d,i) * t1(a,m) * t1(e,j);
                            }
                        }
                    }
                }
            }
        }
    }
}

// X7, Correct
void CCSDT::Fill_X7()
{
    for (int i = 0; i < n_Electrons; i++)
    {
        for (int m = 0; m < n_Electrons; m++)
        {
            X7(i,m) = 0;

            for (int l = 0; l < n_Electrons; l++)
            {
                for (int d = n_Electrons; d < Matrix_Size; d++)
                {
                    for (int e = n_Electrons; e < Matrix_Size; e++)
                    {
                        // CONFIRMED
                        X7(i,m) -= 0.5 * MOs(l,m)(d,e) * tau3(d,e)(l,i);
                       // X7(i,m) -= MOs(l,m)(d,e) * t1(d,l) * t1(e,i);
                    }

                    // CONFIRMED
                    X7(i,m) -= MOs(l,m)(d,i) * t1(d,l);
                }
            }
        }
    }
}

// X8, Correct
void CCSDT::Fill_X8()
{
    for (int a = n_Electrons; a < Matrix_Size; a++)
    {
        for (int e = n_Electrons; e < Matrix_Size; e++)
        {
            X8(a,e) = 0;

            for (int l = 0; l < n_Electrons; l++)
            {
                for (int d = n_Electrons; d < Matrix_Size; d++)
                {
                    X8(a,e) += MOs(l,a)(d,e) * t1(d,l);

                    for (int m = 0; m < n_Electrons; m++)
                    {
                        X8(a,e) -= 0.5 * MOs(l,m)(d,e) * tau3(d,a)(l,m);
                        //X8(a,e) -= MOs(l,m)(d,e) * t1(d,l) * t1(a,m);
                    }
                }

                X8(a,e) -= fs(l,e) * t1(a,l);
            }

        }
    }
}

// X3, weird
void CCSDT::Fill_X3()
{
    for (int j = 0; j < n_Electrons; j++)
    {
        for (int i = 0; i < n_Electrons; i++)
        {
            for (int m = 0; m < n_Electrons; m++)
            {
                for (int l = 0; l < n_Electrons; l++)
                {
                    X3(i,j)(l,m) = 0.5 * MOs(i,j)(l,m);

                    if (Method_Nr > 5)
                    {
                        for (int d = n_Electrons; d < Matrix_Size; d++)
                        {
                            for (int e = n_Electrons; d < Matrix_Size; d++)
                            {
                                // 0.5?
                                X3(i,j)(l,m) +=  0.5 * MOs(l,m)(d,e) * tau3(d,e)(i,j);
                            }
                        }

                        for (int d = n_Electrons; d < Matrix_Size; d++)
                        {
                            X3(i,j)(l,m) += MOs(l,m)(d,j) * t1(d,i);
                            X3(i,j)(l,m) -= MOs(l,m)(d,i) * t1(d,j);
                        }
                    }
                }
            }
        }
    }
}

// X4, weird
void CCSDT::Fill_X4()
{
    for (int b = n_Electrons; b < Matrix_Size; b++)
    {
        for (int a = n_Electrons; a < Matrix_Size; a++)
        {
            for (int e = n_Electrons; e < Matrix_Size; e++)
            {
                for (int d = n_Electrons; d < Matrix_Size; d++)
                {
                    X4(a,b)(d,e) = 0.5 * MOs(a,b)(d,e);

                    if (Method_Nr > 5)
                    {
                        for (int m = 0; m < n_Electrons; m++)
                        {
                            for (int l = 0; l < n_Electrons; l++)
                            {
                                // 0.5?
                                X4(a,b)(d,e) += 0.5 * MOs(l,m)(d,e) * tau3(a,b)(l,m);
                            }
                        }

                        for (int l = 0; l < n_Electrons; l++)
                        {
                            X4(a,b)(d,e) += MOs(l,b)(d,e) * t1(a,l);
                            X4(a,b)(d,e) -= MOs(l,a)(d,e) * t1(b,l);
                        }
                    }
                }
            }
        }
    }
}



// X6, Correct!
void CCSDT::Fill_X6()
{
    for (int a = n_Electrons; a < Matrix_Size; a++)
    {
        for (int l = 0; l < n_Electrons; l++)
        {
            for (int i = 0; i < n_Electrons; i++)
            {
                for (int d = n_Electrons; d < Matrix_Size; d++)
                {
                    // P(i/jk) P(a/bc)

                    // Checked
                    X6(a,l)(i,d) = MOs(a,l)(i,d);

                    if (Method_Nr > 5)
                    {
                        for (int e = n_Electrons; e < Matrix_Size; e++)
                        {
                            for (int m = 0; m < n_Electrons; m++)
                            {
                                X6(a,l)(i,d) += MOs(m,l)(e,d) * t2(a,e)(i,m);
                                X6(a,l)(i,d) -= MOs(m,l)(e,d) * t1(e,i) * t1(a,m);
                            }

                            X6(a,l)(i,d) += MOs(a,l)(e,d) * t1(e,i);
                        }

                        for (int m = 0; m < n_Electrons; m++)
                        {
                            // Checked
                            X6(a,l)(i,d) -= MOs(m,l)(i,d) * t1(a,m);
                        }
                    }
                }
            }
        }
    }
}


// Correct
void CCSDT::Fill_X12()
{
    for (int a = n_Electrons; a < Matrix_Size; a++)
    {
        for (int b = n_Electrons; b < Matrix_Size; b++)
        {
            for (int i = 0; i < n_Electrons; i++)
            {
                for (int d = n_Electrons; d < Matrix_Size; d++)
                {
                    X12(a,b)(i,d) = 0;

                    for (int l = 0; l < n_Electrons; l++)
                    {
                        for (int e = n_Electrons; e < Matrix_Size; e++)
                        {
                            X12(a,b)(i,d) += MOs(l,b)(e,d) * t2(a,e)(i,l);
                        }
                    }

                    if (Method_Nr > 3)
                    {
                        for (int l = 0; l < n_Electrons; l++)
                        {
                            X12(a,b)(i,d) -= MOs(a,l)(i,d) * t1(b,l);

                            for (int e = n_Electrons; e < Matrix_Size; e++)
                            {
                                X12(a,b)(i,d) -= MOs(l,b)(e,d) * t1(e,i) * t1(a,l);

                                for (int m = 0; m < n_Electrons; m++)
                                {
                                    X12(a,b)(i,d) -= MOs(l,m)(e,d) * t1(b,m) * t2(a,e)(i,l);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

// Correct
void CCSDT::Fill_X13()
{
    for (int a = n_Electrons; a < Matrix_Size; a++)
    {
        for (int i =0 ; i < n_Electrons; i++)
        {
            for (int j =0; j < n_Electrons; j++)
            {
                for (int l =0 ; l < n_Electrons; l++)
                {
                    X13(a,l)(i,j) = 0;

                    for (int m = 0; m< n_Electrons; m++)
                    {
                        for (int d = n_Electrons; d < Matrix_Size; d++)
                        {
                            X13(a,l)(i,j) -= MOs(m,l)(d,j) * t2(a,d)(i,m);
                        }
                    }

                    if (Method_Nr > 3)
                    {
                        for (int m = 0; m< n_Electrons; m++)
                        {
                            for (int d = n_Electrons; d < Matrix_Size; d++)
                            {
                                X13(a,l)(i,j) += MOs(m,l)(d,j) * t1(d,i) * t1(a,m);

                                for (int e = n_Electrons; e < Matrix_Size; e++)
                                {
                                    X13(a,l)(i,j) -= MOs(m,l)(d,e) * t2(a,d)(i,m) * t1(e,j);
                                }
                            }
                        }

                        for (int d = n_Electrons; d < Matrix_Size; d++)
                        {
                            X13(a,l)(i,j) -= MOs(a,l)(i,d) * t1(d,j);
                        }
                    }
                }
            }
        }
    }
}

























// Check
void CCSDT::Fill_W1()
{
    for (int i = 0; i < n_Electrons; i++)
    {
        for (int j = 0; j < n_Electrons; j++)
        {
            for (int k = 0; k < n_Electrons; k++)
            {
                for (int l = 0; l < n_Electrons; l++)
                {
                    W1(i,j)(k,l) = MOs(i,j)(k,l);
                    for (int c = n_Electrons; c < Matrix_Size; c++)
                    {
                        W1(i,j)(k,l) -= MOs(j,c)(k,l) * t1(c,i);
                        W1(i,j)(k,l) += MOs(i,c)(k,l) * t1(c,j);

                        for (int d = n_Electrons; d < Matrix_Size; d++)
                        {
                            W1(i,j)(k,l) += 0.5 * MOs(k,l)(c,d) * tau3(c,d)(i,j);
                        }
                    }
                }
            }
        }
    }
}

// Check
void CCSDT::Fill_W2()
{
    for (int a = n_Electrons; a < Matrix_Size; a++)
    {
        for (int i = 0; i < n_Electrons; i++)
        {
            for (int j = 0; j < n_Electrons; j++)
            {
                for (int k = 0; k < n_Electrons; k++)
                {
                    W2(a,k)(i,j) = MOs(i,j)(a,k);

                    for (int c = n_Electrons; c < Matrix_Size; c++)
                    {
                        W2(a,k)(i,j) += MOs(i,c)(a,k) * t1(c,j);
                        W2(a,k)(i,j) -= MOs(j,c)(a,k) * t1(c,i);

                        for (int d = n_Electrons; d < Matrix_Size; d++)
                        {
                            W2(a,k)(i,j) += 0.5 * MOs(a,k)(c,d) * tau3(c,d)(i,j);
                        }
                    }
                }
            }
        }
    }
}

// Check
void CCSDT::Fill_W3()
{
    for (int c = n_Electrons; c < Matrix_Size; c++)
    {
        for (int k = 0; k < n_Electrons; k++)
        {
            for (int l = 0; l < n_Electrons; l++)
            {
                for (int i = 0; i < n_Electrons; i++)
                {
                    W3(c,i)(k,l) = MOs(c,i)(k,l);

                    for (int d = n_Electrons; d < Matrix_Size; d++)
                    {
                        W3(c,i)(k,l) += MOs(c,d)(k,l) * t1(d,i);
                    }
                }
            }
        }
    }
}

// Check
void CCSDT::Fill_W4()
{
    for (int a = n_Electrons; a < Matrix_Size; a++)
    {
        for (int k = 0; k < n_Electrons; k++)
        {
            for (int i = 0; i <n_Electrons; i++)
            {
                for (int c = n_Electrons; c < Matrix_Size; c++)
                {
                    W4(a,k)(i,c) = MOs(a,k)(i,c);

                    for (int d = n_Electrons; d < Matrix_Size; d++)
                    {
                        W4(a,k)(i,c) += MOs(k,a)(c,d) * t1(d,i);
                    }

                    for (int l = 0; l < n_Electrons; l++)
                    {
                        W4(a,k)(i,c) -= t1(a,l) * W3(c,i)(k,l);

                        for (int d = n_Electrons; d < Matrix_Size; d++)
                        {
                            W4(a,k)(i,c) += 0.5 * MOs(k,l)(c,d) * t2(a,d)(i,l);
                        }
                    }
                }
            }
        }
    }
}

// Check
void CCSDT::Fill_F1()
{
    for (int k = 0; k < n_Electrons; k++)
    {
        for (int c = n_Electrons; c < Matrix_Size; c++)
        {
            D1(c,k) = fs(c,k);

            for (int l = 0; l < n_Electrons; l++)
            {
                for (int d = n_Electrons; d < Matrix_Size; d++)
                {
                    D1(c,k) += MOs(k,l)(c,d) * t1(d,l);
                }
            }
        }
    }
}

// Check
void CCSDT::Fill_F2()
{
    for (int i = 0; i < n_Electrons; i++)
    {
        for (int k = 0; k < n_Electrons; k++)
        {
            D2(k,i) = (1 - EqualFunc(k,i)) * fs(k,i);

            for (int c = n_Electrons; c < Matrix_Size; c++)
            {
                D2(k,i) += t1(c,i) * D1(c,k);

                for (int l = 0; l < n_Electrons; l++)
                {
                    D2(k,i) += MOs(i,c)(k,l) * t1(c,l);

                    for (int d = n_Electrons; d < Matrix_Size; d++)
                    {
                        D2(k,i) += 0.5 * MOs(k,l)(c,d) * t2(c,d)(i,l);
                    }
                }
            }
        }
    }
}

// Check
void CCSDT::Fill_F3()
{
    for (int a = n_Electrons; a < Matrix_Size; a++)
    {
        for (int c = n_Electrons; c < Matrix_Size; c++)
        {
            D3(a,c) = (1 - EqualFunc(a,c)) * fs(a,c);

            for (int k = 0; k < n_Electrons; k++)
            {
                D3(a,c) -= t1(a,k) * D1(c,k);

                for (int d = n_Electrons; d < Matrix_Size; d++)
                {
                    D3(a,c) += MOs(a,k)(c,d) * t1(d,k);

                    for (int l = 0; l < n_Electrons; l++)
                    {
                        D3(a,c) -= 0.5 * MOs(k,l)(c,d) * t2(a,d)(k,l);
                    }
                }
            }
        }
    }
}

// Check
void CCSDT::Fill_tau3()
{
    for (int a = n_Electrons; a < Matrix_Size; a++)
    {
        for (int b = n_Electrons; b < Matrix_Size; b++)
        {
            for (int i = 0; i < n_Electrons; i++)
            {
                for (int j = 0; j < n_Electrons; j++)
                {
                    tau3(a,b)(i,j) = t2(a,b)(i,j) + (t1(a,i) * t1(b,j)) - (t1(a,j) * t1(b,i));
                }
            }
        }
    }
}






