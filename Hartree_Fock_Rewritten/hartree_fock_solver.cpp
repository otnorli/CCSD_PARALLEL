#include "hartree_fock_solver.h"

Hartree_Fock_Solver::Hartree_Fock_Solver(int n_N, vec zz, mat rr, string B_s)
{
    n_Nuclei = n_N;
    Z = zz;
    R = rr;
    Basis_Set = B_s;
}

double Hartree_Fock_Solver::get_Energy(int n_S)
{
    // Dummie variabler vi kommer til å trenge underveis
    int i,j,k,t;
    int E_counter1 = 0;
    int E_counter2 = 0;
    int n_Steps = n_S;
    double Energy;
    double NucleiRepulsion;
    int orb1=-1, orb2, orb3, orb4, at1, at2, at3, at4, m,n,o,p,q;
    double sum;

    n_Electrons = 0;
    for (int i=0; i<n_Nuclei; i++)
    {
        n_Electrons += Z(i); // Antallet elektroner
    }

    // Fyller opp initialbetingelsene. Alt dette er input og vi fyller det inn med to klasser foreløpig, siden det er masse tall bare
    matrix_size_setter matset(Z, Basis_Set, n_Nuclei);
    Matrix_Size = matset.Set_Matrix_Size();
    Fill_Alpha Fyll(n_Nuclei, Z, Basis_Set, Matrix_Size);
    alpha = Fyll.Fyll_Opp_Alpha();
    c = Fyll.Fyll_Opp_c();
    n_Basis = Fyll.Fyll_Opp_Nr_Basis_Functions();
    Number_Of_Orbitals = Fyll.Fyll_Opp_Antall_Orbitaler();
    Potenser = Fyll.Fyll_Opp_Potenser();

    // Deklarerer matrisene med data vi kommer til å lagre under beregningene. Kunne lagret mindre men velger dette foreløpig
    EK = zeros(Matrix_Size, Matrix_Size); // Skal inneholde kinetisk energi og elektron-kjerne interaksjon
    O = zeros(Matrix_Size, Matrix_Size); // Skal inneholde overlap
    C = zeros(Matrix_Size);  // Skal inneholde c1, c2, c3 osv til total bølgefunksjon
    P = zeros(Matrix_Size, Matrix_Size); // Skal inneholde electron density, n
    double ****Q; // Q skal holde elektron-elektron interaksjon leddene
    Q = new double***[Matrix_Size];

    for (i=0; i<Matrix_Size; i++)
    {
        Q[i] = new double**[Matrix_Size];
        for (j=0; j<Matrix_Size; j++)
        {
            Q[i][j] = new double*[Matrix_Size];
            for (k=0; k<Matrix_Size; k++)
            {
               Q[i][j][k] = new double[Matrix_Size];
            }
        }
    }

    // Initialiserer Q
    orb1 = -1;
    for (i=0; i<n_Nuclei; i++)
    {
    for (at1=0; at1<Number_Of_Orbitals(i); at1++)
    {
    orb1+=1;
    orb3 = -1;
    for (k=0; k<n_Nuclei; k++)
    {
    for (at3=0; at3 < Number_Of_Orbitals(k); at3++)
    {
        orb3+=1;
        orb2 = -1;
        for (j=0; j<n_Nuclei; j++)
        {
        for (at2 = 0; at2 < Number_Of_Orbitals(j); at2++)
        {
                orb2+=1;
                orb4=-1;
                for (m=0; m<n_Nuclei; m++)
                {
                    for (at4=0; at4<Number_Of_Orbitals(m); at4++)
                    {
                        orb4+=1;
                        Q[orb1][orb3][orb2][orb4] = 0;
                    }
                }
        }
        }

    }
    }
    }
    }

    // Alle beregningene skjer i dette objektet
    Hartree_Integrals HartInt(n_Nuclei, n_Basis, R, Z, alpha, Matrix_Size, c, Number_Of_Orbitals, Potenser, n_Electrons);

    // Gjør litt pre-work, for å effektivisere beregningene
    HartInt.Fill_E_ij();
    cout << "Initialization done!" << endl;

    // Setter opp overlapps matrisen
    O = HartInt.Overlap_Matrix();
    cout << "O done!" << endl;

    // Setter opp elektron-kjerne interaksjoner og kinetisk energiledd
    EK = -0.5*HartInt.Kinetic_Energy()- HartInt.Nuclei_Electron_Interaction();
    cout << "EK done!" << endl;
    cout << EK;

    // Setter opp kjerne-kjerne interaksjoner
    NucleiRepulsion = HartInt.Nuclei_Nuclei_Interaction();
    cout << "NucleiRepulsion done!" << endl;

    // Setter opp V
    eig_sym(eigenvalues_O, eigenvektors_O, O);
    V = eigenvektors_O * diagmat(1.0/sqrt(eigenvalues_O));




    orb1 = -1;
    for (i=0; i<n_Nuclei; i++) // Looper over alle atomene først dawg
    {
    for (at1 = 0; at1 < Number_Of_Orbitals(i); at1++)
    {
        orb1 += 1;
        orb2 = -1;
        for (k=0; k<n_Nuclei; k++)
        {
        for (at3 = 0; at3 < Number_Of_Orbitals(k); at3++)
        {
            orb2 += 1;

            for (n=0; n<n_Basis(orb1); n++)
            {
            for (p=0; p<n_Basis(orb2); p++)
            {
                E_counter2 = 0;
                orb3 = -1;

                for (j=0; j<n_Nuclei; j++)
                {
                for (at2 = 0; at2 < Number_Of_Orbitals(j); at2++)
                {
                    orb3 += 1;
                    orb4 = -1;

                for (m=0; m<n_Nuclei; m++)
                {
                for (at4 = 0; at4 < Number_Of_Orbitals(m); at4++)
                {
                    orb4 += 1;
                    sum = 0;

                    // Må loope over alle mulige comboer av basiser

                    for (o=0; o<n_Basis(orb3); o++)
                    {
                    for (q=0; q<n_Basis(orb4); q++)
                    {
                        Q[orb1][orb3][orb2][orb4] += c(orb1,n)*c(orb2,p)*c(orb3,o)*c(orb4,q)*
                             HartInt.Electron_Electron_Interaction_Single(orb1, orb3, orb2, orb4,
                                                                             i, j, k, m,
                                                                             n, o, p, q,
                                                                             E_counter1, E_counter2);
                        E_counter2 += 3;
                    }}
                }}
                }}
                E_counter1 += 3;
            }
            }
        }
        }
    }
    }
    cout << "Q done!" << endl;

    for (t=0; t<n_Steps; t++)
    {
        // Setter opp D for hvert steg. D = J-K
        D = zeros(Matrix_Size, Matrix_Size);
        for (i=0; i<Matrix_Size; i++)
        {
            for (j=0; j<Matrix_Size; j++)
            {
                for (k=0; k<Matrix_Size; k++)
                {
                    for (m=0; m<Matrix_Size; m++)
                    {
                        D(i,j) += (Q[i][k][j][m]-0.5*Q[i][k][m][j])*P(m,k);
                    }
                }
            }
        }


        // Setter opp F hvert tidssteg, F er fock operator, F = h1 + D
        F = EK+D;
        F = V.t() * F * V;
        eig_sym(eigenvalues_F, eigenvektors_F, F);

        // Oppdaterer C
        C = V * eigenvektors_F.cols(0, n_Electrons/2-1);

        // Normaliserer C
        /*for (i=0; i<(n_Electrons/2-1); i++)
        {
            C.col(i) = Normalize(C.col(i), O);
        }*/
        C = Normalize(C, O);

        // Finn P, P er density matrix
        P = 2*C*C.t();

        // Finner energien
        Energy = 0.0;

        for (i=0; i<Matrix_Size; i++)
        {
            for (j=0; j<Matrix_Size; j++)
            {
                Energy += EK(i,j)*P(i,j);
            }
        }

        for (i=0; i<Matrix_Size; i++)
        {
            for (j=0; j<Matrix_Size; j++)
            {
                for (k=0; k<Matrix_Size; k++)
                {
                    for (m=0; m<Matrix_Size; m++)
                    {
                        Energy += 0.5*P(i,j)*P(m,k)*(Q[i][k][j][m]-0.5*Q[i][k][m][j]);
                    }
                }
            }
        }

        Energy += NucleiRepulsion;

        cout << "Energi: " << Energy << " Steg: " << t+1 << endl;
    }

    cout << "C = " << endl << abs(C) << endl;

    // Sletter variabler siden vi er ferdige
    for (i=0; i<Matrix_Size; i++)
    {
        for (j=0; j<Matrix_Size; j++)
        {
           for (k=0; k<Matrix_Size; k++)
           {
                delete [] Q[i][j][k];
           }
           delete [] Q[i][j];
        }
        delete [] Q[i];
    }
    delete [] Q;
    return Energy;
}

vec Hartree_Fock_Solver::Normalize(vec Vektorn, mat Matrisen)
{
    // Denne funksjonen normaliserer c1, c2, c3... i Vektorn mhp egenvektorene i Matrisen
    int i,j;
    double normaliserings_faktor = 0.0;

    for (i=0; i<Matrix_Size; i++)
    {
        for (j=0; j<Matrix_Size; j++)
        {
           normaliserings_faktor += Vektorn(i) * Matrisen(i,j) * Vektorn(j);
        }
    }

    return Vektorn/sqrt(normaliserings_faktor);
}
