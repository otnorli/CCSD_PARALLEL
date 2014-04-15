#include "hartree_fock_solver.h"

Hartree_Fock_Solver::Hartree_Fock_Solver(int n_N, vec zz, mat rr, string B_s, int n_elec, bool pstf, int ran, int siz)
{
    n_Nuclei = n_N;
    Z = zz;
    R = rr;
    Basis_Set = B_s;
    n_Electrons = n_elec;
    print_stuff = pstf;
    rank = ran;
    size = siz;
}

void Hartree_Fock_Solver::Map_MPI()
{
    // Complete MPI here
    // Variable Fock_MPI scales like n^3, where n is number of Contracted GTOs
    // This is necasary because we want indexes (i) and (j) to be seperated, and also we need to swap k and m,
    // but we can remove one of the indexes due to row and col use later

    WORK_PER_NODE = zeros(size);
    Fock_MPI = (double**) malloc(size * sizeof(double *));

    int index_counter;
    int INDEX_CHECK;

    for (int K = 0 ; K < size; K++)
    {
        index_counter = 0;

        for (int i = 0; i < Matrix_Size; i++)
        {
            for (int k = 0; k < Matrix_Size; k++)
            {
                INDEX_CHECK = i+k;
                if (INDEX_CHECK % size == K)
                {
                    for (int j = 0; j < Matrix_Size; j++)
                    {
                        index_counter += 1;
                    }
                }
            }
        }

        WORK_PER_NODE(K) = index_counter;
        Fock_MPI[K] = (double*) malloc(index_counter * sizeof(double));
    }
}

double Hartree_Fock_Solver::Calc_Integrals_On_The_Fly(int orb1, int orb3, int orb2, int orb4)
{
    // This functions calculates our two electron integrals on the fly
    // Input is indexes, where the indexes orb1, orb2 etc range from 0 -> Number of Contracted GTOs

    int i,j,k,m;

    // Here we find which atom the index belongs to, this is needed in calculations of 2 electron integrals
    i = Calc_Which_Atom_We_Are_Dealing_With(orb1);
    j = Calc_Which_Atom_We_Are_Dealing_With(orb3);
    k = Calc_Which_Atom_We_Are_Dealing_With(orb2);
    m = Calc_Which_Atom_We_Are_Dealing_With(orb4);

    // Here we calculate the two electron integrals, we have already stored E_ij^t so they are quite fast.
    // Symmetry considerations must be applied elsewhere.

    int E_counter1, E_counter2; // These ensures we get the right E_ij^t
    int n,p,o,q;
    double temp = 0;
    E_counter1 = E_index(orb1,orb2);
    for (n=0; n<n_Basis(orb1); n++)
    {
        for (p=0; p<n_Basis(orb2); p++)
        {
            E_counter2 = E_index(orb3, orb4);
            for (o=0; o<n_Basis(orb3); o++)
            {
                for (q=0; q<n_Basis(orb4); q++)
                {
                    temp += c(orb1,n)*c(orb2,p)*c(orb3,o)*c(orb4,q)*
                            HartInt.Electron_Electron_Interaction_Single(orb1, orb3, orb2, orb4,
                                                                            i, j, k, m,
                                                                            n, o, p, q,
                                                                            E_counter1, E_counter2);


                    E_counter2 += 3;
                }
            }
            E_counter1 += 3;
        }
    }
    return temp;
}

int Hartree_Fock_Solver::Calc_Which_Atom_We_Are_Dealing_With(int orb1)
{
    int atom = -1;
    for (int i = 0; i < n_Nuclei; i++)
    {
        for (int j = 0; j < Number_Of_Orbitals(i); j++)
        {
            atom++;
            if (atom == orb1)
            {
                return i;
            }
        }
    }
}

mat Hartree_Fock_Solver::Return_Field_Q(int i, int j)
{
    return field_Q(i,j);
}

mat Hartree_Fock_Solver::ReturnC()
{
    return C;
}

mat Hartree_Fock_Solver::ReturnAlpha()
{
    return alpha;
}

int Hartree_Fock_Solver::ReturnMatrixSize()
{
    return Matrix_Size;
}

mat Hartree_Fock_Solver::ReturnPotenser()
{
    return Potenser;
}

vec Hartree_Fock_Solver::ReturnNumberOfOrbitals()
{
    return Number_Of_Orbitals;
}

vector<cube> Hartree_Fock_Solver::ReturnQ()
{
    return return_Q;
}

//vec Hartree_Fock_Solver::Return_Indexed_Q()
//{
    //return Stored_Indexed_Q;
//    return 0;
//}

mat Hartree_Fock_Solver::return_eigval_F()
{
    return eigenvalues_F;
}

mat Hartree_Fock_Solver::returnC()
{
    return C;
}

int Hartree_Fock_Solver::Get_Integral_Index(int a, int b, int i, int j)
{
    int ab, ij;

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

void Hartree_Fock_Solver::Normalize_small_c()
{
    int i = -1;
    //int angmom;
    double pi_i_3 = acos(-1.0)*acos(-1.0)*acos(-1.0);

    for (int atom_nr=0; atom_nr < n_Nuclei; atom_nr++)
    {
        for (int basis_nr=0; basis_nr < Number_Of_Orbitals(atom_nr); basis_nr++)
        {
            i += 1; // i går fra [0, matrix size), skriver det slik for enkelt å plukke ut hvem atom vi er på

            // s orbital normalization
            if (Potenser(i,0) == 0 && Potenser(i,1) == 0 && Potenser(i,2) == 0)
            {
                for (int j = 0; j < n_Basis(i); j++)
                {
                    c(i,j) *= pow(2*alpha(i,j)/acos(-1.0), 0.75);
                }
            }

            // p orbital normalization
            if (Potenser(i,0) == 1 && Potenser(i,1) == 0 && Potenser(i,2) == 0)
            {
                for (int j = 0; j < n_Basis(i); j++)
                {
                    c(i,j) *= pow(128*alpha(i,j)*alpha(i,j)*alpha(i,j)*alpha(i,j)*alpha(i,j)/pi_i_3, 0.25);
                }
            }

            if (Potenser(i,0) == 0 && Potenser(i,1) == 1 && Potenser(i,2) == 0)
            {
                for (int j = 0; j < n_Basis(i); j++)
                {
                    c(i,j) *= pow(128*alpha(i,j)*alpha(i,j)*alpha(i,j)*alpha(i,j)*alpha(i,j)/pi_i_3, 0.25);
                }
            }

            if (Potenser(i,0) == 0 && Potenser(i,1) == 0 && Potenser(i,2) == 1)
            {
                for (int j = 0; j < n_Basis(i); j++)
                {
                    c(i,j) *= pow(128*alpha(i,j)*alpha(i,j)*alpha(i,j)*alpha(i,j)*alpha(i,j)/pi_i_3, 0.25);
                }
            }

            // d orbital normalization
            if (Potenser(i,0) == 2 && Potenser(i,1) == 0 && Potenser(i,2) == 0)
            {
                for (int j = 0; j < n_Basis(i); j++)
                {
                    c(i,j) *= pow((2048 * alpha(i,j) * alpha(i,j) *
                                  alpha(i,j) * alpha(i,j) * alpha(i,j) *
                                   alpha(i,j) * alpha(i,j)) / (9*pi_i_3), 0.25);
                }
            }

            if (Potenser(i,0) == 0 && Potenser(i,1) == 2 && Potenser(i,2) == 0)
            {
                for (int j = 0; j < n_Basis(i); j++)
                {
                    c(i,j) *= pow((2048 * alpha(i,j) * alpha(i,j) *
                                   alpha(i,j) * alpha(i,j) * alpha(i,j) *
                                   alpha(i,j) * alpha(i,j)) / (9*pi_i_3), 0.25);
                }
            }

            if (Potenser(i,0) == 0 && Potenser(i,1) == 0 && Potenser(i,2) == 2)
            {
                for (int j = 0; j < n_Basis(i); j++)
                {
                    c(i,j) *= pow((2048 * alpha(i,j) * alpha(i,j) *
                                   alpha(i,j) * alpha(i,j) * alpha(i,j) *
                                   alpha(i,j) * alpha(i,j)) / (9*pi_i_3), 0.25);
                }
            }

            if (Potenser(i,0) == 1 && Potenser(i,1) == 1 && Potenser(i,2) == 0)
            {
                for (int j = 0; j < n_Basis(i); j++)
                {
                    c(i,j) *= pow((2048 * alpha(i,j) * alpha(i,j) *
                                   alpha(i,j) * alpha(i,j) * alpha(i,j) *
                                   alpha(i,j) * alpha(i,j)) / (pi_i_3), 0.25);
                }
            }

            if (Potenser(i,0) == 1 && Potenser(i,1) == 0 && Potenser(i,2) == 1)
            {
                for (int j = 0; j < n_Basis(i); j++)
                {
                    c(i,j) *= pow((2048 * alpha(i,j) * alpha(i,j) *
                                   alpha(i,j) * alpha(i,j) * alpha(i,j) *
                                   alpha(i,j) * alpha(i,j)) / (pi_i_3), 0.25);
                }
            }

            if (Potenser(i,0) == 0 && Potenser(i,1) == 1 && Potenser(i,2) == 1)
            {
                for (int j = 0; j < n_Basis(i); j++)
                {
                    c(i,j) *= pow((2048 * alpha(i,j) * alpha(i,j) *
                                   alpha(i,j) * alpha(i,j) * alpha(i,j) *
                                   alpha(i,j) * alpha(i,j)) / (pi_i_3), 0.25);
                }
            }

            // Angulærmoment = 3

            if (Potenser(i,0) == 3 && Potenser(i,1) == 0 && Potenser(i,2) == 0)
            {
                for (int j = 0; j < n_Basis(i); j++)
                {
                    //c(i,j) *= pow(2*alpha(i,j)/acos(-1.0), 0.75) * sqrt(64*alpha(i,j)*alpha(i,j)*alpha(i,j)/9);
                    c(i,j) *= pow((32768 * alpha(i,j) * alpha(i,j) *
                                  alpha(i,j) * alpha(i,j) * alpha(i,j) *
                                  alpha(i,j) * alpha(i,j) * alpha(i,j) *
                                  alpha(i,j)/pi_i_3), 0.25)
                            * sqrt(1/15.0);
                }
            }

            if (Potenser(i,0) == 0 && Potenser(i,1) == 3 && Potenser(i,2) == 0)
            {
                for (int j = 0; j < n_Basis(i); j++)
                {
                    c(i,j) *= pow(32768* alpha(i,j) * alpha(i,j) *
                                  alpha(i,j) * alpha(i,j) * alpha(i,j) *
                                  alpha(i,j) * alpha(i,j) * alpha(i,j) *
                                  alpha(i,j)/pi_i_3, 0.25) * sqrt(1/15.0);
                }
            }

            if (Potenser(i,0) == 0 && Potenser(i,1) == 0 && Potenser(i,2) == 3)
            {
                for (int j = 0; j < n_Basis(i); j++)
                {
                    c(i,j) *= pow(32768 * alpha(i,j) * alpha(i,j) *
                                  alpha(i,j) * alpha(i,j) * alpha(i,j) *
                                  alpha(i,j) * alpha(i,j) * alpha(i,j) *
                                  alpha(i,j)/pi_i_3, 0.25)* sqrt(1/15.0);
                }
            }

            if (Potenser(i,0) == 2 && Potenser(i,1) == 1 && Potenser(i,2) == 0)
            {
                for (int j = 0; j < n_Basis(i); j++)
                {
                    c(i,j) *= pow(32768* alpha(i,j) * alpha(i,j) *
                                  alpha(i,j) * alpha(i,j) * alpha(i,j) *
                                  alpha(i,j) * alpha(i,j) * alpha(i,j) *
                                  alpha(i,j)/pi_i_3, 0.25)* sqrt(1/3.0);
                }
            }

            if (Potenser(i,0) == 2 && Potenser(i,1) == 0 && Potenser(i,2) == 1)
            {
                for (int j = 0; j < n_Basis(i); j++)
                {
                    c(i,j) *= pow(32768* alpha(i,j) * alpha(i,j) *
                                  alpha(i,j) * alpha(i,j) * alpha(i,j) *
                                  alpha(i,j) * alpha(i,j) * alpha(i,j) *
                                  alpha(i,j)/pi_i_3, 0.25)* sqrt(1/3.0);
                }
            }

            if (Potenser(i,0) == 1 && Potenser(i,1) == 2 && Potenser(i,2) == 0)
            {
                for (int j = 0; j < n_Basis(i); j++)
                {
                    c(i,j) *= pow(32768* alpha(i,j) * alpha(i,j) *
                                  alpha(i,j) * alpha(i,j) * alpha(i,j) *
                                  alpha(i,j) * alpha(i,j) * alpha(i,j) *
                                  alpha(i,j)/pi_i_3, 0.25)* sqrt(1/3.0);
                }
            }

            if (Potenser(i,0) == 0 && Potenser(i,1) == 2 && Potenser(i,2) == 1)
            {
                for (int j = 0; j < n_Basis(i); j++)
                {
                    c(i,j) *= pow(32768* alpha(i,j) * alpha(i,j) *
                                  alpha(i,j) * alpha(i,j) * alpha(i,j) *
                                  alpha(i,j) * alpha(i,j) * alpha(i,j) *
                                  alpha(i,j)/pi_i_3, 0.25)* sqrt(1/3.0);
                }
            }

            if (Potenser(i,0) == 1 && Potenser(i,1) == 0 && Potenser(i,2) == 2)
            {
                for (int j = 0; j < n_Basis(i); j++)
                {
                    c(i,j) *= pow(32768* alpha(i,j) * alpha(i,j) *
                                  alpha(i,j) * alpha(i,j) * alpha(i,j) *
                                  alpha(i,j) * alpha(i,j) * alpha(i,j) *
                                  alpha(i,j)/pi_i_3, 0.25)* sqrt(1/3.0);
                }
            }

            if (Potenser(i,0) == 0 && Potenser(i,1) == 1 && Potenser(i,2) == 2)
            {
                for (int j = 0; j < n_Basis(i); j++)
                {
                    c(i,j) *= pow(32768* alpha(i,j) * alpha(i,j) *
                                  alpha(i,j) * alpha(i,j) * alpha(i,j) *
                                  alpha(i,j) * alpha(i,j) * alpha(i,j) *
                                  alpha(i,j)/pi_i_3, 0.25)* sqrt(1/3.0);
                }
            }

            if (Potenser(i,0) == 1 && Potenser(i,1) == 1 && Potenser(i,2) == 1)
            {
                for (int j = 0; j < n_Basis(i); j++)
                {
                    c(i,j) *= pow(32768 * alpha(i,j) * alpha(i,j) *
                                  alpha(i,j) * alpha(i,j) * alpha(i,j) *
                                  alpha(i,j) * alpha(i,j) * alpha(i,j) *
                                  alpha(i,j)/pi_i_3, 0.25);
                }
            }
        }
    }

}

double Hartree_Fock_Solver::Calc_Density()
{
    for (int i = 0; i < Matrix_Size; i++)
    {
        for (int j = 0; j < Matrix_Size; j++)
        {
            P(i,j) = 0.0;
            for (int k = 0; k < ((int) n_Electrons/2); k++)
            {
                P(i,j) += 2*C(i,k)*C(j,k);
            }
        }
    }
    return 0;
}

double Hartree_Fock_Solver::Calc_Energy()
{
    // Effective energy calculations
    Single_E_Energy = accu(EK % P);
    Two_E_Energy = 0.5*accu(Energy_Fock_Matrix % P) - 0.5*Single_E_Energy;
    return Single_E_Energy+Two_E_Energy;
}

double Hartree_Fock_Solver::Make_Fock_Matrix()
{
    F = EK;

    // This must run in parallel, since we only want to store a few of the integrals on each node
    // In fact we are as of now not using the symmetries in the Atomic Orbitals

    int index_counter = 0;
    int work;
    int INDEX_CHECK;
    //vec Vec1 = zeros(Matrix_Size);
    //vec Vec2 = zeros(Matrix_Size);

    for (int i=0; i<Matrix_Size; i++)
    {
        for (int k=0; k<Matrix_Size; k++)
        {
            INDEX_CHECK = i+k;

            if (INDEX_CHECK % size == rank)
            {
                for (int j=0; j< Matrix_Size; j++)
                {
                    // Use this for on the fly calculations, currently we store n^4 / M^2, where M = number processors
                    //for (int m = 0; m < Matrix_Size; m++)
                    //{
                        //Vec1(m) = Calc_Integrals_On_The_Fly(i,k,j,m);
                        //Vec2(m) = Calc_Integrals_On_The_Fly(i,k,m,j);
                    //}

                    Fock_MPI[rank][index_counter] = accu(P.col(k) % (field_Q(i,k).row(j) - 0.5*field_Q(i,k).col(j).t()));
                    //Fock_MPI[rank][index_counter] = accu(P.col(k) % ( Vec1.t() - 0.5 * Vec2.t())); // This is for calculations on the fly :-D slightly ineffective
                    index_counter += 1;
                }
            }
        }
    }

    for (int X = 0; X < size; X++)
    {
        work = WORK_PER_NODE(X);
        MPI_Bcast(Fock_MPI[X], work, MPI_DOUBLE, X, MPI_COMM_WORLD);
    }

    for (int X = 0; X < size; X++)
    {
        index_counter = 0;
        for (int i=0; i<Matrix_Size; i++)
        {
            for (int k=0; k<Matrix_Size; k++)
            {
                INDEX_CHECK = i+k;

                if (INDEX_CHECK % size == X)
                {
                    for (int j=0; j<Matrix_Size; j++)
                    {
                        F(i,j) += Fock_MPI[X][index_counter];
                        index_counter += 1;
                    }
                }
            }
        }
    }



    // Not parallel anymore
    Energy_Fock_Matrix = F;
    DIIS();
    F = V.t() * F * V;

    return 0;
}

double Hartree_Fock_Solver::get_Energy(double toler)
{
    Initialize_DIIS();
    // Gogo HF
    int hf_limit = 300;

    // Dummie variabler vi kommer til å trenge underveis
    int i,j,k;
    double Energy;
    double NucleiRepulsion;
    int orb1=-1, orb2, orb3, orb4, at1, at2, at3, at4, m,n,o,p,q;
    double sum;
    double tolerance = toler;

    // Fyller opp initialbetingelsene. Alt dette er input og vi fyller det inn med to klasser foreløpig, siden det er masse tall bare
    matrix_size_setter matset(Z, Basis_Set, n_Nuclei);
    Matrix_Size = matset.Set_Matrix_Size();

    if (print_stuff == true)
    {
        cout << "Number of orbitals for CCSD: " << 2*Matrix_Size << endl;
    }

    Fill_Alpha Fyll(n_Nuclei, Z, Basis_Set, Matrix_Size, matset.Return_Max_Bas_Func());
    alpha = Fyll.Fyll_Opp_Alpha();
    c = Fyll.Fyll_Opp_c();
    n_Basis = Fyll.Fyll_Opp_Nr_Basis_Functions();
    Number_Of_Orbitals = Fyll.Fyll_Opp_Antall_Orbitaler();
    Potenser = Fyll.Fyll_Opp_Potenser();
    //Stored_Indexed_Q = zeros(Get_Integral_Index(Matrix_Size-1, Matrix_Size-1, Matrix_Size-1, Matrix_Size-1) + 1); // Compact storage


    // For parallel implementation we do this:

    field_Q.set_size(Matrix_Size, Matrix_Size);
    for (int i = 0; i < Matrix_Size; i++)
    {
        for (int j = 0; j < Matrix_Size; j++)
        {
            if ((i+j)%size == rank)
            {
                field_Q(i,j) = zeros(Matrix_Size, Matrix_Size);
            }
        }
    }


    // Deklarerer matrisene med data vi kommer til å lagre under beregningene. Kunne lagret mindre men velger dette foreløpig
    EK = zeros(Matrix_Size, Matrix_Size); // Skal inneholde kinetisk energi og elektron-kjerne interaksjon
    O = zeros(Matrix_Size, Matrix_Size); // Skal inneholde overlap
    // C = zeros(Matrix_Size, Matrix_Size);  // Skal inneholde c1, c2, c3 osv til total bølgefunksjon
    C = zeros(Matrix_Size, Matrix_Size);
    F = zeros(Matrix_Size, Matrix_Size);
    P = zeros(Matrix_Size, Matrix_Size); // Skal inneholde electron density, n
    //int ijkl = get_Q_index(Matrix_Size, Matrix_Size, Matrix_Size, Matrix_Size);
    //Indexed_Q = zeros(ijkl);

    Normalize_small_c();
    Map_MPI();

    // Alle beregningene på integraler skjer i dette objektet, som er et globalt objekt
    //Hartree_Integrals HartInt();
    HartInt.Set_Input(alpha, Matrix_Size, c, Number_Of_Orbitals, Potenser, n_Basis, n_Nuclei, R, Z, n_Electrons);

    // Gjør litt pre-work, for å effektivisere beregningene
    HartInt.Fill_E_ij();

    if (print_stuff)
    {
        cout << "Initialization done!" << endl;
    }

    // Setter opp overlapps matrisen
    O = HartInt.Overlap_Matrix();

    // Setter opp elektron-kjerne interaksjoner og kinetisk energiledd
    EK = -0.5*HartInt.Kinetic_Energy()- HartInt.Nuclei_Electron_Interaction();

    // Setter opp kjerne-kjerne interaksjoner
    NucleiRepulsion = HartInt.Nuclei_Nuclei_Interaction();
    Nuclei_Rep_Energy = NucleiRepulsion;

    // Setter opp V
    V = zeros(Matrix_Size, Matrix_Size);
    eig_sym(eigenvalues_O, eigenvektors_O, O);
    V = eigenvektors_O * diagmat(1.0/sqrt(eigenvalues_O));
    if (print_stuff)
    {
        cout << "h1, S and NucleiRepulsion done!" << endl;
    }

    // Setter opp electron-electron interaksjoner
    //cube temp_cube = zeros(Matrix_Size, Matrix_Size, Matrix_Size);
    //for (int i = 0; i < Matrix_Size; i++)
    //{
    //    return_Q.push_back(temp_cube);
    //}


    // Calc Q not on the fly here, which can be used for small systems
    orb1 = -1;
    E_index = HartInt.get_E_index();
    double temp;

    for (i=0; i<n_Nuclei; i++) // Looper over alle atomene først
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

            //if (orb1 <= orb2) // symetri test 1
           // {

                orb3 = -1;
                for (j=0; j<n_Nuclei; j++)
                {
                for (at2 = 0; at2 < Number_Of_Orbitals(j); at2++)
                {
                orb3 += 1;

                if ((orb1+orb3)%size == rank)
                {
                if (orb1 <= orb3)
                {
                    orb4 = -1;

                    for (m=0; m<n_Nuclei; m++)
                    {
                    for (at4 = 0; at4 < Number_Of_Orbitals(m); at4++)
                    {
                    orb4 += 1;

                    //if (orb3 <= orb4) // symetri test 3
                    //{

                        //temp = 0;

                        //if (Stored_Indexed_Q(Get_Integral_Index(orb1, orb2, orb3, orb4)) == 0)
                        //{
                        temp = Calc_Integrals_On_The_Fly(orb1, orb3, orb2, orb4);




                            // * These symmetries must hold here and also cause our symmetries later in the MOs (integ in CCSD):
                            // *
                            // * kvantefysform !_! orb_nr form
                            // *
                            // *
                            // * <uv|pq>      !_! <13|24>
                            // * <vu|pq>      !_! <31|24>
                            // * <uv|qp>      !_! <13|42>
                            // * <vu|qp>      !_! <31|42>
                            // *
                            // * <qp|uv>      !_! <42|13>
                            // * <pq|vu>      !_! <24|31>
                            // * <qp|vu>      !_! <42|31>
                            // * <pq|uv>      !_! <24|13>
                             //*
                            // * In this program they are implemented with a flip of indicies to accomodate physics/chemical notation
                            // * which is (physics: <12|34>) = chemist: ([13|24])
                            // *
                             //* The flip is done to ensure we can use the stored E(...) in the same order as prior, such that E_counter is correct


                            //Stored_Indexed_Q(Get_Integral_Index(orb1, orb2, orb3, orb4)) = temp; // <-- PHYSICS NOTATION
                        //}



                        field_Q(orb1,orb3)(orb2,orb4) = temp; // <-- CHEMIST NOTATION, FLIPPED INDICIES ORB2 AND ORB3!
                        //field_Q(orb1,orb4)(orb2,orb3) = temp;

                        //field_Q(orb2,orb3)(orb1,orb4) = temp;
                        //field_Q(orb2,orb4)(orb1,orb3) = temp;

                        //field_Q(orb3,orb2)(orb4,orb1) = temp;
                        field_Q(orb3,orb1)(orb4,orb2) = temp;

                        //field_Q(orb4,orb2)(orb3,orb1) = temp;
                        //field_Q(orb4,orb1)(orb3,orb2) = temp;
                      // } // if test 3
                    }
                    }
                } // if test 2

                }
                }
            //} // if test
                }
        }
        }
    }
    }

    if (print_stuff)
    {
        cout << "Q done!" << endl;
    }


    // Kjører Hartree Fock
    bool run_HF = true;
    double hf_new=0, hf_old=0;
    steppp = 0;
    while (run_HF==true)
    {
        hf_old = hf_new;
        sum = Make_Fock_Matrix();
        //DIIS(); // DIIS() inserted into the Make_Fock_Matrix() function
        eig_sym(eigenvalues_F, eigenvektors_F, F);

        // Oppdaterer C
        C = V * eigenvektors_F;

        // Normaliserer C
        for (i = 0; i < Matrix_Size; i++)
        //for (i=0; i<((int) n_Electrons/2); i++)
        {
            C.col(i) = Normalize(C.col(i), O);
        }

        // Finn P, P er density matrix
        sum = Calc_Density();

        // Finner HF energien
        hf_new = 0;
        for (int i = 0; i < ((int) n_Electrons/2); i++)
        {
            hf_new += eigenvalues_F(i);
        }

        steppp += 1;
        if ((sqrt((hf_new - hf_old)*(hf_new - hf_old)) < tolerance) || (steppp > hf_limit))
        {
            run_HF = false;
        }

        Energy = Calc_Energy() + NucleiRepulsion;
        if (print_stuff == true)
        {
            cout << "Energi: " << Energy << " Steg: " << steppp << endl;
        }
    }

    Energy = Calc_Energy() + NucleiRepulsion;
    if (print_stuff)
    {
        cout << "Energy = " << Energy << endl;
        cout << "One-electron energy = " << Single_E_Energy << endl;
        cout << "Two-electron energy = " << Two_E_Energy << endl;
        cout << "Repulsion energy = " << Nuclei_Rep_Energy << endl;
        //cout << "C = " << endl << (C) << endl;
        cout << "Hartree Fock completed" << endl;
    }
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

void Hartree_Fock_Solver::Set_New_R(mat rr)
{
    R = rr;
}

void Hartree_Fock_Solver::DIIS()
{
    double temp_trace;
    // Method to speed up SCF, the iteration part
    // Algorithm based on: 
    // 1) http://onlinelibrary.wiley.com/doi/10.1002/jcc.540030413/pdf
    // 2) http://ac.els-cdn.com/0009261480803964/1-s2.0-0009261480803964-main.pdf?_tid=eee36e34-8dae-11e3-
    // a80b-00000aacb35e&acdnat=1391527002_8853c38053bb3713b2080b7d3c9f1c54

    // Articles refered to somewhat during this function so it should be easy to understand. 
    // Articles by same person and very similar content.

    // Only start after the first iteration, since P is needed and not calculated before our secound round in this function
    if (steppp > 0)
    {
        // Calculate error, this term is 0 when we have converged
        delta_p = F*P*O - O*P*F; // Error, from equation 4 in article (1)

        // Add the new matrixes
        Stored_Error.push_back(delta_p);
        Stored_F.push_back(F);

        if (steppp > number_elements_DIIS)
        {
            // This for loop ensures we only store the pre defined number of elements (this easy algo is actually predefined +1)
            for (int i = 0; i < number_elements_DIIS; i++)
            {
                Stored_Error.at(i) = Stored_Error.at(i+1);
                Stored_F.at(i) = Stored_F.at(i+1);
            }

            // Delete terms we dont need
            Stored_Error.resize(number_elements_DIIS);
            Stored_F.resize(number_elements_DIIS);

            // Make B matrix, which is <i|j>
            mat mat1, mat2;
            for (int i = 0; i < number_elements_DIIS; i++)
            {
                for (int j = 0; j < number_elements_DIIS; j++)
                {
                    mat1 = Stored_Error.at(i);
                    mat2 = Stored_Error.at(j);
                    temp_trace = trace(mat1.t() * mat2);
                    DIIS_B(i,j) = temp_trace;
                }
            }

            // Solve this DIIS_B * DIIS_C = DIIS_Z matrix problem with easy armadillo solve() function
            DIIS_c = solve(DIIS_B, DIIS_Z);

            // Construct F based on a pre defined past number of terms, one term outside loop to make sure no += problems with un-defined sized matrix
            F = DIIS_c.at(0) * Stored_F.at(0);
            for (int i = 1; i < number_elements_DIIS; i++)
            {
                F += DIIS_c.at(i) * Stored_F.at(i);
            }
        }
    }
}

void Hartree_Fock_Solver::Initialize_DIIS()
{
    // Initialize DIIS method, define some stuff here for speedup
    // Please see function DIIS() for info on the method
    number_elements_DIIS = 3; // This term is pre defined, but i choose 3 and not to take it as input (<-- !)
    DIIS_B = zeros(number_elements_DIIS+1, number_elements_DIIS+1);
    DIIS_Z = zeros(number_elements_DIIS+1);
    DIIS_Z(number_elements_DIIS) = -1.0;
    for (int i = 0; i < number_elements_DIIS; i++)
    {
        DIIS_B(i,number_elements_DIIS) = -1.0;
        DIIS_B(number_elements_DIIS,i) = -1.0;
    }
}

void Hartree_Fock_Solver::Delete_Everything()
{
    // Clear everything, this is kinda crazeh
    //Stored_Indexed_Q.clear();
    Stored_Error.clear();
    Stored_F.clear();
    EK.clear();
    O.clear();
    DIIS_B.clear();
    DIIS_Z.clear();
    DIIS_c.clear();
    R.clear();
    n_Basis.clear();
    c.clear();
    Number_Of_Orbitals.clear();
    D.clear();
    V.clear();
    F.clear();
    eigenvalues_F.clear();
    eigenvalues_O.clear();
    eigenvektors_O.clear();
    eigenvektors_F.clear();
    P.clear();
    C.clear();
    V.clear();
    alpha.clear();
    Potenser.clear();
    HartInt.Delete_Everything();
}

void Hartree_Fock_Solver::Unrestricted_Fock_Matrix()
{/*
    // Unfinished implementation of unrestricted HF
    F_up = EK;
    F_down = EK;
    for (int i=0; i<Matrix_Size; i++)
    {
        for (int j=0; j<Matrix_Size; j++)
        {
            for (int k=0; k<Matrix_Size; k++)
            {
                for (int m=0; m<Matrix_Size; m++)
                {
                    F_up(i,j) += P_up(m,k)*(Stored_Indexed_Q(Get_Integral_Index(i,j,k,m)) -
                                            Stored_Indexed_Q(Get_Integral_Index(i,m,k,j)))
                            + P_down(m,k) * Stored_Indexed_Q(Get_Integral_Index(i,j,k,m));
                    F_down(i,j) += P_down(m,k)*(Stored_Indexed_Q(Get_Integral_Index(i,j,k,m)) -
                                            Stored_Indexed_Q(Get_Integral_Index(i,m,k,j)))
                            + P_up(m,k) * Stored_Indexed_Q(Get_Integral_Index(i,j,k,m));
                }
            }
        }
    }
    F_up = V.t() * F_up * V;
    F_down = V.t() * F_down * V;
    */
}

void Hartree_Fock_Solver::Unrestricted_P_Matrix()
{
    // Unfinished implementation of unrestricted HF
}

void Hartree_Fock_Solver::Unrestricted_Energy()
{
    // Unfinished implementation of unrestricted HF
}
