#include "hartree_integrals.h"

Hartree_Integrals::Hartree_Integrals()
{

}

void Hartree_Integrals::Set_Input(mat alp, int n_tot_bf, mat cc, vec n_O, mat pot, vec n_B, int n_N, mat pos, vec zZz, int n_E)
{
    alpha = alp;
    Potenser = pot;
    Matrix_Size = n_tot_bf;
    c = cc;
    n_Orbitals = n_O;
    n_Basis = n_B;
    R = pos;
    n_Nuclei = n_N;
    Z = zZz;
    n_Electrons = n_E;
}

void Hartree_Integrals::Set_R_ijk(double p, int t, int u, int v, rowvec R1, rowvec R2)
{
    int t_max;
    int nn,i,j,k;

    int tt = t;
    int uu = u;
    int vv = v;

    t_max = t+u+v;

    double Boys_arg;
    rowvec Rcp(3);

    Rcp = R1-R2;
    Boys_arg = p*dot(Rcp, Rcp);
    Boys_arg = Boys(Boys_arg, 0); // This is for loading the boys array, to avoid calculating same numbers over and over

    for (nn=0; nn<(t_max+1);nn++)
    {
        R_ijk.at(nn)(0,0,0) = pow(-2*p, nn) * F_Boys(nn);
    }

    for (i=0; i<tt; i++)
    {
        for (nn=0; nn<(t_max-i); nn++)
        {
            R_ijk.at(nn)(i+1,0,0) = Rcp(0) * R_ijk.at(nn+1)(i,0,0);

            if (i > 0)
            {
                R_ijk.at(nn)(i+1,0,0) += i * R_ijk.at(nn+1)(i-1,0,0);
            }
        }
    }

    for (i=0; i<tt+1; i++)
    {
        for (j=0; j<uu; j++)
        {
            for (nn=0; nn<(t_max-i-j); nn++)
            {
                R_ijk.at(nn)(i, j+1, 0) = Rcp(1) * R_ijk.at(nn+1)(i,j,0);

                if (j>0)
                {
                    R_ijk.at(nn)(i, j+1, 0) += j * R_ijk.at(nn+1)(i,j-1,0);
                }
            }
        }
    }

    for (i=0; i<tt+1; i++)
    {
        for (j=0; j<uu+1; j++)
        {
            for (k=0; k<vv; k++)
            {
                for (nn=0; nn<(t_max-i-j-k); nn++)
                {
                    R_ijk.at(nn)(i, j, k+1) = Rcp(2) * R_ijk.at(nn+1)(i,j,k);

                    if (k>0)
                    {
                        R_ijk.at(nn)(i, j, k+1) += k*R_ijk.at(nn+1)(i,j,k-1);
                    }
                }
            }
        }
    }
}

void Hartree_Integrals::Delete_Everything()
{
    n_Basis.clear();
    n_Orbitals.clear();
    Z.clear();
    R.clear();
    c.clear();
    alpha.clear();
    Potenser.clear();
    Ek.clear();
    O.clear();
    E_ij.clear();
    F_Boys.clear();
    R_ijk.clear();
}

void Hartree_Integrals::Fill_E_ij()
{
    // This is long function, but only called one time during the program.
    // This is the genious of E_counter during the program. It finds the correct E_ij always, even in 2-e calculations
    int max_potens = max(max(Potenser));
    int potensgrenser = (4*max_potens+2) * 3 + 1;
    cube temp_cubee = zeros(4*max_potens+3, 4*max_potens+3, 4*max_potens+3);

    for (int nnn=0; nnn < potensgrenser; nnn++)
    {
        R_ijk.push_back(temp_cubee);
    }

    int i,j,k,m,p,q,dim,t;
    int I,J;

    double argz, Xp, Xpb, Xpa, Xab, R1, R2;
    double alpha1, alpha2, pterm;

    int A, B;

    int index_counter1 = -1, index_counter2;
    
    E_index = zeros(Matrix_Size, Matrix_Size);
    int E_counter = 0;
    index_counter1 = -1;
    for (i=0; i<n_Nuclei; i++)
    {
        for (p=0; p<n_Orbitals(i); p++)
        {
            index_counter1 += 1;
            index_counter2 = -1;
            for (j=0; j<n_Nuclei; j++)
            {
                for (q=0; q<n_Orbitals(j); q++)
                {
                    index_counter2 += 1;

                    //if (index_counter1 <= index_counter2)
                        // Not gonna need all, because of symetri.
                        // This if test make sure we only calculate half the matrix + the diagonal
                        // Which is all we will need
                    //{

                    E_index(index_counter1, index_counter2) = E_counter;

                    for (k=0; k<n_Basis(index_counter1); k++)
                    {
                        for (m=0; m<n_Basis(index_counter2); m++)
                        {
                            alpha1 = alpha(index_counter1, k);
                            alpha2 = alpha(index_counter2, m);
                            argz = alpha1+alpha2;
                            pterm = 0.5/argz;

                            for (dim=0; dim<3; dim++)
                            {
                                A = Potenser(index_counter1, dim);
                                B = Potenser(index_counter2, dim);

                                R1 = R(i,dim);
                                R2 = R(j,dim);

                                cube Emat = zeros(A+3, B+3, (A+B+6)); //(A+3)*2+1);

                                // Kalkulasjonen for en bestemt alpha1, alpha2, X1, X2
                                Xp = (alpha1*R1 + alpha2*R2)/argz;
                                Xpb = Xp - R2;
                                Xpa = Xp - R1;
                                Xab = R1-R2;

                                Emat(0,0,0) = exp(-alpha1*alpha2*Xab*Xab/argz);

                                for (I=0; I< (A+2); I++)
                                {
                                    // t = 0
                                    Emat(I+1,0,0) = Xpa*Emat(I,0,0) + Emat(I,0,1);
                                    for (t=1; t<=(I+1); t++)
                                    {
                                        Emat(I+1,0,t) = Emat(I,0,t-1)*pterm + Xpa*Emat(I,0,t) + (t+1)*Emat(I,0,t+1);
                                    }
                                }

                                for (I=0; I<= A+2; I++)
                                {
                                    for (J=0; J<B+2; J++)
                                    {
                                        Emat(I, J+1, 0) = Xpb * Emat(I,J,0) + Emat(I,J,1);

                                        for (t=1; t < (I+J+2); t++)
                                        {
                                            Emat(I,J+1,t) = Emat(I,J,t-1) *pterm + Xpb*Emat(I,J,t) + (t+1)*Emat(I,J,t+1);
                                        }
                                    }
                                }

                                E_ij.push_back(Emat);
                                E_counter += 1;
                            }
                        }
                    }
                    //}
                }
            }
        }
    }
}

mat Hartree_Integrals::get_E_index()
{
    return E_index;
}

int Hartree_Integrals::Double_Factorial(int TALL)
{
    if (TALL == 0)
    {
        return 1;
    }

    else
    {

        int i=1;
        int output=1;

        int limit=1;
        for (i=2; i<=TALL; i++)
        {
            limit *= i;
        }

        for (i=2; i<=limit; i++)
        {
            output *= i;
        }

        return output;
    }
}

void Hartree_Integrals::Set_Boys_Start(int N)
{
    Boys_N = N;

    // Double factorial
    Boys_Start = 1;
    double i = 1;

    while(i < 2*N - 1){
        i += 2;
        Boys_Start *= i;
    }

    //Boys_Start *= sqrt(M_PI)/pow(2,N+1);
}

int Hartree_Integrals::Factorial(int N)
{
    int value=1;
    for (int i = 2; i < (N+1); i++)
    {
        value *= i;
    }
    return value;
}

mat Hartree_Integrals::Overlap_Matrix()
{
    int k,m;
    int E_counter = 0;
    double sum;
    mat Overlap_Matrise = zeros(Matrix_Size, Matrix_Size);

    for (int index_counter1=0; index_counter1 < Matrix_Size; index_counter1++)
    {
        for (int index_counter2=index_counter1; index_counter2 < Matrix_Size; index_counter2++)
        //for (int index_counter2=0; index_counter2 < Matrix_Size; index_counter2++)
        {
            E_counter = E_index(index_counter1, index_counter2);
            sum = 0;
            for (k=0; k<n_Basis(index_counter1); k++)
            {
                for (m=0; m<n_Basis(index_counter2); m++)
                {
                    sum += c(index_counter1, k) * c(index_counter2, m)
                           * Overlap_Integral_Single(index_counter1, index_counter2, k, m, E_counter);
                    E_counter += 3;
                }
            }

            Overlap_Matrise(index_counter1, index_counter2) = sum;
            Overlap_Matrise(index_counter2, index_counter1) = sum;
        }
    }
    return Overlap_Matrise;
}

double Hartree_Integrals::Overlap_Integral_Single(int ind1, int ind2, int bas1, int bas2, int E_counter)
{
    // Returnerer overlappen mellom to basisfunksjoner
    int A=0, B=0, C=0, D=0, E=0, F=0;
    double value;
    double x1, x2;

    x1 = alpha(ind1,bas1);
    x2 = alpha(ind2,bas2);

    A = Potenser(ind1,0);
    B = Potenser(ind1,1);
    C = Potenser(ind1,2);
    D = Potenser(ind2,0);
    E = Potenser(ind2,1);
    F = Potenser(ind2,2);

    // Prefaktoren
    //value = pow(acos(-1.0)/(x1+x2), 1.5);
    value = pow(acos(-1.0)/(x1+x2), 1.5);

    // E_ij termene i x,y,z retning hver for seg
    value *= E_ij.at(E_counter)(A,D,0);
    value *= E_ij.at(E_counter+1)(B,E,0);
    value *= E_ij.at(E_counter+2)(C,F,0);

    return value;
}

mat Hartree_Integrals::Kinetic_Energy()
{
    // Finner kinetic energy termen
    // |psi> er orbitalen vår
    //
    // Denne funksjonen løser <psi|Delta^2|psi>
    //
    // Kombinert med funksjonen NucleiElectronInteraction() løser de to
    // <psi|h1|psi>

    mat Ek = zeros(Matrix_Size, Matrix_Size);
    int k,m;
    double sum;
    int E_counter=0;

    for (int index_counter1=0; index_counter1 < Matrix_Size; index_counter1++)
    {
        for (int index_counter2=index_counter1; index_counter2 < Matrix_Size; index_counter2++)
        {
            sum = 0;
            E_counter = E_index(index_counter1, index_counter2);
            // Må nå beregne <i,j|delta^2|i,j> for alle kombinasjoner av basisfunskjoner for de to orbitalene q og p
            for (k=0; k< n_Basis(index_counter1); k++)
            {
                for (m=0; m<n_Basis(index_counter2); m++)
                {
                    // Må regne summen over i og j <i,j|i,j> for alle basisene som utrykker de to orbitalene vi nå ser på
                    sum += c(index_counter1,k) * c(index_counter2,m) *
                            Kinetic_Energy_Single(index_counter1, index_counter2, k, m, E_counter);
                    E_counter += 3;
                }
            }
            Ek(index_counter1, index_counter2) = sum;
            Ek(index_counter2, index_counter1) = sum;
        }
    }

    return Ek;
}

double Hartree_Integrals::Kinetic_Energy_Single(int ind1, int ind2, int bas1, int bas2, int E_counter)
{
    // Her regner vi ut <ABC|(d/dr)^2|DEF>

    double x1,x2;
    x1 = alpha(ind1, bas1);
    x2 = alpha(ind2, bas2);

    double value=0;
    int A, B, C, D, E, F;

    /*
    A = Set_A_B(bas1);
    B = Set_A_B(bas1-A);
    C = bas1-A-B;

    D = Set_A_B(bas2);
    E = Set_A_B(bas2-D);
    F = bas2-D-E;
    */

    A = Potenser(ind1,0);
    D = Potenser(ind2,0);
    B = Potenser(ind1,1);
    E = Potenser(ind2,1);
    C = Potenser(ind1,2);
    F = Potenser(ind2,2);

    double j_current; // Dette er T_ij leddet, brukes i x y og z

    // Kjører en generell formel for beregning av integralet

    // x_retningz
    j_current = 4*x2*x2*E_ij.at(E_counter)(A,D+2,0) - 2*x2*(2*D+1)*E_ij.at(E_counter)(A,D,0);
    if ((D>1))
    {
        j_current += D*(D-1)*E_ij.at(E_counter)(A,D-2,0);
    }
    value += j_current * E_ij.at(E_counter+1)(B,E,0) * E_ij.at(E_counter+2)(C,F,0);

    // y retningz
    j_current = 4*x2*x2*E_ij.at(E_counter+1)(B, E+2, 0) - 2*x2*(2*E+1)*E_ij.at(E_counter+1)(B, E, 0);
    if (E>1)
    {
        j_current += E*(E-1)*E_ij.at(E_counter+1)(B, E-2, 0);
    }
    value += j_current * E_ij.at(E_counter)(A,D,0)*E_ij.at(E_counter+2)(C,F,0);

    // z retningz
    j_current = 4*x2*x2*E_ij.at(E_counter+2)(C,F+2,0) - 2*x2*(2*F+1)*E_ij.at(E_counter+2)(C,F,0);
    if (F>1)
    {
        j_current += F*(F-1)*E_ij.at(E_counter+2)(C,F-2,0);
    }
    value += j_current * E_ij.at(E_counter)(A,D,0) * E_ij.at(E_counter+1)(B,E,0);

    // Const term
    value *= pow(acos(-1.0)/(x1+x2), 1.5);
    return value;
}

mat Hartree_Integrals::Nuclei_Electron_Interaction()
{
    int i,j,k,m;
    int q,p;
    mat Ne_E_Interaction = zeros(Matrix_Size, Matrix_Size);
    int index_counter1 = -1, index_counter2 = -1;
    double sum;
    int E_counter = 0;

    for (i=0; i<n_Nuclei; i++) // Looper over det første atomet
    {
        for (p=0; p<n_Orbitals(i); p++) // 1 og 1 orbital fra første atom
        {
            index_counter1 += 1;
            index_counter2 = -1; // Regner 1 og 1 rad

            for (j=0; j<n_Nuclei;j++) // Looper over det andre atomet
            {
                for (q=0; q < n_Orbitals(j); q++) // Looper over 1 og 1 orbital fra andre atom
                {
                    index_counter2 += 1;
                    sum = 0;

                    if (index_counter1 <= index_counter2)
                    {

                        E_counter = E_index(index_counter1, index_counter2);

                        // Må nå beregne for alle kombinasjoner av basisfunskjoner for de to orbitalene q og p
                        for (k=0; k< n_Basis(index_counter1); k++)
                        {
                            for (m=0; m<n_Basis(index_counter2); m++)
                            {
                                sum += c(index_counter1, k)*c(index_counter2, m)*
                                        Nuclei_Electron_Interaction_Single(index_counter1, index_counter2, i, j, k, m, E_counter);
                                E_counter += 3;
                            }
                        }

                        Ne_E_Interaction(index_counter1, index_counter2) = sum;
                        Ne_E_Interaction(index_counter2, index_counter1) = sum;
                    }
                }
            }
        }
    }
    return Ne_E_Interaction;
}

double Hartree_Integrals::Nuclei_Electron_Interaction_Single(int ind1, int ind2, int a1, int a2, int bas1, int bas2, int E_counter)
{
    double x1,x2;
    x1 = alpha(ind1, bas1);
    x2 = alpha(ind2, bas2);

    int A, B, C, D, E, F;
    double value=0;
    int t,u,v, nn;
    double temp_sum;
    rowvec Rc;
    double p;

    A = Potenser(ind1,0);
    B = Potenser(ind1,1);
    C = Potenser(ind1,2);
    D = Potenser(ind2,0);
    E = Potenser(ind2,1);
    F = Potenser(ind2,2);

    p = x1+x2;
    Rc = (x1*R.row(a1)+x2*R.row(a2))/p;

    for (t=0; t<(A+D+1); t++)
    {
        for (u=0; u<(B+E+1); u++)
        {
            for (v=0; v<(C+F+1); v++)
            {
                temp_sum = 0;
                for (nn=0; nn<n_Nuclei;nn++)
                {
                    temp_sum += Z(nn)*Nuclei_Electron_Interaction_Single_1d(p,t,u,v,Rc,R.row(nn),0);
                }

                value += E_ij.at(E_counter)(A, D, t) * E_ij.at(E_counter+1)(B, E, u) * E_ij.at(E_counter+2)(C, F, v) * temp_sum;
            }
        }
    }

    value *= 2*acos(-1.0)/p;
    return value;
}

double Hartree_Integrals::Nuclei_Electron_Interaction_Single_1d(double p, int t, int u, int v, rowvec R1, rowvec R2, double n)
{
    int t_max;
    int nn,i,j,k;

    int tt = t;
    int uu = u;
    int vv = v;

    t_max = t+u+v;

    //cube Rinit = zeros<cube>(t_max+1, t_max+1, t_max+1);
    cube Rinit = zeros<cube>(tt+1, uu+1, vv+1);
    vector<cube> Rcur;
    Rcur.clear();

    double Boys_arg;
    rowvec Rcp(3);

    Rcp = R1-R2;
    Boys_arg = p*dot(Rcp, Rcp);

    Rcur.clear();
    for (nn=0; nn<(t_max+1);nn++)
    {
        Rinit(0,0,0) = pow(-2*p, nn) * Boys(Boys_arg, nn);
        Rcur.push_back(Rinit);
    }

    //for (i=0; i<t_max; i++)
    for (i=0; i<tt; i++)
    {
        for (nn=0; nn<(t_max-i); nn++)
        {
            Rcur.at(nn)(i+1,0,0) = Rcp(0) * Rcur.at(nn+1)(i,0,0);

            if (i > 0)
            {
                Rcur.at(nn)(i+1,0,0) += i * Rcur.at(nn+1)(i-1,0,0);
            }
        }
    }

    //for (i=0; i<(t_max+1); i++)
    for (i=0; i<tt+1; i++)
    {
        //for (j=0; j<(t_max); j++)
        for (j=0; j<uu; j++)
        {
            for (nn=0; nn<(t_max-i-j); nn++)
            {
                Rcur.at(nn)(i, j+1, 0) = Rcp(1) * Rcur.at(nn+1)(i,j,0);

                if (j>0)
                {
                    Rcur.at(nn)(i, j+1, 0) += j * Rcur.at(nn+1)(i,j-1,0);
                }
            }
        }
    }

    //for (i=0; i < (t_max+1); i++)
    for (i=0; i<tt+1; i++)
    {
        //for (j=0; j<(t_max+1); j++)
        for (j=0; j<uu+1; j++)
        {
            //for (k=0; k<(t_max); k++)
            for (k=0; k<vv; k++)
            {
                for (nn=0; nn<(t_max-i-j-k); nn++)
                {
                    Rcur.at(nn)(i, j, k+1) = Rcp(2) * Rcur.at(nn+1)(i,j,k);

                    if (k>0)
                    {
                        Rcur.at(nn)(i, j, k+1) += k*Rcur.at(nn+1)(i,j,k-1);
                    }
                }
            }
        }
    }
    return Rcur.at(n)(t,u,v);
}

double Hartree_Integrals::Boys(double x, double n)
{
   double F;
   //int N = (2*n+1)*15;//Boys_N;
   //int N = Boys_N;
   //N = 50;
   //int N = 30;
   int N;

   double exp_term = exp(-x);

   if (x > 50)
   {
       N = 12;
       F_Boys = zeros(N+1);
       Set_Boys_Start(N);

       F = Boys_Start / pow(2.0, N+1) * sqrt(M_PI/pow(x, 2*N+1));
   }

   else
   {
       N = 12;
       F_Boys = zeros(N+1);
       F = 0;
       double sum=0;
       int M;
       for (int j=0; j<100; j++)
       {
           sum = pow(2*x, j);
           M = 2*N+1;
           while (M < (2*N+2+2*j))
           {
               sum /= M;
               M += 2;
           }
           F += sum;
       }
       F *= exp_term;
   }

   F_Boys(N) = F;
   for (int i = N; i > n; i--)
   {
       F_Boys(i-1) = (2*x*F_Boys(i) + exp_term)/(2*i-1);
   }

   return F_Boys(n);
}

double Hartree_Integrals::Nuclei_Nuclei_Interaction()
{
    // Finner nuclei-nuclei interaksjons termen
    // |psi> er bølgefunksjonen vår
    //
    // Bidraget er fra <psi|(Z_j Z_i)/(|R_j-R_i|)|psi>
    //
    // Verdt å merke seg at R og Z er konstant igjennom hele programmet
    // Så dette blir et konstant ledd
    // Pot = sum_i sum_(j<i) Z_j Z_i / |R_j - R_i|

    int i,j;
    double Pot = 0;
    for (i=0; i<n_Nuclei; i++)
    {
        for (j=0; j<i; j++)
        {
            Pot += Z(i) * Z(j) / sqrt(dot(R.row(i)-R.row(j), R.row(i)-R.row(j)));
        }
    }
    return Pot;
}

double Hartree_Integrals::Electron_Electron_Interaction_Single(int ind1, int ind2, int ind3, int ind4,
                                                               int a1, int a2, int a3, int a4,
                                                               int bas1, int bas2, int bas3, int bas4,
                                                               int E_counter1, int E_counter2)
{
    double x1,x2,x3,x4;
    x1 = alpha(ind1, bas1);
    x2 = alpha(ind2, bas2);
    x3 = alpha(ind3, bas3);
    x4 = alpha(ind4, bas4);

    int A1, B1, C1, D1, E1, F1;
    int A2, B2, C2, D2, E2, F2;
    rowvec RP, RQ;

    A1 = Potenser(ind1,0);
    B1 = Potenser(ind1,1);
    C1 = Potenser(ind1,2);
    D1 = Potenser(ind3,0);
    E1 = Potenser(ind3,1);
    F1 = Potenser(ind3,2);
    A2 = Potenser(ind2,0);
    B2 = Potenser(ind2,1);
    C2 = Potenser(ind2,2);
    D2 = Potenser(ind4,0);
    E2 = Potenser(ind4,1);
    F2 = Potenser(ind4,2);

    double value=0, P, Q, tempvalue;

    P = x1+x3;
    Q = x2+x4;
    RP = (x1*R.row(a1) + x3*R.row(a3))/P;
    RQ = (x2*R.row(a2) + x4*R.row(a4))/Q;

    int t1, u1, v1;
    int t2, u2, v2;

    Set_R_ijk(P*Q/(P+Q), (A1+D1+A2+D2+2), (B1+E1+B2+E2+2), (C1+F1+C2+F2+2), RP, RQ);

    for (t1=0; t1<(A1+D1+1); t1++)
    {
        for (u1=0; u1<(B1+E1+1); u1++)
        {
            for (v1=0; v1<(C1+F1+1); v1++)
            {
                tempvalue = 0;
                for (t2=0; t2<(A2+D2+1); t2++)
                {
                    for (u2=0; u2<(B2+E2+1); u2++)
                    {
                        for (v2=0; v2<(C2+F2+1); v2++)
                        {
                            tempvalue += pow(-1.0, (t2+u2+v2))
                                    * R_ijk.at(0)((t1+t2), (u1+u2), (v1+v2))
                                    * E_ij.at(E_counter2)(A2, D2, t2)
                                    * E_ij.at(E_counter2+1)(B2, E2, u2)
                                    * E_ij.at(E_counter2+2)(C2, F2, v2);
                        }
                    }
                }
                value += tempvalue * E_ij.at(E_counter1)(A1, D1, t1)* E_ij.at(E_counter1+1)(B1, E1, u1)* E_ij.at(E_counter1+2)(C1, F1, v1);
            }
        }
    }

    value *= 2*pow(acos(-1.0), 2.5) / (P*Q*sqrt(P+Q));
    return value;
}
