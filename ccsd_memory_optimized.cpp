#include "ccsd_memory_optimized.h"

CCSD_Memory_optimized::CCSD_Memory_optimized(int n_N, vec zz, mat rr, string B_s, int n_Elec, int ran, int siz, Hartree_Fock_Solver *Hartfock)
{
    n_Nuclei = n_N;
    Z = zz;
    R = rr;
    Basis_Set = B_s;
    n_Electrons = n_Elec;
    rank = ran;
    size = siz;
    HartFock = Hartfock;
}

double CCSD_Memory_optimized::CCSD(double toler, bool print_stuff)
{
    /*
     *
     *
     * *Optimization thoughts for programmer:
     *
     *
     *          A MAJOR CHANGE IN PARALLEL IMPLEMENTATION COMMING!
     *
     *          WE SHOULD DEFINE T2_NEW AND PART1_MPI AS TWO ONE DIMENSIONAL ARRAYS OF SIZE = MAX NR OF WORK OF ANY OF THE NODES
     *          THEN EACH OF THE NODES HAVE TWO ARRAYS, ONE OF WHICH IS THEIR WORK, THE OTHER CONTAINS THE WORK RECIEVED
     *          THEN THERE IS ONE COMMUNICATION AND THIS IS MAPPED
     *
     *          => THEN THE TEMP MPI ARRAY's CONTENT IS REPLACED BY NEXT NODES WORK AND THIS IS DISTRIBUTED
     *
     *          THIS REDUCE LEADING TERM STORAGE REQUIRED TO EXACTLY (1 + 1/16) * n_v^2 n_o^2
     *
     *          IF THIS IS NOT GOOD ENOUGH WE MUST STORE SOME MOs ON HARDDISK - CURRENTLY THEY ARE DISTRIBUTED
     *
     *
     *
     *          Current optimizations should then (after the rest) be aimed at the NON ITERATIVE part of CCSD,
     *          HF is fast enough for systems that CCSD can handle, iterations will be very fast, pre-work can and will be further optimized, see Prepear_MOs function
     *
     *
     *
     *
     *          W_4 and W_2 is called depending upon a and b, meaning we can utilize this and only calculate the a and b that we need for each node
     *          in the already parallel t2_new calculation. May cause some aditional calculations, but it is possible
     *          This would also reduce memory needs significantly for what is now our largest array by far
     *          Could also map W_4 out afterwards, since Part1_MPI stores this
     *          Alternatively combine t2_MPI with Part1_MPI and save memory
     *
     *          First solution saves more memeory but more calculation expencive
     *
     *          W4(0,1) er 0 i nedre halvdel, trenger ikke lagres
     *
     *          Division by 2 uneccasary, but OPENMP might be used if we include it :O combine OPENMP and MPI
     *          Can replace i/2 with I++ and I etc as done in meny functions
     *
     *
     *
     * */

    // General starting stuff
    bool continue_ccsd = true;
    double convergance_check;
    int Itermax = 1000;
    iter = 1;
    E_old = 0;
    E_new = 0;

    // Initializion
    double E_HF;

    //Hartree_Fock_Solver HartFock(n_Nuclei, Z, R, Basis_Set, n_Electrons, print_stuff, rank, size);
    E_HF = HartFock->get_Energy(toler); // Calc hartree fock energy
    Matrix_Size = 2*HartFock->ReturnMatrixSize(); // This will be twice as big as in Hatree Fock, because of spin up and down
    unocc_orb = Matrix_Size - n_Electrons; // Number of unocupied orbitals
    mat temp_matrix;
    Speed_Elec = (int) n_Electrons/2;
    Speed_Occ = (int) unocc_orb/2;
    Zero_Matrix = zeros(Speed_Occ, Speed_Occ); // Trying to decleare this here and we will only have to link it later

    // Transform to MO basis, only initialize the MOs we need
    c = HartFock->ReturnC();
    //Integrals = HartFock->Return_Indexed_Q();
    Prepear_AOs();
    fs = Fill_FS(HartFock->return_eigval_F());
    HartFock = NULL; // Remove this memory,  not needed anymore
    R.clear();
    Z.clear();
    //HartFock->Delete_Everything();

    // First we must do some initial mapping to determine how this is gonna work
    Map_T2_For_MPI(); // This initializes the T2 new which will be sent through MPI, and finds what parts of MO9 the different threads will need
    Map_Part1_For_MPI(); // This function initializes everything else that MPI will need, this function consists of pretty much everything in one

    // We now pull out the parts of the MO basis we need, we split up the this 4D array which is a
    // n^4 array, n^4 = (n_u + n_o)^4, where n_u = unoccupied orbitals in HF basis and n_o is
    // occupied orbitals in HF basis. To the right of each term is the nr of elements stored inside it

    ///////////////////////////////////
    // SIZE OF INTEGRAL (MO) STORAGE //
    ///////////////////////////////////

    // We store the following integrals (not updated):  // Solution to storage problems:
    // MO9          // 1/32 n_u^2*(n_o+n_u)^2           // Distributed
    // MO10         // 1/16 n_u*n_o*(n_o+n_u)^2         // Distributed
    // MO6          // 1/32 n_o^2*(n_o+n_u)^2           // Possible to distribute in parallel implementation
    // MOLeftovers  // 1/32 n_o^2 n_u^2                 // Hard to distribute in parallel implementation, but smallest part of MOs, but possible

    // Define variables used in the CCSD method:

    // Not so optimized 4D array
    tau3.set_size(n_Electrons, n_Electrons); // 1/4 n_o^2 n_u^2 // Use this for hyperspeed, can be removed if memoryproblems

    // Optimized 4D arrays
    t2.set_size(unocc_orb, n_Electrons);    // 1/4 n_o^2 n_u^2     // Must be stored!, check size again
    // t2_new is replaced by a double**     // 1/16 n_o^2 n_u^2
    W_1.set_size(n_Electrons, n_Electrons); // 1/4 n_o^4        // Relatively small
    W_2.set_size(n_Electrons, n_Electrons); // 1/4 n_o^3 n_u    // Relatively small
    W_3.set_size(n_Electrons, n_Electrons); // 3/8 n_o^3 n_u    // Relatively small
    W_4.set_size(unocc_orb, n_Electrons);   // 1/2 n_o^2 n_u^2    // Biggest variable

    // In total we end up storing the following leading terms of MOs (updated 28.04.13, somewhat different now, mostly smaller):
    // n_u^4         1/32           =  0.03
    // n_u^3 n_o     (3/32 + 1/2)   =  0.59
    // n_u^2 n_o^2   (9/32  + 7/2)  =  3.78
    // n_o^3 n_u     (5/8 + 3/32)   =  0.72
    // n_o^4         1/32           =  0.03

    // Optimized 2D arrays
    integ2_2D = zeros(unocc_orb, Speed_Elec); // This is used as mapping of MOs
    integ3_2D = zeros(unocc_orb, Speed_Occ); // This is used as mapping of MOs
    integ4_2D = zeros(unocc_orb, Speed_Occ); // This is used as mapping of MOs
    integ5_2D = zeros(unocc_orb, Speed_Occ); // This is used as mapping of MOs
    integ6_2D = zeros(unocc_orb, Speed_Elec); // This is used as mapping of MOs
    integ7_2D = zeros(unocc_orb, Speed_Elec); // This is used as mapping of MOs
    integ8_2D = zeros(n_Electrons, Speed_Elec); // This is used as mapping of MOs
    integ9_2D = zeros(unocc_orb, Speed_Elec); // This is used as mapping of MOs
    integ10_2D = zeros(unocc_orb, Speed_Elec); // This is used as mapping of MOs
    D3 = zeros(unocc_orb, Speed_Occ); // F3 in text
    D2 = zeros(n_Electrons, Speed_Elec); // F2 in text
    T_1 = zeros(unocc_orb, Speed_Elec); // Compact storage for T1 amplitudes
    T_1_new = zeros(unocc_orb, Speed_Elec); // Compact storage for T1 amplitudes
    D1 = zeros(unocc_orb, Speed_Elec); // F1 in text
    FS_AI = zeros(unocc_orb, Speed_Elec); // Compact storage even for Fock eigenvalues
    FS_AB = zeros(unocc_orb, Speed_Occ); // Compact storage even for Fock eigenvalues
    FS_IJ = zeros(n_Electrons, Speed_Elec); // Compact storage even for Fock eigenvalues
    DEN_AI = zeros(unocc_orb, Speed_Elec); // This will also contain DEN_ABIJ the 4D denominator but stored as 2D for memoryconcerns
    tau1 = zeros(n_Electrons, Speed_Elec); // 2D mapping of tau

    // Initialize our 4D arrays
    for (int i = 0; i < n_Electrons; i++)
    {
        for (int j = i+1; j < n_Electrons; j++)
        {
            temp_matrix = zeros(unocc_orb, Speed_Occ);
            tau3(i,j) = temp_matrix;

            temp_matrix = zeros(n_Electrons, Speed_Elec);
            W_1(i,j) = temp_matrix;

            temp_matrix = zeros(unocc_orb, Speed_Elec);
            W_2(i,j) = temp_matrix;
        }
    }

    for (int i = 0; i < n_Electrons; i++)
    {
        for (int j = 0; j < n_Electrons; j++)
        {
            if ((i+j)%2 == 1)
            {
                temp_matrix = zeros(Speed_Occ, Speed_Elec);
                W_3(i,j) = temp_matrix;
            }

            else
            {
                temp_matrix = zeros(unocc_orb, Speed_Elec);
                W_3(i,j) = temp_matrix;
            }
        }
    }

    for (int i = 0; i < unocc_orb; i++)
    {
        for (int j = 0; j < n_Electrons; j++)
        {
            temp_matrix = zeros(unocc_orb, Speed_Elec);
            t2(i,j) = temp_matrix;
            W_4(i,j) = temp_matrix;
        }
    }

    // Coupled Cluster is an iterative process, meaning we make an initial guess of t1 and t2 and calculate the energy
    // Then update t1 and t2 and check for convergance in the energy

    if (print_stuff == true)
    {
        cout << "No more memory needed from here!" << endl;
    }

    // Find denominator matrix, this should not change during our calculations

    for (int a = 0; a < unocc_orb; a++)
    {
        for (int i = 0; i < n_Electrons; i++)
        {
            DEN_AI(a/2,i/2) = (fs(i,i) - fs(a+n_Electrons,a+n_Electrons));
            i++;
        }
        a++;
    }

    for (int a = 0+1; a < unocc_orb; a++)
    {
        for (int i = 1; i < n_Electrons; i++)
        {
            DEN_AI(a/2+Speed_Occ,i/2) = (fs(i,i) - fs(a+n_Electrons,a+n_Electrons));
            i++;
        }
        a++;
    }

    // Set up (1 - delta_pq) * fs, these calculations are the same throughout our entire calculation and
    // does not need to be recalculated. Since it is 2 dimensional arrays we can store them
    // Storing in this way reduce the number of elements stored by a factor of n_Electrons * Unoccupied_Orbitals
    // And then later we reduce this by a factor of 1/2. These terms are as far as i know always 0 in restricted HF basis CCSD with closed shells :O
    // The reason is that eigenvalues of fock matrix (fs) is on diagonal matrix form and we have 1 - delta()
    // This also accomodates matrix matrix multiplication operations without any mapping

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

    // The initial guess of t1 and t2 is made here, t1 is guessed to be 0
    // Compact form of T2

    for (int b = 0; b < unocc_orb; b++)
    {
        for (int j = 0; j < n_Electrons; j++)
        {
            Fill_integ2_2D(b,j);
            for (int a = 0; a < unocc_orb; a++)
            {
                for (int i = 0; i< n_Electrons; i++)
                {
                    if ((a+i+b+j)%2 == 0)
                    {
                        t2(a,i)(b/2+b%2*Speed_Occ,j/2) = integ2_2D(a/2+a%2*Speed_Occ,i/2) / (DEN_AI(a/2+a%2*Speed_Occ,i/2) + DEN_AI(b/2+b%2*Speed_Occ,j/2));
                    }
                }
            }
        }
    }

    // Here starts the iteration, most calculations are taken into functions, even tho they are only called at one place
    // This is to get a more compact algo which is in my opinion more easy to read

    // Define something to measure time
    clock_t start;
    clock_t slutt;
    double time_measured;

    Fill_tau();
    E_new = Calc_Energy(); // Starting energy, with t1=0 guess

    if (print_stuff == true)
    {
        cout << "Energi: " << E_new << " Steg: " << iter << endl;
    }

    while (continue_ccsd == true)
    {
        start = clock();
        E_old = E_new;

        // Update intermediates and find new amplitudes
        // This is crazy implementation. Its insane..
        Fill_F1(); // Initialize everything
        Fill_W1_and_W3(); // Calculate everything
        Distribute_Part1(); // Distribute everything
        Fill_F2(); // These are the F1 * something else part
        Fill_F3(); // These are the F1 * something else part

        // These will go in parallel. Also wide functions to minimize communication.
        // Part two of parallel implementation will be here
        //Fill_W2(); // Moved to Fill_W1_and_W3() function
        Fill_W4_MPI();
        Fill_W4();

        //Fill_t1_new(); // Dont do this anymore
        // This is part 2 of parallel implementation
        Fill_t2_new(); // This function takes care of T1 and T2 in the parallel version of CCSD.
        // This is to ensure we minimize communication

        // Update amplitudes
        Map_T_new();
        T_1 = T_1_new / DEN_AI; // This division is not in parallel, but it is n^2
        Fill_tau(); // NOT IN PARALLEL, n^4

        // Calc energy
        E_new = Calc_Energy(); // NOT IN PARALLEL, n^4
        iter += 1;

        slutt = clock();

        // Output energy
        if (print_stuff == true)
        {
            time_measured = (double) (slutt-start)/CLOCKS_PER_SEC;
            cout << "Energi: " << E_new << " Steg: " << iter << " Time/iter: " << time_measured << endl;
        }

        // Check convergance
        convergance_check = sqrt((E_new - E_old) * (E_new - E_old));
        if (convergance_check < toler || iter >= Itermax)
        {
            continue_ccsd = false;
        }
    }

    // Iterations finished, cout result and return

    if (print_stuff == true)
    {
        cout << "CCSD correction [au] = " << E_new << endl;
    }

    return E_new+E_HF;
}

int CCSD_Memory_optimized::EqualFunc(int a, int b)
{
    // Wondering if a == b? Want 1 returned if this is true? And 0 if they are different?
    // No problem! Simply follow these easy to use and simple step by step instructions:

    // Step 1: Write EqualFunc
    // Step 2: Write thereafter (a, b), if these are the variables you wish to check
    // Step 3: Reconsider if you didnt really want to check (a/2, b/2), or (a%2, b%2)

    // This method is also eggstremely effective in problemsolving with you spouse, child or cellmate
    // Not recomended for use during flue season or pregnancy.

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
    // Do you have a bunch of integrals and amplitudes already calculated, but no way of finding the energy in CCSD?
    // No problem! Simply follow these step by step instructions

    // Step 1: Ensure you do have tau = t2(ij,ab) + t1(ai)t1(bj) - t1(aj)t1(bi) stored
    // Step 2: Call Calc_Energy()

    double E1 = 0;
    for (int i = 0; i < n_Electrons; i++)
    {
        for (int j = i+1; j < n_Electrons; j++)
        {
            Fill_integ4_2D(i,j);
            E1 += accu(tau3(i,j) % integ4_2D);
        }
    }
    E1 *= 0.5;
    E1 += accu(FS_AI % T_1);
    return E1;
}

void CCSD_Memory_optimized::Fill_W1_and_W3()
{
    // Matrix symmetric with W_1(i,j) = W_1(j,i), the terms on the "diagonal" of i==j will be equal to 0 (<-- !)
    // We actually dont even need to store anything in W_1(j,i) since it is accessed symmetricly later on also
    // This halves our storage requirements. Also compress matrix. This in total reduce storage by >75%

    for (int i = 0; i < n_Electrons; i++)
    {
        for (int j = i+1; j < n_Electrons; j++)
        {
            Fill_integ8_2D(i, j);
            W_1(i,j) = integ8_2D;
        }
    }

    // Optimized W_1, W_3 and parts of F1, F2, F3, T1

    int K, L, E, I, J;

    int index_counter;
    int INDEX_CHECK;

    index_counter = 0;
    K = 0;
    for (int k = 0; k < n_Electrons; k++)
    {
        L = 0;
        for (int l = 1; l < k; l++)
        {
            INDEX_CHECK = K*Speed_Elec + Speed_Elec-L;

            if (INDEX_CHECK % size == rank)
            {
                Fill_integ6_2D(l,k);
                Fill_integ4_2D(l,k);

                for (int i = 0; i < n_Electrons; i++)
                {
                    Part1_MPI[rank][index_counter] = accu(integ6_2D(span(Speed_Occ, Speed_Occ + Speed_Occ-1), i/2) % T_1(span(Speed_Occ, Speed_Occ+Speed_Occ-1), l/2));
                    i++;
                    index_counter += 1;
                }

                I = 0;
                for (int i = 0; i < n_Electrons; i++)
                {
                    E = 0;
                    for (int e = 1; e < unocc_orb; e++)
                    {
                        Part1_MPI[rank][index_counter] = -integ6_2D(E+Speed_Occ, I) + accu(integ4_2D(span(0, Speed_Occ-1), E) % T_1(span(0, Speed_Occ-1), I));
                        e++;
                        E++;
                        index_counter += 1;
                    }
                    i++;
                    I++;
                }

                I = 0;
                for (int i = 1; i < n_Electrons; i++)
                {
                    E = 0;
                    for (int e = 0; e < unocc_orb; e++)
                    {
                        Part1_MPI[rank][index_counter] = -integ6_2D(E, I) + accu(integ4_2D(span(Speed_Occ, Speed_Occ + Speed_Occ-1), E) % T_1(span(Speed_Occ, Speed_Occ+Speed_Occ-1), I));
                        e++;
                        index_counter += 1;
                        E++;
                    }
                    i++;
                    I++;
                }

                I = 0;
                for (int i = 0; i < n_Electrons; i++)
                {
                    J = I;
                    for (int j = i+1; j < n_Electrons; j++)
                    {
                        Part1_MPI[rank][index_counter] = accu(integ6_2D(span(Speed_Occ, Speed_Occ-1+Speed_Occ), I) % T_1(span(Speed_Occ, Speed_Occ-1+Speed_Occ), J))
                                - accu(integ6_2D(span(0, Speed_Occ-1),J) % T_1(span(0, Speed_Occ-1), I))
                                - 0.5*accu(integ4_2D % tau3(i,j));
                        j++;
                        J++;
                        index_counter += 1;
                    }
                    i++;
                    I++;
                }

                I = 0;
                for (int i = 1; i < n_Electrons; i++)
                {
                    J = I+1;
                    for (int j = i+1; j < n_Electrons; j++)
                    {

                        Part1_MPI[rank][index_counter] = accu(integ6_2D(span(0, Speed_Occ-1), I) % T_1(span(0, Speed_Occ-1), J))
                                - accu(integ6_2D(span(Speed_Occ, Speed_Occ - 1 + Speed_Occ),J) % T_1(span(Speed_Occ, Speed_Occ -1 + Speed_Occ), I))
                                - 0.5*accu(integ4_2D % tau3(i,j));
                        j++;
                        J++;
                        index_counter += 1;
                    }
                    i++;
                    I++;
                }
            }
            l++;
            L++;
        }

        for (int l = k+1; l < n_Electrons; l++)
        {
            INDEX_CHECK = K*Speed_Elec + Speed_Elec-L;

            if (INDEX_CHECK % size == rank)
            {

                Fill_integ6_2D(k,l);
                for (int i = 0; i < n_Electrons; i++)
                {
                    Part1_MPI[rank][index_counter] = accu(integ6_2D(span(Speed_Occ, Speed_Occ + Speed_Occ-1), i/2) % T_1(span(Speed_Occ, Speed_Occ+Speed_Occ-1), l/2));
                    i++;
                    index_counter += 1;
                }
            }
            l++;
            L++;
        }
        k++;
        K++;
    }

    K = 0;
    for (int k = 1; k < n_Electrons; k++)
    {
        L = 0;
        for (int l = 0; l < k; l++)
        {

            INDEX_CHECK = K*Speed_Elec + Speed_Elec-L;

            if (INDEX_CHECK % size == rank)
            {

                Fill_integ6_2D(l,k);
                Fill_integ4_2D(l,k);

                for (int i = 1; i < n_Electrons; i++)
                {
                    Part1_MPI[rank][index_counter] = accu(integ6_2D(span(0, Speed_Occ-1), i/2) % T_1(span(0, Speed_Occ-1), l/2));
                    i++;
                    index_counter += 1;
                }

                I = 0;
                for (int i = 0; i < n_Electrons; i++)
                {
                    E = 0;
                    for (int e = 1; e < unocc_orb; e++)
                    {
                        Part1_MPI[rank][index_counter] = -integ6_2D(E+Speed_Occ, I) + accu(integ4_2D(span(0, Speed_Occ-1), E) % T_1(span(0, Speed_Occ-1), I));
                        e++;
                        E++;
                        index_counter += 1;
                    }
                    i++;
                    I++;
                }

                I = 0;
                for (int i = 1; i < n_Electrons; i++)
                {
                    E = 0;
                    for (int e = 0; e < unocc_orb; e++)
                    {
                        Part1_MPI[rank][index_counter] = -integ6_2D(E, I) + accu(integ4_2D(span(Speed_Occ, Speed_Occ + Speed_Occ-1), E) % T_1(span(Speed_Occ, Speed_Occ+Speed_Occ-1), I));
                        e++;
                        E++;
                        index_counter += 1;
                    }
                    i++;
                    I++;
                }

                I = 0;
                for (int i = 0; i < n_Electrons; i++)
                {
                    J = I;
                    for (int j = i+1; j < n_Electrons; j++)
                    {
                        Part1_MPI[rank][index_counter] = accu(integ6_2D(span(Speed_Occ, Speed_Occ-1+Speed_Occ), I) % T_1(span(Speed_Occ, Speed_Occ-1+Speed_Occ), J))
                                - accu(integ6_2D(span(0, Speed_Occ-1),J) % T_1(span(0, Speed_Occ-1), I))
                                - 0.5*accu(integ4_2D % tau3(i,j));
                        index_counter += 1;
                        j++;
                        J++;
                    }
                    i++;
                    I++;
                }

                I = 0;
                for (int i = 1; i < n_Electrons; i++)
                {
                    J = I+1;
                    for (int j = i+1; j < n_Electrons; j++)
                    {
                        Part1_MPI[rank][index_counter] = accu(integ6_2D(span(0, Speed_Occ-1), I) % T_1(span(0, Speed_Occ-1), J))
                                - accu(integ6_2D(span(Speed_Occ, Speed_Occ - 1 + Speed_Occ),J) % T_1(span(Speed_Occ, Speed_Occ -1 + Speed_Occ), I))
                                - 0.5*accu(integ4_2D % tau3(i,j));
                        index_counter += 1;
                        j++;
                        J++;
                    }
                    i++;
                    I++;
                }
            }

            L++;
            l++;
        }

        for (int l = k+1; l < n_Electrons; l++)
        {
            INDEX_CHECK = K*Speed_Elec + Speed_Elec-L;

            if (INDEX_CHECK % size == rank)
            {

                Fill_integ6_2D(k,l);

                for (int i = 1; i < n_Electrons; i++)
                {
                    Part1_MPI[rank][index_counter] = accu(integ6_2D(span(0, Speed_Occ-1), i/2) % T_1(span(0, Speed_Occ-1), l/2));
                    i++;
                    index_counter += 1;
                }
            }
            l++;
            L++;
        }
        k++;
        K++;
    }

    K = 0;
    for (int k = 0; k < n_Electrons; k++)
    {
        L = 0;
        for (int l = 0; l < k+1; l++)
        {
            INDEX_CHECK = K*Speed_Elec + Speed_Elec-L;

            if (INDEX_CHECK % size == rank)
            {

                Fill_integ6_2D(l,k);
                Fill_integ4_2D(l,k);

                for (int i = 0; i < n_Electrons; i++)
                {
                    Part1_MPI[rank][index_counter] = accu(integ6_2D(span(0, Speed_Occ-1), i/2) % T_1(span(0, Speed_Occ-1), l/2));
                    i++;
                    index_counter += 1;
                }

                I = 0;
                for (int i = 0; i < n_Electrons; i++)
                {
                    E = 0;
                    for (int e = 0; e < unocc_orb; e++)
                    {
                        Part1_MPI[rank][index_counter] = -integ6_2D(E, I) + accu(integ4_2D(span(0, Speed_Occ-1), E) % T_1(span(0, Speed_Occ-1), I));
                        index_counter += 1;
                        e++;
                        E++;
                    }
                    i++;
                    I++;
                }

                I = 0;
                for (int i = 0; i < n_Electrons; i++)
                {
                    J = I+1;
                    for (int j = i+2; j < n_Electrons; j++)
                    {
                        Part1_MPI[rank][index_counter] = accu(integ6_2D(span(0, Speed_Occ-1),I) % T_1(span(0, Speed_Occ-1), J))
                                - accu(integ6_2D(span(0, Speed_Occ-1),J) % T_1(span(0, Speed_Occ-1), I))
                                - 0.5*accu(integ4_2D % tau3(i,j));
                        index_counter += 1;
                        j++;
                        J++;
                    }
                    i++;
                    I++;
                }
            }

            l++;
            L++;
        }

        for (int l = k+2; l < n_Electrons; l++)
        {
            INDEX_CHECK = K*Speed_Elec + Speed_Elec-L;

            if (INDEX_CHECK % size == rank)
            {
                Fill_integ6_2D(k,l);
                for (int i = 0; i < n_Electrons; i++)
                {
                    Part1_MPI[rank][index_counter] = accu(integ6_2D(span(0, Speed_Occ-1), i/2) % T_1(span(0, Speed_Occ-1), l/2));
                    i++;
                    index_counter += 1;
                }
            }
            l++;
            L++;
        }
        k++;
        K++;
    }

    K = 0;
    for (int k = 1; k < n_Electrons; k++)
    {
        L = 0;
        for (int l = 1; l < k+1; l++)
        {
            INDEX_CHECK = K*Speed_Elec + Speed_Elec-L;

            if (INDEX_CHECK % size == rank)
            {

                Fill_integ6_2D(l,k);
                Fill_integ4_2D(l,k);

                for (int i = 0; i < n_Electrons; i++)
                {
                    Part1_MPI[rank][index_counter] = accu(integ6_2D(span(Speed_Occ, Speed_Occ + Speed_Occ-1), i/2) % T_1(span(Speed_Occ, Speed_Occ+Speed_Occ-1), l/2));
                    i++;
                    index_counter += 1;
                }

                I = 0;
                for (int i = 1; i < n_Electrons; i++)
                {
                    E = 0;
                    for (int e = 1; e < unocc_orb; e++)
                    {
                        Part1_MPI[rank][index_counter] = -integ6_2D(E+Speed_Occ, I) + accu(integ4_2D(span(Speed_Occ, Speed_Occ + Speed_Occ-1), E) % T_1(span(Speed_Occ, Speed_Occ+Speed_Occ-1), I));
                        e++;
                        E++;
                        index_counter += 1;
                    }
                    i++;
                    I++;
                }

                I = 0;
                for (int i = 1; i < n_Electrons; i++)
                {
                    J = I+1;
                    for (int j = i+2; j < n_Electrons; j++)
                    {
                        Part1_MPI[rank][index_counter] = accu(integ6_2D(span(Speed_Occ, Speed_Occ-1+Speed_Occ), I) % T_1(span(Speed_Occ, Speed_Occ-1+Speed_Occ), J))
                                -  accu(integ6_2D(span(Speed_Occ, Speed_Occ - 1 + Speed_Occ),J) % T_1(span(Speed_Occ, Speed_Occ -1 + Speed_Occ), I))
                                 - 0.5*accu(integ4_2D % tau3(i,j));
                        j++;
                        J++;
                        index_counter += 1;
                    }
                    i++;
                    I++;
                }
            }

            l++;
            L++;
        }

        for (int l = k+2; l < n_Electrons; l++)
        {

            INDEX_CHECK = K*Speed_Elec + Speed_Elec-L;

            if (INDEX_CHECK % size == rank)
            {
                Fill_integ6_2D(k,l);
                for (int i = 0; i < n_Electrons; i++)
                {
                    Part1_MPI[rank][index_counter] = accu(integ6_2D(span(Speed_Occ, Speed_Occ + Speed_Occ-1), i/2) % T_1(span(Speed_Occ, Speed_Occ+Speed_Occ-1), l/2));
                    i++;
                    index_counter += 1;
                }
            }
            l++;
            L++;
        }
        k++;
        K++;
    }







    ////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////
    /// Now add W_2 and D_2 and D_3 and D1
    ////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////

    for (int i = 0; i < n_Electrons; i++)
    {
        for (int j = i+1; j < n_Electrons; j++)
        {
            Fill_integ6_2D(i,j);
            for (int k = 0; k < unocc_orb; k++)
            {
                for (int l = 0; l < Speed_Elec; l++)
                {
                    W_2(i,j)(k,l) = integ6_2D(k,l);
                }
            }
        }
    }

    int A, M;

    // Optimized W_2 version
    A = 0;
    for (int a = 0; a < unocc_orb; a++)
    {
        M = 0;
        for (int m = 0; m < n_Electrons; m++)
        {

            INDEX_CHECK = a/2*Speed_Elec + Speed_Elec - m/2;
            if (INDEX_CHECK % size == rank)
            {

                Fill_integ7_2D(a, m);
                Fill_integ5_2D(a, m);

                for (int e = 0; e < unocc_orb; e++)
                {
                    Part1_MPI[rank][index_counter] = accu(integ5_2D(e/2, span()) % T_1(span(0, Speed_Occ-1), m/2).t());
                    index_counter += 1;
                    e++;
                }

                I = 0;
                for (int i = 0; i < n_Electrons; i++)
                {
                    J = I+1;
                    for (int j = i+2; j < n_Electrons; j++)
                    {
                        Part1_MPI[rank][index_counter] = accu(integ7_2D(span(0, Speed_Occ-1), J) % T_1(span(0, Speed_Occ-1),I))
                                - accu(integ7_2D(span(0, Speed_Occ-1), I) % T_1(span(0, Speed_Occ-1),J))
                                + 0.5*accu(integ5_2D % tau3(i,j));
                        index_counter += 1;
                        j++;
                        J++;
                    }
                    i++;
                    I++;
                }
            }
            m++;
            M++;
        }
        a++;
        A++;
    }

    A = 0;
    for (int a = 1; a < unocc_orb; a++)
    {
        M = 0;
        for (int m = 1; m < n_Electrons; m++)
        {
            INDEX_CHECK = a/2*Speed_Elec + Speed_Elec - m/2;
            if (INDEX_CHECK % size == rank)
            {

                Fill_integ7_2D(a, m);
                Fill_integ5_2D(a, m);

                for (int e = 1; e < unocc_orb; e++)
                {
                    Part1_MPI[rank][index_counter] = accu(integ5_2D(e/2+Speed_Occ, span()) % T_1(span(Speed_Occ, unocc_orb-1), m/2).t());
                    e++;
                    index_counter += 1;
                }

                I = 0;
                for (int i = 0; i < n_Electrons; i++)
                {
                    J = I+1;
                    for (int j = i+2; j < n_Electrons; j++)
                    {
                        Part1_MPI[rank][index_counter] = accu(integ7_2D(span(0, Speed_Occ-1), J) % T_1(span(0, Speed_Occ-1),I))
                                - accu(integ7_2D(span(0, Speed_Occ-1), I) % T_1(span(0, Speed_Occ-1),J))
                                + 0.5*accu(integ5_2D % tau3(i,j));
                        index_counter += 1;
                        j++;
                        J++;
                    }
                    i++;
                    I++;
                }

                I = 0;
                for (int i = 1; i < n_Electrons; i++)
                {
                    J = I+1;
                    for (int j = i+2; j < n_Electrons; j++)
                    {
                        Part1_MPI[rank][index_counter] = accu(integ7_2D(span(Speed_Occ, Speed_Occ+Speed_Occ-1), J) % T_1(span(Speed_Occ, Speed_Occ+Speed_Occ-1),I))
                                - accu(integ7_2D(span(Speed_Occ, Speed_Occ+Speed_Occ-1), I) % T_1(span(Speed_Occ, Speed_Occ+Speed_Occ-1),J))
                                + 0.5*accu(integ5_2D % tau3(i,j));
                        index_counter += 1;
                        j++;
                        J++;
                    }
                    i++;
                    I++;
                }
            }

            M++;
            m++;
        }
        a++;
        A++;
    }

    A = 0;
    for (int a = 1; a < unocc_orb; a++)
    {
        M = 0;
        for (int m = 0; m < n_Electrons; m++)
        {
            INDEX_CHECK = a/2*Speed_Elec + Speed_Elec - m/2;
            if (INDEX_CHECK % size == rank)
            {

                Fill_integ7_2D(a, m);
                Fill_integ5_2D(a, m);


                for (int e = 1; e < unocc_orb; e++)
                {
                    Part1_MPI[rank][index_counter] = accu(integ5_2D(e/2+Speed_Occ, span()) % T_1(span(0, Speed_Occ-1), m/2).t());
                    e++;
                    index_counter += 1;
                }

                I = 0;
                for (int i = 0; i < n_Electrons; i++)
                {
                    J = I;
                    for (int j = i+1; j < n_Electrons; j++)
                    {
                        Part1_MPI[rank][index_counter] = accu(integ7_2D(span(0, Speed_Occ-1), J) % T_1(span(0, Speed_Occ-1),I))
                                - accu(integ7_2D(span(Speed_Occ, Speed_Occ+Speed_Occ-1), I) % T_1(span(Speed_Occ, Speed_Occ+Speed_Occ-1),J))
                                + 0.5*accu(integ5_2D % tau3(i,j));
                        index_counter += 1;
                        j++;
                        J++;
                    }
                    i++;
                    I++;
                }

                I = 0;
                for (int i = 1; i < n_Electrons; i++)
                {
                    J = I+1;
                    for (int j = i+1; j < n_Electrons; j++)
                    {
                        Part1_MPI[rank][index_counter] = accu(integ7_2D(span(Speed_Occ, Speed_Occ+Speed_Occ-1), J) % T_1(span(Speed_Occ, Speed_Occ+Speed_Occ-1),I))
                                - accu(integ7_2D(span(0, Speed_Occ-1), I) % T_1(span(0, Speed_Occ-1),J))
                                + 0.5*accu(integ5_2D % tau3(i,j));
                        index_counter += 1;
                        j++;
                        J++;
                    }
                    i++;
                    I++;
                }
            }
            m++;
            M++;
        }
        a++;
        A++;
    }

    A = 0;
    for (int a = 0; a < unocc_orb; a++)
    {
        M = 0;
        for (int m = 1; m < n_Electrons; m++)
        {

            INDEX_CHECK = a/2*Speed_Elec + Speed_Elec - m/2;
            if (INDEX_CHECK % size == rank)
            {

                Fill_integ5_2D(a,m);
                Fill_integ7_2D(a,m);

                for (int e = 0; e < unocc_orb; e++)
                {
                    Part1_MPI[rank][index_counter] = accu(integ5_2D(e/2, span()) % T_1(span(Speed_Occ, unocc_orb-1), m/2).t());
                    index_counter += 1;
                    e++;
                }

                I = 0;
                for (int i = 0; i < n_Electrons; i++)
                {
                    J = I;
                    for (int j = i+1; j < n_Electrons; j++)
                    {
                        Part1_MPI[rank][index_counter] = accu(integ7_2D(span(0, Speed_Occ-1), J) % T_1(span(0, Speed_Occ-1),I))
                                - accu(integ7_2D(span(Speed_Occ, Speed_Occ+Speed_Occ-1), I) % T_1(span(Speed_Occ, Speed_Occ+Speed_Occ-1),J))
                                + 0.5*accu(integ5_2D % tau3(i,j));
                        index_counter += 1;
                        j++;
                        J++;
                    }
                    i++;
                    I++;
                }

                I = 0;
                for (int i = 1; i < n_Electrons; i++)
                {
                    J = I+1;
                    for (int j = i+1; j < n_Electrons; j++)
                    {
                        Part1_MPI[rank][index_counter] = accu(integ7_2D(span(Speed_Occ, Speed_Occ+Speed_Occ-1), J) % T_1(span(Speed_Occ, Speed_Occ+Speed_Occ-1),I))
                                - accu(integ7_2D(span(0, Speed_Occ-1), I) % T_1(span(0, Speed_Occ-1),J))
                                + 0.5*accu(integ5_2D % tau3(i,j));
                        index_counter += 1;
                        j++;
                        J++;
                    }
                    i++;
                    I++;
                }
            }
            m++;
            M++;
        }
        a++;
        A++;
    }

    for (int e = 0; e < unocc_orb; e++)
    {
        for (int m = 0; m < n_Electrons; m++)
        {
            INDEX_CHECK = e/2*Speed_Elec + Speed_Elec - m/2;
            if (INDEX_CHECK % size == rank)
            {
                Fill_integ2_2D(e,m);
                for (int i = 0; i < n_Electrons; i++)
                {
                    Part1_MPI[rank][index_counter] = 0.5*accu(integ2_2D % t2.at(e,i));
                    index_counter += 1;
                    i++;
                }
            }
            m++;
        }
    }

    for (int e = 0; e < unocc_orb; e++)
    {
        for (int m = 1; m < n_Electrons; m++)
        {
            INDEX_CHECK = e/2*Speed_Elec + Speed_Elec - m/2;
            if (INDEX_CHECK % size == rank)
            {
                Fill_integ2_2D(e,m);
                for (int i = 1; i < n_Electrons; i++)
                {
                    Part1_MPI[rank][index_counter] = 0.5*accu(integ2_2D % t2.at(e,i));
                    index_counter += 1;
                    i++;
                }
            }
            m++;
        }
    }


    for (int e = 0; e < unocc_orb; e++)
    {
        for (int m = 0; m < n_Electrons; m++)
        {
            INDEX_CHECK = e/2 * Speed_Elec + Speed_Elec - m/2;
            if (INDEX_CHECK % size == rank)
            {
                Fill_integ2_2D(e,m);

                // F1 term
                Part1_MPI[rank][index_counter] = accu(integ2_2D % T_1);
                index_counter += 1;

                for (int a = 0; a < unocc_orb; a++)
                {
                    Part1_MPI[rank][index_counter] = 0.5*accu(integ2_2D % t2.at(a,m));
                    index_counter += 1;
                    a++;
                }
            }
            m++;
        }
        e++;
    }

    for (int e = 0; e < unocc_orb; e++)
    {
        for (int m = 1; m < n_Electrons; m++)
        {
            INDEX_CHECK = e/2 * Speed_Elec + Speed_Elec - m/2;
            if (INDEX_CHECK % size == rank)
            {
                Fill_integ2_2D(e,m);
                for (int a = 0; a < unocc_orb; a++)
                {
                    Part1_MPI[rank][index_counter] = 0.5*accu(integ2_2D(span(Speed_Occ, unocc_orb-1), span()) % t2.at(a,m)(span(Speed_Occ, unocc_orb-1), span()));
                    index_counter += 1;
                    a++;
                }
            }
            m++;
        }
        e++;
    }

    for (int e = 1; e < unocc_orb; e++)
    {
        for (int m = 0; m < n_Electrons; m++)
        {
            INDEX_CHECK = e/2 * Speed_Elec + Speed_Elec - m/2;
            if (INDEX_CHECK % size == rank)
            {
                Fill_integ2_2D(e,m);
                for (int a = 1; a < unocc_orb; a++)
                {
                    Part1_MPI[rank][index_counter] = 0.5*accu(integ2_2D(span(0,Speed_Occ-1), span()) % t2.at(a,m)(span(0,Speed_Occ-1), span()));
                    index_counter += 1;
                    a++;
                }
            }
            m++;
        }
        e++;
    }

    for (int e = 1; e < unocc_orb; e++)
    {
        for (int m = 1; m < n_Electrons; m++)
        {
            INDEX_CHECK = e/2 * Speed_Elec + Speed_Elec - m/2;
            if (INDEX_CHECK % size == rank)
            {
                Fill_integ2_2D(e,m);

                // F1 term
                Part1_MPI[rank][index_counter] = accu(integ2_2D % T_1);
                index_counter += 1;

                for (int a = 1; a < unocc_orb; a++)
                {
                    Part1_MPI[rank][index_counter] = 0.5*accu(integ2_2D % t2.at(a,m));
                    index_counter += 1;
                    a++;
                }
            }
            m++;
        }
        e++;
    }




    for (int a = 0; a < unocc_orb; a++)
    {
        for (int k = 0; k < n_Electrons; k++)
        {
            INDEX_CHECK = a/2 * Speed_Elec + Speed_Elec - k/2;
            if (INDEX_CHECK % size == rank)
            {
                Fill_integ5_2D(a, k);
                for (int i = 0; i < k; i++)
                {
                    Part1_MPI[rank][index_counter] = 0.5*accu(integ5_2D % tau3(i,k));
                    index_counter += 1;
                    i++;
                }

                // i will be less than k when k = odd number
                // We want i to be an even number since a is even number, hence we start at i = k+1
                for (int i = (k+1+(k+1)%2); i < n_Electrons; i++)
                {
                    Part1_MPI[rank][index_counter] = -0.5*accu(integ5_2D % tau3(k,i));
                    index_counter += 1;
                    i++;
                }
            }
        }
        a++;
    }

    for (int a = 1; a < unocc_orb; a++)
    {
        for (int k = 0; k < n_Electrons; k++)
        {
            INDEX_CHECK = a/2 * Speed_Elec + Speed_Elec - k/2;
            if (INDEX_CHECK % size == rank)
            {
                Fill_integ5_2D(a, k);
                for (int i = 1; i < k; i++)
                {
                    Part1_MPI[rank][index_counter] = 0.5*accu(integ5_2D % tau3(i,k));
                    index_counter += 1;
                    i++;
                }

                // i will be less than k when k = even number
                // We want i to be an odd number since a is an odd number, hence we start at i = k+1
                for (int i = k+1+(k)%2; i < n_Electrons; i++)
                {
                    Part1_MPI[rank][index_counter] = -0.5*accu(integ5_2D % tau3(k,i));
                    index_counter += 1;
                    i++;
                }
            }
        }
        a++;
    }


    // Communicate
    int work;
    for (int X = 0; X < size; X++)
    {
        work = WORK_EACH_NODE_Part1(X);
        MPI_Bcast(Part1_MPI[X], work, MPI_DOUBLE, X, MPI_COMM_WORLD);
    }
}

void CCSD_Memory_optimized::Fill_W2()
{
    // Matrix symmetric with W_2(i,j) = -W_2(j,i) and always 0 on the diagonal where i == j (<-- !)
    // Since symmetry is used later on also W_2(j,i) can be whatever, we do not need to access this ever
    // This halves our storage requirements, also compress matrix = reduce storage 75% + nothing on diagonal

    // This part is done in Fill_W1_And_W3 because faster with MPI this way
}

void CCSD_Memory_optimized::Fill_W4()
{
    //int M, E;

    // We only need to calculate the terms of W_4 of which t2 amplitudes is not equal to 0 later on
    // since W_4 only appear one place, in which it is multiplied by t2 on a accumulative basis
    // We can in fact skip a bunch of calculations using this

    int index_counter, INDEX_CHECK;

    for (int K = 0; K < size; K++)
    {
        index_counter = 0;
        for (int a = 0; a < unocc_orb; a++)
        {
            for (int m = 0; m < n_Electrons; m++)
            {
                INDEX_CHECK = a/2 * Speed_Elec + Speed_Elec - m/2;
                if (INDEX_CHECK % size == K)
                {
                    for (int e = 0; e < unocc_orb; e++)
                    {
                        for(int i = 0; i < n_Electrons; i++)
                        {
                            W_4(a,i)(e/2,m/2) = Part1_MPI[K][index_counter];
                            index_counter += 1;
                            i++;
                        }
                        e++;
                    }
                }
                m++;
            }
            a++;
        }

        for (int a = 1; a < unocc_orb; a++)
        {
            for (int m = 0; m < n_Electrons; m++)
            {
                INDEX_CHECK = a/2 * Speed_Elec + Speed_Elec - m/2;
                if (INDEX_CHECK % size == K)
                {
                    for (int e = 0; e < unocc_orb; e++)
                    {
                        for (int i = 1; i < n_Electrons; i++)
                        {
                            W_4(a,i)(e/2,m/2) = Part1_MPI[K][index_counter];
                            index_counter += 1;
                            i++;
                        }
                        e++;
                    }

                    for (int e = 1; e < unocc_orb; e++)
                    {
                        for (int i = 0; i < n_Electrons; i++)
                        {
                            W_4(a,i)(e/2+Speed_Occ,m/2) = Part1_MPI[K][index_counter];
                            index_counter += 1;
                            i++;
                        }
                        e++;
                    }
                }
                m++;
            }
            a++;
        }

        for (int a = 0; a < unocc_orb; a++)
        {
            for (int m = 1; m < n_Electrons; m++)
            {
                INDEX_CHECK = a/2 * Speed_Elec + Speed_Elec - m/2;
                if (INDEX_CHECK % size == K)
                {
                    for (int e = 0; e < unocc_orb; e++)
                    {
                        for (int i = 1; i < n_Electrons; i++)
                        {
                            W_4(a,i)(e/2,m/2) = Part1_MPI[K][index_counter];
                            i++;
                            index_counter += 1;
                        }
                        e++;
                    }

                    for (int e = 1; e < unocc_orb; e++)
                    {
                        for (int i = 0; i < n_Electrons; i++)
                        {
                            W_4(a,i)(e/2+Speed_Occ, m/2) = Part1_MPI[K][index_counter];
                            i++;
                            index_counter += 1;
                        }
                        e++;
                    }
                }
                m++;
            }
            a++;
        }

        for (int a = 1; a < unocc_orb; a++)
        {
            for (int m = 1; m < n_Electrons; m++)
            {
                INDEX_CHECK = a/2 * Speed_Elec + Speed_Elec - m/2;
                if (INDEX_CHECK % size == K)
                {
                    for (int e = 0; e < unocc_orb; e++)
                    {
                        for (int i = 0; i < n_Electrons; i++)
                        {
                            W_4(a,i)(e/2,m/2) = Part1_MPI[K][index_counter];
                            i++;
                            index_counter += 1;
                        }
                        e++;
                    }

                    for (int e = 1; e < unocc_orb; e++)
                    {
                        for (int i = 1; i < n_Electrons; i++)
                        {
                            W_4(a,i)(e/2+Speed_Occ,m/2) = Part1_MPI[K][index_counter];
                            i++;
                            index_counter += 1;
                        }
                        e++;
                    }
                }
                m++;
            }
            a++;
        }
    }
}

void CCSD_Memory_optimized::Fill_F1()
{
    // Initialize F1, F2 and F3
    D1 = FS_AI;
    D2 = FS_IJ;
    D3 = FS_AB;

    // Also initialize T_1_new here, since it will be easier
    T_1_new = FS_AI;

    // Run F1 in parallel in the W_1_and_W_3 function :O
}

void CCSD_Memory_optimized::Fill_F2()
{
    int I,M;

    // Do most of F2 in parallel, except for n^3 terms

    M = 0;
    for (int m = 0; m < n_Electrons; m++)
    {
        I = 0;
        for (int i = 0; i < n_Electrons; i++)
        {
            D2.at(M,I) += accu(D1(span(0, Speed_Occ-1), m/2) % T_1(span(0, Speed_Occ-1), i/2));
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
            i++;
            I++;
        }
        m++;
        M++;
    }
}

void CCSD_Memory_optimized::Fill_F3()
{
    // F3 , or D3 as it is called in program, is VERY SPECIAL VARIABLE
    // this is extremely important feature described here
    // D3 is the only variable that is NOT stored as "compact storage"
    // It is simply compacted together without any flipping of indexes, just removed all zeroes simply
    // The reason we do not need to do anything special in terms of indexes here is that
    // this variable is ONLY CALLED AS .row, IN OTHER WORDS ONLY ROWS OF F3
    // This means we do not need to flip indexes, and you will see in T2 function that
    // indeed F3 is called as a function of (a) and not (a/2+a%2*Speed_Occ) as is the compact storage way
    // This is extremely important because you may think this is a bug, but its not
    // This is actually faster in this implementation because we can simply call index a or b, and not have the division of 2


    // Most of F3 is done in parallel, except for n^3 terms

    int E;
    for (int a = 0; a < unocc_orb; a++)
    {
        E = 0;
        for (int e = 0; e < unocc_orb; e++)
        {
            D3.at(a,E) -= accu(D1.row(e/2) % T_1.row(a/2));
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
            e++;
            E++;
        }
        a++;
    }
}

void CCSD_Memory_optimized::Fill_tau()
{
    int A, B, I, J;
    I = 0;
    for (int i = 0; i < n_Electrons; i++)
    {
        J = I+1;
        for (int j = i+2; j < n_Electrons; j++)
        {
            A = 0;
            for (int a = 0; a < unocc_orb; a++)
            {
                B = A+1;
                for (int b = a+2; b < unocc_orb; b++)
                {
                    tau3(i,j)(A,B) = t2.at(a,i)(B,J)
                            + T_1(A, I) * T_1(B, J)
                            - T_1(A, J) * T_1(B, I);
                    tau3(i,j)(B,A) = -tau3(i,j)(A,B);
                    b++;
                    B++;
                }
                a++;
                A++;
            }

            A = 0;
            for (int a = 1; a < unocc_orb; a++)
            {
                B = A+1;
                for (int b = a+2; b < unocc_orb; b++)
                {
                    tau3(i,j)(A+Speed_Occ,B) = t2.at(a,i)(B+Speed_Occ, J);
                    tau3(i,j)(B+Speed_Occ,A) = -tau3(i,j)(A+Speed_Occ, B);
                    b++;
                    B++;
                }
                a++;
                A++;
            }
            j++;
            J++;
        }
        i++;
        I++;
    }

    I = 0;
    for (int i = 1; i < n_Electrons; i++)
    {
        J = I+1;
        for (int j = i+2; j < n_Electrons; j++)
        {
            A = 0;
            for (int a = 0; a < unocc_orb; a++)
            {
                B = A+1;
                for (int b = a+2; b < unocc_orb; b++)
                {
                    tau3(i,j)(A,B) = t2.at(a,i)(B,J);
                    tau3(i,j)(B,A) = -tau3(i,j)(A,B);
                    b++;
                    B++;
                }
                a++;
                A++;
            }

            A = 0;
            for (int a = 1; a < unocc_orb; a++)
            {
                B = A+1;
                for (int b = a+2; b < unocc_orb; b++)
                {
                    tau3(i,j)(A+Speed_Occ,B) = t2.at(a,i)(B+Speed_Occ,J)
                            + T_1(A+Speed_Occ, I) * T_1(B+Speed_Occ, J)
                            - T_1(A+Speed_Occ, J) * T_1(B+Speed_Occ, I);
                    tau3(i,j)(B+Speed_Occ,A) = -tau3(i,j)(A+Speed_Occ,B);
                    b++;
                    B++;
                }
                a++;
                A++;
            }
            J++;
            j++;
        }
        I++;
        i++;
    }

    I = 0;
    for (int i = 0; i < n_Electrons; i++)
    {
        J = I;
        for (int j = i+1; j < n_Electrons; j++)
        {
            A = 0;
            for (int a = 0; a < unocc_orb; a++)
            {
                B = A;
                for (int b = a+1; b < unocc_orb; b++)
                {
                    tau3(i,j)(A, B) = t2.at(a,i)(B+Speed_Occ,J)
                            + T_1(A, I) * T_1(B+Speed_Occ, J);
                    tau3(i,j)(B+Speed_Occ,A) = -tau3(i,j)(A,B);
                    b++;
                    B++;
                }
                a++;
                A++;
            }

            A = 0;
            for (int a = 1; a < unocc_orb; a++)
            {
                B = A+1;
                for (int b = a+1; b < unocc_orb; b++)
                {
                    tau3(i,j)(A+Speed_Occ,B) = t2.at(a,i)(B,J)
                            - T_1(A+Speed_Occ, J) * T_1(B, I);
                    tau3(i,j)(B,A) = -tau3(i,j)(A+Speed_Occ,B);
                    b++;
                    B++;
                }
                a++;
                A++;
            }
            J++;
            j++;
        }
        i++;
        I++;
    }

    I = 0;
    for (int i = 1; i < n_Electrons; i++)
    {
        J = I+1;
        for (int j = i+1; j< n_Electrons; j++)
        {
            A = 0;
            for (int a = 0; a < unocc_orb; a++)
            {
                B = A;
                for (int b = a+1; b < unocc_orb; b++)
                {
                    tau3(i,j)(A,B) = t2.at(a,i)(B+Speed_Occ,J)
                            - T_1(A, J) * T_1(B+Speed_Occ, I);
                    tau3(i,j)(B+Speed_Occ,A) = -tau3(i,j)(A,B);
                    b++;
                    B++;
                }
                a++;
                A++;
            }

            A = 0;
            for (int a = 1; a < unocc_orb; a++)
            {
                B = A+1;
                for (int b = a+1; b < unocc_orb; b++)
                {
                    tau3(i,j)(A+Speed_Occ,B) = t2.at(a,i)(B,J)
                            + T_1(A+Speed_Occ, I) * T_1(B, J);
                    tau3(i,j)(B,A) = -tau3(i,j)(A+Speed_Occ,B);
                    b++;
                    B++;
                }
                a++;
                A++;
            }
            J++;
            j++;
        }
        i++;
        I++;
    }
}

void CCSD_Memory_optimized::Fill_2D_tau(int a, int b)
{
    // 2D mapping of tau so we dont need double storage, implemented soon and will be used in the Fill_t2_new() function

    int A = a/2;
    int B = b/2;
    int I,J;

    if (a%2 == 0 && b%2 == 0)
    {
        I = 0;
        for (int i = 0; i < n_Electrons; i++)
        {
            J = I+1;
            for (int j = i+2; j < n_Electrons; j++)
            {
                tau1(I, J) = tau3(i,j)(A,B);
                tau1(J, I) = -tau1(I,J);
                j++;
                J++;
            }
            i++;
            I++;
        }

        I = 0;
        for (int i = 1; i < n_Electrons; i++)
        {
            J = I+1;
            for (int j = i+2; j < n_Electrons; j++)
            {
                tau1(I+Speed_Elec, J) = 0; //tau3(i,j)(a/2, b/2);
                tau1(J+Speed_Elec, I) = 0; //-tau1(i/2+Speed_Elec,j/2);
                j++;
                J++;
            }
            i++;
            I++;
        }
    }

    else if (a%2 == 0 && b%2 == 1)
    {
        I = 0;
        for (int i = 0; i < n_Electrons; i++)
        {
            J = I;
            for (int j = i+1; j < n_Electrons; j++)
            {
                tau1(I, J) = tau3(i,j)(A,B);
                tau1(J+Speed_Elec,i/2) = -tau1(I,J);
                j++;
                J++;
            }
            i++;
            I++;
        }

        I = 0;
        for (int i = 1; i < n_Electrons; i++)
        {
            J = I+1;
            for (int j = i+1; j < n_Electrons; j++)
            {
                tau1(I+Speed_Elec, J) = tau3(i,j)(A, B);
                tau1(J,I) = -tau1(I+Speed_Elec, J);
                j++;
                J++;
            }
            i++;
            I++;
        }
    }

    else if (a%2 == 1 && b%2 == 0)
    {
        I = 0;
        for (int i = 0; i < n_Electrons; i++)
        {
            J = I;
            for (int j = i+1; j < n_Electrons; j++)
            {
                tau1(I, J) = tau3(i,j)(A+Speed_Occ,B);
                tau1(J+Speed_Elec,I) = -tau1(I,J);
                j++;
                J++;
            }
            i++;
            I++;
        }

        I = 0;
        for (int i = 1; i < n_Electrons; i++)
        {
            J = I+1;
            for (int j = i+1; j < n_Electrons; j++)
            {
                tau1(I+Speed_Elec, J) = tau3(i,j)(A+Speed_Occ, B);
                tau1(J,I) = -tau1(I+Speed_Elec,J);
                j++;
                J++;
            }
            i++;
            I++;
        }
    }

    else
    {
        I = 0;
        for (int i = 0; i < n_Electrons; i++)
        {
            J = I+1;
            for (int j = i+2; j < n_Electrons; j++)
            {
                tau1(I, J) = tau3(i,j)(A+Speed_Occ,B);
                tau1(J, I) = -tau1(I,J);
                j++;
                J++;
            }
            i++;
            I++;
        }

        I = 0;
        for (int i = 1; i < n_Electrons; i++)
        {
            J = I+1;
            for (int j = i+2; j < n_Electrons; j++)
            {
                tau1(I+Speed_Elec, J) = tau3(i,j)(A+Speed_Occ, B);
                tau1(J+Speed_Elec, I) = -tau1(I+Speed_Elec,J);
                j++;
                J++;
            }
            i++;
            I++;
        }
    }
}

void CCSD_Memory_optimized::Fill_t1_new()
{
    // This is moved to Fill_t2_new() in parallel implementation :-O
    // Parts go to Fill_W1_and_W3 due to memory distribution
}

void CCSD_Memory_optimized::Fill_t2_new()
{
    // Benchmark H2O STO-3G basis set: -0.0501273 au

    // Optimized Version of T2 calculations
    // Runs in parallel
    // THIS IS THE ONLY PLACE MO9 IS ACCESSED - MEMORY SPREAD ON NODES IN a-b GRID
    // T2-new Stored in an array of size number of nodes times number of calcs for node
    // This is for communication optimization
    // Dont store more than one symetric term for communication minimization
    // Use all symmetries and skip meny calculations where one term will be 0 but the other not or both 0 etc
    // Try not to skip calculations on terms that are both not 0 :-D

    int index_counter;
    int A,B;
    double temp;
    int INDEX_CHECK;
    int work;

        index_counter = 0;

        for (int a = 0; a < unocc_orb; a++) // a even
        {
            A = a/2;
            for (int b = a+2; b < unocc_orb; b++) // b even
            {
                B = Speed_Occ - b/2;
                INDEX_CHECK = A*Speed_Occ+B;

                if (INDEX_CHECK % size == rank)
                {
                    Fill_integ3_2D(a, b);
                    Fill_integ9_2D(a, b);
                    Fill_2D_tau(a, b);

                    for (int i = 0; i < n_Electrons; i++) // i even
                    {
                        for (int j = i+2; j < n_Electrons; j++) // j even
                        {
                            // I_ab^ij, hopefully correct
                            T2_MPI[rank][index_counter] = -MOLeftovers(a/2, b/2)(j/2, i/2) + MOLeftovers(a/2, b/2)(i/2, j/2) // integ2.at(b,j)(a/2,i/2)

                            // P(ab) P(ij) [W_4] t_jk^bc, ARMADILLO
                                    + accu(t2.at(b,j) % W_4.at(a,i))
                                    - accu(t2.at(b,i) % W_4.at(a,j))
                                    - accu(t2.at(a,j) % W_4.at(b,i))
                                    + accu(t2.at(a,i) % W_4.at(b,j))

                            // - P(ab) [] t_k^b, ARMADILLO
                                    - accu(W_2(i,j)(a/2, span()) % T_1.row(b/2))
                                    + accu(W_2(i,j)(b/2, span()) % T_1.row(a/2))

                            // 1/2 [W_1] tau_kl^ab, ARMADILLO
                                    + 0.5*accu(W_1(i,j)(span(0, Speed_Elec-1), span()) % tau1(span(0, Speed_Elec-1), span())) // Half matrix = 0, skip this

                            // P(ij) t_jk^ab [F_2] + P(ab) t_ij^bc [F_3]_c^a, ARMADILLO
                                    - accu(t2.at(b,i)(span(0, Speed_Occ-1), j/2) % D3.row(a).t()) // Notice a is called, and not a/2+a%2*Speed_Occ
                                    + accu(t2.at(a,i)(span(0, Speed_Occ-1), j/2) % D3.row(b).t()) // This is because it is a special variable, descibed in Fill_F3 function

                                    + accu(t2.at(a,j)(b/2, span()) % D2.row(i/2))
                                    - accu(t2.at(a,i)(b/2, span()) % D2.row(j/2))

                            // P(ij) I_ab^cj t_i^c, ARMADILLO
                                    - accu(integ9_2D(span(0, Speed_Occ-1), i/2) % T_1(span(0, Speed_Occ-1), j/2))
                                    + accu(integ9_2D(span(0, Speed_Occ-1), j/2) % T_1(span(0, Speed_Occ-1), i/2))

                            // 0.5 I_ab^cd tau
                                    + 0.5*accu(integ3_2D(span(0, Speed_Occ-1), span()) % tau3(i,j)(span(0, Speed_Occ-1), span())); // Half matrix = 0, skip this


                            index_counter += 1;
                            j++;
                        }
                        i++;
                    }

                /* // Forbidden!
                for (int i = 1; i < n_Electrons; i++) // i odd
                {
                    for (int j = i+2; j < n_Electrons; j++) // j odd
                    {
                        j++;
                    }
                    i++;
                }
                */

                }
                b++;
            }
            a++;
        }

        for (int a = 1; a < unocc_orb; a++) // a odd
        {
            A = a/2;
            for (int b = a+2; b < unocc_orb; b++) // b odd
            {
                B = Speed_Occ - b/2;
                INDEX_CHECK = A*Speed_Occ+B;

                if (INDEX_CHECK % size == rank)
                {
                    Fill_integ3_2D(a, b);
                    Fill_integ9_2D(a, b);
                    Fill_2D_tau(a, b);
/*
                    for (int i = 0; i < n_Electrons; i++) // i even
                    {
                        for (int j = i+2; j < n_Electrons; j++) // j even
                        {

                            // This state illigal according to UN resolution 203

                            // I_ab^ij
                            T2_MPI[rank][index_counter] = 0 ;// integ2.at(b,j)(a/2+Speed_Occ,i/2)

                            // P(ab) P(ij) [W_4] t_jk^bc, ARMADILLO
                                    + accu(t2.at(b,j) % W_4.at(a,i))
                                    - accu(t2.at(b,i) % W_4.at(a,j))
                                    - accu(t2.at(a,j) % W_4.at(b,i))
                                    + accu(t2.at(a,i) % W_4.at(b,j))

                            // - P(ab) [] t_k^b, ARMADILLO
                                    - accu(W_2(i,j)(a/2+Speed_Occ, span()) % T_1.row(b/2+Speed_Occ))
                                    + accu(W_2(i,j)(b/2+Speed_Occ, span()) % T_1.row(a/2+Speed_Occ))

                            // 1/2 [W_1] tau_kl^ab, ARMADILLO
                                    + 0.5*accu(W_1(i,j)(span(0, Speed_Elec-1), span()) % tau1(span(0, Speed_Elec-1), span()))

                            // P(ij) t_jk^ab [F_2] + P(ab) t_ij^bc [F_3]_c^a, ARMADILLO
                                    - accu(t2.at(b,i)(span(0, Speed_Occ-1), j/2) % D3.row(a).t()) // Notice a is called, and not a/2+a%2*Speed_Occ
                                    + accu(t2.at(a,i)(span(0, Speed_Occ-1), j/2) % D3.row(b).t()) // This is because it is a special variable, descibed in Fill_F3 function

                                    + accu(t2.at(a,j)(b/2+Speed_Occ,span()) % D2.row(i/2))
                                    - accu(t2.at(a,i)(b/2+Speed_Occ, span()) % D2.row(j/2))

                            // P(ij) I_ab^cj t_i^c, ARMADILLO
                                    - accu(integ9_2D(span(0, Speed_Occ-1), i/2) % T_1(span(0, Speed_Occ-1), j/2))
                                   + accu(integ9_2D(span(0, Speed_Occ-1), j/2) % T_1(span(0, Speed_Occ-1), i/2))

                            // 0.5 I_ab^cd tau
                                    + 0.5*accu(integ3_2D(span(Speed_Occ, Speed_Occ+Speed_Occ-1), span()) % tau3(i,j)(span(Speed_Occ, Speed_Occ+Speed_Occ-1), span())); // Half matrix = 0, skip this

                            index_counter += 1;
                            j++;
                        }
                        i++;
                    }
*/
                    for (int i = 1; i < n_Electrons; i++) // i odd
                    {
                        for (int j = i+2; j < n_Electrons; j++) // j odd
                        {
                            // I_ab^ij
                            T2_MPI[rank][index_counter] = MOLeftovers(a/2, b/2)(i/2, j/2) - MOLeftovers(a/2, b/2)(j/2,i/2) // integ2.at(b,j)(a/2+Speed_Occ,i/2)

                            // P(ab) P(ij) [W_4] t_jk^bc, ARMADILLO
                                    + accu(t2.at(b,j) % W_4.at(a,i))
                                    - accu(t2.at(b,i) % W_4.at(a,j))
                                    - accu(t2.at(a,j) % W_4.at(b,i))
                                    + accu(t2.at(a,i) % W_4.at(b,j))

                            // - P(ab) [] t_k^b, ARMADILLO
                                    - accu(W_2(i,j)(a/2+Speed_Occ, span()) % T_1.row(b/2+Speed_Occ))
                                    + accu(W_2(i,j)(b/2+Speed_Occ, span()) % T_1.row(a/2+Speed_Occ))

                            // 1/2 [W_1] tau_kl^ab, ARMADILLO
                                    + 0.5*accu(W_1(i,j)(span(Speed_Elec, Speed_Elec+Speed_Elec-1), span()) % tau1(span(Speed_Elec, Speed_Elec+Speed_Elec-1), span()))

                            // P(ij) t_jk^ab [F_2] + P(ab) t_ij^bc [F_3]_c^a, ARMADILLO
                                    - accu(t2.at(b,i)(span(Speed_Occ, Speed_Occ+Speed_Occ-1), j/2) % D3.row(a).t()) // Notice a is called, and not a/2+a%2*Speed_Occ
                                    + accu(t2.at(a,i)(span(Speed_Occ, Speed_Occ+Speed_Occ-1), j/2) % D3.row(b).t()) // This is because it is a special variable, descibed in Fill_F3 function

                                    + accu(t2.at(a,j)(b/2+Speed_Occ, span()) % D2.row(i/2+Speed_Elec))
                                    - accu(t2.at(a,i)(b/2+Speed_Occ, span()) % D2.row(j/2+Speed_Elec))

                            // P(ij) I_ab^cj t_i^c, ARMADILLO
                                    - accu(integ9_2D(span(Speed_Occ, Speed_Occ+Speed_Occ-1), i/2) % T_1(span(Speed_Occ, Speed_Occ+Speed_Occ-1), j/2))
                                    + accu(integ9_2D(span(Speed_Occ, Speed_Occ+Speed_Occ-1), j/2) % T_1(span(Speed_Occ, Speed_Occ+Speed_Occ-1), i/2))

                            // 0.5 I_ab^cd tau
                                    + 0.5*accu(integ3_2D(span(Speed_Occ, Speed_Occ+Speed_Occ-1), span()) % tau3(i,j)(span(Speed_Occ, Speed_Occ+Speed_Occ-1), span()));

                            index_counter += 1;
                            j++;
                        }
                        i++;
                    }

                }
                b++;
            }
            a++;
        }

        for (int a = 0; a < unocc_orb; a++) // a even
        {
            A = a/2;
            for (int b = a+1; b < unocc_orb; b++) // b odd
            {
                B = Speed_Occ - b/2;
                INDEX_CHECK = A*Speed_Occ+B;

                if (INDEX_CHECK % size == rank)
                {
                    Fill_integ3_2D(a, b);
                    Fill_integ9_2D(a, b);
                    Fill_2D_tau(a, b);

                    for (int i = 0; i < n_Electrons; i++) // i even
                    {
                        for (int j = i+1; j < n_Electrons; j++) // j odd
                        {
                            // I_ab^ij
                            T2_MPI[rank][index_counter] = MOLeftovers(a/2, b/2)(i/2, j/2) // integ2.at(b,j)(a/2,i/2)

                            // P(ab) P(ij) [W_4] t_jk^bc, ARMADILLO
                                    - accu(t2.at(a,j)(span(Speed_Occ, Speed_Occ+Speed_Occ-1), span()) % W_4.at(b,i)(span(Speed_Occ, Speed_Occ+Speed_Occ-1), span()))
                                    - accu(t2.at(b,i)(span(0, Speed_Occ-1), span()) % W_4.at(a,j)(span(0, Speed_Occ-1), span()))
                                    + accu(t2.at(a,i) % W_4.at(b,j))
                                    + accu(t2.at(b,j) % W_4.at(a,i))

                            // P(ij) t_jk^ab [F_2] + P(ab) t_ij^bc [F_3]_c^a, ARMADILLO
                                    - accu(t2.at(b,i)(span(0, Speed_Occ-1), j/2) % D3.row(a).t()) // Notice a is called, and not a/2+a%2*Speed_Occ
                                    + accu(t2.at(a,i)(span(Speed_Occ, unocc_orb-1), j/2) % D3.row(b).t()) // This is because it is a special variable, descibed in Fill_F3 function

                                    + accu(t2.at(a,j)(b/2+Speed_Occ, span()) % D2.row(i/2))
                                    - accu(t2.at(a,i)(b/2+Speed_Occ, span()) % D2.row(j/2+Speed_Elec))

                            // 1/2 [W_1] tau_kl^ab, ARMADILLO
                                    + 0.5*accu(W_1(i,j) % tau1)

                            // - P(ab) [W_2] t_k^b, ARMADILLO
                                    - accu(W_2(i,j)(a/2, span()) % T_1.row(b/2+Speed_Occ))
                                    + accu(W_2(i,j)(b/2+Speed_Occ, span()) % T_1.row(a/2))

                            // P(ij) I_ab^cj t_i^c, ARMADILLO
                                    - accu(integ9_2D(span(Speed_Occ, Speed_Occ+Speed_Occ-1), i/2) % T_1(span(Speed_Occ, Speed_Occ+Speed_Occ-1), j/2))
                                    + accu(integ9_2D(span(0, Speed_Occ-1), j/2) % T_1(span(0, Speed_Occ-1), i/2))

                            // 0.5 I_ab^cd tau
                                    + 0.5*accu(integ3_2D % tau3(i,j));

                            index_counter += 1;
                            j++;
                        }
                        i++;
                    }

                    for (int i = 1; i < n_Electrons; i++) // i odd
                    {
                        for (int j = i+1; j < n_Electrons; j++) // j even
                        {
                            // I_ab^ij
                            T2_MPI[rank][index_counter] = -MOLeftovers(a/2, b/2)(j/2, i/2) // integ2.at(b,j)(a/2,i/2)

                            // P(ab) P(ij) [W_4] t_jk^bc, ARMADILLO
                                    - accu(t2.at(a,j) % W_4.at(b,i))
                                    - accu(t2.at(b,i) % W_4.at(a,j))
                                    + accu(t2.at(a,i)(span(Speed_Occ, Speed_Occ+Speed_Occ-1), span()) % W_4.at(b,j)(span(Speed_Occ, Speed_Occ+Speed_Occ-1), span()))
                                    + accu(t2.at(b,j)(span(0, Speed_Occ-1), span()) % W_4.at(a,i)(span(0, Speed_Occ-1), span()))

                            // P(ij) t_jk^ab [F_2] + P(ab) t_ij^bc [F_3]_c^a, ARMADILLO
                                    - accu(t2.at(b,i)(span(0, Speed_Occ-1), j/2) % D3.row(a).t()) // Notice a is called, and not a/2+a%2*Speed_Occ
                                    + accu(t2.at(a,i)(span(Speed_Occ, unocc_orb-1), j/2) % D3.row(b).t()) // This is because it is a special variable, descibed in Fill_F3 function

                                    + accu(t2.at(a,j)(b/2+Speed_Occ, span()) % D2.row(i/2+Speed_Elec))
                                    - accu(t2.at(a,i)(b/2+Speed_Occ, span()) % D2.row(j/2))

                            // 1/2 [W_1] tau_kl^ab, ARMADILLO
                                    + 0.5*accu(W_1(i,j) % tau1)

                            // - P(ab) [W_2] t_k^b, ARMADILLO
                                    - accu(W_2(i,j)(a/2, span()) % T_1.row(b/2+Speed_Occ))
                                    + accu(W_2(i,j)(b/2+Speed_Occ, span()) % T_1.row(a/2))

                            // P(ij) I_ab^cj t_i^c, ARMADILLO
                                    - accu(integ9_2D(span(0, Speed_Occ-1), i/2) % T_1(span(0, Speed_Occ-1), j/2))
                                    + accu(integ9_2D(span(Speed_Occ, Speed_Occ+Speed_Occ-1), j/2) % T_1(span(Speed_Occ, Speed_Occ+Speed_Occ-1), i/2))

                            // 0.5 I_ab^cd tau
                                    + 0.5*accu(integ3_2D % tau3(i,j));

                            index_counter += 1;
                            j++;
                        }
                        i++;
                    }
                }
                b++;
            }
            a++;
        }

        for (int a = 1; a < unocc_orb; a++) // a odd
        {
            A = a/2;
            for (int b = a+1; b < unocc_orb; b++) // b even
            {
                B = Speed_Occ - b/2;
                INDEX_CHECK = A*Speed_Occ+B;

                if (INDEX_CHECK % size == rank)
                {
                    Fill_integ3_2D(a, b);
                    Fill_integ9_2D(a, b);
                    Fill_2D_tau(a, b);

                    for (int i = 0; i < n_Electrons; i++) // i even
                    {
                        for (int j = i+1; j < n_Electrons; j++) // j odd
                        {
                            // I_ab^ij
                            T2_MPI[rank][index_counter] = -MOLeftovers(a/2, b/2)(j/2, i/2) // integ2.at(b,j)(a/2+Speed_Occ,i/2)

                            // P(ab) P(ij) [W_4] t_jk^bc, ARMADILLO
                                    - accu(t2.at(a,j) % W_4.at(b,i))
                                    - accu(t2.at(b,i) % W_4.at(a,j))
                                    + accu(t2.at(a,i)(span(0, Speed_Occ-1), span()) % W_4.at(b,j)(span(0, Speed_Occ-1), span()))
                                    + accu(t2.at(b,j)(span(Speed_Occ, Speed_Occ+Speed_Occ-1), span()) % W_4.at(a,i)(span(Speed_Occ, Speed_Occ+Speed_Occ-1), span()))

                            // P(ij) t_jk^ab [F_2] + P(ab) t_ij^bc [F_3]_c^a, ARMADILLO
                                    + accu(t2.at(a,j)(b/2,span()) % D2.row(i/2))
                                    - accu(t2.at(a,i)(b/2,span()) % D2.row(j/2+Speed_Elec))

                                    - accu(t2.at(b,i)(span(Speed_Occ, unocc_orb-1), j/2) % D3.row(a).t()) // Notice a is called, and not a/2+a%2*Speed_Occ
                                    + accu(t2.at(a,i)(span(0, Speed_Occ-1), j/2) % D3.row(b).t()) // This is because it is a special variable, descibed in Fill_F3 function

                            // 1/2 [W_1] tau_kl^ab, ARMADILLO
                                    + 0.5*accu(W_1(i,j) % tau1)

                            // - P(ab) [W_2] t_k^b, ARMADILLO
                                    - accu(W_2(i,j)(a/2+Speed_Occ, span()) % T_1.row(b/2))
                                    + accu(W_2(i,j)(b/2, span()) % T_1.row(a/2+Speed_Occ))

                            // P(ij) I_ab^cj t_i^c, ARMADILLO
                                    - accu(integ9_2D(span(Speed_Occ, Speed_Occ+Speed_Occ-1), i/2) % T_1(span(Speed_Occ, Speed_Occ+Speed_Occ-1), j/2))
                                    + accu(integ9_2D(span(0, Speed_Occ-1), j/2) % T_1(span(0, Speed_Occ-1), i/2))

                            // 0.5 I_ab^cd tau
                                    + 0.5*accu(integ3_2D % tau3(i,j));

                            index_counter += 1;
                            j++;
                        }
                        i++;
                    }

                    for (int i = 1; i < n_Electrons; i++) // i odd
                    {
                        for (int j = i+1; j < n_Electrons; j++) // j even
                        {
                            // I_ab^ij
                            T2_MPI[rank][index_counter] = MOLeftovers(a/2, b/2)(i/2, j/2) // integ2.at(b,j)(a/2+Speed_Occ,i/2)

                            // P(ab) P(ij) [W_4] t_jk^bc, ARMADILLO
                                    - accu(t2.at(a,j)(span(0, Speed_Occ-1), span()) % W_4.at(b,i)(span(0, Speed_Occ-1), span()))
                                    - accu(t2.at(b,i)(span(Speed_Occ, Speed_Occ+Speed_Occ-1), span()) % W_4.at(a,j)(span(Speed_Occ, Speed_Occ+Speed_Occ-1), span()))
                                    + accu(t2.at(a,i) % W_4.at(b,j))
                                    + accu(t2.at(b,j) % W_4.at(a,i))

                            // P(ij) t_jk^ab [F_2] + P(ab) t_ij^bc [F_3]_c^a, ARMADILLO
                                    + accu(t2.at(a,j)(b/2, span()) % D2.row(i/2+Speed_Elec))
                                    - accu(t2.at(a,i)(b/2, span()) % D2.row(j/2))

                                    - accu(t2.at(b,i)(span(Speed_Occ, unocc_orb-1), j/2) % D3.row(a).t()) // Notice a is called, and not a/2+a%2*Speed_Occ
                                    + accu(t2.at(a,i)(span(0, Speed_Occ-1), j/2) % D3.row(b).t()) // This is because it is a special variable, descibed in Fill_F3 function

                            // 1/2 [W_1] tau_kl^ab, ARMADILLO
                                    + 0.5*accu(W_1(i,j) % tau1)

                            // - P(ab) [W_2] t_k^b, ARMADILLO
                                    - accu(W_2(i,j)(a/2+Speed_Occ, span()) % T_1.row(b/2))
                                    + accu(W_2(i,j)(b/2, span()) % T_1.row(a/2+Speed_Occ))

                            // P(ij) I_ab^cj t_i^c, ARMADILLO
                                    - accu(integ9_2D(span(0, Speed_Occ-1), i/2) % T_1(span(0, Speed_Occ-1), j/2))
                                    + accu(integ9_2D(span(Speed_Occ, Speed_Occ+Speed_Occ-1), j/2) % T_1(span(Speed_Occ, Speed_Occ+Speed_Occ-1), i/2))

                            // 0.5 I_ab^cd tau
                                    + 0.5*accu(integ3_2D % tau3(i,j));

                            index_counter += 1;
                            j++;
                        }
                        i++;
                    }
                }
                b++;
            }
            a++;
        }

        for (int a = 0; a < unocc_orb; a++)
        {
            for (int i = 0; i < n_Electrons; i++)
            {
                INDEX_CHECK = a/2 * Speed_Elec + Speed_Elec - i/2;
                if (INDEX_CHECK % size == rank)
                {
                    Fill_integ10_2D(a, i);

                    // - 1/2 I_kl^ci t_kl^ca - 1/2 I_kl^cd t_kl^ca t_i^d
                    temp = 0;
                    for (int k = 0; k < n_Electrons; k++)
                    {
                        if (k%2 == 1)
                        {
                            temp += accu(W_3(i,k) % t2(a,k)(span(Speed_Occ, Speed_Occ+Speed_Occ-1), span()));
                        }

                        else
                        {
                            temp += accu(W_3(i,k) % t2(a,k));

                        }
                    }

                    T2_MPI[rank][index_counter] = 0.5*temp

                            // t_ik^ac [F_1]_c^k
                            + accu(D1 % t2(a,i))

                            // - t_k^a [F_2]_i^k
                            - accu(D2(span(i%2*Speed_Elec, i%2*Speed_Elec+Speed_Elec-1), i/2) % T_1.row(a/2+a%2*Speed_Occ).t())

                            // f_ac t_i^c
                            + accu(FS_AB(a/2, span()) % T_1(span(0, Speed_Occ-1), i/2).t())

                            // I_ka^ci t_k^c
                            - accu(integ10_2D % T_1);
                    index_counter += 1;
                }

                i++;
            }
            a++;
        }

        for (int a = 1; a < unocc_orb; a++)
        {
            for (int i = 1; i < n_Electrons; i++)
            {
                INDEX_CHECK = a/2 * Speed_Elec + Speed_Elec - i/2;
                if (INDEX_CHECK % size == rank)
                {
                    Fill_integ10_2D(a, i);

                    // - 1/2 I_kl^ci t_kl^ca - 1/2 I_kl^cd t_kl^ca t_i^d
                    temp = 0;
                    for (int k = 0; k < n_Electrons; k++)
                    {
                        if (k%2 == 1)
                        {
                            temp += accu(W_3(i,k) % t2(a,k));
                        }

                        else
                        {
                            temp += accu(W_3(i,k) % t2(a,k)(span(0, Speed_Occ-1), span()));
                        }
                    }

                    T2_MPI[rank][index_counter] = 0.5*temp

                            // t_ik^ac [F_1]_c^k
                            + accu(D1 % t2(a,i))

                            // - t_k^a [F_2]_i^k
                            - accu(D2(span(i%2*Speed_Elec, i%2*Speed_Elec+Speed_Elec-1), i/2) % T_1.row(a/2+a%2*Speed_Occ).t())

                            // f_ac t_i^c
                            + accu(FS_AB(a/2+Speed_Occ, span()) % T_1(span(Speed_Occ, unocc_orb-1), i/2).t())

                            - accu(integ10_2D % T_1);
                    index_counter += 1;
                }
                i++;
            }
            a++;
        }

    for (int K = 0; K < size; K++)
    {
        work = WORK_EACH_NODE(K);
        MPI_Bcast(T2_MPI[K], work, MPI_DOUBLE, K, MPI_COMM_WORLD);
    }
}

void CCSD_Memory_optimized::Map_T_new()
{
    // Mapping for parallel implementation!
    // Needs to resemble the CALCULATE T2 FUNCTION EXACTLY!!!!
    //                       BECAUSE WE NEED TO USE INDEX COUNTER
    //                      this is to make the transfer extreme easy

    double temp;
    int A,B, INDEX_CHECK;
    int index_counter;

    for (int K = 0; K < size; K++)
    {   // This loop is here in transformation to MPI only

        index_counter = 0;

        for (int a = 0; a < unocc_orb; a++) // a even
        {
            A = a/2;
            for (int b = a+2; b < unocc_orb; b++) // b even
            {
                B = Speed_Occ - b/2;
                INDEX_CHECK = A*Speed_Occ+B;

                if (INDEX_CHECK % size == K)
                {
                    for (int i = 0; i < n_Electrons; i++) // i even
                    {
                        for (int j = i+2; j < n_Electrons; j++) // j even
                        {
                            temp = T2_MPI[K][index_counter] / (DEN_AI(a/2+a%2*Speed_Occ,i/2)+DEN_AI(b/2+b%2*Speed_Occ, j/2));

                            t2(a,i)(b/2+b%2*Speed_Occ, j/2) = temp;
                            t2(b,i)(a/2+a%2*Speed_Occ, j/2) = -temp;
                            t2(a,j)(b/2+b%2*Speed_Occ, i/2) = -temp;
                            t2(b,j)(a/2+a%2*Speed_Occ, i/2) = temp;

                            j++;
                            index_counter += 1;
                        }
                        i++;
                    }
                }
                b++;
            }
            a++;
        }

        for (int a = 1; a < unocc_orb; a++) // a odd
        {
            A = a/2;
            for (int b = a+2; b < unocc_orb; b++) // b odd
            {
                B = Speed_Occ - b/2;
                INDEX_CHECK = A*Speed_Occ+B;

                if (INDEX_CHECK % size == K)
                {
                    /*
                    for (int i = 0; i < n_Electrons; i++) // i even
                    {
                        for (int j = i+2; j < n_Electrons; j++) // j even
                        {
                            temp = T2_MPI[K][index_counter] / (DEN_AI(a/2+a%2*Speed_Occ,i/2)+DEN_AI(b/2+b%2*Speed_Occ, j/2));

                            t2(a,i)(b/2+b%2*Speed_Occ, j/2) = temp;
                            t2(b,i)(a/2+a%2*Speed_Occ, j/2) = -temp;
                            t2(a,j)(b/2+b%2*Speed_Occ, i/2) = -temp;
                            t2(b,j)(a/2+a%2*Speed_Occ, i/2) = temp;

                            j++;
                            index_counter += 1;
                        }
                        i++;
                    }
*/
                    for (int i = 1; i < n_Electrons; i++) // i odd
                    {
                        for (int j = i+2; j < n_Electrons; j++) // j odd
                        {
                            temp = T2_MPI[K][index_counter] / (DEN_AI(a/2+a%2*Speed_Occ,i/2)+DEN_AI(b/2+b%2*Speed_Occ, j/2));

                            t2(a,i)(b/2+b%2*Speed_Occ, j/2) = temp;
                            t2(b,i)(a/2+a%2*Speed_Occ, j/2) = -temp;
                            t2(a,j)(b/2+b%2*Speed_Occ, i/2) = -temp;
                            t2(b,j)(a/2+a%2*Speed_Occ, i/2) = temp;

                            j++;
                            index_counter += 1;
                        }
                        i++;
                    }

                }
                b++;
            }
            a++;
        }

        for (int a = 0; a < unocc_orb; a++) // a even
        {
            A = a/2;
            for (int b = a+1; b < unocc_orb; b++) // b odd
            {
                B = Speed_Occ - b/2;
                INDEX_CHECK = A*Speed_Occ+B;

                if (INDEX_CHECK % size == K)
                {
                    for (int i = 0; i < n_Electrons; i++) // i even
                    {
                        for (int j = i+1; j < n_Electrons; j++) // j odd
                        {
                            temp = T2_MPI[K][index_counter] / (DEN_AI(a/2+a%2*Speed_Occ,i/2)+DEN_AI(b/2+b%2*Speed_Occ, j/2));

                            t2(a,i)(b/2+b%2*Speed_Occ, j/2) = temp;
                            t2(b,i)(a/2+a%2*Speed_Occ, j/2) = -temp;
                            t2(a,j)(b/2+b%2*Speed_Occ, i/2) = -temp;
                            t2(b,j)(a/2+a%2*Speed_Occ, i/2) = temp;

                            j++;
                            index_counter += 1;
                        }
                        i++;
                    }

                    for (int i = 1; i < n_Electrons; i++) // i odd
                    {
                        for (int j = i+1; j < n_Electrons; j++) // j even
                        {
                            temp = T2_MPI[K][index_counter] / (DEN_AI(a/2+a%2*Speed_Occ,i/2)+DEN_AI(b/2+b%2*Speed_Occ, j/2));

                            t2(a,i)(b/2+b%2*Speed_Occ, j/2) = temp;
                            t2(b,i)(a/2+a%2*Speed_Occ, j/2) = -temp;
                            t2(a,j)(b/2+b%2*Speed_Occ, i/2) = -temp;
                            t2(b,j)(a/2+a%2*Speed_Occ, i/2) = temp;

                            j++;
                            index_counter += 1;
                        }
                        i++;
                    }
                }
                b++;
            }
            a++;
        }

        for (int a = 1; a < unocc_orb; a++) // a odd
        {
            A = a/2;
            for (int b = a+1; b < unocc_orb; b++) // b even
            {
                B = Speed_Occ - b/2;
                INDEX_CHECK = A*Speed_Occ+B;

                if (INDEX_CHECK % size == K)
                {
                    for (int i = 0; i < n_Electrons; i++) // i even
                    {
                        for (int j = i+1; j < n_Electrons; j++) // j odd
                        {
                            temp = T2_MPI[K][index_counter] / (DEN_AI(a/2+a%2*Speed_Occ,i/2)+DEN_AI(b/2+b%2*Speed_Occ, j/2));

                            t2(a,i)(b/2+b%2*Speed_Occ, j/2) = temp;
                            t2(b,i)(a/2+a%2*Speed_Occ, j/2) = -temp;
                            t2(a,j)(b/2+b%2*Speed_Occ, i/2) = -temp;
                            t2(b,j)(a/2+a%2*Speed_Occ, i/2) = temp;

                            j++;
                            index_counter += 1;
                        }
                        i++;
                    }

                    for (int i = 1; i < n_Electrons; i++) // i odd
                    {
                        for (int j = i+1; j < n_Electrons; j++) // j even
                        {
                            temp = T2_MPI[K][index_counter] / (DEN_AI(a/2+a%2*Speed_Occ,i/2)+DEN_AI(b/2+b%2*Speed_Occ, j/2));

                            t2(a,i)(b/2+b%2*Speed_Occ, j/2) = temp;
                            t2(b,i)(a/2+a%2*Speed_Occ, j/2) = -temp;
                            t2(a,j)(b/2+b%2*Speed_Occ, i/2) = -temp;
                            t2(b,j)(a/2+a%2*Speed_Occ, i/2) = temp;

                            j++;
                            index_counter += 1;
                        }
                        i++;
                    }
                }
                b++;
            }
            a++;
        }

        for (int a = 0; a < unocc_orb; a++)
        {
            for (int i = 0; i < n_Electrons; i++)
            {
                INDEX_CHECK = a/2 * Speed_Elec + Speed_Elec - i/2;

                if (INDEX_CHECK % size == K)
                {
                    T_1_new(a/2, i/2) += T2_MPI[K][index_counter];
                    index_counter += 1;
                }
                i++;
            }
            a++;
        }

        for (int a = 1; a < unocc_orb; a++)
        {
            for (int i = 1; i < n_Electrons; i++)
            {
                INDEX_CHECK = a/2 * Speed_Elec + Speed_Elec - i/2;
                if (INDEX_CHECK % size == K)
                {
                    T_1_new(a/2+Speed_Occ, i/2) += T2_MPI[K][index_counter];
                    index_counter += 1;
                }
                i++;
            }
            a++;
        }
    }
}

void CCSD_Memory_optimized::Prepear_AOs()
{
    int matsize = Matrix_Size/2;

    int A,B, INDEX_CHECK;
    int SENDER;
    mat ttt;

    MO8.set_size(Speed_Elec, Speed_Elec);
    ttt = zeros(Speed_Elec, Speed_Elec);
    for (int a = 0; a < Speed_Elec; a++)
    {
        for (int b = a; b < Speed_Elec; b++)
        {

            MO8(a,b) = ttt;
        }
    }

    MO6_1.set_size(Speed_Elec, Speed_Elec);
    MO6_2.set_size(Speed_Elec, Speed_Elec);
    ttt = zeros(matsize, matsize);
    for (int a = 0; a < Speed_Elec; a++)
    {
        for (int b = a; b < Speed_Elec; b++)
        {
            MO6_1(a,b) = zeros(Speed_Occ, Speed_Elec);
            MO6_2(a,b) = zeros(Speed_Elec, Speed_Occ);
        }
    }

    MO4.set_size(Speed_Elec, Speed_Elec);
    ttt = zeros(Speed_Occ, Speed_Occ);
    for (int a = 0; a < Speed_Elec; a++)
    {
        for (int b = a; b < Speed_Elec; b++)
        {
            MO4(a,b) = ttt;
        }
    }

    MO9.set_size(Speed_Occ, Speed_Occ);
    ttt = zeros(matsize, matsize);
    for (int a = 0; a < Speed_Occ; a++)
    {
        for (int b = a; b < Speed_Occ; b++)
        {
            B = Speed_Occ - b;
            INDEX_CHECK = a*Speed_Occ+B;
            if (INDEX_CHECK % size == rank)
            {
                MO9(a,b) = ttt;
            }
        }
    }

    MO2.set_size(Speed_Occ, Speed_Elec);
    MO10.set_size(Speed_Occ, Speed_Elec);
    MO3.set_size(Speed_Occ, Speed_Elec);
    ttt = zeros(matsize, matsize);
    for (int a = 0; a < Speed_Occ; a++)
    {
        for (int b = 0; b < Speed_Elec; b++)
        {
            INDEX_CHECK = a * Speed_Elec + Speed_Elec - b;
            if (INDEX_CHECK % size == rank)
            {
                MO10(a,b) = ttt;
                MO2(a,b) = zeros(Speed_Occ, Speed_Elec);
                MO3(a,b) = zeros(Speed_Elec, Speed_Occ);
            }
        }
    }

    ttt = zeros(Speed_Elec, Speed_Elec);
    MOLeftovers.set_size(Speed_Occ, Speed_Occ);
    for (int a = 0; a < Speed_Occ; a++)
    {
        for (int b = 0; b < Speed_Occ; b++)
        {
            MOLeftovers(a,b) = ttt;
        }
    }

    // Now we are ready to construct MOs, however there is a problem.
    // Hartree Fock runs in parallel, we do not store all AOs on all threads
    // This is a major problem, and must be solved with communication
    // This is done here by checking first who has the information we need, and bribing this node, to get the information out of him
    // Because of this, our CCSD is now locked to our HF implementation, meaning we cannot use any other Hartree Fock implementation
    // unless we change the lines of code in this function

    int index_counter;
    int send_size = matsize * matsize;
    double *send_matrix_MPI = (double*) malloc(matsize*matsize*sizeof(double));
    mat recieve_matrix = zeros(matsize, matsize);
    cube send_cube = zeros(matsize, matsize, matsize);

    for (int a = 0; a < Speed_Elec; a++)
    {
        compact_mo = zeros(matsize, matsize, matsize);
        for (int k = 0; k < matsize; k++)
        {
            for (int i = 0; i < matsize; i++)
            {
                SENDER = (i+k)%size;
                if (SENDER == rank)
                {
                    recieve_matrix = HartFock->Return_Field_Q(i, k);
                    index_counter = 0;
                    for (int j = 0; j < matsize; j++)
                    {
                        for (int l = 0; l < matsize; l++)
                        {
                            send_matrix_MPI[index_counter] = recieve_matrix(j,l);
                            index_counter += 1;
                        }
                    }
                }

                // THIS COMMUNICATION MUST BE MOVED TO AFTER THE TRANSFORMATION AO -> MO!
                // THIS WILL BE VERY BIG SPEEDUP
                MPI_Bcast(send_matrix_MPI, send_size, MPI_DOUBLE, SENDER, MPI_COMM_WORLD);

                index_counter = 0;
                for (int j = 0; j < matsize; j++)
                {
                    for (int l = 0; l < matsize; l++)
                    {
                        recieve_matrix(j,l) = send_matrix_MPI[index_counter];
                        index_counter += 1;
                    }
                }

                for (int j = 0; j < matsize; j++)
                {
                    for (int l = 0; l < matsize; l++)
                    {
                        // We need to get these integrals here somehow
                        compact_mo(j,k,l) += c(i,a) * recieve_matrix(j,l); // This solution has the AOs distributed among threads and uses communication
                                //HartFock->Calc_Integrals_On_The_Fly(i,k,j,l); // This solution calculates integrals on the fly
                        //Integrals(Return_Integral_Index(i,j,k,l)); // This solution stored integrals as compounded index
                    }
                }
            }
        }

        compact_mo2 = zeros(matsize, matsize, matsize);
        for (int b = 0; b < matsize; b++)
        {
            for (int i = 0; i < matsize; i++)
            {
                for (int j = 0; j < matsize; j++)
                {
                    for (int k = 0; k < matsize; k++)
                    {
                        compact_mo2(b,j,k) += c(i,b) * compact_mo(i,j,k);
                    }
                }
            }
        }

        compact_mo = zeros(matsize, matsize, matsize);
        for (int b = 0; b < matsize; b++)
        {
            for (int g = 0; g < matsize; g++)
            {
                for (int i = 0; i < matsize; i++)
                {
                    for (int j = 0; j < matsize; j++)
                    {
                        compact_mo(b,g,j) += c(i,g) * compact_mo2(b,i,j);
                    }
                }
            }
        }

        compact_mo2 = zeros(matsize, matsize, matsize);

        for (int b = 0; b < matsize; b++)
        {
            for (int g = 0; g < matsize; g++)
            {
                for (int h = 0; h < matsize; h++)
                {
                    for (int i = 0; i < matsize; i++)
                    {
                        compact_mo2(b,g,h) += c(i,h) * compact_mo(b,g,i);
                    }
                }
            }
        }


        // HERE THE COMMUNICATION SHOULD BE PLACED!

            for (int i = a; i < Speed_Elec; i++)
            {
                    for (int b = Speed_Elec; b < matsize; b++)
                    {
                        for (int j = 0; j < Speed_Elec; j++)
                        {
                            MO6_1(a,i)(b-Speed_Elec,j) = compact_mo2(b,i,j);
                        }
                    }
            }

            for (int i = a; i < Speed_Elec; i++)
            {
                    for (int b = 0; b < Speed_Elec; b++)
                    {
                        for (int j = Speed_Elec; j < matsize; j++)
                        {
                            MO6_2(a,i)(b,j-Speed_Elec) = compact_mo2(b,i,j);
                        }
                    }
            }

            for (int i = a; i < Speed_Elec; i++)
            {
                    for (int b = 0; b < Speed_Elec; b++)
                    {
                        for (int j = 0; j < Speed_Elec; j++)
                        {
                            MO8(a,i)(b,j) = compact_mo2(b,i,j);
                        }
                    }
            }

            // MO4 will not be distributed, but it is relatively small
            for (int i = a; i < Speed_Elec; i++)
            {
                    for (int b = Speed_Elec; b < matsize; b++)
                    {
                        for (int j = Speed_Elec; j < matsize; j++)
                        {
                            MO4(a,i)(b-Speed_Elec,j-Speed_Elec) = compact_mo2(b,i,j);
                        }
                    }
            }

    }

    for (int a = Speed_Elec; a < matsize; a++)
    {
        compact_mo = zeros(matsize, matsize, matsize);
        for (int k = 0; k < matsize; k++)
        {
            for (int i = 0; i < matsize; i++)
            {

                SENDER = (i+k)%size;
                if (SENDER == rank)
                {
                    recieve_matrix = HartFock->Return_Field_Q(i, k);
                    index_counter = 0;
                    for (int j = 0; j < matsize; j++)
                    {
                        for (int l = 0; l < matsize; l++)
                        {
                            send_matrix_MPI[index_counter] = recieve_matrix(j,l);
                            index_counter += 1;
                        }
                    }
                }

                MPI_Bcast(send_matrix_MPI, send_size, MPI_DOUBLE, SENDER, MPI_COMM_WORLD);

                index_counter = 0;
                for (int j = 0; j < matsize; j++)
                {
                    for (int l = 0; l < matsize; l++)
                    {
                        recieve_matrix(j,l) = send_matrix_MPI[index_counter];
                        index_counter += 1;
                    }
                }


                for (int j = 0; j < matsize; j++)
                {
                    for (int l = 0; l < matsize; l++)
                    {
                        compact_mo(j,k,l) += c(i,a) * recieve_matrix(j,l); // Communication and distributed memory
                        //HartFock->Calc_Integrals_On_The_Fly(i,k,j,l); // on the fly calculation, not implemented in parallel jet. Should not be needed
                        //Integrals(Return_Integral_Index(i,j,k,l)); // compounded index storage
                    }
                }
            }
        }

        compact_mo2 = zeros(matsize, matsize, matsize);
        for (int b = 0; b < matsize; b++)
        {
            for (int i = 0; i < matsize; i++)
            {
                for (int j = 0; j < matsize; j++)
                {
                    for (int k = 0; k < matsize; k++)
                    {
                        compact_mo2(b,j,k) += c(i,b) * compact_mo(i,j,k);
                    }
                }
            }
        }

        compact_mo = zeros(matsize, matsize, matsize);
        for (int b = 0; b < matsize; b++)
        {
            for (int g = 0; g < matsize; g++)
            {
                for (int i = 0; i < matsize; i++)
                {
                    for (int j = 0; j < matsize; j++)
                    {
                        compact_mo(b,g,j) += c(i,g) * compact_mo2(b,i,j);
                    }
                }
            }
        }

        compact_mo2 = zeros(matsize, matsize, matsize);

        for (int b = 0; b < matsize; b++)
        {
            for (int g = 0; g < matsize; g++)
            {
                for (int h = 0; h < matsize; h++)
                {
                    for (int i = 0; i < matsize; i++)
                    {
                        compact_mo2(b,g,h) += c(i,h) * compact_mo(b,g,i);
                    }
                }
            }
        }

            A = a - Speed_Elec;

            for (int b = 0; b < Speed_Elec; b++)
            {
                // MOLeftovers one is not distributed, but it is relatively small
                for (int i = Speed_Elec; i < matsize; i++)
                {
                    for (int j = 0; j < Speed_Elec; j++)
                    {
                        MOLeftovers(a-Speed_Elec,i-Speed_Elec)(b,j) = compact_mo2(b,i,j);
                    }
                }
            }

            for (int i = a; i < matsize; i++) // (0 -> Speed_Occ) is equal to (a -> matsize) - Speed_Elec
            {
                B = Speed_Occ - (i - Speed_Elec);
                INDEX_CHECK = A*Speed_Occ+B;

                if (INDEX_CHECK % size == rank) // Only the threads who NEED the memory will store it, this is memory distribution :-D
                {
                    for (int b = 0; b < matsize; b++)
                    {
                        for (int j = 0; j < matsize; j++)
                        {
                            MO9(a-Speed_Elec,i-Speed_Elec)(b,j) = compact_mo2(b,i,j);
                        }
                    }
                }
            }

            for (int j = 0; j < Speed_Elec; j++) // <-- j
            {
                INDEX_CHECK = A*Speed_Elec + Speed_Elec - j;
                if (INDEX_CHECK % size == rank)
                {
                    for (int b = Speed_Elec; b < matsize; b++)
                    {
                        for (int i = 0; i < Speed_Elec; i++)
                        {
                            MO2(a-Speed_Elec,j)(b-Speed_Elec,i) = compact_mo2(b,i,j); // <-- different indexing vs MO10 & MO3
                        }
                    }
                }
            }

            for (int i = 0; i < Speed_Elec; i++)
            {
                INDEX_CHECK = A*Speed_Elec + Speed_Elec - i;
                if (INDEX_CHECK % size == rank)
                {
                    for (int b = Speed_Elec; b < matsize; b++)
                    {
                        for (int j = 0; j < Speed_Elec; j++)
                        {
                            MO3(a-Speed_Elec,i)(j,b-Speed_Elec) = compact_mo2(i,j,b); // <-- different indexing vs MO10 & MO2
                        }
                    }
                }
            }

            for (int i = 0; i < Speed_Elec; i++) // <-- i
            {
                INDEX_CHECK = A*Speed_Elec + Speed_Elec - i;
                if (INDEX_CHECK % size == rank)
                {
                    for (int b = 0; b < matsize; b++)
                    {
                        for (int j = 0; j < matsize; j++)
                        {
                            MO10(a-Speed_Elec,i)(b,j) = compact_mo2(b,i,j); // <-- different indexing vs MO2 & MO3
                        }
                    }
                }
            }

    }

    compact_mo = zeros(1,1,1); // Remove this memory, .clear didnt work, needs more armadillo knowledge :-O
    compact_mo2 = zeros(1,1,1); // Remove this memory, .clear didnt work, just put it to 1x1x1 for now :-O

    //Integrals.clear();
    c.clear();
}

long int CCSD_Memory_optimized::Return_Integral_Index(int a, int b, int i, int j)
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

mat CCSD_Memory_optimized::Fill_FS(vec eigval)
{
    mat fss = zeros(Matrix_Size, Matrix_Size);

    for (int i = 0; i < Matrix_Size; i++)
    {
        fss(i,i) = eigval(i/2); // Diagonal matrise
    }

    return fss;
}

void CCSD_Memory_optimized::Fill_integ2_2D(int a, int i)
{
    int B,J;

    if (a%2 == 0 && i%2==0)
    {
        B = 0;
        for (int b = 0; b < unocc_orb; b++)
        {
            J = 0;
            for (int j = 0; j < n_Electrons; j++)
            {
                integ2_2D(B,J) = -MOLeftovers(a/2, b/2)(j/2, i/2)+MOLeftovers(a/2, b/2)(i/2, j/2);
                j++;
                J++;
            }
            b++;
            B++;
        }

        B = 0;
        for (int b = 1; b < unocc_orb; b++)
        {
            J = 0;
            for (int j = 1; j < n_Electrons; j++)
            {
                integ2_2D(B+Speed_Occ,J) = MOLeftovers(a/2, b/2)(i/2, j/2);
                j++;
                J++;
            }
            b++;
            B++;
        }

    }

    else if (a%2 == 0 && i%2==1)
    {
        B = 0;
        for (int b = 0; b < unocc_orb; b++)
        {
            J = 0;
            for (int j = 0; j < n_Electrons; j++)
            {
                integ2_2D(B,J) = 0;
                j++;
                J++;
            }
            b++;
            B++;
        }

        B = 0;
        for (int b = 1; b < unocc_orb; b++)
        {
            J = 0;
            for (int j = 1; j < n_Electrons; j++)
            {
                integ2_2D(B+Speed_Occ,J) = -MOLeftovers(a/2, b/2)(j/2, i/2);
                j++;
                J++;
            }
            b++;
            B++;
        }

    }

    else if (a%2 == 1 && i%2==0)
    {
        B = 0;
        for (int b = 0; b < unocc_orb; b++)
        {
            J = 0;
            for (int j = 0; j < n_Electrons; j++)
            {
                integ2_2D(B,J) = -MOLeftovers(a/2, b/2)(j/2, i/2);
                j++;
                J++;
            }
            b++;
            B++;
        }

        B = 0;
        for (int b = 1; b < unocc_orb; b++)
        {
            J = 0;
            for (int j = 1; j < n_Electrons; j++)
            {
                integ2_2D(B+Speed_Occ,J) = 0;
                j++;
                J++;
            }
            b++;
            B++;
        }

    }

    else
    {
        B = 0;
        for (int b = 0; b < unocc_orb; b++)
        {
            J = 0;
            for (int j = 0; j < n_Electrons; j++)
            {
                integ2_2D(B,J) = MOLeftovers(a/2, b/2)(i/2, j/2);
                j++;
                J++;
            }
            b++;
            B++;
        }

        B = 0;
        for (int b = 1; b < unocc_orb; b++)
        {
            J = 0;
            for (int j = 1; j < n_Electrons; j++)
            {
                integ2_2D(B+Speed_Occ,J) = MOLeftovers(a/2, b/2)(i/2, j/2) - MOLeftovers(a/2, b/2)(j/2,i/2);
                j++;
                J++;
            }
            b++;
            B++;
        }

    }
}

void CCSD_Memory_optimized::Fill_integ3_2D(int a, int i)
{
    int B,J;
    double val1, val2;

    int A = a/2;
    int I = i/2;

    if (a%2 == 0 && i%2 == 0)
    {
        B = 0;
        for (int b = 0; b < unocc_orb; b++)
        {
            J = 0;
            for (int j = 0; j < unocc_orb; j++)
            {
                integ3_2D(B,J) = MO9(A, I)(B+Speed_Elec, (J+Speed_Elec)) - MO9(A, I)(J+Speed_Elec, (B+Speed_Elec));
                j++;
                J++;
            }
            b++;
            B++;
        }

        // Lower half of matrix = 0, apparantly slower
        //integ3_2D(span(Speed_Occ, unocc_orb-1), span()) = Zero_Matrix;

        B = 0;
        for (int b = 1; b < unocc_orb; b++)
        {
            J = 0;
            for (int j = 1; j < unocc_orb; j++)
            {
                integ3_2D(B+Speed_Occ,J) = 0;
                j++;
                J++;
            }
            b++;
            B++;
        }


    }

    else if (a%2 == 1 && i%2 == 1)
    {
        // Upper half of matrix = 0
        //integ3_2D(span(0, Speed_Occ-1), span()) = Zero_Matrix;


        B = 0;
        for (int b = 0; b < unocc_orb; b++)
        {
            J = 0;
            for (int j = 0; j < unocc_orb; j++)
            {
                val1 = EqualFunc(a%2, b%2) * EqualFunc(i%2, j%2) *
                        MO9(A, I)(B+Speed_Elec, (J+Speed_Elec));

                val2 = EqualFunc(a%2, j%2) * EqualFunc(i%2, b%2) *
                        MO9(A, I)(j/2+Speed_Elec, (B+Speed_Elec));

                integ3_2D(B,J) = val1-val2;
                j++;
                J++;
            }
            b++;
            B++;
        }


        B = 0;
        for (int b = 1; b < unocc_orb; b++)
        {
            J = 0;
            for (int j = 1; j < unocc_orb; j++)
            {
                integ3_2D(B+Speed_Occ,J) =  MO9(A, I)(B+Speed_Elec, (J+Speed_Elec)) - MO9(A, I)(J+Speed_Elec, (B+Speed_Elec));
                j++;
                J++;
            }
            b++;
            B++;
        }
    }

    else if (a%2 == 0 && i%2 == 1)
    {
        B = 0;
        for (int b = 0; b < unocc_orb; b++)
        {
            J = 0;
            for (int j = 1; j < unocc_orb; j++)
            {
                integ3_2D(B,J) = MO9(A, I)(B+Speed_Elec, (J+Speed_Elec));
                j++;
                J++;
            }
            b++;
            B++;
        }

        B = 0;
        for (int b = 1; b < unocc_orb; b++)
        {
            J = 0;
            for (int j = 0; j < unocc_orb; j++)
            {
                integ3_2D(B+Speed_Occ,J) = -MO9(A, I)(J+Speed_Elec, (B+Speed_Elec));
                j++;
                J++;
            }
            b++;
            B++;
        }
    }

    else
    {
        B = 0;
        for (int b = 0; b < unocc_orb; b++)
        {
            J = 0;
            for (int j = 1; j < unocc_orb; j++)
            {
                integ3_2D(B,J) = -MO9(A, I)(J+Speed_Elec, (B+Speed_Elec));
                j++;
                J++;
            }
            b++;
            B++;
        }

        B = 0;
        for (int b = 1; b < unocc_orb; b++)
        {
            J = 0;
            for (int j = 0; j < unocc_orb; j++)
            {
                integ3_2D(B+Speed_Occ,J) = MO9(A, I)(B+Speed_Elec, (J+Speed_Elec));
                j++;
                J++;
            }
            b++;
            B++;
        }
    }
}

void CCSD_Memory_optimized::Fill_integ4_2D(int a, int i)
{
    int B,J;
    int A = a/2;
    int I = i/2;

    if (a%2 == 0 && i%2 == 0)
    {
        B = 0;
        for (int b = 0; b < unocc_orb; b++)
        {
            J = 0;
            for (int j = 0; j < unocc_orb; j++)
            {
                integ4_2D(B,J) = MO4(A, I)(B, J) - MO4(A, I)(J, B);
                j++;
                J++;
            }
            b++;
            B++;
        }


       // integ4_2D(span(Speed_Occ, unocc_orb-1), span()) = Zero_Matrix;


        B = 0;
        for (int b = 1; b < unocc_orb; b++)
        {
            J = 0;
            for (int j = 1; j < unocc_orb; j++)
            {
                integ4_2D(B+Speed_Occ,J) = 0;
                j++;
                J++;
            }
            b++;
            B++;
        }

    }

    else if (a%2 == 1 && i%2 == 0)
    {
        B = 0;
        for (int b = 0; b < unocc_orb; b++)
        {
            J = 0;
            for (int j = 0; j < unocc_orb; j++)
            {
                integ4_2D(B,J) = - MO4(A, I)(J, B);
                j++;
                J++;
            }
            b++;
            B++;
        }

        B = 0;
        for (int b = 1; b < unocc_orb; b++)
        {
            J = 0;
            for (int j = 1; j < unocc_orb; j++)
            {
                integ4_2D(B+Speed_Occ,J) = MO4(A, I)(B, J);
                j++;
                J++;
            }
            b++;
            B++;
        }
    }

    else if (a%2 == 1 && i%2 == 1)
    {
        //integ4_2D(span(0, Speed_Occ-1),span()) = Zero_Matrix;

        B = 0;
        for (int b = 0; b < unocc_orb; b++)
        {
            J = 0;
            for (int j = 0; j < unocc_orb; j++)
            {
                integ4_2D(B,J) = 0;
                j++;
                J++;
            }
            b++;
            B++;
        }

        B = 0;
        for (int b = 1; b < unocc_orb; b++)
        {
            J = 0;
            for (int j = 1; j < unocc_orb; j++)
            {
                integ4_2D(B+Speed_Occ,J) = MO4(A, I)(B, J) - MO4(A, I)(J, B);
                j++;
                J++;
            }
            b++;
            B++;
        }
    }

    else
    {
        B = 0;
        for (int b = 0; b < unocc_orb; b++)
        {
            J = 0;
            for (int j = 0; j < unocc_orb; j++)
            {
                integ4_2D(B,J) = MO4(A, I)(B, J);
                j++;
                J++;
            }
            b++;
            B++;
        }

        B = 0;
        for (int b = 1; b < unocc_orb; b++)
        {
            J = 0;
            for (int j = 1; j < unocc_orb; j++)
            {
                integ4_2D(B+Speed_Occ,J) = - MO4(A, I)(J, B);
                j++;
                J++;
            }
            b++;
            B++;
        }
    }
}

void CCSD_Memory_optimized::Fill_integ5_2D(int a, int i)
{
    int B,J;
    if (a%2 == 0 && i%2 == 0)
    {
        B = 0;
        for (int b = 0; b < unocc_orb; b++)
        {
            J = 0;
            for (int j = 0; j < unocc_orb; j++)
            {
                integ5_2D(B,J) =  MO10(a/2, i/2)(b/2+Speed_Elec, j/2+Speed_Elec)-MO10(a/2, i/2)(j/2+Speed_Elec, b/2+Speed_Elec);
                j++;
                J++;
            }
            b++;
            B++;
        }

        B = 0;
        for (int b = 1; b < unocc_orb; b++)
        {
            J = 0;
            for (int j = 1; j < unocc_orb; j++)
            {
                integ5_2D(B+Speed_Occ,J) = 0;
                j++;
                J++;
            }
            b++;
            B++;
        }
    }

    else if (a%2 == 1 && i%2 == 1)
    {
        B = 0;
        for (int b = 0; b < unocc_orb; b++)
        {
            J = 0;
            for (int j = 0; j < unocc_orb; j++)
            {
                integ5_2D(B,J) = 0;
                j++;
                J++;
            }
            b++;
            B++;
        }

        B = 0;
        for (int b = 1; b < unocc_orb; b++)
        {
            J = 0;
            for (int j = 1; j < unocc_orb; j++)
            {
                integ5_2D(B+Speed_Occ,J) = MO10(a/2, i/2)(b/2+Speed_Elec, j/2+Speed_Elec)-MO10(a/2, i/2)(j/2+Speed_Elec, b/2+Speed_Elec);
                j++;
                J++;
            }
            b++;
            B++;
        }
    }

    else if (a%2 == 1 && i%2 == 0)
    {
        B = 0;
        for (int b = 0; b < unocc_orb; b++)
        {
            J = 0;
            for (int j = 1; j < unocc_orb; j++)
            {
                integ5_2D(B,J) = -MO10(a/2, i/2)(j/2+Speed_Elec, b/2+Speed_Elec);
                j++;
                J++;
            }
            b++;
            B++;
        }

        B = 0;
        for (int b = 1; b < unocc_orb; b++)
        {
            J = 0;
            for (int j = 0; j < unocc_orb; j++)
            {
                integ5_2D(B+Speed_Occ,J) = MO10(a/2, i/2)(b/2+Speed_Elec, j/2+Speed_Elec);
                j++;
                J++;
            }
            b++;
            B++;
        }
    }

    else
    {
        B = 0;
        for (int b = 0; b < unocc_orb; b++)
        {
            J = 0;
            for (int j = 1; j < unocc_orb; j++)
            {
                integ5_2D(B,J) = MO10(a/2, i/2)(b/2+Speed_Elec, j/2+Speed_Elec);
                j++;
                J++;
            }
            b++;
            B++;
        }

        B = 0;
        for (int b = 1; b < unocc_orb; b++)
        {
            J = 0;
            for (int j = 0; j < unocc_orb; j++)
            {
                integ5_2D(B+Speed_Occ,J) = -MO10(a/2, i/2)(j/2+Speed_Elec, b/2+Speed_Elec);
                j++;
                J++;
            }
            b++;
            B++;
        }
    }
}

void CCSD_Memory_optimized::Fill_integ6_2D(int a, int i)
{
    int B,J;
    int A = a/2;
    int I = i/2;

    if (a%2 == 0 && i%2 == 0)
    {
        B = 0;
        for (int b = 0; b < unocc_orb; b++)
        {
            J = 0;
            for (int j = 0; j < n_Electrons; j++)
            {
                integ6_2D(B,J) = MO6_1(A, I)(B, J) - MO6_2(A, I)(J, B);
                j++;
                J++;
            }
            b++;
            B++;
        }

        B = 0;
        for (int b = 1; b < unocc_orb; b++)
        {
            J = 0;
            for (int j = 1; j < n_Electrons; j++)
            {
                integ6_2D(B+Speed_Occ,J) = 0;
                j++;
                J++;
            }
            b++;
            B++;
        }
    }

    else if (a%2 == 1 && i%2 == 0)
    {
        B = 0;
        for (int b = 0; b < unocc_orb; b++)
        {
            J = 0;
            for (int j = 0; j < n_Electrons; j++)
            {
                integ6_2D(B,J) = -MO6_2(A, I)(J, B);
                j++;
                J++;
            }
            b++;
            B++;
        }

        B = 0;
        for (int b = 1; b < unocc_orb; b++)
        {
            J = 0;
            for (int j = 1; j < n_Electrons; j++)
            {
                integ6_2D(B+Speed_Occ,J) = MO6_1(A, I)(B, J);
                j++;
                J++;
            }
            b++;
            B++;
        }
    }

    else if (a%2 == 1 && i%2 == 1)
    {
        B = 0;
        for (int b = 0; b < unocc_orb; b++)
        {
            J = 0;
            for (int j = 0; j < n_Electrons; j++)
            {
                integ6_2D(B,J) = 0;
                j++;
                J++;
            }
            b++;
            B++;
        }

        B = 0;
        for (int b = 1; b < unocc_orb; b++)
        {
            J = 0;
            for (int j = 1; j < n_Electrons; j++)
            {
                integ6_2D(B+Speed_Occ,J) = MO6_1(A, I)(B, J)-MO6_2(A, I)(J, B);
                j++;
                J++;
            }
            b++;
            B++;
        }
    }


    else
    {
        B = 0;
        for (int b = 0; b < unocc_orb; b++)
        {
            J = 0;
            for (int j = 0; j < n_Electrons; j++)
            {
                integ6_2D(B,J) = MO6_1(A, I)(B, J);
                j++;
                J++;
            }
            b++;
            B++;
        }

        B = 0;
        for (int b = 1; b < unocc_orb; b++)
        {
            J = 0;
            for (int j = 1; j < n_Electrons; j++)
            {
                integ6_2D(B+Speed_Occ,J) = -MO6_2(A, I)(J, B);
                j++;
                J++;
            }
            b++;
            B++;
        }
    }
}

void CCSD_Memory_optimized::Fill_integ7_2D(int a, int i)
{
    int B,J;
    if (a%2 == 1 && i%2 == 1)
    {
        B = 0;
        for (int b = 0; b < unocc_orb; b++)
        {
            J = 0;
            for (int j = 0; j < n_Electrons; j++)
            {
                integ7_2D(B+Speed_Occ,J) =  MO10(a/2, i/2)(b/2+Speed_Elec, j/2)-MO10(a/2, i/2)(j/2, b/2+Speed_Elec);
                j++;
                J++;
            }
            b++;
            B++;
        }

        B = 0;
        for (int b = 1; b < unocc_orb; b++)
        {
            J = 0;
            for (int j = 1; j < n_Electrons; j++)
            {
                integ7_2D(B,J) = 0;
                j++;
                J++;
            }
            b++;
            B++;
        }
    }

    else if (a%2 == 1 && i%2 == 0)
    {
        B = 0;
        for (int b = 0; b < unocc_orb; b++)
        {
            J = 0;
            for (int j = 1; j < n_Electrons; j++)
            {
                integ7_2D(B,J) =  -MO10(a/2, i/2)(j/2, b/2+Speed_Elec);
                j++;
                J++;
            }
            b++;
            B++;
        }

        B = 0;
        for (int b = 1; b < unocc_orb; b++)
        {
            J = 0;
            for (int j = 0; j < n_Electrons; j++)
            {
                integ7_2D(B+Speed_Occ,J) = MO10(a/2, i/2)(b/2+Speed_Elec, j/2);
                j++;
                J++;
            }
            b++;
            B++;
        }
    }

    else if (a%2==0 && i%2==1)
    {
        B = 0;
        for (int b = 0; b < unocc_orb; b++)
        {
            J = 0;
            for (int j = 1; j < n_Electrons; j++)
            {
                integ7_2D(B+Speed_Occ,J) =  -MO10(a/2, i/2)(j/2, b/2+Speed_Elec);
                j++;
                J++;
            }
            b++;
            B++;
        }

        B = 0;
        for (int b = 1; b < unocc_orb; b++)
        {
            J = 0;
            for (int j = 0; j < n_Electrons; j++)
            {
                integ7_2D(B,J) = MO10(a/2, i/2)(b/2+Speed_Elec, j/2);
                j++;
                J++;
            }
            b++;
            B++;
        }
    }

    else if (a%2 == 0 && i%2 == 0)
    {
        B = 0;
        for (int b = 0; b < unocc_orb; b++)
        {
            J = 0;
            for (int j = 1; j < n_Electrons; j++)
            {
                integ7_2D(B,J) =  MO10(a/2, i/2)(b/2+Speed_Elec, j/2)-MO10(a/2, i/2)(j/2, b/2+Speed_Elec);
                j++;
                J++;
            }
            b++;
            B++;
        }

        B = 0;
        for (int b = 1; b < unocc_orb; b++)
        {
            J = 0;
            for (int j = 0; j < n_Electrons; j++)
            {
                integ7_2D(B+Speed_Occ,J) = 0;
                j++;
                J++;
            }
            b++;
            B++;
        }
    }

    else
    {
        B = 0;
        for (int b = 0; b < unocc_orb; b++)
        {
            J = 0;
            for (int j = 1; j < n_Electrons; j++)
            {
                integ7_2D(B,J) =  -MO10(a/2, i/2)(j/2, b/2+Speed_Elec);
                j++;
                J++;
            }
            b++;
            B++;
        }

        B = 0;
        for (int b = 1; b < unocc_orb; b++)
        {
            J = 0;
            for (int j = 0; j < n_Electrons; j++)
            {
                integ7_2D(B+Speed_Occ,J) = MO10(a/2, i/2)(b/2+Speed_Elec, j/2);
                j++;
                J++;
            }
            b++;
            B++;
        }
    }
}

void CCSD_Memory_optimized::Fill_integ8_2D(int a, int i)
{
    int B,J;
    int A = a/2;
    int I = i/2;
    if (a%2 == 0 && i%2 == 0)
    {
        B = 0;
        for (int b = 0; b < n_Electrons; b++)
        {
            J = 0;
            for (int j = 0; j < n_Electrons; j++)
            {
                integ8_2D(B,J) = MO8(A, I)(B, J) - MO8(A, I)(J, B);
                j++;
                J++;
            }
            b++;
            B++;
        }

        B = 0;
        for (int b = 1; b < n_Electrons; b++)
        {
            J = 0;
            for (int j = 1; j < n_Electrons; j++)
            {
                integ8_2D(B+Speed_Elec,J) = 0;
                j++;
                J++;
            }
            b++;
            B++;
        }
    }

    else if(a%2 == 1 && i%2 == 1)
    {
        B = 0;
        for (int b = 0; b < n_Electrons; b++)
        {
            J = 0;
            for (int j = 0; j < n_Electrons; j++)
            {
                integ8_2D(B,J) = 0;
                j++;
                J++;
            }
            b++;
            B++;
        }

        B = 0;
        for (int b = 1; b < n_Electrons; b++)
        {
            J = 0;
            for (int j = 1; j < n_Electrons; j++)
            {
                integ8_2D(B+Speed_Elec,J) = MO8(A, I)(B, J) - MO8(A, I)(J, B);
                j++;
                J++;
            }
            b++;
            B++;
        }
    }

    else if (a%2 == 0 && i%2 == 1)
    {
        B = 0;
        for (int b = 0; b < n_Electrons; b++)
        {
            J = 0;
            for (int j = 1; j < n_Electrons; j++)
            {
                integ8_2D(B,J) = MO8(A, I)(B, J);
                j++;
                J++;
            }
            b++;
            B++;
        }

        B = 0;
        for (int b = 1; b < n_Electrons; b++)
        {
            J = 0;
            for (int j = 0; j < n_Electrons; j++)
            {
                integ8_2D(B+Speed_Elec,J) = -MO8(A, I)(J, B);
                j++;
                J++;
            }
            b++;
            B++;
        }
    }

    else
    {
        B = 0;
        for (int b = 0; b < n_Electrons; b++)
        {
            J = 0;
            for (int j = 1; j < n_Electrons; j++)
            {
                integ8_2D(B,J) = -MO8(A, I)(J, B);
                j++;
                J++;
            }
            b++;
            B++;
        }

        B = 0;
        for (int b = 1; b < n_Electrons; b++)
        {
            J = 0;
            for (int j = 0; j < n_Electrons; j++)
            {
                integ8_2D(B+Speed_Elec,J) = MO8(A, I)(B, J);
                j++;
                J++;
            }
            b++;
            B++;
        }
    }

}

void CCSD_Memory_optimized::Fill_integ9_2D(int a, int i)
{
    int B,J;
    double val1, val2;

    if ((a+i)%2 == 0)
    {
        B = 0;
        for (int b = 0; b < unocc_orb; b++)
        {
            J = 0;
            for (int j = 0; j < n_Electrons; j++)
            {
                val1 = EqualFunc(a%2, b%2) * EqualFunc(i%2, j%2) *
                        MO9(a/2, i/2)(b/2+Speed_Elec, (j/2));

                val2 = EqualFunc(a%2, j%2) * EqualFunc(i%2, b%2) *
                        MO9(a/2, i/2)(j/2, (b/2+Speed_Elec));

                integ9_2D(B,J) = val1-val2;
                j++;
                J++;
            }
            b++;
            B++;
        }

        B = 0;
        for (int b = 1; b < unocc_orb; b++)
        {
            J = 0;
            for (int j = 1; j < n_Electrons; j++)
            {
                val1 = EqualFunc(a%2, b%2) * EqualFunc(i%2, j%2) *
                        MO9(a/2, i/2)(b/2+Speed_Elec, (j/2));

                val2 = EqualFunc(a%2, j%2) * EqualFunc(i%2, b%2) *
                        MO9(a/2, i/2)(j/2, (b/2+Speed_Elec));

                integ9_2D(B+Speed_Occ,J) = val1-val2;
                j++;
                J++;
            }
            b++;
            B++;
        }
    }

    else
    {
        B = 0;
        for (int b = 0; b < unocc_orb; b++)
        {
            J = 0;
            for (int j = 1; j < n_Electrons; j++)
            {
                val1 = EqualFunc(a%2, b%2) * EqualFunc(i%2, j%2) *
                        MO9(a/2, i/2)(b/2+Speed_Elec, (j/2));

                val2 = EqualFunc(a%2, j%2) * EqualFunc(i%2, b%2) *
                        MO9(a/2, i/2)(j/2, (b/2+Speed_Elec));

                integ9_2D(B,J) = val1-val2;
                j++;
                J++;
            }
            b++;
            B++;
        }

        B = 0;
        for (int b = 1; b < unocc_orb; b++)
        {
            J = 0;
            for (int j = 0; j < n_Electrons; j++)
            {
                val1 = EqualFunc(a%2, b%2) * EqualFunc(i%2, j%2) *
                        MO9(a/2, i/2)(b/2+Speed_Elec, (j/2));

                val2 = EqualFunc(a%2, j%2) * EqualFunc(i%2, b%2) *
                        MO9(a/2, i/2)(j/2, (b/2+Speed_Elec));

                integ9_2D(B+Speed_Occ,J) = val1-val2;
                j++;
                J++;
            }
            b++;
            B++;
        }
    }
}

void CCSD_Memory_optimized::Fill_integ10_2D(int a, int i)
{
    double temp1, temp2;

    int B,J;
    if (a%2 == 0)
    {
        B = 0;
        for (int b = 0; b < unocc_orb; b++)
        {
            J = 0;
            for (int j = 0; j < n_Electrons; j++)
            {
                integ10_2D(B,J) = MO2(a/2, i/2)(b/2, j/2) - MO3(a/2, i/2)(j/2, b/2);
                j++;
                J++;
            }
            b++;
            B++;
        }

        B = 0;
        for (int b = 1; b < unocc_orb; b++)
        {
            J = 0;
            for (int j = 1; j < n_Electrons; j++)
            {
                integ10_2D(B+Speed_Occ,J) = - MO3(a/2, i/2)(j/2, b/2);
                j++;
                J++;
            }
            b++;
            B++;
        }
    }

    else
    {
        B = 0;
        for (int b = 0; b < unocc_orb; b++)
        {
            J = 0;
            for (int j = 1; j < n_Electrons; j++)
            {
                integ10_2D(B,J) = - MO3(a/2, i/2)(j/2, b/2);
                j++;
                J++;
            }
            b++;
            B++;
        }

        B = 0;
        for (int b = 1; b < unocc_orb; b++)
        {
            J = 0;
            for (int j = 0; j < n_Electrons; j++)
            {
                integ10_2D(B+Speed_Occ,J) = MO2(a/2, i/2)(b/2, j/2) - MO3(a/2, i/2)(j/2, b/2);
                j++;
                J++;
            }
            b++;
            B++;
        }
    }
}

void CCSD_Memory_optimized::Map_T2_For_MPI()
{
    // This function is designed to map out T2 new in a way such that the parallel implementation is effective
    // Also to ensure we do not store symmetric terms and pass symmetric terms through parallel

    int number_counter;
    int A,B;
    int INDEX_CHECK;

    T2_MPI = (double**) malloc(size * sizeof(double *));
    WORK_EACH_NODE = zeros(size);

    // Not superb work distribution, but this way ensures somewhat of a communication minimization
    // Communication is only one broadcast per thread per iteration this way,
    // since all work done by one node is stored in one array, and this is broadcasted in one MPI_Bcast call

    for (int K = 0; K < size; K++)
    {
        number_counter = 0;

        for (int a = 0; a < unocc_orb; a++)
        {
            A = a/2;
            for (int b = a+2; b < unocc_orb; b++)
            {

                B = Speed_Occ - b/2;
                INDEX_CHECK = A*Speed_Occ+B;

                if (INDEX_CHECK % size == K)
                {

                    for (int i = 0; i < n_Electrons; i++)
                    {
                        for (int j = i+2; j < n_Electrons; j++)
                        {

                            number_counter += 1;
                            j++;
                        }
                        i++;
                    }
                }

                b++;
            }
            a++;
        }

        for (int a = 1; a < unocc_orb; a++)
        {
            A = a/2;

            for (int b = a+2; b < unocc_orb; b++)
            {
                B = Speed_Occ - b/2;
                INDEX_CHECK = A*Speed_Occ+B;

                if (INDEX_CHECK % size == K)
                {
/*
                    for (int i = 0; i < n_Electrons; i++)
                    {
                        for (int j = i+2; j < n_Electrons; j++)
                        {

                            number_counter += 1;
                            j++;
                        }
                        i++;
                    }*/

                    for (int i = 1; i < n_Electrons; i++)
                    {
                        for (int j = i+2; j < n_Electrons; j++)
                        {
                            number_counter += 1;
                            j++;
                        }
                        i++;
                    }
                }

                b++;
            }
            a++;
        }

        for (int a = 1; a < unocc_orb; a++)
        {
            A = a/2;


            for (int b = a+1; b < unocc_orb; b++)
            {
                B = Speed_Occ - b/2;
                INDEX_CHECK = A*Speed_Occ+B;

                if (INDEX_CHECK % size == K)
                {

                    for (int i = 1; i < n_Electrons; i++)
                    {
                        for (int j = i+1; j < n_Electrons; j++)
                        {

                            number_counter += 1;
                            j++;
                        }
                        i++;
                    }


                    for (int i = 0; i < n_Electrons; i++)
                    {
                        for (int j = i+1; j < n_Electrons; j++)
                        {

                            number_counter += 1;
                            j++;
                        }
                        i++;
                    }
                }
                b++;
            }
            a++;
        }

        for (int a = 0; a < unocc_orb; a++)
        {
            A = a/2;


            for (int b = a+1; b < unocc_orb; b++)
            {
                B = Speed_Occ - b/2;
                INDEX_CHECK = A*Speed_Occ+B;

                if (INDEX_CHECK % size == K)
                {

                    for (int i = 1; i < n_Electrons; i++)
                    {
                        for (int j = i+1; j < n_Electrons; j++)
                        {

                            number_counter += 1;
                            j++;
                        }
                        i++;
                    }


                    for (int i = 0; i < n_Electrons; i++)
                    {
                        for (int j = i+1; j < n_Electrons; j++)
                        {

                            number_counter += 1;
                            j++;
                        }
                        i++;
                    }
                }
                b++;
            }
            a++;
        }

        for (int a = 0; a < unocc_orb; a++)
        {
            for (int i = 0; i < n_Electrons; i++)
            {
                INDEX_CHECK = a/2 * Speed_Elec + Speed_Elec - i/2;
                if (INDEX_CHECK % size == K)
                {
                    number_counter += 1;
                }

                i++;
            }
            a++;
        }

        for (int a = 1; a < unocc_orb; a++)
        {
            for (int i = 1; i < n_Electrons; i++)
            {
                INDEX_CHECK = a/2 * Speed_Elec + Speed_Elec - i/2;
                if (INDEX_CHECK % size == K)
                {
                    number_counter += 1;
                }
                i++;
            }
            a++;
        }




        if (rank==0)
        {
            cout << "Work : " << number_counter << " for node " << K << endl;
        }

        WORK_EACH_NODE(K) = number_counter;

        // Determine the size of nodes T2 calculations
        T2_MPI[K] = (double*) malloc(number_counter*sizeof(double));

    }
}

void CCSD_Memory_optimized::Fill_W4_MPI()
{

    int index_counter = 0;
    int INDEX_CHECK;

    for (int a = 0; a < unocc_orb; a++)
    {
        for (int m = 0; m < n_Electrons; m++)
        {
            INDEX_CHECK = a/2 * Speed_Elec + Speed_Elec - m/2;
            if (INDEX_CHECK % size == rank)
            {
                Fill_integ7_2D(a,m);
                Fill_integ5_2D(a,m);

                for (int e = 0; e < unocc_orb; e++)
                {
                    Fill_integ2_2D_even_even(e, m);
                    for(int i = 0; i < n_Electrons; i++)
                    {
                        Part1_MPI[rank][index_counter] = -integ7_2D(e/2,i/2) - accu(W_3.at(i,m)(e/2,span()) % T_1.row(a/2))
                                + accu(integ5_2D(span(0, Speed_Occ-1),e/2) % T_1(span(0, Speed_Occ-1),i/2))
                                + 0.5*accu(integ2_2D % t2.at(a,i));
                        index_counter += 1;
                        i++;
                    }
                    e++;
                }
            }
            m++;
        }
        a++;
    }

    for (int a = 1; a < unocc_orb; a++)
    {
        for (int m = 0; m < n_Electrons; m++)
        {
            INDEX_CHECK = a/2 *Speed_Elec + Speed_Elec - m/2;
            if (INDEX_CHECK % size == rank)
            {

                Fill_integ7_2D(a,m);
                Fill_integ5_2D(a,m);

                for (int e = 0; e < unocc_orb; e++)
                {
                    Fill_integ2_2D_even_even(e, m);
                    for (int i = 1; i < n_Electrons; i++)
                    {
                        Part1_MPI[rank][index_counter] = -integ7_2D(e/2,i/2) - accu(W_3.at(i,m)(e/2,span()) % T_1.row(a/2+Speed_Occ))
                                + accu(integ5_2D(span(Speed_Occ, Speed_Occ+Speed_Occ-1),e/2) % T_1(span(Speed_Occ, Speed_Occ+Speed_Occ-1),i/2))
                                + 0.5*accu(integ2_2D % t2.at(a,i));;
                        index_counter += 1;
                        i++;
                    }
                    e++;
                }

                for (int e = 1; e < unocc_orb; e++)
                {
                    Fill_integ2_2D_odd_even(e, m);
                    for (int i = 0; i < n_Electrons; i++)
                    {
                        Part1_MPI[rank][index_counter] = -integ7_2D(e/2+Speed_Occ,i/2) - accu(W_3.at(i,m)(e/2+Speed_Occ,span()) % T_1.row(a/2+Speed_Occ))
                                + accu(integ5_2D(span(0, Speed_Occ-1),e/2) % T_1(span(0, Speed_Occ-1),i/2))
                                + 0.5*accu(integ2_2D % t2.at(a,i));
                        index_counter += 1;
                        i++;
                    }
                    e++;
                }
            }
            m++;
        }
        a++;
    }

    for (int a = 0; a < unocc_orb; a++)
    {
        for (int m = 1; m < n_Electrons; m++)
        {
            INDEX_CHECK = a/2 *Speed_Elec + Speed_Elec - m/2;
            if (INDEX_CHECK % size == rank)
            {

                Fill_integ7_2D(a,m);
                Fill_integ5_2D(a,m);

                for (int e = 0; e < unocc_orb; e++)
                {
                    Fill_integ2_2D_even_odd(e, m);
                    for (int i = 1; i < n_Electrons; i++)
                    {
                        Part1_MPI[rank][index_counter] = -integ7_2D(e/2,i/2) - accu(W_3.at(i,m)(e/2,span()) % T_1.row(a/2))
                                + accu(integ5_2D(span(Speed_Occ, Speed_Occ+Speed_Occ-1),e/2) % T_1(span(Speed_Occ, Speed_Occ+Speed_Occ-1),i/2))
                                + 0.5*accu(integ2_2D(span(Speed_Occ, Speed_Occ+Speed_Occ-1), span()) % t2.at(a,i)(span(Speed_Occ, Speed_Occ+Speed_Occ-1), span()));
                        index_counter += 1;
                        i++;
                    }
                    e++;
                }

                for (int e = 1; e < unocc_orb; e++)
                {
                    Fill_integ2_2D_odd_odd(e, m);
                    for (int i = 0; i < n_Electrons; i++)
                    {
                        Part1_MPI[rank][index_counter] = -integ7_2D(e/2+Speed_Occ,i/2) - accu(W_3.at(i,m)(e/2,span()) % T_1.row(a/2+a%2*Speed_Occ))
                                + accu(integ5_2D(span(0, Speed_Occ-1),e/2) % T_1(span(0, Speed_Occ-1),i/2))
                                + 0.5*accu(integ2_2D % t2.at(a,i));
                        index_counter += 1;
                        i++;
                    }
                    e++;
                }
            }
            m++;
        }
        a++;
    }

    for (int a = 1; a < unocc_orb; a++)
    {
        for (int m = 1; m < n_Electrons; m++)
        {
            INDEX_CHECK = a/2 *Speed_Elec + Speed_Elec - m/2;
            if (INDEX_CHECK % size == rank)
            {

            Fill_integ7_2D(a,m);
            Fill_integ5_2D(a,m);
                for (int e = 0; e < unocc_orb; e++)
                {
                    Fill_integ2_2D_even_odd(e, m);
                    for (int i = 0; i < n_Electrons; i++)
                    {
                        Part1_MPI[rank][index_counter] = -integ7_2D(e/2,i/2) - accu(W_3.at(i,m)(e/2,span()) % T_1.row(a/2+Speed_Occ))
                                + 0.5*accu(integ2_2D(span(0, Speed_Occ-1), span()) % t2.at(a,i)(span(0,Speed_Occ-1), span()));
                        index_counter += 1;
                        i++;
                    }
                    e++;
                }

                for (int e = 1; e < unocc_orb; e++)
                {
                    Fill_integ2_2D_odd_odd(e, m);
                    for (int i = 1; i < n_Electrons; i++)
                    {
                        Part1_MPI[rank][index_counter] = -integ7_2D(e/2+Speed_Occ,i/2) - accu(W_3.at(i,m)(e/2+Speed_Occ,span()) % T_1.row(a/2+Speed_Occ))
                                + accu(integ5_2D(span(Speed_Occ, Speed_Occ+Speed_Occ-1),e/2) % T_1(span(Speed_Occ, Speed_Occ+Speed_Occ-1),i/2))
                                + 0.5*accu(integ2_2D % t2.at(a,i));
                        index_counter += 1;
                        i++;
                    }
                    e++;
                }
            }
            m++;
        }
        a++;
    }


    // Communicate
    int work;
    for (int X = 0; X < size; X++)
    {
        work = WORK_EACH_NODE_Part2(X);
        MPI_Bcast(Part1_MPI[X], work, MPI_DOUBLE, X, MPI_COMM_WORLD);
    }


}

void CCSD_Memory_optimized::Map_Part1_For_MPI()
{
    // Mapping of how much work goes where for part 1 of MPI
    // Mapping is done for communication minimization and this will only be called way before iterations start
    // How much work per node is actually calculated before any work is done, and will in future be implemented in a script

    // ------------------------------------------------------------------------------------------------------------------------
    // This way you can know how your scaling will be PRIOR to initiating a supercomputer run which will take meny hours
    // This will be awsome
    // ------------------------------------------------------------------------------------------------------------------------

    // This function will need some aditional terms, for initialization of W_1 and W_2 :O

    Part1_MPI = (double**)malloc(size*sizeof(double*));
    WORK_EACH_NODE_Part1 = zeros(size);
    WORK_EACH_NODE_Part2 = zeros(size);
    int index_counter;
    int index_counter2;
    int INDEX_CHECK;

    for (int K = 0; K < size; K++)
    {
        index_counter = 0;
        for (int k = 0; k < n_Electrons; k++)
        {
            for (int l = 1; l < k; l++)
            {
                INDEX_CHECK = k/2*Speed_Elec + Speed_Elec-l/2;

                if (INDEX_CHECK % size == K)
                {


                    for (int i = 0; i < n_Electrons; i++)
                    {
                        index_counter += 1;
                        i++;
                    }

                    for (int i = 0; i < n_Electrons; i++)
                    {
                        for (int e = 1; e < unocc_orb; e++)
                        {
                            index_counter += 1;
                            e++;
                        }
                        i++;
                    }

                    for (int i = 1; i < n_Electrons; i++)
                    {
                        for (int e = 0; e < unocc_orb; e++)
                        {
                            index_counter += 1;
                            e++;
                        }
                        i++;
                    }

                    for (int i = 0; i < n_Electrons; i++)
                    {
                        for (int j = i+1; j < n_Electrons; j++)
                        {
                            index_counter += 1;
                            j++;
                        }
                        i++;
                    }

                    for (int i = 1; i < n_Electrons; i++)
                    {
                        for (int j = i+1; j < n_Electrons; j++)
                        {
                            index_counter += 1;
                            j++;
                        }
                        i++;
                    }

                }

                l++;
            }

            for (int l = k+1; l < n_Electrons; l++)
            {
                INDEX_CHECK = k/2*Speed_Elec + Speed_Elec-l/2;

                if (INDEX_CHECK % size == K)
                {
                    for (int i = 0; i < n_Electrons; i++)
                    {
                        index_counter += 1;
                        i++;
                    }
                }

                l++;
            }
            k++;
        }

        for (int k = 1; k < n_Electrons; k++)
        {
            for (int l = 0; l < k; l++)
            {
                INDEX_CHECK = k/2*Speed_Elec + Speed_Elec-l/2;

                if (INDEX_CHECK % size == K)
                {

                    for (int i = 1; i < n_Electrons; i++)
                    {
                        index_counter += 1;
                        i++;
                    }

                    for (int i = 0; i < n_Electrons; i++)
                    {
                        for (int e = 1; e < unocc_orb; e++)
                        {
                            index_counter += 1;
                            e++;
                        }
                        i++;
                    }

                    for (int i = 1; i < n_Electrons; i++)
                    {
                        for (int e = 0; e < unocc_orb; e++)
                        {
                            index_counter += 1;
                            e++;
                        }
                        i++;
                    }

                    for (int i = 0; i < n_Electrons; i++)
                    {
                        for (int j = i+1; j < n_Electrons; j++)
                        {
                            index_counter += 1;
                            j++;
                        }
                        i++;
                    }

                    for (int i = 1; i < n_Electrons; i++)
                    {
                        for (int j = i+1; j < n_Electrons; j++)
                        {
                            index_counter += 1;
                            j++;
                        }
                        i++;
                    }
                }
                l++;
            }

            for (int l = k+1; l < n_Electrons; l++)
            {
                INDEX_CHECK = k/2*Speed_Elec + Speed_Elec-l/2;

                if (INDEX_CHECK % size == K)
                {
                    for (int i = 1; i < n_Electrons; i++)
                    {
                        index_counter += 1;
                        i++;
                    }
                }
                l++;
            }
            k++;
        }

        for (int k = 0; k < n_Electrons; k++)
        {
            for (int l = 0; l < k+1; l++)
            {
                INDEX_CHECK = k/2*Speed_Elec + Speed_Elec-l/2;

                if (INDEX_CHECK % size == K)
                {
                    for (int i = 0; i < n_Electrons; i++)
                    {
                        index_counter += 1;
                        i++;
                    }

                    for (int i = 0; i < n_Electrons; i++)
                    {
                        for (int e = 0; e < unocc_orb; e++)
                        {
                            index_counter += 1;
                            e++;
                        }
                        i++;
                    }

                    for (int i = 0; i < n_Electrons; i++)
                    {
                        for (int j = i+2; j < n_Electrons; j++)
                        {
                            index_counter += 1;
                            j++;
                        }
                        i++;
                    }
                }
                l++;
            }

            for (int l = k+2; l < n_Electrons; l++)
            {
                INDEX_CHECK = k/2*Speed_Elec + Speed_Elec-l/2;

                if (INDEX_CHECK % size == K)
                {
                    for (int i = 0; i < n_Electrons; i++)
                    {
                        index_counter += 1;
                        i++;
                    }
                }

                l++;
            }
            k++;
        }

        for (int k = 1; k < n_Electrons; k++)
        {
            for (int l = 1; l < k+1; l++)
            {
                INDEX_CHECK = k/2*Speed_Elec + Speed_Elec-l/2;

                if (INDEX_CHECK % size == K)
                {
                    for (int i = 0; i < n_Electrons; i++)
                    {
                        index_counter += 1;
                        i++;
                    }

                    for (int i = 1; i < n_Electrons; i++)
                    {
                        for (int e = 1; e < unocc_orb; e++)
                        {
                            index_counter += 1;
                            e++;
                        }
                        i++;
                    }

                    for (int i = 1; i < n_Electrons; i++)
                    {
                        for (int j = i+2; j < n_Electrons; j++)
                        {
                            index_counter += 1;
                            j++;
                        }
                        i++;
                    }
                }
                l++;
            }

            for (int l = k+2; l < n_Electrons; l++)
            {
                INDEX_CHECK = k/2*Speed_Elec + Speed_Elec-l/2;

                if (INDEX_CHECK % size == K)
                {
                    for (int i = 0; i < n_Electrons; i++)
                    {
                        index_counter += 1;
                        i++;
                    }
                }
                l++;
            }
            k++;
        }






        ///////////////////////////////////////
        ///////////////////////////////////////
        ///////////////////////////////////////
        ///////////////////////////////////////
        ///////////////////////////////////////
        ///////////////////////////////////////
        ///////////////////////////////////////
        ///////////////////////////////////////
        ///////////////////////////////////////


            for (int a = 0; a < unocc_orb; a++)
            {
                for (int m = 0; m < n_Electrons; m++)
                {
                    INDEX_CHECK = a/2*Speed_Elec + Speed_Elec - m/2;
                    if (INDEX_CHECK % size == K)
                    {

                        for (int e = 0; e < unocc_orb; e++)
                        {
                            index_counter += 1;
                            e++;
                        }

                        for (int i = 0; i < n_Electrons; i++)
                        {
                            for (int j = i+2; j < n_Electrons; j++)
                            {
                                index_counter += 1;
                                j++;
                            }
                            i++;
                        }
                    }
                    m++;
                }
                a++;
            }

            for (int a = 1; a < unocc_orb; a++)
            {
                for (int m = 1; m < n_Electrons; m++)
                {
                    INDEX_CHECK = a/2*Speed_Elec + Speed_Elec - m/2;
                    if (INDEX_CHECK % size == K)
                    {
                        for (int e = 1; e < unocc_orb; e++)
                        {
                            index_counter += 1;
                            e++;
                        }

                        for (int i = 0; i < n_Electrons; i++)
                        {
                            for (int j = i+2; j < n_Electrons; j++)
                            {
                                index_counter += 1;
                                j++;
                            }
                            i++;
                        }

                        for (int i = 1; i < n_Electrons; i++)
                        {
                            for (int j = i+2; j < n_Electrons; j++)
                            {
                                index_counter += 1;
                                j++;
                            }
                            i++;
                        }
                    }
                    m++;
                }
                a++;
            }

            for (int a = 1; a < unocc_orb; a++)
            {
                for (int m = 0; m < n_Electrons; m++)
                {
                    INDEX_CHECK = a/2*Speed_Elec + Speed_Elec - m/2;
                    if (INDEX_CHECK % size == K)
                    {
                        for (int e = 1; e < unocc_orb; e++)
                        {
                            index_counter += 1;
                            e++;
                        }

                        for (int i = 0; i < n_Electrons; i++)
                        {
                            for (int j = i+1; j < n_Electrons; j++)
                            {
                                index_counter += 1;
                                j++;
                            }
                            i++;
                        }

                        for (int i = 1; i < n_Electrons; i++)
                        {
                            for (int j = i+1; j < n_Electrons; j++)
                            {
                                index_counter += 1;
                                j++;
                            }
                            i++;
                        }
                    }
                    m++;
                }
                a++;
            }

            for (int a = 0; a < unocc_orb; a++)
            {
                for (int m = 1; m < n_Electrons; m++)
                {
                    INDEX_CHECK = a/2*Speed_Elec + Speed_Elec - m/2;
                    if (INDEX_CHECK % size == K)
                    {

                        for (int e = 0; e < unocc_orb; e++)
                        {
                            index_counter += 1;
                            e++;
                        }


                        for (int i = 0; i < n_Electrons; i++)
                        {
                            for (int j = i+1; j < n_Electrons; j++)
                            {
                                index_counter += 1;
                                j++;
                            }
                            i++;
                        }

                        for (int i = 1; i < n_Electrons; i++)
                        {
                            for (int j = i+1; j < n_Electrons; j++)
                            {
                                index_counter += 1;
                                j++;
                            }
                            i++;
                        }
                    }
                    m++;
                }
                a++;
            }


            for (int e = 0; e < unocc_orb; e++)
            {
                for (int m = 0; m < n_Electrons; m++)
                {
                    INDEX_CHECK = e/2*Speed_Elec + Speed_Elec - m/2;
                    if (INDEX_CHECK % size == K)
                    {
                        for (int i = 0; i < n_Electrons; i++)
                        {
                            index_counter += 1;
                            i++;
                        }
                    }
                    m++;
                }
            }

            for (int e = 0; e < unocc_orb; e++)
            {
                for (int m = 1; m < n_Electrons; m++)
                {
                    INDEX_CHECK = e/2*Speed_Elec + Speed_Elec - m/2;
                    if (INDEX_CHECK % size == K)
                    {
                        for (int i = 1; i < n_Electrons; i++)
                        {
                            index_counter += 1;
                            i++;
                        }
                    }
                    m++;
                }
            }

            for (int e = 0; e < unocc_orb; e++)
            {
                for (int m = 0; m < n_Electrons; m++)
                {
                    INDEX_CHECK = e/2 * Speed_Elec + Speed_Elec - m/2;
                    if (INDEX_CHECK % size == K)
                    {
                        index_counter += 1;
                        for (int a = 0; a < unocc_orb; a++)
                        {
                            index_counter += 1;
                            a++;
                        }
                    }
                    m++;
                }
                e++;
            }

            for (int e = 0; e < unocc_orb; e++)
            {
                for (int m = 1; m < n_Electrons; m++)
                {
                    INDEX_CHECK = e/2 * Speed_Elec + Speed_Elec - m/2;
                    if (INDEX_CHECK % size == K)
                    {
                        for (int a = 0; a < unocc_orb; a++)
                        {
                            index_counter += 1;
                            a++;
                        }
                    }
                    m++;
                }
                e++;
            }

            for (int e = 1; e < unocc_orb; e++)
            {
                for (int m = 0; m < n_Electrons; m++)
                {
                    INDEX_CHECK = e/2 * Speed_Elec + Speed_Elec - m/2;
                    if (INDEX_CHECK % size == K)
                    {
                        for (int a = 1; a < unocc_orb; a++)
                        {
                            index_counter += 1;
                            a++;
                        }
                    }
                    m++;
                }
                e++;
            }

            for (int e = 1; e < unocc_orb; e++)
            {
                for (int m = 1; m < n_Electrons; m++)
                {
                    INDEX_CHECK = e/2 * Speed_Elec + Speed_Elec - m/2;
                    if (INDEX_CHECK % size == K)
                    {
                        index_counter += 1;
                        for (int a = 1; a < unocc_orb; a++)
                        {
                            index_counter += 1;
                            a++;
                        }
                    }
                    m++;
                }
                e++;
            }



        ///////////////////////////////////////
        ///////////////////////////////////////
        ///////////////////////////////////////
        ///////////////////////////////////////
        ////////// T_1_new here ///////////////
        ///////////////////////////////////////
        ///////////////////////////////////////
        ///////////////////////////////////////
        ///////////////////////////////////////


            for (int a = 0; a < unocc_orb; a++)
            {
                for (int k = 0; k < n_Electrons; k++)
                {
                    INDEX_CHECK = a/2 * Speed_Elec + Speed_Elec - k/2;
                    if (INDEX_CHECK % size == K)
                    {
                        for (int i = 0; i < k; i++)
                        {
                            index_counter += 1;
                            i++;
                        }

                        // i will be less than k when k = odd number
                        // We want i to be an even number since a is even number, hence we start at i = k+1
                        for (int i = (k+1+(k+1)%2); i < n_Electrons; i++)
                        {
                            index_counter += 1;
                            i++;
                        }
                    }
                }
                a++;
            }

            for (int a = 1; a < unocc_orb; a++)
            {
                for (int k = 0; k < n_Electrons; k++)
                {
                    INDEX_CHECK = a/2 * Speed_Elec + Speed_Elec - k/2;
                    if (INDEX_CHECK % size == K)
                    {
                        for (int i = 1; i < k; i++)
                        {
                            index_counter += 1;
                            i++;
                        }

                        for (int i = k+1+(k)%2; i < n_Electrons; i++)
                        {
                            index_counter += 1;
                            i++;
                        }
                    }
                }
                a++;
            }





            if (rank == 0)
            {
                cout << "Work: " << index_counter << " for node: " << K << endl;
            }








            ///////////////////////////////////////
            ///////////////////////////////////////
            ///////////////////////////////////////
            ///////////////////////////////////////
            ////////// W_4 here ///////////////////
            ///////////////////////////////////////
            ///////////////////////////////////////
            ///////////////////////////////////////
            ///////////////////////////////////////

            index_counter2 = 0;

            for (int a = 0; a < unocc_orb; a++)
            {
                for (int m = 0; m < n_Electrons; m++)
                {
                    INDEX_CHECK = a/2 * Speed_Elec + Speed_Elec - m/2;
                    if (INDEX_CHECK % size == K)
                    {
                        for (int e = 0; e < unocc_orb; e++)
                        {
                            for(int i = 0; i < n_Electrons; i++)
                            {
                                index_counter2 += 1;
                                i++;
                            }
                            e++;
                        }
                    }
                    m++;
                }
                a++;
            }

            for (int a = 1; a < unocc_orb; a++)
            {
                for (int m = 0; m < n_Electrons; m++)
                {
                    INDEX_CHECK = a/2 * Speed_Elec + Speed_Elec - m/2;
                    if (INDEX_CHECK % size == K)
                    {

                        for (int e = 0; e < unocc_orb; e++)
                        {
                            for (int i = 1; i < n_Electrons; i++)
                            {
                                index_counter2 += 1;
                                i++;
                            }
                            e++;
                        }

                        for (int e = 1; e < unocc_orb; e++)
                        {
                            for (int i = 0; i < n_Electrons; i++)
                            {
                                index_counter2 += 1;
                                i++;
                            }
                            e++;
                        }
                    }
                    m++;
                }
                a++;
            }

            for (int a = 0; a < unocc_orb; a++)
            {
                for (int m = 1; m < n_Electrons; m++)
                {
                    INDEX_CHECK = a/2 * Speed_Elec + Speed_Elec - m/2;
                    if (INDEX_CHECK % size == K)
                    {
                        for (int e = 0; e < unocc_orb; e++)
                        {
                            for (int i = 1; i < n_Electrons; i++)
                            {
                                index_counter2 += 1;
                                i++;
                            }
                            e++;
                        }

                        for (int e = 1; e < unocc_orb; e++)
                        {
                            for (int i = 0; i < n_Electrons; i++)
                            {
                                index_counter2 += 1;
                                i++;
                            }
                            e++;
                        }
                    }
                    m++;
                }
                a++;
            }

            for (int a = 1; a < unocc_orb; a++)
            {
                for (int m = 1; m < n_Electrons; m++)
                {
                    INDEX_CHECK = a/2 * Speed_Elec + Speed_Elec - m/2;
                    if (INDEX_CHECK % size == K)
                    {
                        for (int e = 0; e < unocc_orb; e++)
                        {
                            for (int i = 0; i < n_Electrons; i++)
                            {
                                index_counter2 += 1;
                                i++;
                            }
                            e++;
                        }

                        for (int e = 1; e < unocc_orb; e++)
                        {
                            for (int i = 1; i < n_Electrons; i++)
                            {
                                index_counter2 += 1;
                                i++;
                            }
                            e++;
                        }
                    }
                    m++;
                }
                a++;
            }

            if (rank == 0)
            {
                //cout << "Work per node part 2: " << index_counter2 << " for rank " << K << endl;
            }


        WORK_EACH_NODE_Part1(K) = index_counter;
        WORK_EACH_NODE_Part2(K) = index_counter2;

        if (index_counter2 > index_counter)
        {
            index_counter = index_counter2;
        }

        Part1_MPI[K] = (double*)malloc(index_counter * sizeof(double));
    }
}

void CCSD_Memory_optimized::Distribute_Part1()
{
    int K, L, E, I, J;
    double temp;

    int index_counter;
    int INDEX_CHECK;

    for (int X = 0; X < size; X++)
    {
        index_counter = 0;

        K = 0;
        for (int k = 0; k < n_Electrons; k++)
        {
            L = 0;
            for (int l = 1; l < k; l++)
            {

                INDEX_CHECK = K*Speed_Elec+Speed_Elec-L;
                if (INDEX_CHECK % size == X)
                {

                    for (int i = 0; i < n_Electrons; i++)
                    {
                        D2(K,i/2) += Part1_MPI[X][index_counter];
                        i++;
                        index_counter += 1;
                    }

                    I = 0;
                    for (int i = 0; i < n_Electrons; i++)
                    {
                        E = 0;
                        for (int e = 1; e < unocc_orb; e++)
                        {
                            W_3(i,k)(E+Speed_Occ,L) = Part1_MPI[X][index_counter];
                            e++;
                            E++;
                            index_counter += 1;
                        }
                        i++;
                        I++;
                    }

                    I = 0;
                    for (int i = 1; i < n_Electrons; i++)
                    {
                        E = 0;
                        for (int e = 0; e < unocc_orb; e++)
                        {
                            W_3(i,k)(E,L) = Part1_MPI[X][index_counter];
                            e++;
                            E++;
                            index_counter += 1;
                        }
                        i++;
                        I++;
                    }

                    I = 0;
                    for (int i = 0; i < n_Electrons; i++)
                    {
                        J = I;
                        for (int j = i+1; j < n_Electrons; j++)
                        {
                            temp = Part1_MPI[X][index_counter];

                            W_1(i,j)(K,L) += temp;
                            W_1(i,j)(L+Speed_Elec,K) -= temp;
                            index_counter += 1;
                            j++;
                            J++;
                        }
                        i++;
                        I++;
                    }

                    I = 0;
                    for (int i = 1; i < n_Electrons; i++)
                    {
                        J = I+1;
                        for (int j = i+1; j < n_Electrons; j++)
                        {
                            temp = Part1_MPI[X][index_counter];

                            W_1(i,j)(K,L) += temp;
                            W_1(i,j)(L+Speed_Elec,K) -= temp;
                            index_counter += 1;
                            j++;
                            J++;
                        }
                        i++;
                        I++;
                    }
                }
                l++;
                L++;
            }

            for (int l = k+1; l < n_Electrons; l++)
            {
                INDEX_CHECK = K*Speed_Elec+Speed_Elec-L;
                if (INDEX_CHECK % size == X)
                {

                    for (int i = 0; i < n_Electrons; i++)
                    {
                        D2.at(K,i/2) -= Part1_MPI[X][index_counter];
                        i++;
                        index_counter += 1;
                    }

                    I = 0;
                    for (int i = 0; i < n_Electrons; i++)
                    {
                        E = 0;
                        for (int e = 1; e < unocc_orb; e++)
                        {
                            W_3(i,k)(E+Speed_Occ,L) = -W_3(i,l)(E,K);
                            e++;
                            E++;
                        }
                        i++;
                        I++;
                    }

                    // NEW
                    I = 0;
                    for (int i = 1; i < n_Electrons; i++)
                    {
                        E = 0;
                        for (int e = 0; e < unocc_orb; e++)
                        {
                            W_3(i,k)(E,L) = -W_3(i,l)(E,K);
                            e++;
                            E++;
                        }
                        i++;
                        I++;
                    }
                }

                l++;
                L++;
            }
            k++;
            K++;
        }

        K = 0;
        for (int k = 1; k < n_Electrons; k++)
        {
            L = 0;
            for (int l = 0; l < k; l++)
            {

                INDEX_CHECK = K*Speed_Elec+Speed_Elec-L;
                if (INDEX_CHECK % size == X)
                {

                    for (int i = 1; i < n_Electrons; i++)
                    {
                        D2(K+Speed_Elec,i/2) += Part1_MPI[X][index_counter];
                        i++;
                        index_counter += 1;
                    }

                    I = 0;
                    for (int i = 0; i < n_Electrons; i++)
                    {
                        E = 0;
                        for (int e = 1; e < unocc_orb; e++)
                        {
                            W_3(i,k)(E,L) = Part1_MPI[X][index_counter];
                            e++;
                            E++;
                            index_counter += 1;
                        }
                        i++;
                        I++;
                    }

                    I = 0;
                    for (int i = 1; i < n_Electrons; i++)
                    {
                        E = 0;
                        for (int e = 0; e < unocc_orb; e++)
                        {
                            W_3(i,k)(E,L) = Part1_MPI[X][index_counter];
                            e++;
                            E++;
                            index_counter += 1;
                        }
                        i++;
                        I++;
                    }

                    I = 0;
                    for (int i = 0; i < n_Electrons; i++)
                    {
                        J = I;
                        for (int j = i+1; j < n_Electrons; j++)
                        {
                            temp = Part1_MPI[X][index_counter];

                            W_1(i,j)(K+Speed_Elec,L) += temp;
                            W_1(i,j)(L,K) -= temp;
                            j++;
                            J++;
                            index_counter += 1;
                        }
                        i++;
                        I++;
                    }

                    I = 0;
                    for (int i = 1; i < n_Electrons; i++)
                    {
                        J = I+1;
                        for (int j = i+1; j < n_Electrons; j++)
                        {
                            temp = Part1_MPI[X][index_counter];

                            W_1(i,j)(K+Speed_Elec,L) += temp;
                            W_1(i,j)(L,K) -= temp;
                            j++;
                            index_counter += 1;
                            J++;
                        }
                        i++;
                        I++;
                    }
                }
                L++;
                l++;
            }

            for (int l = k+1; l < n_Electrons; l++)
            {

                INDEX_CHECK = K*Speed_Elec+Speed_Elec-L;
                if (INDEX_CHECK % size == X)
                {

                    for (int i = 1; i < n_Electrons; i++)
                    {
                        D2.at(K+Speed_Elec,i/2) -= Part1_MPI[X][index_counter];
                        i++;
                        index_counter += 1;
                    }

                    I = 0;
                    for (int i = 1; i < n_Electrons; i++)
                    {
                        E = 0;
                        for (int e = 0; e < unocc_orb; e++)
                        {
                            W_3(i,k)(E,L) = -W_3(i,l)(E,K);
                            e++;
                            E++;
                        }
                        i++;
                        I++;
                    }

                    I = 0;
                    for (int i = 0; i < n_Electrons; i++)
                    {
                        E = 0;
                        for (int e = 1; e < unocc_orb; e++)
                        {
                            W_3(i,k)(E,L) = -W_3(i,l)(E+Speed_Occ,K);
                            e++;
                            E++;
                        }
                        i++;
                        I++;
                    }
                }

                l++;
                L++;
            }
            k++;
            K++;
        }

        K = 0;
        for (int k = 0; k < n_Electrons; k++)
        {
            L = 0;
            for (int l = 0; l < k+1; l++)
            {

                INDEX_CHECK = K*Speed_Elec+Speed_Elec-L;
                if (INDEX_CHECK % size == X)
                {

                    for (int i = 0; i < n_Electrons; i++)
                    {
                        D2.at(K,i/2) += Part1_MPI[X][index_counter];
                        i++;
                        index_counter += 1;
                    }

                    I = 0;
                    for (int i = 0; i < n_Electrons; i++)
                    {
                        E = 0;
                        for (int e = 0; e < unocc_orb; e++)
                        {
                            W_3(i,k)(E,L) = Part1_MPI[X][index_counter];
                            W_3(i,l)(E,K) = -W_3(i,k)(E,L);
                            index_counter += 1;
                            e++;
                            E++;
                        }
                        i++;
                        I++;
                    }

                    I = 0;
                    for (int i = 0; i < n_Electrons; i++)
                    {
                        J = I+1;
                        for (int j = i+2; j < n_Electrons; j++)
                        {
                            temp = Part1_MPI[X][index_counter];

                            W_1(i,j)(K,L) += temp;
                            W_1(i,j)(L,K) -= temp;
                            index_counter += 1;
                            j++;
                            J++;
                        }
                        i++;
                        I++;
                    }
                }

                l++;
                L++;
            }

            for (int l = k+2; l < n_Electrons; l++)
            {
                INDEX_CHECK = K*Speed_Elec+Speed_Elec-L;
                if (INDEX_CHECK % size == X)
                {
                    for (int i = 0; i < n_Electrons; i++)
                    {
                        D2.at(K,i/2) -= Part1_MPI[X][index_counter];
                        i++;
                        index_counter += 1;
                    }
                }

                l++;
                L++;
            }
            k++;
            K++;
        }

        K = 0;
        for (int k = 1; k < n_Electrons; k++)
        {
            L = 0;
            for (int l = 1; l < k+1; l++)
            {
                INDEX_CHECK = K*Speed_Elec+Speed_Elec-L;
                if (INDEX_CHECK % size == X)
                {
                    for (int i = 0; i < n_Electrons; i++)
                    {
                        D2.at(K+Speed_Elec,i/2) += Part1_MPI[X][index_counter];
                        i++;
                        index_counter += 1;
                    }

                    I = 0;
                    for (int i = 1; i < n_Electrons; i++)
                    {
                        E = 0;
                        for (int e = 1; e < unocc_orb; e++)
                        {
                            W_3(i,k)(E+Speed_Occ,L) = Part1_MPI[X][index_counter];
                            W_3(i,l)(E+Speed_Occ,K) = -W_3(i,k)(E+Speed_Occ,L);
                            e++;
                            E++;
                            index_counter += 1;
                        }
                        i++;
                        I++;
                    }

                    I = 0;
                    for (int i = 1; i < n_Electrons; i++)
                    {
                        J = I+1;
                        for (int j = i+2; j < n_Electrons; j++)
                        {
                            temp = Part1_MPI[X][index_counter];

                            W_1(i,j)(K+Speed_Elec,L) += temp;
                            W_1(i,j)(L+Speed_Elec,K) -= temp;

                            index_counter += 1;
                            j++;
                            J++;
                        }
                        i++;
                        I++;
                    }
                }
                l++;
                L++;
            }

            for (int l = k+2; l < n_Electrons; l++)
            {
                INDEX_CHECK = K*Speed_Elec+Speed_Elec-L;
                if (INDEX_CHECK % size == X)
                {
                    for (int i = 0; i < n_Electrons; i++)
                    {
                        D2.at(K+Speed_Elec,i/2) -= Part1_MPI[X][index_counter];
                        i++;
                        index_counter += 1;
                    }
                }

                l++;
                L++;
            }
            k++;
            K++;
        }







        ////////////////////////////////////
        ////////////////////////////////////
        ///7////////////////////////////////
        /// W_2 PART HERE! /////////////////
        ////////////////////////////////////
        ////////////////////////////////////
        ////////////////////////////////////

        int A, M;

        // Optimized W_2 version
        A = 0;
        for (int a = 0; a < unocc_orb; a++)
        {
            M = 0;
            for (int m = 0; m < n_Electrons; m++)
            {
                INDEX_CHECK = a/2*Speed_Elec + Speed_Elec - m/2;
                if (INDEX_CHECK % size == X)
                {
                    for (int e = 0; e < unocc_orb; e++)
                    {
                        D3(a, e/2) += Part1_MPI[X][index_counter];
                        e++;
                        index_counter += 1;
                    }

                    I = 0;
                    for (int i = 0; i < n_Electrons; i++)
                    {
                        J = I+1;
                        for (int j = i+2; j < n_Electrons; j++)
                        {
                            W_2(i,j)(A,M) += Part1_MPI[X][index_counter];
                            index_counter += 1;
                            j++;
                            J++;
                        }
                        i++;
                        I++;
                    }
                }

                m++;
                M++;
            }
            a++;
            A++;
        }

        A = 0;
        for (int a = 1; a < unocc_orb; a++)
        {
            M = 0;
            for (int m = 1; m < n_Electrons; m++)
            {
                INDEX_CHECK = a/2*Speed_Elec + Speed_Elec - m/2;
                if (INDEX_CHECK % size == X)
                {
                    for (int e = 1; e < unocc_orb; e++)
                    {
                        D3(a, e/2) += Part1_MPI[X][index_counter];
                        index_counter += 1;
                        e++;
                    }

                    I = 0;
                    for (int i = 0; i < n_Electrons; i++)
                    {
                        J = I+1;
                        for (int j = i+2; j < n_Electrons; j++)
                        {

                            W_2(i,j)(A+Speed_Occ,M) += Part1_MPI[X][index_counter];
                            index_counter += 1;
                            j++;
                            J++;
                        }
                        i++;
                        I++;
                    }

                    I = 0;
                    for (int i = 1; i < n_Electrons; i++)
                    {
                        J = I+1;
                        for (int j = i+2; j < n_Electrons; j++)
                        {

                            W_2(i,j)(A+Speed_Occ,M) += Part1_MPI[X][index_counter];
                            index_counter += 1;
                            j++;
                            J++;
                        }
                        i++;
                        I++;
                    }
                }
                M++;
                m++;
            }
            a++;
            A++;
        }

        A = 0;
        for (int a = 1; a < unocc_orb; a++)
        {
            M = 0;
            for (int m = 0; m < n_Electrons; m++)
            {
                INDEX_CHECK = a/2*Speed_Elec + Speed_Elec - m/2;

                if (INDEX_CHECK % size == X)
                {
                    for (int e = 1; e < unocc_orb; e++)
                    {
                        D3(a, e/2) += Part1_MPI[X][index_counter];
                        index_counter += 1;
                        e++;
                    }

                    I = 0;
                    for (int i = 0; i < n_Electrons; i++)
                    {
                        J = I;
                        for (int j = i+1; j < n_Electrons; j++)
                        {

                            W_2(i,j)(A+Speed_Occ,M) += Part1_MPI[X][index_counter];
                            index_counter += 1;
                            j++;
                            J++;
                        }
                        i++;
                        I++;
                    }

                    I = 0;
                    for (int i = 1; i < n_Electrons; i++)
                    {
                        J = I+1;
                        for (int j = i+1; j < n_Electrons; j++)
                        {

                            W_2(i,j)(A+Speed_Occ,M) += Part1_MPI[X][index_counter];
                            index_counter += 1;
                            j++;
                            J++;
                        }
                        i++;
                        I++;
                    }
                }

                m++;
                M++;
            }
            a++;
            A++;
        }

        A = 0;
        for (int a = 0; a < unocc_orb; a++)
        {
            M = 0;
            for (int m = 1; m < n_Electrons; m++)
            {
                INDEX_CHECK = a/2*Speed_Elec + Speed_Elec - m/2;
                if (INDEX_CHECK % size == X)
                {
                    for (int e = 0; e < unocc_orb; e++)
                    {
                        D3(a,e/2) += Part1_MPI[X][index_counter];
                        index_counter += 1;
                        e++;
                    }

                    I = 0;
                    for (int i = 0; i < n_Electrons; i++)
                    {
                        J = I;
                        for (int j = i+1; j < n_Electrons; j++)
                        {

                            W_2(i,j)(A,M) += Part1_MPI[X][index_counter];
                            index_counter += 1;
                            j++;
                            J++;
                        }
                        i++;
                        I++;
                    }

                    I = 0;
                    for (int i = 1; i < n_Electrons; i++)
                    {
                        J = I+1;
                        for (int j = i+1; j < n_Electrons; j++)
                        {

                            W_2(i,j)(A,M) += Part1_MPI[X][index_counter];
                            index_counter += 1;
                            j++;
                            J++;
                        }
                        i++;
                        I++;
                    }
                }
                m++;
                M++;
            }
            a++;
            A++;
        }


        ////////////////////////////////////
        ////////////////////////////////////
        ///7////////////////////////////////
        /// W_2 PART ENDING HERE! //////////
        ////////////////////////////////////
        ////////////////////////////////////
        ////////////////////////////////////


        ////////////////////////////////////
        ////////////////////////////////////
        ////////////////////////////////////
        ////////////////////////////////////
        ////////////////////////////////////
        /// Add F2 and F3 parts here ///////
        ////////////////////////////////////
        ////////////////////////////////////
        ////////////////////////////////////
        ////////////////////////////////////
        ////////////////////////////////////


        for (int e = 0; e < unocc_orb; e++)
        {
            M = 0;
            for (int m = 0; m < n_Electrons; m++)
            {
                INDEX_CHECK = e/2*Speed_Elec + Speed_Elec - M;
                if (INDEX_CHECK % size == X)
                {
                    I = 0;
                    for (int i = 0; i < n_Electrons; i++)
                    {
                        D2(M,I) += Part1_MPI[X][index_counter];
                        index_counter += 1;
                        i++;
                        I++;
                    }
                }
                m++;
                M++;
            }
        }

        for (int e = 0; e < unocc_orb; e++)
        {
            M = 0;
            for (int m = 1; m < n_Electrons; m++)
            {
                INDEX_CHECK = e/2*Speed_Elec + Speed_Elec - M;
                if (INDEX_CHECK % size == X)
                {
                    I = 0;
                    for (int i = 1; i < n_Electrons; i++)
                    {
                        D2(M+Speed_Elec,I) += Part1_MPI[X][index_counter];
                        index_counter += 1;
                        i++;
                        I++;
                    }
                }
                m++;
                M++;
            }
        }

        for (int e = 0; e < unocc_orb; e++)
        {
            for (int m = 0; m < n_Electrons; m++)
            {
                INDEX_CHECK = e/2 * Speed_Elec + Speed_Elec - m/2;
                if (INDEX_CHECK % size == X)
                {
                    D1(e/2, m/2) += Part1_MPI[X][index_counter];
                    index_counter += 1;

                    for (int a = 0; a < unocc_orb; a++)
                    {
                        D3.at(a,e/2) -= Part1_MPI[X][index_counter];
                        index_counter += 1;
                        a++;
                    }
                }
                m++;
            }
            e++;
        }

        for (int e = 0; e < unocc_orb; e++)
        {
            for (int m = 1; m < n_Electrons; m++)
            {
                INDEX_CHECK = e/2 * Speed_Elec + Speed_Elec - m/2;
                if (INDEX_CHECK % size == X)
                {
                    for (int a = 0; a < unocc_orb; a++)
                    {
                        D3.at(a,e/2) -= Part1_MPI[X][index_counter];
                        index_counter += 1;
                        a++;
                    }
                }
                m++;
            }
            e++;
        }

        for (int e = 1; e < unocc_orb; e++)
        {
            for (int m = 0; m < n_Electrons; m++)
            {
                INDEX_CHECK = e/2 * Speed_Elec + Speed_Elec - m/2;
                if (INDEX_CHECK % size == X)
                {
                    for (int a = 1; a < unocc_orb; a++)
                    {
                        D3.at(a,e/2) -= Part1_MPI[X][index_counter];
                        index_counter += 1;
                        a++;
                    }
                }
                m++;
            }
            e++;
        }

        for (int e = 1; e < unocc_orb; e++)
        {
            for (int m = 1; m < n_Electrons; m++)
            {
                INDEX_CHECK = e/2 * Speed_Elec + Speed_Elec - m/2;
                if (INDEX_CHECK % size == X)
                {
                    D1(e/2+Speed_Occ, m/2) += Part1_MPI[X][index_counter];
                    index_counter += 1;

                    for (int a = 1; a < unocc_orb; a++)
                    {
                        D3.at(a,e/2) -= Part1_MPI[X][index_counter];
                        index_counter += 1;
                        a++;
                    }
                }
                m++;
            }
            e++;
        }

        ////////////////////////////////////
        ////////////////////////////////////
        ///7////////////////////////////////
        /// T_1_new here! //////////////////
        ////////////////////////////////////
        ////////////////////////////////////
        ////////////////////////////////////

        for (int a = 0; a < unocc_orb; a++)
        {
            for (int k = 0; k < n_Electrons; k++)
            {
                INDEX_CHECK = a/2 * Speed_Elec + Speed_Elec - k/2;
                if (INDEX_CHECK % size == X)
                {
                    for (int i = 0; i < k; i++)
                    {
                        T_1_new(a/2, i/2) += Part1_MPI[X][index_counter];
                        index_counter += 1;
                        i++;
                    }

                    // i will be less than k when k = odd number
                    // We want i to be an even number since a is even number, hence we start at i = k+1
                    for (int i = (k+1+(k+1)%2); i < n_Electrons; i++)
                    {
                        T_1_new(a/2, i/2) += Part1_MPI[X][index_counter];
                        index_counter += 1;
                        i++;
                    }
                }
            }
            a++;
        }

        for (int a = 1; a < unocc_orb; a++)
        {
            for (int k = 0; k < n_Electrons; k++)
            {
                INDEX_CHECK = a/2 * Speed_Elec + Speed_Elec - k/2;
                if (INDEX_CHECK % size == X)
                {
                    for (int i = 1; i < k; i++)
                    {
                        T_1_new(a/2+Speed_Occ, i/2) += Part1_MPI[X][index_counter];
                        index_counter += 1;
                        i++;
                    }

                    // i will be less than k when k = even number
                    // We want i to be an odd number since a is an odd number, hence we start at i = k+1
                    for (int i = k+1+(k)%2; i < n_Electrons; i++)
                    {
                        T_1_new(a/2+Speed_Occ, i/2) += Part1_MPI[X][index_counter];
                        index_counter += 1;
                        i++;
                    }
                }
            }
            a++;
        }


    }
}

void CCSD_Memory_optimized::Fill_integ2_2D_even_even(int a, int i)
{
    int B,J;
        B = 0;
        for (int b = 0; b < unocc_orb; b++)
        {
            J = 0;
            for (int j = 0; j < n_Electrons; j++)
            {
                integ2_2D(B,J) = -MOLeftovers(a/2, b/2)(j/2, i/2)+MOLeftovers(a/2, b/2)(i/2, j/2);
                j++;
                J++;
            }
            b++;
            B++;
        }
        B = 0;
        for (int b = 1; b < unocc_orb; b++)
        {
            J = 0;
            for (int j = 1; j < n_Electrons; j++)
            {
                integ2_2D(B+Speed_Occ,J) = MOLeftovers(a/2, b/2)(i/2, j/2);
                j++;
                J++;
            }
            b++;
            B++;
        }

}

void CCSD_Memory_optimized::Fill_integ2_2D_even_odd(int a, int i)
{

    int B,J;
        B = 0;
        for (int b = 0; b < unocc_orb; b++)
        {
            J = 0;
            for (int j = 0; j < n_Electrons; j++)
            {
                integ2_2D(B,J) = 0;
                j++;
                J++;
            }
            b++;
            B++;
        }

        B = 0;
        for (int b = 1; b < unocc_orb; b++)
        {
            J = 0;
            for (int j = 1; j < n_Electrons; j++)
            {
                integ2_2D(B+Speed_Occ,J) = -MOLeftovers(a/2, b/2)(j/2, i/2);
                j++;
                J++;
            }
            b++;
            B++;
        }

}

void CCSD_Memory_optimized::Fill_integ2_2D_odd_even(int a, int i)
{
    int B,J;
        B = 0;
        for (int b = 0; b < unocc_orb; b++)
        {
            J = 0;
            for (int j = 0; j < n_Electrons; j++)
            {
                integ2_2D(B,J) = -MOLeftovers(a/2, b/2)(j/2, i/2);
                j++;
                J++;
            }
            b++;
            B++;
        }
        B = 0;
        for (int b = 1; b < unocc_orb; b++)
        {
            J = 0;
            for (int j = 1; j < n_Electrons; j++)
            {
                integ2_2D(B+Speed_Occ,J) = 0;
                j++;
                J++;
            }
            b++;
            B++;
        }

}

void CCSD_Memory_optimized::Fill_integ2_2D_odd_odd(int a, int i)
{
    int B,J;
        B = 0;
        for (int b = 0; b < unocc_orb; b++)
        {
            J = 0;
            for (int j = 0; j < n_Electrons; j++)
            {
                integ2_2D(B,J) = MOLeftovers(a/2, b/2)(i/2, j/2);
                j++;
                J++;
            }
            b++;
            B++;
        }

        B = 0;
        for (int b = 1; b < unocc_orb; b++)
        {
            J = 0;
            for (int j = 1; j < n_Electrons; j++)
            {
                integ2_2D(B+Speed_Occ,J) = MOLeftovers(a/2, b/2)(i/2, j/2) - MOLeftovers(a/2, b/2)(j/2,i/2);
                j++;
                J++;
            }
            b++;
            B++;
        }

}

// Big program
