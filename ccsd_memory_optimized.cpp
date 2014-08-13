#include "ccsd_memory_optimized.h"

CCSD_Memory_optimized::CCSD_Memory_optimized(int n_N, vec zz, mat rr, string B_s, int n_Elec, int ran, int siz, Hartree_Fock_Solver *Hartfock, bool frez)
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

double CCSD_Memory_optimized::CCSD(double toler, bool print_stuff)
{
    /*
     *
     *          Optimization thoughts for programmer:
     *
     *
     *          Work distribution part 1 not optimal, this is common:
     *
                Node : 0 work 1: 11412
                Node : 1 work 1: 10680
                Node : 2 work 1: 10312
                Node : 3 work 1: 9580
                Node : 4 work 1: 11780
                Node : 5 work 1: 11412
                Node : 6 work 1: 10680
                Node : 7 work 1: 10312

                Work distribution part 2 and 3 totally optimal.

     *
     * */

    // General starting stuff. Initialize a few values
    bool continue_ccsd = true;
    double convergance_check;
    int Itermax = 1000;
    iter = 1;
    E_old = 0;
    E_new = 0;
    double E_HF;

    // Calculate the Hartree Fock Energy. Also find the coefficients C and number of basis functions
    E_HF = HartFock->get_Energy(toler, 0); // Calc hartree fock energy
    Matrix_Size = 2*HartFock->ReturnMatrixSize(); // This will be twice as big as in Hatree Fock, because of spin up and down
    unocc_orb = Matrix_Size - n_Electrons; // Number of unocupied orbitals
    mat temp_matrix;
    Speed_Elec = (int) n_Electrons/2; // Number of Correlated Electrons
    Speed_Occ = (int) unocc_orb/2; // Number of unoccopied Spacial Orbitals
    Zero_Matrix = zeros(Speed_Occ, Speed_Occ);
    jump = size*2; // Used in parallel implementation to avoid if tests

    // Figure out how many orbitals to freeze, not used jet
    int how_meny_orbs_to_freze = Freeze_How_Meny_Orbitals();

    if (print_stuff == true)
    {
        cout << "Number of Electrons for CCSD: " << n_Electrons - 2 * how_meny_orbs_to_freze << endl;
        cout << "Number of Virtuals for CCSD: " << unocc_orb << endl;
    }

    // Get coefficients
    c = HartFock->ReturnC();

    // Perform AO->MO transformation
    Prepear_AOs(how_meny_orbs_to_freze);

    // Get Fock Eigenvalues
    fs = Fill_FS(HartFock->return_eigval_F(), how_meny_orbs_to_freze);
    HartFock = NULL; // Remove this memory,  not needed anymore
    R.clear();
    Z.clear();

    // First we must do some initial mapping to determine how
    // the Parallel Implementation will work.
    // Need some displacements and worksizes for use in global communication later
    // These numbers are constant throughout iterations
    // Also allocate these arrays
    ccsd_non_iterative_part ccsd_init(unocc_orb, n_Electrons, rank, size);
    WORK_EACH_NODE = ccsd_init.return_Work_T2();
    WORK_EACH_NODE_Part1 = ccsd_init.return_Work_P1();
    WORK_EACH_NODE_Part2 = ccsd_init.return_Work_P2();
    Global_Displacement_1 = ccsd_init.ret_Global_Disp1();
    Global_Displacement_2 = ccsd_init.ret_Global_Disp2();
    Global_Worksize_1 = ccsd_init.ret_Global_Work1();
    Global_Worksize_2 = ccsd_init.ret_Global_Work2();
    Local_Displacement1 = Speed_Elec * Speed_Occ;
    Where_To_Start_Part2 = ccsd_init.Return_Start_Part2_Pos();

    Global_Displacement_1_1 = (int**)malloc(size*sizeof(int*));
    Global_Worksize_1_1 = (int**)malloc(size*sizeof(int*));
    Global_Displacement_2_1 = (int**)malloc(size*sizeof(int*));
    Global_Worksize_2_1 = (int**)malloc(size*sizeof(int*));

    for (int i = 0; i < size; i++)
    {
        Global_Displacement_1_1[i] = (int*)malloc(size*sizeof(int));
        Global_Worksize_1_1[i] = (int*)malloc(size*sizeof(int));

        Global_Displacement_2_1[i] = (int*)malloc(size*sizeof(int));
        Global_Worksize_2_1[i] = (int*)malloc(size*sizeof(int));

        for (int j = 0; j < size; j++)
        {
            Global_Worksize_1_1[i][j] = Global_Worksize_1[j][i];
            Global_Worksize_2_1[i][j] = Global_Worksize_2[j][i];
        }

        Global_Displacement_1_1[i][0] = 0;
        Global_Displacement_2_1[i][0] = 0;
        for (int j = 1; j < size; j++)
        {
            Global_Displacement_1_1[i][j] = Global_Displacement_1_1[i][j-1] + Global_Worksize_1_1[i][j-1];
            Global_Displacement_2_1[i][j] = Global_Displacement_2_1[i][j-1] + Global_Worksize_2_1[i][j-1];
        }

        for (int j = 0; j < size; j++)
        {
            // Number of bits, not number of doubles, will be required later in use of MPI_Alltoallw
            Global_Displacement_1_1[i][j] = Global_Displacement_1_1[i][j] * sizeof(double);
            Global_Displacement_1[i][j] = Global_Displacement_1[i][j] * sizeof(double);

            Global_Displacement_2_1[i][j] = Global_Displacement_2_1[i][j] * sizeof(double);
            Global_Displacement_2[i][j] = Global_Displacement_2[i][j] * sizeof(double);

        }
    }

    Index_Swapping_W_4 = (double****)malloc(unocc_orb * sizeof(double***));
    for (int i = 0; i < unocc_orb; i++)
    {
        Index_Swapping_W_4[i] = (double***)malloc(n_Electrons * sizeof(double**));
        for (int j = 0; j < n_Electrons; j++)
        {
            // INDEX_CHECK, distribute
            if ((i/2 * Speed_Elec + j/2)%size == rank)
            {
                Index_Swapping_W_4[i][j] = (double**)malloc(unocc_orb * sizeof(double*));
                for (int k = 0; k < unocc_orb; k++)
                {
                    Index_Swapping_W_4[i][j][k] = (double*)malloc(n_Electrons * sizeof(double));
                }
            }
        }
    }

    Work_Each_Node_T2_Parallel = (int*)malloc(size * sizeof(int));
    Displacement_Each_Node_T2_Parallel = (int*)malloc(size * sizeof(int));

    Displacement_Each_Node_part1_Parallel = (int*)malloc(size*sizeof(int));
    Work_Each_Node_part1_Parallel = (int*)malloc(size*sizeof(int));

    for (int i = 0; i < size; i++)
    {
        Work_Each_Node_T2_Parallel[i] = WORK_EACH_NODE(i);
        Work_Each_Node_part1_Parallel[i] = WORK_EACH_NODE_Part1(i);
    }

    Displacement_Each_Node_T2_Parallel[0] = 0;
    Displacement_Each_Node_part1_Parallel[0] = 0;
    for (int i = 1; i < size; i++)
    {
        Displacement_Each_Node_part1_Parallel[i] = Displacement_Each_Node_part1_Parallel[i-1] + Work_Each_Node_part1_Parallel[i-1];
        Displacement_Each_Node_T2_Parallel[i] = Displacement_Each_Node_T2_Parallel[i-1] + Work_Each_Node_T2_Parallel[i-1];
    }

    // Figure out exactly how large arrays need to be. We will reuse arrays
    // so we only need the largest value

    int max_work_part_1_mpi = max(max(WORK_EACH_NODE_Part1));
    int max_work_part_2_mpi = max(max(WORK_EACH_NODE_Part2));
    int max_work_t2_amplitudes = max(max(WORK_EACH_NODE));

    if (max_work_part_1_mpi < max_work_part_2_mpi)
    {
        max_work_part_1_mpi = max_work_part_2_mpi;
    }

    if (max_work_t2_amplitudes < max_work_part_1_mpi)
    {
        max_work_t2_amplitudes = max_work_part_1_mpi;
    }

    for (int i = 0; i < size; i++)
    {

        max_work_part_1_mpi = 0;

        for (int j = 0; j < size; j++)
        {
            max_work_part_1_mpi += Global_Worksize_1[i][j];
        }

        if (max_work_t2_amplitudes < max_work_part_1_mpi)
        {
            max_work_t2_amplitudes = max_work_part_1_mpi;
        }

        max_work_part_1_mpi = 0;

        for (int j = 0; j < size; j++)
        {
            max_work_part_1_mpi += Global_Worksize_1_1[i][j];
        }

        if (max_work_t2_amplitudes < max_work_part_1_mpi)
        {
            max_work_t2_amplitudes = max_work_part_1_mpi;
        }

        max_work_part_1_mpi = 0;
        for (int j = 0; j < size; j++)
        {
            max_work_part_1_mpi += Global_Worksize_2[i][j];
        }

        if (max_work_t2_amplitudes < max_work_part_1_mpi)
        {
            max_work_t2_amplitudes = max_work_part_1_mpi;
        }

        max_work_part_1_mpi = 0;
        for (int j = 0; j < size; j++)
        {
            max_work_part_1_mpi += Global_Worksize_2_1[i][j];
        }

        if (max_work_t2_amplitudes < max_work_part_1_mpi)
        {
            max_work_t2_amplitudes = max_work_part_1_mpi;
        }
    }

    if (max_work_t2_amplitudes < size)
    {
        max_work_t2_amplitudes = size;
    }

    double temp_t2 = 0;
    double temp_t22 = 0;

    for (int i = 0; i < size; i++)
    {
        temp_t2 += WORK_EACH_NODE(i);
        temp_t22 += WORK_EACH_NODE_Part1(i);
    }

    if (temp_t2 < temp_t22)
    {
        temp_t2 = temp_t22;
    }

    if (max_work_t2_amplitudes < temp_t2)
    {
        max_work_t2_amplitudes = temp_t2;
    }


    mpi_types_array = (MPI_Datatype*)malloc(size*sizeof(MPI_Datatype));
    for (int i = 0; i < size; i++)
    {
        mpi_types_array[i] = MPI_DOUBLE;
    }

    // Allocate the arrays we will use
    SHARED_INFO_MPI = (double*)malloc(max_work_t2_amplitudes * sizeof(double));
    MY_OWN_MPI = (double*)malloc(max_work_t2_amplitudes * sizeof(double));
    Nr_Parallel_Operations = 3;

    // Define variables used in the CCSD method:

    tau3.set_size(n_Electrons, n_Electrons); // 1/4 n_o^2 n_u^2 // Use this for hyperspeed, can be removed if memoryproblems

    // Optimized 4D arrays
    t2.set_size(unocc_orb, n_Electrons);    // 1/2 n_o^2 n_u^2     // Must be stored somehow, possible with memory distribution possibly if MO storage changed
    W_1.set_size(n_Electrons, n_Electrons); // 1/4 n_o^4        // Relatively small
    W_2.set_size(n_Electrons, n_Electrons); // 1/4 n_o^3 n_u    // Relatively small
    W_3.set_size(n_Electrons, n_Electrons); // 3/8 n_o^3 n_u    // Relatively small

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

    // Memory distribution arrays
    W_4.set_size(unocc_orb, n_Electrons);   // distributed...
    W_5.set_size(unocc_orb, unocc_orb);

    // Initialize our 4D arrays
    for (int i = 0; i < n_Electrons; i++)
    {
        temp_matrix = zeros(unocc_orb, Speed_Elec);

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
            // If this if test passes the W_3 values will be zero in half the matrix due
            // to spin. Can ignore storage of these values.
            // This is the only array were we utilize compact storage to its fullest.
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

            // Memory distributed
            if ((i/2 * Speed_Elec + j/2)%size == rank)
            {
                temp_matrix = zeros(unocc_orb, Speed_Elec);
                W_4(i,j) = temp_matrix;
            }
        }
    }

    for (int i = 0; i < unocc_orb; i++)
    {
        for (int j = i+1; j < unocc_orb; j++)
        {
            // Memory distributed
            if ((i/2 * Speed_Occ - Calc_sum_a_n(i/2) + j/2)%size == rank)
            {
                temp_matrix = zeros(n_Electrons, Speed_Elec);
                W_5(i,j) = temp_matrix; // This array can be removed in future implementation
            }
        }
    }

    // Coupled Cluster is an iterative process, meaning we make an initial guess of t1 and t2 and calculate the energy
    // Then update t1 and t2 and check for convergance in the energy

    if (print_stuff == true)
    {
        // If this is printed the system will work.
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
                    // Compact storage in its original non-simplified way
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

    // Measure time
    clock_t start;
    clock_t slutt;
    double time_measured;

    // Initial energy calculated
    Fill_tau();
    E_new = Calc_Energy(); // Starting energy, with t1=0 guess

    if (print_stuff == true)
    {
        cout << "Energi: " << E_new << " Steg: " << iter << endl;
    }

    // Synchronize threads before we start iterations, for accurate time per iteration measurement
    MPI_Barrier(MPI_COMM_WORLD);

    // Start iterations
    while (continue_ccsd == true)
    {
        start = clock();
        E_old = E_new;

        // Part 1 of parallel implementation, see text
        Fill_F1(); // Initialize everything in part 1
        Fill_W1_and_W3(); // Calculate everything of part 1
        Distribute_Part1(); // Distribute everything of part 1
        Fill_F2(); // These are the F1 * something else part, n^3 scaling
        Fill_F3(); // These are the F1 * something else part, n^3 scaling

        // Part 2 of parallel implementation, see text
        Fill_W4_MPI(); // Find W4
        Fill_W4(); // Redistribute W_4
        Fill_W5(); // Redistribute [P(ij) P(ab) W_4*T2], to keep everything distributed in memory

        // This is part 3 of parallel implementation
        Fill_t2_new(); // Find new amplitudes

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

    // Iterations finished, print result and return

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
            j++;
        }
        i++;
    }

    for (int i = 0; i < n_Electrons; i++)
    {
        for (int j = i+2; j < n_Electrons; j++)
        {
            Fill_integ4_2D(i,j);
            E1 += accu(tau3(i,j)(span(0, Speed_Occ-1), span()) % integ4_2D(span(0, Speed_Occ-1), span()));
            j++;
        }
        i++;
    }

    for (int i = 1; i < n_Electrons; i++)
    {
        for (int j = i+2; j < n_Electrons; j++)
        {
            Fill_integ4_2D(i,j);
            E1 += accu(tau3(i,j)(span(Speed_Occ, Speed_Occ+Speed_Occ-1), span()) % integ4_2D(span(Speed_Occ, Speed_Occ+Speed_Occ-1), span()));
            j++;
        }
        i++;
    }

    for (int i = 1; i < n_Electrons; i++)
    {
        for (int j = i+1; j < n_Electrons; j++)
        {
            Fill_integ4_2D(i,j);
            E1 += accu(tau3(i,j) % integ4_2D);
            j++;
        }
        i++;
    }

    E1 *= 0.5;
    //E1 += accu(FS_AI % T_1);
    return E1;
}

void CCSD_Memory_optimized::Fill_W1_and_W3()
{
    // Matrix symmetric with W_1(i,j) = W_1(j,i), the terms on the "diagonal" of i==j
    // will be equal to 0 (<-- !)
    // We actually dont even need to store anything in W_1(j,i)
    // since it is accessed symmetricly later on also
    // This halves our storage requirements. Also compress matrix with compact storage.
    // This in total reduce storage by >75%

    for (int i = 0; i < n_Electrons; i++)
    {
        for (int j = i+1; j < n_Electrons; j++)
        {
            // initialize W1
            Fill_integ8_2D(i, j);
            W_1(i,j) = integ8_2D;
        }
    }

    // Optimized W_1, W_3 and parts of F1, F2, F3, T1

    int K, L, E, I, J;
    int KK;

    int index_counter;
    int INDEX_CHECK;

    ////////////////////////////////////////////
    ////////////////////////////////////////////
    ////////////////////////////////////////////
    ////////////////////////////////////////////
    ///                                      ///
    ///     W1 Calculation                   ///
    ///                                      ///
    ///                                      ///
    ////////////////////////////////////////////
    ////////////////////////////////////////////
    ////////////////////////////////////////////
    ////////////////////////////////////////////


    index_counter = 0;
    K = 0;
    for (int k = 0; k < n_Electrons; k++) // k is always even number here, k++ later on also for k+2
    {
        KK = K*Speed_Elec;
        L = 0;
        for (int l = 1; l < k; l++) // here l is always odd number and smaller than k
        {
            INDEX_CHECK = KK+L+Local_Displacement1;

            if (INDEX_CHECK % size == rank)
            {
                Fill_integ6_2D(l,k);
                Fill_integ4_2D(l,k);

                for (int i = 0; i < n_Electrons; i++)
                {
                    MY_OWN_MPI[index_counter] = accu(integ6_2D(span(Speed_Occ, Speed_Occ + Speed_Occ-1), i/2) % T_1(span(Speed_Occ, Speed_Occ+Speed_Occ-1), l/2));
                    i++;
                    index_counter += 1;
                }

                I = 0;
                // Only certain values for i and e are allowed when k and l are even and odd
                for (int i = 0; i < n_Electrons; i++)
                {
                    E = 0;
                    for (int e = 1; e < unocc_orb; e++)
                    {
                        MY_OWN_MPI[index_counter] = -integ6_2D(E+Speed_Occ, I) + accu(integ4_2D(span(0, Speed_Occ-1), E) % T_1(span(0, Speed_Occ-1), I));
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
                        MY_OWN_MPI[index_counter] = -integ6_2D(E, I) + accu(integ4_2D(span(Speed_Occ, Speed_Occ + Speed_Occ-1), E) % T_1(span(Speed_Occ, Speed_Occ+Speed_Occ-1), I));
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
                        MY_OWN_MPI[index_counter] = accu(integ6_2D(span(Speed_Occ, Speed_Occ-1+Speed_Occ), I) % T_1(span(Speed_Occ, Speed_Occ-1+Speed_Occ), J))
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

                        MY_OWN_MPI[index_counter] = accu(integ6_2D(span(0, Speed_Occ-1), I) % T_1(span(0, Speed_Occ-1), J))
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
            INDEX_CHECK = KK+L+Local_Displacement1;

            if (INDEX_CHECK % size == rank)
            {

                Fill_integ6_2D(k,l);
                for (int i = 0; i < n_Electrons; i++)
                {
                    MY_OWN_MPI[index_counter] = accu(integ6_2D(span(Speed_Occ, Speed_Occ + Speed_Occ-1), i/2) % T_1(span(Speed_Occ, Speed_Occ+Speed_Occ-1), l/2));
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
        KK = K*Speed_Elec;
        L = 0;
        for (int l = 0; l < k; l++)
        {

            INDEX_CHECK = KK+L+Local_Displacement1;

            if (INDEX_CHECK % size == rank)
            {

                Fill_integ6_2D(l,k);
                Fill_integ4_2D(l,k);

                I = 0;
                for (int i = 0; i < n_Electrons; i++)
                {
                    E = 0;
                    for (int e = 1; e < unocc_orb; e++)
                    {
                        MY_OWN_MPI[index_counter] = -integ6_2D(E+Speed_Occ, I) + accu(integ4_2D(span(0, Speed_Occ-1), E) % T_1(span(0, Speed_Occ-1), I));
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
                        MY_OWN_MPI[index_counter] = -integ6_2D(E, I) + accu(integ4_2D(span(Speed_Occ, Speed_Occ + Speed_Occ-1), E) % T_1(span(Speed_Occ, Speed_Occ+Speed_Occ-1), I));
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
                        MY_OWN_MPI[index_counter] = accu(integ6_2D(span(Speed_Occ, Speed_Occ-1+Speed_Occ), I) % T_1(span(Speed_Occ, Speed_Occ-1+Speed_Occ), J))
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
                        MY_OWN_MPI[index_counter] = accu(integ6_2D(span(0, Speed_Occ-1), I) % T_1(span(0, Speed_Occ-1), J))
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
        k++;
        K++;
    }

    K = 0;
    for (int k = 0; k < n_Electrons; k++)
    {
        KK = K*Speed_Elec;
        L = 0;
        for (int l = 0; l < k+1; l++)
        {
            INDEX_CHECK = KK+L+Local_Displacement1;

            if (INDEX_CHECK % size == rank)
            {

                Fill_integ6_2D(l,k);
                Fill_integ4_2D(l,k);

                for (int i = 0; i < n_Electrons; i++)
                {
                    MY_OWN_MPI[index_counter] = accu(integ6_2D(span(0, Speed_Occ-1), i/2) % T_1(span(0, Speed_Occ-1), l/2));
                    i++;
                    index_counter += 1;
                }


                I = 0;
                for (int i = 0; i < n_Electrons; i++)
                {
                    E = 0;
                    for (int e = 0; e < unocc_orb; e++)
                    {
                        MY_OWN_MPI[index_counter] = -integ6_2D(E, I) + accu(integ4_2D(span(0, Speed_Occ-1), E) % T_1(span(0, Speed_Occ-1), I));
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
                        MY_OWN_MPI[index_counter] = accu(integ6_2D(span(0, Speed_Occ-1),I) % T_1(span(0, Speed_Occ-1), J))
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
            INDEX_CHECK = KK+L+Local_Displacement1;

            if (INDEX_CHECK % size == rank)
            {
                Fill_integ6_2D(k,l);
                for (int i = 0; i < n_Electrons; i++)
                {
                    MY_OWN_MPI[index_counter] = accu(integ6_2D(span(0, Speed_Occ-1), i/2) % T_1(span(0, Speed_Occ-1), l/2));
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
        KK = K * Speed_Elec;
        L = 0;
        for (int l = 1; l < k+1; l++)
        {
            INDEX_CHECK = KK+L+Local_Displacement1;

            if (INDEX_CHECK % size == rank)
            {
                Fill_integ6_2D(l,k);
                Fill_integ4_2D(l,k);

                I = 0;
                for (int i = 1; i < n_Electrons; i++)
                {
                    E = 0;
                    for (int e = 1; e < unocc_orb; e++)
                    {
                        MY_OWN_MPI[index_counter] = -integ6_2D(E+Speed_Occ, I) + accu(integ4_2D(span(Speed_Occ, Speed_Occ + Speed_Occ-1), E) % T_1(span(Speed_Occ, Speed_Occ+Speed_Occ-1), I));
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
                        MY_OWN_MPI[index_counter] = accu(integ6_2D(span(Speed_Occ, Speed_Occ-1+Speed_Occ), I) % T_1(span(Speed_Occ, Speed_Occ-1+Speed_Occ), J))
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
        k++;
        K++;
    }







    ////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////
    /// Now add W_2 and D_2 and D_3 and D1                  ////////////////////////
    ////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////

    for (int i = 0; i < n_Electrons; i++)
    {
        for (int j = i+1; j < n_Electrons; j++)
        {
            // Initialize W2
            Fill_integ6_2D(i,j);
            W_2(i,j) = integ6_2D;

        }
    }

    int A, M;

    // Optimized W_2 version
    A = 0;
    for (int a = 0; a < unocc_orb; a++) // Only even a
    {
        M = 0;
        for (int m = 0; m < n_Electrons; m++) // only even m, m++ in for loop and m++ later also
        {

            INDEX_CHECK = a/2 * Speed_Elec + m/2;
            if (INDEX_CHECK % size == rank)
            {
                Fill_integ7_2D(a, m);
                Fill_integ5_2D(a, m);
                // Only even values of e allowed
                for (int e = 0; e < unocc_orb; e++)
                {
                    MY_OWN_MPI[index_counter] = accu(integ5_2D(e/2, span()) % T_1(span(0, Speed_Occ-1), m/2).t());
                    index_counter += 1;
                    e++;
                }

                // only even values for i and j allowed.
                // odd i and odd j is skipped,
                // odd i and even j is skipped,
                // even i and odd j is skipped
                // => Most calculations skipped and not needed
                // Also the calculations that are performed use span() to dodge multiplication by 0
                // example 0.5*accu(integ5_2D(span(0, Speed_Occ-1), span())  some lines down

                I = 0;
                for (int i = 0; i < n_Electrons; i++)
                {
                    J = I+1;
                    for (int j = i+2; j < n_Electrons; j++)
                    {
                        MY_OWN_MPI[index_counter] = accu(integ7_2D(span(0, Speed_Occ-1), J) % T_1(span(0, Speed_Occ-1),I))
                                - accu(integ7_2D(span(0, Speed_Occ-1), I) % T_1(span(0, Speed_Occ-1),J))
                                + 0.5*accu(integ5_2D(span(0, Speed_Occ-1), span()) % tau3(i,j)(span(0, Speed_Occ-1), span()));
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

    // odd a
    A = 0;
    for (int a = 1; a < unocc_orb; a++)
    {
        // odd m
        M = 0;
        for (int m = 1; m < n_Electrons; m++)
        {
            INDEX_CHECK = a/2 * Speed_Elec + m/2;
            if (INDEX_CHECK % size == rank)
            {

                Fill_integ7_2D(a, m);
                Fill_integ5_2D(a, m);

                // even i and even j calculated,
                // odd i and odd j calculated
                // however:
                // even i and odd j ignored
                // odd i and even j ignored
                // Half calculations just ignored because it results in 0,
                // Also calculation that is performed is reduced massively

                I = 0;
                for (int i = 0; i < n_Electrons; i++)
                {
                    J = I+1;
                    for (int j = i+2; j < n_Electrons; j++)
                    {
                        MY_OWN_MPI[index_counter] = accu(integ7_2D(span(0, Speed_Occ-1), J) % T_1(span(0, Speed_Occ-1),I))
                                - accu(integ7_2D(span(0, Speed_Occ-1), I) % T_1(span(0, Speed_Occ-1),J))
                                + 0.5*accu(integ5_2D(span(0, Speed_Occ-1),span()) % tau3(i,j)(span(0, Speed_Occ-1), span()));
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
                        MY_OWN_MPI[index_counter] = accu(integ7_2D(span(Speed_Occ, Speed_Occ+Speed_Occ-1), J) % T_1(span(Speed_Occ, Speed_Occ+Speed_Occ-1),I))
                                - accu(integ7_2D(span(Speed_Occ, Speed_Occ+Speed_Occ-1), I) % T_1(span(Speed_Occ, Speed_Occ+Speed_Occ-1),J))
                                + 0.5*accu(integ5_2D(span(Speed_Occ, Speed_Occ+Speed_Occ-1), span()) % tau3(i,j)(span(Speed_Occ, Speed_Occ+Speed_Occ-1), span()));
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
    // odd a
    for (int a = 1; a < unocc_orb; a++)
    {
        M = 0;
        // even m
        for (int m = 0; m < n_Electrons; m++)
        {
            INDEX_CHECK = a/2 * Speed_Elec  + m/2;
            if (INDEX_CHECK % size == rank)
            {

                Fill_integ7_2D(a, m);
                Fill_integ5_2D(a, m);

                // only certain values allowed for i and j. most calculations ignored because not needed
                I = 0;
                for (int i = 0; i < n_Electrons; i++)
                {
                    J = I;
                    for (int j = i+1; j < n_Electrons; j++)
                    {
                        MY_OWN_MPI[index_counter] = accu(integ7_2D(span(0, Speed_Occ-1), J) % T_1(span(0, Speed_Occ-1),I))
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
                        MY_OWN_MPI[index_counter] = accu(integ7_2D(span(Speed_Occ, Speed_Occ+Speed_Occ-1), J) % T_1(span(Speed_Occ, Speed_Occ+Speed_Occ-1),I))
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

            INDEX_CHECK = a/2 * Speed_Elec + m/2;
            if (INDEX_CHECK % size == rank)
            {

                Fill_integ5_2D(a,m);
                Fill_integ7_2D(a,m);

                for (int e = 0; e < unocc_orb; e++)
                {
                    MY_OWN_MPI[index_counter] = accu(integ5_2D(e/2, span()) % T_1(span(Speed_Occ, unocc_orb-1), m/2).t());
                    index_counter += 1;
                    e++;
                }

                I = 0;
                for (int i = 0; i < n_Electrons; i++)
                {
                    J = I;
                    for (int j = i+1; j < n_Electrons; j++)
                    {
                        MY_OWN_MPI[index_counter] = accu(integ7_2D(span(0, Speed_Occ-1), J) % T_1(span(0, Speed_Occ-1),I))
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
                        MY_OWN_MPI[index_counter] = accu(integ7_2D(span(Speed_Occ, Speed_Occ+Speed_Occ-1), J) % T_1(span(Speed_Occ, Speed_Occ+Speed_Occ-1),I))
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

    // Seperation noted here, F2 and F3 parts start here

    for (int e = 0; e < unocc_orb; e++)
    {
        for (int m = 0; m < n_Electrons; m++)
        {
            INDEX_CHECK = e/2* Speed_Elec + m/2;
            if (INDEX_CHECK % size == rank)
            {
                Fill_integ2_2D_even_even(e,m);
                for (int i = 0; i < n_Electrons; i++)
                {
                    MY_OWN_MPI[index_counter] = 0.5*accu(integ2_2D % t2.at(e,i));
                    index_counter += 1;
                    i++;
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
            INDEX_CHECK = e/2* Speed_Elec + m/2;
            if (INDEX_CHECK % size == rank)
            {
                Fill_integ2_2D_odd_even(e,m);
                for (int i = 0; i < n_Electrons; i++)
                {
                    MY_OWN_MPI[index_counter] = 0.5*accu(integ2_2D(span(0, Speed_Occ-1), span()) % t2.at(e,i)(span(0, Speed_Occ-1), span()));
                    index_counter += 1;
                    i++;
                }
            }
            m++;
        }
        e++;
    }

    // Seperation noted here

    for (int e = 0; e < unocc_orb; e++)
    {
        for (int m = 0; m < n_Electrons; m++)
        {
            INDEX_CHECK = e/2* Speed_Elec+ m/2;
            if (INDEX_CHECK % size == rank)
            {
                Fill_integ2_2D_even_even(e,m);

                // F1 term
                MY_OWN_MPI[index_counter] = accu(integ2_2D % T_1);
                index_counter += 1;

                for (int a = 0; a < unocc_orb; a++)
                {
                    MY_OWN_MPI[index_counter] = 0.5*accu(integ2_2D % t2.at(a,m));
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
            INDEX_CHECK = e/2 * Speed_Elec+ m/2;
            if (INDEX_CHECK % size == rank)
            {
                Fill_integ2_2D_even_odd(e,m);
                for (int a = 0; a < unocc_orb; a++)
                {
                    MY_OWN_MPI[index_counter] = 0.5*accu(integ2_2D(span(Speed_Occ, unocc_orb-1), span()) % t2.at(a,m)(span(Speed_Occ, unocc_orb-1), span()));
                    index_counter += 1;
                    a++;
                }
            }
            m++;
        }
        e++;
    }

    // part of T1_new here
    // This part of T1 amplitudes calculated here because parts of MOs are distributed
    // in a way that fits in good here.

    for (int a = 0; a < unocc_orb; a++)
    {
        for (int k = 0; k < n_Electrons; k++)
        {
            INDEX_CHECK = a/2 * Speed_Elec+ k/2;
            if (INDEX_CHECK % size == rank)
            {
                Fill_integ5_2D(a, k);
                for (int i = 0; i < k; i++)
                {
                    MY_OWN_MPI[index_counter] = 0.5*accu(integ5_2D % tau3(i,k));
                    index_counter += 1;
                    i++;
                }

                // i will be less than k when k = odd number
                // We want i to be an even number since a is even number, hence we start at i = k+1
                for (int i = (k+1+(k+1)%2); i < n_Electrons; i++)
                {
                    MY_OWN_MPI[index_counter] = -0.5*accu(integ5_2D % tau3(k,i));
                    index_counter += 1;
                    i++;
                }
            }
        }
        a++;
    }
}

void CCSD_Memory_optimized::Fill_W2()
{
    // Matrix symmetric with W_2(i,j) = -W_2(j,i) and always 0 on the diagonal where i == j (<-- !)
    // When i = j the amplitude itself is illigal, so no need even if it is not zero.. which it is
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

    // This function is a mapping of W4 into armadillo fields, W4 is calculated elsewere

    int index_counter;
    //int AA,I,INDEX_CHECK;

    // Scatter the information, all to all scatter through MPI_Alltoallw call
    MPI_Alltoallw(MY_OWN_MPI, Global_Worksize_2[rank], Global_Displacement_2[rank],             mpi_types_array, // This last one is filled with MPI_Doubles
                  SHARED_INFO_MPI, Global_Worksize_2_1[rank], Global_Displacement_2_1[rank],    mpi_types_array,
                  MPI_COMM_WORLD);

    index_counter = 0;

    for (int J = 0; J < size ; J++)
    {
            // Arrange this so one node gets its information in sequence
            for (int a = 0; a < unocc_orb; a++)
            {
                for(int i = Where_To_Start_Part2(rank,a); i < n_Electrons; i+=jump)
                {
                        for (int m = Where_To_Start_Part2(J,a); m < n_Electrons; m+=jump)
                        {
                          //  INDEX_CHECK = AA + m/2;
                          //  if (INDEX_CHECK % size == J)
                          //  {
                                for (int e = 0; e < unocc_orb; e++)
                                {

                                    W_4(a,i)(e/2,m/2) = SHARED_INFO_MPI[index_counter];
                                    index_counter += 1;
                                    e++;
                                }
                           // }
                            //m++;
                        }

                        for (int m = Where_To_Start_Part2(J,a)+1; m < n_Electrons; m+=jump)
                        {
                          //  INDEX_CHECK = AA + m/2;
                          //  if (INDEX_CHECK % size == J)
                          //  {
                                for (int e = 1; e < unocc_orb; e++)
                                {
                                    W_4(a,i)(e/2+Speed_Occ,m/2) = SHARED_INFO_MPI[index_counter];
                                    index_counter += 1;
                                    e++;
                                }
                           // }
                           // m++;
                        }

                }
                a++;
            }

            for (int a = 0; a < unocc_orb; a++)
            {
                //AA = a/2*Speed_Elec;
                for (int i = Where_To_Start_Part2(rank,a)+1; i < n_Electrons; i+=jump)
                {
                        for (int m = Where_To_Start_Part2(J,a)+1; m < n_Electrons; m+=jump)
                        {
                            //INDEX_CHECK = AA + m/2;
                           // if (INDEX_CHECK % size == J)
                           // {
                                for (int e = 0; e < unocc_orb; e++)
                                {
                                    W_4(a,i)(e/2,m/2) = SHARED_INFO_MPI[index_counter];
                                    index_counter += 1;
                                    e++;
                                }
                            //}
                           // m++;
                        }

                }
                a++;
            }

            for (int a = 1; a < unocc_orb; a++)
            {
                for (int i = Where_To_Start_Part2(rank,a-1); i < n_Electrons; i+=jump)
                {
                        for (int m = Where_To_Start_Part2(J,a-1); m < n_Electrons; m+=jump)
                        {
                           // INDEX_CHECK = AA + m/2;
                          //  if (INDEX_CHECK % size == J)
                           // {
                                for (int e = 1; e < unocc_orb; e++)
                                {
                                    W_4(a,i)(e/2+Speed_Occ,m/2) = SHARED_INFO_MPI[index_counter];
                                    index_counter += 1;
                                    e++;
                                }
                           // }
                            //m++;
                        }

                }
                a++;
            }

            for (int a = 1; a < unocc_orb; a++)
            {
                for (int i = Where_To_Start_Part2(rank,a-1)+1; i < n_Electrons; i+=jump)
                {
                        for (int m = Where_To_Start_Part2(J,a-1); m < n_Electrons; m+=jump)
                        {
                            //INDEX_CHECK = AA + m/2;
                            //if (INDEX_CHECK % size == J)
                            //{
                                for (int e = 0; e < unocc_orb; e++)
                                {
                                    W_4(a,i)(e/2,m/2) = SHARED_INFO_MPI[index_counter];
                                    index_counter += 1;
                                    e++;
                                }
                            //}
                            //m++;
                        }

                        for (int m = Where_To_Start_Part2(J,a-1)+1; m < n_Electrons; m+=jump)
                        {
                           // INDEX_CHECK = AA + m/2;
                           // if (INDEX_CHECK % size == J)
                           // {
                                for (int e = 1; e < unocc_orb; e++)
                                {
                                    W_4(a,i)(e/2+Speed_Occ,m/2) = SHARED_INFO_MPI[index_counter];
                                    index_counter += 1;
                                    e++;
                                }
                           // }
                            //m++;
                        }

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

    // Run F1 in parallel in the W_1_and_W_3 function
}

void CCSD_Memory_optimized::Fill_F2()
{
    int I,M;

    // spin up up or down down wont matter
    for (int e = 1; e < unocc_orb; e++)
    {
        for (int m = 1; m < n_Electrons; m++)
        {
            D1(e/2+Speed_Occ,m/2) = D1(e/2, m/2);
        }
    }

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

    // spin up up or down down wont matter, upper and lower half of matrix will be the same
    // because of compact storage
    for (int k = 0; k < n_Electrons; k++)
    {
        for (int i = 0; i < n_Electrons; i++)
        {
            D2(k/2+Speed_Elec,i/2) = D2(k/2,i/2);
            i++;
        }
        k++;
    }
}

void CCSD_Memory_optimized::Fill_F3()
{
    // F3 is not stored as compact storage normally is, since it is only accesed row wise later
    // This means we dont need to simplify that column access part, so just store it normal without zeroes


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

    // spin up up or down down wont matter, compact storage turns this into
    // upper matrix and lower matrix which are identical
    for (int a = 1; a < unocc_orb; a++)
    {
        E = 0;
        for (int e = 1; e < unocc_orb; e++)
        {
            D3.at(a,E) = D3(a-1,E);// - accu(D1.row(e/2+Speed_Occ) % T_1.row(a/2+Speed_Occ));
            e++;
            E++;
        }
        a++;
    }
}

void CCSD_Memory_optimized::Fill_tau()
{
    // Calculate tau , not done in parallel
    // Parallel here may be more effective, but it is good for now.
    // n^4 scaling, which is small since largest calculations are n^6
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
    // 2D mapping of tau so we dont need double storage,
    // This function enables us to map into very effective external math library use.
    // mapping is n^4 operation

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
    // This is moved to Fill_t2_new() in parallel implementation
    // Parts go to Fill_W1_and_W3 due to memory distribution
}

void CCSD_Memory_optimized::Fill_t2_new()
{
    // Benchmark H2O STO-3G basis set: -0.0501273 au, 20 iteration

    // Optimized Version of T2 calculations
    // Runs in parallel
    // THIS IS THE ONLY PLACE MO9 IS ACCESSED - MEMORY SPREAD ON NODES IN a-b GRID
    // T2-new Stored in an array of size number of calcs for node
    // This is for communication optimization
    // Dont store more than one symetric term for communication minimization
    // Use all symmetries and skip meny calculations where one term will be 0 but the other not or both 0 etc
    // Try not to skip calculations on terms that are both not 0 :-D

    // Function is now prepeared for reading MOs from disk in an extremely optimized way

    int index_counter;
    int A,B, AA;
    int sum_a_n;
    double temp; // Used for optimization in T1 calculations
    int INDEX_CHECK;
    //int work;




        index_counter = 0;
        sum_a_n = 0;
        for (int a = 0; a < unocc_orb; a++)
        {
            sum_a_n += a/2;
            A = a/2;
            AA = A* Speed_Occ - sum_a_n;

            // If we want to read from file and distribute other arrays:
            // Read in from file here a 1 dimensional array that contains single bar integrals for a fixed index "a".
            // We can use the same ones for even a and odd a and for all b,i,j, to ensure we keep our speed advantage

            // MPI_File_read(...)

            // At this point calculations start and
            // a is an even number

            for (int b = a+2; b < unocc_orb; b++) // b even
            {
                B = b/2;
                INDEX_CHECK = AA+B;

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
                            MY_OWN_MPI[index_counter] = (-MOLeftovers(a/2, b/2)(j/2, i/2) + MOLeftovers(a/2, b/2)(i/2, j/2) // integ2.at(b,j)(a/2,i/2)

                            // P(ab) P(ij) [W_4] t_jk^bc, ARMADILLO
/*
                                    + accu(t2.at(b,j) % W_4.at(a,i))
                                    - accu(t2.at(b,i) % W_4.at(a,j))
                                    - accu(t2.at(a,j) % W_4.at(b,i))
                                    + accu(t2.at(a,i) % W_4.at(b,j))
*/

                                                         + W_5(a,b)(i/2,j/2)
                                                         - W_5(a,b)(j/2,i/2)
                                                      //   - W_5(b,a)(i/2,j/2)
                                                      //   + W_5(b,a)(j/2,i/2)

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
                                    + 0.5*accu(integ3_2D(span(0, Speed_Occ-1), span()) % tau3(i,j)(span(0, Speed_Occ-1), span()))) // Half matrix = 0, skip this

                            // Denominator
                                     / (DEN_AI(a/2,i/2)+DEN_AI(b/2, j/2));


                            index_counter++;
                            j++;
                        }
                        i++;
                    }
                }
                b++;
            }

            for (int b = a+1; b < unocc_orb; b++) // b odd
            {
                B = b/2;
                INDEX_CHECK = AA+B;

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
                            MY_OWN_MPI[index_counter] = (MOLeftovers(a/2, b/2)(i/2, j/2) // integ2.at(b,j)(a/2,i/2)

                            // P(ab) P(ij) [W_4] t_jk^bc, ARMADILLO
/*
                                    - accu(t2.at(a,j)(span(Speed_Occ, Speed_Occ+Speed_Occ-1), span()) % W_4.at(b,i)(span(Speed_Occ, Speed_Occ+Speed_Occ-1), span()))
                                    - accu(t2.at(b,i)(span(0, Speed_Occ-1), span()) % W_4.at(a,j)(span(0, Speed_Occ-1), span()))
                                    + accu(t2.at(a,i) % W_4.at(b,j))
                                    + accu(t2.at(b,j) % W_4.at(a,i))
*/

                                                         + W_5(a,b)(i/2,j/2)
                                                       //  - W_5(b,a)(i/2,j/2)
                                                         - W_5(a,b)(j/2+Speed_Elec, i/2)
                                                       //  + W_5(b,a)(j/2+Speed_Elec, i/2)

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
                                    + 0.5*accu(integ3_2D % tau3(i,j)))

                            // Denominator
                                     / (DEN_AI(a/2,i/2)+DEN_AI(b/2+Speed_Occ, j/2));

                            index_counter++;
                            j++;
                        }
                        i++;
                    }

                    for (int i = 1; i < n_Electrons; i++) // i odd
                    {
                        for (int j = i+1; j < n_Electrons; j++) // j even
                        {
                            // I_ab^ij
                            MY_OWN_MPI[index_counter] = (-MOLeftovers(a/2, b/2)(j/2, i/2) // integ2.at(b,j)(a/2,i/2)

                            // P(ab) P(ij) [W_4] t_jk^bc, ARMADILLO
/*
                                    - accu(t2.at(a,j) % W_4.at(b,i))
                                    - accu(t2.at(b,i) % W_4.at(a,j))
                                    + accu(t2.at(a,i)(span(Speed_Occ, Speed_Occ+Speed_Occ-1), span()) % W_4.at(b,j)(span(Speed_Occ, Speed_Occ+Speed_Occ-1), span()))
                                    + accu(t2.at(b,j)(span(0, Speed_Occ-1), span()) % W_4.at(a,i)(span(0, Speed_Occ-1), span()))
*/

                                                         + W_5(a,b)(i/2+Speed_Elec,j/2)
                                                      //   - W_5(b,a)(i/2+Speed_Elec,j/2)
                                                         - W_5(a,b)(j/2, i/2)
                                                      //   + W_5(b,a)(j/2, i/2)

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
                                    + 0.5*accu(integ3_2D % tau3(i,j)))

                            // Denominator
                                    / (DEN_AI(a/2,i/2)+DEN_AI(b/2+Speed_Occ, j/2));

                            index_counter++;
                            j++;
                        }
                        i++;
                    }
                }
                b++;
            }



            a++;
            // At this point a is an odd number, but a/2 is the same value as before so we can use the same single bar integrals potentially

            for (int b = a+2; b < unocc_orb; b++) // b odd
            {
                B = b/2;
                INDEX_CHECK = AA+B;

                if (INDEX_CHECK % size == rank)
                {
                    Fill_integ3_2D(a, b);
                    Fill_integ9_2D(a, b);
                    Fill_2D_tau(a, b);

                    for (int i = 1; i < n_Electrons; i++) // i odd
                    {
                        for (int j = i+2; j < n_Electrons; j++) // j odd
                        {
                            // I_ab^ij
                            MY_OWN_MPI[index_counter] = (MOLeftovers(a/2, b/2)(i/2, j/2) - MOLeftovers(a/2, b/2)(j/2,i/2) // integ2.at(b,j)(a/2+Speed_Occ,i/2)

                            // P(ab) P(ij) [W_4] t_jk^bc, ARMADILLO
/*
                                    + accu(t2.at(b,j) % W_4.at(a,i))
                                    - accu(t2.at(b,i) % W_4.at(a,j))
                                    - accu(t2.at(a,j) % W_4.at(b,i))
                                    + accu(t2.at(a,i) % W_4.at(b,j))
*/

                                                         + W_5(a,b)(i/2+Speed_Elec, j/2)
                                                         - W_5(a,b)(j/2+Speed_Elec, i/2)
                                                       //  - W_5(b,a)(i/2+Speed_Elec, j/2)
                                                       //  + W_5(b,a)(j/2+Speed_Elec, i/2)

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
                                    + 0.5*accu(integ3_2D(span(Speed_Occ, Speed_Occ+Speed_Occ-1), span()) % tau3(i,j)(span(Speed_Occ, Speed_Occ+Speed_Occ-1), span())))

                            // Denominator
                                     / (DEN_AI(a/2+Speed_Occ,i/2)+DEN_AI(b/2+Speed_Occ, j/2));

                            index_counter++;
                            j++;
                        }
                        i++;
                    }

                }
                b++;
            }

            for (int b = a+1; b < unocc_orb; b++) // b even
            {
                B = b/2;
                INDEX_CHECK = AA+B;

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
                            MY_OWN_MPI[index_counter] = (-MOLeftovers(a/2, b/2)(j/2, i/2) // integ2.at(b,j)(a/2+Speed_Occ,i/2)

                            // P(ab) P(ij) [W_4] t_jk^bc, ARMADILLO
/*
                                    - accu(t2.at(a,j) % W_4.at(b,i))
                                    - accu(t2.at(b,i) % W_4.at(a,j))
                                    + accu(t2.at(a,i)(span(0, Speed_Occ-1), span()) % W_4.at(b,j)(span(0, Speed_Occ-1), span()))
                                    + accu(t2.at(b,j)(span(Speed_Occ, Speed_Occ+Speed_Occ-1), span()) % W_4.at(a,i)(span(Speed_Occ, Speed_Occ+Speed_Occ-1), span()))
*/

                                                         + W_5(a,b)(i/2,j/2)
                                                        // - W_5(b,a)(i/2,j/2)
                                                         - W_5(a,b)(j/2+Speed_Elec,i/2)
                                                        // + W_5(b,a)(j/2+Speed_Elec,i/2)

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
                                    + 0.5*accu(integ3_2D % tau3(i,j)))

                            // Denominator
                                     / (DEN_AI(a/2+Speed_Occ,i/2)+DEN_AI(b/2, j/2));

                            index_counter++;
                            j++;
                        }
                        i++;
                    }

                    for (int i = 1; i < n_Electrons; i++) // i odd
                    {
                        for (int j = i+1; j < n_Electrons; j++) // j even
                        {
                            // I_ab^ij
                            MY_OWN_MPI[index_counter] = (MOLeftovers(a/2, b/2)(i/2, j/2) // integ2.at(b,j)(a/2+Speed_Occ,i/2)

                            // P(ab) P(ij) [W_4] t_jk^bc, ARMADILLO
/*
                                    - accu(t2.at(a,j)(span(0, Speed_Occ-1), span()) % W_4.at(b,i)(span(0, Speed_Occ-1), span()))
                                    - accu(t2.at(b,i)(span(Speed_Occ, Speed_Occ+Speed_Occ-1), span()) % W_4.at(a,j)(span(Speed_Occ, Speed_Occ+Speed_Occ-1), span()))
                                    + accu(t2.at(a,i) % W_4.at(b,j))
                                    + accu(t2.at(b,j) % W_4.at(a,i))

                                                     */
                                                         + W_5(a,b)(i/2+Speed_Elec,j/2)
                                                         - W_5(a,b)(j/2, i/2)
                                                         //- W_5(b,a)(i/2+Speed_Elec, j/2)
                                                         //+ W_5(b,a)(j/2, i/2)


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
                                    + 0.5*accu(integ3_2D % tau3(i,j)))

                            // Denominator
                                     / (DEN_AI(a/2+Speed_Occ,i/2)+DEN_AI(b/2, j/2));

                            index_counter++;
                            j++;
                        }
                        i++;
                    }
                }
                b++;
            }
        }



        // T1 part
        // We only need  to calculate either up up or down down T1 amplitudes,
        // Since they will be the same due to spin restriction
        for (int a = 0; a < unocc_orb; a++)
        {
            for (int i = 0; i < n_Electrons; i++)
            {
                INDEX_CHECK = a/2 * Speed_Elec+ i/2;
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

                    MY_OWN_MPI[index_counter] = (0.5*temp

                            // t_ik^ac [F_1]_c^k
                            + accu(D1 % t2(a,i))

                            // - t_k^a [F_2]_i^k
                            - accu(D2(span(i%2*Speed_Elec, i%2*Speed_Elec+Speed_Elec-1), i/2) % T_1.row(a/2+a%2*Speed_Occ).t())

                            // f_ac t_i^c
                            //+ accu(FS_AB(a/2, span()) % T_1(span(0, Speed_Occ-1), i/2).t())

                            // I_ka^ci t_k^c
                            - accu(integ10_2D % T_1));
                    index_counter += 1;
                }

                i++;
            }
            a++;
        }
}

void CCSD_Memory_optimized::Map_T_new()
{
    // Mapping for parallel implementation!
    // Will resemble the CALCULATE T2 FUNCTION EXACTLY!!!!
    //                       and we can then use index counter
    //                      this is to make the transfer extreme easy and efficient

    double temp;
    int A,B,I,J, INDEX_CHECK;
    int index_counter;
    int work;
    int sum_a_n;
    int AA;


    // All-to-all communication for effective MPI use.
    MPI_Allgatherv(MY_OWN_MPI, WORK_EACH_NODE(rank), MPI_DOUBLE,
                   SHARED_INFO_MPI, Work_Each_Node_T2_Parallel, Displacement_Each_Node_T2_Parallel,
                   MPI_DOUBLE, MPI_COMM_WORLD);


    index_counter = 0;

    for (int K = 0; K < size; K++)
    {   // This loop is here to pull out the correct numbers in the correct sequence,
        // Since we just gather all numbers in a quite random way and then used allgatherv
        // However since index counter counts through the values calculated by each CPU,
        // we can pull them out just were they should go
        // We copy pasted the function above and changed instead of calculating we map into t2 amplitudes

        sum_a_n = 0;

        for (int a = 0; a < unocc_orb; a++)
        {
            sum_a_n += a/2;
            A = a/2;
            AA = A * Speed_Occ - sum_a_n;


            // At this point a is an even number

            for (int b = a+2; b < unocc_orb; b++) // b even
            {
                B = b/2;
                INDEX_CHECK = AA+B;

                if (INDEX_CHECK % size == K)
                {
                    for (int i = 0; i < n_Electrons; i++) // i even
                    {
                        I = i/2;
                        for (int j = i+2; j < n_Electrons; j++) // j even
                        {
                            J = j/2;
                            temp = SHARED_INFO_MPI[index_counter];// / (DEN_AI(a/2+a%2*Speed_Occ,i/2)+DEN_AI(b/2+b%2*Speed_Occ, j/2));

                            t2(a,i)(B, J) = temp;
                            t2(b,i)(A, J) = -temp;
                            t2(a,j)(B, I) = -temp;
                            t2(b,j)(A, I) = temp;

                            index_counter++;
                            j++;
                        }
                        i++;
                    }
                }
                b++;
            }

            for (int b = a+1; b < unocc_orb; b++) // b odd
            {
                B = b/2;
                INDEX_CHECK = AA+B;

                if (INDEX_CHECK % size == K)
                {
                    for (int i = 0; i < n_Electrons; i++) // i even
                    {
                        I = i/2;
                        for (int j = i+1; j < n_Electrons; j++) // j odd
                        {
                            J = j/2;
                            temp = SHARED_INFO_MPI[index_counter];// / (DEN_AI(a/2+a%2*Speed_Occ,i/2)+DEN_AI(b/2+b%2*Speed_Occ, j/2));

                            t2(a,i)(B+Speed_Occ, J) = temp;
                            t2(b,i)(A, J) = -temp;
                            t2(a,j)(B+Speed_Occ, I) = -temp;
                            t2(b,j)(A, I) = temp;

                            index_counter++;
                            j++;
                        }
                        i++;
                    }

                    for (int i = 1; i < n_Electrons; i++) // i odd
                    {
                        I = i/2;
                        for (int j = i+1; j < n_Electrons; j++) // j even
                        {
                            J = j/2;

                            temp = SHARED_INFO_MPI[index_counter];// / (DEN_AI(a/2+a%2*Speed_Occ,i/2)+DEN_AI(b/2+b%2*Speed_Occ, j/2));

                            t2(a,i)(B+Speed_Occ, J) = temp;
                            t2(b,i)(A, J) = -temp;
                            t2(a,j)(B+Speed_Occ, I) = -temp;
                            t2(b,j)(A, I) = temp;

                            index_counter++;
                            j++;
                        }
                        i++;
                    }
                }
                b++;
            }



            a++;
            // At this point a is an odd number, but a/2 is the same value as before
            for (int b = a+2; b < unocc_orb; b++) // b odd
            {
                B = b/2;
                INDEX_CHECK = AA+B;

                if (INDEX_CHECK % size == K)
                {
                    for (int i = 1; i < n_Electrons; i++) // i odd
                    {
                        I = i/2;
                        for (int j = i+2; j < n_Electrons; j++) // j odd
                        {
                            J = j/2;

                            temp = SHARED_INFO_MPI[index_counter];// / (DEN_AI(a/2+a%2*Speed_Occ,i/2)+DEN_AI(b/2+b%2*Speed_Occ, j/2));

                            t2(a,i)(B+Speed_Occ, J) = temp;
                            t2(b,i)(A+Speed_Occ, J) = -temp;
                            t2(a,j)(B+Speed_Occ, I) = -temp;
                            t2(b,j)(A+Speed_Occ, I) = temp;

                            index_counter++;
                            j++;
                        }
                        i++;
                    }

                }
                b++;
            }

            for (int b = a+1; b < unocc_orb; b++) // b even
            {
                B = b/2;
                INDEX_CHECK = AA+B;

                if (INDEX_CHECK % size == K)
                {
                    for (int i = 0; i < n_Electrons; i++) // i even
                    {
                        I = i/2;
                        for (int j = i+1; j < n_Electrons; j++) // j odd
                        {
                            J = j/2;
                            temp = SHARED_INFO_MPI[index_counter];// / (DEN_AI(a/2+a%2*Speed_Occ,i/2)+DEN_AI(b/2+b%2*Speed_Occ, j/2));

                            t2(a,i)(B, J) = temp;
                            t2(b,i)(A+Speed_Occ, J) = -temp;
                            t2(a,j)(B, I) = -temp;
                            t2(b,j)(A+Speed_Occ, I) = temp;

                            index_counter++;
                            j++;
                        }
                        i++;
                    }

                    for (int i = 1; i < n_Electrons; i++) // i odd
                    {
                        I = i/2;
                        for (int j = i+1; j < n_Electrons; j++) // j even
                        {
                            J = j/2;

                            temp = SHARED_INFO_MPI[index_counter];// / (DEN_AI(a/2+a%2*Speed_Occ,i/2)+DEN_AI(b/2+b%2*Speed_Occ, j/2));

                            t2(a,i)(B, J) = temp; // -t2(a,i-1)(b/2,j/2); // temp;
                            t2(b,i)(A+Speed_Occ, J) = -temp; // -t2(b,i-1)(A+Speed_Occ, j/2); // -temp;
                            t2(a,j)(B, I) = -temp; // -t2(a,j-1)(b/2, i/2); //-temp;
                            t2(b,j)(A+Speed_Occ, I) = temp; // -t2(b,j-1)(A+Speed_Occ, i/2); //temp;

                            index_counter++;
                            j++;
                        }
                        i++;
                    }
                }
                b++;
            }
        }




        // T1 part

        for (int a = 0; a < unocc_orb; a++)
        {
            A = a/2;
            for (int i = 0; i < n_Electrons; i++)
            {
                INDEX_CHECK = A *Speed_Elec+i/2;

                if (INDEX_CHECK % size == K)
                {
                    T_1_new(A, i/2) += SHARED_INFO_MPI[index_counter];
                    index_counter++;
                }
                i++;
            }
            a++;
        }
    }


    for (int a = 1; a < unocc_orb; a++)
    {
        A = a/2;
        for (int i = 1; i < n_Electrons; i++)
        {
            T_1_new(A+Speed_Occ, i/2) = T_1_new(A,i/2);
        }
    }
}

void CCSD_Memory_optimized::Prepear_AOs(int nr_freeze)
{

    // Freeze Core = Ignorer det første antallet frozen core, ikke lagre de bare ....
    //


    // Effective AOtoMO transformation

    int matsize = Matrix_Size/2;
    int matsize_frozen = matsize - nr_freeze;
    int Speed_Elec_Frozen = Speed_Elec - nr_freeze;

    int A,B, INDEX_CHECK, INDEX_CHECK3;
    mat ttt;

    // Allocate arrays needed to store single bar four index integrals, or MOs.

    MO8.set_size(Speed_Elec_Frozen, Speed_Elec_Frozen);
    ttt = zeros(Speed_Elec_Frozen, Speed_Elec_Frozen);
    for (int a = 0; a < Speed_Elec_Frozen; a++)
    {
        for (int b = a; b < Speed_Elec_Frozen; b++)
        {
            MO8(a,b) = ttt;
        }
    }

    MO6_1.set_size(Speed_Elec_Frozen, Speed_Elec_Frozen);
    MO6_2.set_size(Speed_Elec_Frozen, Speed_Elec_Frozen);

    for (int a = 0; a < Speed_Elec_Frozen; a++)
    {
        for (int b = a; b < Speed_Elec_Frozen; b++)
        {
            MO6_1(a,b) = zeros(Speed_Occ, Speed_Elec_Frozen);
            MO6_2(a,b) = zeros(Speed_Elec_Frozen, Speed_Occ);
        }
    }

    MO4.set_size(Speed_Elec_Frozen, Speed_Elec_Frozen);
    ttt = zeros(Speed_Occ, Speed_Occ);
    for (int a = 0; a < Speed_Elec_Frozen; a++)
    {
        for (int b = a; b < Speed_Elec_Frozen; b++)
        {
            MO4(a,b) = ttt;
        }
    }

    int index_A;
    MO9.set_size(Speed_Occ, Speed_Occ);
    ttt = zeros(matsize_frozen, matsize_frozen);
    for (int a = 0; a < Speed_Occ; a++)
    {
        for (int b = a; b < Speed_Occ; b++)
        {
            index_A = Calc_sum_a_n(a);
            INDEX_CHECK = a*Speed_Occ - index_A + b;
            if (INDEX_CHECK % size == rank)
            {
                MO9(a,b) = ttt;
            }
        }
    }

    MO2.set_size(Speed_Occ, Speed_Elec_Frozen);
    MO10.set_size(Speed_Occ, Speed_Elec_Frozen);
    MO3.set_size(Speed_Occ, Speed_Elec_Frozen);
    ttt = zeros(matsize_frozen, matsize_frozen);
    for (int a = 0; a < Speed_Occ; a++)
    {
        for (int b = 0; b < Speed_Elec_Frozen; b++)
        {
            // Distribute in memory
            INDEX_CHECK = a+ b;
            if ((a*Speed_Elec_Frozen+b)%size == rank)
            {
                MO10(a,b) = ttt;
                MO2(a,b) = zeros(Speed_Occ, Speed_Elec_Frozen);
                MO3(a,b) = zeros(Speed_Elec_Frozen, Speed_Occ);
            }
        }
    }

    ttt = zeros(Speed_Elec_Frozen, Speed_Elec_Frozen);
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

    // Map out size of arrays for MPI use
    mat work_to_do = MO_Grid_A_B(nr_freeze);
    mat work_do_2 = MO_Grid_I_J(nr_freeze);
    mat where_to_start2 = zeros(matsize, size);
    mat where_to_start = zeros(matsize, size);
    int dist;
    int G;

    // Figure out were each CPU starts calculations
    for (int K = 0; K < size; K++)
    {
        for (int i = 0; i < matsize; i++)
        {
            G = 0;
            dist = 0;
            for (int NN = 0; NN < i; NN++)
            {
                dist += NN;
            }

            while ((i+G+dist)%size != K)
            {
                G++;
            }
            where_to_start2(i,K) = G;
        }

        for (int i = 0; i < matsize; i++)
        {
            G = 0;
            while ((i+G)%size != K)
            {
                G++;
            }
            where_to_start(i,K) = G;
        }
    }

    int cube_send_size = max(max(work_do_2));
    int nr2 = max(max(work_to_do));
    if (nr2 > cube_send_size)
    {
        cube_send_size = nr2;
    }

    if (cube_send_size == 0)
    {
        cube_send_size = 1; // just in case
    }
    int ccube_send_size = matsize*matsize*matsize;

    vector<mat> compact_mo5;
    mat temp_mat = zeros(matsize, matsize);

    vec local_dist1 = zeros(matsize), local_dist2 = zeros(matsize);

    for (int a = nr_freeze; a < matsize; a++)
    {
        A = a - Speed_Elec;

        if (a >= Speed_Elec)
        {
            // Figure out local distribution

            for (int i = a; i < matsize*size; i++) // (0 -> Speed_Occ) is equal to (a -> matsize) - Speed_Elec
            {
                B = (i - Speed_Elec);
                index_A = Calc_sum_a_n(A);
                INDEX_CHECK = A*Speed_Occ - index_A + B;

                if (INDEX_CHECK % size == rank)
                {
                    local_dist1(a) = i;
                    i = matsize*size+1;
                }
            }

            for (int j = nr_freeze; j < matsize*size; j++) // <-- j Where_To_Start_Part2(rank,a
            {
                if ((A*Speed_Elec_Frozen+j-nr_freeze)%size == rank)
                {
                    local_dist2(a) = j;
                    j = matsize*size+1;
                }
            }
        }
    }

    // Use these as one dimensional intermediates
    vec compact_mo10 = zeros(matsize);
    vec compact_mo11 = zeros(matsize);

    // Use this as n^3 sized intermediate to store two quarter transformed values
    for (int i = 0; i < matsize; i++)
    {
        compact_mo5.push_back(temp_mat);
    }

    // Allocate arrays needed for MPI all-to-all function use
    // Displacement and nubmer of bytes for communication
    int** recieve_array = (int**)malloc(matsize*sizeof(int*));
    int** displacement_array = (int**)malloc(matsize*sizeof(int*));
    int** recieve_array2 = (int**)malloc(matsize*sizeof(int*));
    int** displacement_array2 = (int**)malloc(matsize*sizeof(int*));

    for (int i = 0; i < matsize; i++)
    {
        recieve_array2[i] = (int*)malloc(size*sizeof(int));
        displacement_array2[i] = (int*)malloc(size*sizeof(int));
        recieve_array[i] = (int*)malloc(size*sizeof(int));
        displacement_array[i] = (int*)malloc(size*sizeof(int));
    }

    for (int a = 0; a < matsize; a++)
    {
        for (int i = 0; i < size; i++)
        {
            recieve_array[a][i] = work_do_2(i,a);
        }

        displacement_array[a][0] = 0;
        for (int i = 1; i < size; i++)
        {
            displacement_array[a][i] = displacement_array[a][i-1] + recieve_array[a][i-1];
        }
    }

    for (int i = 0; i < size; i++)
    {
        for (int a = 0; a < matsize; a++)
        {
            recieve_array2[a][i] = work_to_do(i,a);
        }
    }

    for (int a = 0; a < matsize; a++)
    {
        displacement_array2[a][0] = 0;
        for (int i = 1; i < size; i++)
        {
            displacement_array2[a][i] = displacement_array2[a][i-1] + recieve_array2[a][i-1];
        }
    }

    mat recieve_matrix = zeros(matsize, matsize); // Recieve information from Hatree Fock object
    double *send_cube_MPI = (double*) malloc(cube_send_size*sizeof(double)); // Use this to store calculations performed on a node
    double *recieve_cube_MPI = (double*) malloc(ccube_send_size*sizeof(double)); // Use this to send and recieve information between nodes

    // We use .t() on c. No reason for this, just preferance
    c = c.t();


    /*
    // Method 2, not used alternative algorithm... ignore this if
    // you are not trying to optimize further and wish to apply
    // principles stated in text

    mat calc1 = zeros(matsize,matsize);
    field<mat> temp1;
    temp1.set_size(matsize,matsize);
    field<mat> temp2;
    temp2.set_size(matsize,matsize);

    for (int i = 0; i < matsize; i++)
    {
        for (int j = 0; j <= i; j++)
        {
            temp1(i,j) = calc1;
            temp2(i,j) = calc1;
        }
    }
    */


    clock_t start1; // Measure time for AO -> MO transfomation
    clock_t slutt1;
    double time_ao_mo;
    start1 = clock();

    for (int a = nr_freeze; a < matsize; a++)
    {
        INDEX_CHECK3 = 0;
        for (int j = 0; j < matsize; j++)
        {
            for (int l = where_to_start2(j,rank); l <= j; l+=size)
            {
                // Get AOs
                recieve_matrix = HartFock->Return_Field_Q(j, l);

                // First quarter transformation
                compact_mo10 = c.row(a) * recieve_matrix;

                // Secound quarter transformation
                compact_mo11(span(0,a)) = c(span(0,a),span()) * compact_mo10.t();

                // Make it one dimensional array for MPI communication
                for (int b = 0; b <= a; b++)
                {
                    send_cube_MPI[INDEX_CHECK3] = compact_mo11(b);
                    INDEX_CHECK3 += 1;
                }
            }
        }

        // Communicate data with all-to-all communication
        MPI_Allgatherv(send_cube_MPI, work_do_2(rank,a), MPI_DOUBLE, recieve_cube_MPI,
                       recieve_array[a], displacement_array[a],
                       MPI_DOUBLE, MPI_COMM_WORLD);

        INDEX_CHECK3 = 0;
        for (int K = 0; K < size; K++)
        {
            for (int i = 0; i < matsize; i++)
            {
                for (int j = where_to_start2(i,K); j <= i; j+=size)
                {
                    for (int k = 0; k <= a; k++)
                    {
                        // Flip indexes at same time as mapping, to facilitate matrix multiplications
                        compact_mo5.at(k)(j,i) = recieve_cube_MPI[INDEX_CHECK3];
                        compact_mo5.at(k)(i,j) = recieve_cube_MPI[INDEX_CHECK3]; // Symmetry
                        INDEX_CHECK3 += 1;
                    }
                }
            }
        }

        // Continue in parallel, new grid now of g and b
        INDEX_CHECK3 = 0;
        for (int b = 0; b <= a; b++)
        {
            for (int g = where_to_start(b, rank); g < matsize; g+=size)
            {
                // Third quarter transform
                compact_mo10 = c.row(g) * compact_mo5.at(b);

                // Fourth quarter transform
                compact_mo11(span(0,g)) = c(span(0,g), span()) * compact_mo10.t();

                // Make one dimensioal array
                for (int h = 0; h <= g; h++)
                {
                    send_cube_MPI[INDEX_CHECK3] = compact_mo11(h);
                    INDEX_CHECK3 += 1;
                }
            }
        }

        // Communicate data with all-to-all communication
        MPI_Allgatherv(send_cube_MPI, work_to_do(rank,a), MPI_DOUBLE, recieve_cube_MPI,
                       recieve_array2[a], displacement_array2[a],
                       MPI_DOUBLE, MPI_COMM_WORLD);

        INDEX_CHECK3 = 0;
        for (int K = 0; K < size; K++)
        {
            for (int i = 0; i <= a; i++)
            {
                for (int k = where_to_start(i,K); k < matsize; k+=size)
                {
                    for (int j = 0; j <= k; j++)
                    {
                        compact_mo5.at(i)(j,k) = recieve_cube_MPI[INDEX_CHECK3];
                        compact_mo5.at(i)(k,j) = recieve_cube_MPI[INDEX_CHECK3];
                        INDEX_CHECK3 += 1;
                    }
                }
            }
        }



        // Lets distribute once again, now even a new grid (third grid contained in this function)
        // This grid is optimized for CCSD calculations later, and have nothing to do
        // with the AOtoMO transformation.
        // We just store the values we need as we need them for later calculations.
        if (a >= nr_freeze && a < Speed_Elec)
        {
            // Using a-b symmetry
            for (int b = nr_freeze; b <= a; b++)
            {
                for (int i = a; i < Speed_Elec; i++)
                {
                    for (int j = Speed_Elec; j < matsize; j++)
                    {
                        MO6_2(a-nr_freeze, i-nr_freeze)(b-nr_freeze,j-Speed_Elec) = compact_mo5.at(b)(i,j);
                    }
                }
            }

            // Using a-b symmetry
            for (int b = nr_freeze; b <= a; b++)
            {
                for (int i = b; i < Speed_Elec; i++)
                {
                    for (int j = Speed_Elec; j < matsize; j++)
                    {
                        MO6_2(b-nr_freeze, i-nr_freeze)(a-nr_freeze,j-Speed_Elec) = compact_mo5.at(b)(i,j);
                    }
                }
            }

            // Using a-b symmetry
            for (int i = a; i < Speed_Elec; i++)
            {
                for (int b = nr_freeze; b <= a; b++)
                {
                    for (int j = nr_freeze; j < Speed_Elec; j++)
                    {
                        MO8(a-nr_freeze,i-nr_freeze)(b-nr_freeze,j-nr_freeze) = compact_mo5.at(b)(i,j);
                    }
                }
            }

            // Using a-b symmetry
            for (int b = nr_freeze; b <= a; b++)
            {
                for (int i = b; i < Speed_Elec; i++)
                {
                    for (int j = nr_freeze; j < Speed_Elec; j++)
                    {
                        MO8(b-nr_freeze,i-nr_freeze)(a-nr_freeze,j-nr_freeze) = compact_mo5.at(b)(i,j);
                    }
                }
            }
        }


        else if (a >= Speed_Elec)
        {
            // Using a-b symmetry
            for (int b = nr_freeze; b < Speed_Elec; b++)
            {
                for (int i = Speed_Elec; i < matsize; i++)
                {
                    for (int j = nr_freeze; j < Speed_Elec; j++)
                    {
                        MOLeftovers(a-Speed_Elec,i-Speed_Elec)(b-nr_freeze,j-nr_freeze) = compact_mo5.at(b)(i,j);
                    }
                }
            }

            // Using a-b symmetry
            for (int b = nr_freeze; b < Speed_Elec; b++)
            {
                for (int i = b; i < Speed_Elec; i++)
                {
                    for (int j = nr_freeze; j < Speed_Elec; j++)
                    {
                        MO6_1(b-nr_freeze,i-nr_freeze)(a-Speed_Elec,j-nr_freeze) = compact_mo5.at(b)(i,j); //compact_mo2(b,i,j);
                    }
                }
            }

            // Using a-b symmetry
            for (int b = nr_freeze; b < Speed_Elec; b++)
            {
                for (int i = b; i < Speed_Elec; i++)
                {
                    for (int j = Speed_Elec; j < matsize; j++)
                    {
                        MO4(b-nr_freeze,i-nr_freeze)(a-Speed_Elec,j-Speed_Elec) = compact_mo5.at(b)(i,j); //compact_mo2(b,i,j);
                    }
                }
            }

            A = a - Speed_Elec;

            // Using a-b symmetry
            for (int i = local_dist1(a); i < matsize; i+=size) // (0 -> Speed_Occ) is equal to (a -> matsize) - Speed_Elec
            {
                for (int b = nr_freeze; b <= a; b++)
                {
                    for (int j = nr_freeze; j < matsize; j++)
                    {
                        MO9(a-Speed_Elec,i-Speed_Elec)(b-nr_freeze,j-nr_freeze) = compact_mo5.at(b)(i,j); // compact_mo2(b,i,j);
                    }
                }
            }

            // Using a-b symmetry
            for (int b = Speed_Elec; b <= a; b++)
            {
                for (int i = local_dist1(b); i < matsize; i+=size) // (0 -> Speed_Occ) is equal to (a -> matsize) - Speed_Elec
                {
                    for (int j = nr_freeze; j < matsize; j++)
                    {
                        MO9(b-Speed_Elec,i-Speed_Elec)(a-nr_freeze,j-nr_freeze) = compact_mo5.at(b)(i,j);
                    }
                }
            }

            // Using a-b symmetry
            for (int j = local_dist2(a); j < Speed_Elec; j+=size)
            {
                for (int b = Speed_Elec; b <= a; b++)
                {
                    for (int i = nr_freeze; i < Speed_Elec; i++)
                    {
                        MO2(a-Speed_Elec,j-nr_freeze)(b-Speed_Elec,i-nr_freeze) = compact_mo5.at(b)(i,j);
                    }
                }
            }

            // Using a-b symmetry
            for (int b = Speed_Elec; b <= a; b++)
            {
                for (int j = local_dist2(b); j < Speed_Elec; j+=size)
                {
                    for (int i = nr_freeze; i < Speed_Elec; i++)
                    {
                        MO2(b-Speed_Elec, j-nr_freeze)(a-Speed_Elec,i-nr_freeze) = compact_mo5.at(b)(i,j);
                    }
                }
            }

            // Using a-b symmetry
            for (int i = local_dist2(a); i < Speed_Elec; i+=size)
            {
                for (int b = Speed_Elec; b < matsize; b++)
                {
                    for (int j = nr_freeze; j < Speed_Elec; j++)
                    {
                        MO3(a-Speed_Elec,i-nr_freeze)(j-nr_freeze,b-Speed_Elec) = compact_mo5.at(i)(j,b);
                    }
                }
            }

            // Using a-b symmetry
            for (int i = local_dist2(a); i < Speed_Elec; i+=size) // <-- i
            {
                for (int b = nr_freeze; b <= a; b++)
                {
                    for (int j = nr_freeze; j < matsize; j++)
                    {
                        MO10(a-Speed_Elec,i-nr_freeze)(b-nr_freeze,j-nr_freeze) = compact_mo5.at(b)(i,j);
                    }
                }
            }

            // Using a-b symmetry
            for (int b = Speed_Elec; b <= a; b++)
            {
                for (int i = local_dist2(b); i < Speed_Elec; i+=size) // <-- i
                {
                    for (int j = nr_freeze; j < matsize; j++)
                    {
                        MO10(b-Speed_Elec,i-nr_freeze)(a-nr_freeze,j-nr_freeze) = compact_mo5.at(b)(i,j);
                    }
                }
            }
        }
    }




/*
    // Method 2:
    // Ignore this if you are not trying to optimize further.


    // QT1 & QT2
    for (int j = 0; j < matsize; j++)
    {
        for (int k = 0; k <= j; k++)
        {
            recieve_matrix = HartFock->Return_Field_Q(j, k);

            // QT1
            calc1 = c * recieve_matrix;
            for (int a = 0; a < matsize; a++)
            {
                // QT2
                temp1(j,k)(a,span(0,a)) = c(span(0,a),span()) * calc1.row(a);
            }
        }
    }

    // Communication
    for (int a = 0; a < matsize; a++)
    {
        for (int b = 0; b <= a; b++)
        {
            for (int j = 0; j < matsize; j++)
            {
                for (int k = 0; k <= j; k++)
                {
                    temp2(a,b)(j,k) = temp1(j,k)(a,b);
                    temp2(a,b)(k,j) = temp1(j,k)(a,b);
                }
            }
        }
    }

    // QT3 & QT4
    for (int a = 0; a < matsize; a++)
    {
        for (int b = 0; b <= a; b++)
        {
            // QT3
            calc1 = c * temp2(a,b);
            for (int g = 0; g < matsize; g++)
            {
                // QT4
                temp1(a,b)(g, span(0,g)) = c(span(0,g),span()) * calc1.row(g);
            }
        }
    }

    // New distribution
    for (int a = 0; a < matsize; a++)
    {
        for (int b = 0; b <= a; b++)
        {
            for (int j = 0; j < matsize; j++)
            {
                for (int k = 0;k <= j; k++)
                {
                    compact_mo5.at(b)(j,k) = temp1(a,b)(j,k);
                    compact_mo5.at(b)(k,j) = temp1(a,b)(j,k);
                }
            }
        }


        if (a >= nr_freeze && a < Speed_Elec)
        {
            // Using a-b symmetry
            for (int b = nr_freeze; b <= a; b++)
            {
                for (int i = a; i < Speed_Elec; i++)
                {
                    for (int j = Speed_Elec; j < matsize; j++)
                    {
                        MO6_2(a-nr_freeze, i-nr_freeze)(b-nr_freeze,j-Speed_Elec) = compact_mo5.at(b)(i,j);
                    }
                }
            }

            // Using a-b symmetry
            for (int b = nr_freeze; b <= a; b++)
            {
                for (int i = b; i < Speed_Elec; i++)
                {
                    for (int j = Speed_Elec; j < matsize; j++)
                    {
                        MO6_2(b-nr_freeze, i-nr_freeze)(a-nr_freeze,j-Speed_Elec) = compact_mo5.at(b)(i,j);
                    }
                }
            }

            // Using a-b symmetry
            for (int i = a; i < Speed_Elec; i++)
            {
                for (int b = nr_freeze; b <= a; b++)
                {
                    for (int j = nr_freeze; j < Speed_Elec; j++)
                    {
                        MO8(a-nr_freeze,i-nr_freeze)(b-nr_freeze,j-nr_freeze) = compact_mo5.at(b)(i,j);
                    }
                }
            }

            // Using a-b symmetry
            for (int b = nr_freeze; b <= a; b++)
            {
                for (int i = b; i < Speed_Elec; i++)
                {
                    for (int j = nr_freeze; j < Speed_Elec; j++)
                    {
                        MO8(b-nr_freeze,i-nr_freeze)(a-nr_freeze,j-nr_freeze) = compact_mo5.at(b)(i,j);
                    }
                }
            }
        }

        else if (a >= Speed_Elec)
        {
            // Using a-b symmetry
            for (int b = nr_freeze; b < Speed_Elec; b++)
            {
                for (int i = Speed_Elec; i < matsize; i++)
                {
                    for (int j = nr_freeze; j < Speed_Elec; j++)
                    {
                        MOLeftovers(a-Speed_Elec,i-Speed_Elec)(b-nr_freeze,j-nr_freeze) = compact_mo5.at(b)(i,j);
                    }
                }
            }

            // Using a-b symmetry
            for (int b = nr_freeze; b < Speed_Elec; b++)
            {
                for (int i = b; i < Speed_Elec; i++)
                {
                    for (int j = nr_freeze; j < Speed_Elec; j++)
                    {
                        MO6_1(b-nr_freeze,i-nr_freeze)(a-Speed_Elec,j-nr_freeze) = compact_mo5.at(b)(i,j); //compact_mo2(b,i,j);
                    }
                }
            }

            // Using a-b symmetry
            for (int b = nr_freeze; b < Speed_Elec; b++)
            {
                for (int i = b; i < Speed_Elec; i++)
                {
                    for (int j = Speed_Elec; j < matsize; j++)
                    {
                        MO4(b-nr_freeze,i-nr_freeze)(a-Speed_Elec,j-Speed_Elec) = compact_mo5.at(b)(i,j); //compact_mo2(b,i,j);
                    }
                }
            }

            A = a - Speed_Elec;

            // Using a-b symmetry
            for (int i = local_dist1(a); i < matsize; i+=size) // (0 -> Speed_Occ) is equal to (a -> matsize) - Speed_Elec
            {
                for (int b = nr_freeze; b <= a; b++)
                {
                    for (int j = nr_freeze; j < matsize; j++)
                    {
                        MO9(a-Speed_Elec,i-Speed_Elec)(b-nr_freeze,j-nr_freeze) = compact_mo5.at(b)(i,j); // compact_mo2(b,i,j);
                    }
                }
            }

            // Using a-b symmetry
            for (int b = Speed_Elec; b <= a; b++)
            {
                for (int i = local_dist1(b); i < matsize; i+=size) // (0 -> Speed_Occ) is equal to (a -> matsize) - Speed_Elec
                {
                    for (int j = nr_freeze; j < matsize; j++)
                    {
                        MO9(b-Speed_Elec,i-Speed_Elec)(a-nr_freeze,j-nr_freeze) = compact_mo5.at(b)(i,j);
                    }
                }
            }

            // Using a-b symmetry
            for (int j = local_dist2(a); j < Speed_Elec; j+=size)
            {
                for (int b = Speed_Elec; b <= a; b++)
                {
                    for (int i = nr_freeze; i < Speed_Elec; i++)
                    {
                        MO2(a-Speed_Elec,j-nr_freeze)(b-Speed_Elec,i-nr_freeze) = compact_mo5.at(b)(i,j);
                    }
                }
            }

            // Using a-b symmetry
            for (int b = Speed_Elec; b <= a; b++)
            {
                for (int j = local_dist2(b); j < Speed_Elec; j+=size)
                {
                    for (int i = nr_freeze; i < Speed_Elec; i++)
                    {
                        MO2(b-Speed_Elec, j-nr_freeze)(a-Speed_Elec,i-nr_freeze) = compact_mo5.at(b)(i,j);
                    }
                }
            }

            // Using a-b symmetry
            for (int i = local_dist2(a); i < Speed_Elec; i+=size)
            {
                for (int b = Speed_Elec; b < matsize; b++)
                {
                    for (int j = nr_freeze; j < Speed_Elec; j++)
                    {
                        MO3(a-Speed_Elec,i-nr_freeze)(j-nr_freeze,b-Speed_Elec) = compact_mo5.at(i)(j,b);
                    }
                }
            }

            // Using a-b symmetry
            for (int i = local_dist2(a); i < Speed_Elec; i+=size) // <-- i
            {
                for (int b = nr_freeze; b <= a; b++)
                {
                    for (int j = nr_freeze; j < matsize; j++)
                    {
                        MO10(a-Speed_Elec,i-nr_freeze)(b-nr_freeze,j-nr_freeze) = compact_mo5.at(b)(i,j);
                    }
                }
            }

            // Using a-b symmetry
            for (int b = Speed_Elec; b <= a; b++)
            {
                for (int i = local_dist2(b); i < Speed_Elec; i+=size) // <-- i
                {
                    for (int j = nr_freeze; j < matsize; j++)
                    {
                        MO10(b-Speed_Elec,i-nr_freeze)(a-nr_freeze,j-nr_freeze) = compact_mo5.at(b)(i,j);
                    }
                }
            }
        }

    }

*/


    slutt1 = clock();
    time_ao_mo = (double) (slutt1-start1)/CLOCKS_PER_SEC;

    if (rank == 0)
    {
        cout << "Time for AO->MO transfomation: " <<  time_ao_mo << endl;
    }

    c.clear();

    // Ignore the electrons that are frozen in future calculations
    // This is how to freeze core in CCSD... but we have not finished implementation
    // of frozen core orbitals.
    n_Electrons -= 2*nr_freeze;
    Speed_Elec -= nr_freeze;
}

long int CCSD_Memory_optimized::Return_Integral_Index(int a, int b, int i, int j)
{
    // Function not in use anymore
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

mat CCSD_Memory_optimized::Fill_FS(vec eigval, int nr_freeze)
{
    // Fill fock eigenvalues into a diagonal matrix
    int frozen = 2 * nr_freeze;
    mat fss = zeros(Matrix_Size-frozen, Matrix_Size-frozen);

    for (int i = frozen; i < Matrix_Size; i++)
    {
        fss(i-frozen,i-frozen) = eigval(i/2); // Diagonal matrise
    }

    return fss;
}

void CCSD_Memory_optimized::Fill_integ2_2D(int a, int i)
{
    // 2D Mapping for External Math Library use
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

        // lower half 0

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
    // 2D Mapping for External Math Library use
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
                integ3_2D(B,J) = MO9(A, I)(B+Speed_Elec, (J+Speed_Elec)) - MO9(A, I)(J+Speed_Elec, (B+Speed_Elec));
                j++;
                J++;
            }
            b++;
            B++;
        }

        // Lower half of matrix = 0, apparantly slower
        //integ3_2D(span(Speed_Occ, unocc_orb-1), span()) = Zero_Matrix;
/*
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
*/

    }

    else if (a%2 == 1 && i%2 == 1)
    {
        // Upper half of matrix = 0
        //integ3_2D(span(0, Speed_Occ-1), span()) = Zero_Matrix;

/*
        B = 0;
        for (int b = 0; b < unocc_orb; b++)
        {
            J = 0;
            for (int j = 0; j < unocc_orb; j++)
            {
                integ3_2D(B,J) = 0;
                j++;
                J++;
            }
            b++;
            B++;
        }
*/

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
    // 2D Mapping for External Math Library use
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

        // upper half 0
        /*
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
        */
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

        // upper half 0
        /*
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
        */

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
    // 2D Mapping for External Math Library use
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

        // lower half 0
        /*
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
        */
    }

    else if (a%2 == 1 && i%2 == 1)
    {
        // upper half 0
        /*
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
        */

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
    // 2D Mapping for External Math Library use
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

        // lower half 0
        /*
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
        */
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
        // lower half 0
        /*
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
        */

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
    // 2D Mapping for External Math Library use
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

        // lower half 0
        /*
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
        */
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

        // lower half 0
        /*
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
        */
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
    // 2D Mapping for External Math Library use
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

        // lower half 0
        /*
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
        */
    }

    else if(a%2 == 1 && i%2 == 1)
    {
        // upper half 0
        /*
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
        */

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
    // 2D Mapping for External Math Library use
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
                integ9_2D(B,J) = MO9(A, I)(B+Speed_Elec, J)-MO9(A, I)(J, (B+Speed_Elec));
                j++;
                J++;
            }
            b++;
            B++;
        }

        // Lower half 0
        /*
        B = 0;
        for (int b = 1; b < unocc_orb; b++)
        {
            J = 0;
            for (int j = 1; j < n_Electrons; j++)
            {
                integ9_2D(B+Speed_Occ,J) = 0;
                j++;
                J++;
            }
            b++;
            B++;
        }
        */
    }

    else if(a%2 == 1 && i%2 == 1)
    {
        // Upper half 0
        /*
        B = 0;
        for (int b = 0; b < unocc_orb; b++)
        {
            J = 0;
            for (int j = 0; j < n_Electrons; j++)
            {
                integ9_2D(B,J) = 0;
                j++;
                J++;
            }
            b++;
            B++;
        }
        */

        B = 0;
        for (int b = 1; b < unocc_orb; b++)
        {
            J = 0;
            for (int j = 1; j < n_Electrons; j++)
            {
                integ9_2D(B+Speed_Occ,J) = MO9(A, I)(B+Speed_Elec, J)-MO9(A, I)(J, (B+Speed_Elec));
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
            for (int j = 1; j < n_Electrons; j++)
            {
                integ9_2D(B,J) = MO9(A, I)(B+Speed_Elec, J);
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
                integ9_2D(B+Speed_Occ,J) = -MO9(A, I)(J, (B+Speed_Elec));
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
                integ9_2D(B,J) = -MO9(A, I)(J, (B+Speed_Elec));
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
                integ9_2D(B+Speed_Occ,J) = MO9(A, I)(B+Speed_Elec, J);
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
    // 2D Mapping for External Math Library use
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

void CCSD_Memory_optimized::Fill_W4_MPI()
{
    // Calculate W4 as given in the text.
    // Run this calculation in parallel

    int index_counter = 0;
   // int INDEX_CHECK;
   // int A, I, AA;
    double temp;

    for (int a = 0; a < unocc_orb; a++)
    {
       // A = a/2*Speed_Elec;
        for (int m = Where_To_Start_Part2(rank,a); m < n_Electrons; m+=jump)
        {
            //INDEX_CHECK = A + m/2;
            //if (INDEX_CHECK % size == rank)
            //{
                Fill_integ7_2D(a,m);
                Fill_integ5_2D(a,m);

                for (int e = 0; e < unocc_orb; e++)
                {
                    Fill_integ2_2D_even_even(e, m);
                    for(int i = 0; i < n_Electrons; i++)
                    {
                        Index_Swapping_W_4[a][m][e][i] = -integ7_2D(e/2,i/2) - accu(W_3.at(i,m)(e/2,span()) % T_1.row(a/2))
                                + accu(integ5_2D(span(0, Speed_Occ-1),e/2) % T_1(span(0, Speed_Occ-1),i/2))
                                + 0.5*accu(integ2_2D % t2.at(a,i));
                        i++;
                    }
                    e++;
                }
            //}
            //m++;
        }
        a++;
    }

    for (int a = 1; a < unocc_orb; a++)
    {
       // A = a/2*Speed_Elec;
        for (int m = Where_To_Start_Part2(rank,a-1); m < n_Electrons; m+=jump)
        {
           // INDEX_CHECK = A + m/2;
           // if (INDEX_CHECK % size == rank)
           // {

                Fill_integ7_2D(a,m);
                Fill_integ5_2D(a,m);

                for (int e = 0; e < unocc_orb; e++)
                {
                    Fill_integ2_2D_even_even(e, m);
                    for (int i = 1; i < n_Electrons; i++)
                    {
                        Index_Swapping_W_4[a][m][e][i] = temp - integ7_2D(e/2,i/2) + 0.5*accu(integ2_2D % t2.at(a,i))
                                - accu(W_3.at(i,m)(e/2,span()) % T_1.row(a/2+Speed_Occ))
                                + accu(integ5_2D(span(Speed_Occ, Speed_Occ+Speed_Occ-1),e/2) % T_1(span(Speed_Occ, Speed_Occ+Speed_Occ-1),i/2));
                        i++;
                    }
                    e++;
                }

                for (int e = 1; e < unocc_orb; e++)
                {
                    Fill_integ2_2D_odd_even(e, m);
                    for (int i = 0; i < n_Electrons; i++)
                    {
                        Index_Swapping_W_4[a][m][e][i] = -integ7_2D(e/2+Speed_Occ,i/2) - accu(W_3.at(i,m)(e/2+Speed_Occ,span()) % T_1.row(a/2+Speed_Occ))
                                + accu(integ5_2D(span(0, Speed_Occ-1),e/2) % T_1(span(0, Speed_Occ-1),i/2))
                                + 0.5*accu(integ2_2D(span(0, Speed_Occ-1), span()) % t2.at(a,i)(span(0, Speed_Occ-1), span()));
                        i++;
                    }
                    e++;
                }
            //}
            //m++;
        }
        a++;
    }

    for (int a = 0; a < unocc_orb; a++)
    {
        //A = a/2*Speed_Elec; // m is odd number
        for (int m = Where_To_Start_Part2(rank,a)+1; m < n_Electrons; m+=jump)
        {
            //INDEX_CHECK = A + m/2;
            //if (INDEX_CHECK % size == rank)
           // {

                Fill_integ7_2D(a,m);
                Fill_integ5_2D(a,m);

                for (int e = 0; e < unocc_orb; e++)
                {
                    Fill_integ2_2D_even_odd(e, m);
                    for (int i = 1; i < n_Electrons; i++)
                    {
                        Index_Swapping_W_4[a][m][e][i] = -integ7_2D(e/2,i/2) - accu(W_3.at(i,m)(e/2,span()) % T_1.row(a/2))
                                + accu(integ5_2D(span(Speed_Occ, Speed_Occ+Speed_Occ-1),e/2) % T_1(span(Speed_Occ, Speed_Occ+Speed_Occ-1),i/2))
                                + 0.5*accu(integ2_2D(span(Speed_Occ, Speed_Occ+Speed_Occ-1), span()) % t2.at(a,i)(span(Speed_Occ, Speed_Occ+Speed_Occ-1), span()));
                        i++;
                    }
                    e++;
                }

                for (int e = 1; e < unocc_orb; e++)
                {
                    Fill_integ2_2D_odd_odd(e, m);
                    for (int i = 0; i < n_Electrons; i++)
                    {
                        Index_Swapping_W_4[a][m][e][i] = -integ7_2D(e/2+Speed_Occ,i/2) - accu(W_3.at(i,m)(e/2,span()) % T_1.row(a/2+a%2*Speed_Occ))
                                + accu(integ5_2D(span(0, Speed_Occ-1),e/2) % T_1(span(0, Speed_Occ-1),i/2))
                                + 0.5*accu(integ2_2D % t2.at(a,i));
                        i++;
                    }
                    e++;
                }
            //}
            //m++;
        }
        a++;
    }

    for (int a = 1; a < unocc_orb; a++)
    {
        //A = a/2*Speed_Elec;
        for (int m = Where_To_Start_Part2(rank,a-1)+1; m < n_Electrons; m+=jump)
        {
            //INDEX_CHECK = A + m/2;
            //if (INDEX_CHECK % size == rank)
            {

                Fill_integ7_2D(a,m);
                Fill_integ5_2D(a,m);

                for (int e = 1; e < unocc_orb; e++)
                {
                    Fill_integ2_2D_odd_odd(e, m);
                    for (int i = 1; i < n_Electrons; i++)
                    {
                        Index_Swapping_W_4[a][m][e][i] = -integ7_2D(e/2+Speed_Occ,i/2) - accu(W_3.at(i,m)(e/2+Speed_Occ,span()) % T_1.row(a/2+Speed_Occ))
                                + accu(integ5_2D(span(Speed_Occ, Speed_Occ+Speed_Occ-1),e/2) % T_1(span(Speed_Occ, Speed_Occ+Speed_Occ-1),i/2))
                                + 0.5*accu(integ2_2D % t2.at(a,i));
                        i++;
                    }
                    e++;
                }
            }
            //m++;
        }
        a++;
    }




    // 1D mapping to make sure we can use Scatter later
    // (done like this because we are also remapping to a different grid)
    index_counter = 0;

    for (int K = 0; K < size; K++)
    {
        // Arrange this so one node gets its information in sequence

        for (int a = 0; a < unocc_orb; a++)
        {
           // AA = a/2 * Speed_Elec;
            for(int i = Where_To_Start_Part2(K,a); i < n_Electrons; i+=jump)
            {
              //  I = i/2;
              //  if ((AA+I)%size == K)
               // {
                    for (int m =Where_To_Start_Part2(rank,a); m < n_Electrons; m+=jump)
                    {
                        //INDEX_CHECK = AA + m/2;
                        //if (INDEX_CHECK % size == rank)
                        //{
                            for (int e = 0; e < unocc_orb; e++)
                            {
                                MY_OWN_MPI[index_counter] = Index_Swapping_W_4[a][m][e][i];
                                index_counter += 1;
                                e++;
                            }
                        //}
                       // m++;
                    }

                    for (int m = Where_To_Start_Part2(rank,a)+1; m < n_Electrons; m+=jump)
                    {
                       // INDEX_CHECK = AA + m/2;
                       // if (INDEX_CHECK % size == rank)
                       // {
                            for (int e = 1; e < unocc_orb; e++)
                            {
                                MY_OWN_MPI[index_counter] = Index_Swapping_W_4[a][m][e][i];
                                index_counter += 1;
                                e++;
                            }
                       // }
                        //m++;
                    }
                //}
                //i++;
            }
            a++;
        }

        for (int a = 0; a < unocc_orb; a++)
        {
           // AA = a/2 * Speed_Elec;
            for (int i = Where_To_Start_Part2(K,a)+1; i < n_Electrons; i+=jump)
            {
               // I = i/2;
               // if ((AA+I)%size == K)
               // {
                    for (int m = Where_To_Start_Part2(rank,a)+1; m < n_Electrons; m+=jump)
                    {
                        //INDEX_CHECK = AA + m/2;
                        //if (INDEX_CHECK % size == rank)
                        //{
                            for (int e = 0; e < unocc_orb; e++)
                            {
                                MY_OWN_MPI[index_counter] = Index_Swapping_W_4[a][m][e][i];
                                index_counter += 1;
                                e++;
                            }
                        //}
                        //m++;
                    }
               // }
                //i++;
            }
            a++;
        }

        for (int a = 1; a < unocc_orb; a++)
        {
           // AA = a/2 * Speed_Elec;
            for (int i = Where_To_Start_Part2(K,a-1); i < n_Electrons; i+=jump)
            {
              //  I = i/2;
              //  if ((AA+I)%size == K)
              //  {
                    for (int m = Where_To_Start_Part2(rank,a-1); m < n_Electrons; m+=jump)
                    {
                        //INDEX_CHECK = AA + m/2;
                        //if (INDEX_CHECK % size == rank)
                        //{
                            for (int e = 1; e < unocc_orb; e++)
                            {
                                MY_OWN_MPI[index_counter] = Index_Swapping_W_4[a][m][e][i];
                                index_counter += 1;
                                e++;
                            }
                        //}
                        //m++;
                    }
//                }
               // i++;
            }
            a++;
        }

        for (int a = 1; a < unocc_orb; a++)
        {
          //  AA = a/2 * Speed_Elec;
            for (int i = Where_To_Start_Part2(K,a-1)+1; i < n_Electrons; i+=jump)
            {
               // I = i/2;
                //if ((AA+I)%size == K)
               // {
                    for (int m = Where_To_Start_Part2(rank,a-1); m < n_Electrons; m+=jump)
                    {
                      //  INDEX_CHECK = AA + m/2;
                      //  if (INDEX_CHECK % size == rank)
                      //  {
                            for (int e = 0; e < unocc_orb; e++)
                            {
                                MY_OWN_MPI[index_counter] = Index_Swapping_W_4[a][m][e][i];
                                index_counter += 1;
                                e++;
                            }
                       // }
                        //m++;
                    }

                    for (int m = Where_To_Start_Part2(rank,a-1)+1; m < n_Electrons; m+=jump)
                    {
                       // INDEX_CHECK = AA + m/2;
                       // if (INDEX_CHECK % size == rank)
                       // {
                            for (int e = 1; e < unocc_orb; e++)
                            {
                                MY_OWN_MPI[index_counter] = Index_Swapping_W_4[a][m][e][i];
                                index_counter += 1;
                                e++;
                            }
                       // }
                       // m++;
                    }
                //}
                //i++;
            }
            a++;
        }
    }
}

mat CCSD_Memory_optimized::MO_Grid_A_B(int frez)
{

    // Determine workload for each node in AO->MO transformation procedure, secound grid
    int matsize = Matrix_Size/2;
    mat WorkLoad = zeros(size, matsize);
    int index_counter;

    for (int K = 0; K < size; K++)
    {
        for (int a = 0; a < matsize; a++)
        {
            index_counter = 0;
            for (int i = 0; i <= a; i++)
            {
                for (int j = 0; j < matsize; j++)
                {
                    if ((i+j)%size==K)
                    {
                        for (int k = 0; k <= j; k++)
                        {
                            index_counter += 1;
                        }
                    }
                }
            }
            WorkLoad(K,a) = index_counter;
        }
    }

    return WorkLoad;
}

mat CCSD_Memory_optimized::MO_Grid_I_J(int frez)
{
    // Determine workload for each node in AO->MO transformation procedure, first grid
    int matsize = Matrix_Size/2;
    mat WorkLoad = zeros(size, matsize);
    int index_counter;
    int grid;

    for (int K = 0; K < size; K++)
    {
        for (int a = 0; a < matsize; a++)
        {
            grid = 0;
            index_counter = 0;
            for (int i = 0; i < matsize; i++)
            {
                for (int j = 0; j <= i; j++)
                {
                    if (grid%size == K)
                    {
                        for (int k = 0; k <= a; k++)
                        {
                            index_counter += 1;
                        }
                    }
                    grid += 1;
                }
            }
            WorkLoad(K,a) = index_counter;
        }
    }

    return WorkLoad;
}

void CCSD_Memory_optimized::Where_To_Start_On_What()
{
    // Function no longer in use

    // This functions goal is to handle the W_4 calculation effectively
    // W_4(a,i)(e,m)
    // We wish to transform ourselves from one grid to another, without any additional calculations
    // This new grid is over a and i
    // And we also wish to replace MPI_Bcast with MPI_Scatterv by sorting our one dimensioanl distributed array
    // into an array where rows a and cols i are consecutive
/*
    int B;
    Displacement_W_4 = (int*)malloc(size * sizeof(int));
    Start_Pos = zeros(unocc_orb, size);
    Double_Size = 2 * size; // Twice the size of size! Will be used later for speedup

    // This function defines first where to start (i) at a fixed index (a)

    for (int a = 0; a < unocc_orb; a++)
    {
        B = 0;
        while ((a+B)%size != rank)
        {
            B++;
        }
        Start_Pos(a, rank) = B;
    }
    */
}

int CCSD_Memory_optimized::Calc_sum_a_n(int a)
{
    // Figure out how to get optimal work distribution in situation were
    // a matrix of jobs to distribute is symmetric and also jobs on diagonal
    // can be ignored and dodged and not performed

    int sum = 0;
    for (int i = 0; i <= a; i++)
    {
        sum += i;
    }
    return sum;
}

int CCSD_Memory_optimized::Freeze_How_Meny_Orbitals()
{
    int val = 0;

    // This functions returns val, which is how meny orbitals we freeze.
    // Freeze in this setting means ignore the contribution to the wavefunction from excitations of these orbitals
    // Freeze core not implemented jet

    if (freeze_core == true)
    {
        for (int i = 0; i < n_Electrons; i++)
        {
            // For H and He we dont do any core

            // For atoms nr 3 - 10 we ignore one contracted GTO
            if (Z(i) > 2 && Z(i) < 11)
            {
                val += 1; // Give values in number of SPACIAL orbitals to freeze
            }

            // For atoms nr 11 - ??? we ignore ??? contracted GTOs

            // -0.2473 // Test calculations, not to bad
            // -0.269453
        }
    }

    return val;
}

void CCSD_Memory_optimized::Distribute_Part1()
{
    // Distribute calculated terms in part 1 of calculation
    // Puts the values in correct arrays for external math library use
    int K, L, E, I, J;
    double temp;

    int index_counter;
    int INDEX_CHECK;

    int KK;

    // all to all optimized MPI communication
    MPI_Allgatherv(MY_OWN_MPI, Work_Each_Node_part1_Parallel[rank], MPI_DOUBLE,
                   SHARED_INFO_MPI, Work_Each_Node_part1_Parallel, Displacement_Each_Node_part1_Parallel,
                   MPI_DOUBLE, MPI_COMM_WORLD);
    index_counter = 0;


    // Loop over CPUs, since all CPUs have contributed in this calculation
    for (int X = 0; X < size; X++)
    {
        K = 0;
        for (int k = 0; k < n_Electrons; k++)
        {
            KK = K*Speed_Elec;
            L = 0;
            for (int l = 1; l < k; l++)
            {
                // Check which CPU sent the value we need
                INDEX_CHECK = KK+L+Local_Displacement1;
                if (INDEX_CHECK % size == X)
                {

                    for (int i = 0; i < n_Electrons; i++)
                    {
                        D2(K,i/2) += SHARED_INFO_MPI[index_counter];
                        i++;
                        index_counter += 1;
                    }

                    I = 0;
                    for (int i = 0; i < n_Electrons; i++)
                    {
                        E = 0;
                        for (int e = 1; e < unocc_orb; e++)
                        {
                            W_3(i,k)(E+Speed_Occ,L) = SHARED_INFO_MPI[index_counter];
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
                            W_3(i,k)(E,L) = SHARED_INFO_MPI[index_counter];
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
                            temp = SHARED_INFO_MPI[index_counter];

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
                            temp = SHARED_INFO_MPI[index_counter];

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
                INDEX_CHECK = KK+L+Local_Displacement1;
                if (INDEX_CHECK % size == X)
                {

                    for (int i = 0; i < n_Electrons; i++)
                    {
                        D2.at(K,i/2) -= SHARED_INFO_MPI[index_counter];
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
            KK = K*Speed_Elec;
            L = 0;
            for (int l = 0; l < k; l++)
            {

                INDEX_CHECK = KK+L+Local_Displacement1;
                if (INDEX_CHECK % size == X)
                {
                    I = 0;
                    for (int i = 0; i < n_Electrons; i++)
                    {
                        E = 0;
                        for (int e = 1; e < unocc_orb; e++)
                        {
                            W_3(i,k)(E,L) = SHARED_INFO_MPI[index_counter];
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
                            W_3(i,k)(E,L) = SHARED_INFO_MPI[index_counter];
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
                            temp = SHARED_INFO_MPI[index_counter];

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
                            temp = SHARED_INFO_MPI[index_counter];

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
            k++;
            K++;
        }

        K = 0;
        for (int k = 0; k < n_Electrons; k++)
        {
            KK = K *Speed_Elec;
            L = 0;
            for (int l = 0; l < k+1; l++)
            {

                INDEX_CHECK = KK+L+Local_Displacement1;
                if (INDEX_CHECK % size == X)
                {

                    for (int i = 0; i < n_Electrons; i++)
                    {
                        D2.at(K,i/2) += SHARED_INFO_MPI[index_counter];
                        i++;
                        index_counter += 1;
                    }

                    I = 0;
                    for (int i = 0; i < n_Electrons; i++)
                    {
                        E = 0;
                        for (int e = 0; e < unocc_orb; e++)
                        {
                            W_3(i,k)(E,L) = SHARED_INFO_MPI[index_counter];
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
                            temp = SHARED_INFO_MPI[index_counter];

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
                INDEX_CHECK = KK+L+Local_Displacement1;
                if (INDEX_CHECK % size == X)
                {
                    for (int i = 0; i < n_Electrons; i++)
                    {
                        D2.at(K,i/2) -= SHARED_INFO_MPI[index_counter];
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
            KK = K*Speed_Elec;
            L = 0;
            for (int l = 1; l < k+1; l++)
            {
                INDEX_CHECK = KK+L+Local_Displacement1;
                if (INDEX_CHECK % size == X)
                {/*
                    for (int i = 0; i < n_Electrons; i++)
                    {
                        D2.at(K+Speed_Elec,i/2) += SHARED_INFO_MPI[index_counter];
                        i++;
                        index_counter += 1;
                    }
                    */

                    I = 0;
                    for (int i = 1; i < n_Electrons; i++)
                    {
                        E = 0;
                        for (int e = 1; e < unocc_orb; e++)
                        {
                            W_3(i,k)(E+Speed_Occ,L) = SHARED_INFO_MPI[index_counter];
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
                            temp = SHARED_INFO_MPI[index_counter];

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
                INDEX_CHECK = a/2* Speed_Elec +m/2;
                if (INDEX_CHECK % size == X)
                {
                    for (int e = 0; e < unocc_orb; e++)
                    {
                        D3(a, e/2) += SHARED_INFO_MPI[index_counter];
                        e++;
                        index_counter += 1;
                    }

                    I = 0;
                    for (int i = 0; i < n_Electrons; i++)
                    {
                        J = I+1;
                        for (int j = i+2; j < n_Electrons; j++)
                        {
                            W_2(i,j)(A,M) += SHARED_INFO_MPI[index_counter];
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
                INDEX_CHECK = a/2* Speed_Elec +m/2;
                if (INDEX_CHECK % size == X)
                {/*
                    for (int e = 1; e < unocc_orb; e++)
                    {
                        D3(a, e/2) += SHARED_INFO_MPI[index_counter];
                        index_counter += 1;
                        e++;
                    }
                    */

                    I = 0;
                    for (int i = 0; i < n_Electrons; i++)
                    {
                        J = I+1;
                        for (int j = i+2; j < n_Electrons; j++)
                        {

                            W_2(i,j)(A+Speed_Occ,M) += SHARED_INFO_MPI[index_counter];
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

                            W_2(i,j)(A+Speed_Occ,M) += SHARED_INFO_MPI[index_counter];
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
                INDEX_CHECK = a/2* Speed_Elec + m/2;

                if (INDEX_CHECK % size == X)
                {/*
                    for (int e = 1; e < unocc_orb; e++)
                    {
                        D3(a, e/2) += SHARED_INFO_MPI[index_counter];
                        index_counter += 1;
                        e++;
                    }
                    */

                    I = 0;
                    for (int i = 0; i < n_Electrons; i++)
                    {
                        J = I;
                        for (int j = i+1; j < n_Electrons; j++)
                        {

                            W_2(i,j)(A+Speed_Occ,M) += SHARED_INFO_MPI[index_counter];
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

                            W_2(i,j)(A+Speed_Occ,M) += SHARED_INFO_MPI[index_counter];
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
                INDEX_CHECK = a/2* Speed_Elec + m/2;
                if (INDEX_CHECK % size == X)
                {
                    for (int e = 0; e < unocc_orb; e++)
                    {
                        D3(a,e/2) += SHARED_INFO_MPI[index_counter];
                        index_counter += 1;
                        e++;
                    }

                    I = 0;
                    for (int i = 0; i < n_Electrons; i++)
                    {
                        J = I;
                        for (int j = i+1; j < n_Electrons; j++)
                        {

                            W_2(i,j)(A,M) += SHARED_INFO_MPI[index_counter];
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

                            W_2(i,j)(A,M) += SHARED_INFO_MPI[index_counter];
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
                INDEX_CHECK = e/2 * Speed_Elec+ M;
                if (INDEX_CHECK % size == X)
                {
                    I = 0;
                    for (int i = 0; i < n_Electrons; i++)
                    {
                        D2(M,I) += SHARED_INFO_MPI[index_counter];
                        index_counter += 1;
                        i++;
                        I++;
                    }
                }
                m++;
                M++;
            }
            e++;
        }

        for (int e = 1; e < unocc_orb; e++)
        {
            M = 0;
            for (int m = 0; m < n_Electrons; m++)
            {
                INDEX_CHECK = e/2* Speed_Elec+ M;
                if (INDEX_CHECK % size == X)
                {
                    I = 0;
                    for (int i = 0; i < n_Electrons; i++)
                    {
                        D2(M,I) += SHARED_INFO_MPI[index_counter];
                        index_counter += 1;
                        i++;
                        I++;
                    }
                }
                m++;
                M++;
            }
            e++;
        }

        // Seperation noted here

        for (int e = 0; e < unocc_orb; e++)
        {
            for (int m = 0; m < n_Electrons; m++)
            {
                INDEX_CHECK = e/2* Speed_Elec + m/2;
                if (INDEX_CHECK % size == X)
                {
                    D1(e/2, m/2) += SHARED_INFO_MPI[index_counter];
                    index_counter += 1;

                    for (int a = 0; a < unocc_orb; a++)
                    {
                        D3.at(a,e/2) -= SHARED_INFO_MPI[index_counter];
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
                INDEX_CHECK = e/2* Speed_Elec + m/2;
                if (INDEX_CHECK % size == X)
                {
                    for (int a = 0; a < unocc_orb; a++)
                    {
                        D3.at(a,e/2) -= SHARED_INFO_MPI[index_counter];
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
                INDEX_CHECK = a/2 * Speed_Elec+ k/2;
                if (INDEX_CHECK % size == X)
                {
                    for (int i = 0; i < k; i++)
                    {
                        T_1_new(a/2, i/2) += SHARED_INFO_MPI[index_counter];
                        index_counter += 1;
                        i++;
                    }

                    // i will be less than k when k = odd number
                    // We want i to be an even number since a is even number, hence we start at i = k+1
                    for (int i = (k+1+(k+1)%2); i < n_Electrons; i++)
                    {
                        T_1_new(a/2, i/2) += SHARED_INFO_MPI[index_counter];
                        index_counter += 1;
                        i++;
                    }
                }
            }
            a++;
        }
    }

    // Map out remaining symmetries
    for (int k = 0; k < n_Electrons; k++)
    {
        K = k/2;
        for (int l = k+1; l < n_Electrons; l++)
        {
            L = l/2;

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

            l++;
        }
        k++;
    }

    for (int k = 1; k < n_Electrons; k++)
    {
        K = k/2;
        for (int l = k+1; l < n_Electrons; l++)
        {
            L = l/2;

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

            l++;
        }
        k++;
    }
}

void CCSD_Memory_optimized::Fill_integ2_2D_even_even(int a, int i)
{
    // 2D Mapping for External Math Library use
    int B,J;
    int A,I;
    A = a/2;
    I = i/2;

        B = 0;
        for (int b = 0; b < unocc_orb; b++)
        {
            J = 0;
            for (int j = 0; j < n_Electrons; j++)
            {
                integ2_2D(B,J) = -MOLeftovers(A, B)(J, I)+MOLeftovers(A, B)(I, J);
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
                integ2_2D(B+Speed_Occ,J) = MOLeftovers(A, B)(I, J);
                j++;
                J++;
            }
            b++;
            B++;
        }

}

void CCSD_Memory_optimized::Fill_integ2_2D_even_odd(int a, int i)
{

    // 2D Mapping for External Math Library use
    int B,J;
    int A,I;
    A = a/2;
    I = i/2;


    // upper half 0
    /*
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
        */

        B = 0;
        for (int b = 1; b < unocc_orb; b++)
        {
            J = 0;
            for (int j = 1; j < n_Electrons; j++)
            {
                integ2_2D(B+Speed_Occ,J) = -MOLeftovers(A, B)(J, I);
                j++;
                J++;
            }
            b++;
            B++;
        }

}

void CCSD_Memory_optimized::Fill_integ2_2D_odd_even(int a, int i)
{
    // 2D Mapping for External Math Library use
    int B,J;
    int A,I;
    A = a/2;
    I = i/2;
        B = 0;
        for (int b = 0; b < unocc_orb; b++)
        {
            J = 0;
            for (int j = 0; j < n_Electrons; j++) // ?
            {
                integ2_2D(B,J) = -MOLeftovers(A, B)(J, I);
                j++;
                J++;
            }
            b++;
            B++;
        }

        // lower half 0
        /*
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
        */
}

void CCSD_Memory_optimized::Fill_integ2_2D_odd_odd(int a, int i)
{
    // 2D Mapping for External Math Library use
    int B,J;
    int A,I;
    A = a/2;
    I = i/2;
        B = 0;
        for (int b = 0; b < unocc_orb; b++)
        {
            J = 0;
            for (int j = 0; j < n_Electrons; j++)
            {
                integ2_2D(B,J) = MOLeftovers(A, B)(I, J);
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
                integ2_2D(B+Speed_Occ,J) = MOLeftovers(A, B)(I, J) - MOLeftovers(A, B)(J, I);
                j++;
                J++;
            }
            b++;
            B++;
        }

}

void CCSD_Memory_optimized::Fill_W5()
{
    // Calculate P(ij) P(ab) W4 * t2
    // In parallel and all out memory distributed
    // Everything is distributed in memory except T2 amplitudes,
    // befoer and after calculations

    //int INDEX_CHECK;
    int index_counter;
    int A,B;//,I;
    int AA;

    int sum_a_n;
    index_counter = 0;

    for (int K =0 ; K < size; K++)
    {

        // All terms where a > b

        // even-even-even-even
        sum_a_n = 0;
        for (int a = 0; a < unocc_orb; a++)
        {
            sum_a_n += a/2;
            A = a/2 * Speed_Occ - sum_a_n;
            AA = a/2 * Speed_Elec;

            for (int b = a+2; b < unocc_orb; b++)
            {
                B = b/2;
                if ((A+B)%size == K)
                {
                    for (int i = Where_To_Start_Part2(rank,a); i < n_Electrons; i+=jump)
                    {
                       // if ((AA+i/2)%size == rank)
                       // {
                            for (int j = 0; j < n_Electrons; j++)
                            {
                                MY_OWN_MPI[index_counter] = accu(t2(b,j) % W_4(a,i));
                                index_counter += 1;
                                j++;
                            }
                       // }
                       // i++;
                    }
                }
                b++;
            }
            a++;
        }

        sum_a_n = 0;
        // even-odd -(even-odd / odd-even)
        for (int a = 0; a < unocc_orb; a++)
        {
            sum_a_n += a/2;
            A = a/2 * Speed_Occ - sum_a_n;
            AA = a/2 * Speed_Elec;

            for (int b = a+1; b < unocc_orb; b++)
            {
                B = b/2;

                if ((A+B)%size == K)
                {
                    for (int i = Where_To_Start_Part2(rank,a); i < n_Electrons; i+=jump)
                    {
                       // if ((AA+i/2)%size == rank)
                       // {
                            for (int j = 1; j < n_Electrons; j++)
                            {
                                MY_OWN_MPI[index_counter] = accu(t2.at(b,j) % W_4.at(a,i));
                                index_counter += 1;
                                j++;
                            }
                       // }
                       // i++;
                    }

                    for (int i = Where_To_Start_Part2(rank,a)+1; i < n_Electrons; i+=jump)
                    {
                       // if ((AA+i/2)%size == rank)
                       // {
                            for (int j = 0; j < n_Electrons; j++)
                            {
                                MY_OWN_MPI[index_counter] = accu(t2.at(b,j)(span(0, Speed_Occ-1), span()) % W_4.at(a,i)(span(0, Speed_Occ-1), span()));
                                index_counter += 1;
                                j++;
                            }
                       // }
                       // i++;
                    }
                }
                b++;
            }
            a++;
        }


        sum_a_n = 0;
        // odd-odd-odd-odd
        for (int a = 1; a < unocc_orb; a++)
        {
            sum_a_n += a/2;
            A = a/2 * Speed_Occ - sum_a_n;
            AA = a/2 * Speed_Elec;

            for (int b = a+2; b < unocc_orb; b++)
            {
                B = b/2;

                if ((A+B)%size == K)
                {

                    for (int i = Where_To_Start_Part2(rank,a-1)+1; i < n_Electrons; i+=jump)
                    {
                       // if ((AA+i/2)%size == rank)
                       // {
                            for (int j = 1; j < n_Electrons; j++)
                            {
                                MY_OWN_MPI[index_counter] = accu(t2.at(b,j) % W_4.at(a,i));
                                index_counter += 1;
                                j++;
                            }
                       // }
                       // i++;
                    }
                }

                b++;
            }
            a++;
        }

        sum_a_n = 0;
        // odd-even - (odd-even / even-odd)
        for (int a = 1; a < unocc_orb; a++)
        {
            sum_a_n += a/2;
            AA = a/2 * Speed_Elec;
            A = a/2 * Speed_Occ - sum_a_n;

            for (int b = a+1; b < unocc_orb; b++)
            {
                B = b/2;

                if ((A+B)%size == K)
                {

                    for (int i = Where_To_Start_Part2(rank,a-1)+1; i < n_Electrons; i+=jump)
                    {
                       // if ((AA+i/2)%size == rank)
                       // {
                            for (int j = 0; j < n_Electrons; j++)
                            {
                                MY_OWN_MPI[index_counter] = accu(t2.at(b,j) % W_4.at(a,i));
                                index_counter += 1;
                                j++;
                            }
                       // }
                        //i++;
                    }

                    for (int i = Where_To_Start_Part2(rank,a-1); i < n_Electrons; i+=jump)
                    {
                       // if ((AA+i/2)%size == rank)
                      //  {
                            for (int j = 1; j < n_Electrons; j++)
                            {
                                MY_OWN_MPI[index_counter] = accu(t2.at(b,j)(span(Speed_Occ, Speed_Occ+Speed_Occ-1), span()) % W_4.at(a,i)(span(Speed_Occ, Speed_Occ+Speed_Occ-1), span()));
                                index_counter += 1;
                                j++;
                            }
                      //  }
                       // i++;
                    }
                }
                b++;
            }
            a++;
        }

        // Terms where b > a, might be gridded differently than previous b < a but we still need it

        // even-even-even-even
        for (int a = 0; a < unocc_orb; a++)
        {
            A = a/2;
            AA = a/2 * Speed_Elec;

            for (int b = 0; b < a; b++)
            {
                B = b/2 * Speed_Occ - Calc_sum_a_n(b/2);
                if ((A+B)%size == K)
                {
                    for (int i = Where_To_Start_Part2(rank,a); i < n_Electrons; i+=jump)
                    {
                       // if ((AA+i/2)%size == rank)
                       // {
                            for (int j = 0; j < n_Electrons; j++)
                            {
                                MY_OWN_MPI[index_counter] = accu(t2(b,j) % W_4(a,i));
                                index_counter += 1;
                                j++;
                            }
                       // }
                       // i++;
                    }
                }
                b++;
            }
            a++;
        }

        // even-odd -(even-odd / odd-even)
        for (int a = 0; a < unocc_orb; a++)
        {
            A = a/2;
            AA = a/2 * Speed_Elec;

            for (int b = 1; b < a; b++)
            {
                B = b/2 * Speed_Occ - Calc_sum_a_n(b/2);

                if ((A+B)%size == K)
                {
                    for (int i = Where_To_Start_Part2(rank,a); i < n_Electrons; i+=jump)
                    {
                      //  if ((AA+i/2)%size == rank)
                      //  {
                            for (int j = 1; j < n_Electrons; j++)
                            {
                                MY_OWN_MPI[index_counter] = accu(t2.at(b,j) % W_4.at(a,i));
                                index_counter += 1;
                                j++;
                            }
                      // }
                      //  i++;
                    }

                    for (int i = Where_To_Start_Part2(rank,a)+1; i < n_Electrons; i+=jump)
                    {
                        //if ((AA+i/2)%size == rank)
                        //{
                            for (int j = 0; j < n_Electrons; j++)
                            {
                                MY_OWN_MPI[index_counter] = accu(t2.at(b,j)(span(0, Speed_Occ-1), span()) % W_4.at(a,i)(span(0, Speed_Occ-1), span()));
                                index_counter += 1;
                                j++;
                            }
                        //}
                        //i++;
                    }
                }
                b++;
            }
            a++;
        }


        // odd-odd-odd-odd
        for (int a = 1; a < unocc_orb; a++)
        {
            A = a/2;
            AA = a/2 * Speed_Elec;

            for (int b = 1; b < a; b++)
            {
                B = b/2 * Speed_Occ - Calc_sum_a_n(b/2);

                if ((A+B)%size == K)
                {

                    for (int i = Where_To_Start_Part2(rank,a-1)+1; i < n_Electrons; i+=jump)
                    {
                       // if ((AA+i/2)%size == rank)
                       // {
                            for (int j = 1; j < n_Electrons; j++)
                            {
                                MY_OWN_MPI[index_counter] = accu(t2.at(b,j) % W_4.at(a,i));
                                index_counter += 1;
                                j++;
                            }
                      //  }
                       // i++;
                    }
                }
                b++;
            }
            a++;
        }

        // odd-even - (odd-even / even-odd)
        for (int a = 1; a < unocc_orb; a++)
        {
            A = a/2;
            AA = a/2 * Speed_Elec;

            for (int b = 0; b < a; b++)
            {
                B = b/2 * Speed_Occ - Calc_sum_a_n(b/2);

                if ((A+B)%size == K)
                {

                    for (int i = Where_To_Start_Part2(rank,a-1)+1; i < n_Electrons; i+=jump)
                    {
                       // if ((AA+i/2)%size == rank)
                       // {
                            for (int j = 0; j < n_Electrons; j++)
                            {
                                MY_OWN_MPI[index_counter] = accu(t2.at(b,j) % W_4.at(a,i));
                                index_counter += 1;
                                j++;
                            }
                        //}
                        //i++;
                    }

                    for (int i = Where_To_Start_Part2(rank,a-1); i < n_Electrons; i+=jump)
                    {
                      //  if ((AA+i/2)%size == rank)
                      //  {
                            for (int j = 1; j < n_Electrons; j++)
                            {
                                MY_OWN_MPI[index_counter] = accu(t2.at(b,j)(span(Speed_Occ, Speed_Occ+Speed_Occ-1), span()) % W_4.at(a,i)(span(Speed_Occ, Speed_Occ+Speed_Occ-1), span()));
                                index_counter += 1;
                                j++;
                            }
                      //  }
                      //  i++;
                    }
                }
                b++;
            }
            a++;
        }
    }

    // Zero out the terms in W_5 local at each node
    for (int i = 0; i < unocc_orb; i++)
    {
        for (int j = 0; j < unocc_orb; j++)
        {
            if ((i/2 * Speed_Occ - Calc_sum_a_n(i/2) + j/2)%size == rank)
            {
                W_5(i,j) = zeros(n_Electrons, Speed_Elec);
            }
        }
    }

    // Scatter with some all to all personalised communication, Displacement must be given in number of bits, not number of doubles
    MPI_Alltoallw(MY_OWN_MPI, Global_Worksize_1[rank], Global_Displacement_1[rank],             mpi_types_array,
                  SHARED_INFO_MPI, Global_Worksize_1_1[rank], Global_Displacement_1_1[rank],    mpi_types_array,
                  MPI_COMM_WORLD);

    index_counter = 0;

    // Scatter this stuff, each node will need its own a and b parts
    for (int J = 0; J < size ; J++) // Since work is ID'ed 0, 1, 2, 3, 4, etc and distributed very good in the case where there are more cores than jobs we can use "size" and modify this value to ensure those nodes with no calculations dont call this function to avoid bugs
    {
        /*
        // Can use scatter here now and remove for K loop
        MPI_Scatterv(MY_OWN_MPI, Global_Worksize_1[J], Global_Displacement_1[J], MPI_DOUBLE, SHARED_INFO_MPI, Global_Worksize_1[J][rank], MPI_DOUBLE, J, MPI_COMM_WORLD);
        index_counter = 0;
        */

        sum_a_n = 0;
        // even-even-even-even
        for (int a = 0; a < unocc_orb; a++)
        {
            sum_a_n += a/2;
            A = a/2 * Speed_Occ - sum_a_n;
            AA = a/2 * Speed_Elec;

            for (int b = a+2; b < unocc_orb; b++)
            {
                B = b/2;
                if ((A+B)%size == rank)
                {
                    for (int i = Where_To_Start_Part2(J,a); i < n_Electrons; i+=jump)
                    {
                      //  if ((AA+i/2)%size == J)
                      //  {
                            for (int j = 0; j < n_Electrons; j++)
                            {
                                W_5(a,b)(i/2,j/2) += SHARED_INFO_MPI[index_counter];
                                index_counter += 1;
                                j++;
                            }
                       // }
                       // i++;
                    }
                }
                b++;
            }
            a++;
        }

        sum_a_n = 0;
        //even-odd -(even-odd / odd-even)
        for (int a = 0; a < unocc_orb; a++)
        {
            sum_a_n += a/2;
            A = a/2 * Speed_Occ - sum_a_n;
            AA = a/2 * Speed_Elec;
            for (int b = a+1; b < unocc_orb; b++)
            {
                B = b/2;

                if ((A+B)%size == rank)
                {
                    for (int i = Where_To_Start_Part2(J,a); i < n_Electrons; i+=jump)
                    {
                      //  if ((AA+i/2)%size == J)
                      //  {
                            for (int j = 1; j < n_Electrons; j++)
                            {
                                W_5(a,b)(i/2,j/2) += SHARED_INFO_MPI[index_counter];
                                index_counter += 1;
                                j++;
                            }
                      // }
                      //  i++;
                    }

                    for (int i = Where_To_Start_Part2(J,a)+1; i < n_Electrons; i+=jump)
                    {
                      //  if ((AA+i/2)%size == J)
                      //  {
                            for (int j = 0; j < n_Electrons; j++)
                            {
                                W_5(a,b)(i/2+Speed_Elec,j/2) += SHARED_INFO_MPI[index_counter];
                                index_counter += 1;
                                j++;
                            }
                      //  }
                       // i++;
                    }
                }
                b++;
            }
            a++;
        }

        sum_a_n = 0;
        // odd-odd-odd-odd
        for (int a = 1; a < unocc_orb; a++)
        {
            sum_a_n += a/2;
            A = a/2 * Speed_Occ - sum_a_n;
            AA = a/2 * Speed_Elec;
            for (int b = a+2; b < unocc_orb; b++)
            {
                B = b/2;

                if ((A+B)%size == rank)
                {
                    for (int i = Where_To_Start_Part2(J,a-1)+1; i < n_Electrons; i+=jump)
                    {
                       // if ((AA+i/2)%size == J)
                       // {
                            for (int j = 1; j < n_Electrons; j++)
                            {
                                W_5(a,b)(i/2+Speed_Elec,j/2) += SHARED_INFO_MPI[index_counter];
                                index_counter += 1;
                                j++;
                            }
                      //  }
                       // i++;
                    }
                }
                b++;
            }
            a++;
        }

        sum_a_n = 0;
        // odd-even - (odd-even / even-odd)
        for (int a = 1; a < unocc_orb; a++)
        {
            sum_a_n += a/2;
            A = a/2 * Speed_Occ - sum_a_n;
            AA = a/2 * Speed_Elec;
            for (int b = a+1; b < unocc_orb; b++)
            {
                B = b/2;

                if ((A+B)%size == rank)
                {

                    for (int i = Where_To_Start_Part2(J,a-1)+1; i < n_Electrons; i+=jump)
                    {
                       // if ((AA+i/2)%size == J)
                       // {
                            for (int j = 0; j < n_Electrons; j++)
                            {
                                W_5(a,b)(i/2+Speed_Elec,j/2) += SHARED_INFO_MPI[index_counter];
                                index_counter += 1;
                                j++;
                            }
                       // }
                        //i++;
                    }

                    for (int i = Where_To_Start_Part2(J,a-1); i < n_Electrons; i+=jump)
                    {
                       // if ((AA+i/2)%size == J)
                      //  {
                            for (int j = 1; j < n_Electrons; j++)
                            {
                                W_5(a,b)(i/2,j/2) += SHARED_INFO_MPI[index_counter];
                                index_counter += 1;
                                j++;
                            }
                       // }
                       // i++;
                    }
                }
                b++;
            }
            a++;
        }

        // Terms where b > a, might be gridded differently

        // even-even-even-even
        for (int a = 0; a < unocc_orb; a++)
        {
            A = a/2;
            AA = a/2 * Speed_Elec;

            for (int b = 0; b < a; b++)
            {
                B = b/2 * Speed_Occ - Calc_sum_a_n(b/2);
                if ((A+B)%size == rank)
                {
                    for (int i = Where_To_Start_Part2(J,a); i < n_Electrons; i+=jump)
                    {
                       // if ((AA+i/2)%size == J)
                       // {
                            for (int j = 0; j < n_Electrons; j++)
                            {
                                W_5(b,a)(i/2,j/2) -= SHARED_INFO_MPI[index_counter];
                                index_counter += 1;
                                j++;
                            }
                       // }
                       // i++;
                    }
                }
                b++;
            }
            a++;
        }

        // even-odd -(even-odd / odd-even)
        for (int a = 0; a < unocc_orb; a++)
        {
            A = a/2;
            AA = a/2 * Speed_Elec;

            for (int b = 1; b < a; b++)
            {
                B = b/2 * Speed_Occ - Calc_sum_a_n(b/2);

                if ((A+B)%size == rank)
                {
                    for (int i = Where_To_Start_Part2(J,a); i < n_Electrons; i+=jump)
                    {
                       // if ((AA+i/2)%size == J)
                       // {
                            for (int j = 1; j < n_Electrons; j++)
                            {
                                W_5(b,a)(i/2,j/2) -= SHARED_INFO_MPI[index_counter];
                                index_counter += 1;
                                j++;
                            }
                        //}
                       // i++;
                    }

                    for (int i = Where_To_Start_Part2(J,a)+1; i < n_Electrons; i+=jump)
                    {
                      //  if ((AA+i/2)%size == J)
                      //  {
                            for (int j = 0; j < n_Electrons; j++)
                            {
                                W_5(b,a)(i/2+Speed_Elec,j/2) -= SHARED_INFO_MPI[index_counter];
                                index_counter += 1;
                                j++;
                            }
                      //  }
                       // i++;
                    }
                }
                b++;
            }
            a++;
        }

        // odd-odd-odd-odd
        for (int a = 1; a < unocc_orb; a++)
        {
            AA = a/2 * Speed_Elec;
            A = a/2;

            for (int b = 1; b < a; b++)
            {
                B = b/2 * Speed_Occ - Calc_sum_a_n(b/2);

                if ((A+B)%size == rank)
                {

                    for (int i = Where_To_Start_Part2(J,a-1)+1; i < n_Electrons; i+=jump)
                    {
                        //if ((AA+i/2)%size == J)
                        //{
                            for (int j = 1; j < n_Electrons; j++)
                            {
                                W_5(b,a)(i/2+Speed_Elec,j/2) -= SHARED_INFO_MPI[index_counter];
                                index_counter += 1;
                                j++;
                            }
                        //}
                        //i++;
                    }
                }

                b++;
            }
            a++;
        }

        // odd-even - (odd-even / even-odd)
        for (int a = 1; a < unocc_orb; a++)
        {
            A = a/2;
            AA = a/2 * Speed_Elec;

            for (int b = 0; b < a; b++)
            {
                B = b/2 * Speed_Occ - Calc_sum_a_n(b/2);

                if ((A+B)%size == rank)
                {

                    for (int i = Where_To_Start_Part2(J,a-1)+1; i < n_Electrons; i+=jump)
                    {
                       // if ((AA+i/2)%size == J)
                       // {
                            for (int j = 0; j < n_Electrons; j++)
                            {
                                W_5(b,a)(i/2+Speed_Elec,j/2) -= SHARED_INFO_MPI[index_counter];
                                index_counter += 1;
                                j++;
                            }
                       // }
                      //  i++;
                    }

                    for (int i = Where_To_Start_Part2(J,a-1); i < n_Electrons; i+=jump)
                    {
                       // if ((AA+i/2)%size == J)
                       // {
                            for (int j = 1; j < n_Electrons; j++)
                            {
                                W_5(b,a)(i/2,j/2) -= SHARED_INFO_MPI[index_counter];
                                index_counter += 1;
                                j++;
                            }
                       // }
                      //  i++;
                    }
                }
                b++;
            }
            a++;
        }
    }
}

// Big program
// Please see text if you want to optimize further.
//
