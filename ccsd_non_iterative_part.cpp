#include "ccsd_non_iterative_part.h"

ccsd_non_iterative_part::ccsd_non_iterative_part(int unc, int elec, int ran, int siz)
{
    unocc_orb = unc;
    n_Electrons = elec;
    rank = ran;
    size = siz;

}

int ccsd_non_iterative_part::Calc_sum_a_n(int a)
{
    int sum = 0;
    for (int i = 0; i <= a; i++)
    {
        sum += i;
    }
    return sum;
}

int ccsd_non_iterative_part::ret_max_terms()
{
    return max_terms;
}

void ccsd_non_iterative_part::Map_T2_For_MPI()
{
    // This function is designed to map out T2 new in a way such that the parallel implementation is effective
    // Also to ensure we do not store symmetric terms and pass symmetric terms through parallel

    if(rank==0)
    cout << endl << endl;
    int number_counter;
    int A,B;
    int INDEX_CHECK;

    WORK_EACH_NODE = zeros(size);

    // We are gridding here over A/2 and B/2, which is unocc x unocc large.
    // Lets say we have 4 distinct A/2 values, meaning A and B goes from 0 to 7 ( => A/2 and B/2 goes from 0 to 3)
    // However we will access this symetricly, for a 4x4 example we would have:

    /*
    a = 0, 1, 2, 3

    b = 1,
    b = 2, 2
    b = 3, 3, 3

    Meaning the a grid over a and b would be marked like this, X = calculations will be done, O = no calculations to be done

    (a,b)

    (0,0) (1,0) (2,0) (3,0)
    (0,1) (1,1) ...
    ....

    X X X X
    O X X X
    O O X X
    O O O X

    This is best achieved from the distribution where X is calculated like this:

    X = A * U + B - sum(a_n)

    where sum(a_n) is calculated like this:
    A = 0 : sum(a_n) = 0
    A = 1 : sum(a_n) = 1 + 0
    A = 2 : sum(a_n) = 2 + 1 + 0
    A = 3 : sum(a_n) = 3 + 2 + 1 + 0
    A = 4 : sum(a_n) = 4 + 3 + 2 + 1 + 0


    If we set the value of X to be : X = B + A * unocc_orb + A this is achieved



    */

    int sum_a_n;

    for (int K = 0; K < size; K++)
    {
        number_counter = 0;

        sum_a_n = 0;

        for (int a = 0; a < unocc_orb; a++)
        {
            sum_a_n += a/2;

            A = a/2 * unocc_orb/2 - sum_a_n;

            for (int b = a+2; b < unocc_orb; b++) // b even
            {
                B = b/2;
                INDEX_CHECK = A+B;

                if (INDEX_CHECK % size == K)
                {
                    for (int i = 0; i < n_Electrons; i++) // i even
                    {
                        for (int j = i+2; j < n_Electrons; j++) // j even
                        {
                            number_counter += 1;
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
                INDEX_CHECK = A+B;

                if (INDEX_CHECK % size == K)
                {
                    for (int i = 0; i < n_Electrons; i++) // i even
                    {
                        for (int j = i+1; j < n_Electrons; j++) // j odd
                        {
                            number_counter += 1;
                            j++;
                        }
                        i++;
                    }

                    for (int i = 1; i < n_Electrons; i++) // i odd
                    {
                        for (int j = i+1; j < n_Electrons; j++) // j even
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
            // At this point a is an odd number, but a/2 is the same value as before so we can use the same single bar integrals potentially

            for (int b = a+2; b < unocc_orb; b++) // b odd
            {
                B = b/2;
                INDEX_CHECK = A+B;

                if (INDEX_CHECK % size == K)
                {
                    for (int i = 1; i < n_Electrons; i++) // i odd
                    {
                        for (int j = i+2; j < n_Electrons; j++) // j odd
                        {
                            number_counter += 1;
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
                INDEX_CHECK = A+B;

                if (INDEX_CHECK % size == K)
                {
                    for (int i = 0; i < n_Electrons; i++) // i even
                    {
                        for (int j = i+1; j < n_Electrons; j++) // j odd
                        {
                            number_counter += 1;
                            j++;
                        }
                        i++;
                    }

                    for (int i = 1; i < n_Electrons; i++) // i odd
                    {
                        for (int j = i+1; j < n_Electrons; j++) // j even
                        {
                            number_counter += 1;
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
            for (int i = 0; i < n_Electrons; i++)
            {
                INDEX_CHECK = a/2 * n_Electrons/2+ i/2;
                if (INDEX_CHECK % size == K)
                {
                    number_counter += 1;
                }

                i++;
            }
            a++;
        }
/*
        for (int a = 1; a < unocc_orb; a++)
        {
            for (int i = 1; i < n_Electrons; i++)
            {
                INDEX_CHECK = a/2 * n_Electrons/2+ i/2;
                if (INDEX_CHECK % size == K)
                {
                    number_counter += 1;
                }
                i++;
            }
            a++;
        }
*/



        if (rank==0)
        {
            cout << "Work t2 : " << number_counter << " for node " << K << endl;
        }

        WORK_EACH_NODE(K) = number_counter;

    }



    if(rank==0)
    cout << endl << endl;
}

vec ccsd_non_iterative_part::return_Work_T2()
{
    Map_T2_For_MPI(); // first map
    return WORK_EACH_NODE; // Then return
}

void ccsd_non_iterative_part::Map_Part1_For_MPI()
{
    // Mapping of how much work goes where for part 1 of MPI
    // Mapping is done for communication minimization and this will only be called way before iterations start
    // How much work per node is actually calculated before any work is done, and will in future be implemented in a script

    // ------------------------------------------------------------------------------------------------------------------------
    // This way you can know how your scaling will be PRIOR to initiating a supercomputer run which will take meny hours
    // This will be awsome
    // ------------------------------------------------------------------------------------------------------------------------

    // This function will need some aditional terms, for initialization of W_1 and W_2 :O

    WORK_EACH_NODE_Part1 = zeros(size);
    WORK_EACH_NODE_Part2 = zeros(size);

    int index_counter;
    int index_counter2;
    int INDEX_CHECK;
    int A,I, B;

    int** Displacement_W_4;
    int** Nr_Elems_W_4;
    int** Displacement_W_4_2;
    int** Nr_Elems_W_4_2;

    // Use this in W_4 only
    Displacement_W_4 = (int**)malloc(size*sizeof(int*));
    Nr_Elems_W_4 = (int**)malloc(size*sizeof(int*));
    Global_Displacement_1 = (int**)malloc(size*sizeof(int*));
    Global_Displacement_2 = (int**)malloc(size*sizeof(int*));
    Global_Worksize_1 = (int**)malloc(size*sizeof(int*));
    Global_Worksize_2 = (int**)malloc(size*sizeof(int*));

    for (int i = 0; i < size; i++)
    {
        Displacement_W_4[i] = (int*)malloc(size * sizeof(int));
        Nr_Elems_W_4[i] = (int*)malloc(size*sizeof(int));

        for (int j = 0; j < size; j++)
        {
            Displacement_W_4[i][j] = 0;
            Nr_Elems_W_4[i][j] = 0;
        }
    }

    Displacement_W_4_2 = (int**)malloc(size*sizeof(int*));
    Nr_Elems_W_4_2 = (int**)malloc(size*sizeof(int*));


    for (int i = 0; i < size; i++)
    {
        Displacement_W_4_2[i] = (int*)malloc(size * sizeof(int));
        Nr_Elems_W_4_2[i] = (int*)malloc(size*sizeof(int));

        for (int j = 0; j < size; j++)
        {
            Displacement_W_4_2[i][j] = 0;
            Nr_Elems_W_4_2[i][j] = 0;
        }
    }

    int AA;

    for (int J = 0; J < size; J++)
    {
        for (int K = 0; K < size; K++)
        {
            // Arrange this so one node gets its information in sequence
            for (int a = 0; a < unocc_orb; a++)
            {
                AA = a/2*n_Electrons/2;
                A = a/2;
                for(int i = 0; i < n_Electrons; i++)
                {
                    I = i/2;
                    if ((AA+I)%size == K)
                    {
                        for (int m = 0; m < n_Electrons; m++)
                        {
                            INDEX_CHECK = AA + m/2;
                            if (INDEX_CHECK % size == J)
                            {
                                for (int e = 0; e < unocc_orb; e++)
                                {
                                    Nr_Elems_W_4_2[J][K] += 1;
                                    e++;
                                }
                            }
                            m++;
                        }

                        for (int m = 1; m < n_Electrons; m++)
                        {
                            INDEX_CHECK = AA + m/2;
                            if (INDEX_CHECK % size == J)
                            {
                                for (int e = 1; e < unocc_orb; e++)
                                {
                                    Nr_Elems_W_4_2[J][K] += 1;
                                    e++;
                                }
                            }
                            m++;
                        }
                    }
                    i++;
                }
                a++;
            }

            for (int a = 0; a < unocc_orb; a++)
            {
                AA = a/2*n_Electrons/2;
                A = a/2;
                for (int i = 1; i < n_Electrons; i++)
                {
                    I = i/2;
                    if ((AA+I)%size == K)
                    {
                        for (int m = 1; m < n_Electrons; m++)
                        {
                            INDEX_CHECK = AA + m/2;
                            if (INDEX_CHECK % size == J)
                            {
                                for (int e = 0; e < unocc_orb; e++)
                                {
                                    Nr_Elems_W_4_2[J][K] += 1;
                                    e++;
                                }
                            }
                            m++;
                        }
                    }
                    i++;
                }
                a++;
            }

            for (int a = 1; a < unocc_orb; a++)
            {
                A = a/2;
                AA = a/2*n_Electrons/2;
                for (int i = 0; i < n_Electrons; i++)
                {
                    I = i/2;
                    if ((AA+I)%size == K)
                    {
                        for (int m = 0; m < n_Electrons; m++)
                        {
                            INDEX_CHECK = AA + m/2;
                            if (INDEX_CHECK % size == J)
                            {
                                for (int e = 1; e < unocc_orb; e++)
                                {
                                    Nr_Elems_W_4_2[J][K] += 1;
                                    e++;
                                }
                            }
                            m++;
                        }
                    }
                    i++;
                }
                a++;
            }

            for (int a = 1; a < unocc_orb; a++)
            {
                A = a/2;
                AA = a/2*n_Electrons/2;
                for (int i = 1; i < n_Electrons; i++)
                {
                    I = i/2;
                    if ((AA+I)%size == K)
                    {
                        for (int m = 0; m < n_Electrons; m++)
                        {
                            INDEX_CHECK = AA + m/2;
                            if (INDEX_CHECK % size == J)
                            {
                                for (int e = 0; e < unocc_orb; e++)
                                {
                                    Nr_Elems_W_4_2[J][K] += 1;
                                    e++;
                                }
                            }
                            m++;
                        }

                        for (int m = 1; m < n_Electrons; m++)
                        {
                            INDEX_CHECK = AA + m/2;
                            if (INDEX_CHECK % size == J)
                            {
                                for (int e = 1; e < unocc_orb; e++)
                                {
                                    Nr_Elems_W_4_2[J][K] += 1;
                                    e++;
                                }
                            }
                            m++;
                        }
                    }
                    i++;
                }
                a++;
            }
        }


    }







    int sum_a_n;
    int Speed_Occ = unocc_orb/2;

    for (int J = 0; J < size; J++)
    {
    for (int K = 0; K < size; K++)
    {
        sum_a_n = 0;
        for (int a = 0; a < unocc_orb; a++)
        {
            sum_a_n += a/2;
            A = a/2 * unocc_orb/2 - sum_a_n;
            AA = a/2 * n_Electrons/2;
            for (int b = a+2; b < unocc_orb; b++)
            {
                B = b/2;
                if ((A+B)%size == K)
                {
                    for (int i = 0; i < n_Electrons; i++)
                    {
                        if ((AA+i/2)%size == J)
                        {
                            for (int j = 0; j < n_Electrons; j++)
                            {
                                Nr_Elems_W_4[J][K]+= 1;
                                j++;
                            }
                        }
                        i++;
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
            A = a/2 * unocc_orb/2 - sum_a_n;
            AA = a/2 * n_Electrons/2;

            for (int b = a+1; b < unocc_orb; b++)
            {
                B = b/2;

                if ((A+B)%size == K)
                {
                    for (int i = 0; i < n_Electrons; i++)
                    {
                        if ((AA+i/2)%size == J)
                        {
                            for (int j = 1; j < n_Electrons; j++)
                            {
                                Nr_Elems_W_4[J][K]+= 1;
                                j++;
                            }
                        }
                        i++;
                    }

                    for (int i = 1; i < n_Electrons; i++)
                    {
                        if ((AA+i/2)%size == J)
                        {
                            for (int j = 0; j < n_Electrons; j++)
                            {
                                Nr_Elems_W_4[J][K]+= 1;
                                j++;
                            }
                        }
                        i++;
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
            A = a/2 * unocc_orb/2 - sum_a_n;
            AA = a/2 * n_Electrons/2;

            for (int b = a+2; b < unocc_orb; b++)
            {
                B = b/2;

                if ((A+B)%size == K)
                {

                    for (int i = 1; i < n_Electrons; i++)
                    {
                        if ((AA+i/2)%size == J)
                        {
                            for (int j = 1; j < n_Electrons; j++)
                            {
                                Nr_Elems_W_4[J][K]+= 1;
                                j++;
                            }
                        }
                        i++;
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
            A = a/2 * unocc_orb/2 - sum_a_n;
            AA = a/2 * n_Electrons/2;

            for (int b = a+1; b < unocc_orb; b++)
            {
                B = b/2;

                if ((A+B)%size == K)
                {

                    for (int i = 1; i < n_Electrons; i++)
                    {
                        if ((AA+i/2)%size == J)
                        {
                            for (int j = 0; j < n_Electrons; j++)
                            {
                                Nr_Elems_W_4[J][K]+= 1;
                                j++;
                            }
                        }
                        i++;
                    }

                    for (int i = 0; i < n_Electrons; i++)
                    {
                        if ((AA+i/2)%size == J)
                        {
                            for (int j = 1; j < n_Electrons; j++)
                            {
                                Nr_Elems_W_4[J][K]+= 1;
                                j++;
                            }
                        }
                        i++;
                    }
                }
                b++;
            }
            a++;
        }










        // even-even-even-even
        for (int a = 0; a < unocc_orb; a++)
        {
            A = a/2;
            AA = a/2 * n_Electrons/2;

            for (int b = 0; b < a; b++)
            {
                B = b/2 * Speed_Occ - Calc_sum_a_n(b/2);
                if ((A+B)%size == K)
                {
                    for (int i = 0; i < n_Electrons; i++)
                    {
                        if ((AA+i/2)%size == J)
                        {
                            for (int j = 0; j < n_Electrons; j++)
                            {
                                Nr_Elems_W_4[J][K]+= 1;
                                j++;
                            }
                        }
                        i++;
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
            AA = a/2 * n_Electrons/2;

            for (int b = 1; b < a; b++)
            {
                B = b/2 * Speed_Occ - Calc_sum_a_n(b/2);

                if ((A+B)%size == K)
                {
                    for (int i = 0; i < n_Electrons; i++)
                    {
                        if ((AA+i/2)%size == J)
                        {
                            for (int j = 1; j < n_Electrons; j++)
                            {
                                Nr_Elems_W_4[J][K]+= 1;
                                j++;
                            }
                        }
                        i++;
                    }

                    for (int i = 1; i < n_Electrons; i++)
                    {
                        if ((AA+i/2)%size == J)
                        {
                            for (int j = 0; j < n_Electrons; j++)
                            {
                                Nr_Elems_W_4[J][K]+= 1;
                                j++;
                            }
                        }
                        i++;
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
            AA = a/2 * n_Electrons/2;

            for (int b = 1; b < a; b++)
            {
                B = b/2 * Speed_Occ - Calc_sum_a_n(b/2);

                if ((A+B)%size == K)
                {

                    for (int i = 1; i < n_Electrons; i++)
                    {
                        if ((AA+i/2)%size == J)
                        {
                            for (int j = 1; j < n_Electrons; j++)
                            {
                                Nr_Elems_W_4[J][K]+= 1;
                                j++;
                            }
                        }
                        i++;
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
            AA = a/2 * n_Electrons/2;

            for (int b = 0; b < a; b++)
            {
                B = b/2 * Speed_Occ - Calc_sum_a_n(b/2);

                if ((A+B)%size == K)
                {

                    for (int i = 1; i < n_Electrons; i++)
                    {
                        if ((AA+i/2)%size == J)
                        {
                            for (int j = 0; j < n_Electrons; j++)
                            {
                                Nr_Elems_W_4[J][K]+= 1;
                                j++;
                            }
                        }
                        i++;
                    }

                    for (int i = 0; i < n_Electrons; i++)
                    {
                        if ((AA+i/2)%size == J)
                        {
                            for (int j = 1; j < n_Electrons; j++)
                            {
                                Nr_Elems_W_4[J][K]+= 1;
                                j++;
                            }
                        }
                        i++;
                    }
                }
                b++;
            }
            a++;
        }





    }
    }

    // Now we know how to communication is going to work
    // To ensure we do not get bugs we must ensure that only the nodes that are going to recieve any information is a part of the communication group
    // This means we cannot use MPI_COMM_WORLD

    for (int K = 0; K < size; K++)
    {
        Displacement_W_4[K][0] = 0;
        Displacement_W_4_2[K][0] = 0;
        for (int J = 1; J < size; J++)
        {
            Displacement_W_4[K][J] = Nr_Elems_W_4[K][J-1] + Displacement_W_4[K][J-1];
            Displacement_W_4_2[K][J] = Nr_Elems_W_4_2[K][J-1] + Displacement_W_4_2[K][J-1];
        }
    }

    // Figure our who is going to be a part of this communication, if someone gets 0 bytes then it will cause bug
   // MPI_WORLD_W4_GROUPS_1 = (MPI_Comm*)malloc(size*sizeof(MPI_Comm));
   // MPI_WORLD_W4_GROUPS_2 = (MPI_Comm*)malloc(size*sizeof(MPI_Comm));

    max_terms = 0;

    for (int K = 0; K < size; K++)
    {
        index_counter = 0;
        index_counter2 = 0;

        for (int J = 0; J < size; J++)
        {
            if (Nr_Elems_W_4[K][J]!= 0)
            {
                index_counter += 1;
            }

            if (Nr_Elems_W_4_2[K][J] != 0)
            {
                index_counter2 += 1;
            }
        }

        // Figure out what goes into this communication group
       //int* local_group = (int*)malloc(index_counter*sizeof(int));
        //int* local_group2 = (int*)malloc(index_counter2*sizeof(int));

        // Define  displacement and amount of doubles for communication later
        Global_Displacement_1[K] = (int*)malloc(size*sizeof(int));
        Global_Displacement_2[K] = (int*)malloc(size*sizeof(int));
        Global_Worksize_1[K] = (int*)malloc(size*sizeof(int));
        Global_Worksize_2[K] = (int*)malloc(size*sizeof(int));

        // Copy the nubmers we need
        for (int J = 0; J < size; J++)
        {
                Global_Worksize_1[K][J] = Nr_Elems_W_4[K][J];

                if (Global_Worksize_1[K][J] > max_terms)
                {
                    max_terms = Global_Worksize_1[K][J];
                }

                Global_Worksize_2[K][J] = Nr_Elems_W_4_2[K][J];

                if (Global_Worksize_2[K][J] > max_terms)
                {
                    max_terms = Global_Worksize_1[K][J];
                }
        }

        // Copy the displacements we need
        Global_Displacement_1[K][0] = 0;
        for (int J = 1; J < size; J++)
        {
            Global_Displacement_1[K][J] = Global_Displacement_1[K][J-1] + Global_Worksize_1[K][J-1];
        }

        Global_Displacement_2[K][0] = 0;
        for (int J = 1; J < size; J++)
        {
            Global_Displacement_2[K][J] = Global_Displacement_2[K][J-1] + Global_Worksize_2[K][J-1];
        }
/*
        index_counter = 0;
        index_counter2 = 0;

        // Define the ranks we need
        for (int J = 0; J < size; J++)
        {
            if (Nr_Elems_W_4[K][J]!= 0)
            {
                local_group[index_counter] = J;
                index_counter += 1;
            }

            if (Nr_Elems_W_4_2[K][J] != 0)
            {
                local_group2[index_counter2] = J;
                index_counter2 += 1;
            }
        }

        MPI_Group list_of_ranks, list_of_ranks2;
        MPI_Group_incl(MPI_COMM_WORLD, index_counter, local_group, list_of_ranks);
        MPI_Group_incl(MPI_COMM_WORLD, index_counter, local_group2, list_of_ranks2);

        // Define the new communicator group
        MPI_Comm_create(MPI_COMM_WORLD, list_of_ranks, &MPI_WORLD_W4_GROUPS_1[K]);
        MPI_Comm_create(MPI_COMM_WORLD, list_of_ranks2, &MPI_WORLD_W4_GROUPS_2[K]);

        // Remove the communication group we have already defined
        free(local_group);
        free(local_group2);
        */
    }




    // Free some memory up
    for (int i = 0; i < size; i++)
    {
        free(Nr_Elems_W_4[i]);
        free(Nr_Elems_W_4_2[i]);
        free(Displacement_W_4[i]);
        free(Displacement_W_4_2[i]);
    }

    free(Nr_Elems_W_4);
    free(Nr_Elems_W_4_2);
    free(Displacement_W_4);
    free(Displacement_W_4_2);

    Local_Displacement_1 = unocc_orb/2 * n_Electrons/2;
    //Local_Displacement_2 = 2*Local_Displacement_1;


    int KK;
    int L, E;
    int Speed_Elec = (int) n_Electrons/2;

    for (int K = 0; K < size; K++)
    {
        index_counter = 0;





        ///////////// W1 part


        for (int k = 0; k < n_Electrons; k++)
        {
            KK = k/2*Speed_Elec;
            L = 0;
            for (int l = 1; l < k; l++)
            {
                INDEX_CHECK = KK+L+Local_Displacement_1;
                if (INDEX_CHECK % size == K)
                {
                    for (int i = 0; i < n_Electrons; i++)
                    {
                        i++;
                        index_counter += 1;
                    }

                    for (int i = 0; i < n_Electrons; i++)
                    {
                        E = 0;
                        for (int e = 1; e < unocc_orb; e++)
                        {
                            e++;
                            index_counter += 1;
                        }
                        i++;
                    }

                    for (int i = 1; i < n_Electrons; i++)
                    {
                        for (int e = 0; e < unocc_orb; e++)
                        {
                            e++;
                            index_counter += 1;
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
                L++;
            }

            for (int l = k+1; l < n_Electrons; l++)
            {
                INDEX_CHECK = KK+L+Local_Displacement_1;
                if (INDEX_CHECK % size == K)
                {

                    for (int i = 0; i < n_Electrons; i++)
                    {
                        i++;
                        index_counter += 1;
                    }
                }
                l++;
                L++;
            }
            k++;
        }



        for (int k = 1; k < n_Electrons; k++)
        {
            KK = k/2*Speed_Elec;
            L = 0;
            for (int l = 0; l < k; l++)
            {

                INDEX_CHECK = KK+L+Local_Displacement_1;
                if (INDEX_CHECK % size == K)
                {
                    for (int i = 0; i < n_Electrons; i++)
                    {
                        for (int e = 1; e < unocc_orb; e++)
                        {
                            e++;
                            index_counter += 1;
                        }
                        i++;
                    }

                    for (int i = 1; i < n_Electrons; i++)
                    {
                        for (int e = 0; e < unocc_orb; e++)
                        {
                            e++;
                            index_counter += 1;
                        }
                        i++;
                    }

                    for (int i = 0; i < n_Electrons; i++)
                    {
                        for (int j = i+1; j < n_Electrons; j++)
                        {
                            j++;
                            index_counter += 1;
                        }
                        i++;
                    }

                    for (int i = 1; i < n_Electrons; i++)
                    {
                        for (int j = i+1; j < n_Electrons; j++)
                        {
                            j++;
                            index_counter += 1;
                        }
                        i++;
                    }

                }
                L++;
                l++;
            }
            k++;
        }

        for (int k = 0; k < n_Electrons; k++)
        {
            KK = k/2 *Speed_Elec;
            L = 0;
            for (int l = 0; l < k+1; l++)
            {
                INDEX_CHECK = KK+L+Local_Displacement_1;
                if (INDEX_CHECK % size == K)
                {
                    for (int i = 0; i < n_Electrons; i++)
                    {
                        i++;
                        index_counter += 1;
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

                    I = 0;
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
                L++;
            }

            for (int l = k+2; l < n_Electrons; l++)
            {
                INDEX_CHECK = KK+L+Local_Displacement_1;
                if (INDEX_CHECK % size == K)
                {
                    for (int i = 0; i < n_Electrons; i++)
                    {
                        i++;
                        index_counter += 1;
                    }
                }

                l++;
                L++;
            }
            k++;
        }

        for (int k = 1; k < n_Electrons; k++)
        {
            KK = k/2*Speed_Elec;
            L = 0;
            for (int l = 1; l < k+1; l++)
            {
                INDEX_CHECK = KK+L+Local_Displacement_1;
                if (INDEX_CHECK % size == K)
                {
                    for (int i = 1; i < n_Electrons; i++)
                    {
                        for (int e = 1; e < unocc_orb; e++)
                        {
                            e++;
                            index_counter += 1;
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
                L++;
            }
            k++;
        }


        ///////////////////////////////////////
        ///////////////////////////////////////
        ///////////////////////////////////////
        ///////////////////////////////////////
        //////////// W2 part    ///////////////
        ///////////////////////////////////////
        ///////////////////////////////////////
        ///////////////////////////////////////
        ///////////////////////////////////////


        int A, M;

        // Optimized W_2 version
        A = 0;
        for (int a = 0; a < unocc_orb; a++)
        {
            M = 0;
            for (int m = 0; m < n_Electrons; m++)
            {
                INDEX_CHECK = a/2* Speed_Elec +m/2;
                if (INDEX_CHECK % size == K)
                {
                    for (int e = 0; e < unocc_orb; e++)
                    {
                        e++;
                        index_counter += 1;
                    }

                    I = 0;
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
                if (INDEX_CHECK % size == K)
                {

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

                if (INDEX_CHECK % size == K)
                {
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
                M++;
            }
            a++;
            A++;
        }


        ///////////////////////////////////////
        ///////////////////////////////////////
        ///////////////////////////////////////
        ///////////////////////////////////////
        //////////// F2 and F3  ///////////////
        ///////////////////////////////////////
        ///////////////////////////////////////
        ///////////////////////////////////////
        ///////////////////////////////////////


        for (int e = 0; e < unocc_orb; e++)
        {
            M = 0;
            for (int m = 0; m < n_Electrons; m++)
            {
                INDEX_CHECK = e/2 * Speed_Elec+ M;
                if (INDEX_CHECK % size == K)
                {
                    for (int i = 0; i < n_Electrons; i++)
                    {
                        index_counter += 1;
                        i++;
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
                if (INDEX_CHECK % size == K)
                {
                    for (int i = 0; i < n_Electrons; i++)
                    {
                        index_counter += 1;
                        i++;

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
                INDEX_CHECK = e/2* Speed_Elec + m/2;
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

        for (int a = 0; a < unocc_orb; a++)
        {
            for (int k = 0; k < n_Electrons; k++)
            {
                INDEX_CHECK = a/2 * Speed_Elec+ k/2;
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

        WORK_EACH_NODE_Part1(K) = index_counter;
        if (rank == 0)
        {
            cout << "Node : " << K << " work 1: " << index_counter << endl;
        }
    }

    if(rank==0)
    cout << endl << endl;

    for (int K = 0; K < size; K++)
    {








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

            for (int i = 0; i < size; i++)
            {
                index_counter2 += Global_Worksize_2[K][i];
            }


/*
            for (int a = 0; a < unocc_orb; a++)
            {
                for (int m = 0; m < n_Electrons; m++)
                {
                    INDEX_CHECK = a/2 +m/2;
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
                    INDEX_CHECK = a/2+ m/2;
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
                    INDEX_CHECK = a/2 +m/2;
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
                    INDEX_CHECK = a/2+ m/2;
                    if (INDEX_CHECK % size == K)
                    {
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


            */
            if (rank == 0)
            {
                cout << "Node : " << K << " work 2: " << index_counter2 << endl;
            }



        WORK_EACH_NODE_Part2(K) = index_counter2;
    }

    if(rank==0)
    cout << endl << endl;
}

vec ccsd_non_iterative_part::return_Work_P1()
{
    Map_Part1_For_MPI();
    return WORK_EACH_NODE_Part1;
}

vec ccsd_non_iterative_part::return_Work_P2()
{
    return WORK_EACH_NODE_Part2;
}

int **ccsd_non_iterative_part::ret_Global_Disp1()
{
    return Global_Displacement_1;
}

int **ccsd_non_iterative_part::ret_Global_Disp2()
{
    return Global_Displacement_2;
}

int **ccsd_non_iterative_part::ret_Global_Work1()
{
    return Global_Worksize_1;
}

int **ccsd_non_iterative_part::ret_Global_Work2()
{
    return Global_Worksize_2;
}

void ccsd_non_iterative_part::Figure_Out_Where_To_Start()
{
    Where_To_Start_Part1 = zeros(size, unocc_orb);
    Where_To_Start_Part2 = zeros(size, unocc_orb);
    Where_To_start_T2 = zeros(size, unocc_orb);

    int A, B, E;
    int I;
    int INDEX_CHECK;
    int speed_elec = n_Electrons/2;
    int sum_a_n;

    for (int K = 0; K < size; K++)
    {
        for (int i = 0; i < unocc_orb; i++)
        {
            I = i/2 * speed_elec;

            for (int j = 0; j < 2*size; j++)
            {
                INDEX_CHECK = I + j/2;

                if (INDEX_CHECK % size == K)
                {
                    E = j;
                    j += 2*size;
                }
                j++;
            }

            Where_To_Start_Part2(K, i) = E;
            i++;
        }
    }

    for (int K = 0; K < size; K++)
    {
        sum_a_n = 0;
        for (int i = 0; i < unocc_orb; i++)
        {
            sum_a_n += i/2;

            I = i/2 * speed_elec - sum_a_n;

            for (int j = 0; j < 2*size; j++)
            {
                INDEX_CHECK = I + j/2;

                if (INDEX_CHECK % size == K)
                {
                    E = j;
                    j += 2*size;
                }
                j++;
            }

            Where_To_start_T2(K, i) = E;
            i++;
        }
    }

}

mat ccsd_non_iterative_part::Return_Start_Part2_Pos()
{
    Figure_Out_Where_To_Start();
    return Where_To_Start_Part2;
}

mat ccsd_non_iterative_part::Return_Start_T2_Pos()
{
    return Where_To_start_T2;
}


