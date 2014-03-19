#include "coupled_cluster_integrals.h"

coupled_cluster_integrals::coupled_cluster_integrals(int mat_siz, mat ccc, int n_E)
{
    Matrix_Size = mat_siz;
    c = ccc;
    n_Electrons = n_E;
}

void coupled_cluster_integrals::Mount_Indexed_Integrals(vec intsrr)
{
    Indexed_integrals = intsrr;
}

void coupled_cluster_integrals::Return_Integrals(vec Integrals, int Version)
{
    // This is ment to be an efficient and good way to calculate the double integrals we will need

    // Transformation from Atomic Orbital (AO) to Molecule Orbital (MO)
    // Sometimes this procedure is refered to as AO2MO or AOtoMO

    vector<cube> temp;
    vector <cube> temp2;
    vector <cube> temp3;
    cube temp_cube = zeros(Matrix_Size, Matrix_Size, Matrix_Size);
    mat field_filler_matrix2;


    for (int i = 0; i < Matrix_Size; i++)
    {
        temp.push_back(temp_cube);
        temp2.push_back(temp_cube);
        temp3.push_back(temp_cube);
        temp4.push_back(temp_cube);
    }

    for (int a = 0; a < Matrix_Size; a++)
    {
        for (int i = 0; i < Matrix_Size; i++)
        {
            for (int j = 0; j < Matrix_Size; j++)
            {
                for (int k = 0; k < Matrix_Size; k++)
                {
                    for (int l = 0; l < Matrix_Size; l++)
                    {
                        temp.at(a)(j,k,l) += Integrals(Return_Integral_Index(i,j,k,l)) * c(i,a);
                    }
                }
            }
        }
    }
    Integrals.clear();

    for (int a = 0; a < Matrix_Size; a++)
    {
        for (int b = 0; b < Matrix_Size; b++)
        {
            for (int i = 0; i < Matrix_Size; i++)
            {
                for (int j = 0; j < Matrix_Size; j++)
                {
                    for (int k = 0; k < Matrix_Size; k++)
                    {
                        temp2.at(a)(b,j,k) += c(i,b) * temp.at(a)(i,j,k);
                    }
                }
            }
        }
    }
    temp.clear();

    for (int a = 0; a < Matrix_Size; a++)
    {
        for (int b = 0; b < Matrix_Size; b++)
        {
            for (int g = 0; g < Matrix_Size; g++)
            {
                for (int i = 0; i < Matrix_Size; i++)
                {
                    for (int j = 0; j < Matrix_Size; j++)
                    {
                        temp3.at(a)(b,g,j) += c(i,g) * temp2.at(a)(b,i,j);
                    }
                }
            }
        }
    }

    temp2.clear();

    for (int a = 0; a < Matrix_Size; a++)
    {
        for (int b = 0; b < Matrix_Size; b++)
        {
            for (int g = 0; g < Matrix_Size; g++)
            {
                for (int h = 0; h < Matrix_Size; h++)
                {
                    for (int i = 0; i < Matrix_Size; i++)
                    {
                        temp4.at(a)(b,g,h) += c(i,h) * temp3.at(a)(b,g,i);
                    }
                }
            }
        }
    }

    temp3.clear();

    //vector<cube> value;
    double val1, val2;

    int unocc_orb = 2*Matrix_Size - n_Electrons;
    int Speed_Occ = unocc_orb/2;

    if (Version == 2)
    {

        Rear_cube.set_size(unocc_orb, n_Electrons);
        field_filler_matrix2 = zeros(unocc_orb, n_Electrons/2);
        for (int i = 0; i < n_Electrons; i++)
        {
            for (int a = 0; a < unocc_orb; a++)
            {
                Rear_cube(a,i) = field_filler_matrix2;
            }
        }

        // Using properties of the integrals we can store only half the elements


        // Compact form of I(a,b,i,j), dont store any zeroes
        int B, J;
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
                            val1 = EqualFunc((a+n_Electrons)%2, i%2) * EqualFunc((b+n_Electrons)%2, j%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) (i/2), (int) ((b+n_Electrons)/2), (int) (j/2));

                            val2 = EqualFunc((a+n_Electrons)%2, j%2) * EqualFunc((b+n_Electrons)%2, i%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) (j/2), (int) ((b+n_Electrons)/2), (int) (i/2));

                            Rear_cube(a,i)(B,J) = val1-val2; //temp_cube3(a+n_Electrons,b+n_Electrons)(i,j);
                            j++;
                            J++;
                        }
                    }

                    else
                    {
                        J = 0;
                        for (int j = 1; j < n_Electrons; j++)
                        {
                            val1 = EqualFunc((a+n_Electrons)%2, i%2) * EqualFunc((b+n_Electrons)%2, j%2) *
                                     temp4.at((int) ((a+n_Electrons)/2))((int) (i/2), (int) ((b+n_Electrons)/2), (int) (j/2));


                            val2 = EqualFunc((a+n_Electrons)%2, j%2) * EqualFunc((b+n_Electrons)%2, i%2) *
                                     temp4.at((int) ((a+n_Electrons)/2))((int) (j/2), (int) ((b+n_Electrons)/2), (int) (i/2));

                             Rear_cube(a,i)(B,J) = val1-val2;
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
                            val1 = EqualFunc((a+n_Electrons)%2, i%2) * EqualFunc((b+n_Electrons)%2, j%2) *
                                     temp4.at((int) ((a+n_Electrons)/2))((int) (i/2), (int) ((b+n_Electrons)/2), (int) (j/2));

                            val2 = EqualFunc((a+n_Electrons)%2, j%2) * EqualFunc((b+n_Electrons)%2, i%2) *
                                     temp4.at((int) ((a+n_Electrons)/2))((int) (j/2), (int) ((b+n_Electrons)/2), (int) (i/2));

                            Rear_cube(a,i)(B+Speed_Occ,J) = val1-val2;
                            j++;
                            J++;
                        }
                    }

                    else
                    {
                        J = 0;
                        for (int j = 1; j < n_Electrons; j++)
                        {
                            val1 = EqualFunc((a+n_Electrons)%2, i%2) * EqualFunc((b+n_Electrons)%2, j%2) *
                                     temp4.at((int) ((a+n_Electrons)/2))((int) (i/2), (int) ((b+n_Electrons)/2), (int) (j/2));

                            val2 = EqualFunc((a+n_Electrons)%2, j%2) * EqualFunc((b+n_Electrons)%2, i%2) *
                                     temp4.at((int) ((a+n_Electrons)/2))((int) (j/2), (int) ((b+n_Electrons)/2), (int) (i/2));

                            Rear_cube(a,i)(B+Speed_Occ,J) = val1-val2;
                            j++;
                            J++;
                        }
                    }
                }
                b++;
                B++;
            }
        }


        // Compact form of I(a,b,c,d), dont store any zeroes
        Rear_cube2.set_size(unocc_orb, unocc_orb);
        field_filler_matrix2 = zeros(unocc_orb, Speed_Occ);
        for (int i = 0; i < unocc_orb; i++)
        {
            for (int a = 0; a < unocc_orb; a++)
            {
                Rear_cube2(a,i) = field_filler_matrix2;
            }
        }

        for (int a = 0; a < unocc_orb; a++)
        {
            B = 0;
            for (int b = 0; b < unocc_orb; b++)
            {
                for (int i = 0; i < unocc_orb; i++)
                {
                    if ((a+b+i)%2 == 0)
                    {
                        J = 0;
                        for (int j = 0; j < unocc_orb; j++)
                        {
                            val1 = EqualFunc((a+n_Electrons)%2, (b+n_Electrons)%2) * EqualFunc((i+n_Electrons)%2, (j+n_Electrons)%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) ((b+n_Electrons)/2), (int) ((i+n_Electrons)/2), (int) ((j+n_Electrons)/2));

                            val2 = EqualFunc((a+n_Electrons)%2, (j+n_Electrons)%2) * EqualFunc((i+n_Electrons)%2, (b+n_Electrons)%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) ((j+n_Electrons)/2), (int) ((i+n_Electrons)/2), (int) ((b+n_Electrons)/2));

                            Rear_cube2(a,i)(B,J) = val1 - val2;
                            j++;
                            J++;
                        }
                    }

                    else
                    {
                        J = 0;
                        for (int j = 1; j < unocc_orb; j++)
                        {
                            val1 = EqualFunc((a+n_Electrons)%2, (b+n_Electrons)%2) * EqualFunc((i+n_Electrons)%2, (j+n_Electrons)%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) ((b+n_Electrons)/2), (int) ((i+n_Electrons)/2), (int) ((j+n_Electrons)/2));

                            val2 = EqualFunc((a+n_Electrons)%2, (j+n_Electrons)%2) * EqualFunc((i+n_Electrons)%2, (b+n_Electrons)%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) ((j+n_Electrons)/2), (int) ((i+n_Electrons)/2), (int) ((b+n_Electrons)/2));

                            Rear_cube2(a,i)(B,J) = val1 - val2;
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
                for (int i = 0; i < unocc_orb; i++)
                {
                    if ((a+b+i)%2 == 0)
                    {
                        J = 0;
                        for (int j = 0; j < unocc_orb; j++)
                        {
                            val1 = EqualFunc((a+n_Electrons)%2, (b+n_Electrons)%2) * EqualFunc((i+n_Electrons)%2, (j+n_Electrons)%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) ((b+n_Electrons)/2), (int) ((i+n_Electrons)/2), (int) ((j+n_Electrons)/2));

                            val2 = EqualFunc((a+n_Electrons)%2, (j+n_Electrons)%2) * EqualFunc((i+n_Electrons)%2, (b+n_Electrons)%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) ((j+n_Electrons)/2), (int) ((i+n_Electrons)/2), (int) ((b+n_Electrons)/2));

                            Rear_cube2(a,i)(B+Speed_Occ,J) = val1 - val2;
                            j++;
                            J++;
                        }
                    }

                    else
                    {
                        J = 0;
                        for (int j = 1; j < unocc_orb; j++)
                        {
                            val1 = EqualFunc((a+n_Electrons)%2, (b+n_Electrons)%2) * EqualFunc((i+n_Electrons)%2, (j+n_Electrons)%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) ((b+n_Electrons)/2), (int) ((i+n_Electrons)/2), (int) ((j+n_Electrons)/2));

                            val2 = EqualFunc((a+n_Electrons)%2, (j+n_Electrons)%2) * EqualFunc((i+n_Electrons)%2, (b+n_Electrons)%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) ((j+n_Electrons)/2), (int) ((i+n_Electrons)/2), (int) ((b+n_Electrons)/2));

                            Rear_cube2(a,i)(B+Speed_Occ,J) = val1 - val2;
                            j++;
                            J++;
                        }
                    }
                }
                b++;
                B++;
            }
        }

        // Compact form of I(i,j,a,b), dont store any zeroes
        Rear_cube3.set_size(n_Electrons, n_Electrons);
        field_filler_matrix2 = zeros(unocc_orb, Speed_Occ);
        for (int i = 0; i < n_Electrons; i++)
        {
            for (int a = 0; a < n_Electrons; a++)
            {
                Rear_cube3(a,i) = field_filler_matrix2;
            }
        }

        for (int a = 0; a < n_Electrons; a++)
        {
            B = 0;
            for (int b = 0; b < unocc_orb; b++)
            {
                for (int i = 0; i < n_Electrons; i++)
                {
                    if ((a+b+i)%2 == 0)
                    {
                        J = 0;
                        for (int j = 0; j < unocc_orb; j++)
                        {
                            val1 = EqualFunc(a%2, (b+n_Electrons)%2) * EqualFunc(i%2, (j+n_Electrons)%2) *
                                    temp4.at((int) (a/2))((int) ((b+n_Electrons)/2), (int) (i/2), (int) ((j+n_Electrons)/2));

                            val2 = EqualFunc(a%2, (j+n_Electrons)%2) * EqualFunc(i%2, (b+n_Electrons)%2) *
                                    temp4.at((int) (a/2))((int) ((j+n_Electrons)/2), (int) (i/2), (int) ((b+n_Electrons)/2));

                            Rear_cube3(a,i)(B,J) = val1-val2;
                            j++;
                            J++;
                        }
                    }

                    else
                    {
                        J = 0;
                        for (int j = 1; j < unocc_orb; j++)
                        {
                            val1 = EqualFunc(a%2, (b+n_Electrons)%2) * EqualFunc(i%2, (j+n_Electrons)%2) *
                                    temp4.at((int) (a/2))((int) ((b+n_Electrons)/2), (int) (i/2), (int) ((j+n_Electrons)/2));

                            val2 = EqualFunc(a%2, (j+n_Electrons)%2) * EqualFunc(i%2, (b+n_Electrons)%2) *
                                    temp4.at((int) (a/2))((int) ((j+n_Electrons)/2), (int) (i/2), (int) ((b+n_Electrons)/2));

                            Rear_cube3(a,i)(B,J) = val1-val2;
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
                    if ((a+b+i)%2 == 0)
                    {
                        J = 0;
                        for (int j = 0; j < unocc_orb; j++)
                        {
                            val1 = EqualFunc(a%2, (b+n_Electrons)%2) * EqualFunc(i%2, (j+n_Electrons)%2) *
                                    temp4.at((int) (a/2))((int) ((b+n_Electrons)/2), (int) (i/2), (int) ((j+n_Electrons)/2));

                            val2 = EqualFunc(a%2, (j+n_Electrons)%2) * EqualFunc(i%2, (b+n_Electrons)%2) *
                                    temp4.at((int) (a/2))((int) ((j+n_Electrons)/2), (int) (i/2), (int) ((b+n_Electrons)/2));

                            Rear_cube3(a,i)(B+Speed_Occ,J) = val1-val2;
                            j++;
                            J++;
                        }
                    }

                    else
                    {
                        J = 0;
                        for (int j = 1; j < unocc_orb; j++)
                        {
                            val1 = EqualFunc(a%2, (b+n_Electrons)%2) * EqualFunc(i%2, (j+n_Electrons)%2) *
                                    temp4.at((int) (a/2))((int) ((b+n_Electrons)/2), (int) (i/2), (int) ((j+n_Electrons)/2));

                            val2 = EqualFunc(a%2, (j+n_Electrons)%2) * EqualFunc(i%2, (b+n_Electrons)%2) *
                                    temp4.at((int) (a/2))((int) ((j+n_Electrons)/2), (int) (i/2), (int) ((b+n_Electrons)/2));

                            Rear_cube3(a,i)(B+Speed_Occ,J) = val1-val2;
                            j++;
                            J++;
                        }
                    }
                }
                b++;
                B++;
            }
        }


        // Compact form of I(a,i,b,c), dont store any zeroes
        Rear_cube4.set_size(unocc_orb, n_Electrons);
        field_filler_matrix2 = zeros(unocc_orb, Speed_Occ);
        for (int i = 0; i < n_Electrons; i++)
        {
            for (int a = 0; a < unocc_orb; a++)
            {
                Rear_cube4(a,i) = field_filler_matrix2;
            }
        }

        for (int a = 0; a < unocc_orb; a++)
        {
            B = 0;
            for (int b = 0; b < unocc_orb; b++)
            {
                for (int i = 0; i < n_Electrons; i++)
                {
                    if ((a+b+i)%2 == 0)
                    {
                        J = 0;
                        for (int j = 0; j < unocc_orb; j++)
                        {
                            val1 = EqualFunc((a+n_Electrons)%2, (b+n_Electrons)%2) * EqualFunc(i%2, (j+n_Electrons)%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) ((b+n_Electrons)/2), (int) (i/2), (int) ((j+n_Electrons)/2));

                            val2 = EqualFunc((a+n_Electrons)%2, (j+n_Electrons)%2) * EqualFunc(i%2, (b+n_Electrons)%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) ((j+n_Electrons)/2), (int) (i/2), (int) ((b+n_Electrons)/2));

                            Rear_cube4(a,i)(B,J) = val1-val2;
                            j++;
                            J++;
                        }
                    }

                    else
                    {
                        J = 0;
                        for (int j = 1; j < unocc_orb; j++)
                        {
                            val1 = EqualFunc((a+n_Electrons)%2, (b+n_Electrons)%2) * EqualFunc(i%2, (j+n_Electrons)%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) ((b+n_Electrons)/2), (int) (i/2), (int) ((j+n_Electrons)/2));

                            val2 = EqualFunc((a+n_Electrons)%2, (j+n_Electrons)%2) * EqualFunc(i%2, (b+n_Electrons)%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) ((j+n_Electrons)/2), (int) (i/2), (int) ((b+n_Electrons)/2));

                            Rear_cube4(a,i)(B,J) = val1-val2;
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
                    if ((a+b+i)%2 == 0)
                    {
                        J = 0;
                        for (int j = 0; j < unocc_orb; j++)
                        {
                            val1 = EqualFunc((a+n_Electrons)%2, (b+n_Electrons)%2) * EqualFunc(i%2, (j+n_Electrons)%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) ((b+n_Electrons)/2), (int) (i/2), (int) ((j+n_Electrons)/2));

                            val2 = EqualFunc((a+n_Electrons)%2, (j+n_Electrons)%2) * EqualFunc(i%2, (b+n_Electrons)%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) ((j+n_Electrons)/2), (int) (i/2), (int) ((b+n_Electrons)/2));

                            Rear_cube4(a,i)(B+Speed_Occ,J) = val1-val2;
                            j++;
                            J++;
                        }
                    }

                    else
                    {
                        J = 0;
                        for (int j = 1; j < unocc_orb; j++)
                        {
                            val1 = EqualFunc((a+n_Electrons)%2, (b+n_Electrons)%2) * EqualFunc(i%2, (j+n_Electrons)%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) ((b+n_Electrons)/2), (int) (i/2), (int) ((j+n_Electrons)/2));

                            val2 = EqualFunc((a+n_Electrons)%2, (j+n_Electrons)%2) * EqualFunc(i%2, (b+n_Electrons)%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) ((j+n_Electrons)/2), (int) (i/2), (int) ((b+n_Electrons)/2));

                            Rear_cube4(a,i)(B+Speed_Occ,J) = val1-val2;
                            j++;
                            J++;
                        }
                    }
                }
                b++;
                B++;
            }
        }

        // Compact form of I(i,j,a,k), dont store any zeroes
        Rear_cube6.set_size(n_Electrons, n_Electrons);
        field_filler_matrix2 = zeros(unocc_orb, n_Electrons/2);
        for (int i = 0; i < n_Electrons; i++)
        {
            for (int a = 0; a < n_Electrons; a++)
            {
                Rear_cube6(a,i) = field_filler_matrix2;
            }
        }

        for (int a = 0; a < n_Electrons; a++)
        {
            B = 0;
            for (int b = 0; b < unocc_orb; b++)
            {
                for (int i = 0; i < n_Electrons; i++)
                {
                    if ((a+b+i)%2 == 0)
                    {
                        J = 0;
                        for (int j = 0; j < n_Electrons; j++)
                        {
                            val1 = EqualFunc(a%2, (b+n_Electrons)%2) * EqualFunc(i%2, j%2) *
                                    temp4.at((int) (a/2))((int) ((b+n_Electrons)/2), (int) (i/2), (int) (j/2));

                            val2 = EqualFunc(a%2, j%2) * EqualFunc(i%2, (b+n_Electrons)%2) *
                                    temp4.at((int) (a/2))((int) (j/2), (int) (i/2), (int) ((b+n_Electrons)/2));

                            Rear_cube6(a,i)(B,J) = val1-val2;
                            j++;
                            J++;
                        }
                    }

                    else
                    {
                        J = 0;
                        for (int j = 1; j < n_Electrons; j++)
                        {
                            val1 = EqualFunc(a%2, (b+n_Electrons)%2) * EqualFunc(i%2, j%2) *
                                    temp4.at((int) (a/2))((int) ((b+n_Electrons)/2), (int) (i/2), (int) (j/2));

                            val2 = EqualFunc(a%2, j%2) * EqualFunc(i%2, (b+n_Electrons)%2) *
                                    temp4.at((int) (a/2))((int) (j/2), (int) (i/2), (int) ((b+n_Electrons)/2));

                            Rear_cube6(a,i)(B,J) = val1-val2;
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
                    if ((a+b+i)%2 == 0)
                    {
                        J = 0;
                        for (int j = 0; j < n_Electrons; j++)
                        {
                            val1 = EqualFunc(a%2, (b+n_Electrons)%2) * EqualFunc(i%2, j%2) *
                                    temp4.at((int) (a/2))((int) ((b+n_Electrons)/2), (int) (i/2), (int) (j/2));

                            val2 = EqualFunc(a%2, j%2) * EqualFunc(i%2, (b+n_Electrons)%2) *
                                    temp4.at((int) (a/2))((int) (j/2), (int) (i/2), (int) ((b+n_Electrons)/2));

                            Rear_cube6(a,i)(B+Speed_Occ,J) = val1-val2;
                            j++;
                            J++;
                        }
                    }

                    else
                    {
                        J = 0;
                        for (int j = 1; j < n_Electrons; j++)
                        {
                            val1 = EqualFunc(a%2, (b+n_Electrons)%2) * EqualFunc(i%2, j%2) *
                                    temp4.at((int) (a/2))((int) ((b+n_Electrons)/2), (int) (i/2), (int) (j/2));

                            val2 = EqualFunc(a%2, j%2) * EqualFunc(i%2, (b+n_Electrons)%2) *
                                    temp4.at((int) (a/2))((int) (j/2), (int) (i/2), (int) ((b+n_Electrons)/2));

                            Rear_cube6(a,i)(B+Speed_Occ,J) = val1-val2;
                            j++;
                            J++;
                        }
                    }
                }
                b++;
                B++;
            }
        }

        // Compact form of I(a,i,b,j), dont store any zeroes
        Rear_cube7.set_size(unocc_orb, n_Electrons);
        field_filler_matrix2 = zeros(unocc_orb, n_Electrons/2);
        for (int i = 0; i < n_Electrons; i++)
        {
            for (int a = 0; a < unocc_orb; a++)
            {
                Rear_cube7(a,i) = field_filler_matrix2;
            }
        }

        for (int a = 0; a < unocc_orb; a++)
        {
            B = 0;
            for (int b = 0; b < unocc_orb; b++)
            {
                for (int i = 0; i < n_Electrons; i++)
                {
                    if ((a+b+i)%2 == 0)
                    {
                        J = 0;
                        for (int j = 0; j < n_Electrons; j++)
                        {
                            val1 = EqualFunc((a+n_Electrons)%2, (b+n_Electrons)%2) * EqualFunc(i%2, j%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) ((b+n_Electrons)/2), (int) (i/2), (int) (j/2));

                            val2 = EqualFunc((a+n_Electrons)%2, j%2) * EqualFunc(i%2, (b+n_Electrons)%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) (j/2), (int) (i/2), (int) ((b+n_Electrons)/2));

                            Rear_cube7(a,i)(B,J) = val1-val2;
                            j++;
                            J++;
                        }
                    }

                    else
                    {
                        J = 0;
                        for (int j = 1; j < n_Electrons; j++)
                        {
                            val1 = EqualFunc((a+n_Electrons)%2, (b+n_Electrons)%2) * EqualFunc(i%2, j%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) ((b+n_Electrons)/2), (int) (i/2), (int) (j/2));

                            val2 = EqualFunc((a+n_Electrons)%2, j%2) * EqualFunc(i%2, (b+n_Electrons)%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) (j/2), (int) (i/2), (int) ((b+n_Electrons)/2));

                            Rear_cube7(a,i)(B,J) = val1-val2;
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
                    if ((a+b+i)%2 == 0)
                    {
                        J = 0;
                        for (int j = 0; j < n_Electrons; j++)
                        {
                            val1 = EqualFunc((a+n_Electrons)%2, (b+n_Electrons)%2) * EqualFunc(i%2, j%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) ((b+n_Electrons)/2), (int) (i/2), (int) (j/2));

                            val2 = EqualFunc((a+n_Electrons)%2, j%2) * EqualFunc(i%2, (b+n_Electrons)%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) (j/2), (int) (i/2), (int) ((b+n_Electrons)/2));

                            Rear_cube7(a,i)(B+Speed_Occ,J) = val1-val2;
                            j++;
                            J++;
                        }
                    }

                    else
                    {
                        J = 0;
                        for (int j = 1; j < n_Electrons; j++)
                        {
                            val1 = EqualFunc((a+n_Electrons)%2, (b+n_Electrons)%2) * EqualFunc(i%2, j%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) ((b+n_Electrons)/2), (int) (i/2), (int) (j/2));

                            val2 = EqualFunc((a+n_Electrons)%2, j%2) * EqualFunc(i%2, (b+n_Electrons)%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) (j/2), (int) (i/2), (int) ((b+n_Electrons)/2));

                            Rear_cube7(a,i)(B+Speed_Occ,J) = val1-val2;
                            j++;
                            J++;
                        }
                    }
                }
                b++;
                B++;
            }
        }

        // Compact form of I(i,j,k,l), dont store any zeroes
        Rear_cube8.set_size(n_Electrons, n_Electrons);
        field_filler_matrix2 = zeros(n_Electrons, n_Electrons/2);
        for (int i = 0; i < n_Electrons; i++)
        {
            for (int a = 0; a < n_Electrons; a++)
            {
                Rear_cube8(a,i) = field_filler_matrix2;
            }
        }

        for (int a = 0; a < n_Electrons; a++)
        {
            B = 0;
            for (int b = 0; b < n_Electrons; b++)
            {
                for (int i = 0; i < n_Electrons; i++)
                {
                    if ((a+b+i)%2 == 0)
                    {
                        J = 0;
                        for (int j = 0; j < n_Electrons; j++)
                        {
                            val1 = EqualFunc(a%2, b%2) * EqualFunc(i%2, j%2) *
                                    temp4.at((int) (a/2))((int) (b/2), (int) (i/2), (int) (j/2));

                            val2 = EqualFunc(a%2, j%2) * EqualFunc(i%2, b%2) *
                                    temp4.at((int) (a/2))((int) (j/2), (int) (i/2), (int) (b/2));

                            Rear_cube8(a,i)(B,J) = val1-val2;
                            j++;
                            J++;
                        }
                    }

                    else
                    {
                        J = 0;
                        for (int j = 1; j < n_Electrons; j++)
                        {
                            val1 = EqualFunc(a%2, b%2) * EqualFunc(i%2, j%2) *
                                    temp4.at((int) (a/2))((int) (b/2), (int) (i/2), (int) (j/2));

                            val2 = EqualFunc(a%2, j%2) * EqualFunc(i%2, b%2) *
                                    temp4.at((int) (a/2))((int) (j/2), (int) (i/2), (int) (b/2));

                            Rear_cube8(a,i)(B,J) = val1-val2;
                            j++;
                            J++;
                        }
                    }
                }
                b++;
                B++;
            }

            B = 0;
            for (int b = 1; b < n_Electrons; b++)
            {
                for (int i = 0; i < n_Electrons; i++)
                {
                    if ((a+b+i)%2 == 0)
                    {
                        J = 0;
                        for (int j = 0; j < n_Electrons; j++)
                        {
                            val1 = EqualFunc(a%2, b%2) * EqualFunc(i%2, j%2) *
                                    temp4.at((int) (a/2))((int) (b/2), (int) (i/2), (int) (j/2));

                            val2 = EqualFunc(a%2, j%2) * EqualFunc(i%2, b%2) *
                                    temp4.at((int) (a/2))((int) (j/2), (int) (i/2), (int) (b/2));

                            Rear_cube8(a,i)(B+n_Electrons/2,J) = val1-val2;
                            j++;
                            J++;
                        }
                    }

                    else
                    {
                        J = 0;
                        for (int j = 1; j < n_Electrons; j++)
                        {
                            val1 = EqualFunc(a%2, b%2) * EqualFunc(i%2, j%2) *
                                    temp4.at((int) (a/2))((int) (b/2), (int) (i/2), (int) (j/2));

                            val2 = EqualFunc(a%2, j%2) * EqualFunc(i%2, b%2) *
                                    temp4.at((int) (a/2))((int) (j/2), (int) (i/2), (int) (b/2));

                            Rear_cube8(a,i)(B+n_Electrons/2,J) = val1-val2;
                            j++;
                            J++;
                        }
                    }
                }
                b++;
                B++;
            }
        }

        // Compact form of I(a,b,c,i), dont store any zeroes
        Rear_cube9.set_size(unocc_orb, unocc_orb);
        field_filler_matrix2 = zeros(unocc_orb, n_Electrons/2);
        for (int i = 0; i < unocc_orb; i++)
        {
            for (int a = 0; a < unocc_orb; a++)
            {
                Rear_cube9(a,i) = field_filler_matrix2;
            }
        }

        for (int a = 0; a < unocc_orb; a++)
        {
            B = 0;
            for (int b = 0; b < unocc_orb; b++)
            {
                for (int i = 0; i < unocc_orb; i++)
                {
                    if ((a+b+i)%2 == 0)
                    {
                        J = 0;
                        for (int j = 0; j < n_Electrons; j++)
                        {
                            val1 = EqualFunc((a+n_Electrons)%2, (b+n_Electrons)%2) * EqualFunc((i+n_Electrons)%2, j%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) ((b+n_Electrons)/2), (int) ((i+n_Electrons)/2), (int) (j/2));

                            val2 = EqualFunc((a+n_Electrons)%2, j%2) * EqualFunc((i+n_Electrons)%2, (b+n_Electrons)%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) (j/2), (int) ((i+n_Electrons)/2), (int) ((b+n_Electrons)/2));

                            Rear_cube9(a,i)(B,J) = val1-val2;
                            j++;
                            J++;
                        }
                    }

                    else
                    {
                        J = 0;
                        for (int j = 1; j < n_Electrons; j++)
                        {
                            val1 = EqualFunc((a+n_Electrons)%2, (b+n_Electrons)%2) * EqualFunc((i+n_Electrons)%2, j%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) ((b+n_Electrons)/2), (int) ((i+n_Electrons)/2), (int) (j/2));

                            val2 = EqualFunc((a+n_Electrons)%2, j%2) * EqualFunc((i+n_Electrons)%2, (b+n_Electrons)%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) (j/2), (int) ((i+n_Electrons)/2), (int) ((b+n_Electrons)/2));

                            Rear_cube9(a,i)(B,J) = val1-val2;
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
                for (int i = 0; i < unocc_orb; i++)
                {
                    if ((a+b+i)%2 == 0)
                    {
                        J = 0;
                        for (int j = 0; j < n_Electrons; j++)
                        {
                            val1 = EqualFunc((a+n_Electrons)%2, (b+n_Electrons)%2) * EqualFunc((i+n_Electrons)%2, j%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) ((b+n_Electrons)/2), (int) ((i+n_Electrons)/2), (int) (j/2));

                            val2 = EqualFunc((a+n_Electrons)%2, j%2) * EqualFunc((i+n_Electrons)%2, (b+n_Electrons)%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) (j/2), (int) ((i+n_Electrons)/2), (int) ((b+n_Electrons)/2));

                            Rear_cube9(a,i)(B+Speed_Occ,J) = val1-val2;
                            j++;
                            J++;
                        }
                    }

                    else
                    {
                        J = 0;
                        for (int j = 1; j < n_Electrons; j++)
                        {
                            val1 = EqualFunc((a+n_Electrons)%2, (b+n_Electrons)%2) * EqualFunc((i+n_Electrons)%2, j%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) ((b+n_Electrons)/2), (int) ((i+n_Electrons)/2), (int) (j/2));

                            val2 = EqualFunc((a+n_Electrons)%2, j%2) * EqualFunc((i+n_Electrons)%2, (b+n_Electrons)%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) (j/2), (int) ((i+n_Electrons)/2), (int) ((b+n_Electrons)/2));

                            Rear_cube9(a,i)(B+Speed_Occ,J) = val1-val2;
                            j++;
                            J++;
                        }
                    }
                }
                b++;
                B++;
            }
        }

        // Compact form of I(a,b,c,i), dont store any zeroes
        Rear_cube10.set_size(unocc_orb, n_Electrons);
        field_filler_matrix2 = zeros(unocc_orb, n_Electrons/2);
        for (int i = 0; i < n_Electrons; i++)
        {
            for (int a = 0; a < unocc_orb; a++)
            {
                Rear_cube10(a,i) = field_filler_matrix2;
            }
        }

        for (int a = 0; a < unocc_orb; a++)
        {
            B = 0;
            for (int b = 0; b < unocc_orb; b++)
            {
                for (int i = 0; i < n_Electrons; i++)
                {
                    if ((a+b+i)%2 == 0)
                    {
                        J = 0;
                        for (int j = 0; j < n_Electrons; j++)
                        {
                            val1 = EqualFunc((a+n_Electrons)%2, (b+n_Electrons)%2) * EqualFunc(j%2, i%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) ((b+n_Electrons)/2), (int) (j/2), (int) (i/2));

                            val2 = EqualFunc((a+n_Electrons)%2, i%2) * EqualFunc(j%2, (b+n_Electrons)%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) (i/2), (int) (j/2), (int) ((b+n_Electrons)/2));

                            Rear_cube10(a,i)(B,J) = val1-val2;
                            j++;
                            J++;
                        }
                    }

                    else
                    {
                        J = 0;
                        for (int j = 1; j < n_Electrons; j++)
                        {
                            val1 = EqualFunc((a+n_Electrons)%2, (b+n_Electrons)%2) * EqualFunc(j%2, i%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) ((b+n_Electrons)/2), (int) (j/2), (int) (i/2));

                            val2 = EqualFunc((a+n_Electrons)%2, i%2) * EqualFunc(j%2, (b+n_Electrons)%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) (i/2), (int) (j/2), (int) ((b+n_Electrons)/2));

                            Rear_cube10(a,i)(B,J) = val1-val2;
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
                    if ((a+b+i)%2 == 0)
                    {
                        J = 0;
                        for (int j = 0; j < n_Electrons; j++)
                        {
                            val1 = EqualFunc((a+n_Electrons)%2, (b+n_Electrons)%2) * EqualFunc(j%2, i%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) ((b+n_Electrons)/2), (int) (j/2), (int) (i/2));

                            val2 = EqualFunc((a+n_Electrons)%2, i%2) * EqualFunc(j%2, (b+n_Electrons)%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) (i/2), (int) (j/2), (int) ((b+n_Electrons)/2));

                            Rear_cube10(a,i)(B+Speed_Occ,J) = val1-val2;
                            j++;
                            J++;
                        }
                    }

                    else
                    {
                        J = 0;
                        for (int j = 1; j < n_Electrons; j++)
                        {
                            val1 = EqualFunc((a+n_Electrons)%2, (b+n_Electrons)%2) * EqualFunc(j%2, i%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) ((b+n_Electrons)/2), (int) (j/2), (int) (i/2));

                            val2 = EqualFunc((a+n_Electrons)%2, i%2) * EqualFunc(j%2, (b+n_Electrons)%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) (i/2), (int) (j/2), (int) ((b+n_Electrons)/2));

                            Rear_cube10(a,i)(B+Speed_Occ,J) = val1-val2;
                            j++;
                            J++;
                        }
                    }
                }
                b++;
                B++;
            }
        }

        temp4.clear();
    }

    if (Version==1)
    {
        temp_cube3.set_size(2*Matrix_Size, 2*Matrix_Size);
        field_filler_matrix2 = zeros(2*Matrix_Size, 2*Matrix_Size);

        for (int i = 0; i < 2*Matrix_Size; i++)
        {
            for (int j = 0; j < 2*Matrix_Size; j++)
            {
                temp_cube3(i,j) = field_filler_matrix2;
            }
        }

        temp_cube33.set_size(2*Matrix_Size-n_Electrons, n_Electrons);
        field_filler_matrix2 = zeros(2*Matrix_Size-n_Electrons, n_Electrons);

        for (int i = 0; i < (2*Matrix_Size-n_Electrons); i++)
        {
            for (int j = 0; j < n_Electrons; j++)
            {
                temp_cube33(i,j) = field_filler_matrix2;
            }
        }

        for (int i = 0; i < (2*Matrix_Size); i++)
        {
            for (int j = 0; j < (2*Matrix_Size); j++)
            {
                for (int k = 0; k < (2*Matrix_Size); k++)
                {
                    for (int l = 0; l < (2*Matrix_Size); l++)
                    {
                        val1 = EqualFunc(i%2, k%2) * EqualFunc(j%2, l%2) *
                                temp4.at((int) (i/2))((int) (k/2), (int) (j/2), (int) (l/2));

                        val2 = EqualFunc(i%2, l%2) * EqualFunc(j%2, k%2) *
                                temp4.at((int) (i/2))((int) (l/2), (int) (j/2), (int) (k/2));

                        temp_cube3.at(i,j)(k,l) = val1-val2;
                    }
                }
            }
        }

        int temp_i;
        int temp_j;
        for (int i = n_Electrons; i < (2*Matrix_Size); i++)
        {
            temp_i = i - n_Electrons;
            for (int j = n_Electrons; j < (2*Matrix_Size); j++)
            {
                temp_j = j - n_Electrons;
                for (int k = 0; k < n_Electrons; k++)
                {
                    for (int l = 0; l < n_Electrons; l++)
                    {
                        temp_cube33(temp_i,k)(temp_j,l) = temp_cube3(i,j)(k,l);
                    }
                }
            }
        }
        temp4.clear();
    }
}

void coupled_cluster_integrals::Prepear_Integrals()
{
    // This is ment to be an efficient and good way to calculate the double integrals we will need

    // Transformation from Atomic Orbital (AO) to Molecule Orbital (MO)
    // Sometimes this procedure is refered to as AO2MO or AOtoMO
    cube temp_cube = zeros(Matrix_Size, Matrix_Size, Matrix_Size);
    mat field_filler_matrix2;

    vector <cube> temp;
    vector <cube> temp2;
    vector <cube> temp3;

    for (int i = 0; i < Matrix_Size; i++)
    {
        temp.push_back(temp_cube);
    }

    for (int a = 0; a < Matrix_Size; a++)
    {
        for (int i = 0; i < Matrix_Size; i++)
        {
            for (int j = 0; j < Matrix_Size; j++)
            {
                for (int k = 0; k < Matrix_Size; k++)
                {
                    for (int l = 0; l < Matrix_Size; l++)
                    {
                        temp.at(a)(j,k,l) += c(i,a) * Indexed_integrals(Return_Integral_Index(i,j,k,l));//Integrals.at(i)(k,j,l) * c(i,a);
                    }
                }
            }
        }
    }
    Integrals.clear();

    for (int i = 0; i < Matrix_Size; i++)
    {
        temp2.push_back(temp_cube);
    }

    for (int a = 0; a < Matrix_Size; a++)
    {
        for (int b = 0; b < Matrix_Size; b++)
        {
            for (int i = 0; i < Matrix_Size; i++)
            {
                for (int j = 0; j < Matrix_Size; j++)
                {
                    for (int k = 0; k < Matrix_Size; k++)
                    {
                        temp2.at(a)(b,j,k) += c(i,b) * temp.at(a)(i,j,k);
                    }
                }
            }
        }
    }

    temp.clear();

    for (int i = 0; i < Matrix_Size; i++)
    {
        temp3.push_back(temp_cube);
    }

    for (int a = 0; a < Matrix_Size; a++)
    {
        for (int b = 0; b < Matrix_Size; b++)
        {
            for (int g = 0; g < Matrix_Size; g++)
            {
                for (int i = 0; i < Matrix_Size; i++)
                {
                    for (int j = 0; j < Matrix_Size; j++)
                    {
                        temp3.at(a)(b,g,j) += c(i,g) * temp2.at(a)(b,i,j);
                    }
                }
            }
        }
    }

    temp2.clear();

    for (int i = 0; i < Matrix_Size; i++)
    {
        temp4.push_back(temp_cube);
    }

    for (int a = 0; a < Matrix_Size; a++)
    {
        for (int b = 0; b < Matrix_Size; b++)
        {
            for (int g = 0; g < Matrix_Size; g++)
            {
                for (int h = 0; h < Matrix_Size; h++)
                {
                    for (int i = 0; i < Matrix_Size; i++)
                    {
                        temp4.at(a)(b,g,h) += c(i,h) * temp3.at(a)(b,g,i);
                    }
                }
            }
        }
    }

    temp3.clear();

    double val1, val2;

    int unocc_orb = 2*Matrix_Size - n_Electrons;
    int Speed_Occ = unocc_orb/2;

        Rear_cube.set_size(unocc_orb, n_Electrons);
        field_filler_matrix2 = zeros(unocc_orb, n_Electrons/2);
        for (int i = 0; i < n_Electrons; i++)
        {
            for (int a = 0; a < unocc_orb; a++)
            {
                Rear_cube(a,i) = field_filler_matrix2;
            }
        }

        // Using properties of the integrals we can store only half the elements


        // Compact form of I(a,b,i,j), dont store any zeroes
        int B, J;
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
                            val1 = EqualFunc((a+n_Electrons)%2, i%2) * EqualFunc((b+n_Electrons)%2, j%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) (i/2), (int) ((b+n_Electrons)/2), (int) (j/2));

                            val2 = EqualFunc((a+n_Electrons)%2, j%2) * EqualFunc((b+n_Electrons)%2, i%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) (j/2), (int) ((b+n_Electrons)/2), (int) (i/2));

                            Rear_cube(a,i)(B,J) = val1-val2; //temp_cube3(a+n_Electrons,b+n_Electrons)(i,j);
                            j++;
                            J++;
                        }
                    }

                    else
                    {
                        J = 0;
                        for (int j = 1; j < n_Electrons; j++)
                        {
                            val1 = EqualFunc((a+n_Electrons)%2, i%2) * EqualFunc((b+n_Electrons)%2, j%2) *
                                     temp4.at((int) ((a+n_Electrons)/2))((int) (i/2), (int) ((b+n_Electrons)/2), (int) (j/2));


                            val2 = EqualFunc((a+n_Electrons)%2, j%2) * EqualFunc((b+n_Electrons)%2, i%2) *
                                     temp4.at((int) ((a+n_Electrons)/2))((int) (j/2), (int) ((b+n_Electrons)/2), (int) (i/2));

                             Rear_cube(a,i)(B,J) = val1-val2;
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
                            val1 = EqualFunc((a+n_Electrons)%2, i%2) * EqualFunc((b+n_Electrons)%2, j%2) *
                                     temp4.at((int) ((a+n_Electrons)/2))((int) (i/2), (int) ((b+n_Electrons)/2), (int) (j/2));

                            val2 = EqualFunc((a+n_Electrons)%2, j%2) * EqualFunc((b+n_Electrons)%2, i%2) *
                                     temp4.at((int) ((a+n_Electrons)/2))((int) (j/2), (int) ((b+n_Electrons)/2), (int) (i/2));

                            Rear_cube(a,i)(B+Speed_Occ,J) = val1-val2;
                            j++;
                            J++;
                        }
                    }

                    else
                    {
                        J = 0;
                        for (int j = 1; j < n_Electrons; j++)
                        {
                            val1 = EqualFunc((a+n_Electrons)%2, i%2) * EqualFunc((b+n_Electrons)%2, j%2) *
                                     temp4.at((int) ((a+n_Electrons)/2))((int) (i/2), (int) ((b+n_Electrons)/2), (int) (j/2));

                            val2 = EqualFunc((a+n_Electrons)%2, j%2) * EqualFunc((b+n_Electrons)%2, i%2) *
                                     temp4.at((int) ((a+n_Electrons)/2))((int) (j/2), (int) ((b+n_Electrons)/2), (int) (i/2));

                            Rear_cube(a,i)(B+Speed_Occ,J) = val1-val2;
                            j++;
                            J++;
                        }
                    }
                }
                b++;
                B++;
            }
        }


        // Compact form of I(a,b,c,d), dont store any zeroes
        Rear_cube2.set_size(unocc_orb, unocc_orb);
        field_filler_matrix2 = zeros(unocc_orb, Speed_Occ);
        for (int i = 0; i < unocc_orb; i++)
        {
            for (int a = i+1; a < unocc_orb; a++)
            {
                Rear_cube2(i,a) = field_filler_matrix2;
            }
        }

        for (int a = 0; a < unocc_orb; a++)
        {
            B = 0;
            for (int b = 0; b < unocc_orb; b++)
            {
                for (int i = a+1; i < unocc_orb; i++)
                {
                    if ((a+b+i)%2 == 0)
                    {
                        J = 0;
                        for (int j = 0; j < unocc_orb; j++)
                        {
                            val1 = EqualFunc((a+n_Electrons)%2, (b+n_Electrons)%2) * EqualFunc((i+n_Electrons)%2, (j+n_Electrons)%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) ((b+n_Electrons)/2), (int) ((i+n_Electrons)/2), (int) ((j+n_Electrons)/2));

                            val2 = EqualFunc((a+n_Electrons)%2, (j+n_Electrons)%2) * EqualFunc((i+n_Electrons)%2, (b+n_Electrons)%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) ((j+n_Electrons)/2), (int) ((i+n_Electrons)/2), (int) ((b+n_Electrons)/2));

                            Rear_cube2(a,i)(B,J) = val1 - val2;
                            j++;
                            J++;
                        }
                    }

                    else
                    {
                        J = 0;
                        for (int j = 1; j < unocc_orb; j++)
                        {
                            val1 = EqualFunc((a+n_Electrons)%2, (b+n_Electrons)%2) * EqualFunc((i+n_Electrons)%2, (j+n_Electrons)%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) ((b+n_Electrons)/2), (int) ((i+n_Electrons)/2), (int) ((j+n_Electrons)/2));

                            val2 = EqualFunc((a+n_Electrons)%2, (j+n_Electrons)%2) * EqualFunc((i+n_Electrons)%2, (b+n_Electrons)%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) ((j+n_Electrons)/2), (int) ((i+n_Electrons)/2), (int) ((b+n_Electrons)/2));

                            Rear_cube2(a,i)(B,J) = val1 - val2;
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
                for (int i = a+1; i < unocc_orb; i++)
                {
                    if ((a+b+i)%2 == 0)
                    {
                        J = 0;
                        for (int j = 0; j < unocc_orb; j++)
                        {
                            val1 = EqualFunc((a+n_Electrons)%2, (b+n_Electrons)%2) * EqualFunc((i+n_Electrons)%2, (j+n_Electrons)%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) ((b+n_Electrons)/2), (int) ((i+n_Electrons)/2), (int) ((j+n_Electrons)/2));

                            val2 = EqualFunc((a+n_Electrons)%2, (j+n_Electrons)%2) * EqualFunc((i+n_Electrons)%2, (b+n_Electrons)%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) ((j+n_Electrons)/2), (int) ((i+n_Electrons)/2), (int) ((b+n_Electrons)/2));

                            Rear_cube2(a,i)(B+Speed_Occ,J) = val1 - val2;
                            j++;
                            J++;
                        }
                    }

                    else
                    {
                        J = 0;
                        for (int j = 1; j < unocc_orb; j++)
                        {
                            val1 = EqualFunc((a+n_Electrons)%2, (b+n_Electrons)%2) * EqualFunc((i+n_Electrons)%2, (j+n_Electrons)%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) ((b+n_Electrons)/2), (int) ((i+n_Electrons)/2), (int) ((j+n_Electrons)/2));

                            val2 = EqualFunc((a+n_Electrons)%2, (j+n_Electrons)%2) * EqualFunc((i+n_Electrons)%2, (b+n_Electrons)%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) ((j+n_Electrons)/2), (int) ((i+n_Electrons)/2), (int) ((b+n_Electrons)/2));

                            Rear_cube2(a,i)(B+Speed_Occ,J) = val1 - val2;
                            j++;
                            J++;
                        }
                    }
                }
                b++;
                B++;
            }
        }

        // Compact form of I(i,j,a,b), dont store any zeroes
        Rear_cube3.set_size(n_Electrons, n_Electrons);
        field_filler_matrix2 = zeros(unocc_orb, Speed_Occ);
        for (int i = 0; i < n_Electrons; i++)
        {
            for (int a = 0; a < n_Electrons; a++)
            {
                Rear_cube3(a,i) = field_filler_matrix2;
            }
        }

        for (int a = 0; a < n_Electrons; a++)
        {
            B = 0;
            for (int b = 0; b < unocc_orb; b++)
            {
                for (int i = 0; i < n_Electrons; i++)
                {
                    if ((a+b+i)%2 == 0)
                    {
                        J = 0;
                        for (int j = 0; j < unocc_orb; j++)
                        {
                            val1 = EqualFunc(a%2, (b+n_Electrons)%2) * EqualFunc(i%2, (j+n_Electrons)%2) *
                                    temp4.at((int) (a/2))((int) ((b+n_Electrons)/2), (int) (i/2), (int) ((j+n_Electrons)/2));

                            val2 = EqualFunc(a%2, (j+n_Electrons)%2) * EqualFunc(i%2, (b+n_Electrons)%2) *
                                    temp4.at((int) (a/2))((int) ((j+n_Electrons)/2), (int) (i/2), (int) ((b+n_Electrons)/2));

                            Rear_cube3(a,i)(B,J) = val1-val2;
                            j++;
                            J++;
                        }
                    }

                    else
                    {
                        J = 0;
                        for (int j = 1; j < unocc_orb; j++)
                        {
                            val1 = EqualFunc(a%2, (b+n_Electrons)%2) * EqualFunc(i%2, (j+n_Electrons)%2) *
                                    temp4.at((int) (a/2))((int) ((b+n_Electrons)/2), (int) (i/2), (int) ((j+n_Electrons)/2));

                            val2 = EqualFunc(a%2, (j+n_Electrons)%2) * EqualFunc(i%2, (b+n_Electrons)%2) *
                                    temp4.at((int) (a/2))((int) ((j+n_Electrons)/2), (int) (i/2), (int) ((b+n_Electrons)/2));

                            Rear_cube3(a,i)(B,J) = val1-val2;
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
                    if ((a+b+i)%2 == 0)
                    {
                        J = 0;
                        for (int j = 0; j < unocc_orb; j++)
                        {
                            val1 = EqualFunc(a%2, (b+n_Electrons)%2) * EqualFunc(i%2, (j+n_Electrons)%2) *
                                    temp4.at((int) (a/2))((int) ((b+n_Electrons)/2), (int) (i/2), (int) ((j+n_Electrons)/2));

                            val2 = EqualFunc(a%2, (j+n_Electrons)%2) * EqualFunc(i%2, (b+n_Electrons)%2) *
                                    temp4.at((int) (a/2))((int) ((j+n_Electrons)/2), (int) (i/2), (int) ((b+n_Electrons)/2));

                            Rear_cube3(a,i)(B+Speed_Occ,J) = val1-val2;
                            j++;
                            J++;
                        }
                    }

                    else
                    {
                        J = 0;
                        for (int j = 1; j < unocc_orb; j++)
                        {
                            val1 = EqualFunc(a%2, (b+n_Electrons)%2) * EqualFunc(i%2, (j+n_Electrons)%2) *
                                    temp4.at((int) (a/2))((int) ((b+n_Electrons)/2), (int) (i/2), (int) ((j+n_Electrons)/2));

                            val2 = EqualFunc(a%2, (j+n_Electrons)%2) * EqualFunc(i%2, (b+n_Electrons)%2) *
                                    temp4.at((int) (a/2))((int) ((j+n_Electrons)/2), (int) (i/2), (int) ((b+n_Electrons)/2));

                            Rear_cube3(a,i)(B+Speed_Occ,J) = val1-val2;
                            j++;
                            J++;
                        }
                    }
                }
                b++;
                B++;
            }
        }


        // Compact form of I(a,i,b,c), dont store any zeroes
        Rear_cube4.set_size(unocc_orb, n_Electrons);
        field_filler_matrix2 = zeros(unocc_orb, Speed_Occ);
        for (int i = 0; i < n_Electrons; i++)
        {
            for (int a = 0; a < unocc_orb; a++)
            {
                Rear_cube4(a,i) = field_filler_matrix2;
            }
        }

        for (int a = 0; a < unocc_orb; a++)
        {
            B = 0;
            for (int b = 0; b < unocc_orb; b++)
            {
                for (int i = 0; i < n_Electrons; i++)
                {
                    if ((a+b+i)%2 == 0)
                    {
                        J = 0;
                        for (int j = 0; j < unocc_orb; j++)
                        {
                            val1 = EqualFunc((a+n_Electrons)%2, (b+n_Electrons)%2) * EqualFunc(i%2, (j+n_Electrons)%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) ((b+n_Electrons)/2), (int) (i/2), (int) ((j+n_Electrons)/2));

                            val2 = EqualFunc((a+n_Electrons)%2, (j+n_Electrons)%2) * EqualFunc(i%2, (b+n_Electrons)%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) ((j+n_Electrons)/2), (int) (i/2), (int) ((b+n_Electrons)/2));

                            Rear_cube4(a,i)(B,J) = val1-val2;
                            j++;
                            J++;
                        }
                    }

                    else
                    {
                        J = 0;
                        for (int j = 1; j < unocc_orb; j++)
                        {
                            val1 = EqualFunc((a+n_Electrons)%2, (b+n_Electrons)%2) * EqualFunc(i%2, (j+n_Electrons)%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) ((b+n_Electrons)/2), (int) (i/2), (int) ((j+n_Electrons)/2));

                            val2 = EqualFunc((a+n_Electrons)%2, (j+n_Electrons)%2) * EqualFunc(i%2, (b+n_Electrons)%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) ((j+n_Electrons)/2), (int) (i/2), (int) ((b+n_Electrons)/2));

                            Rear_cube4(a,i)(B,J) = val1-val2;
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
                    if ((a+b+i)%2 == 0)
                    {
                        J = 0;
                        for (int j = 0; j < unocc_orb; j++)
                        {
                            val1 = EqualFunc((a+n_Electrons)%2, (b+n_Electrons)%2) * EqualFunc(i%2, (j+n_Electrons)%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) ((b+n_Electrons)/2), (int) (i/2), (int) ((j+n_Electrons)/2));

                            val2 = EqualFunc((a+n_Electrons)%2, (j+n_Electrons)%2) * EqualFunc(i%2, (b+n_Electrons)%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) ((j+n_Electrons)/2), (int) (i/2), (int) ((b+n_Electrons)/2));

                            Rear_cube4(a,i)(B+Speed_Occ,J) = val1-val2;
                            j++;
                            J++;
                        }
                    }

                    else
                    {
                        J = 0;
                        for (int j = 1; j < unocc_orb; j++)
                        {
                            val1 = EqualFunc((a+n_Electrons)%2, (b+n_Electrons)%2) * EqualFunc(i%2, (j+n_Electrons)%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) ((b+n_Electrons)/2), (int) (i/2), (int) ((j+n_Electrons)/2));

                            val2 = EqualFunc((a+n_Electrons)%2, (j+n_Electrons)%2) * EqualFunc(i%2, (b+n_Electrons)%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) ((j+n_Electrons)/2), (int) (i/2), (int) ((b+n_Electrons)/2));

                            Rear_cube4(a,i)(B+Speed_Occ,J) = val1-val2;
                            j++;
                            J++;
                        }
                    }
                }
                b++;
                B++;
            }
        }

        // Compact form of I(i,j,a,k), dont store any zeroes
        Rear_cube6.set_size(n_Electrons, n_Electrons);
        field_filler_matrix2 = zeros(unocc_orb, n_Electrons/2);
        for (int i = 0; i < n_Electrons; i++)
        {
            for (int a = 0; a < n_Electrons; a++)
            {
                Rear_cube6(a,i) = field_filler_matrix2;
            }
        }

        for (int a = 0; a < n_Electrons; a++)
        {
            B = 0;
            for (int b = 0; b < unocc_orb; b++)
            {
                for (int i = 0; i < n_Electrons; i++)
                {
                    if ((a+b+i)%2 == 0)
                    {
                        J = 0;
                        for (int j = 0; j < n_Electrons; j++)
                        {
                            val1 = EqualFunc(a%2, (b+n_Electrons)%2) * EqualFunc(i%2, j%2) *
                                    temp4.at((int) (a/2))((int) ((b+n_Electrons)/2), (int) (i/2), (int) (j/2));

                            val2 = EqualFunc(a%2, j%2) * EqualFunc(i%2, (b+n_Electrons)%2) *
                                    temp4.at((int) (a/2))((int) (j/2), (int) (i/2), (int) ((b+n_Electrons)/2));

                            Rear_cube6(a,i)(B,J) = val1-val2;
                            j++;
                            J++;
                        }
                    }

                    else
                    {
                        J = 0;
                        for (int j = 1; j < n_Electrons; j++)
                        {
                            val1 = EqualFunc(a%2, (b+n_Electrons)%2) * EqualFunc(i%2, j%2) *
                                    temp4.at((int) (a/2))((int) ((b+n_Electrons)/2), (int) (i/2), (int) (j/2));

                            val2 = EqualFunc(a%2, j%2) * EqualFunc(i%2, (b+n_Electrons)%2) *
                                    temp4.at((int) (a/2))((int) (j/2), (int) (i/2), (int) ((b+n_Electrons)/2));

                            Rear_cube6(a,i)(B,J) = val1-val2;
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
                    if ((a+b+i)%2 == 0)
                    {
                        J = 0;
                        for (int j = 0; j < n_Electrons; j++)
                        {
                            val1 = EqualFunc(a%2, (b+n_Electrons)%2) * EqualFunc(i%2, j%2) *
                                    temp4.at((int) (a/2))((int) ((b+n_Electrons)/2), (int) (i/2), (int) (j/2));

                            val2 = EqualFunc(a%2, j%2) * EqualFunc(i%2, (b+n_Electrons)%2) *
                                    temp4.at((int) (a/2))((int) (j/2), (int) (i/2), (int) ((b+n_Electrons)/2));

                            Rear_cube6(a,i)(B+Speed_Occ,J) = val1-val2;
                            j++;
                            J++;
                        }
                    }

                    else
                    {
                        J = 0;
                        for (int j = 1; j < n_Electrons; j++)
                        {
                            val1 = EqualFunc(a%2, (b+n_Electrons)%2) * EqualFunc(i%2, j%2) *
                                    temp4.at((int) (a/2))((int) ((b+n_Electrons)/2), (int) (i/2), (int) (j/2));

                            val2 = EqualFunc(a%2, j%2) * EqualFunc(i%2, (b+n_Electrons)%2) *
                                    temp4.at((int) (a/2))((int) (j/2), (int) (i/2), (int) ((b+n_Electrons)/2));

                            Rear_cube6(a,i)(B+Speed_Occ,J) = val1-val2;
                            j++;
                            J++;
                        }
                    }
                }
                b++;
                B++;
            }
        }

        // Compact form of I(a,i,b,j), dont store any zeroes
        Rear_cube7.set_size(unocc_orb, n_Electrons);
        field_filler_matrix2 = zeros(unocc_orb, n_Electrons/2);
        for (int i = 0; i < n_Electrons; i++)
        {
            for (int a = 0; a < unocc_orb; a++)
            {
                Rear_cube7(a,i) = field_filler_matrix2;
            }
        }

        for (int a = 0; a < unocc_orb; a++)
        {
            B = 0;
            for (int b = 0; b < unocc_orb; b++)
            {
                for (int i = 0; i < n_Electrons; i++)
                {
                    if ((a+b+i)%2 == 0)
                    {
                        J = 0;
                        for (int j = 0; j < n_Electrons; j++)
                        {
                            val1 = EqualFunc((a+n_Electrons)%2, (b+n_Electrons)%2) * EqualFunc(i%2, j%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) ((b+n_Electrons)/2), (int) (i/2), (int) (j/2));

                            val2 = EqualFunc((a+n_Electrons)%2, j%2) * EqualFunc(i%2, (b+n_Electrons)%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) (j/2), (int) (i/2), (int) ((b+n_Electrons)/2));

                            Rear_cube7(a,i)(B,J) = val1-val2;
                            j++;
                            J++;
                        }
                    }

                    else
                    {
                        J = 0;
                        for (int j = 1; j < n_Electrons; j++)
                        {
                            val1 = EqualFunc((a+n_Electrons)%2, (b+n_Electrons)%2) * EqualFunc(i%2, j%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) ((b+n_Electrons)/2), (int) (i/2), (int) (j/2));

                            val2 = EqualFunc((a+n_Electrons)%2, j%2) * EqualFunc(i%2, (b+n_Electrons)%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) (j/2), (int) (i/2), (int) ((b+n_Electrons)/2));

                            Rear_cube7(a,i)(B,J) = val1-val2;
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
                    if ((a+b+i)%2 == 0)
                    {
                        J = 0;
                        for (int j = 0; j < n_Electrons; j++)
                        {
                            val1 = EqualFunc((a+n_Electrons)%2, (b+n_Electrons)%2) * EqualFunc(i%2, j%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) ((b+n_Electrons)/2), (int) (i/2), (int) (j/2));

                            val2 = EqualFunc((a+n_Electrons)%2, j%2) * EqualFunc(i%2, (b+n_Electrons)%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) (j/2), (int) (i/2), (int) ((b+n_Electrons)/2));

                            Rear_cube7(a,i)(B+Speed_Occ,J) = val1-val2;
                            j++;
                            J++;
                        }
                    }

                    else
                    {
                        J = 0;
                        for (int j = 1; j < n_Electrons; j++)
                        {
                            val1 = EqualFunc((a+n_Electrons)%2, (b+n_Electrons)%2) * EqualFunc(i%2, j%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) ((b+n_Electrons)/2), (int) (i/2), (int) (j/2));

                            val2 = EqualFunc((a+n_Electrons)%2, j%2) * EqualFunc(i%2, (b+n_Electrons)%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) (j/2), (int) (i/2), (int) ((b+n_Electrons)/2));

                            Rear_cube7(a,i)(B+Speed_Occ,J) = val1-val2;
                            j++;
                            J++;
                        }
                    }
                }
                b++;
                B++;
            }
        }

        // Compact form of I(i,j,k,l), dont store any zeroes
        Rear_cube8.set_size(n_Electrons, n_Electrons);
        field_filler_matrix2 = zeros(n_Electrons, n_Electrons/2);
        for (int i = 0; i < n_Electrons; i++)
        {
            for (int a = i+1; a < n_Electrons; a++)
            {
                Rear_cube8(i,a) = field_filler_matrix2;
            }
        }

        for (int a = 0; a < n_Electrons; a++)
        {
            B = 0;
            for (int b = 0; b < n_Electrons; b++)
            {
                for (int i = a+1; i < n_Electrons; i++)
                {
                    if ((a+b+i)%2 == 0)
                    {
                        J = 0;
                        for (int j = 0; j < n_Electrons; j++)
                        {
                            val1 = EqualFunc(a%2, b%2) * EqualFunc(i%2, j%2) *
                                    temp4.at((int) (a/2))((int) (b/2), (int) (i/2), (int) (j/2));

                            val2 = EqualFunc(a%2, j%2) * EqualFunc(i%2, b%2) *
                                    temp4.at((int) (a/2))((int) (j/2), (int) (i/2), (int) (b/2));

                            Rear_cube8(a,i)(B,J) = val1-val2;
                            j++;
                            J++;
                        }
                    }

                    else
                    {
                        J = 0;
                        for (int j = 1; j < n_Electrons; j++)
                        {
                            val1 = EqualFunc(a%2, b%2) * EqualFunc(i%2, j%2) *
                                    temp4.at((int) (a/2))((int) (b/2), (int) (i/2), (int) (j/2));

                            val2 = EqualFunc(a%2, j%2) * EqualFunc(i%2, b%2) *
                                    temp4.at((int) (a/2))((int) (j/2), (int) (i/2), (int) (b/2));

                            Rear_cube8(a,i)(B,J) = val1-val2;
                            j++;
                            J++;
                        }
                    }
                }
                b++;
                B++;
            }

            B = 0;
            for (int b = 1; b < n_Electrons; b++)
            {
                for (int i = a+1; i < n_Electrons; i++)
                {
                    if ((a+b+i)%2 == 0)
                    {
                        J = 0;
                        for (int j = 0; j < n_Electrons; j++)
                        {
                            val1 = EqualFunc(a%2, b%2) * EqualFunc(i%2, j%2) *
                                    temp4.at((int) (a/2))((int) (b/2), (int) (i/2), (int) (j/2));

                            val2 = EqualFunc(a%2, j%2) * EqualFunc(i%2, b%2) *
                                    temp4.at((int) (a/2))((int) (j/2), (int) (i/2), (int) (b/2));

                            Rear_cube8(a,i)(B+n_Electrons/2,J) = val1-val2;
                            j++;
                            J++;
                        }
                    }

                    else
                    {
                        J = 0;
                        for (int j = 1; j < n_Electrons; j++)
                        {
                            val1 = EqualFunc(a%2, b%2) * EqualFunc(i%2, j%2) *
                                    temp4.at((int) (a/2))((int) (b/2), (int) (i/2), (int) (j/2));

                            val2 = EqualFunc(a%2, j%2) * EqualFunc(i%2, b%2) *
                                    temp4.at((int) (a/2))((int) (j/2), (int) (i/2), (int) (b/2));

                            Rear_cube8(a,i)(B+n_Electrons/2,J) = val1-val2;
                            j++;
                            J++;
                        }
                    }
                }
                b++;
                B++;
            }
        }

        // Compact form of I(a,b,c,i), dont store any zeroes
        Rear_cube9.set_size(unocc_orb, unocc_orb);
        field_filler_matrix2 = zeros(unocc_orb, n_Electrons/2);
        for (int i = 0; i < unocc_orb; i++)
        {
            for (int a = i+1; a < unocc_orb; a++)
            {
                Rear_cube9(i,a) = field_filler_matrix2;
            }
        }

        for (int a = 0; a < unocc_orb; a++)
        {
            B = 0;
            for (int b = 0; b < unocc_orb; b++)
            {
                for (int i = a+1; i < unocc_orb; i++)
                {
                    if ((a+b+i)%2 == 0)
                    {
                        J = 0;
                        for (int j = 0; j < n_Electrons; j++)
                        {
                            val1 = EqualFunc((a+n_Electrons)%2, (b+n_Electrons)%2) * EqualFunc((i+n_Electrons)%2, j%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) ((b+n_Electrons)/2), (int) ((i+n_Electrons)/2), (int) (j/2));

                            val2 = EqualFunc((a+n_Electrons)%2, j%2) * EqualFunc((i+n_Electrons)%2, (b+n_Electrons)%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) (j/2), (int) ((i+n_Electrons)/2), (int) ((b+n_Electrons)/2));

                            Rear_cube9(a,i)(B,J) = val1-val2;
                            j++;
                            J++;
                        }
                    }

                    else
                    {
                        J = 0;
                        for (int j = 1; j < n_Electrons; j++)
                        {
                            val1 = EqualFunc((a+n_Electrons)%2, (b+n_Electrons)%2) * EqualFunc((i+n_Electrons)%2, j%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) ((b+n_Electrons)/2), (int) ((i+n_Electrons)/2), (int) (j/2));

                            val2 = EqualFunc((a+n_Electrons)%2, j%2) * EqualFunc((i+n_Electrons)%2, (b+n_Electrons)%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) (j/2), (int) ((i+n_Electrons)/2), (int) ((b+n_Electrons)/2));

                            Rear_cube9(a,i)(B,J) = val1-val2;
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
                for (int i = a+1; i < unocc_orb; i++)
                {
                    if ((a+b+i)%2 == 0)
                    {
                        J = 0;
                        for (int j = 0; j < n_Electrons; j++)
                        {
                            val1 = EqualFunc((a+n_Electrons)%2, (b+n_Electrons)%2) * EqualFunc((i+n_Electrons)%2, j%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) ((b+n_Electrons)/2), (int) ((i+n_Electrons)/2), (int) (j/2));

                            val2 = EqualFunc((a+n_Electrons)%2, j%2) * EqualFunc((i+n_Electrons)%2, (b+n_Electrons)%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) (j/2), (int) ((i+n_Electrons)/2), (int) ((b+n_Electrons)/2));

                            Rear_cube9(a,i)(B+Speed_Occ,J) = val1-val2;
                            j++;
                            J++;
                        }
                    }

                    else
                    {
                        J = 0;
                        for (int j = 1; j < n_Electrons; j++)
                        {
                            val1 = EqualFunc((a+n_Electrons)%2, (b+n_Electrons)%2) * EqualFunc((i+n_Electrons)%2, j%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) ((b+n_Electrons)/2), (int) ((i+n_Electrons)/2), (int) (j/2));

                            val2 = EqualFunc((a+n_Electrons)%2, j%2) * EqualFunc((i+n_Electrons)%2, (b+n_Electrons)%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) (j/2), (int) ((i+n_Electrons)/2), (int) ((b+n_Electrons)/2));

                            Rear_cube9(a,i)(B+Speed_Occ,J) = val1-val2;
                            j++;
                            J++;
                        }
                    }
                }
                b++;
                B++;
            }
        }

        // Compact form of I(a,b,c,i), dont store any zeroes
        Rear_cube10.set_size(unocc_orb, n_Electrons);
        field_filler_matrix2 = zeros(unocc_orb, n_Electrons/2);
        for (int i = 0; i < n_Electrons; i++)
        {
            for (int a = 0; a < unocc_orb; a++)
            {
                Rear_cube10(a,i) = field_filler_matrix2;
            }
        }

        for (int a = 0; a < unocc_orb; a++)
        {
            B = 0;
            for (int b = 0; b < unocc_orb; b++)
            {
                for (int i = 0; i < n_Electrons; i++)
                {
                    if ((a+b+i)%2 == 0)
                    {
                        J = 0;
                        for (int j = 0; j < n_Electrons; j++)
                        {
                            val1 = EqualFunc((a+n_Electrons)%2, (b+n_Electrons)%2) * EqualFunc(j%2, i%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) ((b+n_Electrons)/2), (int) (j/2), (int) (i/2));

                            val2 = EqualFunc((a+n_Electrons)%2, i%2) * EqualFunc(j%2, (b+n_Electrons)%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) (i/2), (int) (j/2), (int) ((b+n_Electrons)/2));

                            Rear_cube10(a,i)(B,J) = val1-val2;
                            j++;
                            J++;
                        }
                    }

                    else
                    {
                        J = 0;
                        for (int j = 1; j < n_Electrons; j++)
                        {
                            val1 = EqualFunc((a+n_Electrons)%2, (b+n_Electrons)%2) * EqualFunc(j%2, i%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) ((b+n_Electrons)/2), (int) (j/2), (int) (i/2));

                            val2 = EqualFunc((a+n_Electrons)%2, i%2) * EqualFunc(j%2, (b+n_Electrons)%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) (i/2), (int) (j/2), (int) ((b+n_Electrons)/2));

                            Rear_cube10(a,i)(B,J) = val1-val2;
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
                    if ((a+b+i)%2 == 0)
                    {
                        J = 0;
                        for (int j = 0; j < n_Electrons; j++)
                        {
                            val1 = EqualFunc((a+n_Electrons)%2, (b+n_Electrons)%2) * EqualFunc(j%2, i%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) ((b+n_Electrons)/2), (int) (j/2), (int) (i/2));

                            val2 = EqualFunc((a+n_Electrons)%2, i%2) * EqualFunc(j%2, (b+n_Electrons)%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) (i/2), (int) (j/2), (int) ((b+n_Electrons)/2));

                            Rear_cube10(a,i)(B+Speed_Occ,J) = val1-val2;
                            j++;
                            J++;
                        }
                    }

                    else
                    {
                        J = 0;
                        for (int j = 1; j < n_Electrons; j++)
                        {
                            val1 = EqualFunc((a+n_Electrons)%2, (b+n_Electrons)%2) * EqualFunc(j%2, i%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) ((b+n_Electrons)/2), (int) (j/2), (int) (i/2));

                            val2 = EqualFunc((a+n_Electrons)%2, i%2) * EqualFunc(j%2, (b+n_Electrons)%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) (i/2), (int) (j/2), (int) ((b+n_Electrons)/2));

                            Rear_cube10(a,i)(B+Speed_Occ,J) = val1-val2;
                            j++;
                            J++;
                        }
                    }
                }
                b++;
                B++;
            }
        }

        temp4.clear();

}

void coupled_cluster_integrals::Prepear_Integralsv2()
{
    cube temp_cube = zeros(Matrix_Size, Matrix_Size, Matrix_Size);
    mat field_filler_matrix2;

    vector <cube> temp;
    vector <cube> temp2;
    vector <cube> temp3;

    for (int i = 0; i < Matrix_Size; i++)
    {
        temp.push_back(temp_cube);
    }

    for (int a = 0; a < Matrix_Size; a++)
    {
        for (int i = 0; i < Matrix_Size; i++)
        {
            for (int j = 0; j < Matrix_Size; j++)
            {
                for (int k = 0; k < Matrix_Size; k++)
                {
                    for (int l = 0; l < Matrix_Size; l++)
                    {
                        temp.at(a)(j,k,l) += c(i,a) * Indexed_integrals(Return_Integral_Index(i,j,k,l));//Integrals.at(i)(k,j,l) * c(i,a);
                    }
                }
            }
        }
    }
    Integrals.clear();

    for (int i = 0; i < Matrix_Size; i++)
    {
        temp2.push_back(temp_cube);
    }

    for (int a = 0; a < Matrix_Size; a++)
    {
        for (int b = 0; b < Matrix_Size; b++)
        {
            for (int i = 0; i < Matrix_Size; i++)
            {
                for (int j = 0; j < Matrix_Size; j++)
                {
                    for (int k = 0; k < Matrix_Size; k++)
                    {
                        temp2.at(a)(b,j,k) += c(i,b) * temp.at(a)(i,j,k);
                    }
                }
            }
        }
    }

    temp.clear();

    for (int i = 0; i < Matrix_Size; i++)
    {
        temp3.push_back(temp_cube);
    }

    for (int a = 0; a < Matrix_Size; a++)
    {
        for (int b = 0; b < Matrix_Size; b++)
        {
            for (int g = 0; g < Matrix_Size; g++)
            {
                for (int i = 0; i < Matrix_Size; i++)
                {
                    for (int j = 0; j < Matrix_Size; j++)
                    {
                        temp3.at(a)(b,g,j) += c(i,g) * temp2.at(a)(b,i,j);
                    }
                }
            }
        }
    }

    temp2.clear();

    for (int i = 0; i < Matrix_Size; i++)
    {
        temp4.push_back(temp_cube);
    }

    for (int a = 0; a < Matrix_Size; a++)
    {
        for (int b = 0; b < Matrix_Size; b++)
        {
            for (int g = 0; g < Matrix_Size; g++)
            {
                for (int h = 0; h < Matrix_Size; h++)
                {
                    for (int i = 0; i < Matrix_Size; i++)
                    {
                        temp4.at(a)(b,g,h) += c(i,h) * temp3.at(a)(b,g,i);
                    }
                }
            }
        }
    }

    temp3.clear();

    double val1, val2;

    int unocc_orb = 2*Matrix_Size - n_Electrons;
    int Speed_Occ = unocc_orb/2;

        Rear_cube.set_size(unocc_orb, n_Electrons);
        field_filler_matrix2 = zeros(unocc_orb, n_Electrons);
        for (int i = 0; i < n_Electrons; i++)
        {
            for (int a = 0; a < unocc_orb; a++)
            {
                Rear_cube(a,i) = field_filler_matrix2;
            }
        }

        // Using properties of the integrals we can store only half the elements


        // Compact form of I(a,b,i,j), dont store any zeroes
        int B, J;
        for (int a = 0; a < unocc_orb; a++)
        {
            for (int b = 0; b < unocc_orb; b++)
            {
                for (int i = 0; i < n_Electrons; i++)
                {
                    if ((a-b+i)%2 == 0)
                    {
                        for (int j = 0; j < n_Electrons; j++)
                        {
                            val1 = EqualFunc((a+n_Electrons)%2, i%2) * EqualFunc((b+n_Electrons)%2, j%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) (i/2), (int) ((b+n_Electrons)/2), (int) (j/2));

                            val2 = EqualFunc((a+n_Electrons)%2, j%2) * EqualFunc((b+n_Electrons)%2, i%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) (j/2), (int) ((b+n_Electrons)/2), (int) (i/2));

                            Rear_cube(a,i)(b,j) = val1-val2;
                            j++;
                        }
                    }

                    else
                    {
                        for (int j = 1; j < n_Electrons; j++)
                        {
                            val1 = EqualFunc((a+n_Electrons)%2, i%2) * EqualFunc((b+n_Electrons)%2, j%2) *
                                     temp4.at((int) ((a+n_Electrons)/2))((int) (i/2), (int) ((b+n_Electrons)/2), (int) (j/2));


                            val2 = EqualFunc((a+n_Electrons)%2, j%2) * EqualFunc((b+n_Electrons)%2, i%2) *
                                     temp4.at((int) ((a+n_Electrons)/2))((int) (j/2), (int) ((b+n_Electrons)/2), (int) (i/2));

                             Rear_cube(a,i)(b,j) = val1-val2;
                            j++;
                        }
                    }
                }
                b++;
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
                            val1 = EqualFunc((a+n_Electrons)%2, i%2) * EqualFunc((b+n_Electrons)%2, j%2) *
                                     temp4.at((int) ((a+n_Electrons)/2))((int) (i/2), (int) ((b+n_Electrons)/2), (int) (j/2));

                            val2 = EqualFunc((a+n_Electrons)%2, j%2) * EqualFunc((b+n_Electrons)%2, i%2) *
                                     temp4.at((int) ((a+n_Electrons)/2))((int) (j/2), (int) ((b+n_Electrons)/2), (int) (i/2));

                            Rear_cube(a,i)(b,j) = val1-val2;
                            j++;
                        }
                    }

                    else
                    {
                        J = 0;
                        for (int j = 1; j < n_Electrons; j++)
                        {
                            val1 = EqualFunc((a+n_Electrons)%2, i%2) * EqualFunc((b+n_Electrons)%2, j%2) *
                                     temp4.at((int) ((a+n_Electrons)/2))((int) (i/2), (int) ((b+n_Electrons)/2), (int) (j/2));

                            val2 = EqualFunc((a+n_Electrons)%2, j%2) * EqualFunc((b+n_Electrons)%2, i%2) *
                                     temp4.at((int) ((a+n_Electrons)/2))((int) (j/2), (int) ((b+n_Electrons)/2), (int) (i/2));

                            Rear_cube(a,i)(b,j) = val1-val2;
                            j++;
                        }
                    }
                }
                b++;
            }
        }


        // Compact form of I(a,b,c,d), dont store any zeroes
        Rear_cube2.set_size(unocc_orb, unocc_orb);
        field_filler_matrix2 = zeros(unocc_orb, unocc_orb);
        for (int i = 0; i < unocc_orb; i++)
        {
            for (int a = i+1; a < unocc_orb; a++)
            {
                Rear_cube2(i,a) = field_filler_matrix2;
            }
        }

        for (int a = 0; a < unocc_orb; a++)
        {
            for (int b = 0; b < unocc_orb; b++)
            {
                for (int i = a+1; i < unocc_orb; i++)
                {
                    if ((a+b+i)%2 == 0)
                    {
                        for (int j = 0; j < unocc_orb; j++)
                        {
                            val1 = EqualFunc((a+n_Electrons)%2, (b+n_Electrons)%2) * EqualFunc((i+n_Electrons)%2, (j+n_Electrons)%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) ((b+n_Electrons)/2), (int) ((i+n_Electrons)/2), (int) ((j+n_Electrons)/2));

                            val2 = EqualFunc((a+n_Electrons)%2, (j+n_Electrons)%2) * EqualFunc((i+n_Electrons)%2, (b+n_Electrons)%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) ((j+n_Electrons)/2), (int) ((i+n_Electrons)/2), (int) ((b+n_Electrons)/2));

                            Rear_cube2(a,i)(b,j) = val1 - val2;
                            j++;
                        }
                    }

                    else
                    {
                        J = 0;
                        for (int j = 1; j < unocc_orb; j++)
                        {
                            val1 = EqualFunc((a+n_Electrons)%2, (b+n_Electrons)%2) * EqualFunc((i+n_Electrons)%2, (j+n_Electrons)%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) ((b+n_Electrons)/2), (int) ((i+n_Electrons)/2), (int) ((j+n_Electrons)/2));

                            val2 = EqualFunc((a+n_Electrons)%2, (j+n_Electrons)%2) * EqualFunc((i+n_Electrons)%2, (b+n_Electrons)%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) ((j+n_Electrons)/2), (int) ((i+n_Electrons)/2), (int) ((b+n_Electrons)/2));

                            Rear_cube2(a,i)(b,j) = val1 - val2;
                            j++;
                        }
                    }
                }
                b++;
            }

            for (int b = 1; b < unocc_orb; b++)
            {
                for (int i = a+1; i < unocc_orb; i++)
                {
                    if ((a+b+i)%2 == 0)
                    {
                        for (int j = 0; j < unocc_orb; j++)
                        {
                            val1 = EqualFunc((a+n_Electrons)%2, (b+n_Electrons)%2) * EqualFunc((i+n_Electrons)%2, (j+n_Electrons)%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) ((b+n_Electrons)/2), (int) ((i+n_Electrons)/2), (int) ((j+n_Electrons)/2));

                            val2 = EqualFunc((a+n_Electrons)%2, (j+n_Electrons)%2) * EqualFunc((i+n_Electrons)%2, (b+n_Electrons)%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) ((j+n_Electrons)/2), (int) ((i+n_Electrons)/2), (int) ((b+n_Electrons)/2));

                            Rear_cube2(a,i)(b,j) = val1 - val2;
                            j++;
                        }
                    }

                    else
                    {
                        for (int j = 1; j < unocc_orb; j++)
                        {
                            val1 = EqualFunc((a+n_Electrons)%2, (b+n_Electrons)%2) * EqualFunc((i+n_Electrons)%2, (j+n_Electrons)%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) ((b+n_Electrons)/2), (int) ((i+n_Electrons)/2), (int) ((j+n_Electrons)/2));

                            val2 = EqualFunc((a+n_Electrons)%2, (j+n_Electrons)%2) * EqualFunc((i+n_Electrons)%2, (b+n_Electrons)%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) ((j+n_Electrons)/2), (int) ((i+n_Electrons)/2), (int) ((b+n_Electrons)/2));

                            Rear_cube2(a,i)(b,j) = val1 - val2;
                            j++;
                        }
                    }
                }
                b++;
            }
        }

        // Compact form of I(i,j,a,b), dont store any zeroes
        Rear_cube3.set_size(n_Electrons, n_Electrons);
        field_filler_matrix2 = zeros(unocc_orb, unocc_orb);
        for (int i = 0; i < n_Electrons; i++)
        {
            for (int a = 0; a < n_Electrons; a++)
            {
                Rear_cube3(a,i) = field_filler_matrix2;
            }
        }

        for (int a = 0; a < n_Electrons; a++)
        {
            for (int b = 0; b < unocc_orb; b++)
            {
                for (int i = 0; i < n_Electrons; i++)
                {
                    if ((a+b+i)%2 == 0)
                    {
                        for (int j = 0; j < unocc_orb; j++)
                        {
                            val1 = EqualFunc(a%2, (b+n_Electrons)%2) * EqualFunc(i%2, (j+n_Electrons)%2) *
                                    temp4.at((int) (a/2))((int) ((b+n_Electrons)/2), (int) (i/2), (int) ((j+n_Electrons)/2));

                            val2 = EqualFunc(a%2, (j+n_Electrons)%2) * EqualFunc(i%2, (b+n_Electrons)%2) *
                                    temp4.at((int) (a/2))((int) ((j+n_Electrons)/2), (int) (i/2), (int) ((b+n_Electrons)/2));

                            Rear_cube3(a,i)(b,j) = val1-val2;
                            j++;
                        }
                    }

                    else
                    {
                        J = 0;
                        for (int j = 1; j < unocc_orb; j++)
                        {
                            val1 = EqualFunc(a%2, (b+n_Electrons)%2) * EqualFunc(i%2, (j+n_Electrons)%2) *
                                    temp4.at((int) (a/2))((int) ((b+n_Electrons)/2), (int) (i/2), (int) ((j+n_Electrons)/2));

                            val2 = EqualFunc(a%2, (j+n_Electrons)%2) * EqualFunc(i%2, (b+n_Electrons)%2) *
                                    temp4.at((int) (a/2))((int) ((j+n_Electrons)/2), (int) (i/2), (int) ((b+n_Electrons)/2));

                            Rear_cube3(a,i)(b,j) = val1-val2;
                            j++;
                        }
                    }
                }
                b++;
            }

            for (int b = 1; b < unocc_orb; b++)
            {
                for (int i = 0; i < n_Electrons; i++)
                {
                    if ((a+b+i)%2 == 0)
                    {
                        for (int j = 0; j < unocc_orb; j++)
                        {
                            val1 = EqualFunc(a%2, (b+n_Electrons)%2) * EqualFunc(i%2, (j+n_Electrons)%2) *
                                    temp4.at((int) (a/2))((int) ((b+n_Electrons)/2), (int) (i/2), (int) ((j+n_Electrons)/2));

                            val2 = EqualFunc(a%2, (j+n_Electrons)%2) * EqualFunc(i%2, (b+n_Electrons)%2) *
                                    temp4.at((int) (a/2))((int) ((j+n_Electrons)/2), (int) (i/2), (int) ((b+n_Electrons)/2));

                            Rear_cube3(a,i)(b,j) = val1-val2;
                            j++;
                        }
                    }

                    else
                    {
                        J = 0;
                        for (int j = 1; j < unocc_orb; j++)
                        {
                            val1 = EqualFunc(a%2, (b+n_Electrons)%2) * EqualFunc(i%2, (j+n_Electrons)%2) *
                                    temp4.at((int) (a/2))((int) ((b+n_Electrons)/2), (int) (i/2), (int) ((j+n_Electrons)/2));

                            val2 = EqualFunc(a%2, (j+n_Electrons)%2) * EqualFunc(i%2, (b+n_Electrons)%2) *
                                    temp4.at((int) (a/2))((int) ((j+n_Electrons)/2), (int) (i/2), (int) ((b+n_Electrons)/2));

                            Rear_cube3(a,i)(b,j) = val1-val2;
                            j++;
                        }
                    }
                }
                b++;
            }
        }


        // Compact form of I(a,i,b,c), dont store any zeroes
        Rear_cube4.set_size(unocc_orb, n_Electrons);
        field_filler_matrix2 = zeros(unocc_orb, unocc_orb);
        for (int i = 0; i < n_Electrons; i++)
        {
            for (int a = 0; a < unocc_orb; a++)
            {
                Rear_cube4(a,i) = field_filler_matrix2;
            }
        }

        for (int a = 0; a < unocc_orb; a++)
        {
            B = 0;
            for (int b = 0; b < unocc_orb; b++)
            {
                for (int i = 0; i < n_Electrons; i++)
                {
                    if ((a+b+i)%2 == 0)
                    {
                        J = 0;
                        for (int j = 0; j < unocc_orb; j++)
                        {
                            val1 = EqualFunc((a+n_Electrons)%2, (b+n_Electrons)%2) * EqualFunc(i%2, (j+n_Electrons)%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) ((b+n_Electrons)/2), (int) (i/2), (int) ((j+n_Electrons)/2));

                            val2 = EqualFunc((a+n_Electrons)%2, (j+n_Electrons)%2) * EqualFunc(i%2, (b+n_Electrons)%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) ((j+n_Electrons)/2), (int) (i/2), (int) ((b+n_Electrons)/2));

                            Rear_cube4(a,i)(b,j) = val1-val2;
                            j++;
                            J++;
                        }
                    }

                    else
                    {
                        J = 0;
                        for (int j = 1; j < unocc_orb; j++)
                        {
                            val1 = EqualFunc((a+n_Electrons)%2, (b+n_Electrons)%2) * EqualFunc(i%2, (j+n_Electrons)%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) ((b+n_Electrons)/2), (int) (i/2), (int) ((j+n_Electrons)/2));

                            val2 = EqualFunc((a+n_Electrons)%2, (j+n_Electrons)%2) * EqualFunc(i%2, (b+n_Electrons)%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) ((j+n_Electrons)/2), (int) (i/2), (int) ((b+n_Electrons)/2));

                            Rear_cube4(a,i)(b,j) = val1-val2;
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
                    if ((a+b+i)%2 == 0)
                    {
                        J = 0;
                        for (int j = 0; j < unocc_orb; j++)
                        {
                            val1 = EqualFunc((a+n_Electrons)%2, (b+n_Electrons)%2) * EqualFunc(i%2, (j+n_Electrons)%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) ((b+n_Electrons)/2), (int) (i/2), (int) ((j+n_Electrons)/2));

                            val2 = EqualFunc((a+n_Electrons)%2, (j+n_Electrons)%2) * EqualFunc(i%2, (b+n_Electrons)%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) ((j+n_Electrons)/2), (int) (i/2), (int) ((b+n_Electrons)/2));

                            Rear_cube4(a,i)(b,j) = val1-val2;
                            j++;
                            J++;
                        }
                    }

                    else
                    {
                        J = 0;
                        for (int j = 1; j < unocc_orb; j++)
                        {
                            val1 = EqualFunc((a+n_Electrons)%2, (b+n_Electrons)%2) * EqualFunc(i%2, (j+n_Electrons)%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) ((b+n_Electrons)/2), (int) (i/2), (int) ((j+n_Electrons)/2));

                            val2 = EqualFunc((a+n_Electrons)%2, (j+n_Electrons)%2) * EqualFunc(i%2, (b+n_Electrons)%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) ((j+n_Electrons)/2), (int) (i/2), (int) ((b+n_Electrons)/2));

                            Rear_cube4(a,i)(b,j) = val1-val2;
                            j++;
                            J++;
                        }
                    }
                }
                b++;
                B++;
            }
        }

        // Compact form of I(i,j,a,k), dont store any zeroes
        Rear_cube6.set_size(n_Electrons, n_Electrons);
        field_filler_matrix2 = zeros(unocc_orb, n_Electrons);
        for (int i = 0; i < n_Electrons; i++)
        {
            for (int a = 0; a < n_Electrons; a++)
            {
                Rear_cube6(a,i) = field_filler_matrix2;
            }
        }

        for (int a = 0; a < n_Electrons; a++)
        {
            B = 0;
            for (int b = 0; b < unocc_orb; b++)
            {
                for (int i = 0; i < n_Electrons; i++)
                {
                    if ((a+b+i)%2 == 0)
                    {
                        J = 0;
                        for (int j = 0; j < n_Electrons; j++)
                        {
                            val1 = EqualFunc(a%2, (b+n_Electrons)%2) * EqualFunc(i%2, j%2) *
                                    temp4.at((int) (a/2))((int) ((b+n_Electrons)/2), (int) (i/2), (int) (j/2));

                            val2 = EqualFunc(a%2, j%2) * EqualFunc(i%2, (b+n_Electrons)%2) *
                                    temp4.at((int) (a/2))((int) (j/2), (int) (i/2), (int) ((b+n_Electrons)/2));

                            Rear_cube6(a,i)(b,j) = val1-val2;
                            j++;
                            J++;
                        }
                    }

                    else
                    {
                        J = 0;
                        for (int j = 1; j < n_Electrons; j++)
                        {
                            val1 = EqualFunc(a%2, (b+n_Electrons)%2) * EqualFunc(i%2, j%2) *
                                    temp4.at((int) (a/2))((int) ((b+n_Electrons)/2), (int) (i/2), (int) (j/2));

                            val2 = EqualFunc(a%2, j%2) * EqualFunc(i%2, (b+n_Electrons)%2) *
                                    temp4.at((int) (a/2))((int) (j/2), (int) (i/2), (int) ((b+n_Electrons)/2));

                            Rear_cube6(a,i)(b,j) = val1-val2;
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
                    if ((a+b+i)%2 == 0)
                    {
                        J = 0;
                        for (int j = 0; j < n_Electrons; j++)
                        {
                            val1 = EqualFunc(a%2, (b+n_Electrons)%2) * EqualFunc(i%2, j%2) *
                                    temp4.at((int) (a/2))((int) ((b+n_Electrons)/2), (int) (i/2), (int) (j/2));

                            val2 = EqualFunc(a%2, j%2) * EqualFunc(i%2, (b+n_Electrons)%2) *
                                    temp4.at((int) (a/2))((int) (j/2), (int) (i/2), (int) ((b+n_Electrons)/2));

                            Rear_cube6(a,i)(b,j) = val1-val2;
                            j++;
                            J++;
                        }
                    }

                    else
                    {
                        J = 0;
                        for (int j = 1; j < n_Electrons; j++)
                        {
                            val1 = EqualFunc(a%2, (b+n_Electrons)%2) * EqualFunc(i%2, j%2) *
                                    temp4.at((int) (a/2))((int) ((b+n_Electrons)/2), (int) (i/2), (int) (j/2));

                            val2 = EqualFunc(a%2, j%2) * EqualFunc(i%2, (b+n_Electrons)%2) *
                                    temp4.at((int) (a/2))((int) (j/2), (int) (i/2), (int) ((b+n_Electrons)/2));

                            Rear_cube6(a,i)(b,j) = val1-val2;
                            j++;
                            J++;
                        }
                    }
                }
                b++;
                B++;
            }
        }

        // Compact form of I(a,i,b,j), dont store any zeroes
        Rear_cube7.set_size(unocc_orb, n_Electrons);
        field_filler_matrix2 = zeros(unocc_orb, n_Electrons);
        for (int i = 0; i < n_Electrons; i++)
        {
            for (int a = 0; a < unocc_orb; a++)
            {
                Rear_cube7(a,i) = field_filler_matrix2;
            }
        }

        for (int a = 0; a < unocc_orb; a++)
        {
            B = 0;
            for (int b = 0; b < unocc_orb; b++)
            {
                for (int i = 0; i < n_Electrons; i++)
                {
                    if ((a+b+i)%2 == 0)
                    {
                        J = 0;
                        for (int j = 0; j < n_Electrons; j++)
                        {
                            val1 = EqualFunc((a+n_Electrons)%2, (b+n_Electrons)%2) * EqualFunc(i%2, j%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) ((b+n_Electrons)/2), (int) (i/2), (int) (j/2));

                            val2 = EqualFunc((a+n_Electrons)%2, j%2) * EqualFunc(i%2, (b+n_Electrons)%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) (j/2), (int) (i/2), (int) ((b+n_Electrons)/2));

                            Rear_cube7(a,i)(b,j) = val1-val2;
                            j++;
                            J++;
                        }
                    }

                    else
                    {
                        J = 0;
                        for (int j = 1; j < n_Electrons; j++)
                        {
                            val1 = EqualFunc((a+n_Electrons)%2, (b+n_Electrons)%2) * EqualFunc(i%2, j%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) ((b+n_Electrons)/2), (int) (i/2), (int) (j/2));

                            val2 = EqualFunc((a+n_Electrons)%2, j%2) * EqualFunc(i%2, (b+n_Electrons)%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) (j/2), (int) (i/2), (int) ((b+n_Electrons)/2));

                            Rear_cube7(a,i)(b,j) = val1-val2;
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
                    if ((a+b+i)%2 == 0)
                    {
                        J = 0;
                        for (int j = 0; j < n_Electrons; j++)
                        {
                            val1 = EqualFunc((a+n_Electrons)%2, (b+n_Electrons)%2) * EqualFunc(i%2, j%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) ((b+n_Electrons)/2), (int) (i/2), (int) (j/2));

                            val2 = EqualFunc((a+n_Electrons)%2, j%2) * EqualFunc(i%2, (b+n_Electrons)%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) (j/2), (int) (i/2), (int) ((b+n_Electrons)/2));

                            Rear_cube7(a,i)(b,j) = val1-val2;
                            j++;
                            J++;
                        }
                    }

                    else
                    {
                        J = 0;
                        for (int j = 1; j < n_Electrons; j++)
                        {
                            val1 = EqualFunc((a+n_Electrons)%2, (b+n_Electrons)%2) * EqualFunc(i%2, j%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) ((b+n_Electrons)/2), (int) (i/2), (int) (j/2));

                            val2 = EqualFunc((a+n_Electrons)%2, j%2) * EqualFunc(i%2, (b+n_Electrons)%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) (j/2), (int) (i/2), (int) ((b+n_Electrons)/2));

                            Rear_cube7(a,i)(b,j) = val1-val2;
                            j++;
                            J++;
                        }
                    }
                }
                b++;
                B++;
            }
        }

        // Compact form of I(i,j,k,l), dont store any zeroes
        Rear_cube8.set_size(n_Electrons, n_Electrons);
        field_filler_matrix2 = zeros(n_Electrons, n_Electrons);
        for (int i = 0; i < n_Electrons; i++)
        {
            for (int a = i+1; a < n_Electrons; a++)
            {
                Rear_cube8(i,a) = field_filler_matrix2;
            }
        }

        for (int a = 0; a < n_Electrons; a++)
        {
            B = 0;
            for (int b = 0; b < n_Electrons; b++)
            {
                for (int i = a+1; i < n_Electrons; i++)
                {
                    if ((a+b+i)%2 == 0)
                    {
                        J = 0;
                        for (int j = 0; j < n_Electrons; j++)
                        {
                            val1 = EqualFunc(a%2, b%2) * EqualFunc(i%2, j%2) *
                                    temp4.at((int) (a/2))((int) (b/2), (int) (i/2), (int) (j/2));

                            val2 = EqualFunc(a%2, j%2) * EqualFunc(i%2, b%2) *
                                    temp4.at((int) (a/2))((int) (j/2), (int) (i/2), (int) (b/2));

                            Rear_cube8(a,i)(b,j) = val1-val2;
                            j++;
                            J++;
                        }
                    }

                    else
                    {
                        J = 0;
                        for (int j = 1; j < n_Electrons; j++)
                        {
                            val1 = EqualFunc(a%2, b%2) * EqualFunc(i%2, j%2) *
                                    temp4.at((int) (a/2))((int) (b/2), (int) (i/2), (int) (j/2));

                            val2 = EqualFunc(a%2, j%2) * EqualFunc(i%2, b%2) *
                                    temp4.at((int) (a/2))((int) (j/2), (int) (i/2), (int) (b/2));

                            Rear_cube8(a,i)(b,j) = val1-val2;
                            j++;
                            J++;
                        }
                    }
                }
                b++;
                B++;
            }

            B = 0;
            for (int b = 1; b < n_Electrons; b++)
            {
                for (int i = a+1; i < n_Electrons; i++)
                {
                    if ((a+b+i)%2 == 0)
                    {
                        J = 0;
                        for (int j = 0; j < n_Electrons; j++)
                        {
                            val1 = EqualFunc(a%2, b%2) * EqualFunc(i%2, j%2) *
                                    temp4.at((int) (a/2))((int) (b/2), (int) (i/2), (int) (j/2));

                            val2 = EqualFunc(a%2, j%2) * EqualFunc(i%2, b%2) *
                                    temp4.at((int) (a/2))((int) (j/2), (int) (i/2), (int) (b/2));

                            Rear_cube8(a,i)(b,j) = val1-val2;
                            j++;
                            J++;
                        }
                    }

                    else
                    {
                        J = 0;
                        for (int j = 1; j < n_Electrons; j++)
                        {
                            val1 = EqualFunc(a%2, b%2) * EqualFunc(i%2, j%2) *
                                    temp4.at((int) (a/2))((int) (b/2), (int) (i/2), (int) (j/2));

                            val2 = EqualFunc(a%2, j%2) * EqualFunc(i%2, b%2) *
                                    temp4.at((int) (a/2))((int) (j/2), (int) (i/2), (int) (b/2));

                            Rear_cube8(a,i)(b,j) = val1-val2;
                            j++;
                            J++;
                        }
                    }
                }
                b++;
                B++;
            }
        }

        // Compact form of I(a,b,c,i), dont store any zeroes
        Rear_cube9.set_size(unocc_orb, unocc_orb);
        field_filler_matrix2 = zeros(unocc_orb, n_Electrons);
        for (int i = 0; i < unocc_orb; i++)
        {
            for (int a = i+1; a < unocc_orb; a++)
            {
                Rear_cube9(i,a) = field_filler_matrix2;
            }
        }

        for (int a = 0; a < unocc_orb; a++)
        {
            B = 0;
            for (int b = 0; b < unocc_orb; b++)
            {
                for (int i = a+1; i < unocc_orb; i++)
                {
                    if ((a+b+i)%2 == 0)
                    {
                        J = 0;
                        for (int j = 0; j < n_Electrons; j++)
                        {
                            val1 = EqualFunc((a+n_Electrons)%2, (b+n_Electrons)%2) * EqualFunc((i+n_Electrons)%2, j%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) ((b+n_Electrons)/2), (int) ((i+n_Electrons)/2), (int) (j/2));

                            val2 = EqualFunc((a+n_Electrons)%2, j%2) * EqualFunc((i+n_Electrons)%2, (b+n_Electrons)%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) (j/2), (int) ((i+n_Electrons)/2), (int) ((b+n_Electrons)/2));

                            Rear_cube9(a,i)(b,j) = val1-val2;
                            j++;
                            J++;
                        }
                    }

                    else
                    {
                        J = 0;
                        for (int j = 1; j < n_Electrons; j++)
                        {
                            val1 = EqualFunc((a+n_Electrons)%2, (b+n_Electrons)%2) * EqualFunc((i+n_Electrons)%2, j%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) ((b+n_Electrons)/2), (int) ((i+n_Electrons)/2), (int) (j/2));

                            val2 = EqualFunc((a+n_Electrons)%2, j%2) * EqualFunc((i+n_Electrons)%2, (b+n_Electrons)%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) (j/2), (int) ((i+n_Electrons)/2), (int) ((b+n_Electrons)/2));

                            Rear_cube9(a,i)(b,j) = val1-val2;
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
                for (int i = a+1; i < unocc_orb; i++)
                {
                    if ((a+b+i)%2 == 0)
                    {
                        J = 0;
                        for (int j = 0; j < n_Electrons; j++)
                        {
                            val1 = EqualFunc((a+n_Electrons)%2, (b+n_Electrons)%2) * EqualFunc((i+n_Electrons)%2, j%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) ((b+n_Electrons)/2), (int) ((i+n_Electrons)/2), (int) (j/2));

                            val2 = EqualFunc((a+n_Electrons)%2, j%2) * EqualFunc((i+n_Electrons)%2, (b+n_Electrons)%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) (j/2), (int) ((i+n_Electrons)/2), (int) ((b+n_Electrons)/2));

                            Rear_cube9(a,i)(b,j) = val1-val2;
                            j++;
                            J++;
                        }
                    }

                    else
                    {
                        J = 0;
                        for (int j = 1; j < n_Electrons; j++)
                        {
                            val1 = EqualFunc((a+n_Electrons)%2, (b+n_Electrons)%2) * EqualFunc((i+n_Electrons)%2, j%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) ((b+n_Electrons)/2), (int) ((i+n_Electrons)/2), (int) (j/2));

                            val2 = EqualFunc((a+n_Electrons)%2, j%2) * EqualFunc((i+n_Electrons)%2, (b+n_Electrons)%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) (j/2), (int) ((i+n_Electrons)/2), (int) ((b+n_Electrons)/2));

                            Rear_cube9(a,i)(b,j) = val1-val2;
                            j++;
                            J++;
                        }
                    }
                }
                b++;
                B++;
            }
        }

        // Compact form of I(a,b,c,i), dont store any zeroes
        Rear_cube10.set_size(unocc_orb, n_Electrons);
        field_filler_matrix2 = zeros(unocc_orb, n_Electrons);
        for (int i = 0; i < n_Electrons; i++)
        {
            for (int a = 0; a < unocc_orb; a++)
            {
                Rear_cube10(a,i) = field_filler_matrix2;
            }
        }

        for (int a = 0; a < unocc_orb; a++)
        {
            B = 0;
            for (int b = 0; b < unocc_orb; b++)
            {
                for (int i = 0; i < n_Electrons; i++)
                {
                    if ((a+b+i)%2 == 0)
                    {
                        J = 0;
                        for (int j = 0; j < n_Electrons; j++)
                        {
                            val1 = EqualFunc((a+n_Electrons)%2, (b+n_Electrons)%2) * EqualFunc(j%2, i%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) ((b+n_Electrons)/2), (int) (j/2), (int) (i/2));

                            val2 = EqualFunc((a+n_Electrons)%2, i%2) * EqualFunc(j%2, (b+n_Electrons)%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) (i/2), (int) (j/2), (int) ((b+n_Electrons)/2));

                            Rear_cube10(a,i)(b,j) = val1-val2;
                            j++;
                            J++;
                        }
                    }

                    else
                    {
                        J = 0;
                        for (int j = 1; j < n_Electrons; j++)
                        {
                            val1 = EqualFunc((a+n_Electrons)%2, (b+n_Electrons)%2) * EqualFunc(j%2, i%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) ((b+n_Electrons)/2), (int) (j/2), (int) (i/2));

                            val2 = EqualFunc((a+n_Electrons)%2, i%2) * EqualFunc(j%2, (b+n_Electrons)%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) (i/2), (int) (j/2), (int) ((b+n_Electrons)/2));

                            Rear_cube10(a,i)(b,j) = val1-val2;
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
                    if ((a+b+i)%2 == 0)
                    {
                        J = 0;
                        for (int j = 0; j < n_Electrons; j++)
                        {
                            val1 = EqualFunc((a+n_Electrons)%2, (b+n_Electrons)%2) * EqualFunc(j%2, i%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) ((b+n_Electrons)/2), (int) (j/2), (int) (i/2));

                            val2 = EqualFunc((a+n_Electrons)%2, i%2) * EqualFunc(j%2, (b+n_Electrons)%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) (i/2), (int) (j/2), (int) ((b+n_Electrons)/2));

                            Rear_cube10(a,i)(b,j) = val1-val2;
                            j++;
                            J++;
                        }
                    }

                    else
                    {
                        J = 0;
                        for (int j = 1; j < n_Electrons; j++)
                        {
                            val1 = EqualFunc((a+n_Electrons)%2, (b+n_Electrons)%2) * EqualFunc(j%2, i%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) ((b+n_Electrons)/2), (int) (j/2), (int) (i/2));

                            val2 = EqualFunc((a+n_Electrons)%2, i%2) * EqualFunc(j%2, (b+n_Electrons)%2) *
                                    temp4.at((int) ((a+n_Electrons)/2))((int) (i/2), (int) (j/2), (int) ((b+n_Electrons)/2));

                            Rear_cube10(a,i)(b,j) = val1-val2;
                            j++;
                            J++;
                        }
                    }
                }
                b++;
                B++;
            }
        }

        temp4.clear();
}

int coupled_cluster_integrals::Return_Integral_Index(int a, int b, int i, int j)
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

field<mat> coupled_cluster_integrals::Return_Integrals_Part_2()
{
    return temp_cube33;
}

field<mat> coupled_cluster_integrals::Return_V2_integrals()
{
    return temp_cube3;
}

void coupled_cluster_integrals::Delete_Everything()
{
    Rear_cube.set_size(0,0);
    Rear_cube2.set_size(0,0);
    Rear_cube3.set_size(0,0);
    Rear_cube4.set_size(0,0);
    Rear_cube6.set_size(0,0);
    Rear_cube7.set_size(0,0);
    Rear_cube8.set_size(0,0);
    Rear_cube9.set_size(0,0);
    Rear_cube10.set_size(0,0);
}

double coupled_cluster_integrals::EqualFunc(int a, int b)
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

mat coupled_cluster_integrals::Return_FS(vec eigval)
{
    mat fs = zeros(2*Matrix_Size, 2*Matrix_Size);

    for (int i = 0; i < 2*Matrix_Size; i++)
    {
        fs(i,i) = eigval((int) i/2); // Diagonal matrise
    }

    return fs;
}

field<mat> coupled_cluster_integrals::Return_Rearranged_Integrals2()
{
    return Rear_cube;
}

field<mat> coupled_cluster_integrals::Return_Rearranged_Integrals3()
{
    return Rear_cube2;
}

field<mat> coupled_cluster_integrals::Return_Rearranged_Integrals4()
{
    return Rear_cube3;
}

field<mat> coupled_cluster_integrals::Return_Rearranged_Integrals5()
{
    return Rear_cube4;
}

field<mat> coupled_cluster_integrals::Return_Rearranged_Integrlas6()
{
    return Rear_cube6;
}

field<mat> coupled_cluster_integrals::Return_Rearranged_Integrals7()
{
    return Rear_cube7;
}

field<mat> coupled_cluster_integrals::Return_Rearranged_Integrals8()
{
    return Rear_cube8;
}

field<mat> coupled_cluster_integrals::Return_Rearranged_Integrals9()
{
    return Rear_cube9;
}

field<mat> coupled_cluster_integrals::Return_Rearranged_Integrals10()
{
    return Rear_cube10;
}
