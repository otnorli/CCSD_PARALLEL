#include "matrix_size_setter.h"

matrix_size_setter::matrix_size_setter(vec zZz, string sys_string, int n_N)
{
    System_String = sys_string;
    Z = zZz;
    n_Nuclei = n_N;
}

int matrix_size_setter::Set_Matrix_Size()
{
    int i;
    value = 0;
    max_bas_func = 0;
    int security = 5;

    for (i=0; i<n_Nuclei; i++)
    {
        if (System_String == "aug-cc-pVDZ")
        {
            if (Z(i) == 1)
            {
                value += 9;
                if (max_bas_func < 3)
                {
                    max_bas_func = 3;
                }
            }

            if (Z(i) == 17)
            {
                value += 29;
                if (max_bas_func < 11)
                {
                    max_bas_func= 11;
                }
            }

            if (Z(i) == 8)
            {
                value += 25;

                if (max_bas_func < 8)
                {
                    max_bas_func = 8;
                }
            }

            if (Z(i) == 16)
            {
                value += 27;
                if (max_bas_func < 11)
                {
                    max_bas_func = 11;
                }
            }

            if ((Z(i) == 21) || (Z(i) == 22) || (Z(i) == 26) || (Z(i) == 24))
            {
                value += 69;
                if (max_bas_func < 19)
                {
                    max_bas_func = 19;
                }
            }
        }

        if (System_String == "aug-cc-pVTZ")
        {
            if (Z(i) == 4)
            {
                value += 55;
                if (max_bas_func < 9)
                {
                    max_bas_func = 9;
                }
            }
        }

        if (System_String == "aug-cc-pVQZ")
        {
            if (Z(i) == 1)
            {
                value += 55;
                if (max_bas_func < 3)
                {
                    max_bas_func = 3;
                }
            }
        }

        if (System_String == "dz")
        {
            if (Z(i) == 1)
            {
                value += 2;
                if (max_bas_func < 3)
                {
                    max_bas_func = 3;
                }
            }

            if (Z(i) == 8)
            {
                value += 10;
                if (max_bas_func < 6)
                {
                    max_bas_func = 6;
                }
            }
        }

        if (System_String == "dzp")
        {
            if (Z(i) == 1)
            {
                value += 5;
                if (max_bas_func < 3)
                {
                    max_bas_func = 3;
                }
            }

            if (Z(i) == 5)
            {
                value += 16;
                if (max_bas_func < 5)
                {
                    max_bas_func = 5;
                }
            }

            if (Z(i) == 9)
            {
                value += 16;
                if (max_bas_func < 6)
                {
                    max_bas_func = 6;
                }
            }
        }

        if (System_String == "dzvp")
        {
            if (Z(i) == 1 || Z(i) == 2)
            {
                value += 2;
                if (max_bas_func < 4)
                {
                    max_bas_func = 4;
                }
            }

            if (Z(i) > 2 && Z(i) < 11)
            {
                value += 19;
                if (max_bas_func < 6)
                {
                    max_bas_func = 6;
                }
            }
        }

        if (System_String == "SARC-DKH")
        {
            if (Z(i) == 92)
            {
                value += 190;

                if (max_bas_func < 9)
                {
                    max_bas_func = 9;
                }
            }
        }

        if (System_String == "6-311-2d2p")
        {
            if (Z(i) == 1 || Z(i) == 2)
            {
                value += 10;
                if (max_bas_func < 3)
                {
                    max_bas_func = 3;
                }
            }

            else if (Z(i) > 2 && Z(i) < 11)
            {
                value += 29;
                if (max_bas_func < 6)
                {
                    max_bas_func = 6;
                }
            }

            else if (Z(i) > 10 && Z(i) < 19)
            {
                value += 37;
                if (max_bas_func < 6)
                {
                    max_bas_func = 6;
                }
            }

            else if (Z(i) > 18 && Z(i) < 21)
            {
                value += 51;
                if (max_bas_func < 6)
                {
                    max_bas_func = 6;
                }
            }
        }

        if (System_String == "6-31G")
        {
            if (Z(i) == 1 || Z(i) == 2)
            {
                value += 2;
                if (max_bas_func < 3)
                {
                    max_bas_func = 3;
                }
            }

            if (Z(i) > 2 && Z(i) < 11)
            {
                value += 9;
                if (max_bas_func < 6)
                {
                    max_bas_func = 6;
                }
            }

            if (Z(i) > 10 && Z(i) < 19)
            {
                value += 13;
                if (max_bas_func < 6)
                {
                    max_bas_func = 6;
                }
            }

            if (Z(i) > 18 && Z(i) < 21)
            {
                value += 17;
                if (max_bas_func < 6)
                {
                    max_bas_func = 6;
                }
            }

            if (Z(i) > 20 && Z(i) < 31)
            {
                value +=29;
                if (max_bas_func < 6)
                {
                    max_bas_func = 6;
                }
            }
        }

        if (System_String == "6-311ss")
        {
            if (Z(i) == 1 || Z(i) == 2)
            {
                value += 6;
                if (max_bas_func < 3)
                {
                    max_bas_func = 3;
                }
            }

            if (Z(i) > 2 && Z(i) < 11)
            {
                value += 19;
                if (max_bas_func < 6)
                {
                    max_bas_func = 6;
                }
            }

            if (Z(i) > 10 && Z(i) < 19)
            {
                value += 27;
                if (max_bas_func < 6)
                {
                    max_bas_func = 6;
                }
            }

            if (Z(i) > 18 && Z(i) < 21)
            {
                value += 41;
                if (max_bas_func < 6)
                {
                    max_bas_func = 6;
                }
            }

            if (Z(i) > 30 && Z(i) < 37)
            {
                value += 47;
                if (max_bas_func < 6)
                {
                    max_bas_func = 6;
                }
            }

            if (Z(i) == 53)
            {
                value += 67;
                if (max_bas_func < 6)
                {
                    max_bas_func = 6;
                }
            }
        }

        if (System_String == "cc-pVDZ")
        {
            if (Z(i) == 1 || Z(i) == 2)
            {
                value += 5;
                if (max_bas_func < 4)
                {
                    max_bas_func = 4;
                }
            }

            if (Z(i) > 2 && Z(i) < 11)
            {
                value += 15;
                if (max_bas_func < 9)
                {
                    max_bas_func = 9;
                }
            }
        }

        if (System_String == "STO-3G")
        {
            if (Z(i) == 1)
            {
                value += 1;
                if (max_bas_func < 3)
                {
                    max_bas_func = 3;
                }

            }

            if (Z(i) == 3)
            {
                value += 5;
                if (max_bas_func < 3)
                {
                    max_bas_func = 3;
                }
            }

            if (Z(i) == 4 || Z(i) == 5)
            {
                value += 5;
                if (max_bas_func < 3)
                {
                    max_bas_func = 3;
                }
            }

            if (Z(i) == 6)
            {
                value += 5;
                if (max_bas_func < 3)
                {
                    max_bas_func = 3;
                }
            }


            if (Z(i) == 7)
            {
                value += 5;
                if (max_bas_func < 3)
                {
                    max_bas_func = 3;
                }
            }

            if (Z(i) == 8 || Z(i) == 9)
            {
                value += 5;
                if (max_bas_func < 3)
                {
                    max_bas_func = 3;
                }
            }

            if ((Z(i) >= 16) && (Z(i) <= 35))
            {
                value += 19;
                if (max_bas_func < 3)
                {
                    max_bas_func = 3;
                }
            }
/*
            if (Z(i) == 27)
            {
                value += 23;
                if (max_bas_func < 3)
                {
                    max_bas_func = 3;
                }
            }

            if ((Z(i) >= 28) || (Z(i) <= 34))
            {
                value += 29;
                if (max_bas_func < 3)
                {
                    max_bas_func = 3;
                }
            }
            */

           // if (Z(i) >= 37 || Z(i) <= )


        }

        if (System_String == "cc-pVTZ")
        {
            if (Z(i) == 1)
            {
                value += 15;
                if (max_bas_func < 5)
                {
                    max_bas_func = 5;
                }
            }

            if (Z(i) == 4)
            {
                value += 35;
                if (max_bas_func < 11)
                {
                    max_bas_func = 11;
                }
            }

            if (Z(i) == 8)
            {
                value += 35;
                if (max_bas_func < 10)
                {
                    max_bas_func = 10;
                }
            }

            if (Z(i) == 10)
            {
                value += 35;
                if (max_bas_func < 10)
                {
                    max_bas_func = 10;
                }
            }
        }

        if (System_String == "3-21G")
        {
            if (Z(i) > 0 && Z(i) < 3)
            {
                value += 2;
                if (max_bas_func < 2)
                {
                    max_bas_func = 2;
                }
            }

            if (Z(i) > 2 && Z(i) < 11)
            {
                value += 9;
                if (max_bas_func < 3)
                {
                    max_bas_func = 3;
                }

            }

            if (Z(i) > 10 && Z(i) < 19)
            {
                value += 13;
            }

            if (Z(i) > 18 && Z(i) < 21)
            {
                value += 17;
            }

            if (Z(i) > 20 && Z(i) < 31)
            {
                value += 29;
            }

            if (Z(i) > 30 && Z(i) < 36)
            {
                value += 23;
            }

            // Flere atomer nedover her
        }

        if (System_String == "4-31G")
        {
            if (Z(i) == 1)
            {
                value += 2;
                if (max_bas_func < 3)
                {
                    max_bas_func = 3;
                }
            }


            if (Z(i) == 8)
            {
                value += 9;
                if (max_bas_func < 4)
                {
                    max_bas_func = 4;
                }
            }
        }

        if (System_String == "6-311G")
        {
            if (Z(i) == 1)
            {
                value += 3;
            }
        }

        if (System_String == "Thjiisen")
        {
            if (Z(i) == 1)
            {
                value += 4;
            }
        }

        if (System_String == "6-311-3d3p")
        {
            if (Z(i) == 1)
            {
                value += 19;

                if (max_bas_func < 3)
                {
                    max_bas_func = 3;
                }
            }

            if (Z(i) == 3)
            {
                value += 45;
                if (max_bas_func < 6)
                {
                    max_bas_func = 6;
                }
            }

            if (Z(i) == 8)
            {
                value += 45;
                if (max_bas_func < 6)
                {
                    max_bas_func = 6;
                }
            }
        }

        if (max_bas_func == 0)
        {
            max_bas_func = security;
        }
    }

    return value;
}

int matrix_size_setter::Return_Max_Bas_Func()
{
    return max_bas_func;
}
