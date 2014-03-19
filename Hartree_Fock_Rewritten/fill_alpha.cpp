#include "fill_alpha.h"

Fill_Alpha::Fill_Alpha(int n_N, vec zZz, string b_s, int mat_size)
{
    n_Nuclei = n_N;
    Z = zZz;
    Basis_Set = b_s;
    Matrix_Size = mat_size;
}

mat Fill_Alpha::Fyll_Opp_Alpha()
{
    int max_numb_bas = 5;
    uint i;
    int index_counter = -1;
    n_Orbitals = zeros(n_Nuclei);

    alpha = zeros(Matrix_Size, max_numb_bas);
    c = zeros(Matrix_Size, max_numb_bas);
    Basis_Per_Orbital = zeros(Matrix_Size);
    Potenser = zeros(Matrix_Size, 3);

    if (Basis_Set == "Thjiisen")
    {
        for (i=0; i<n_Nuclei; i++)
        {
            if (Z(i) == 1)
            {
                n_Orbitals(i) = 4;

                index_counter += 1;
                Basis_Per_Orbital(index_counter) = 1;
                alpha(index_counter,0) = 13.00773;
                c(index_counter,0) = 1.0;

                index_counter += 1;
                Basis_Per_Orbital(index_counter) = 1;
                alpha(index_counter,0) = 1.962079;
                c(index_counter,0) = 1.0;

                index_counter += 1;
                Basis_Per_Orbital(index_counter) = 1;
                alpha(index_counter,0) = 0.444529;
                c(index_counter,0) = 1.0;

                index_counter += 1;
                Basis_Per_Orbital(index_counter) = 1;
                alpha(index_counter,0) = 0.1219492;
                c(index_counter,0) = 1.0;
            }
        }
    }


    if (Basis_Set == "3-21G")
    {
        for (i=0; i < n_Nuclei; i++)
        {
            if (Z(i) == 1)
            {
                n_Orbitals(i) = 2;

                index_counter += 1;
                Basis_Per_Orbital(index_counter) = 2;
                alpha(index_counter,1)  = 5.4471780;
                c(index_counter, 1)     = 0.1562850;
                alpha(index_counter, 0) = 0.8245470;
                c(index_counter, 0)     = 0.9046910;

                index_counter += 1;
                Basis_Per_Orbital(index_counter) = 1;
                alpha(index_counter, 0) = 0.1831920;
                c(index_counter, 0)     = 1.0;

            }

            if (Z(i) == 4)
            {
                n_Orbitals(i) = 5;

                index_counter += 1;
                Basis_Per_Orbital(index_counter) = 3;
                alpha(index_counter,0)  = 2.2220500;
                c(index_counter,0)      = 0.6959340;
                alpha(index_counter,1)  = 10.728900;
                c(index_counter,1)      = 0.3660960;
                alpha(index_counter,2)  = 71.887600;
                c(index_counter,2)      = 0.0644263;

                index_counter += 1;
                Basis_Per_Orbital(index_counter) = 2;
                alpha(index_counter,0)  = 0.2688810;
                c(index_counter,1)      = -0.421064;
                alpha(index_counter,1)  = 1.2954800;
                c(index_counter,0)      = 1.2240700;

                index_counter += 1;
                Potenser(index_counter,2) = 1;
                Basis_Per_Orbital(index_counter) = 2;
                alpha(index_counter,0)  = 0.2688810;
                c(index_counter,0)      = 0.8825280;
                alpha(index_counter,1)  = 1.2954800;
                c(index_counter,1)      = 0.2051320;

                index_counter += 1;
                Basis_Per_Orbital(index_counter) = 1;
                alpha(index_counter,0)  = 0.0773500;
                c(index_counter,0)      = 1.0000000;

                index_counter += 1;
                Potenser(index_counter,1) = 1;
                Basis_Per_Orbital(index_counter) = 1;
                alpha(index_counter,0)  = 0.0773500;
                c(index_counter,0)      = 1.0000000;
            }
        }
    }

    if (Basis_Set == "4-31G")
    {
        for (i=0; i<n_Nuclei; i++)
        {
            if (Z(i) == 8)
            {
                n_Orbitals(i) = 5;

                index_counter += 1;
                Basis_Per_Orbital(index_counter) = 4;
                alpha(index_counter,0)  = 7.9786772;
                c(index_counter,0)      = 0.5600108;
                alpha(index_counter,1)  = 29.906408;
                c(index_counter,1)      = 0.4348836;
                alpha(index_counter,2)  = 133.12928;
                c(index_counter,2)      = 0.1228292;
                alpha(index_counter,3)  = 883.27286;
                c(index_counter,4)      = 0.0175506;

                index_counter += 1;
                Basis_Per_Orbital(index_counter) = 3;
                alpha(index_counter,0)  = 1.0709836;
                c(index_counter,0)      = 1.1594979;
                alpha(index_counter,1)  = 3.7800860;
                c(index_counter,1)      = -0.1772865;
                alpha(index_counter,2)  = 16.194447;
                c(index_counter,2)      = -0.1134010;

                index_counter += 1;
                Potenser(index_counter,0) = 1;
                Basis_Per_Orbital(index_counter) = 3;
                alpha(index_counter,0)  = 1.0709836;
                c(index_counter,0)      = 0.7346079;
                alpha(index_counter,1)  = 3.7800860;
                c(index_counter,1)      = 0.3312254;
                alpha(index_counter,2)  = 16.194447;
                c(index_counter,2)      = 0.0685453;

                index_counter += 1;
                Basis_Per_Orbital(index_counter) = 1;
                alpha(index_counter,0) = 0.2838798;
                c(index_counter,0) = 1.0;

                index_counter += 1;
                Potenser(index_counter,1) = 1;
                Basis_Per_Orbital(index_counter) = 1;
                alpha(index_counter,0) = 0.2838798;
                c(index_counter,0) = 1.0;
            }
        }
    }

    // Fra Z = 2 og ut pls oppdater

    if (Basis_Set == "6-311G")
    {
        for (i=0; i< n_Nuclei; i++)
        {

        if (Z(i) == 1)
        {
          // #BASIS SET: (5s) -> [3s]
            n_Orbitals(i) = 3;

          // H, S
            index_counter += 1;
            Basis_Per_Orbital(index_counter) = 3;
            alpha(index_counter, 2) = 33.8650;
            c(index_counter, 2) = 0.0254938;
            alpha(index_counter, 1) = 5.09479;
            c(index_counter, 1) = 0.1903730;
            alpha(index_counter, 0) = 1.15879;
            c(index_counter, 0) =  0.8521610;

          // H, S
            index_counter += 1;
            Basis_Per_Orbital(index_counter) = 1;
            alpha(index_counter,0) = 0.3258400;
            c(index_counter, 0) = 1.0000000;

          // H, S
            index_counter += 1;
            Basis_Per_Orbital(index_counter) = 1;
            alpha(index_counter,0) = 0.1027410;
            c(index_counter, 0) = 1.0000000;
        }

          if (Z(i) == 2)
        {
          // #BASIS SET: (5s) -> [3s]
              n_Orbitals(i) = 3;

          // He, S
              index_counter += 1;
              Basis_Per_Orbital(index_counter) = 3;
              alpha(index_counter,0) = 98.1243; c(index_counter, 0) = 0.0287452;
              alpha(index_counter,1) = 14.7689; c(index_counter, 1) = 0.2080610;
              alpha(index_counter,2) = 3.31883; c(index_counter, 2) = 0.8376350;

          // He, S
              index_counter += 1;
              Basis_Per_Orbital(index_counter) = 1;
              alpha(index_counter,0) = 0.874047; c(index_counter, 0) = 1.000000;

          // He, S
              index_counter += 1;
              Basis_Per_Orbital(index_counter) = 1;
              alpha(index_counter,0) = 0.2445640; c(index_counter, 0) = 1.00000;
        }

          if (Z(i) == 3)
        {
          // #BASIS SET: (11s,5p) -> [4s,3p]
              n_Orbitals(i) = 5;

          // Li, S
              index_counter += 1;
              Basis_Per_Orbital(index_counter) = 6;
              alpha(index_counter,0) = 900.4600000; c(index_counter,0) = 0.00228704;
          alpha(index_counter,1) = 134.4330000; c(index_counter, 1) = 0.0176350;
          alpha(index_counter,2) = 30.4365000; c(index_counter,2) = 0.0873434;
          alpha(index_counter,3) = 8.6263900; c(index_counter,3) =  0.2809770;
          alpha(index_counter,4) = 2.4833200; c(index_counter,4) = 0.6587410;
          alpha(index_counter,5) = 0.3031790; c(index_counter,5) = 0.1187120;

          // Li, SP
          index_counter += 1;
          Basis_Per_Orbital(index_counter) = 3;
          alpha(index_counter,6) = 4.8689000; c(index_counter, 0) = 0.0933293; // 0.0327661
          alpha(index_counter,7) = 0.8569240; c(index_counter, 1) = 0.9430450; // 0.1597920
          alpha(index_counter,8) = 0.2432270; c(index_counter, 2) = -0.00279827; // 0.8856670

          // Li, SP
          index_counter += 1;
          Basis_Per_Orbital(index_counter) = 1;
          alpha(index_counter,9) = 0.0635070; c(index_counter, 0 ) = 1.0;//, 1.0

          // Li, SP
          index_counter += 1;
          Basis_Per_Orbital(index_counter) = 1;
          alpha(index_counter,10) = 0.0243683; c(index_counter, 0 ) = 1.0;//, 1.0
        }

          if (Z(i) == 4)
        {
          // #BASIS SET: (11s,5p) -> [4s,3p]

          // Be, S
              index_counter += 1;
              Basis_Per_Orbital(index_counter) = 6;
          alpha(index_counter,0) = 1682.8000000; c(index_counter, 0) = 0.00228574;
          alpha(index_counter,1) = 251.7150000; c(index_counter, 1) = 0.0175938;
          alpha(index_counter,2) = 57.4116000; c(index_counter, 2) = 0.0863315;
          alpha(index_counter,3) = 16.5171000; c(index_counter, 3) = 0.2818350;
          alpha(index_counter,4) = 4.8536400; c(index_counter, 4) = 0.6405940;
          alpha(index_counter,5) = 0.6268630; c(index_counter, 5) = 0.1444670;

          // Be, SP
          index_counter += 1;
          Basis_Per_Orbital(index_counter) = 3;
          alpha(index_counter,6) = 8.3093800; c(index_counter, 0) = 0.1086210;// 0.0361344
          alpha(index_counter,7) = 1.7407500; c(index_counter, 1) = 0.9273010;// 0.2169580
          alpha(index_counter,8) = 0.4858160; c(index_counter,2) = -0.00297169; // 0.8418390

          // Be, SP
          index_counter += 1;
          Basis_Per_Orbital(index_counter) = 1;
          alpha(index_counter,9) = 0.1636130; c(index_counter,0) = 1.0000000;c(index_counter, 0 ) = 1.0;

          // Be, SP
          index_counter += 1;
          Basis_Per_Orbital(index_counter) = 1;
          alpha(index_counter,10) = 0.0567285; c(index_counter,0) = 1.0000000; c(index_counter, 0 ) = 1.0;
        }

          if (Z(i) == 5)
        {
          // #BASIS SET: (11s,5p) -> [4s,3p]

          // B, S
              index_counter += 1;
              Basis_Per_Orbital(index_counter) = 6;
          alpha(index_counter,0) = 2858.8900000; c(index_counter, 0) = 0.00215375;
          alpha(index_counter,1) =428.1400000; c(index_counter, 1) = 0.0165823;
          alpha(index_counter,2) = 97.5282000; c(index_counter, 2) = 0.0821870;
          alpha(index_counter,3) = 27.9693000; c(index_counter, 3) = 0.2766180;
          alpha(index_counter,4) = 8.2157700; c(index_counter, 4) = 0.6293160;
          alpha(index_counter,5) = 1.1127800; c(index_counter, 5) = 0.1737700;

          // B, SP
          index_counter += 1;
          Basis_Per_Orbital(index_counter) = 3;
          alpha(index_counter,6) = 13.2415000; c(index_counter,0 ) = 0.1174430;// 0.0418100
          alpha(index_counter,7) = 3.0016600; c(index_counter, 1) = 0.9180020;// 0.2365750
          alpha(index_counter,8) = 0.9128560; c(index_counter,2) = -0.00265105;// 0.8162140

          // B, SP
          index_counter += 1;
          Basis_Per_Orbital(index_counter) = 1;
          alpha(index_counter,9) =  0.3154540; c(index_counter, 0 ) = 1.0;//, 1.0000000

          // B, SP
          index_counter += 1;
          Basis_Per_Orbital(index_counter) = 1;
          alpha(index_counter,10) = 0.0988563; c(index_counter, 0 ) = 1.0;//, 1.0000000
        }

          if (Z(i) == 6)
        {
          // #BASIS SET: (11s,5p) -> [4s,3p]

          // C, S
              index_counter += 1;
              Basis_Per_Orbital(index_counter) = 6;
          alpha(index_counter,0) = 4563.2400000; c(index_counter, 0) = 0.00196665;
          alpha(index_counter,1) = 682.0240000; c(index_counter, 1) = 0.0152306;
          alpha(index_counter,2) = 154.9730000; c(index_counter, 2) = 0.0761269;
          alpha(index_counter,3) = 44.4553000; c(index_counter, 3) = 0.2608010;
          alpha(index_counter,4) = 13.0290000; c(index_counter, 4) = 0.6164620;
          alpha(index_counter,5) = 1.8277300; c(index_counter, 5) = 0.2210060;

          // C, SP
          index_counter += 1;
          Basis_Per_Orbital(index_counter) = 3;
          alpha(index_counter,6) = 20.9642000; c(index_counter, 0) = 0.1146600, 0.0402487;
          alpha(index_counter,7) = 4.8033100; c(index_counter, 1) = 0.9199990, 0.2375940;
          alpha(index_counter,8) = 1.4593300; c(index_counter,2) = -0.00303068;//0.8158540

          // C, SP
          index_counter += 1;
          Basis_Per_Orbital(index_counter) = 1;
          alpha(index_counter,9) = 0.4834560; c(index_counter, 0 ) = 1.0;//, 1.0000000

          // C, SP
          index_counter += 1;
          Basis_Per_Orbital(index_counter) = 1;
          alpha(index_counter,10) = 0.1455850; c(index_counter, 0 ) = 1.0;//, 1.0000000
        }

          if (Z(i) == 7)

        {
          // #BASIS SET: (11s,5p) -> [4s,3p]

          // N, S
              index_counter += 1;
              Basis_Per_Orbital(index_counter) = 6;
          alpha(index_counter,0) = 6293.4800000; c(index_counter, 0) = 0.00196979;
          alpha(index_counter,1) = 949.0440000; c(index_counter, 1) = 0.0149613;
          alpha(index_counter,2) = 218.7760000; c(index_counter, 2) = 0.0735006;
          alpha(index_counter,3) = 63.6916000; c(index_counter, 3) = 0.2489370;
          alpha(index_counter,4) = 18.8282000; c(index_counter, 4) = 0.6024600;
          alpha(index_counter,5) = 2.7202300; c(index_counter, 5) = 0.2562020;

          // N, SP
          index_counter += 1;
          Basis_Per_Orbital(index_counter) = 3;
          alpha(index_counter,0) = 30.6331000; c(index_counter, 0) = 0.1119060;// 0.0383119
          alpha(index_counter,1) = 7.0261400; c(index_counter, 1) = 0.9216660;// 0.2374030
          alpha(index_counter,2) = 2.1120500; c(index_counter,2) = -0.00256919; // 0.8175920

          // N, SP
          index_counter += 1;
          Basis_Per_Orbital(index_counter) = 1;
          alpha(index_counter,0) = 0.6840090; c(index_counter, 0 ) = 1.0;//, 1.0000000

          // N, SP
          index_counter += 1;
          Basis_Per_Orbital(index_counter) = 1;
          alpha(index_counter,0) = 0.2008780; c(index_counter, 0 ) = 1.0;//, 1.0000000
        }

        if (Z(i) == 8)
          {
            // #BASIS SET: (11s,5p) -> [4s,3p]

            // O, S
            index_counter += 1;
            Basis_Per_Orbital(index_counter) = 6;
            alpha(index_counter,0) = 8588,50; c(index_counter,0) = 0.00189515;
            alpha(index_counter,1) = 1297.2300000; c(index_counter, 1) = 0.0143859;
            alpha(index_counter,2) = 299.2960000; c(index_counter, 2) = 0.0707320;
            alpha(index_counter,3) = 87.3771000; c(index_counter, 3) = 0.240001;
            alpha(index_counter,4) = 25.6789000; c(index_counter, 4) = 0.594797;
            alpha(index_counter,5) = 3.7400400; c(index_counter, 5) = 0.2808020;

            // O, SP
            index_counter += 1;
            Basis_Per_Orbital(index_counter) = 3;
            alpha(index_counter,0) = 42.1175000; c(index_counter, 0) = 0.1138890;// 0.0365114
            alpha(index_counter,1) = 9.6283700; c(index_counter, 1) = 0.9208110; // 0.2371530;
            alpha(index_counter,2) = 2.8533200; c(index_counter,2) = -0.00327447;// 0.8197020

            // O, SP
            alpha(index_counter,0) = 0.9056610; c(index_counter, 0 ) = 1.0;//, 1.0000000

            // O, SP
            alpha(index_counter,10) = 0.2556110; c(index_counter, 0 ) = 1.0;//, 1.0
          }

        if (Z(i) == 9)
          {
            // #BASIS SET: (11s,5p) -> [4s,3p]

            // F, S
            alpha(index_counter,0) = 11427.1; c(index_counter, 0) = 0.00180093;
            alpha(index_counter,1) = 1722.35; c(index_counter, 1) = 0.01347419;
            alpha(index_counter,2) = 395.746; c(index_counter, 2) = 0.0681334;
            alpha(index_counter,3) = 115.139; c(index_counter, 3) = 0.2333250;
            alpha(index_counter,4) = 33.6026; c(index_counter, 4) = 0.5890860;
            alpha(index_counter,5) = 4.91901; c(index_counter, 5) = 0.2995050;

            // F, SP
            alpha(index_counter,6) = 55.4441; c(index_counter, 0) = 0.1145360; //, 0.0354609
            alpha(index_counter,7) = 12.6323; c(index_counter, 1) = 0.9205120; // 0.2374510
            alpha(index_counter,8) = 3.71756; c(index_counter, 2) = -0.00337804; // 0.8204580

            // F, SP
            alpha(index_counter,9) = 1.16545; c(index_counter, 0 ) = 1.0;//, 1.0

            // F, SP
            alpha(index_counter,10) = 0.3218920; c(index_counter, 0 ) = 1.0;//, 1.0
          }

        if (Z(i) == 10)
          {
            // #BASIS SET: (11s,5p) -> [4s,3p]

            // Ne, Sc(index_counter, 0 ) = 1.0;
            alpha(index_counter,0) = 13995.7; c(index_counter, 0) = 0.00183276;
            alpha(index_counter,1) = 2117.1; c(index_counter, 1) = 0.0138827;
            alpha(index_counter,2) = 490.425; c(index_counter, 2) = 0.0680687;
            alpha(index_counter,3) = 143,833; c(index_counter, 3) = 0.2313280;
            alpha(index_counter,4) = 41.9265; c(index_counter, 4) = 0.58589;
            alpha(index_counter,5) = 6.15684; c(index_counter, 5) = 0.305883;

            // Ne, SP
            alpha(index_counter,6) = 69.1211; c(index_counter, 0) = 0.119149; //, 0.0356574
            alpha(index_counter,7) = 15.835; c(index_counter, 1) = 0.9173750; //, 0.2394770
            alpha(index_counter,8) = 4.67326; c(index_counter, 2) = -0.00405839; //, 0.8184610

            // Ne, SP
            alpha(index_counter,9) = 1.45756; c(index_counter, 0 ) = 1.0;//, 1.0

            // Ne, SP
            alpha(index_counter,10) = 0.397057; c(index_counter, 0 ) = 1.0;//, 1.0
          }

        if (Z(i) == 11)
          {
            // #BASIS SET: (13s,9p) -> [6s,5p]

            // Na, S
            alpha(index_counter,0) = 36166.4; c(index_counter, 0) = 0.001032;
            alpha(index_counter,1) = 5372.58; c(index_counter, 1) = 0.0080710;
            alpha(index_counter,2) = 1213.21; c(index_counter, 2) = 0.0421290;
            alpha(index_counter,3) = 339.623; c(index_counter, 3) = 0.1697890;
            alpha(index_counter,4) = 109.553; c(index_counter, 4) = 0.5146210;
            alpha(index_counter,5) = 38.7773; c(index_counter, 5) = 0.3798170;

            // Na, S
            alpha(index_counter,6) = 38.7773; c(index_counter, 0) = 0.3747620;
            alpha(index_counter,7) = 14.5759; c(index_counter, 1) = 0.5757690;
            alpha(index_counter,8) = 5.27993; c(index_counter, 2) = 0.1129330;

            // Na, S
            alpha(index_counter,9) = 1.82777; c(index_counter, 0 ) = 1.0;

            // Na, S
            alpha(index_counter,10) = 0.619948; c(index_counter, 0 ) = 1.0;

            // Na, S
            alpha(index_counter,11) = 0.05724; c(index_counter, 0 ) = 1.0;

            // Na, S
            alpha(index_counter,12) = 0.0240480; c(index_counter, 0 ) = 1.0;

            // Na, P
            alpha(index_counter,13) = 144.645; c(index_counter, 0) = 0.0114850;
            alpha(index_counter,14) = 33.9074; c(index_counter, 1) = 0.0823830;
            alpha(index_counter,15) = 10.6285; c(index_counter, 2) = 0.3196580;
            alpha(index_counter,16) = 3.82389; c(index_counter, 3) = 0.7012950;

            // Na, P
            alpha(index_counter,17) = 1.44429; c(index_counter, 0) = 0.6385060;
            alpha(index_counter,18) = 0.55261; c(index_counter, 1) = 0.4253650;

            // Na, P
            alpha(index_counter,19) = 0.18872; c(index_counter, 0 ) = 1.0;

            // Na, P
            alpha(index_counter,20) = 0.046501; c(index_counter, 0 ) = 1.0;

            // Na, P
            alpha(index_counter,21) = 0.016285; c(index_counter, 0 ) = 1.0;
          }

        if (Z(i) == 12)
          {
            // #BASIS SET: (13s,9p) -> [6s,5p]

            // Mg, S
            alpha(index_counter,0) = 43866.5; c(index_counter, 0) = 0.0009180;
            alpha(index_counter,1) = 6605.37; c(index_counter, 1) = 0.0070470;
            alpha(index_counter,2) = 1513.26; c(index_counter, 2) = 0.0359410;
            alpha(index_counter,3) = 432.317; c(index_counter, 3) = 0.1414610;
            alpha(index_counter,4) = 142.149; c(index_counter, 4) = 0.4268640;
            alpha(index_counter,5) = 51.3983; c(index_counter, 5) = 0.4979750;

            // Mg, S
            alpha(index_counter,6) = 51.3983; c(index_counter, 0) = 0.2513550;
            alpha(index_counter,7) = 19.9196; c(index_counter, 1) = 0.6186710;
            alpha(index_counter,8) = 8.02474; c(index_counter, 2) = 0.1884170;

            // Mg, S
            alpha(index_counter,9) = 2.50817; c(index_counter, 0 ) = 1.0;

            // Mg, S
            alpha(index_counter,10) = 0.871531; c(index_counter, 0 ) = 1.0;

            // Mg, S
            alpha(index_counter,11) = 0.108188; c(index_counter, 0 ) = 1.0;

            // Mg, S
            alpha(index_counter,12) = 0.040130; c(index_counter, 0 ) = 1.0;

            // Mg, P
            alpha(index_counter,13) = 193.854; c(index_counter, 0) = 0.0101880;
            alpha(index_counter,14) = 45.4420; c(index_counter, 1) = 0.0753600;
            alpha(index_counter,15) = 14.1864; c(index_counter, 2) = 0.3074190;
            alpha(index_counter,16) = 5.05751; c(index_counter, 3) = 0.7175750;

            // Mg, P
            alpha(index_counter,17) = 1.88861; c(index_counter, 0) = 0.6673390;
            alpha(index_counter,18) = 0.722652; c(index_counter, 1) = 0.3946490;

            // Mg, P
            alpha(index_counter,19) = 0.236417; c(index_counter, 0 ) = 1.0;

            // Mg, P
            alpha(index_counter,20) = 0.093358; c(index_counter, 0 ) = 1.0;

            // Mg, P
            alpha(index_counter,21) = 0.034809; c(index_counter, 0 ) = 1.0;
          }

        if (Z(i) == 13)
          {
            // #BASIS SET: (13s,9p) -> [6s,5p]

            // Al, S
            alpha(index_counter,0) = 54866.489; c(index_counter, 0) = 0.0008390;
            alpha(index_counter,1) = 8211.7665; c(index_counter, 1) = 0.0065270;
            alpha(index_counter,2) = 1866.1761; c(index_counter, 2) = 0.0336660;
            alpha(index_counter,3) = 531.12934; c(index_counter, 3) = 0.1329020;
            alpha(index_counter,4) = 175.11797; c(index_counter, 4) = 0.4012660;
            alpha(index_counter,5) = 64.005500; c(index_counter, 5) = 0.5313380;

            // Al, S
            alpha(index_counter,6) = 64.005500; c(index_counter, 0) = 0.2023050;
            alpha(index_counter,7) = 25.292507; c(index_counter, 1) = 0.6247900;
            alpha(index_counter,8) = 10.534910; c(index_counter, 2) = 0.2274390;

            // Al, S
            alpha(index_counter,9) = 3.2067110; c(index_counter, 0 ) = 1.0;

            // Al, S
            alpha(index_counter,10) = 1.152555; c(index_counter, 0 ) = 1.0;

            // Al, S
            alpha(index_counter,11) = 0.176678; c(index_counter, 0 ) = 1.0;

            // Al, S
            alpha(index_counter,12) = 0.065237; c(index_counter, 0 ) = 1.0;

            // Al, P
            alpha(index_counter,13) = 259.28362; c(index_counter, 0) = 0.0094480;
            alpha(index_counter,14) = 61.076870; c(index_counter, 1) = 0.0709740;
            alpha(index_counter,15) = 19.303237; c(index_counter, 2) = 0.2956360;
            alpha(index_counter,16) = 7.0108820; c(index_counter, 3) = 0.7282190;

            // Al, P
            alpha(index_counter,17) = 2.6738650; c(index_counter, 0) = 0.6444670;
            alpha(index_counter,18) = 1.0365960; c(index_counter, 1) = 0.4174130;

            // Al, P
            alpha(index_counter,19) = 0.3168190; c(index_counter, 0 ) = 1.0;

            // Al, P
            alpha(index_counter,20) = 0.1142570; c(index_counter, 0 ) = 1.0;

            // Al, P
            alpha(index_counter,21) = 0.0413970; c(index_counter, 0 ) = 1.0;
          }

          if (Z(i) == 14)
        {
          // #BASIS SET: (13s,9p) -> [6s,5p]

          // Si, S
          alpha(index_counter,0) = 69379.2300000; c(index_counter, 0) = 0.0007570;
          alpha(index_counter,1) = 10354.9400000; c(index_counter, 1) = 0.0059320;
          alpha(index_counter,2) = 2333.8796000; c(index_counter, 2) = 0.0310880;
          alpha(index_counter,3) = 657.1429500; c(index_counter, 3) = 0.1249670;
          alpha(index_counter,4) = 214.3011300; c(index_counter, 4) = 0.3868970;
          alpha(index_counter,5) = 77.6291680; c(index_counter, 5) = 0.5548880;

          // Si, S
          alpha(index_counter,6) = 77.6291680; c(index_counter, 0) = 0.1778810;
          alpha(index_counter,7) = 30.6308070; c(index_counter, 1) = 0.6277650;
          alpha(index_counter,8) = 12.8012950; c(index_counter, 2) = 0.2476230;

          // Si, S
          alpha(index_counter,9) = 3.9268660; c(index_counter, 0 ) = 1.0;

          // Si, S
          alpha(index_counter,10) = 1.4523430; c(index_counter, 0 ) = 1.0;

          // Si, S
          alpha(index_counter,11) = 0.2562340; c(index_counter, 0 ) = 1.0;

          // Si, S
          alpha(index_counter,12) = 0.0942790; c(index_counter, 0 ) = 1.0;

          // Si, P
          alpha(index_counter,13) = 335.4831900; c(index_counter, 0) = 0.0088660;
          alpha(index_counter,14) = 78.9003660; c(index_counter, 1) = 0.0682990;
          alpha(index_counter,15) = 24.9881500; c(index_counter, 2) = 0.2909580;
          alpha(index_counter,16) = 9.2197110; c(index_counter, 3) = 0.7321170;

          // Si, P
          alpha(index_counter,17) = 3.6211400; c(index_counter, 0) = 0.6198790;
          alpha(index_counter,18) = 1.4513100; c(index_counter, 1) = 0.4391480;

          // Si, P
          alpha(index_counter,19) = 0.5049770; c(index_counter, 0 ) = 1.0;

          // Si, P
          alpha(index_counter,20) = 0.1863170; c(index_counter, 0 ) = 1.0;

          // Si, P
          alpha(index_counter,21) = 0.0654320; c(index_counter, 0 ) = 1.0;
        }
        }
        }
    return alpha;
}

mat Fill_Alpha::Fyll_Opp_c()
{
    return c;
}

vec Fill_Alpha::Fyll_Opp_Nr_Basis_Functions()
{
    return Basis_Per_Orbital;
}

vec Fill_Alpha::Fyll_Opp_Antall_Orbitaler()
{
    return n_Orbitals;
}

mat Fill_Alpha::Fyll_Opp_Potenser()
{
    return Potenser;
}
