module Atomdata
using ..AtomicMod:makeatom
#using ..Atomic:AtomicMod


atoms = Dict()

#      name                     name  Z  row col mass nval nsemi  orbitals total_energy    e_orb1, e_orb2
#                                                (amu)                      (Ryd)             (eV)
#
atoms["T"] = makeatom("T",    0, 0, 0, 0.0,  1.0,   0,   [:s], 0.0, [0.0], 0.0);
atoms["X"] = makeatom("X",    0, 0, 0, 0.0,  1.0,   0,   [:s], 0.0, [0.0], 0.0, 0.5);
atoms["Xa"] = makeatom("Xa",    0, 0, 0, 0.0,  1.0,   0,   [:s], 3.0, [3.0* 13.605693122], 0.0, 0.5);

atoms["Is"] = makeatom("Si", 14, 3, 4, 28.09,  4.0,   0,   [:s, :p], -9.15243175, [-10.6594, -3.9477],0.007075865418, 0.5300262399999837);

atoms["H" ] = makeatom("H" , 1, 1, 1, 1.0079,   1.0,   0,   [:s ],     -0.90531367, [-6.3233], 0.0, 0.9194820210526328);
#atoms["H" ] = makeatom("H" , 1, 1, 1, 1.0079,   1.0,   0,   [:s, :p ],     -0.90531367, [-6.3233, 0.3], 0.0, 0.9194820210526328);

atoms["Li"] = makeatom("Li", 3, 2, 1, 6.941,    1.0,   2,   [:s, :p], -14.23107803, [-2.7786, -1.0063], 0.0, 0.3666106694736778);
atoms["Be"] = makeatom("Be", 4, 2, 2, 9.012,    2.0,   2,   [:s, :p], -27.80694399, [-5.5055, -1.9313], 0.0, 0.5967315368420951);
atoms["B" ] = makeatom("B" , 5, 2, 3, 10.81,    3.0,   0,   [:s, :p], -5.84632720, [-9.3213, -3.5217], 0.0, 0.6907072084210543);
atoms["C" ] = makeatom("C" , 6, 2, 4, 12.011,   4.0,   0,   [:s, :p], -10.77360203, [-13.6099, -5.2024],0.003505246254, 0.8299248757894666 );
atoms["N" ] = makeatom("N" , 7, 2, 5, 14.007,   5.0,   0,   [:s, :p], -19.47575537, [-18.4071, -7.0112],0.0, 0.9670884968421092 );
atoms["O" ] = makeatom("O" , 8, 2, 6, 15.999,   6.0,   0,   [:s, :p], -31.71327683, [-23.7413, -8.9545],0.0, 1.100101094736835 );
atoms["F" ] = makeatom("F" , 9, 2, 7, 18.998,   7.0,   0,   [:s, :p], -48.22196210, [-29.6972, -11.0286],0.0, 1.2312884800000015 );
#atoms["F" ] = makeatom("F" , 9, 2, 7, 18.998,   5.0,   0,   [:p], -48.22196210, [-11.0286],0.0, 1.2312884800000015 );

atoms["Na"] = makeatom("Na", 11, 3, 1, 22.98,  1.0,   8,   [:s, :p], -95.06749946, [-2.6455, -0.6556], 0.0, 0.35322667789488926 );
atoms["Mg"] = makeatom("Mg", 12, 3, 2, 24.31,  2.0,   8,   [:s, :p], -124.99127898, [-4.6187, -1.2416], 0.0, 0.49684460631582567);
atoms["Al"] = makeatom("Al", 13, 3, 3, 26.98,  3.0,   0,   [:s, :p],  -6.37406439, [-7.6515, -2.6256], 0.0, 0.45958674526316);

#atoms["Si"] = makeatom("Si", 14, 3, 4, 28.09,  4.0,   0,   [:s, :p, :d], -9.15243181, [-10.7049, -3.9932, -1.0 ],0.0, 0.5459773389473677);
#atoms["Si"] = makeatom("Si", 14, 3, 4, 28.09,  4.0,   0,   [:s, :d, :p], -9.15243181, [-10.7049, -1.0, -3.9932],0.0, 0.5459773389473677);

atoms["Si"] = makeatom("Si", 14, 3, 4, 28.09,  4.0,   0,   [:s, :p], -9.15243181, [-10.7049, -3.9932],0.0, 0.5459773389473677);

#atoms["Is"] = makeatom("Si", 14, 3, 4, 28.09,  4.0,   0,   [:s, :p], -9.15243175, [-10.6594, -3.9477],0.007075865418, 0.5300262399999837);
atoms["P" ] = makeatom("P" , 15, 3, 5, 30.97,  5.0,   0,   [:s, :p], -15.07327848, [-13.8557, -5.4181], 0.0, 0.6280724631578809);
atoms["S" ] = makeatom("S" , 16, 3, 6, 32.06,  6.0,   0,   [:s, :p], -23.70562681, [-17.1360, -6.9201], 0.0, 0.7069649852631505);
atoms["Cl"] = makeatom("Cl", 17, 3, 7, 35.45,  7.0,   0,   [:s, :p], -33.06905525, [-20.5627, -8.5062], 0.0, 0.7840431242104959);

atoms["K" ] = makeatom("K" , 19, 4, 1, 39.098,  1.0,   8,   [:s, :d, :p], -56.84482202 , [-2.1989, -0.1312, -0.6747], 0.0, 0.29296732631577904);
#atoms["K" ] = makeatom("K" , 19, 4, 1, 39.098,  1.0,   8,   [:s, :d], -56.84482202 , [-2.1989, -0.1312], 0.0, 0.29296732631577904);
#atoms["K" ] = makeatom("K" , 19, 4, 1, 39.098,  1.0,   8,   [:s, :d, :p], -56.84482202 , [-2.1989, -0.1312, -0.6747], 0.0, 0.29296732631577904);
#atoms["K" ] = makeatom("K" , 19, 4, 1, 39.098,  1.0,   8,   [:s, :p], -56.84482202 , [-2.1989, -0.6747], 0.0, 0.29296732631577904);
atoms["Ca"] = makeatom("Ca", 20, 4, 2, 40.078,  2.0,   8,   [:s, :d, :p], -74.68919704, [-3.6547, -1.8454, -1.2718], 0.0, 0.39270917052637117);

atoms["Sc"] = makeatom("Sc", 21, 4, 2.05, 44.956,  3.0,   8,   [:s, :d, :p], -93.89389249,  [-4.0850, -3.1428, -1.3456], 0.0, 0.4517225347368296);
atoms["Ti"] = makeatom("Ti", 22, 4, 2.15, 47.867,  4.0,   8,   [:s, :d, :p], -118.63461912, [-4.3470, -4.0656, -1.3310], 0.0, 0.40357462736845867);
atoms["V" ] = makeatom("V" , 23, 4, 2.25, 50.941,  5.0,    8,   [:s, :d, :p], -144.34826807, [-4.3727, -4.2223, -1.2153], 0.0, 0.4090880505264785);
atoms["Cr"] = makeatom("Cr", 24, 4, 2.35, 51.996,  6.0,    8,   [:s, :d, :p], -174.69181247, [-4.3746, -4.2928, -1.0764], 0.0, 0.4194641515788422);
atoms["Mn"] = makeatom("Mn", 25, 4, 2.45, 54.938,  7.0,    8,   [:s, :d, :p], -210.40629067, [-4.3786, -4.3505, -0.9746], 0.0, 0.4316654315786859);
atoms["Fe"] = makeatom("Fe", 26, 4, 2.55, 55.845,  8.0,    8,   [:s, :d, :p], -249.56808439, [-4.3911, -4.4108, -0.8715], 0.0, 0.4458202021055272);
atoms["Co"] = makeatom("Co", 27, 4, 2.65, 58.993,  9.0,    8,   [:s, :d, :p], -297.72249071, [-4.2335, -4.2995, -0.5167], 0.0, 0.4611156042104707);
atoms["Ni"] = makeatom("Ni", 28, 4, 2.75, 58.963,  10.0,   8,   [:s, :d, :p], -342.65159683, [-4.2480, -4.3731, -0.3946], 0.0, 0.4802667873683323);
atoms["Cu"] = makeatom("Cu", 29, 4, 2.85, 63.546,  11.0,   8,   [:s, :d, :p], -403.19763559, [-4.5506, -4.9911, -0.5855], 0.0, 0.5280792084213793);
atoms["Zn"] = makeatom("Zn", 30, 4, 2.95, 65.38,   12.0,   8,   [:s, :d, :p], -461.12493484, [-5.9305, -10.1008, -1.0049], 0.0, 0.5943674189476091);

#atoms["Ga"] = makeatom("Ga", 30, 4, 3.0,  69.72,   3.0,   16,   [:s, :p], -414.75774465, [-8.8812, -2.4866], 0.0, 0.4697093810);

atoms["Ga"] = makeatom("Ga", 30, 4, 3.0,  69.72,   13.0,   6,   [:s,:d, :p], -414.75774465, [-8.8812,-19.0647, -2.4866], 0.0, 0.4697093810);

#atoms["Ge"] = makeatom("Ge", 31, 4, 4.0,  72.63,   4.0,    10,   [:s,:p, :d], -212.84384318, [-11.6643,  -3.8180, -1.0], 0.0, 0.5314231326);
#atoms["Ge"] = makeatom("Ge", 31, 4, 4.0,  72.63,   4.0,    10,   [:s,:d, :p], -212.84384318, [-11.6643, -1.0, -3.8180], 0.0, 0.5314231326);
atoms["Ge"] = makeatom("Ge", 31, 4, 4.0,  72.63,   4.0,    10,   [:s, :p], -212.84384318, [-11.6643, -3.8180], 0.0, 0.5314231326);

atoms["As"] = makeatom("As", 32, 4, 5.0,  74.92,   5.0,    0,   [:s, :p], -39.62501000, [-14.4349, -5.1279], 0.0,  0.5906811873);
atoms["Se"] = makeatom("Se", 33, 4, 6.0,  78.97,   6.0,    0,   [:s, :p],  -43.03651432, [-17.2230, -6.4409], 0.0, 0.6459585094);
atoms["Br"] = makeatom("Br", 34, 4, 7.0,  79.90,   7.0,    0,   [:s, :p], -40.38769338, [-20.0640, -7.7862], 0.0,  0.6992542568);

###

atoms["Rb"] = makeatom("Rb", 37, 5, 1, 85.468,  1.0,   8,   [:s, :d, :p], -53.05281238, [-2.0994, -0.3, -0.5921], 0.0, 0.282637);
#atoms["Rb"] = makeatom("Rb", 37, 5, 1, 85.468,  1.0,   8,   [:s, :p], -53.05281238, [-2.0994, -0.5921], 0.0, 0.282637);
atoms["Sr"] = makeatom("Sr", 38, 5, 2, 87.62,  2.0,   8,   [:s, :d, :p], -69.92604259, [-3.4017, -1.1543, -1.1318], 0.0,  0.365728);
                           
atoms["Y" ] = makeatom("Y" , 39, 5, 2.05, 88.906,  3.0,   8,    [:s, :d, :p], -92.30423365,  [-3.9864, -2.3907, -1.2867], 0.0,  0.5538339);
atoms["Zr"] = makeatom("Zr", 40, 5, 2.15, 91.224,  4.0,   8,    [:s, :d, :p], -98.57421435,  [-4.3347, -3.4817, -1.3201], 0.0,  0.4905697);
atoms["Nb"] = makeatom("Nb", 41, 5, 2.25, 92.906,  5.0,    8,   [:s, :d, :p], -117.19066274 ,  [-4.5053, -4.3124, -1.3159],0.0, 0.3883595);
atoms["Mo"] = makeatom("Mo", 42, 5, 2.35, 95.95,   6.0,    8,   [:s, :d, :p], -138.51328928,  [-4.2469, -4.2022, -1.0704], 0.0, 0.3888637);
atoms["Tc"] = makeatom("Tc", 43, 5, 2.45, 98.0,    7.0,    8,   [:s, :d, :p], -173.39188188,  [-3.9325, -3.9874, -0.7985], 0.0, 0.4606537);
atoms["Ru"] = makeatom("Ru", 44, 5, 2.55, 101.07,  8.0,    8,   [:s, :d, :p], -193.54935366,  [-3.5929, -3.7493, -0.4736], 0.0, 0.6120038);
atoms["Rh"] = makeatom("Rh", 45, 5, 2.65, 102.91,  9.0,    6,   [:s, :d, :p], -170.64161619,  [-3.3079, -3.6154, -0.2568], 0.0, 0.7119646);
atoms["Pd"] = makeatom("Pd", 46, 5, 2.75, 106.42,  10.0,   6,   [:s, :d, :p], -198.31288985,  [-3.2109, -3.9641, -0.1564], 0.0, 0.7660070);
atoms["Ag"] = makeatom("Ag", 47, 5, 2.85, 107.87,  11.0,   8,   [:s, :d, :p], -295.31434496,  [-4.3268, -7.3569, -0.5745], 0.0, 0.5027618);
atoms["Cd"] = makeatom("Cd", 48, 5, 2.95, 112.41,  12.0,   0,   [:s, :d, :p], -114.39800378,  [-5.6202, -11.5954, -1.0602], 0.0, 0.5507178);
                           
#atoms["In"] = makeatom("In", 49, 5, 3.0,  114.82,   3.0,   10,   [:s, :p], -132.97861648 , [-8.1551, -2.4193], 0.0, 0.4292669);
atoms["In"] = makeatom("In", 49, 5, 3.0,  114.82,   13.0,   0,   [:s,:d, :p], -132.97861648 , [-8.1551, -18.4828, -2.4193], 0.0, 0.4292669);
atoms["Sn"] = makeatom("Sn", 50, 5, 4.0,  118.71,   4.0,   10,   [:s, :p], -158.61463088, [-10.4877, -3.5995], 0.0,  0.4766457);
atoms["Sb"] = makeatom("Sb", 51, 5, 5.0,  121.76,   5.0,   10,   [:s, :p], -184.28336729, [-12.7656, -4.7332], 0.0,  0.5212358);
atoms["Te"] = makeatom("Te", 52, 5, 6.0,  127.60,   6.0,    0,   [:s, :p], -26.19141773, [-15.0608, -5.8767], 0.0,   0.5638125);
atoms["I" ] = makeatom("I" , 53, 5, 7.0,  126.90,   7.0,    0,   [:s, :p], -65.01087258 , [-17.3546, -7.0038], 0.0,  0.6039061);

###

atoms["Cs"] = makeatom("Cs", 55, 6, 1, 132.9,  1.0,   8,   [:s, :d, :p ], -62.83174426, [-1.9276, -0.3, -0.5431 ], 0.0, 0.262230);
#atoms["Cs"] = makeatom("Cs", 55, 6, 1, 132.9,  1.0,   8,   [:s, :p ], -62.83174426, [-1.9276, -0.5431 ], 0.0, 0.262230);
atoms["Ba"] = makeatom("Ba", 56, 6, 2, 139.3,  2.0,   8,   [:s, :d, :p], -70.01638694, [-3.0703, -1.8523, -1.0457], 0.0, 0.332378);
                           
#atoms["La"] = makeatom("La", 57, 6, 2.05, 138.9,  3.0,   8,    [:s, :d, :p], -101.82373826,  [-3.0750,-2.2784,-0.6575], 0.0, 0.5); #f orbitals?
atoms["La"] = makeatom("La", 57, 6, 2.05, 138.9,  3.0,   8,    [:s, :d], -101.82373826,  [-3.0750,-2.2784], 0.0, 0.5); #f orbitals?
atoms["Hf"] = makeatom("Hf", 72, 6, 2.15, 178.5,  4.0,   8,    [:s, :d, :p], -157.97552412,  [-5.0187, -2.6052, -1.2994], 0.0, 0.5723412); #lan
atoms["Ta"] = makeatom("Ta", 73, 6, 2.25, 181.0,  5.0,    8,   [:s, :d, :p], -141.00656328,  [-5.3481, -3.4993, -1.2418], 0.0, 0.6127437);
atoms["W" ] = makeatom("W" , 74, 6, 2.35, 183.8,  6.0,    8,   [:s, :d, :p], -158.16868802,  [-5.6245, -4.3914, -1.2092], 0.0, 0.6426146);
atoms["Re"] = makeatom("Re", 75, 6, 2.45, 186.2,  7.0,    8,   [:s, :d, :p], -189.48871022,  [-5.8618, -5.2763, -1.1642], 0.0, 0.5348510);
atoms["Os"] = makeatom("Os", 76, 6, 2.55, 190.2,  8.0,    8,   [:s, :d, :p], -197.71975092,  [-5.9152, -5.8238, -1.0547], 0.0, 0.4740821);
atoms["Ir"] = makeatom("Ir", 77, 6, 2.65, 192.2,  9.0,    6,   [:s, :d, :p], -180.96527345,  [-5.5876, -5.6416, -0.8041], 0.0, 0.4882183);
atoms["Pt"] = makeatom("Pt", 78, 6, 2.75, 195.1,  10.0,   6,   [:s, :d, :p], -210.06169163,  [-5.3377, -5.5187, -0.5975], 0.0, 0.5793362);
atoms["Au"] = makeatom("Au", 79, 6, 2.85, 197.0,  11.0,   0,   [:s, :d, :p], -114.46631173,  [-5.7155, -6.8176, -0.6483], 0.0, 0.5508305);
atoms["Hg"] = makeatom("Hg", 80, 6, 2.95, 200.6,  12.0,   0,   [:s, :d, :p], -107.48183748,  [-9.7901, -6.7757, -0.8974], 0.0, 0.5915364);
                           
atoms["Tl"] = makeatom("Tl", 81, 6, 3.0,  204.4,   13.0,   0,   [:s, :d, :p], -144.53037033, [ -9.4990,-15.3839,  -2.2651], 0.0, 0.4202546); #include d
#atoms["Tl"] = makeatom("Tl", 81, 6, 3.0,  204.4,   13.0,   0,   [:d, :s, :p], -144.53037033, [-15.3839, -9.4990,  -2.2651], 0.0, 0.4202546); #include d
#atoms["Tl"] = makeatom("Tl", 81, 6, 3.0,  204.4,   3.0,   10,   [:s, :p], -144.53037033, [-9.4990, -2.2651], 0.0, 0.4202546);

atoms["Pb"] = makeatom("Pb", 82, 6, 4.0,  207.2,   4.0,   10,   [:s, :p], -163.61141039, [-11.9754, -3.3964], 0.0, 0.4606776);
atoms["Bi"] = makeatom("Bi", 83, 6, 5.0,  209.0,   5.0,   10,   [:s, :p], -184.68674216, [-14.4080, -4.4700], 0.0, 0.4996092);



#U = IE - EA  (ionization energy - electron affinity

#stop printing
nothing

#approximate metallic radius / single covalent radius 
#for use in prerelaxing structures in pm
atom_radius = Dict()

atom_radius[ "X" ] =  200
atom_radius[ "Xa" ] =  200
atom_radius[ "T" ] =  200

atom_radius[ "H" ] =  50
atom_radius[ "He" ] =  50
atom_radius[ "Li" ] =  152
atom_radius[ "Be" ] =  112
atom_radius[ "B" ] =   82
atom_radius[ "C" ] =   77
atom_radius[ "N" ] =   75
atom_radius[ "O" ] =   73
atom_radius[ "F" ] =   71
atom_radius[ "Ne" ] =  80
atom_radius[ "Na" ] =  186
atom_radius[ "Mg" ] =  160
atom_radius[ "Al" ] =  143
atom_radius[ "Si" ] =  111
atom_radius[ "P" ] =   106
atom_radius[ "S" ] =   102
atom_radius[ "Cl" ] =  99
atom_radius[ "Ar" ] =  100
atom_radius[ "K" ] =   227
atom_radius[ "Ca" ] =  197
atom_radius[ "Sc" ] =  162
atom_radius[ "Ti" ] =  147
atom_radius[ "V" ] =   134
atom_radius[ "Cr" ] =  128
atom_radius[ "Mn" ] =  127
atom_radius[ "Fe" ] =  126
atom_radius[ "Co" ] =  125
atom_radius[ "Ni" ] =  124
atom_radius[ "Cu" ] =  128
atom_radius[ "Zn" ] =  134
atom_radius[ "Ga" ] =  135
atom_radius[ "Ge" ] =  122
atom_radius[ "As" ] =  119
atom_radius[ "Se" ] =  116
atom_radius[ "Br" ] =  114
atom_radius[ "Kr" ] =  100
atom_radius[ "Rb" ] =  248
atom_radius[ "Sr" ] =  215
atom_radius[ "Y" ] =  180
atom_radius[ "Zr" ] =  160
atom_radius[ "Nb" ] =  146
atom_radius[ "Mo" ] =  139
atom_radius[ "Tc" ] =  136
atom_radius[ "Ru" ] =  134
atom_radius[ "Rh" ] =  134
atom_radius[ "Pd" ] =  137
atom_radius[ "Ag" ] =  144
atom_radius[ "Cd" ] =  151
atom_radius[ "In" ] =  167
atom_radius[ "Sn" ] =  141
atom_radius[ "Sb" ] =  138
atom_radius[ "Te" ] =  135
atom_radius[ "I" ] =   133
atom_radius[ "Xe" ] =  130
atom_radius[ "Cs" ] =  265
atom_radius[ "Ba" ] =  222
atom_radius[ "La" ] =  187
atom_radius[ "Ce" ] =  181
atom_radius[ "Pr" ] =  181
atom_radius[ "Nd" ] =  181
atom_radius[ "Pm" ] =  181
atom_radius[ "Sm" ] =  181
atom_radius[ "Eu" ] =  181
atom_radius[ "Gd" ] =  181
atom_radius[ "Tb" ] =  181
atom_radius[ "Dy" ] =  181
atom_radius[ "Ho" ] =  181
atom_radius[ "Er" ] =  181
atom_radius[ "Tm" ] =  181
atom_radius[ "Yb" ] =  181
atom_radius[ "Lu" ] =  173
atom_radius[ "Hf" ] =  159
atom_radius[ "Ta" ] =  146
atom_radius[ "W" ] =  139
atom_radius[ "Re" ] =  137
atom_radius[ "Os" ] =  135
atom_radius[ "Ir" ] =  135
atom_radius[ "Pt" ] =  135
atom_radius[ "Au" ] =  144
atom_radius[ "Hg" ] =  151
atom_radius[ "Tl" ] =  170
atom_radius[ "Pb" ] =  147
atom_radius[ "Bi" ] =  146
atom_radius[ "Po" ] =  140
atom_radius[ "At" ] =  140


atom_prefered_oxidation = Dict()


atom_prefered_oxidation[ "H" ] =  [1, -1, 0]
atom_prefered_oxidation[ "He" ] = [0]
atom_prefered_oxidation[ "Li" ] = [1, 0]
atom_prefered_oxidation[ "Be" ] = [2, 0]
atom_prefered_oxidation[ "B" ] =  [3, 0]
atom_prefered_oxidation[ "C" ] =  [4, -4, 0]
atom_prefered_oxidation[ "N" ] =  [-3, 0]
atom_prefered_oxidation[ "O" ] =  [-2, 0]
atom_prefered_oxidation[ "F" ] =  [-1, 0]
atom_prefered_oxidation[ "Ne" ] = [0]
atom_prefered_oxidation[ "Na" ] = [1, 0]
atom_prefered_oxidation[ "Mg" ] = [2, 0]
atom_prefered_oxidation[ "Al" ] = [3, 0]
atom_prefered_oxidation[ "Si" ] = [4, -4, 2, 0]
atom_prefered_oxidation[ "P" ] =  [-3, 3, 5, 0]
atom_prefered_oxidation[ "S" ] =  [-2, 6, 0]
atom_prefered_oxidation[ "Cl" ] = [-1, 0]
atom_prefered_oxidation[ "Ar" ] = [0]
atom_prefered_oxidation[ "K" ] =  [1, 0]
atom_prefered_oxidation[ "Ca" ] = [2, 0]
atom_prefered_oxidation[ "Sc" ] = [3, 2, 0]
atom_prefered_oxidation[ "Ti" ] = [4, 3, 2, 0]
atom_prefered_oxidation[ "V" ] =  [5, 3, 0]
atom_prefered_oxidation[ "Cr" ] = [4, 3, 6, 2, 0]
atom_prefered_oxidation[ "Mn" ] = [4, 3, 2, 0]
atom_prefered_oxidation[ "Fe" ] = [4, 3, 2, 0]
atom_prefered_oxidation[ "Co" ] = [3, 4, 2, 0]
atom_prefered_oxidation[ "Ni" ] = [2, 4, 1, 0]
atom_prefered_oxidation[ "Cu" ] = [2, 4, 1, 0]
atom_prefered_oxidation[ "Zn" ] = [2, 0]
atom_prefered_oxidation[ "Ga" ] = [3, 0]
atom_prefered_oxidation[ "Ge" ] = [4, -4, 2, 0]
atom_prefered_oxidation[ "As" ] = [-3, 3, 5, 0]
atom_prefered_oxidation[ "Se" ] = [-2, 2, 4, 6, 0]
atom_prefered_oxidation[ "Br" ] = [-1, 0]
atom_prefered_oxidation[ "Kr" ] = [0]
atom_prefered_oxidation[ "Rb" ] = [1, 0]
atom_prefered_oxidation[ "Sr" ] = [2, 0]
atom_prefered_oxidation[ "Y" ] =  [3, 0]
atom_prefered_oxidation[ "Zr" ] = [4, 2, 0]
atom_prefered_oxidation[ "Nb" ] = [5, 0]
atom_prefered_oxidation[ "Mo" ] = [6, 4, 0]
atom_prefered_oxidation[ "Tc" ] = [4, 3, 2, 0]
atom_prefered_oxidation[ "Ru" ] = [4, 3, 2, 0]
atom_prefered_oxidation[ "Rh" ] = [4, 3, 2, 0]
atom_prefered_oxidation[ "Pd" ] = [4, 2, 0]
atom_prefered_oxidation[ "Ag" ] = [1, 2, 0]
atom_prefered_oxidation[ "Cd" ] = [2,0]
atom_prefered_oxidation[ "In" ] = [3,0]
atom_prefered_oxidation[ "Sn" ] = [4, 2, -4, 0]
atom_prefered_oxidation[ "Sb" ] = [-3, 3, 5, 0]
atom_prefered_oxidation[ "Te" ] = [-2, 6, 0]
atom_prefered_oxidation[ "I" ] =  [-1, 0]
atom_prefered_oxidation[ "Xe" ] = [0]
atom_prefered_oxidation[ "Cs" ] = [1, 0]
atom_prefered_oxidation[ "Ba" ] = [2, 0]
atom_prefered_oxidation[ "La" ] = [3, 0]
atom_prefered_oxidation[ "Ce" ] = [3, 0]
atom_prefered_oxidation[ "Pr" ] = [3, 0]
atom_prefered_oxidation[ "Nd" ] = [3, 0]
atom_prefered_oxidation[ "Pm" ] = [3, 0]
atom_prefered_oxidation[ "Sm" ] = [3, 0]
atom_prefered_oxidation[ "Eu" ] = [2, 3, 0]
atom_prefered_oxidation[ "Gd" ] = [3, 0]
atom_prefered_oxidation[ "Tb" ] = [3, 0]
atom_prefered_oxidation[ "Dy" ] = [3, 0]
atom_prefered_oxidation[ "Ho" ] = [3, 0]
atom_prefered_oxidation[ "Er" ] = [3, 0]
atom_prefered_oxidation[ "Tm" ] = [3, 0]
atom_prefered_oxidation[ "Yb" ] = [3, 0]
atom_prefered_oxidation[ "Lu" ] = [3, 0]
atom_prefered_oxidation[ "Hf" ] = [4, 2, 0]
atom_prefered_oxidation[ "Ta" ] = [5, 3, 0]
atom_prefered_oxidation[ "W" ] =  [6, 4, 0]
atom_prefered_oxidation[ "Re" ] = [4, 3, 6, 0]
atom_prefered_oxidation[ "Os" ] = [4, 3, 0]
atom_prefered_oxidation[ "Ir" ] = [4, 3, 0]
atom_prefered_oxidation[ "Pt" ] = [2, 4, 0]
atom_prefered_oxidation[ "Au" ] = [3, 1, 0]
atom_prefered_oxidation[ "Hg" ] = [2, 1, 0]
atom_prefered_oxidation[ "Tl" ] = [3, 1, 0]
atom_prefered_oxidation[ "Pb" ] = [4, 2, -4, 0]
atom_prefered_oxidation[ "Bi" ] = [3, -3, 0]
atom_prefered_oxidation[ "Po" ] = [-2, 0]
atom_prefered_oxidation[ "At" ] = [-1, 0]


cutoff_dist = Dict()

function get_cutoff(at1, at2)

    if (at1,at2) in keys(cutoff_dist)
        return cutoff_dist[(at1,at2)]
    else
        rad1 = atom_radius[at1] / 0.529 / 100.0
        rad2 = atom_radius[at2] / 0.529 / 100.0
        
        cutoff2X    = (rad1 + rad2) / 2.0 * 7.0
        cutoff_onX  = (rad1 + rad2) / 2.0 * 6.0

        cutoff2X   = max(min(cutoff2X,   19.01), 15.01) #2body
        cutoff_onX = max(min(cutoff_onX, 18.01), 14.51) #onsite

        cutoff_dist[(at1,at2)] = [cutoff2X, cutoff_onX]
        cutoff_dist[(at2,at1)] = [cutoff2X, cutoff_onX]
        
        return [cutoff2X, cutoff_onX]
    end
end

function get_cutoff(at1, at2, at3)
    if (at1,at2,at3) in keys(cutoff_dist)
        return cutoff_dist[(at1,at2,at3)]
    else
        rad1 = atom_radius[at1] / 0.529 / 100.0
        rad2 = atom_radius[at2] / 0.529 / 100.0
        rad3 = atom_radius[at3] / 0.529 / 100.0

        cutoff3bX = minimum([rad1, rad2, rad3])*5.0
#        cutoff3bX    = (rad1 + rad2 + rad3) / 3.0 * 5.0
        cutoff3bX  = max(min(cutoff3bX,  13.51),  10.01) #3body

        cutoff_dist[(at1,at2,at3)] = cutoff3bX
        cutoff_dist[(at1,at3,at2)] = cutoff3bX
        cutoff_dist[(at2,at1,at3)] = cutoff3bX
        cutoff_dist[(at2,at3,at1)] = cutoff3bX
        cutoff_dist[(at3,at1,at2)] = cutoff3bX
        cutoff_dist[(at3,at2,at1)] = cutoff3bX
        
        return cutoff3bX
    end
end        

min_dimer_dist_dict = Dict()

min_dimer_dist_dict["Pd"] =  3.73246
min_dimer_dist_dict["Si"] =  3.33537
min_dimer_dist_dict["C"] =  1.92837
min_dimer_dist_dict["P"] =  2.86897
min_dimer_dist_dict["Nb"] =  3.10416
min_dimer_dist_dict["Ag"] =  3.8093
min_dimer_dist_dict["Ru"] =  2.95446
min_dimer_dist_dict["Sb"] =  3.75256
min_dimer_dist_dict["Cs"] =  4.49968
min_dimer_dist_dict["Be"] =  3.2197
min_dimer_dist_dict["Sr"] =  3.95284
min_dimer_dist_dict["Ta"] =  3.2543
min_dimer_dist_dict["Ga"] =  3.94627
min_dimer_dist_dict["Te"] =  3.89208
min_dimer_dist_dict["Y"] =  4.13392
min_dimer_dist_dict["Ca"] =  5.03584
min_dimer_dist_dict["K"] =  4.29261
min_dimer_dist_dict["Au"] =  3.74318
min_dimer_dist_dict["V"] =  2.59038
min_dimer_dist_dict["Mo"] =  2.88349
min_dimer_dist_dict["Br"] =  3.46016
min_dimer_dist_dict["Zn"] =  3.88444
min_dimer_dist_dict["Sc"] =  3.49735
min_dimer_dist_dict["H"] =  1.16274
min_dimer_dist_dict["Al"] =  3.9947
min_dimer_dist_dict["S"] =  2.87718
min_dimer_dist_dict["Ge"] =  3.53994
min_dimer_dist_dict["Fe"] =  2.68548
min_dimer_dist_dict["Mg"] =  4.609
min_dimer_dist_dict["F"] =  2.12195
min_dimer_dist_dict["La"] =  3.92813
min_dimer_dist_dict["W"] =  3.03842
min_dimer_dist_dict["Li"] =  2.88861
min_dimer_dist_dict["O"] =  1.84197
min_dimer_dist_dict["B"] =  2.47166
min_dimer_dist_dict["Mn"] =  2.47598
min_dimer_dist_dict["Bi"] =  3.96946
min_dimer_dist_dict["Re"] =  3.06043
min_dimer_dist_dict["Hf"] =  3.71755
min_dimer_dist_dict["Cl"] =  3.00634
min_dimer_dist_dict["In"] =  4.58286
min_dimer_dist_dict["Rb"] =  4.50733
min_dimer_dist_dict["Se"] =  3.29481
min_dimer_dist_dict["Ba"] =  4.10738
min_dimer_dist_dict["Rh"] =  3.20676
min_dimer_dist_dict["Ti"] =  2.86562
min_dimer_dist_dict["Zr"] =  3.39044
min_dimer_dist_dict["Cr"] =  2.39985
min_dimer_dist_dict["Cu"] =  3.29852
min_dimer_dist_dict["Tc"] =  2.90788
min_dimer_dist_dict["Tl"] =  4.8252
min_dimer_dist_dict["Hg"] =  4.27212
min_dimer_dist_dict["N"] =  1.67319
min_dimer_dist_dict["Na"] =  3.50215
min_dimer_dist_dict["Os"] =  3.09582
min_dimer_dist_dict["Ni"] =  3.11813
min_dimer_dist_dict["Co"] =  2.91763
min_dimer_dist_dict["Pt"] =  3.47029
min_dimer_dist_dict["Ir"] =  3.25225
min_dimer_dist_dict["As"] =  3.18345
min_dimer_dist_dict["Sn"] =  4.14125
min_dimer_dist_dict["I"] =  4.02336
min_dimer_dist_dict["Cd"] =  4.29162
min_dimer_dist_dict["Pb"] =  4.35765




end #ends module
 
