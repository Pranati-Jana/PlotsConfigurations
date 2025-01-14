DYstitching = {}

# "NWevents" are total number of events in the samples before ANY cuts, all weighted by 'genWeight'.
# "Nevents" are total number of unweighted events in the samples before ANY cuts. (But these aren't needed for stitching)
# "GenW" is average ABSOLUTE value of all 'genWeight'. Ideally, the 'genWeight' is just this value *1 or *(-1), but in 2018 this is only the case for most events, not all of them.
DYstitching["2016"] = {}
DYstitching["2016"]["NWevents"] = {}
DYstitching["2016"]["NWevents"]["incl"] = 1897141621849.12 # _ext2
DYstitching["2016"]["NWevents"]["HT_70-100"] = 9691660.0
DYstitching["2016"]["NWevents"]["HT_100-200"] = 11017086.0 # nominal + _ext1
DYstitching["2016"]["NWevents"]["HT_200-400"] = 9609137.0 # nominal + _ext1
DYstitching["2016"]["NWevents"]["HT_400-600"] = 9725661.0 # nominal + _ext1
DYstitching["2016"]["NWevents"]["HT_600-800"] = 8292957.0
DYstitching["2016"]["NWevents"]["HT_800-1200"] = 2673066.0
DYstitching["2016"]["NWevents"]["HT_1200-2500"] = 596079.0
DYstitching["2016"]["NWevents"]["HT_2500-Inf"] = 399492.0
DYstitching["2016"]["NWevents"]["MLL_100-200"] = 677564380.018
DYstitching["2016"]["NWevents"]["MLL_200-400"] = 7015716.51774
DYstitching["2016"]["NWevents"]["MLL_400-500"] = 368792.490469
DYstitching["2016"]["NWevents"]["MLL_500-700"] = 211573.519843
DYstitching["2016"]["NWevents"]["MLL_700-800"] = 33038.6888613
DYstitching["2016"]["NWevents"]["MLL_800-1000"] = 27995.4083682
DYstitching["2016"]["NWevents"]["MLL_1000-1500"] = 14522.590115
DYstitching["2016"]["NWevents"]["MLL_1500-2000"] = 1989.39385564
DYstitching["2016"]["NWevents"]["MLL_2000-3000"] = 473.865105087
DYstitching["2016"]["NWevents"]["PTLL_50-100"] = 105391679401.22 # nominal + _ext3
DYstitching["2016"]["NWevents"]["PTLL_100-250"] = 16042219725.78 # nominal + _ext1 + _ext2 + _ext5
DYstitching["2016"]["NWevents"]["PTLL_250-400"] = 148903207.691 # nominal + _ext1 + _ext2 + _ext5
DYstitching["2016"]["NWevents"]["PTLL_400-650"] = 1458753.03603 # nominal + _ext1 + _ext2
DYstitching["2016"]["NWevents"]["PTLL_650-Inf"] = 137769.984046 # nominal + _ext1 + _ext2
#DYstitching["2016"]["Nevents"] = {}
#DYstitching["2016"]["Nevents"]["incl"] = 120777245.0 # _ext2
#DYstitching["2016"]["Nevents"]["HT_70-100"] = 9691660.0
#DYstitching["2016"]["Nevents"]["HT_100-200"] = 11017086.0 # nominal + _ext1
#DYstitching["2016"]["Nevents"]["HT_200-400"] = 9609137.0 # nominal + _ext1
#DYstitching["2016"]["Nevents"]["HT_400-600"] = 9725661.0 # nominal + _ext1
#DYstitching["2016"]["Nevents"]["HT_600-800"] = 8292957.0
#DYstitching["2016"]["Nevents"]["HT_800-1200"] = 2673066.0
#DYstitching["2016"]["Nevents"]["HT_1200-2500"] = 596079.0
#DYstitching["2016"]["Nevents"]["HT_2500-Inf"] = 399492.0
#DYstitching["2016"]["Nevents"]["MLL_100-200"] = 1083606.0
#DYstitching["2016"]["Nevents"]["MLL_200-400"] = 298679.0
#DYstitching["2016"]["Nevents"]["MLL_400-500"] = 287262.0
#DYstitching["2016"]["Nevents"]["MLL_500-700"] = 280940.0
#DYstitching["2016"]["Nevents"]["MLL_700-800"] = 276235.0
#DYstitching["2016"]["Nevents"]["MLL_800-1000"] = 271768.0
#DYstitching["2016"]["Nevents"]["MLL_1000-1500"] = 258620.0
#DYstitching["2016"]["Nevents"]["MLL_1500-2000"] = 258625.0
#DYstitching["2016"]["Nevents"]["MLL_2000-3000"] = 255342.0
#DYstitching["2016"]["Nevents"]["PTLL_50-100"] = 115343756.0 # nominal + _ext3
#DYstitching["2016"]["Nevents"]["PTLL_100-250"] = 80736664.0 # nominal + _ext1 + _ext2 + _ext5
#DYstitching["2016"]["Nevents"]["PTLL_250-400"] = 20669993.0 # nominal + _ext1 + _ext2 + _ext5
#DYstitching["2016"]["Nevents"]["PTLL_400-650"] = 1625936.0 # nominal + _ext1 + _ext2
#DYstitching["2016"]["Nevents"]["PTLL_650-Inf"] = 1627882.0 # nominal + _ext1 + _ext2
DYstitching["2016"]["XS"] = {}
DYstitching["2016"]["XS"]["incl"] = 6025.2
DYstitching["2016"]["XS"]["HT_70-100"] = 209.1280411
DYstitching["2016"]["XS"]["HT_100-200"] = 181.4330386
DYstitching["2016"]["XS"]["HT_200-400"] = 50.45414011
DYstitching["2016"]["XS"]["HT_400-600"] = 6.988987742
DYstitching["2016"]["XS"]["HT_600-800"] = 1.682625263
DYstitching["2016"]["XS"]["HT_800-1200"] = 0.7759524256
DYstitching["2016"]["XS"]["HT_1200-2500"] = 0.1863565946
DYstitching["2016"]["XS"]["HT_2500-Inf"] = 0.0043881193
DYstitching["2016"]["XS"]["MLL_100-200"] = 226.6
DYstitching["2016"]["XS"]["MLL_200-400"] = 7.77
DYstitching["2016"]["XS"]["MLL_400-500"] = 0.4065
DYstitching["2016"]["XS"]["MLL_500-700"] = 0.2334
DYstitching["2016"]["XS"]["MLL_700-800"] = 0.03614
DYstitching["2016"]["XS"]["MLL_800-1000"] = 0.03047
DYstitching["2016"]["XS"]["MLL_1000-1500"] = 0.01636
DYstitching["2016"]["XS"]["MLL_1500-2000"] = 0.00218
DYstitching["2016"]["XS"]["MLL_2000-3000"] = 0.0005156
DYstitching["2016"]["XS"]["PTLL_50-100"] = 354.8
DYstitching["2016"]["XS"]["PTLL_100-250"] = 81.22
DYstitching["2016"]["XS"]["PTLL_250-400"] = 2.991
DYstitching["2016"]["XS"]["PTLL_400-650"] = 0.3882
DYstitching["2016"]["XS"]["PTLL_650-Inf"] = 0.03737
DYstitching["2016"]["GenW"] = {}
DYstitching["2016"]["GenW"]["incl"] = 23443.4238281
DYstitching["2016"]["GenW"]["HT_70-100"] = 1.0
DYstitching["2016"]["GenW"]["HT_100-200"] = 1.0
DYstitching["2016"]["GenW"]["HT_200-400"] = 1.0
DYstitching["2016"]["GenW"]["HT_400-600"] = 1.0
DYstitching["2016"]["GenW"]["HT_600-800"] = 1.0
DYstitching["2016"]["GenW"]["HT_800-1200"] = 1.0
DYstitching["2016"]["GenW"]["HT_1200-2500"] = 1.0
DYstitching["2016"]["GenW"]["HT_2500-Inf"] = 1.0
DYstitching["2016"]["GenW"]["MLL_100-200"] = 963.769104004
DYstitching["2016"]["GenW"]["MLL_200-400"] = 40.8800773621
DYstitching["2016"]["GenW"]["MLL_400-500"] = 2.43926501274
DYstitching["2016"]["GenW"]["MLL_500-700"] = 1.46828174591
DYstitching["2016"]["GenW"]["MLL_700-800"] = 0.241346806288
DYstitching["2016"]["GenW"]["MLL_800-1000"] = 0.21275369823
DYstitching["2016"]["GenW"]["MLL_1000-1500"] = 0.121011503041
DYstitching["2016"]["GenW"]["MLL_1500-2000"] = 0.0178087167442
DYstitching["2016"]["GenW"]["MLL_2000-3000"] = 0.00467793131247
DYstitching["2016"]["GenW"]["PTLL_50-100"] = 2507.18735146
DYstitching["2016"]["GenW"]["PTLL_100-250"] = 557.085410719
DYstitching["2016"]["GenW"]["PTLL_250-400"] = 19.6042160161
DYstitching["2016"]["GenW"]["PTLL_400-650"] = 2.32259188016
DYstitching["2016"]["GenW"]["PTLL_650-Inf"] = 0.206063707404

DYstitching["2017"] = {}
DYstitching["2017"]["NWevents"] = {}
DYstitching["2017"]["NWevents"]["incl"] = 3676449948155.45 # nominal + _ext1
DYstitching["2017"]["NWevents"]["HT_70-100"] = 9333543.0
DYstitching["2017"]["NWevents"]["HT_100-200"] = 15124171.0 # _newpmx + _ext1
DYstitching["2017"]["NWevents"]["HT_200-400"] = 11896758.0 # nominal + _ext1
DYstitching["2017"]["NWevents"]["HT_400-600"] = 11294006.0 # _newpmx + _ext1
DYstitching["2017"]["NWevents"]["HT_600-800"] = 8691608.0
DYstitching["2017"]["NWevents"]["HT_800-1200"] = 3089712.0
DYstitching["2017"]["NWevents"]["HT_1200-2500"] = 616923.0
DYstitching["2017"]["NWevents"]["HT_2500-Inf"] = 401334.0
DYstitching["2017"]["NWevents"]["MLL_100-200"] = 10527337145.4
DYstitching["2017"]["NWevents"]["MLL_200-400"] = 80609420.4812
DYstitching["2017"]["NWevents"]["MLL_400-500"] = 762279.972199
DYstitching["2017"]["NWevents"]["MLL_500-700"] = 410961.374568
DYstitching["2017"]["NWevents"]["MLL_700-800"] = 67318.4698278
DYstitching["2017"]["NWevents"]["MLL_800-1000"] = 61595.4142462
DYstitching["2017"]["NWevents"]["MLL_1000-1500"] = 27023.7381562
DYstitching["2017"]["NWevents"]["MLL_1500-2000"] = 3650.43573602
DYstitching["2017"]["NWevents"]["MLL_2000-3000"] = 907.251873148
DYstitching["2017"]["NWevents"]["MLL_3000-Inf"] = 46.752888016
DYstitching["2017"]["NWevents"]["PTLL_0-50"] = 11674084360009.28
DYstitching["2017"]["NWevents"]["PTLL_50-100"] = 89926628391.3
DYstitching["2017"]["NWevents"]["PTLL_100-250"] = 13550028759.0
DYstitching["2017"]["NWevents"]["PTLL_250-400"] = 172400451.16
DYstitching["2017"]["NWevents"]["PTLL_400-650"] = 1960515.32297
DYstitching["2017"]["NWevents"]["PTLL_650-Inf"] = 197651.309224
#DYstitching["2017"]["Nevents"] = {}
#DYstitching["2017"]["Nevents"]["incl"] = 206031656.0 # nominal + _ext1
#DYstitching["2017"]["Nevents"]["HT_70-100"] = 9344037.0
#DYstitching["2017"]["Nevents"]["HT_100-200"] = 15147827.0 # _newpmx + _ext1
#DYstitching["2017"]["Nevents"]["HT_200-400"] = 11929310.0 # nominal + _ext1
#DYstitching["2017"]["Nevents"]["HT_400-600"] = 11343818.0 # _newpmx + _ext1
#DYstitching["2017"]["Nevents"]["HT_600-800"] = 8743640.0
#DYstitching["2017"]["Nevents"]["HT_800-1200"] = 3114980.0
#DYstitching["2017"]["Nevents"]["HT_1200-2500"] = 625517.0
#DYstitching["2017"]["Nevents"]["HT_2500-Inf"] = 419308.0
#DYstitching["2017"]["Nevents"]["MLL_100-200"] = 13787181.0
#DYstitching["2017"]["Nevents"]["MLL_200-400"] = 2893179.0
#DYstitching["2017"]["Nevents"]["MLL_400-500"] = 504582.0
#DYstitching["2017"]["Nevents"]["MLL_500-700"] = 478062.0
#DYstitching["2017"]["Nevents"]["MLL_700-800"] = 503757.0
#DYstitching["2017"]["Nevents"]["MLL_800-1000"] = 540993.0
#DYstitching["2017"]["Nevents"]["MLL_1000-1500"] = 444230.0
#DYstitching["2017"]["Nevents"]["MLL_1500-2000"] = 458896.0
#DYstitching["2017"]["Nevents"]["MLL_2000-3000"] = 502544.0
#DYstitching["2017"]["Nevents"]["MLL_3000-Inf"] = 498600.0
#DYstitching["2017"]["Nevents"]["PTLL_0-50"] = 94420459.0
#DYstitching["2017"]["Nevents"]["PTLL_50-100"] = 90629617.0
#DYstitching["2017"]["Nevents"]["PTLL_100-250"] = 61882913.0
#DYstitching["2017"]["Nevents"]["PTLL_250-400"] = 20970944.0
#DYstitching["2017"]["Nevents"]["PTLL_400-650"] = 1766556.0
#DYstitching["2017"]["Nevents"]["PTLL_650-Inf"] = 1885153.0
DYstitching["2017"]["XS"] = {}
DYstitching["2017"]["XS"]["incl"] = 6189.39
DYstitching["2017"]["XS"]["HT_70-100"] = 196.7442
DYstitching["2017"]["XS"]["HT_100-200"] = 186.5538
DYstitching["2017"]["XS"]["HT_200-400"] = 56.34828
DYstitching["2017"]["XS"]["HT_400-600"] = 8.068944
DYstitching["2017"]["XS"]["HT_600-800"] = 2.018394
DYstitching["2017"]["XS"]["HT_800-1200"] = 0.9324216
DYstitching["2017"]["XS"]["HT_1200-2500"] = 0.2238414
DYstitching["2017"]["XS"]["HT_2500-Inf"] = 0.004015944
DYstitching["2017"]["XS"]["MLL_100-200"] = 247.8
DYstitching["2017"]["XS"]["MLL_200-400"] = 8.502
DYstitching["2017"]["XS"]["MLL_400-500"] = 0.4514
DYstitching["2017"]["XS"]["MLL_500-700"] = 0.2558
DYstitching["2017"]["XS"]["MLL_700-800"] = 0.04023
DYstitching["2017"]["XS"]["MLL_800-1000"] = 0.03406
DYstitching["2017"]["XS"]["MLL_1000-1500"] = 0.01828
DYstitching["2017"]["XS"]["MLL_1500-2000"] = 0.002367
DYstitching["2017"]["XS"]["MLL_2000-3000"] = 0.0005409
DYstitching["2017"]["XS"]["MLL_3000-Inf"] = 0.00003048
DYstitching["2017"]["XS"]["PTLL_0-50"] = 106300.0
DYstitching["2017"]["XS"]["PTLL_50-100"] = 407.9
DYstitching["2017"]["XS"]["PTLL_100-250"] = 96.8
DYstitching["2017"]["XS"]["PTLL_250-400"] = 3.774
DYstitching["2017"]["XS"]["PTLL_400-650"] = 0.5164
DYstitching["2017"]["XS"]["PTLL_650-Inf"] = 0.04796
DYstitching["2017"]["GenW"] = {}
DYstitching["2017"]["GenW"]["incl"] = 26331.2011719
DYstitching["2017"]["GenW"]["HT_70-100"] = 1.0
DYstitching["2017"]["GenW"]["HT_100-200"] = 1.0
DYstitching["2017"]["GenW"]["HT_200-400"] = 1.0
DYstitching["2017"]["GenW"]["HT_400-600"] = 1.0
DYstitching["2017"]["GenW"]["HT_600-800"] = 1.0
DYstitching["2017"]["GenW"]["HT_800-1200"] = 1.0
DYstitching["2017"]["GenW"]["HT_1200-2500"] = 1.0
DYstitching["2017"]["GenW"]["HT_2500-Inf"] = 1.0
DYstitching["2017"]["GenW"]["MLL_100-200"] = 1081.3458252
DYstitching["2017"]["GenW"]["MLL_200-400"] = 43.2147254944
DYstitching["2017"]["GenW"]["MLL_400-500"] = 2.39647388458
DYstitching["2017"]["GenW"]["MLL_500-700"] = 1.38987624645
DYstitching["2017"]["GenW"]["MLL_700-800"] = 0.217704832554
DYstitching["2017"]["GenW"]["MLL_800-1000"] = 0.185504332185
DYstitching["2017"]["GenW"]["MLL_1000-1500"] = 0.0990991294384
DYstitching["2017"]["GenW"]["MLL_1500-2000"] = 0.0131418416277
DYstitching["2017"]["GenW"]["MLL_2000-3000"] = 0.00318269222043
DYstitching["2017"]["GenW"]["MLL_3000-Inf"] = 0.000275370111922
DYstitching["2017"]["GenW"]["PTLL_0-50"] = 125423.01519
DYstitching["2017"]["GenW"]["PTLL_50-100"] = 2720.92335987
DYstitching["2017"]["GenW"]["PTLL_100-250"] = 599.822883664
DYstitching["2017"]["GenW"]["PTLL_250-400"] = 21.2592097076
DYstitching["2017"]["GenW"]["PTLL_400-650"] = 2.71139467184
DYstitching["2017"]["GenW"]["PTLL_650-Inf"] = 0.239168909401

DYstitching["2018"] = {}
DYstitching["2018"]["NWevents"] = {}
DYstitching["2018"]["NWevents"]["incl"] = 3462175967212.74 # nominal + _ext2
DYstitching["2018"]["NWevents"]["HT_70-100"] = 10010341.4025
DYstitching["2018"]["NWevents"]["HT_100-200"] = 11516745.8532
DYstitching["2018"]["NWevents"]["HT_200-400"] = 10840078.364
DYstitching["2018"]["NWevents"]["HT_400-600"] = 18902490.5772 # nominal + _ext1
DYstitching["2018"]["NWevents"]["HT_600-800"] = 8826238.14915
DYstitching["2018"]["NWevents"]["HT_800-1200"] = 3120982.10607
DYstitching["2018"]["NWevents"]["HT_1200-2500"] = 531566.862752
DYstitching["2018"]["NWevents"]["HT_2500-Inf"] = 415517.021934
DYstitching["2018"]["NWevents"]["MLL_100-200"] = 10869479189.4
DYstitching["2018"]["NWevents"]["MLL_200-400"] = 72226325.7221
DYstitching["2018"]["NWevents"]["MLL_400-500"] = 735637.201303
DYstitching["2018"]["NWevents"]["MLL_500-700"] = 416537.715495
DYstitching["2018"]["NWevents"]["MLL_700-800"] = 59976.059433
DYstitching["2018"]["NWevents"]["MLL_800-1000"] = 57465.163306
DYstitching["2018"]["NWevents"]["MLL_1000-1500"] = 30049.9509265
DYstitching["2018"]["NWevents"]["MLL_1500-2000"] = 4330.5045304
DYstitching["2018"]["NWevents"]["MLL_2000-3000"] = 899.311337645
DYstitching["2018"]["NWevents"]["MLL_3000-Inf"] = 51.0181199711
DYstitching["2018"]["NWevents"]["PTLL_0-50"] = 8792621513619.46
DYstitching["2018"]["NWevents"]["PTLL_50-100"] = 119578608420.74
DYstitching["2018"]["NWevents"]["PTLL_100-250"] = 17888382439.0
DYstitching["2018"]["NWevents"]["PTLL_250-400"] = 165953665.353
DYstitching["2018"]["NWevents"]["PTLL_400-650"] = 2070432.0925
DYstitching["2018"]["NWevents"]["PTLL_650-Inf"] = 198926.018914
#DYstitching["2018"]["Nevents"] = {}
#DYstitching["2018"]["Nevents"]["incl"] = 194117151.0 # nominal + _ext2
#DYstitching["2018"]["Nevents"]["HT_70-100"] = 10019684.0
#DYstitching["2018"]["Nevents"]["HT_100-200"] = 11530510.0
#DYstitching["2018"]["Nevents"]["HT_200-400"] = 10860690.0
#DYstitching["2018"]["Nevents"]["HT_400-600"] = 18959400.0 # nominal + _ext1
#DYstitching["2018"]["Nevents"]["HT_600-800"] = 8862104.0
#DYstitching["2018"]["Nevents"]["HT_800-1200"] = 3138129.0
#DYstitching["2018"]["Nevents"]["HT_1200-2500"] = 536416.0
#DYstitching["2018"]["Nevents"]["HT_2500-Inf"] = 427051.0
#DYstitching["2018"]["Nevents"]["MLL_100-200"] = 14255862.0
#DYstitching["2018"]["Nevents"]["MLL_200-400"] = 2587941.0
#DYstitching["2018"]["Nevents"]["MLL_400-500"] = 487372.0
#DYstitching["2018"]["Nevents"]["MLL_500-700"] = 481980.0
#DYstitching["2018"]["Nevents"]["MLL_700-800"] = 446493.0
#DYstitching["2018"]["Nevents"]["MLL_800-1000"] = 504407.0
#DYstitching["2018"]["Nevents"]["MLL_1000-1500"] = 492743.0
#DYstitching["2018"]["Nevents"]["MLL_1500-2000"] = 541461.0
#DYstitching["2018"]["Nevents"]["MLL_2000-3000"] = 489512.0
#DYstitching["2018"]["Nevents"]["MLL_3000-Inf"] = 496959.0
#DYstitching["2018"]["Nevents"]["PTLL_0-50"] = 71114716.0
#DYstitching["2018"]["Nevents"]["PTLL_50-100"] = 120435658.0
#DYstitching["2018"]["Nevents"]["PTLL_100-250"] = 81687671.0
#DYstitching["2018"]["Nevents"]["PTLL_250-400"] = 20199757.0
#DYstitching["2018"]["Nevents"]["PTLL_400-650"] = 1863809.0
#DYstitching["2018"]["Nevents"]["PTLL_650-Inf"] = 1895998.0
DYstitching["2018"]["XS"] = {}
DYstitching["2018"]["XS"]["incl"] = 6189.39
DYstitching["2018"]["XS"]["HT_70-100"] = 196.7442
DYstitching["2018"]["XS"]["HT_100-200"] = 186.5538
DYstitching["2018"]["XS"]["HT_200-400"] = 56.34828
DYstitching["2018"]["XS"]["HT_400-600"] = 8.068944
DYstitching["2018"]["XS"]["HT_600-800"] = 2.018394
DYstitching["2018"]["XS"]["HT_800-1200"] = 0.9324216
DYstitching["2018"]["XS"]["HT_1200-2500"] = 0.2238414
DYstitching["2018"]["XS"]["HT_2500-Inf"] = 0.004015944
DYstitching["2018"]["XS"]["MLL_100-200"] = 247.8
DYstitching["2018"]["XS"]["MLL_200-400"] = 8.502
DYstitching["2018"]["XS"]["MLL_400-500"] = 0.4514
DYstitching["2018"]["XS"]["MLL_500-700"] = 0.2558
DYstitching["2018"]["XS"]["MLL_700-800"] = 0.04023
DYstitching["2018"]["XS"]["MLL_800-1000"] = 0.03406
DYstitching["2018"]["XS"]["MLL_1000-1500"] = 0.01828
DYstitching["2018"]["XS"]["MLL_1500-2000"] = 0.002367
DYstitching["2018"]["XS"]["MLL_2000-3000"] = 0.0005409
DYstitching["2018"]["XS"]["MLL_3000-Inf"] = 0.00003048
DYstitching["2018"]["XS"]["PTLL_0-50"] = 106300.0
DYstitching["2018"]["XS"]["PTLL_50-100"] = 407.9
DYstitching["2018"]["XS"]["PTLL_100-250"] = 96.8
DYstitching["2018"]["XS"]["PTLL_250-400"] = 3.774
DYstitching["2018"]["XS"]["PTLL_400-650"] = 0.5164
DYstitching["2018"]["XS"]["PTLL_650-Inf"] = 0.04796
DYstitching["2018"]["GenW"] = {}
DYstitching["2018"]["GenW"]["incl"] = 26320.0368582
DYstitching["2018"]["GenW"]["HT_70-100"] = 0.999905441832
DYstitching["2018"]["GenW"]["HT_100-200"] = 0.999883631108
DYstitching["2018"]["GenW"]["HT_200-400"] = 0.999837144759
DYstitching["2018"]["GenW"]["HT_400-600"] = 0.999768702394
DYstitching["2018"]["GenW"]["HT_600-800"] = 0.999715788218
DYstitching["2018"]["GenW"]["HT_800-1200"] = 0.999643181447
DYstitching["2018"]["GenW"]["HT_1200-2500"] = 0.999521690168
DYstitching["2018"]["GenW"]["HT_2500-Inf"] = 0.999217616794
DYstitching["2018"]["GenW"]["MLL_100-200"] = 1080.81259633
DYstitching["2018"]["GenW"]["MLL_200-400"] = 43.182890384
DYstitching["2018"]["GenW"]["MLL_400-500"] = 2.39421090155
DYstitching["2018"]["GenW"]["MLL_500-700"] = 1.38846091102
DYstitching["2018"]["GenW"]["MLL_700-800"] = 0.217464477459
DYstitching["2018"]["GenW"]["MLL_800-1000"] = 0.185289886054
DYstitching["2018"]["GenW"]["MLL_1000-1500"] = 0.0989730613257
DYstitching["2018"]["GenW"]["MLL_1500-2000"] = 0.0131235878648
DYstitching["2018"]["GenW"]["MLL_2000-3000"] = 0.00317773212428
DYstitching["2018"]["GenW"]["MLL_3000-Inf"] = 0.000274766467068
DYstitching["2018"]["GenW"]["PTLL_0-50"] = 125423.2262
DYstitching["2018"]["GenW"]["PTLL_50-100"] = 2720.92898209
DYstitching["2018"]["GenW"]["PTLL_100-250"] = 599.823462922
DYstitching["2018"]["GenW"]["PTLL_250-400"] = 21.2591663687
DYstitching["2018"]["GenW"]["PTLL_400-650"] = 2.71138456518
DYstitching["2018"]["GenW"]["PTLL_650-Inf"] = 0.239170265942
