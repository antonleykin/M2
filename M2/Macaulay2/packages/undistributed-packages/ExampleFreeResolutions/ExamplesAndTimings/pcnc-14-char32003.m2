R1 = (ZZ/32003)[a, b, c, d, e, f, g, h, i, j, k, l, m]
I1 = ideal(
    -2879*a^2-7674*a*b-13341*b^2+2373*a*c+7926*b*c-478*c^2-4105*a*d-9188*b*d+14375*c*d+2531*d^2+14476*a*e-7975*b*e+5587*c*e-1528*d*e+3970*a*f-8284*b*f+15235*c*f-11832*a*g+4146*b*g+14201*c*g+5277*a*h-4090*b*h-572*a*i-6405*a*j+b*j-11747*a*k,
    -9457*a^2-7203*a*b+8257*b^2+13018*a*c-3748*b*c-5591*c^2-10836*a*d+3933*b*d-14528*c*d-5717*d^2+9085*a*e-11571*b*e-7112*c*e-1648*d*e+5612*a*f+4638*b*f+10702*c*f+8985*a*g-393*b*g-12581*c*g+8670*a*h+12549*b*h-13877*a*i+7262*a*j-13968*a*k+c*k+1371*a*l-7856*a*m,
    -10633*a^2+6610*a*b-1380*b^2+1568*a*c-13301*b*c-708*c^2-6168*a*d-1297*b*d-12406*c*d+12441*d^2-666*a*e+4975*b*e+8175*c*e-3068*d*e-9319*a*f-883*b*f-14785*c*f-4583*a*g-15653*b*g-3312*c*g-4634*a*h-7327*b*h+c*h+14847*a*i-2500*a*j,
    -556*a^2-7597*a*b-5258*b^2+1134*a*c+14432*b*c-11128*c^2+411*a*d+562*b*d+2518*c*d+5938*d^2+2025*a*e-13714*b*e-3167*c*e+9869*d*e+3798*a*f+4605*b*f-9128*c*f-15634*a*g-9047*b*g+6636*c*g-6428*a*h-1421*b*h+4866*a*i+15723*a*j+14817*a*k+2205*a*l+c*l-6462*b*m,
    -10413*a^2+2827*a*b+9174*b^2-2489*a*c+15823*b*c+14352*c^2-12171*a*d+3098*b*d+7547*c*d-15617*d^2+7994*a*e+13746*b*e-4274*c*e+6526*d*e+7090*a*f+8609*b*f+7121*c*f+10751*a*g-5899*b*g+15130*c*g-7817*a*h+15356*b*h-15357*a*i+c*i-14731*a*j+7629*a*k,
    -15312*a^2+15148*a*b-1745*b^2-6201*a*c-6889*b*c+1357*c^2-13453*a*d+4313*b*d-113*c*d+504*d^2-4911*a*e+3253*b*e+10796*c*e+5983*d*e-15471*a*f+5929*b*f+7211*c*f-711*a*g+10969*b*g-517*c*g-4913*a*h-1468*b*h-11070*a*i-13685*a*j-6405*a*k+b*k+5397*a*l,
    1159*a^2+2856*a*b-5385*b^2-7113*a*c-1324*b*c+7948*c^2-4708*a*d+8022*b*d-12259*c*d-412*d^2-1919*a*e+12914*b*e+7300*c*e+3303*d*e-5866*a*f+2696*b*f-4233*c*f-11064*a*g-9214*b*g+14240*c*g+8958*a*h+15069*b*h+12030*a*i-5868*a*j-2175*a*k-6405*a*l+b*l+13452*a*m,
    -12343*a^2-13480*a*b-15887*b^2+6880*a*c+11177*b*c+2377*c^2-1793*a*d+7280*b*d+15123*c*d-1177*d^2+6325*a*e+2910*b*e+6913*c*e+3246*d*e+6232*a*f+3326*b*f-8463*c*f-362*a*g-7786*b*g+2967*c*g+8334*a*h-1899*b*h-4279*a*i+b*i-15622*a*j,
    15352*a^2-13674*a*b+8685*b^2-9601*a*c+3191*b*c+15826*c^2-470*a*d+5860*b*d-10940*c*d-1783*d^2+5224*a*e-13945*b*e+12903*c*e+1037*d*e+13894*a*f-498*b*f-2394*c*f-1340*a*g-4353*b*g-2940*c*g+12554*a*h+6324*b*h-2784*a*i-13968*a*j+c*j+5199*a*k+4643*a*l,
    2133*a^2+6966*a*b+685*b^2+15995*a*c-7125*b*c-15877*c^2+6471*a*d-6776*b*d-2818*c*d-2553*d^2+5086*a*e+4184*b*e-12665*c*e-5038*d*e-8827*a*f-14426*b*f-3126*c*f+d*f+8165*a*g+9639*b*g-15061*c*g-13450*a*h+8917*b*h+7497*a*i,
    -7941*a^2+12454*a*b+3296*b^2-8486*a*c+1960*b*c+4957*c^2-1195*a*d+1332*b*d+8762*c*d-564*d^2+10486*a*e+7067*b*e-3762*c*e-3373*d*e+14004*a*f+14232*b*f-10403*c*f-2287*a*g-14924*b*g-2055*c*g+d*g-11006*a*h+15322*b*h+15858*a*i-11732*a*j,
    -2930*a^2-13584*a*b-1310*b^2-8161*a*c-12719*b*c-7422*c^2-15714*a*d+4472*b*d+1554*c*d-14628*d^2+7386*a*e+15877*b*e-5843*c*e+4692*d*e-7889*a*f-9487*b*f-979*c*f-9631*a*g-7463*b*g+10021*c*g+2922*a*h+3688*b*h+d*h-5038*a*i+8660*a*j-481*a*k,
    13970*a^2+9163*a*b+11724*b^2+15268*a*c+4895*b*c-8999*c^2-15083*a*d-5017*b*d-8000*c*d-4110*d^2+4032*a*e+3074*b*e-2090*c*e-4086*d*e-11872*a*f+11616*b*f+3135*c*f-2253*a*g+5940*b*g+12828*c*g-3869*a*h+4897*b*h+2569*a*i+d*i-8989*a*j+4040*a*k-781*a*l,
    -8098*a^2+247*a*b+14786*b^2+2928*a*c+8224*b*c+5596*c^2+10188*a*d-12692*b*d-246*c*d-9357*d^2+7597*a*e-10946*b*e+5053*c*e-115*d*e+6071*a*f-7210*b*f-9158*c*f+13766*a*g+13692*b*g+13074*c*g-4736*a*h-4889*b*h-6466*a*i-5850*a*j+d*j+4174*a*k+9007*a*l-2638*a*m,
    8984*a^2+10921*a*b+8720*b^2+11572*a*c-7793*b*c+8312*c^2-12524*a*d-11954*b*d-8818*c*d-14337*d^2+3848*a*e-1362*b*e-11453*c*e-4058*d*e-6768*a*f+13294*b*f-3383*c*f-1444*a*g-13926*b*g-7880*c*g+2515*a*h+2573*b*h-5189*a*i+6472*a*j-10969*a*k+d*k-10067*a*l-9293*b*m,
    -1631*a^2-7375*a*b-4469*b^2-13023*a*c+4026*b*c+9639*c^2-12712*a*d-2042*b*d+7238*c*d-295*d^2+4363*a*e-6916*b*e+1823*c*e+9951*d*e-4128*a*f-1450*b*f-926*c*f+9697*a*g+3560*b*g+5722*c*g-8085*a*h-5579*b*h+2976*a*i+2919*a*j+7663*a*k+14182*a*l+d*l-15129*c*m,
    1170*a^2-14687*a*b-13125*b^2-10928*a*c+12290*b*c+182*c^2-13551*a*d+10074*b*d-482*c*d+12834*d^2+5464*a*e-14748*b*e-2511*c*e+938*d*e+e^2+6756*a*f-14560*b*f+5452*c*f+14885*a*g+3080*b*g+6227*c*g+4172*a*h+14105*b*h-13038*a*i,
    -597*a^2+11739*a*b-12010*b^2+13835*a*c-5519*b*c-9449*c^2+2583*a*d-10826*b*d-7836*c*d+7093*d^2-7032*a*e+6254*b*e-5121*c*e-3381*d*e-8606*a*f+4256*b*f+798*c*f+e*f+3963*a*g-5983*b*g+102*c*g+2469*a*h-15777*b*h+11848*a*i-4807*a*j,
    9277*a^2-14080*a*b-3252*b^2+3656*a*c+10856*b*c+3280*c^2+1179*a*d+4279*b*d+8566*c*d-4466*d^2+15824*a*e+15331*b*e-1898*c*e+457*d*e-9804*a*f-4890*b*f-12325*c*f-10693*a*g+1221*b*g-2987*c*g+e*g+14866*a*h-12442*b*h+15424*a*i+2645*a*j+1978*a*k,
    -13147*a^2+14570*a*b-4208*b^2+9486*a*c+10666*b*c+8186*c^2+3214*a*d-4137*b*d-11792*c*d-4388*d^2-13139*a*e+11799*b*e+11942*c*e-5464*d*e+9779*a*f+15273*b*f+5489*c*f+1069*a*g-3911*b*g-14798*c*g-15715*a*h-4905*b*h+e*h+2160*a*i-2488*a*j-12546*a*k+7004*a*l,
    11567*a^2-440*a*b+14919*b^2-11551*a*c-10954*b*c+5814*c^2-13779*a*d+2712*b*d+15775*c*d+1224*d^2-9777*a*e-4292*b*e-11474*c*e-12296*d*e-4030*a*f+3475*b*f-15106*c*f-14256*a*g-12408*b*g-10123*c*g-15822*a*h+11374*b*h+9137*a*i+e*i-1298*a*j-5004*a*k-1076*a*l+13682*a*m,
    10829*a^2-15045*a*b-9312*b^2+12111*a*c-76*b*c-7894*c^2-13217*a*d-239*b*d-12800*c*d+13358*d^2-14594*a*e-11038*b*e-5292*c*e+13035*d*e+8675*a*f+1374*b*f+6183*c*f-7476*a*g-1297*b*g+10788*c*g+14261*a*h-15529*b*h-7537*a*i+11471*a*j+e*j+3186*a*k-545*a*l-15384*b*m,
    12793*a^2+801*a*b+7529*b^2-8930*a*c+13009*b*c+11261*c^2+3960*a*d+14005*b*d-9937*c*d+1891*d^2+3506*a*e-6726*b*e+9708*c*e-4633*d*e-9449*a*f-12215*b*f+10856*c*f-13940*a*g-15381*b*g+14765*c*g-14164*a*h-12517*b*h+13825*a*i-693*a*j-1624*a*k+e*k+12256*a*l+12046*c*m,
    9656*a^2-1706*a*b-8785*b^2+13823*a*c+4896*b*c+15908*c^2-660*a*d+5924*b*d+13718*c*d+177*d^2-12920*a*e-6014*b*e-8060*c*e-4036*d*e+7720*a*f-1251*b*f-6271*c*f-5272*a*g+9674*b*g-10087*c*g+10895*a*h-8000*b*h+9780*a*i+9913*a*j-9834*a*k-14127*a*l+e*l-15062*d*m,
    -11212*a^2+9628*a*b+13912*b^2-6030*a*c+10649*b*c+14053*c^2-5540*a*d+12578*b*d-9131*c*d-10738*d^2+744*a*e+10566*b*e-9913*c*e-10132*d*e-11421*a*f+13119*b*f-1670*c*f+f^2+7494*a*g+1669*b*g-1800*c*g+1574*a*h+8350*b*h-14969*a*i-8651*a*j-3130*a*k,
    6551*a^2+3996*a*b-7421*b^2-10832*a*c+9558*b*c+5928*c^2+1973*a*d+7145*b*d+15092*c*d-8247*d^2+6657*a*e-7919*b*e-11698*c*e-1707*d*e+6539*a*f-11288*b*f+11969*c*f+12374*a*g-9392*b*g+14279*c*g+f*g+14835*a*h+8455*b*h+5364*a*i+3250*a*j+1416*a*k-2872*a*l,
    198*a^2-11252*a*b-13827*b^2-8060*a*c-3072*b*c+8182*c^2-986*a*d+11184*b*d+6140*c*d+5034*d^2-3758*a*e+11250*b*e-9924*c*e+8025*d*e+10910*a*f+15003*b*f-7695*c*f-8736*a*g+10043*b*g-11151*c*g-3728*a*h+680*b*h+f*h+15627*a*i-6342*a*j+12534*a*k-12654*a*l-13626*a*m,
    6921*a^2+3559*a*b+7923*b^2-13195*a*c-946*b*c+6368*c^2-10144*a*d+13520*b*d-4155*c*d+9393*d^2+8178*a*e-2075*b*e-7125*c*e-6269*d*e+293*a*f+9539*b*f-11515*c*f+13270*a*g+4003*b*g-9478*c*g-12659*a*h-12600*b*h+6057*a*i+f*i+3241*a*j-7676*a*k-8772*a*l+14343*b*m,
    15890*a^2-485*a*b-2865*b^2+1709*a*c+7649*b*c+15154*c^2-573*a*d-502*b*d-11303*c*d+15329*d^2-2894*a*e-6780*b*e-8602*c*e+13135*d*e+372*a*f-193*b*f-3666*c*f-2402*a*g-156*b*g+5977*c*g-8697*a*h-8578*b*h+6753*a*i+15626*a*j+f*j-13840*a*k-13790*a*l-12269*c*m,
    851*a^2+2783*a*b-12697*b^2-1573*a*c-8704*b*c-12070*c^2-8949*a*d+8923*b*d-150*c*d+4499*d^2+3026*a*e-993*b*e+12929*c*e+995*d*e+15051*a*f-2544*b*f+13361*c*f-13123*a*g-11199*b*g+15040*c*g-9740*a*h-353*b*h-6750*a*i-6498*a*j-7315*a*k+f*k-6443*a*l-361*d*m,
    -11688*a^2-14277*a*b-1760*b^2-13199*a*c-14863*b*c-14650*c^2-3696*a*d+411*b*d-9400*c*d-11299*d^2-6002*a*e+4122*b*e+9673*c*e+11827*d*e-7995*a*f-6769*b*f-11706*c*f-2893*a*g-11316*b*g+12566*c*g+3210*a*h-14794*b*h-12645*a*i+597*a*j+9663*a*k+13807*a*l+f*l+7066*e*m,
    -9455*a^2+2631*a*b+13477*b^2+5787*a*c-12305*b*c-2967*c^2-15707*a*d+14313*b*d-11170*c*d-13709*d^2-1703*a*e+5902*b*e-1200*c*e-12710*d*e+2536*a*f+320*b*f+8342*c*f-14259*a*g+6536*b*g+15587*c*g+g^2-8837*a*h+10440*b*h-14407*a*i-11506*a*j+12597*a*k+11074*a*l-11611*a*m,
    626*a^2-9899*a*b+9713*b^2-11122*a*c+4359*b*c+9049*c^2-15883*a*d-1944*b*d+3669*c*d-5480*d^2-8423*a*e-7972*b*e+463*c*e-10943*d*e-14914*a*f+5727*b*f+7290*c*f+5885*a*g-15110*b*g-933*c*g-10491*a*h-4827*b*h+g*h-6833*a*i-1285*a*j-13595*a*k+407*a*l-3544*b*m,
    9360*a^2-6278*a*b+7926*b^2-8195*a*c+7445*b*c+15249*c^2-13368*a*d+3714*b*d-7875*c*d-13228*d^2+14411*a*e+2301*b*e+3480*c*e-9095*d*e-10992*a*f+2846*b*f+2822*c*f+13695*a*g+4222*b*g-9180*c*g-5789*a*h-1995*b*h-447*a*i+g*i+2679*a*j-3745*a*k-9258*a*l+9294*c*m,
    -2225*a^2-1812*a*b-14044*b^2+14528*a*c+1917*b*c-10690*c^2-2011*a*d-7726*b*d+11227*c*d-13911*d^2+13920*a*e+13668*b*e+10186*c*e-8897*d*e+9292*a*f-9641*b*f-5624*c*f+10696*a*g-12338*b*g-8497*c*g-8705*a*h+4415*b*h+13829*a*i+1226*a*j+g*j+3010*a*k-13580*a*l+11786*d*m,
    3754*a^2+2579*a*b+3588*b^2-11716*a*c-7109*b*c+12102*c^2-10661*a*d-8458*b*d-742*c*d+7977*d^2+482*a*e+10557*b*e-12558*c*e-2328*d*e-14266*a*f+6585*b*f+6611*c*f-9981*a*g-13103*b*g+6230*c*g+12012*a*h+3156*b*h+9167*a*i+7885*a*j-8831*a*k+g*k+9213*a*l-2146*e*m,
    14542*a^2-1100*a*b-5928*b^2-12245*a*c+14925*b*c-11371*c^2+7123*a*d-6204*b*d+5271*c*d+7558*d^2+6150*a*e-961*b*e-12149*c*e-3274*d*e-13479*a*f+3411*b*f+4794*c*f+24*a*g+15405*b*g+198*c*g+10229*a*h-5308*b*h+2061*a*i+10285*a*j+6233*a*k+13099*a*l+g*l+14716*f*m,
    -9242*a^2-6040*a*b+5572*b^2-8489*a*c-6472*b*c-1831*c^2+6524*a*d+1218*b*d+11408*c*d-4814*d^2-1156*a*e-3944*b*e+2573*c*e+2146*d*e-5579*a*f-4870*b*f+8968*c*f-597*a*g-14007*b*g+11476*c*g+13748*a*h+14418*b*h+h^2-9819*a*i-8893*a*j+10580*a*k-1020*a*l-1859*c*m,
    6544*a^2+15826*a*b-2938*b^2-858*a*c-9249*b*c-12127*c^2+1545*a*d-5534*b*d-13406*c*d-14884*d^2-2247*a*e-6466*b*e+9438*c*e+12567*d*e-8457*a*f+12374*b*f-4008*c*f+15970*a*g+11204*b*g+13116*c*g+6041*a*h+15382*b*h+6691*a*i+h*i-5456*a*j-1492*a*k-5921*a*l-1653*d*m,
    1603*a^2-2015*a*b-8036*b^2+7638*a*c+15550*b*c+15017*c^2-6594*a*d+829*b*d+4141*c*d+965*d^2+7983*a*e-10976*b*e+6238*c*e-15798*d*e-1573*a*f-9243*b*f+14452*c*f+7676*a*g-3066*b*g-10038*c*g+2810*a*h-2535*b*h-2862*a*i+14409*a*j+h*j-5324*a*k-7934*a*l-13318*e*m,
    -6956*a^2+12754*a*b-915*b^2+3316*a*c-14307*b*c+10575*c^2+61*a*d+9690*b*d+4070*c*d-5681*d^2-3669*a*e+7082*b*e-13438*c*e-9240*d*e+6275*a*f-4006*b*f+9388*c*f-13324*a*g-11868*b*g-10133*c*g-999*a*h+4785*b*h+12180*a*i+6708*a*j+7110*a*k+h*k-1929*a*l+2470*f*m,
    14313*a^2+8199*a*b+6100*b^2+12429*a*c-10724*b*c+251*c^2-5111*a*d-8611*b*d+622*c*d+9242*d^2-3027*a*e-13773*b*e+13159*c*e-15521*d*e-14513*a*f-5873*b*f-12674*c*f+2076*a*g-6694*b*g+1119*c*g+1573*a*h+14895*b*h+12581*a*i+9932*a*j-8413*a*k-15705*a*l+h*l+14147*g*m,
    -1107*a^2+3660*a*b-7901*b^2+9428*a*c+9172*b*c+12685*c^2-14172*a*d-4144*b*d+1527*c*d-640*d^2+12164*a*e+7*b*e-11748*c*e+12121*d*e-7502*a*f+13833*b*f-9993*c*f-15567*a*g+418*b*g+15396*c*g+7320*a*h-253*b*h-3294*a*i+i^2+14025*a*j+10916*a*k-5637*a*l-8712*e*m,
    15463*a^2-6262*a*b-15021*b^2-10426*a*c-12253*b*c-5377*c^2+11043*a*d+6357*b*d-15529*c*d-9891*d^2-3237*a*e+10713*b*e-6797*c*e+11499*d*e-2765*a*f+6214*b*f-1978*c*f+15507*a*g+2271*b*g-12340*c*g-3007*a*h-10280*b*h-1222*a*i+3149*a*j+i*j-12352*a*k+6565*a*l+7952*f*m,
    1002*a^2+13450*a*b+11343*b^2-15339*a*c-12781*b*c-5239*c^2+1239*a*d-13442*b*d+2802*c*d-1727*d^2+14729*a*e+6924*b*e+10802*c*e-4779*d*e+146*a*f+3162*b*f+9245*c*f-1715*a*g-6868*b*g+14906*c*g+10495*a*h+13330*b*h+1336*a*i-14461*a*j-3496*a*k+i*k-10530*a*l+7468*g*m,
    -15935*a^2-11864*a*b-1259*b^2-4370*a*c-2973*b*c+1772*c^2+13588*a*d+12473*b*d+2679*c*d-4649*d^2+11136*a*e+7256*b*e-1548*c*e+10924*d*e+5213*a*f+9250*b*f+4878*c*f-7399*a*g-8553*b*g+3583*c*g+7596*a*h-1315*b*h+3579*a*i-6750*a*j-14956*a*k-7667*a*l+i*l-4690*h*m,
    -1811*a^2+15975*a*b-4110*b^2-3981*a*c+7814*b*c-12586*c^2-4749*a*d-7601*b*d-8044*c*d+5095*d^2+15935*a*e-10576*b*e+9007*c*e+2028*d*e-6326*a*f-4564*b*f+11677*c*f+14808*a*g-3025*b*g-5845*c*g-9061*a*h-12633*b*h-11237*a*i+10934*a*j+j^2+11236*a*k+15022*a*l-13383*g*m,
    1184*a^2+7487*a*b-13592*b^2-10579*a*c+9964*b*c+7923*c^2+14628*a*d+15013*b*d-12043*c*d-11477*d^2+5809*a*e-1350*b*e-6067*c*e+13866*d*e+13291*a*f-10043*b*f-5606*c*f-12060*a*g-13678*b*g-11881*c*g+15799*a*h-4490*b*h-8688*a*i-14110*a*j-3801*a*k+j*k-3286*a*l+5783*h*m,
    11267*a^2+13675*a*b-13080*b^2-12531*a*c+13987*b*c+1737*c^2-11479*a*d-2516*b*d-1660*c*d+312*d^2-3469*a*e-12816*b*e-2480*c*e+6015*d*e-5111*a*f-6501*b*f-1999*c*f+7271*a*g-5224*b*g+6262*c*g+7719*a*h+9535*b*h-11497*a*i-7448*a*j+14448*a*k+4479*a*l+j*l+1103*i*m,
    -14831*a^2-5823*a*b-9965*b^2+6035*a*c+9925*b*c+10836*c^2+1752*a*d-11642*b*d-12553*c*d-12784*d^2-4036*a*e-9739*b*e-6673*c*e-15010*d*e-13269*a*f-8723*b*f+4789*c*f+13839*a*g+10215*b*g-4787*c*g-4639*a*h+11790*b*h+2254*a*i+8175*a*j+5248*a*k+k^2+1255*a*l-15646*i*m,
    -11825*a^2-15973*a*b+2289*b^2-1642*a*c-4923*b*c+15388*c^2-6258*a*d-15849*b*d+10907*c*d-3031*d^2+5435*a*e-7665*b*e+2019*c*e+3427*d*e+13267*a*f-3395*b*f-377*c*f+11879*a*g+10047*b*g-6818*c*g-2218*a*h-7278*b*h-6336*a*i-3884*a*j+2020*a*k-14148*a*l+k*l+2938*j*m,
    2093*a^2-15275*a*b+9226*b^2-4769*a*c-1058*b*c+3637*c^2-13774*a*d-2763*b*d+12439*c*d-11817*d^2+4591*a*e-13451*b*e+1360*c*e-13759*d*e+6900*a*f-13287*b*f-13349*c*f-11774*a*g+13544*b*g+15636*c*g+5786*a*h+14624*b*h+8282*a*i+9729*a*j+2123*a*k-5406*a*l+l^2+10511*k*m
    )

