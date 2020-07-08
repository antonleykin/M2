path=prepend("/Users/jo/Documents/GoodGit/AntonM2/M2/Macaulay2/packages",path)
needsPackage"NumericalDecomposition"
     R = CC[x1,x2,x3,y];
     F = gateSystem polySystem {x1-x2-y};
     pt = point{{2,1,1,3_CC}};
     G={{0,1,2},{3}}--Groups the variables
     P = multiaffineDimension(F,G,pt)
     ambDim P == #G--number of groups
     scs = getSequenceSC P   
     V= gens R
     debug NumericalDecomposition
     describeSCS(V, G,scs)
-*
Prints:
slice in {x1, x2, x3}
slice in {x1, x2, x3}
coarsen:  {{x1, x2, x3}, {y}}
slice in {x1, x2, x3, y}
*-

     R = CC[x1,x2,x3,y];
     V= gens R
     F = gateSystem polySystem {x1-x2-y};
     pt = point{{2,1,1,3_CC}};
     G={{0,1,2},{3}}--Groups the variables
     scs =  getSequenceSC (F,G,pt)
     describeSCS  (V,G,scs)


     V= gens R
     F = gateSystem polySystem {x1-2};
     G= {{0},{1},{2},{3}}
     scs =  getSequenceSC (F,G,pt)
     describeSCS  (V,G,scs)


     R = CC[x_0..x_9];

     V= gens R
     F = gateSystem polySystem {x_0-1};
     G= apply(10,i->{i})
     pt = point{{1}|apply(#gens R-1,i->random CC)};
     #gens R
     #flatten G
     #coordinates pt
     scs =  getSequenceSC (F,G,pt)
     describeSCS  (V,G,scs)


end
restart
load "/Users/jo/Documents/GoodGit/AntonM2/M2/Macaulay2/packages/NumericalDecomposition/TESTground/describeSCSeq.m2"
debug NumericalDecomposition
check NumericalDecomposition
installPackage"NumericalDecomposition"
