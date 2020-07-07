path=prepend("/Users/jo/Documents/GoodGit/AntonM2/M2/Macaulay2/packages",path)
needsPackage"NumericalDecomposition"
     R = CC[x1,x2,x3,y];
     F = gateSystem polySystem {x1-x2-y};
     pt = point{{2,1,1,3_CC}};
     G={{0,1,2},{3}}
     P = multiaffineDimension(F,G,pt)
     G={{0,1,2},{3}}--Groups the variables
     ambDim P == #G--number of groups
     scs = getSequenceSC P   
     describeSCS(V, G,scs)
-*
Prints:
slice in {x1, x2, x3}
slice in {x1, x2, x3}
coarsen:  {{x1, x2, x3}, {y}}
slice in {x1, x2, x3, y}
*-

end
load "/Users/jo/Documents/GoodGit/AntonM2/M2/Macaulay2/packages/NumericalDecomposition/TESTground/describeSCSeq.m2"
debug NumericalDecomposition

