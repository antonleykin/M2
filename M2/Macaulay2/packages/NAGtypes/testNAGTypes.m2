--In this file we test out the NAG types for multiprojective varieties. 

restart
debug needsPackage "NAGtypes"
debug needsPackage "NumericalAlgebraicGeometry"
--load "/Users/jo/Documents/GitHub/M2/M2/Macaulay2/packages/NAGtypes/WSet-abstract.m2"
moveSlicingVariety

R=QQ[x1,x2,x3]

* `` lays foundation
* `WSet-ambient.m2` describes affine and projective `Ambient` and `SlicingVariety`
* `WSet-multiaffine.m2` a stab at mutli- affine/projective witness set/collection

