--errorDepth = 0
needsPackage "NumericalDecomposition"
R = CC[x,y,z]
numericalAffineSpace R
NV = numericalAffineSpace R
comps = components NV
C = first comps
sph = (x^2+y^2+z^2-1); 
--NAGtrace 4
uGeneration {sph*(y-x^2), sph*(z-x^3)}

end
restart
load "sphere-and-twisted-cubic.m2"
