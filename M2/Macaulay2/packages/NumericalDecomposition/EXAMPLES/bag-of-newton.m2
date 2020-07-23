debug needsPackage "NumericalDecomposition"
-- k points close to V(F)
newtonBag = (F,k) -> (
    ret := new MutableList;
    i := 0;
    numIterations := 5;
    while i < k do (
	p := point random(CC^1,CC^(numVariables F));
	for iteration to numIterations do (
	    p = newton(F,p);	   
	    );
       	if p.ErrorBoundEstimate < getDefault(CorrectorTolerance) then (
	    ret#i = p;
	    i = i + 1;
	    )
	);	
    return toList ret
    )

end

load "bag-of-newton.m2"
declareVariable \ {x,y,z,w}
f = 1 + 2*x + 3*y^2+4*z^3+5*w^4
h = 1 + 2*x + 3*y+ 5*z + 7*w + 11*z*y + 13 * x* z + 17*x*w + 19*y*z+23*y*w +29*z*w+31*x*y*z+37*x*y*w+41*x*z*w+43*y*z*w+47*x*y*z*w
F = gateSystem(matrix{{x,y,z,w}},transpose gateMatrix{{f,h}})
k = 10
newtonBag(F,k)
