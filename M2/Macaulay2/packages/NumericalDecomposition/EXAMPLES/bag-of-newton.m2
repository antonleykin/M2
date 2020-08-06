setRandomSeed 0
path = prepend("../../", path)
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

restart
load "bag-of-newton.m2"

declareVariable \ {x,y,z,w}
f = 1 + 2*x + 3*y^2 + 4*z^3 + 5*w^4
g = 1 + 2*x + 3*y+ 5*z + 7*w 
h = g + 11*z*y + 13 * x* z + 17*x*w + 19*y*z+23*y*w +29*z*w+31*x*y*z+37*x*y*w+41*x*z*w+43*y*z*w+47*x*y*z*w
F = gateSystem({x,y,z,w},transpose gateMatrix{{f,h}})
k = 100
xs=newtonBag(F,k)
x0 = first xs
W0 = witnessCurve(F, entries transpose vars F, x0)
populate W0
degree W0
membershipTest(x0, W0)
tally for x in xs list membershipTest(x, W0)
tally for i from 1 to 100 list membershipTest(point random(CC^1,CC^4), W0)


-*
Given: a bag of points
Want: 
(1) create a WitnessCurve from each point, if possible (if a random slice is transverse)
(2) for each point, say on which component(s) it is 

Methods to create: 

membershipTest(Point,WitnessCurve)

*-

uninstallPackage "MonodromySolver"
restart
installPackage "MonodromySolver"
check "MonodromySolver"
