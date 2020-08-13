setRandomSeed 0
path = prepend("../../", path)
debug needsPackage "NumericalDecomposition"
-- k points close to V(F)

newtonBag = method(Options=>{Verbose=>false})
newtonBag (Thing, ZZ, Function,Function) := o ->  (F,k,generatePoint,isGoodPoint) -> (
    ret := new MutableList;
    i := 0;
    numIterations := 5;
    --Statistics are printed when Verbose=>true
    totalSeedsUsed := 0;
    while i < k do (
	totalSeedsUsed=totalSeedsUsed+1;
	--p := point random(CC^1,CC^(numVariables F));
	p := generatePoint();
	for iteration to numIterations do (
	    p = newton(F,p);	   
	    );
    	if isGoodPoint p then (
	    ret#i = p;
	    i = i + 1;
	    )	    
    	);
    if o.Verbose then print("Number of seeds used: "|totalSeedsUsed);	
    if o.Verbose then print("Convergence percentage: "|toString sub(#ret/totalSeedsUsed,RR));	
    return toList ret
    )

newtonBag(Thing, ZZ, Function) := o -> (F,k,generatePoint) ->(
    isGoodPoint :=    p -> p.ErrorBoundEstimate < getDefault(CorrectorTolerance);
    newtonBag(F,k,generatePoint,isGoodPoint,o)
    )

newtonBag(Thing, ZZ) := o ->  (F,k) ->(
    generatePoint := () -> point random(CC^1,CC^(numVariables F));
    newtonBag(F,k,generatePoint,o)
    )



end

restart
load "bag-of-newton.m2"
declareVariable \ {x,y,z,w}

-- one component in A^4
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
--- all true?
elapsedTime tally for x in xs list membershipTest(x, W0)
elapsedTime tally for x in xs list membershipTest(x, W0, Backtrack=>true)
tally for i from 1 to 100 list membershipTest(point random(CC^1,CC^4), W0)

-- two components
f = x*z-y^2
g = y*w-z^2
F = gateSystem({x,y,z,w},transpose gateMatrix{{f,g}})
k = 100
xs=newtonBag(F,k)
x0 = first xs
W0 = witnessCurve(F, entries transpose vars F, x0)
populate W0
degree W0
membershipTest(x0, W0)
tally for x in xs list membershipTest(x, W0)
x1 = first select(1, xs, x -> not membershipTest(x, W0))
W1 = witnessCurve(F, entries transpose vars F, x1)
populate W1
degree W1
tally for x in xs list membershipTest(x, W1)


-- fiber product
-*
this example has 52 IrrComps of dim 2, with
1 of degree 5
15 of degree 20
25 of degree 60
11 of degree 120
*-
restart
load "bag-of-newton.m2"
V=toList vars(X_1..X_5,A,B)
F = gateSystem(V,transpose gateMatrix{for i from 1 to 5 list X_i^5+A*X_i+B})
k = 100
elapsedTime xs=newtonBag(F,k);

witnessCurves = new MutableList 
nWitnessCurves = 0
elapsedTime witnessLabels = apply(xs, x -> (
    xWitnesses := toList positions(witnessCurves, W -> membershipTest(x, W, Backtrack=>true));
    if (length xWitnesses == 0) then (
        << " making witness curve" << endl;
        elapsedTime witnessCurves#nWitnessCurves = witnessCurve(F, {flatten entries transpose vars F}, x);
        << "populating witness curve" << endl;
        elapsedTime populate(witnessCurves#nWitnessCurves, Verbose=>false, NumberOfNodes=>3);
        xWitnesses = {nWitnessCurves};
        nWitnessCurves = nWitnessCurves + 1;
        << "now we've got " << nWitnessCurves << " components" << endl;
    );
    xWitnesses
    )
);
tally witnessLabels
toList apply(witnessCurves, w -> degree w)
first positions(witnessLabels, i -> i == {14})


-*
Given: a bag of points
Want: 
(1) create a WitnessCurve from each point, if possible (if a random slice is transverse)
(2) for each point, say on which component(s) it is 

Methods to create: 

membershipTest(Point,WitnessCurve)

*-


--Example of Newton Bag for isolated points. 

--Case v=4
v=4;
R = CC[X_0, X_1, X_2, X_3, nLittle]
sys = {2*X_0*X_1+2*X_2*X_3,
    X_1^2+2*X_0*X_2+X_3^2,
    2*X_1*X_2+2*X_0*X_3,
    .0625*X_0^2+.75*nLittle-1,
    X_1*X_3-4*nLittle}


--Using the standard method we have a poor convergence percentage to find five solutions.
B = newtonBag(polySystem (sys),5,Verbose=>true);
--Number of seeds used: 3735
--Convergence percentage: .00133869

--Let's generate our points in a structured manner to capture structure of the problem. 
--Specifically, we fix nLittle, and fix the first half of X_1...,X_v rounded down (A)
--The remaining variables can be sovled for in terms of the fixed ones (x0 and B).
generatePoint = ()->(
    chooseN:=random CC;
    x0:=sqrt(    (v)^2-4*chooseN*(v-1));
    if v%2==1 then(
	A := apply(sub((v-1)/2,ZZ),o->random CC	);
	B := apply(A,i->4*chooseN/i) ;
	);
    if v%2==0 then(
	A := apply(sub((v)/2,ZZ)-1,o->random CC	);
	B := {sqrt(4*chooseN)}|apply(A,i->4*chooseN/i) ;
	);
    point(matrix{{x0}|A|B|{chooseN}})
    )

--We see a much better convergence percentage by generating the points in the structured manner.

B = newtonBag(polySystem sys,5,generatePoint,Verbose=>true);netList B
--Number of seeds used: 13
--Convergence percentage: .384615

--Looking at the solutions we find some with non-positive integer last coordinates. 
-* 
ii53 : B = newtonBag(polySystem sys,5,generatePoint,Verbose=>true);netList B

       +----------------------------------------------------------------------------------------------------------------------+
oo54 = |{2, -2*ii, 2, 2*ii, 1}                                                                                                |
       +----------------------------------------------------------------------------------------------------------------------+
       |{4, -5.65967e-14+6.90155e-15*ii, -5.65967e-14+6.90126e-15*ii, -5.65969e-14+6.90128e-15*ii, 2.82395e-14-3.58282e-15*ii}|
       +----------------------------------------------------------------------------------------------------------------------+
       |{2, -2*ii, 2, 2*ii, 1}                                                                                                |
       +----------------------------------------------------------------------------------------------------------------------+
       |{4, 1.35052e-14+1.82267e-14*ii, -1.45392e-14-1.81778e-14*ii, 1.45071e-14+1.71347e-14*ii, 7.03564e-15+7.59842e-15*ii}  |
       +----------------------------------------------------------------------------------------------------------------------+
       |{4, -2.28908e-16-1.06978e-16*ii, -2.28911e-16-1.06978e-16*ii, -2.2891e-16-1.06976e-16*ii, 1.133e-16+5.30256e-17*ii}   |
       +----------------------------------------------------------------------------------------------------------------------+
*-

--If we only care about solutions with a positive integer in the last coordinate,
--then we can set isGoodPoint as follows

isGoodPoint=p->(
    nLittle:=(coordinates p)#-1;
    --Check tolerance
    if p.ErrorBoundEstimate < getDefault(CorrectorTolerance) 
    --check last coordinate is greater than or equal to 1
    then if abs(imaginaryPart(nLittle))<getDefault(CorrectorTolerance) and realPart( nLittle)>1-getDefault(CorrectorTolerance) 
    then (
    	--check last coordinate is a positive integer
	while nLittle>1+getDefault(CorrectorTolerance) do nLittle =nLittle-1; 
	if min{abs(nLittle), abs(nLittle-1)}<getDefault(CorrectorTolerance) 
	then return true
	);
    return false
    )

B = newtonBag(polySystem (sys),5,generatePoint,isGoodPoint,Verbose=>true);
netList B
--Number of seeds used: 21
--Convergence percentage: .238095

-*

       +----------------------+
oo63 = |{2, -2*ii, 2, 2*ii, 1}|
       +----------------------+
       |{2, 2*ii, 2, -2*ii, 1}|
       +----------------------+
       |{2, -2*ii, 2, 2*ii, 1}|
       +----------------------+
       |{2, 2*ii, 2, -2*ii, 1}|
       +----------------------+
       |{2, 2*ii, 2, -2*ii, 1}|
       +----------------------+
*-
