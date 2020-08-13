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


