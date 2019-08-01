needs "step1.m2"
debug needsPackage "MonodromySolver"	
declareVariable \ {x,y,z,w}
f = 1 + 2*x + 3*y^2+4*z^3+5*w^4
h = 1 + 2*x + 3*y+ 5*z + 7*w + 11*z*y + 13 * x* z + 17*x*w + 19*y*z+23*y*w +29*z*w+31*x*y*z+37*x*y*w+41*x*z*w+43*y*z*w+47*x*y*z*w
F = gateSystem(matrix{{x,y,z,w}},transpose gateMatrix{{f,h}})
G = new VariableGroup from {{0},{1},{2},{3}}
pt = point{{4,-2/3,-1.44419, .765304}}
P = multiaffineDimension(F,G,pt)

-- output of step 1), input to step 2)
assert (
    P == convexHull transpose matrix{
    {1,1,0,0},
    {1,0,1,0},
    {1,0,0,1},
    {0,1,1,0},
    {0,1,0,1},
    {0,0,1,1}
    }
    )

needs "step2.m2"
SCseq = getSequenceSC P 

-- output of step 2
assert(SCseq == ({2,3},2,{1,2},{0,1},0))

-- output of step 3 
-*
-- outline for step 3 function
-- output: 
  -- Params: List of parameters
  -- L: GateMatrix representing the slice
makeSlices = (I,G,SCseq) := -> (
*-

needs "step3.m2"
(params, L) = makeSliceSystem(F,G,SCseq)

pp = flatten entries params
assert(
    L === transpose gateMatrix{
    	{pp#0*z+pp#1*w+pp#2,
     	    pp#3*x+pp#4*y+pp#5*z+pp#6*w+pp#7
     	    }
 	}
    )

needs "steps4-6.m2"
masterGS = gateSystem(
    params,
    vars F,
    (gateMatrix F)||L
    )
x0 =point sub(matrix pt,CC)
createSeedPair(masterGS,x0)

errorDepth=2
G=populateWitnessCollection(masterGS,x0,Verbose=>true
--    ,TraceTolerance=>1e-10
--    ,AugmentNodeCount=>2
--    ,AugmentEdgeCount=>3
--    ,NumberOfRepeats=>40
--    ,MaxNumTraceTests=>5
    )
(toList G.Vertices)/(v->#points v.PartialSols)


V = first G.Vertices
monodromyResult = (masterGS, V.BasePoint, V.PartialSols)
assert (length last monodromyResult == 12)
end--

restart 
needs "octahedron-do-steps.m2"
    
-- example
x1 = point{{-2_CC,1,0,0}}
errorDepth = 100
NAGtrace 3
elapsedTime membershipTest(x1,F,monodromyResult)
elapsedTime membershipTest(x1,F,monodromyResult,Backtrack=>true)
x2 = point{{0,0,0,0}}
membershipTest(x2,F,monodromyResult)

