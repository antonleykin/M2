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

-- prepare the homotopy graph
populateWitnessCollection = method(
    Options=>{NumberOfEdges => 2,
       NumberOfRepeats => 15, BatchSize => infinity, Potential => null,
       EdgesSaturated => false, Randomizer => null, SelectEdgeAndDirection =>selectRandomEdgeAndDirection
       , Verbose => false, FilterCondition => null, TargetSolutionCount =>
       null, "new tracking routine" => true,
       NumberOfNodes => 2, MaxNumTraceTests => 3, TraceTolerance => 1e-4})
populateWitnessCollection (GateSystem, Point) := o -> (MasterGS, x0) -> (
    n := # coordinates x0;
    p1 := first createSeedPair(MasterGS,x0);
    G := homotopyGraph MasterGS;
    addNode(G, p1, pointArray{x0});
    -- assume complete graph 
    completeGraphInit(G,p1,first G.Vertices,o.NumberOfNodes,o.NumberOfEdges);
    coreNodes := toList G.Vertices;
    -- two "special nodes"
    addNode(G, point{drop(coordinates p1,-1)|{random CC}}, pointArray{});
    addNode(G, point{drop(coordinates p1,-1)|{random CC}}, pointArray{});
    specialNodes := toList take(G.Vertices,-2);
    for v in coreNodes do for u in specialNodes do for e from 1 to o.NumberOfEdges do addEdge(G,u,v);
    -- annoying things that need to be set before calling the core solver
    USEtrackHomotopy=true;
    setTrackTime(G,0);
    -- calling the core solver
    mutableo := new MutableHashTable from o;
    mutableo#StoppingCriterion = ((n,L) -> (n >= o.NumberOfRepeats));
    remove(mutableo, MaxNumTraceTests);
    remove(mutableo, TraceTolerance);
    numTraceTests := 0;
    traceVal := infinity;
    V := first G.Vertices;
    pencilNodes := {V} | specialNodes;
    pencilCounts := {0,0,0};
    newSols := {{},{},{}};
    partialTraces := map(CC^3,CC^n,0);
    done := false;
    while not done do (
	-- we should really augment G at some point in this loop...
    	coreMonodromySolve(G,
    	    first G.Vertices,
	    new OptionTable from mutableo
    	    );
    	newSols = for i from 0 to 2 list(
	    matrix \ drop(points (pencilNodes#i).PartialSols, pencilCounts#i)
	    );
	pencilCounts = pencilCounts + length \ newSols;
	partialTraces = partialTraces + fold(newSols/sum,(a,b)->a||b);
    	S := first SVD(
	    (partialTraces^{1}-partialTraces^{0})||
	    (partialTraces^{2}-partialTraces^{0})
	    );
	traceVal = min S;
	numTraceTests = numTraceTests + 1;
	if (traceVal < o.TraceTolerance) then (
	    << "trace test succeeds :)" << endl;
	    done = true;
	    )
	else (
	    << "trace test fails :(" << endl;
	    if (o.MaxNumTraceTests < numTraceTests) then done = true;
	    );
	);
    G
    )

end--

restart 
needs "octahedron-do-steps.m2"
MasterGS = gateSystem(
    params,
    vars F,
    (gateMatrix F)||L
    )
x0 =point sub(matrix pt,CC)
createSeedPair(MasterGS,x0)

G=populateWitnessCollection(MasterGS,x0,Verbose=>true)

(toList V.Graph.Vertices)/(v->#points v.PartialSols)
# toList V.Graph.Edges

-- STEP 5: trace test
p1 = V.BasePoint
 -- trace test passes!

    
-- example
MonodromyResult = (MasterGS, p1, p1Sols)
x1 = point{{-2_CC,1,0,0}}
member(x1,F,MonodromyResult)
x2 = point{{0,0,0,0}}
member(x2,F,MonodromyResult)