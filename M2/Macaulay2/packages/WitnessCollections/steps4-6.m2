-- prepare the homotopy graph
populateWitnessCollection = method(
    Options=>{AugmentEdgeCount=>0, AugmentNodeCount=> 1, NumberOfEdges => 1,
       NumberOfRepeats => 10, BatchSize => infinity, Potential => null,
       EdgesSaturated => false, Randomizer => null, SelectEdgeAndDirection =>selectBestEdgeAndDirection,
       Potential => potentialLowerBound,
       Verbose => false, FilterCondition => null, TargetSolutionCount =>
       null, "new tracking routine" => true,
       NumberOfNodes => 2, MaxNumTraceTests => 3, TraceTolerance => 1e-4})
populateWitnessCollection (GateSystem, Point) := o -> (masterGS, x0) -> (
    n := # coordinates x0;
    p1 := first createSeedPair(masterGS,x0);
    G := homotopyGraph(masterGS, Potential=>o.Potential);
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
    -- we need to mangle the options before calling the core solver
    mutableo := new MutableHashTable from o;
    mutableo#StoppingCriterion = ((n,L) -> (n >= o.NumberOfRepeats));
    remove(mutableo, MaxNumTraceTests);
    remove(mutableo, TraceTolerance);
    remove(mutableo, AugmentNodeCount);
    remove(mutableo, AugmentEdgeCount);
    -- trace test initialize
    numTraceTests := 0;
    traceVal := infinity;
    V := first G.Vertices;
    pencilNodes := {V} | specialNodes;
    pencilCounts := {0,0,0};
    newSols := {{},{},{}};
    partialTraces := map(CC^3,CC^n,0);
    done := false;
    nMonodromyPaths := 0;
    while not done do (
	if (numTraceTests > 0) then completeGraphAugment(G,p1,V,o.AugmentEdgeCount,o.NumberOfEdges,o.AugmentNodeCount);
    	nMonodromyPaths = nMonodromyPaths + last coreMonodromySolve(G,
    	    first G.Vertices,
	    new OptionTable from mutableo
    	    ); -- this modifies G
    	newSols = for i from 0 to 2 list(
	    matrix \ drop(points (pencilNodes#i).PartialSols, pencilCounts#i)
	    );
	pencilCounts = pencilCounts + length \ newSols;
	partialTraces = partialTraces + fold(
	    newSols/(sL -> if #sL>0 then sum sL else map(CC^1,CC^n,0))
	    ,(a,b)->a||b);
    	S := first SVD(
	    (partialTraces^{1}-partialTraces^{0})||
	    (partialTraces^{2}-partialTraces^{0})
	    );
	traceVal = min S;
	NotZeroMatrix := (max S > o.TraceTolerance);
	numTraceTests = numTraceTests + 1;
	if (traceVal < o.TraceTolerance and NotZeroMatrix) then (
	    << "trace test succeeds :)" << endl;
	    done = true;
	    )
	else (
	    << "trace test fails :(" << endl;
	    if (o.MaxNumTraceTests < numTraceTests) then done = true;
	    );
	<< "we have tracked " << nMonodromyPaths << " paths so far (includes trace test nodes)" << endl;
	<< " and we ran " << numTraceTests << " trace test(s)" << endl;
	);
    G
    )

membershipTest = method(Options=>{Backtrack=>false})
membershipTest (Point, GateSystem, Sequence) := o -> (x1, F, MonodromyResult) -> (
    (masterGS, p1, p1Sols) := MonodromyResult;
    assert(instance(masterGS,GateSystem) and instance(p1,Point));
    sliceTargParams := first createSeedPair(masterGS, x1);
    specializationParameters := if o.Backtrack then transpose((matrix sliceTargParams)|(matrix p1)) else transpose((matrix p1)|(matrix sliceTargParams));
    MembershipHomotopy := specialize(parametricSegmentHomotopy masterGS, specializationParameters);
    startSols := if o.Backtrack then pointArray{x1} else p1Sols;
    -- in the line below: it would be nice if trackHomotopy could take a PointArray of solutions
    targetSols := pointArray trackHomotopy(MembershipHomotopy,points startSols);
    if o.Backtrack then member(first points targetSols, p1Sols) else member(x1, targetSols)
    )
end--


