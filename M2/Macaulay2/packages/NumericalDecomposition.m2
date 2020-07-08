newPackage(
	"NumericalDecomposition",
    	Version => "1.14",
    	Date => "Aug 2019",
    	Authors => {
	     {Name => "Timothy Duff", Email => "tduff3@gatech.edu"},
	     {Name => "Anton Leykin", Email => "leykin@math.gatech.edu"},
	     {Name => "Jose Rodriguez", Email => "jrodriguez43@wisc.edu"}
	     },
    	Headline => "numerical decomposition of varieties",
	PackageImports => {},
	PackageExports => {"NAGtypes","NumericalAlgebraicGeometry","Polyhedra","MonodromySolver"},
	DebuggingMode => true,
	AuxiliaryFiles => true
    	)
debug MonodromySolver

-- types
export lines ///VariableGroup
WitnessCurve
///
-- functions and options
export lines ///multiaffineDimension
getSequenceSC
makeSliceSystem
Backtrack
MaxNumTraceTests
TraceTolerance
witnessCurve
populate///
-- TYPES --
VariableGroup = new Type of List
WitnessCurve = new Type of HashTable

-- SERVICE FUNCTIONS

-- these should be handled in NumericalAG
createJacobian = method()
createJacobian GateSystem := GS -> (
    F := gateMatrix GS;
    I := vars GS;
    J := diff(I,F);
    (m,n) := (numrows J,numcols J);
    JGS := gateSystem(I, transpose gateMatrix{flatten entries J});
    (m,n,JGS)
    )
evaluateJacobian (Point, ZZ, ZZ, GateSystem) := (pt, m,n,JGS) -> matrix(evaluate(JGS,matrix pt),m,n)

-- TO DO: find a better way to generate unique symbols
protect sP
sliceParam = symbol sP
sliceParamCounter = 0
newSliceParam = () -> (
    ret := inputGate sliceParam_sliceParamCounter;
    sliceParamCounter = sliceParamCounter + 1;
    ret
    )

-- override buggy code in MonodromySolver
createSeedPair (System, Point) := o -> (P, x0) -> (
    G := if instance(P, GateSystem) then P else gateSystem P.PolyMap;
    n := numVariables G;
    m := numParameters G;
    N := numFunctions G;
    I := id_(CC^m);
    A := random(CC^0,CC^N);
    scan(m, i -> A = A || evaluate(G, point I_{i}, x0));
    b := evaluate(G, point matrix 0_(CC^m), x0);
    K := numericalKernel(transpose A, 1e-5) ;
    offset := solve(transpose A,transpose b,ClosestFit=>true);
    p0 := point(K* random(CC^(numcols K), CC^1) - offset);
    (p0, x0)
    )

-- STEP 1 --

--This function returns a polyhedron.

multiaffineDimension = method(TypicalValue=>Thing)
multiaffineDimension(GateSystem, List, Point) := (F,G,pt)->multiaffineDimension(F,new VariableGroup from G,pt)
multiaffineDimension(GateSystem, VariableGroup, Point) := (F,G,pt)->( --(polynomial system, k variable groups, a general witness point)
    (m,n,JF) := createJacobian F; --JF is a m by n matrix
    Jpt := evaluateJacobian(pt,m,n,JF);
    thePartialJacs:= apply(G, vg -> Jpt_vg); -- Jpt^vg takes rows
    --all nonempty subsets of [0,1,..,#G-1]
		if instance(ring matrix pt,InexactFieldFamily) or instance(ring matrix pt,InexactField)
		then useRank :=  numericalRank
		else useRank =  rank;
		intrinsicCodimension :=  useRank  Jpt;
    subsetsIV := drop(subsets(#G),1);
    M := {};
    v := {};
    scan(subsetsIV, Icomplement->(
    	    M = append(M,apply(#G,i->if member(i,Icomplement) then 0 else 1 ));
    	    varsInIcom := flatten G_Icomplement;
	    v = append(v,
		n-intrinsicCodimension +
		-#varsInIcom +useRank Jpt_varsInIcom
	    )));
    --Inequalities.
    M = matrix M;
    v = transpose matrix {v};
--    print (M,v);
    --Equalities
    N := matrix {apply(#G,i->1)};
    w := matrix{{n -intrinsicCodimension}};
    --  print (N,w);
    --M*matrix{{e_1},...,{e_k}} \leq v
    --N*matrix{{e_1},...,{e_k}} = w
    P:=polyhedronFromHData(M,v,N,w);
    return P
    )



-- STEP 2 --
getSequenceSC = method(TypicalValue=>Sequence)
--TODO: have an option Strategy => ZZ
getSequenceSC (Polyhedron) := (P)->(
    SCS := ();
    if isEmpty(P) then error" P is empty. " ;
    -- k is the number of factors and number of variable groups
    k := ambDim P;
    -- bfe is a maximal integer vector such that for each i with bfe_i=!=0 there exists a vector bfe+ei in P.
    bfe := new MutableList from {};
    newP:=P;
    scan(k,i->(
	    --ei is the ith basis vector
    	    ei := for j to k-1 list if i==j then 1 else 0;
	    --project newP to the ith coordinate to get a lattice polytopy in R,
	    -- which is just a list of integers. The largest integer is maxEi
    	    maxEi := max\\flatten\\flatten\flatten\entries \latticePoints affineImage(matrix{ei},newP);
    	    --if maxEi is greater than one then Bertini's theorem applies and we can slice.
	    if maxEi > 1 then (
		bfe#i = maxEi-1;
	    	if bfe#i>0 and not isEmpty(newP)
	    	then (
		    Q :=polyhedronFromHData(-matrix{ei},-matrix {{bfe#i}}); -- newM \leq e_I
    	    	    ---newP on the lhs consists of integer vectors in newP on the rhs that are greater than or equal to bfe (coordinatewise).
  	    	    newP = intersection(newP,Q);
	    	    ))
	    else bfe#i = 0));
    idMatrix:= diagonalMatrix(apply(k,i->1));
    --newP becomes a matroid polytope by shifting by -bfe
    newP=affineImage(idMatrix,newP,-matrix transpose {toList bfe});
    scan(#bfe,i->scan(bfe#i,j -> SCS=append(SCS, i)));
    --print("SCS"=>SCS);
    dimLowerBound := 1;
    numFactors := k;
    scan(k,i->(
    	    ai := for j to k-1 list if k-1-j<=i then 1 else 0;
--    	    print ai;
    	    maxAi := latticePoints affineImage(matrix{ai},newP);
	    maxAi = maxAi/entries/flatten/flatten//flatten//max;
    	    --maxAi is the dimension of the projection to the i+1 last coordinates.
--	    print maxAi;
	    if i>0 then (--then we coarsen the last two factors
		SCS = append(SCS,{k-1-i,k-i});
   	    	numFactors = numFactors - 1);
    	    if maxAi > dimLowerBound then (
--    	    	print("test"=>(k-1-i,numFactors-1)); --these numbers should be the same
	    	SCS = append(SCS, numFactors-1)
		);
	    dimLowerBound = max(maxAi,dimLowerBound)));
    assert(numFactors ==1);
    SCS=append(SCS,0);
    SCS
    )

-- Step 2 without needing Step 1. 
getSequenceSC(GateSystem, List, Point) := (F,G,pt)->( --(polynomial system, k variable groups, a general witness point)
    (m,n,JF) := createJacobian F; --JF is a m by n matrix
    Jpt := evaluateJacobian(pt,m,n,JF);
    thePartialJacs:= apply(G, vg -> Jpt_vg); -- Jpt^vg takes rows
    --all nonempty subsets of [0,1,..,#G-1]
    if instance(ring matrix pt,InexactFieldFamily) or instance(ring matrix pt,InexactField)
    then useRank :=  numericalRank
    else useRank =  rank;
    intrinsicCodimension :=  useRank  Jpt;
    -- We project onto the last coordinate, last two coordinates, etc, and determine the dimension of the image. 
    imageDimensions :=  apply(#G-1,i->(
	Icomplement := for j to #G-1-i list j;
	--print Icomplement;
	varsInIcom := flatten G_Icomplement;
	(n-intrinsicCodimension) - (#varsInIcom -useRank Jpt_varsInIcom)	
	));
    imageDimensions=append(imageDimensions,(n-intrinsicCodimension));
    scs :={};
    D:=0;
    k:=#G-1;
    numSlices:=0;
    scan(imageDimensions,i->(
	    scan(i-numSlices-1, j->(
		    scs = append(scs, k);
		    numSlices = numSlices +1
		    )
		);
	    if k>0 then (
		scs = append(scs,{k-1,k});
		k=k-1)	    
	    )
	);
    toSequence scs
    )


describeSCS = method()
--V is a list of variables
--G is a 2D-List of integers giving the variables by grouping
--SCSeq is a slice coarsening sequence
--Convention choice: for {a,b,c} in SCseq we coarsen by putting the variables in groups a b c all in group a and delete the old groups b c.
describeSCS (List,List,Sequence) := (V,G,SCseq) -> (    
    G = new MutableList from G;
    scan(SCseq,i->(
	    if instance(i,ZZ) 
	    then print ("slice in "|toString apply(G#i,j->V#j))
	    else if instance(i,List)
	    then (
		CG:=apply(i, j -> G#j );
		print ("coarsen:  "|toString apply(CG,X->apply(X,x ->V#x)));
		G#(i#0) = flatten CG;
		scan(drop(i,1), j -> G#j = null);
		G = new MutableList from delete(null,toList G);		
		)
	    )
	)
    )






-- STEP 3 --
makeSliceSystem = method()
makeSliceSystem (GateSystem, VariableGroup, Sequence) := (F, G, SCseq) -> (
    ParamList := new MutableList from {};
    Leqs := new MutableList from {};
    Blocks := new MutableHashTable from for i from 0 to #G-1 list i=>G#i; -- keeps track of identified variable groups as coarsening progresses
    localSliceParamCounter := 0; -- we want to keep track of how many sliceParamList have been generated for each witness collection
    sliceCounter := 0;
    I := vars F;
    for s in SCseq do (
    	if instance(s,ZZ) then (-- make a slice
	    blockSize := #(Blocks#s);
	    -- add in slice to L...
	    Leqs#sliceCounter = (sum for i in Blocks#s list newSliceParam()*I_(0,i))+newSliceParam();
	    sliceCounter = sliceCounter + 1;
	    -- ... and add in parameters to ParamList
	    for j from 0 to blockSize do (
	    	ParamList#(localSliceParamCounter+j) =  sliceParam_(sliceParamCounter-blockSize-1+j);
	    	);
	    localSliceParamCounter = localSliceParamCounter + blockSize + 1;
	    )
    	else if (instance(s,List) and #s == 2) then (
	    -- coarsen varParts
	    (i,j) := (first s, last s);
    	    Blocks#i = Blocks#i | Blocks#j;
	    for k from j to #keys Blocks-2 do Blocks#j = Blocks#(j+1);
	    remove(Blocks,#Blocks-1);
	    )
    	else error "invalid encoding of coarsening / slicing moves";
    	);
    gateSystem(
	gateMatrix{toList ParamList} | parameters F, -- parameters
    	vars F,
    	(gateMatrix F) || transpose gateMatrix{toList Leqs}
    	)
    )

-- STEP 4-6
-- prepare the homotopy graph
populate = method(
    Options=>{AugmentEdgeCount=>0, AugmentNodeCount=> 1, NumberOfEdges => 1,
       NumberOfRepeats => 10, BatchSize => infinity, Potential => null,
       EdgesSaturated => false, Randomizer => (p->p), SelectEdgeAndDirection =>selectBestEdgeAndDirection,
       Potential => potentialLowerBound,
       Verbose => false, FilterCondition => null, TargetSolutionCount =>
       null, "new tracking routine" => true,
       NumberOfNodes => 2, MaxNumTraceTests => 3, TraceTolerance => 1e-4})
populate WitnessCurve := o -> W -> (
    masterGS := W#"system with slices";
    p1 := W.cache#"SpecializationParameters";
    x0 := first W.cache#"WitnessPoints";
    n := # coordinates x0;
    G := homotopyGraph(W#"system with slices", Potential=>o.Potential, Randomizer=>o.Randomizer);
    addNode(G, p1, W.cache#"WitnessPoints");
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
--    G
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

-*
 todo:
 -) isComplete (wrt. some optional tolerance)
 -) store partial trace test values in W
 -) repopulate
 -) enhancements (rational map, etc.)
 -) refinement (Anton's end)
*-
witnessCurve = method(Options=>{})--eventual options might include SC strategy, for instance
witnessCurve (GateSystem, List, Point) := o -> (F, Blocks, pt) -> (
    I := vars F;
    V := toList flatten entries I;
--    assert(set flatten Blocks == set V);
    assert(# flatten Blocks == # V);
    G := new VariableGroup from for b in Blocks list(b/(v->position(V,x->x===v)));
    V/(v->position(Blocks,B->member(v,B)));
    m := numcols parameters F;
    P := multiaffineDimension(F,G,pt);
    SCseq := getSequenceSC P;
    masterGS := makeSliceSystem(F,G,SCseq); -- this should always take the point, have the default of returning a square subsystem, and allow for parameters in F
    new WitnessCurve from {
	-- unnecessary keys (for keeping record)
	"original system" => F, 
	"SC sequence" => SCseq,
	-- necessary keys
	"system with slices" => masterGS,
	"NumVars" => # V,
	"NumParams" => m,
	"NumSliceParams" => (numcols parameters masterGS)-m,
	"NumSlices" => numFunctions masterGS - numFunctions F, -- assumes square
	cache => new CacheTable from {
	    "SpecializationParameters" => first createSeedPair(masterGS, pt),
	    "WitnessPoints" => pointArray {pt}
	    }
	}
    )

points WitnessCurve := W -> points W.cache#"WitnessPoints"

TEST /// -- one factor of dim=2
needsPackage "NumericalDecomposition"
declareVariable \ {x,y}
F = gateSystem(
    gateMatrix{{x,y}},
    gateMatrix{{y^2-x^3-x+1}}
    )
pt = point{{2.00000001,-3+0.0000001*ii}}
assert(norm evaluate(F,pt)<0.0001)
setRandomSeed 0
W = witnessCurve(F, {{x,y}}, pt)
populate W
points W
///

TEST /// -- two factors of dim=1
needsPackage "NumericalDecomposition"
declareVariable \ {x,y}
F = gateSystem(
    gateMatrix{{x,y}},
    gateMatrix{{y^2-x^3-x+1}}
    )
pt = point{{2.00000001,-3+0.0000001*ii}}
assert(norm evaluate(F,pt)<0.0001)
setRandomSeed 0
W = witnessCurve(F, {{x},{y}}, pt)
populate W
points W
///

TEST /// -- one factor of dim=4 with parameters
restart
needsPackage "NumericalDecomposition"
needs "./NumericalDecomposition/nag_supplemental.m2"
vars(x0,x1,y0,y1,ax,ay) -- x0 & y0 are homogenizing variables
fs={x0^3*y1^2+y0^2*(-x1^3-x1*x0^2+x0^3),ax*x1+x0-1,ay*y1+y0-1}
F = gateSystem({ax,ay},{x0,x1,y0,y1},fs)
pt'a = point{{0_CC,0}}
pt = point{{1,2.00000001,1,-3+0.0000001*ii}}
setRandomSeed 0
W = witnessCurve(F, {{x0,y0,x1,y1}}, pt)
points W
errorDepth = 0
populate(W,Verbose=>true)
GS = W#"system with slices"
p0 = W.cache#"SpecializationParameters"
x0 = first W.cache#"WitnessPoints"
assert areEqual(0, norm evaluate(GS,p0,x0))
monodromySolve(GS,p0,{x0},Verbose=>true, NumberOfNodes=>4, NumberOfEdges=>1, Randomizer=>(p->p))
///

-*
-- symbolic 
restart
R=QQ[x0,x1,y1,y0]
ax = random QQ
ay = random QQ
I=ideal(x0^3*y1^2+y0^2*(-x1^3-x1*x0^2+x0^3),ax*x1+x0-1,ay*y1+y0-1, random(1,R)-random QQ)
degree I
*-

-- IN: a list of points
-- OUT: a hash table of 
--      WitnessCurve => list of points (that belong to the component witnessed by the key)
decompose(List, System, List) := (pts,system,varParts) -> (
    F := gateSystem system;
    curves := new MutableHashTable;
    for pt in pts do (
	C := select(1,keys curves,c->"is0n"(pt,c)); -- TODO: implement isOn (replaces membershipTest)
	if C=!=null then curves#C = append(curves#C,pt) else (
	    C = witnessCurve(F,varParts,pt);
	    populate C;
	    curves#C = {pt};
	    )
	); 
    new HashTable from curves
    )

beginDocumentation()
needs "./NumericalDecomposition/Documentation/doc12.m2"
needs "./NumericalDecomposition/Documentation/doc3456.m2"
end

restart
uninstallPackage "NumericalDecomposition"
installPackage "NumericalDecomposition"
installPackage("NumericalDecomposition", RemakeAllDocumentation=>true)
check "NumericalDecomposition"
peek NumericalDecomposition
