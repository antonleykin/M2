--Organization:

--0) imports, global variables, & uncategorized service functions
   -- this contains some syntactic sugar that may be controversial
--1) creating and manipulating GateSystems and related objects
--2) extracting information from gate systems (squaring, seeding, root counts)

-- Conventions that need to be better enforced:
-- inputs should be row matrices
-- outputs should be column matrices
-- parameters come before variables
-- "GS" is better variable name for a GateSystem than G or another single letter

--0) imports, global variables, & uncategorized service functions
needsPackage "MonodromySolver"
needsPackage "Polyhedra"
needsPackage "NumericalAlgebraicGeometry"

-- debugging
-*
0 => no debugging
1 => some debugging
2 => much debugging
*-
DBG = 0 

-- global symbols that let us create slicing parameters on the fly
sliceCounter = 0
SS = symbol SS

-- fold objects with a "||" method in a list (ie matrices along rows)
foldVert = method()
foldVert List := L -> if (#L == 0) then L else fold(L,(a,b)->a||b)

-- fold objects with a "|" method in a list (ie matrices along cols)
foldHor = method()
foldHor List := L -> if (#L == 0) then L else fold(L,(a,b)->a|b)

-- helpers for squareDown
-- orthonormal basis for col(L) using SVD
ONB = L -> (
    (S,U,Vt) := SVD L;
    r := # select(S,s->not areEqual(s,0));
    U_{0..r-1}
    )
-- component of col(L) that is perpendicular to M
perp = method(Options=>{UseSVD=>true})
perp (Matrix, Matrix) := o -> (M, L) -> if areEqual(norm L, 0) then M else (
    Lortho := if o.UseSVD then ONB L else first complexQR L; -- QR seems buggy
    Lperp := M-Lortho*conjugate transpose Lortho * M;
    if o.UseSVD then ONB Lperp else first complexQR Lperp
    )

--1) creating and manipulating GateSystems and related objects

isSquare = GS -> (numVariables GS == numFunctions GS)

-- gates for small determinants
dimensions = method()
dimensions Thing := M -> (numrows M, numcols M)
det2 = M -> M_(0,0)*M_(1,1)-M_(1,0)*M_(0,1)
det3 = M -> M_(0,0)*det2(M^{1,2}_{1,2})-M_(0,1)*det2(M^{1,2}_{0,2})+M_(0,2)*det2(M^{1,2}_{0,1})
laplaceDet = M -> (
    (m, n) := dimensions M;
    assert(m==n);
    if (m==2) then det2 M else if (m==3) then det3 M else sum(m,i-> (
            inds := delete(i, toList(0..m-1));
            (-1)^i*M_(0,i) * laplaceDet(M_inds^inds)
            )
        )
    )
-- convenience functions for minors
minors (GateMatrix, ZZ, Sequence, Boolean) := o ->  (M, k, S, laplace) -> (
    (Sm, Sn) := (first S, last S);
    (m,n) := dimensions M;
    assert(k<=min(m,n));
    assert(all(Sm,s->#s==k));
    assert(all(Sn,s->#s==k));
    gateMatrix{flatten apply(Sm,sm->apply(Sn, sn -> 
	    if (laplace) then laplaceDet submatrix(M,sm,sn)
	    else det submatrix(M,sm,sn)
	    ))}
    )

allMinors = method(Options=>{Laplace=>false})
allMinors (GateMatrix, ZZ) := o -> (M, k) -> (
    (m, n ) := dimensions M;
    s := (subsets(0..m-1,k),subsets(0..n-1,k));
    minors(M, k, s, o.Laplace)
    )

maxMinors = method(Options=>{})
maxMinors GateMatrix := o -> M -> (
    r := min dimensions M;
    useLaplace := (r < 4);
    allMinors(M,r, Laplace=> useLaplace)
    )

-- "join" of two GateSystems (take all functions from both)
GateSystem || GateSystem := (P, Q) -> (
    allVars := unique( (flatten entries vars P) | (flatten entries vars Q) );
    allParams := unique( (flatten entries parameters P) | (flatten entries parameters Q) );
    gateSystem(
	gateMatrix{allParams},
	gateMatrix{allVars},
    	(gateMatrix P)||(gateMatrix Q)
	)
    )

-- sum of two GateSystems
GateSystem + GateSystem := (P, Q) -> (
    if (numFunctions P =!= numFunctions Q) then error "can only add GateSystems of the same shape";
    H := P || Q;
    gateSystem(parameters H, vars H, gateMatrix P + gateMatrix Q)
    )

-- take some functions from the GateMatrix \ GateSystem
GateMatrix ^ BasicList := (M, inds) -> M^(toList inds)
GateSystem ^ BasicList := (P, inds) -> gateSystem(parameters P, vars P, (gateMatrix P)^inds)


-- append some slices to a given GateSystem
sliceSystem = method(Options => {Affine => 0, Homog => 0})
sliceSystem (GateMatrix, GateSystem) := o -> (X, F) -> (
    if (o.Affine <= 0 and o.Homog <= 0) then error("you did not do the slice");
    F || foldVert for i from 0 to o.Affine + o.Homog - 1 list (
	m := if i < o.Affine then numcols X + 1 else numcols X;
	X' := if i < o.Affine then transpose(X | gateMatrix{{1_CC}}) else transpose X;
    	sliceParams := gateMatrix{for i from 0 to m-1 list (
	    	ret := SS_sliceCounter;
	    	sliceCounter = sliceCounter + 1;
	    	ret)};
    	gateSystem(sliceParams, vars F, sliceParams * X')
	)
    )
sliceSystem GateSystem := o -> F -> sliceSystem(vars F, F, o)

-- sub new values for variables
sub (GateMatrix, GateSystem) := (X, F) -> gateSystem(parameters F, X, sub(gateMatrix F, vars F, X))    

-- trace
trace GateMatrix := SumGate => M -> (
    if numrows M =!= numcols M then error "expected a square matrix";
    sum(numrows M, i -> M_(i, i))
    )
-- product
GateMatrix * Gate := (M, a) -> a*M

gateSystem (BasicList, BasicList, GateMatrix) := (P, X, F) -> gateSystem(gateMatrix{toList P}, gateMatrix{toList X}, columnize F)
gateSystem (BasicList, BasicList, BasicList) := (P, X, Fs) -> foldVert apply(toList Fs, f -> gateSystem(P, X, gateMatrix f))
gateSystem (BasicList, BasicList) := (X, Fs) -> foldVert apply(toList Fs, f -> gateSystem({}, X, gateMatrix f))
gateSystem (Thing, Thing, Gate) := (X, P, g) -> gateSystem(X, P, gateMatrix{{g}})
gateSystem (List, Thing) := (X, F) -> gateSystem({}, X, F)

vars IndexedVariable := x -> declareVariable x
vars Symbol := x -> declareVariable x
vars InputGate := x -> x
InputGate .. InputGate := (A, B) -> value \ (A.Name .. B.Name)

gateMatrix VisibleList := L -> if instance(L, GateMatrix) then L else gateMatrix{toList L}
gateMatrix (List, ZZ, ZZ) := (L, m, n) -> (
    if # L =!= m*n then error "wrong number of list entries";
    gateMatrix for i from 0 to m-1 list take(L,{n*i,n*(i+1)-1})
    )
gateMatrix (List, ZZ) := (L, m) -> gateMatrix({L}, 1, m)


--2) extracting information from gate systems (squaring, seeding, root counts, initializing Jacobians)

-- make all parameters new variables
flatten GateSystem := GS -> gateSystem(vars GS | parameters GS, gateMatrix GS)

-- 1xm matrix of random parameter values for a GateSystem
randomParameters = method()
randomParameters (Thing, GateSystem) := (FF, GS) -> random(FF^1, FF^(numParameters GS))
randomParameters GateSystem := GS -> randomParameters(CC_53, GS)

-- reshape a matrix to be a matrix of outputs
columnize = method()
columnize Matrix := M -> transpose matrix{flatten entries M}
columnize GateMatrix := M -> transpose gateMatrix{flatten entries M}

gateMatrix Gate := g -> gateMatrix{{g}}

GateSystem | GateSystem := (A, B) -> A||B

-- switch variables for parameters
switch GateSystem := GS -> gateSystem(vars GS, parameters GS, gateMatrix GS)

evaluate (GateMatrix, GateMatrix, Point) := (GM, varMat, xVals) -> (
    Peval := gateSystem(varMat, transpose gateMatrix{flatten entries GM});
    result := evaluate(Peval, xVals);
    matrix(result, numrows GM, numcols GM)
    )
evaluate (GateMatrix, VisibleList, Point) := (GM, xVars, xVals) -> evaluate(GM, gateMatrix{xVars}, xVals)
evaluate (GateMatrix, VisibleList, Matrix) := (GM, xVars, xVals) -> (
    I := if instance(xVars, GateMatrix) then xVars else gateMatrix{toList xVars};
    evaluate(GM, I, point xVals)
    )
evaluate (Matrix, Thing, Thing) := (M, thing1, thing2) -> evaluate(gateMatrix M, thing1, thing2)


newton (GateSystem, Matrix, Matrix) := (GS, p0, x0) -> (
    D0 := evaluateJacobian(GS, p0, x0);
    F0 := transpose evaluate(GS, p0, x0);
    dx0 := transpose solve(D0, F0, ClosestFit => true); -- does this use SVD?
    x0 - dx0
    )

newtonFixedPointSystem = method()
newtonFixedPointSystem GateSystem := GS -> (
    fs := flatten entries gateMatrix GS;
    F := sum(fs/(f->f^2));
    gateSystem(parameters GS, vars GS, transpose diff(vars GS, gateMatrix{{F}}))
    )

refine (GateSystem, Matrix, Matrix) := o -> (GS, p0, x0) -> (
    absRes := norm evaluate(GS, point p0, point x0);
    maxIters := if instance(o.Iterations, Nothing) then 5 else o.Iterations;
    residualTolerance := if instance(o.ResidualTolerance, Nothing) then 1e-6 else o.ResidualTolerance;
    cur := x0;
    nIters := 0;
    while (absRes > residualTolerance and nIters < maxIters) do (
	cur = newton(GS, p0, cur);
	absRes = norm evaluate(GS, point p0, point cur);
	nIters = nIters + 1;
	if (DBG > 1) then (
	    << nIters << "th iteration" << endl;
	    << cur << endl;
	    );
    	);
    cur
    )

conditionNumber = (nFP, p0, cur) -> (
    S := first SVD evaluateJacobian(nFP, p0, cur);
    (max S)/(min S)
    )

-- very primitive atm, but works ok for square systems
-- todo: implement greedy alternation heuristic for (well-constrained) overdetermined systems
seedNewton = method(Options=>
    -- these are all options for "refine" for now
    {Iterations=>10, ResidualTolerance=>1e-6}
    )
seedNewton (Matrix, Matrix, GateSystem) := o -> (p0, x0, GS) -> (sub(p0,CC), refine(GS, sub(p0,CC), sub(x0,CC), o))
seedNewton (Matrix, GateSystem) := o -> (p0, GS) -> seedNewton(p0, random(CC^1,CC^(numVariables GS)), GS, o)
seedNewton GateSystem := o -> GS -> seedNewton(randomParameters GS, GS, o)

///TEST
restart
-- clean this up
needs "nag_supplemental.m2"
setRandomSeed 0
vars {x,y,p,q,r,s}
G=gateSystem({p,q,r,s}, {x,y}, {p*(2*x*y^2+x*y)+q*(x^2-1), r*(x^2*y-x*y^2)+s*(1-x*y)})
x0=random(CC^1,CC^2)
p0 = sub(matrix{{0,1,0,1}},CC)
DBG=2
evaluateJacobian(G,p0,matrix{{1_CC,1}})
evaluate(G,p0,matrix{{1_CC,1}})
newton(G,p0,newton(G,p0,newton(G,p0,newton(G,p0,newton(G,p0,newton(G,p0,newton(G,p0,newton(G,p0,x0))))))))
(p1, x1)=seedNewton(p0,x0,G)
evaluate(G,p0,x0)
assert(areEqual(
	norm evaluate(G,p1,x1)
	,0))
///

rowSelector = method(Options=>{BlockSize=>1,UseSVD=>true,Verbose=>false})
rowSelector (Point, Point, GateSystem) := o -> (y0, c0, GS) -> (
    (n, m, N) := (numVariables GS, numParameters GS, numFunctions GS);
    blockSize := o.BlockSize;
    numBlocks = ceiling(N/blockSize);
    numIters=0;
    L := matrix{for i from 1 to n list 0_CC}; -- initial "basis" for row space
    r := 0;
    goodRows := {};
    diffIndices := {};
    while (r < n and numIters < numBlocks) do (
    	diffIndices = for j from numIters*blockSize to min((numIters+1)*blockSize, N)-1 list j;
	if o.Verbose then << "processing rows " << first diffIndices << " thru " << last diffIndices << endl;
    	newRows := evaluateJacobian(GS^diffIndices, y0, c0);
    	for j from 0 to numrows newRows - 1 do (
	    tmp := transpose perp(transpose newRows^{j}, transpose L);
	    if not areEqual(0, norm tmp) then (
		if o.Verbose then << "added row " << blockSize*numIters+j << endl;
	    	if areEqual(norm L^{0}, 0) then L = tmp else L = L || tmp;
	    	goodRows = append(goodRows, blockSize*numIters+j);
		);
    	    );
    	r = numericalRank L;
    	numIters = numIters+1;
	);
    if o.Verbose then << "the rows selected are " << goodRows << endl;
    goodRows
    )

-- candidate
squareDown = method(Options=>{BlockSize=>1, Verbose=>false})
squareDown (Point, Point, GateSystem) := o -> (y0, c0, F) -> F^(rowSelector(y0, c0, F, BlockSize => o.BlockSize, Verbose=>o.Verbose))

-- IN: gateSystem
-- OUT: polynomials with random rational numbers substituted for parameters 
expand = method()
expand (Thing, GateSystem) := (FF, GS) -> (
    n := numVariables GS;
    aaa := symbol aaa;
    flatGS := flatten GS;
    nm := numVariables GS;
    R := frac(QQ[aaa_1..aaa_n]);
    RCC := FF[gens R];
    ringEls := flatten entries evaluate(flatGS, vars R | randomParameters(QQ,GS));
    if any(ringEls, r -> first degree denominator r > 0) then error " didnt want to expand rational functions";
    polys := numerator \ ringEls;
    polys/(p -> sub(p,RCC))
    )
expand GateSystem := GS -> expand(CC_53, GS)

-- store precomputed Newton Polytopes as their vertices
newtonPolytopes = GS -> (
    if (not GS#?"Newton Polytopes") then (
    	expGS := expand GS;
	GS#"Newton Polytopes" = vertices \ newtonPolytope \ expGS;
--	assert all(GS#"Newton Polytopes", np->areEqual(np, sort np)); -- this should always be true
	);
    GS#"Newton Polytopes"
    )

computeMixedVolume GateSystem := GS -> (
    assert isSquare GS;
    computeMixedVolume expand GS
    )
sparseSubsystemIndices = GS -> (
    (N, n) := (numFunctions GS, numVariables GS);
    assert(N >= n);
    nPs := newtonPolytopes GS;
    nPClasses := partition(i->nPs#i,0..N-1);
    classSizes := (values nPClasses)/length;
    nClasses := #classSizes;
    nPPartitions := toList \ select(partitions(numVariables GS,max classSizes), p -> #p == nClasses);
    apply(nPPartitions, p -> foldHor apply(values nPClasses,p, (eqInds,pEl) -> take(eqInds, pEl)))
    )
computeMixedVolume (GateSystem, List) := (GS, inds) -> inds/(ind -> computeMixedVolume GS^ind)
computeAllMixedVolumes = method()
computeAllMixedVolumes GateSystem := GS -> (
    inds := sparseSubsystemIndices GS;
    vols := computeMixedVolume(GS, inds);
    (inds, vols)
    )

-- monomial matrix: yVars^F
exp (GateMatrix, Matrix) := (yVars, F) -> gateMatrix{for i from 0 to numrows F-1 list compress product for j from 0 to numcols F - 1 list (
	e := sub(F_(i, j),ZZ);
	if e >0 then yVars_(0,j)^e else if e<0 then inputGate(1)/yVars_(0,j)^(-e) else inputGate(1)
	)
    }
exp (Matrix, Matrix) := (y0, F) -> matrix{for i from 0 to numrows F-1 list product for j from 0 to numcols F - 1 list (
	e := sub(F_(i, j),ZZ);
	y0_(0,j)^e
	)
    }

-- supports = vertices of Newton polytopes
sparseFamily = method()
sparseFamily GateSystem := GS -> (
    assert isSquare GS;
    ww := symbol ww;
    n := numVariables GS;
    X := vars GS;
    vertsGS := newtonPolytopes GS;
    coefficientParameters := apply(n, i -> 
	apply(numcols vertsGS#i, j -> 
	    declareVariable ww_(flatten entries (vertsGS#i)_{j}, i)
	    )
	);
    exps := unique flatten(entries \ transpose \ vertsGS);
    mons := exp(X, sub( matrix exps ,ZZ));
    polys := for cPs in coefficientParameters list sum(cPs, cp -> (
	    cpExp := (cp.Name)#1#0;
	    monIndex := position(exps, e -> e == cpExp);
	    cp * mons_(0, monIndex)
	    )
	);
    gateSystem(gateMatrix{foldHor coefficientParameters}, vars GS, transpose gateMatrix{polys})
    )

degrees GateSystem := GS -> first \ degree \ (expand GS)

///TEST
vars {x,y,p}
G=gateSystem({x,y}, {x^2*y+2*x*y^2+x*y-1, x^2*y-x*y^2-x*y+2})
computeMixedVolume G
computeMixedVolume gateSystem({x,y}, {x^2*y+2*x*y^2+x*y-1, x^2*y-x*y^2-x*y+p})
///

end--
restart
needs "nag_supplemental.m2"
needsPackage "MonodromySolver"
(M, N) = (3, 4)
assert(D==0)
vars(T_(M+1)..T_(M+N))
for j from M+1 to M+N do S_j = inputGate(1)/T_j
-- unknowns: G1,G2, (links) z1,z2 (pivots ) + conjugates (8) and T1..TN for pose requirements
vars(T_(M+1)..T_(M+N))
-- parameters: D1..D(M+N) and conjugates (precision 
vars(d_(1,1)..d_(M+N,2))
for j from 1 to M+N do D_j = matrix{{d_(j,1)},{d_(j,2)}}
