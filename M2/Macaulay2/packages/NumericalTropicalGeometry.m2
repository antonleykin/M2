newPackage(
	"NumericalTropicalGeometry",
	Version => "0.9",
	Date=> "September 2014",
	Authors => {
	    {Name=>"Anton Leykin", Email=>"leykin@math.gatech.edu"},
	    {Name=>"Josephine Yu", Email=>"jyu@math.gatech.edu"}						
	    },	  
	Headline => "Calculates the Tropical Variety of a 1-dimensional Ideal by Numerical Methods",
	PackageExports => {"NumericalAlgebraicGeometry"},
	PackageImports => {"PHCpack"},
	DebuggingMode=> true
	)

export{sampleAmoeba, Factor, Multiplier, Extra, 
    StartSystemAndSolutions, RayIntegerBound,
    tropicalIntersectionMultiplicityTimesRayMultiplicity, Epsilon,
    guessRays,validateRay,
    combineDuplicates, reformatValidatedRays, appendNewData,
    PrincipalVector,GoodPoints,    
    degreeOfCurveInTorus,tropicalDegree}

INFINITY := 1e9; -- threshold for infinite solutions

sampleAmoeba = method(Options=>{Factor=>5_RR, Multiplier=>1000_RR, Extra=>1000_CC, 
	RayIntegerBound=>9, Software=>M2engine})

ShortApproxCollinearVectorAlt = method()
ShortApproxCollinearVectorAlt (List, RR) := (p, multiplier) -> (
     -- print p;
     n := #p;
     if n < 2 then return p; -- check for special case
     pivot := max(p/abs);
     pivotloc :=position(p, s->abs(s)==pivot); 
     M := mutableMatrix(random(RR^(n-1), RR^n)); --generate the upper matrix
     for i from 0 to pivotloc-1 do(
	for j from 0 to n-1 do(
		if i==j then M_(i,j)=-p#(pivotloc) else (if j==pivotloc then M_(i,j)=p#(i) else M_(i,j)=0);
	);
     ); 
     for i from pivotloc to n-2 do(
        for j from 0 to n-1 do(
		if i+1==j then M_(i,j)=-p#(pivotloc) else (if j==pivotloc then M_(i,j)=p#(i+1) else M_(i,j)=0);
	);
     );
--     print M;
     M = matrix(M)*(multiplier/pivot);
     Mzz := mutableMatrix(random(ZZ^(n-1),ZZ^n));
     for i from 0 to n-2 do(
	for j from 0 to n-1 do(
		Mzz_(i,j)=round(M_(i,j));
	);
     ); 
     M = matrix(Mzz) || map ZZ^n;
     -- << "M = " << M << endl;
     L := LLL M;
     -- print L;
     --print toString (Mp,L);
     v := flatten entries L_{0}^(toList(n-1..2*n-2));
     i := position(v, x->x!=0);
     if v#i*p#i<0 then v=-v;
     -- print v;
     return v; 
     ) 

--------------------------------------------------------------------------------------
guessRays = method()
guessRays (Ideal, Number, List, ZZ) := (I,C,numAttemptsList,coord'bound) -> (
    R := ring I;
    n := numgens R;
    RCC := CC[gens R];
    ICC := sub(I,RCC);
    LoL := flatten apply(#numAttemptsList, d->
    	apply(numAttemptsList#d, count->(
    	    	ind := take(random toList (0..(n-1)),d+1);
    	    	w := apply(n,i->if member(i,ind) then (-1)^(random 2) else 0);  
    	    	--try 
	    	sampleAmoeba(ICC,w,C,Software=>PHCPACK) 
	    	--else {}
    	    	))
    	);
    select(unique join toSequence LoL, w->max (w/abs) < coord'bound)
    )

-- if binom=(a,b) then "slice" with (a + Ctb), pick C nonzero randomly?
tropicalIntersectionMultiplicityTimesRayMultiplicity = method(Options=>{Software=>M2engine, Epsilon=>0})
tropicalIntersectionMultiplicityTimesRayMultiplicity(List,Ideal,Sequence,Symbol) := o -> (w,I,binom,xxx) -> (
    R := ring I;
    assert(numgens R == numgens I + 1);
    (a,b) := binom;
    ab := matrix{first exponents a - first exponents b};
    c := 1/(ab*transpose matrix{w})_(0,0);
    assert( c==1 );
    tt := symbol tt;
    Rt := (coefficientRing R) (monoid [tt,gens R]);
    RtoRt := map(Rt,R,drop(Rt_*,1));
    Rt'frac := frac Rt;
    t := first Rt_*;
    tt = first Rt'frac_*;
    degenerate := f -> (
	M := map(frac Rt, Rt, {tt}|apply(numgens R,i->tt^(w#i)*Rt'frac_(i+1)));
	lift(numerator M f, Rt)
	);
    It := ideal apply(
	(RtoRt I)_*
	| 
	{RtoRt a + t * RtoRt b}
	,
	degenerate
	);
    alarm 2; -- number of seconds
    try ( 
	It'sat := ideal mingens saturate(It,t); 
	print "-- saturation";
	print toString It'sat;
	assert(numgens It'sat == numgens It)
	)
    else (
	print "validateRay: saturation did not finish... or not a complete intersection";
	It'sat = ideal apply(It_*, f -> (while f%t == 0 do f = f//t; f));   
	print "-- pseudo-saturation";
	print toString It'sat;
	);
    alarm 0;
    RtCC := CC[gens Rt];
    F := ((map(RtCC, Rt, gens RtCC)) It'sat)_*; 
    assert(numgens Rt == #F+1);
    Fepsilon := F|{RtCC_0-o.Epsilon}; -- target system: Epsilon = 0 by default
    done := false; count := 10;
    while not done do (
	count = count - 1;
	(S,solsS) := solveGenericSystemInTorus(F|{RtCC_0-1}); -- we solve a random system with the same support
	<< "mixed volume = " << #solsS << endl;
	done = all(solsS, s->isFullNumericalRank evaluate (jacobian polySystem S, s)) or (count == 0)
	);    
    solsEpsilon := track(S,Fepsilon,solsS,CorrectorTolerance=>1e-6,EndZoneFactor=>1e-24,Software=>o.Software);
    scan(solsEpsilon, s -> if status s === null then s.SolutionStatus = Regular); -- remove when bertiniTrackHomotpy is fixed !!!
    solsEpsilon = select(solsEpsilon, s->status s =!= Infinity);
    if o.Software === PHCPACK then postProcessPHCpack solsEpsilon;
    << "#solutions after the first step (t=epsilon): " << #solsEpsilon << endl;
    if o.Epsilon == 0 then sols := solsEpsilon
    else ( -- compute the limit when t->0
	F0 := F|{RtCC_0};
	sols = if #solsEpsilon == 0 then solsEpsilon
	else track(Fepsilon,F0,solsEpsilon,CorrectorTolerance=>1e-6,EndZoneFactor=>1e-24,Software=>o.Software);
	); 
    if o.Software === PHCPACK then postProcessPHCpack sols;
    good'pts := positions(sols, s->status s =!= Infinity);
    << "#solutions after the second step (t=0): " << #good'pts << endl;
    if o.Software === M2engine then (
    	print toString S;
    	print toString (solsS_good'pts / coordinates);
	);
    sols_good'pts
    )

-- returns (m, pts, pos) where 
--   m is the _alleged_ multiplicity 
--   pts are _all_ limit points  
--   pos is the list of numbers of solutions that are not on the torus boundary
validateRay = method(Options=>{Software=>M2engine,Epsilon=>0})
validateRay (List,Ideal) := o -> (w,I) -> (
    w = -w; -- min/max switch
    R := ring I;
    i := position(w, a->abs a == 1);
    if i =!= null then (
    	binom := if w#i == 1 then (R_i,1_R) else (1_R,R_i)
	)
    else (
    	i = position(w, a->a!=0);
	j := position(w, a->gcd(a,w#i)==1);
	g := gcdCoefficients(w#i,w#j);
	binom = ( 1_R , 1_R );
	if g#1 > 0 then binom = (binom#0 * R_i^(g#1), binom#1)
	else binom = (binom#0, binom#1 * R_i^(-g#1));
	if g#2 > 0 then binom = (binom#0 * R_j^(g#2), binom#1)
	else binom = (binom#0, binom#1 * R_j^(-g#2))
	);
    sols := tropicalIntersectionMultiplicityTimesRayMultiplicity(w, I, binom, getSymbol "xxx",o);
    pos := positions(sols, s->all(drop(coordinates s,1), x->not areEqual(x,0_CC)));
    (#pos, sols, pos)
    )

-- tropical degree by projectivization
-- in: list of pairs (ray,multiplicity)
-- out: integer
tropicalDegree = method()
tropicalDegree List := L -> sum apply(L, wm -> (
	(w,m) := wm;
	w' := -w|{0};
	a := min w';
	w' = apply(w', x->x-a);
	m*w'
	)) 

-- intersect amoeba with a hyperplane (in log coordinates)
-- (note: paths are tracked with external software ONE AT A TIME; this may change later)
sampleAmoeba (Ideal, List, Number) := o -> (I,w,C) -> ( 
    R := ring I;
    monom1 := product(#w, i->if w#i>0 then R_i else 1);
    monom2 := product(#w, i->if w#i<0 then R_i else 1); {* exp(2*pi*ii*random RR)* *}
    sys := I_*|{monom1 - C*monom2}; 
    if o.Software =!= PHCPACK then (
    	sols := solveSystem(sys,PostProcess=>false,Software=>o.Software);
	) 
    else (
	done := false; count := 1;
	while not done do (
	    count = count - 1;
	    (S,solsS) := solveGenericSystemInTorus (sys);
	    << "mixed volume = " << #solsS << endl;
	    done = all(solsS, s->isFullNumericalRank evaluate (jacobian polySystem S, s)) or (count == 0)
	    );
        sols = if true -- o.Software =!= PHCPACK --!!!
                 then track(S, sys, solsS, InfinityThreshold=>INFINITY, Software=>o.Software)
		 else --track(S, sys, solsS, Software=>BERTINI);
		 track(S, sys, solsS, InfinityThreshold=>INFINITY, Software=>M2engine);
    	scan(sols, s -> if status s === null then s.SolutionStatus = Regular); -- remove when bertiniTrackHomotpy is fixed !!!
	);
    regsols := select(sols, s->status s == Regular);
    << #regsols << " initial solutions" << endl;
    sys2 := I_*|{monom1 - o.Factor*C*monom2};
    regsols2 := apply(regsols, s->(
	    sols2 := if o.Software =!= PHCPACK --!!!
	    then track(sys, sys2, {s}, InfinityThreshold=>INFINITY, Software=>o.Software)
	    else --track(sys, sys2, {s}, Software=>BERTINI);
	    	track(sys, sys2, {s}, InfinityThreshold=>INFINITY, Software=>M2engine);
    	    scan(sols2, s -> if status s === null then s.SolutionStatus = Regular); -- remove when bertiniTrackHomotpy is fixed !!!
    	    select(sols2, s->status s == Regular)
	    ));
    good'paths := positions(regsols2, sols->#sols>0);
    log1 := apply(regsols_good'paths, s->apply(coordinates s, c->log abs c));
    log2 := apply(regsols2_good'paths, s->apply(coordinates first s, c->log abs c));
    if #good'paths>0 then (
	print "--";
    	<< "-- success with " << last sys << endl;
	<< "number of surviving points =  " << #good'paths << endl;  
	);
    diff'log'sols := log2-log1;
    apply(diff'log'sols, diff'log -> ShortApproxCollinearVectorAlt(diff'log, o.Multiplier))
    )

postProcessPHCpack = method()
postProcessPHCpack List := sols -> scan(sols, 
    s->if any(coordinates s, x -> abs x > getDefault InfinityThreshold)
    then s.SolutionStatus=Infinity
    )

degreeOfCurveInTorus = method()
degreeOfCurveInTorus Ideal := I -> (
    R := ring I;
    RCC := CC[gens R];
    ICC := sub(I, RCC);
    F := ICC_* | {random(1,RCC)+1};
    done := false; count := 5;
    while not done and count>0 do (
	count = count - 1;
    	(S,solsS) := solveGenericSystemInTorus F;
    	<< "mixed volume = " << #solsS << endl;
    	done = all(solsS, s->isFullNumericalRank evaluate (jacobian polySystem S, s))
    	);    
    # track(S,F,solsS,Software=>BERTINI)
    )

-- DATA COLLECTION AND PROCESSING FUNCTIONS

-- In: rays and verification data 
-- Out: (rays,mults)
--     rays, a list of sequences of integers
--     mults, a list of integers
combineDuplicates = method()
combineDuplicates List := rays'and'v'all -> (
    rays := unique apply(rays'and'v'all, r->r.PrincipalVector); 
    mults := apply(rays, r->(
	L := select(rays'and'v'all, r'->r'.PrincipalVector==r);
	pts := {};
	scan(L,r'-> (
		gpts := r'.Points_(r'.GoodPoints);
	    	pts' := select(gpts, p'->all(pts,p->not areEqual(p,p',Tolerance=>0.0001)));
	    	if #pts'>0 then (
		    if #pts > 0 then (
			<< "-- extra points for ray " << r << endl;
		    	<< "  appending " << matrix sortSolutions pts' << endl;
			<< "  from      " << matrix sortSolutions gpts << endl;
		    	<< "  to        " << matrix sortSolutions pts << endl;
			);
		    pts =  pts | pts'; 
	    	    );
	    	));
	#pts
	));
    (rays,mults)
    )


-- In:
reformatValidatedRays = method()
reformatValidatedRays (List,List) := (new'rays,v) -> (
    good'rays := positions(v,i->first i > 0);
    new'rays'and'v := apply(good'rays, i->new HashTable from {
    	    PrincipalVector => new'rays#i,
    	    Multiplicity => first v#i,
    	    Points => v#i#1 / coordinates,
    	    GoodPoints => v#i#2 
    	    }
    	); 
    -- things to print -------------------------------------------------------
    scan(good'rays, i->(
    	    << "-- ray of multiplicity " << v#i#0 << endl; 
    	    << toSequence new'rays#i << endl;
    	    --    << "-- leading coeff's of Puiseux series" << endl;
    	    << "  " << VerticalList( (v#i#1)_(v#i#2) / (p->toSequence drop(coordinates p,1)) ) << endl;
    	    ));
    new'rays'and'v
    );

appendNewData = method()
appendNewData (String, List, List) := (filename, new'rays, validated'rays'and'v) -> (
    copyFile(filename, filename|"~");
    out := openOutAppend filename;
    out << "-- BATCH: " << flush; run ("date >> "|filename);
    out << "known'rays = known'rays | " << toExternalString new'rays << endl; 
    out << "rays'and'v = rays'and'v | " << toString validated'rays'and'v << endl;
    out << endl << close;
    )
end;

