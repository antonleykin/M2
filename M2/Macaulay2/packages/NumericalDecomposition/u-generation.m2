uGeneration = method()
uGeneration System??? := NumericalVariety => F -> (
-- solves a system of polynomials equation-by-equation
-- IN:  F = list of polynomials
-- OUT: a NumericalVariety
     o = fillInDefaultOptions o;
     
     checkCCpolynomials F;
     F = toCCpolynomials(F,53);
     
     R := ring F#0;
     V := numericalAffineSpace R; -- current solution components
     for f in F do V = uHypersurfaceSection(V,f,o); -- intersect with hypersurface	 
     V
     )

uHypersurfaceSection = method(Options =>{})
uHypersurfaceSection(NumericalVariety,System???) := o -> (c1,f) -> (
    if DBG>1 then << "-- uHypersurfaceSection: "<< endl << 
    "NumericalVariety: " << c1 << endl <<
    "hypersurface: " << f << endl;
    d := sum degree f;
    R := ring f;
    -- if getDefault Normalize then f = f/sqrt(numgens R * BombieriWeylNormSquared f); 
    c2 := new MutableHashTable; -- new components
    for comp in components c1 do (
	if DBG>2 then << "*** processing component " << peek comp << endl;
	result := uHypersurfaceSection(c1,f);
        -- postprocess into c2
    	);
    numericalVariety flatten apply(keys c2, i->apply(keys c2#i, j->c2#i#j))
    )

uHypersurfaceSection =method()
uHypersurfaceSection(WitnessSet,System???) := (comp,f) -> (
    (cIn,cOut) := splitWitness(comp,f); 
    if cIn =!= null then (
	if DBG>2 then << "( regeneration: " << net cIn << " is contained in V(f) for" << endl <<  
	<< "  f = " << f << " )" << endl;
	scan(
	    partitionViaDeflationSequence( cIn.Points, 
		polySystem(equations polySystem cIn | slice cIn) 
		-- this is overdetermined!
		),
	    pts -> (
		oldW := witnessSet(cIn.Equations, cIn.Slice, pts);
		if DBG>2 then << "   old component " << peek oldW << endl;
		check oldW;    
		insertComponent(oldW,c2);
		)
	    );
	); 
    if cOut =!= null 
    and dim cOut > 0 -- 0-dimensional components outside V(f) discarded
    then (
	s := cOut#Slice;
	-- RM := (randomUnitaryMatrix numcols s)^(toList(0..d-2)); -- pick d-1 random orthogonal row-vectors (this is wrong!!! is there a good way to pick d-1 random hyperplanes???)
	RM := random(CC^(d-1),CC^(numcols s));
	dWS := {cOut} | apply(d-1, i->(
		newSlice := RM^{i} || submatrix'(s,{0},{}); -- replace the first row
		moveSlice(cOut,newSlice,Software=>o.Software)
		));
	slice' := submatrix'(comp#Slice,{0},{});
	local'regular'seq := equations polySystem comp;
	S := polySystem( local'regular'seq
	    | { product flatten apply( dWS, w->sliceEquations(w.Slice^{0},R) ) } -- product of linear factors
	    | sliceEquations(slice',R) );
	T := polySystem( local'regular'seq
	    | {f}
	    | sliceEquations(slice',R) );
	-- deflate if singular
	P := first comp.Points;
	if status P =!= Singular then targetPoints := track(S,T,flatten apply(dWS,points), 
	    NumericalAlgebraicGeometry$gamma=>exp(random(0.,2*pi)*ii),
	    Software=>o.Software)
	else (
	    seq := P.DeflationSequenceMatrices;
	    S' := squareUp(deflate(S, seq), squareUpMatrix P.cache.LiftedSystem); -- square-up using the same matrix
	    T' := squareUp(deflate(T, seq), squareUpMatrix P.cache.LiftedSystem); -- square-up using the same matrix
	    S'sols := flatten apply(dWS,W->apply(W.Points,p->p.cache.LiftedPoint));
	    T'.PolyMap = (map(ring S', ring T', vars ring S')) T'.PolyMap; -- hack!!!: rewrite with trackHomotopy
	    lifted'w' := track(S',T',S'sols, 
		NumericalAlgebraicGeometry$gamma=>exp(random(0.,2*pi)*ii), Software=>o.Software);
	    targetPoints = apply(lifted'w', p->(
		    q := project(p,T.NumberOfVariables);
		    q.cache.SolutionSystem = T;
		    q.cache.LiftedSystem = T';
		    q.cache.LiftedPoint = p;
		    q.cache.SolutionStatus = Singular;
		    q
		    ));
	    );    
	LARGE := 100; ---!!!
	refinedPoints := refine(T, targetPoints, 
	    ErrorTolerance=>DEFAULT.ErrorTolerance*LARGE,
	    ResidualTolerance=>DEFAULT.ResidualTolerance*LARGE,
	    Software=>o.Software);
	regPoints := select(refinedPoints, p->p.cache.SolutionStatus===Regular);
	singPoints := select(refinedPoints, p->p.cache.SolutionStatus===Singular);
	targetPoints = if o.Output == Regular then regPoints else regPoints | solutionsWithMultiplicity singPoints;
	if DBG>2 then << "( regeneration: " << net cOut << " meets V(f) at " 
	<< #targetPoints << " points for" << endl 
	<< "  f = " << f << " )" << endl;
	f' := ideal (equations comp | {f});
	nonJunkPoints := select(targetPoints, p-> not isPointOnAnyComponent(p,c2)); -- this is very slow		    
	scan(partitionViaDeflationSequence(nonJunkPoints,T),
	    pts -> (
		newW := witnessSet(f',slice',selectUnique(pts, Tolerance=>1e-4));--!!!
		if DBG>2 then << "   new component " << peek newW << endl;
		check newW;    
		insertComponent(newW,c2);
		)
	    )
	) 
    ) -- end uHypersurfaceSection(WitnessSet,...)
    

TEST ///
R = CC[x,y,z]
numericalAffineSpace R
NV = numericalAffineSpace R
comps = components NV
C = first comps

sph = (x^2+y^2+z^2-1); 

regeneration {sph*(x-1)*(y-x^2), sph*(y-2)*(z-x^3)}

///
