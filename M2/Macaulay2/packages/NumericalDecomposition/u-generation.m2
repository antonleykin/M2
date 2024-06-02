debug NumericalAlgebraicGeometry;
DBG = 0;

uGeneration = method(Options=>{})
uGeneration List -*System???*- := NumericalVariety => o -> F -> (
-- solves a system of polynomials equation-by-equation
-- IN:  F = list of polynomials
-- OUT: a NumericalVariety
     -- o = fillInDefaultOptions o;
     
     --checkCCpolynomials F;
     --F = toCCpolynomials(F,53);
     
     R := ring F#0;
     V := numericalAffineSpace R; -- current solution components
     for f in F do (
	 print ("before: V" => V);
	 print ("before: f" => f);
	 V = uHypersurfaceSection(V,f,o); -- intersect with hypersurface
	 print ("after: V" => V);
	 );
     V
     )

uHypersurfaceSection = method(Options =>{})
uHypersurfaceSection(NumericalVariety,RingElement-*System???*-) := o -> (c1,f) -> (
    if DBG>1 then << "-- uHypersurfaceSection: "<< endl << 
    "NumericalVariety: " << c1 << endl <<
    "hypersurface: " << f << endl;
    d := sum degree f;
    R := ring f;
    -- if getDefault Normalize then f = f/sqrt(numgens R * BombieriWeylNormSquared f); 
    c2 := new NumericalVariety; -- new components
    for comp in components c1 do (
	if DBG>2 then << "*** processing component " << peek comp << endl;
	result := uHypersurfaceSection(comp,f);
        -- postprocess into c2 ...
	for newComponent in components result do insert(newComponent, c2);
    	);
    c2
    )
uHypersurfaceSection(WitnessSet,RingElement-*System???*-) := o -> (comp,f) -> (
    (cIn,cOut) := splitWitness(comp,f); 
    result := numericalVariety{};
    if cIn =!= null then (
	if DBG>2 then << "( u-generation: " << net cIn << " is contained in V(f) for" << endl <<  
	<< "  f = " << f << " )" << endl;
	insert(cIn,result);
	); 
    if cOut =!= null and dim cOut > 0 -- 0-dimensional components outside V(f) discarded
    then (
	newComponent := hypersurfaceSection(cOut, numericalVariety(
		if cIn===null then {} else {cIn}
		), f); -- todo: replace with u-generation call
	if newComponent =!= null then -- null is returned when the intersection is empty (otherwise it is of dimension one less) 
	-*
	using old code that may return several components
	*-
	-- newComponentsFromOldNumericalVariety := components newComponent;
	-- scan(newComponentsFromOldNumericalVariety, c -> insert(c,result));
	insert(newComponent, result);
	);
    result
    ) -- end uHypersurfaceSection(WitnessSet,...)


--isPointOnAnyComponent = method()
isPointOnAnyComponent(AbstractPoint,NumericalVariety) := (p,V) -> any(components V, c-> isOn(p,c))

--splitWitness = method(TypicalValue=>Sequence, Options =>{Tolerance=>1e-6})
splitWitness (WitnessSet,RingElement) := Sequence => o -> (w,f) -> (
-- splits the witness set into two parts: one contained in {f=0}, the other not
-- IN:  comp = a witness set
--      f = a polynomial
-- OUT: (w1,w2) = two witness sets   
     w1 := {}; w2 := {};
     for x in points w do 
	 if residual(matrix {{f}}, matrix x) < 1e-6 -- o.Tolerance 
	 then w1 = w1 | {x}
	 else w2 = w2 | {x};   
     ( if #w1===0 then null else witnessSet(ideal equations w + ideal f, 
	     -* this is "stretching" the convention that this has to be a complete intersection *-
	     w.Slice, w1), 
       if #w2===0 then null else witnessSet(w.Equations, w.Slice, w2) 
       )
   )

insert(WitnessSet,NumericalVariety) := (comp,V) -> (
    d := dim comp; 
    if not V#?d then V#d = {};
    V#d = V#d | {comp};
    )      


---------------------------
-------- OLD CODE below !!!

hypersurfaceSection(WitnessSet,NumericalVariety,RingElement) := o -> (cOut,varietyToAvoid,f) -> (
    if DBG>1 then << "-- hypersurfaceSection: "<< endl << 
    "equidimensional component: " << cOut << endl <<
    "hypersurface: " << f << endl;
    d := sum degree f;
    R := ring f;
    if DBG>2 then << "*** processing component " << peek cOut << endl;
    if dim cOut == 0 -- 0-dimensional components outside V(f) discarded
    then null else (
	s := cOut#Slice;
	-- RM := (randomUnitaryMatrix numcols s)^(toList(0..d-2)); -- pick d-1 random orthogonal row-vectors (this is wrong!!! is there a good way to pick d-1 random hyperplanes???)
	RM := random(CC^(d-1),CC^(numcols s));
	dWS := {cOut} | apply(d-1, i->(
		newSlice := RM^{i} || submatrix'(s,{0},{}); -- replace the first row
		moveSlice(cOut,newSlice,Software=>o.Software)
		));
	slice' := submatrix'(cOut#Slice,{0},{});
	local'regular'seq := equations polySystem cOut;
	    

	S := polySystem( local'regular'seq
	    | { product flatten apply( dWS, w->sliceEquations(w.Slice^{0},R) ) } -- product of linear factors
	    | sliceEquations(slice',R) );
	T := polySystem( local'regular'seq
	    | {f}
	    | sliceEquations(slice',R) );
    
	P := first cOut.Points;

	targetPoints := if any(cOut.Points, p->status p === Singular) then (
	    << "warning: singular points encountered in the INPUT by hypersurfaceSection(WitnessSet,...): " << endl;
	    -- do nothing: discard
	    {}
	    )
	else trackHomotopy(segmentHomotopy(S,T,NumericalAlgebraicGeometry$gamma=>exp(random(0.,2*pi)*ii)),
		flatten apply(dWS,points), Software=>o.Software);
	-*
	LARGE := 100; ---!!!
	refinedPoints := refine(T, targetPoints, 
	    ErrorTolerance=>DEFAULT.ErrorTolerance*LARGE,
	    ResidualTolerance=>DEFAULT.ResidualTolerance*LARGE,
	    Software=>o.Software);
	*-
	refinedPoints := targetPoints;
	regPoints := select(refinedPoints, p->p.cache.SolutionStatus===Regular);
	singPoints := select(refinedPoints, p->p.cache.SolutionStatus===Singular);
	print (#regPoints, #singPoints);
	if #singPoints > 0 then
 	  << "warning: singular points encountered in the COMPUTATION by hypersurfaceSection(WitnessSet,...): " << endl;
	targetPoints = regPoints;
	if DBG>2 then << "( regeneration: " << net cOut << " meets V(f) at " 
	<< #targetPoints << " points for" << endl 
	<< "  f = " << f << " )" << endl;
	f' := ideal (equations cOut | {f});
	nonJunkPoints := select(targetPoints, p-> not isPointOnAnyComponent(p,varietyToAvoid)); -- this is very slow		    
	junkPoints := select(targetPoints, p-> isPointOnAnyComponent(p,varietyToAvoid)); -- this is very slow		    
	if #junkPoints > 0 and DBG>2 then << "   #junk points = " << #junkPoints << endl;
	newW := witnessSet(f',slice',selectUnique(nonJunkPoints, Tolerance=>1e-4));--!!!
	if DBG>2 then << "   new component " << peek newW << endl;
	check newW;
	newW
	)
    )

TEST ///
restart
--errorDepth = 0
needsPackage "NumericalDecomposition"
R = CC[x,y,z]
numericalAffineSpace R
NV = numericalAffineSpace R
comps = components NV
C = first comps
sph = (x^2+y^2+z^2-1); 
--NAGtrace 4
uGeneration {sph*(y-x^2), sph*(z-x^3)}
///


TEST ///
restart
needsPackage "NumericalDecomposition"
R = CC[x,y,z]
numericalAffineSpace R
NV = numericalAffineSpace R
comps = components NV
C = first comps
sph = (x^2+y^2+z^2-1); 
uGeneration {sph*(x-1)*(y-x^2), sph*(y-2)*(z-x^3)}
///


end--
restart
needsPackage "NumericalDecomposition"
n = 6
R = CC[x_1..x_n,y_1..y_n]
M = genericMatrix(R,n,2)
F = for i from 1 to n-1 list det M^{i-1,i}
V = uGeneration F
decompose \ components V
