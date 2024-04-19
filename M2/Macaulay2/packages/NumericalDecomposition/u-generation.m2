DBG = 3;

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
     for f in F do V = uHypersurfaceSection(V,f,o); -- intersect with hypersurface	 
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
	newComponent := hypersurfaceSection(numericalVariety {cOut}, f); -- todo: replace with u-generation call
	if newComponent =!= null then -- null is returned when the intersection is empty (otherwise it is of dimension one less) 
	-*
	using old code that may return several components
	*-
	newComponentsFromOldNumericalVariety := components newComponent;
	scan(newComponentsFromOldNumericalVariety, c -> insert(c,result));
--	insert(newComponent, result); -- todo: reinstate this syntax
	);
    result
    ) -- end uHypersurfaceSection(WitnessSet,...)


isPointOnAnyComponent = method()
isPointOnAnyComponent(AbstractPoint,HashTable) := (p,H) -> any(keys H, d -> any(keys H#d, k -> isOn(p,H#d#k)))

splitWitness = method(TypicalValue=>Sequence, Options =>{Tolerance=>1e-6})
splitWitness (WitnessSet,RingElement) := Sequence => o -> (w,f) -> (
-- splits the witness set into two parts: one contained in {f=0}, the other not
-- IN:  comp = a witness set
--      f = a polynomial
-- OUT: (w1,w2) = two witness sets   
     w1 := {}; w2 := {};
     for x in points w do 
	 if residual(matrix {{f}}, matrix x) < o.Tolerance 
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


TEST ///
restart
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
