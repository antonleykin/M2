
getSequenceSC = method(TypicalValue=>Thing)
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
  	    	    --print latticePoints newP;
	    	    ))
	    else bfe#i = 0));
    --print ("bfe"=>bfe);
    idMatrix:= diagonalMatrix(apply(k,i->1));
    --newP is now a matroid polytope by shifting by -bfe
    newP = affineImage(idMatrix,newP,-matrix transpose {toList bfe});
    scan(bfe,i->scan(bfe#i,j -> SCS=append(SCS, i)));
--    print("SCS"=>SCS);
    numFactors := k;
    coarsenPositions := {};
    previousA := 0;
    scan(k,i->(
    	    ai := for j to k-1 list if j<=i then 1 else 0;
--    	    print ai;
    	    --maxAi is the dimension of the projection to the first i+1 coordinates.  
    	    maxAi := max\\flatten\\flatten\flatten\entries\latticePoints affineImage(matrix{ai},newP);
    	    if previousA<maxAi
	    then coarsenPositions = append(coarsenPositions,i);
	    previousA = maxAi;
	    ));
    leadVG :=first coarsenPositions;	
    coarsenPositions = drop(coarsenPositions,1); 
    --print coarsenPositions;
    scan(#coarsenPositions,shift->(
	    SCS = append(SCS,{leadVG,coarsenPositions_shift-shift});
	    SCS = append(SCS,leadVG)));
    scan(k-1-#coarsenPositions,i->SCS = append(SCS,{0,1}));
    SCS=append(SCS,0);
    return  SCS
    )

end
restart
needs "step1.m2"
needs "step2-Alternate.m2"

declareVariable \ {x,y,z,w}
f = 1 + 2*x + 3*y^2+4*z^3+5*w^4
h = 1 + 2*x + 3*y+ 5*z + 7*w + 11*z*y + 13 * x* z + 17*x*w + 19*y*z+23*y*w +29*z*w+31*x*y*z+37*x*y*w+41*x*z*w+43*y*z*w+47*x*y*z*w
F = gateSystem(matrix{{x,y,z,w}},transpose gateMatrix{{f,h}})
G = new VariableGroup from {{0},{1},{2},{3}}
pt = point{{4,-2/3,-1.44419, .765304}}
P = multiaffineDimension(F,G,pt)
latticePoints P 
getSequenceSC (P)

declareVariable \ {d,x,y,z,w};
f = 1 + 2*x + 3*y^2+4*z^3+5*w^4;
h = 1 + 2*x + 3*y+ 5*z + 7*w + 11*z*y + 13 * x* z + 17*x*w + 19*y*z+23*y*w +29*z*w+31*x*y*z+37*x*y*w+41*x*z*w+43*y*z*w+47*x*y*z*w;
f3 = d;
F = gateSystem(matrix{{d,x,y,z,w}},transpose gateMatrix{{f,h,f3}});
G = new VariableGroup from {{0},{1},{2},{3},{4}};
pt = point{{0,4,-2/3,-1.44419, .765304}};
P = multiaffineDimension(F,G,pt);
latticePoints P ;
getSequenceSC (P)

X = {x,y,z,w,d}
declareVariable \ X
f = 1 + 2*x + 3*y^2+4*z^3+5*w^4
h = 1 + 2*x + 3*y+ 5*z + 7*w + 11*z*y + 13 * x* z + 17*x*w + 19*y*z+23*y*w +29*z*w+31*x*y*z+37*x*y*w+41*x*z*w+43*y*z*w+47*x*y*z*w
f3 = d
F = gateSystem(matrix{X},transpose gateMatrix{{f,h,f3}})
G = new VariableGroup from {{0},{1},{2},{3},{4}}
pt = point{{0,4,-2/3,-1.44419, .765304}}
P = multiaffineDimension(F,G,pt)
latticePoints P 
getSequenceSC (P)









