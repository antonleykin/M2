
SequenceSC = new Type of MutableList

getSequenceSC = method(TypicalValue=>Thing)
getSequenceSC (Polyhedron) := (P)->(
    SCS := {};
    if isEmpty(P) then error"P is empty." ;
    k:=ambDim P;
    bfe:=new MutableList from {};
    newP:=P;
    scan(k,i->(
	    --ith basis vector
    	    ei := for j to k-1 list if i==j then 1 else 0;
	    --project to first coordinate
    	    maxEi := latticePoints affineImage(matrix{ei},newP);
	    maxEi = maxEi/entries/flatten/flatten//flatten//max;
    	    --if maxEi is greater than one then Bertini's theorem applies and we can slice.
	    if maxEi > 1 then (
		bfe#i = maxEi-1;
	    	if bfe#i>0 and not isEmpty(newP) 
	    	then (
		    Q :=polyhedronFromHData(-matrix{ei},-matrix {{bfe#i}}); -- newM \leq e_I  
  	    	    newP = intersection(newP,Q);
  	    	    print latticePoints newP;
	    	    ))
	    else bfe#i = 0));
    print ("bfe"=>bfe);
    idMatrix:= diagonalMatrix(apply(k,i->1));
    newP=affineImage(idMatrix,newP,-matrix transpose {toList bfe});
    scan(bfe,i->scan(bfe#i,j -> SCS=append(SCS, i)));
    print 1;
    dimLowerBound = 1;
    scan(k,i->(
    	    ai := for j to k-1 list if j<=i then 1 else 0;
    	    maxAi := latticePoints affineImage(matrix{ai},newP);
	    maxAi = maxAi/entries/flatten/flatten//flatten//max;
    	    print ("i"=>i);
	    print maxAi;
	    print SCS;
    	    if maxAi>dimLowerBound then SCS=append(SCS,{i-1,i});
	    dimLowerBound = max(maxAi,dimLowerBound)));	    
    SCS=append(SCS,0);
    return toSequence SCS
    )

end
restart
needs "multiaffineDimension.m2"
needs "step2.m2"

declareVariable \ {x,y,z,w}
f = 1 + 2*x + 3*y^2+4*z^3+5*w^4
h = 1 + 2*x + 3*y+ 5*z + 7*w + 11*z*y + 13 * x* z + 17*x*w + 19*y*z+23*y*w +29*z*w+31*x*y*z+37*x*y*w+41*x*z*w+43*y*z*w+47*x*y*z*w
F = gateSystem(matrix{{x,y,z,w}},transpose gateMatrix{{f,h}})
G = new VariableGroup from {{0},{1},{2},{3}}
pt = point{{4,-2/3,-1.44419, .765304}}
P = multiaffineDimension(F,G,pt)
latticePoints P 
getSequenceSC (P)


declareVariable \ {d,x,y,z,w}
f = 1 + 2*x + 3*y^2+4*z^3+5*w^4
h = 1 + 2*x + 3*y+ 5*z + 7*w + 11*z*y + 13 * x* z + 17*x*w + 19*y*z+23*y*w +29*z*w+31*x*y*z+37*x*y*w+41*x*z*w+43*y*z*w+47*x*y*z*w
f3 = d
F = gateSystem(matrix{{d,x,y,z,w}},transpose gateMatrix{{f,h,f3}})
G = new VariableGroup from {{0},{1},{2},{3},{4}}
pt = point{{0,4,-2/3,-1.44419, .765304}}
P = multiaffineDimension(F,G,pt)
latticePoints P 
getSequenceSC (P)
