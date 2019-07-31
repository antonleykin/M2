needs "common.m2"
-- NEW
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
    (gateMatrix{toList ParamList}, transpose gateMatrix{toList Leqs})
    )
end

restart
needs "step3.m2"

-- (mostly) COPIED FROM OCTAHEDRON
declareVariable \ {x,y,z,w}
f = 1 + 2*x + 3*y^2+4*z^3+5*w^4
h = 1 + 2*x + 3*y+ 5*z + 7*w + 11*z*y + 13 * x* z + 17*x*w + 19*y*z+23*y*w +29*z*w+31*x*y*z+37*x*y*w+41*x*z*w+43*y*z*w+47*x*y*z*w
F = gateSystem(
    matrix{{x,y,z,w}},
    transpose gateMatrix{{f,h}}
    )
G = new VariableGroup from {{0},{1},{2},{3}}
pt = point{{4,-2/3,-1.44419, .765304}}
evaluate(F,pt)

(m,n,JF)=createJacobian F
Jpt = evaluateJacobian(pt,m,n,JF)

P = convexHull transpose matrix{
    {1,1,0,0},
    {1,0,1,0},
    {1,0,0,1},
    {0,1,1,0},
    {0,1,0,1},
    {0,0,1,1}
    }

SCseq =({2,3},2,{0,1},{0,1},0)


(Params, L) = makeSliceSystem(F,G,SCseq)
MasterGS = gateSystem(
    Params,
    vars F,
    (gateMatrix F) || L
    )