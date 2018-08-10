-- multi-affine and multi-projective witness sets 

MultiSlicingVariety = new Type of SlicingVariety
MultiSlicingVariety.synonym = "multi-affine slice"
multiSlicingVariety = method()
multiSlicingVariety(Ambient,List) := (A,MM) -> new MultiSlicingVariety from {"ambient"=>A, "maps"=>MM}
codim MultiSlicingVariety :=  {} >> o -> S -> rank@@target\S#"maps"
net MultiSlicingVariety := S -> net "slice of codim " | net codim S
map MultiSlicingVariety := o -> S -> transpose matrix{S#"maps"/transpose}

randomSlicingVariety(MultiAffineSpace,List) := (A,K) -> ( -- K = list of codimensions 
    if not (all(K,k->class k===ZZ and k>=0) and sum K>0) then error" codimensions are not nonnegative integers.";
    R := ring A;
    N := dim A;
    M := apply(#N,i->(
	    n := N#i; 
	    k := K#i;
	    if k==0 then map(R^0,R^1,0)
	    else (variables(i,A) * random(R^n,R^k) - matrix {toList (k:1_R)})
	    ));
    multiSlicingVariety( A, M/rationalMap )
    )

MultiAffineWSet = new Type of WSet
multiAffineWSet = method(TypicalValue=>MultiAffineWSet)
multiAffineWSet (PolySystem,SlicingVariety,List) := (F,S,pts) -> (
      R := ring F;
      -- assert isHomogeneous F.PolyMap; --!!! 
      new MultiAffineWSet from {
      	  "ambient" => multiAffineSpace R,
	  "equations" => F, 
      	  "slice" => S,
      	  "points" => pts
      	  }
      )
dim MultiAffineWSet := W -> codim W#"slice"
degree MultiAffineWSet := W -> # points W 

net MultiAffineWSet := W -> net "multiAffineWSet(dim=" | net dim W | ",deg=" | net degree W | ")" 

slicingVariety MultiAffineWSet := W -> W#"slice" 

-- this returns points "downstairs"
points MultiAffineWSet := W -> W#"points"

----------------------------------------------------------
-- Witness collection for a multi-projective variety 
-- should store a collection of multi-affine witness sets

WCollection = new Type of MutableHashTable
WCollection.synonym = "multi-projective witness set collection"

net WCollection := W -> if #W#"witnesses">0 then
"[dim="| dim W | " deg="| degree W |" in " | net ambient W | "]" else
"uninitialized WCollection"

wCollection = method(TypicalValue=>WCollection, 
    Options=>{Tolerance=>0.000001} -- used to determine whether the point is "at infinity": see toChart
    )
wCollection(Ambient,PolySystem) := o -> (A,F) ->
  new WCollection from {
      "ambient" => A,
      "equations" => F,
      "witnesses" => null,
      Tolerance => o.Tolerance
      }


--Need to decide if we want to return the partial dimension set or report an error. 
dim WCollection := W ->  if W#"witnesses"=!=null then keys W#"witnesses" else error "WCollection not initialized"  
codim WCollection := {} >> o -> W ->  apply(dim W, d -> dim ambient W - d)
ambient WCollection := W -> W#"ambient"
WCollection _ Sequence := (W,d) -> if member(d,witnessKeys W) then W#"witnesses"#d else null
WCollection _ List := (W,d) ->  (
    r := select(witnessKeys W, k -> first k == d);
    if #r>0 then W_(first r) 
    else null
    ) 

-- degree WCollection := W -> #W.Points
-- points WCollection := W -> W.Points

addWSet = method()
addWSet (WCollection, MultiSlicingVariety, List) := (W,S,pts) ->(
    if #pts == 0 then error "refusing to add an empty witness set"; 
    if W#"witnesses"===null (	
      	-- do computation of the dim set
	-- (described in the "polymatroid" literature?)
      	);
    W#"witnesses"#(codim S) = append(
    	W#"witnesses"#(codim S),
    	multiAffineWSet(W#"equations",S,pts)
	)    	
    )

--MultiprojectiveNAGtypes
--First we define multi-affine witness sets and collections. 
--A witness set consists 
TEST ///--multiaffine example
restart
debug needsPackage "NAGtypes"
debug needsPackage "NumericalAlgebraicGeometry"
errorDepth = 2
A = multiAffineSpace(CC_53,{1,1},x)
use ring A 
F = polySystem {x_(1,1)^2-x_(0,1)}
S=randomSlicingVariety(A,{1,0})
peek S
pts = {point{{1,1}},point{{1,-1}}}
W = wCollection(A,F) 
peek W
peek  W#"equations"
addWSet(W,S,pts)--Here we should also be able to do addWSet(WCollection,WSet)
addWSet(W,S,pts)
W#"witnesses"#{1,0}
peek W
first pairs( W#"witnesses")
dim W

-- how to perturb a slice?
-- how to "populate" a w.set? (suppose you have a partial w.set)
-- how to populate a neighboring dimension w.set ("hop" or "hopscotch"?) 
-- how to do a trace test?
-- how to check irreducibility?

///



TEST ///--multiprojective example
restart
debug needsPackage "NAGtypes"
debug needsPackage "NumericalAlgebraicGeometry"
errorDepth = 2
A = multiProjectiveSpace(CC_53,{1,1},x)
use ring A 
-- multi=homogenized parabola y-z^2=0 where y=x_(0,0) and z=x_(1,0)
F = polySystem {x_(0,0)*x_(1,1)^2-x_(1,0)^2*x_(0,1)}
H = {x_(1,1)-1,x_(0,1)-1}
S = multiSlicingVariety(A, {rationalMap matrix{{x_(0,0)-1}}, rationalMap map((ring A)^0,(ring A)^1,0)})
pts = {point{{1,1,1,1}},point{{1,1,-1,1}}}
W = wCollection(A,F) 
peek W
peek  W#"equations"

addWSet(W,H,S,pts)
peek W
first pairs( W#"witnesses")

dim W
W_{1,0}
assert(W_{0,1} === null)
///