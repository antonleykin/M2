-- -*- coding: utf-8 -*-
-- licensed under GPL v2 or any later version
newPackage(
     "WitnessCollections",
     Version => "1.12",
     Date => "Summer 2018",
     Headline => "witness collections representing multiaffine varieties",
     HomePage => "",
     AuxiliaryFiles => false,
     Authors => {
	  {Name => "Jose Rodriguez", Email => "jose@madison"},
	  {Name => "Anton Leykin", Email => "leykin@math.gatech.edu"}
	  },
     PackageImports => {"NumericalAlgebraicGeometry", "Polyhedra"},
     DebuggingMode => false 
     --DebuggingMode => true 
     )

-------------------------------------------------------------------------
-- To do: 
-- how to perturb a slice?
-- how to "populate" a w.set? (suppose you have a partial w.set)
-- how to populate a neighboring dimension w.set ("hop" or "hopscotch"?) 
-- how to do a trace test?
-- how to check irreducibility?
-------------------------------------------------------------------------


debug needsPackage "NAGtypes"

--Input: PolySystem and a point. We assume it is a general point, 
---I.e., it only lives on one irreducible component.
--Output: dimension set of a variety that is an irreducible component of the PolySystem containing that point. 
--Let's also assume the variety is d dimensional as an affine space or maybe not. 

dimensionPolytope = method()

dimensionPolytope(Point,MultiAffineWSet):= (pt,MAWS) ->dimensionPolytope(pt,ambient MAWS,equations MAWS)  
dimensionPolytope(Point,WCollection):= (pt,WC) ->dimensionPolytope(pt,ambient WC,equations WC)  
dimensionPolytope(Point,Ambient,PolySystem):= (pt,A,F) ->( 
    jacF := jacobian F; --This ia a matrix. 
    thePartialJacs := apply(variables(A),v->sub(diff (v,F.PolyMap),matrix pt));
    theJacs := apply(subsets thePartialJacs,
	I -> if #I>0 
    	then matrix{I}
	);--drop one for the empty set
--    print 1;
  theRanks:=drop(apply(theJacs,j->if j=!=null 	      then (numericalRank j)),1);--drop 1 for the empty set
  theVarGroups:= apply(  drop(subsets(#dim A),1),I->apply(#dim A,i->if member(i,I) then 1 else 0 ));
--    print 2;
  M:=matrix drop(theVarGroups,-1);
--    print M;
  v:=transpose matrix {drop(theRanks,-1)};
--    print v;
  N:=matrix {theVarGroups_-1};
--    print N;
  w:=transpose matrix{{theRanks_-1}};
--    print w;
--  print (M,v,N,w);
  polyhedronFromHData(M,v,N,w)
  --apply(latticePoints P,p->flatten entries (transpose matrix {dim A}-p))
    ) 




---Input: A witness collection of a multi affine variety--
--it is really a partial witness collection. 
---(Possibly irreducible and non equidimensional in the polytope sense; 
--we assume nothing).
--Output: Describes a union of components--
--where each component is equidimensional in the dimension polytope sense.
--(These do not have to be irreducible)


memberEqualEqualOut:=(x,L)->(
  theOut:=null;
  scan(L,i->if i==x then (theOut=i; break )) ;
  return theOut );

dimensionPartition = method()
dimensionPartition MultiAffineWSet := WS -> ( 
  P:=new MutableHashTable;
  apply(points WS,p->(
    dp:=dimensionPolytope(p,WS);
    theKey:=memberEqualEqualOut(dp,keys P);
    if null=== theKey then P#dp={p}
    else P#theKey=append(P#theKey,p)));
  print apply(keys P,i->vertices i);
  print  netList pairs P;
  theValues:=apply(keys P,
    dimpoly->multiAffineWSet(equations WS,slicingVariety WS, P#dimpoly,
      "dimension polytope"=>dimpoly
      )  );
  return new MutableHashTable from transpose {keys P,theValues}
  )
dimensionPartition (WCollection,List) := (WC,sliceType) -> ( 
  allSlice:=apply(WC#"witnesses"#sliceType,i->dimensionPartition i);
  P:=new MutableHashTable;
  print ( allSlice/pairs);
  apply(allSlice,gws->apply(keys gws,dp->(
    theKey:=memberEqualEqualOut(dp,keys P);
    if null=== theKey then P#dp={gws#dp}
    else P#theKey=append(P#theKey,gws#dp))));
  return P
  )

dimensionPartition WCollection := WC -> ( 
  Q:=new MutableHashTable;
  allTypes:=apply(WC#"witnesses"//keys,st->(--st=sliceType
    gws:=dimensionPartition(WC,st);
    apply(keys gws,dp->(
      theKey:=memberEqualEqualOut(dp,keys Q);
      if null=== theKey 
      then (Q#dp=wCollection(ambient WC,equations WC);
	  theKey=dp);
      apply(gws#dp,ws->addWSet(Q#theKey,slicingVariety ws, points ws)
	  )))));
  return Q
  )


TEST ///--dimension partition example
restart
debug needsPackage "NAGtypes"
debug needsPackage "WitnessCollections"
needsPackage "Polyhedra"
errorDepth = 2
A = multiAffineSpace(CC_53,{1,1},symbol x)
use ring A 
x=first flatten entries first variables A
y=first flatten entries last variables A
F = polySystem {x*y*(x^2-y)}
S1 = randomSlicingVariety(A,{1,0})
S2 = randomSlicingVariety(A,{0,1})
S1#"maps"=replace(0,matrix {{x-3_CC}},S1#"maps") 
S2#"maps"=replace(1,matrix {{y-4_CC}},S2#"maps") 
pts1 = {point{sub(matrix{{3,0}},CC)},point{sub(matrix{{3,9}},CC)}}
pts2 = {point{sub(matrix{{0,4}},CC)},point{sub(matrix{{2,4}},CC)},point{sub(matrix{{-2,4}},CC)}}
W1 = multiAffineWSet(F,S1,pts1)
W2 = multiAffineWSet(F,S2,pts2)
keys W1
assert (dim W1==codim S1)
assert(dim W2=={0,1})
assert(degree W1==2)
assert(degree W2==3)
thePartition= dimensionPartition(W2)
peek thePartition
thePartition#(first keys thePartition)
dim first keys thePartition
dim last keys thePartition
assert(#keys thePartition==2)
--
WC=wCollection(A,F)
addWSet(WC,S1,pts1)
addWSet(WC,S2,pts2)
addWSet(WC,S1,pts1)
addWSet(WC,S2,pts2)
keys WC
P=dimensionPartition( WC,{0,1})
P#(first keys P)
keys P
P=dimensionPartition( WC)
P#(last keys P)#"witnesses"//peek

///




TEST ///--multiaffine example
restart
debug needsPackage "NAGtypes"
debug needsPackage "WitnessCollections"
needsPackage "Polyhedra"
errorDepth = 2
A = multiAffineSpace(CC_53,{3,2},symbol x)
use ring A 
variables A
F = polySystem {x_(1,1)^2-x_(0,1)}
S = randomSlicingVariety(A,{1,0})
assert(codim S == {1,0})
assert(dim A == codim S + dim S)
P = point{apply(dim ring A,i->1_CC)}
W = wCollection(A,F) 
peek W
assert( latticePoints dimensionPolytope(P,W) / entries // sort ==  {{{0}, {1}}, {{1}, {0}}} )

latticePoints dimensionPolytope(
    point random(CC^1,CC^(dim ring A))
    ,W)

codimG = {1,2}
G = randomSlicingVariety(A,codimG)
latticePoints dimensionPolytope(first pts,wCollection(A,polySystem map G)) 


A = multiAffineSpace(CC_53,{2,3},symbol x)
use ring A 
F = polySystem {random({1,1},ring A),random({1,0},ring A)}

///





TEST ///
---Affine funtf 3v5
restart
debug needsPackage "NAGtypes"
debug needsPackage "WitnessCollections"

R=CC[{w1v1, w1v2, w1v3}]**CC[{w2v1, w2v2, w2v3}]**CC[{w3v1, w3v2, w3v3}]**CC[{w4v1, w4v2, w4v3}]**CC[{w5v1, w5v2, w5v3}]
--R=CC[w1v1]**CC[w1v2]**CC[ w1v3]**CC[w2v1]**CC[ w2v2]**CC[ w2v3]**CC[w3v1]**CC[w3v2]**CC[ w3v3]**CC[w4v1]**CC[w4v2]**CC[ w4v3]**CC[w5v1]**CC[ w5v2]**CC[ w5v3]
--R=CC[{w1v1, w2v1, w3v1, w4v1, w5v1}]**CC[ {w1v2, w2v2, w3v2, w4v2, w5v2}] ** CC[{w1v3, w2v3, w3v3, w4v3, w5v3}]

jade0 = w1v1^2+w2v1^2+w3v1^2+w4v1^2+w5v1^2-1 ; 
jade1 = w1v1*w1v2+w2v1*w2v2+w3v1*w3v2+w4v1*w4v2+w5v1*w5v2 ; 
jade2 = w1v1*w1v3+w2v1*w2v3+w3v1*w3v3+w4v1*w4v3+w5v1*w5v3 ; 
jade3 = w1v2*w1v3+w2v2*w2v3+w3v2*w3v3+w4v2*w4v3+w5v2*w5v3 ; 
jade4 = w1v1^2+w2v1^2+w3v1^2+w4v1^2+w5v1^2-w1v2^2-w2v2^2-w3v2^2-w4v2^2-w5v2^2 ; 
jade5 = w1v2^2+w2v2^2+w3v2^2+w4v2^2+w5v2^2-w1v3^2-w2v3^2-w3v3^2-w4v3^2-w5v3^2 ; 
jade6 = w1v1^2-w2v1^2+w1v2^2-w2v2^2+w1v3^2-w2v3^2 ; 
jade7 = w2v1^2-w3v1^2+w2v2^2-w3v2^2+w2v3^2-w3v3^2 ; 
jade8 = w3v1^2-w4v1^2+w3v2^2-w4v2^2+w3v3^2-w4v3^2 ; 
jade9 = w4v1^2-w5v1^2+w4v2^2-w5v2^2+w4v3^2-w5v3^2 ; 
A = multiAffineSpace(R)
F=polySystem {jade0, jade1, jade2, jade3, jade4, jade5, jade6, jade7, jade8, jade9 }; 
pt=point{flatten {{.878350739656677-.655267722694227*ii,
.992311073747755+.624621562952355*ii,
-.0756482712266543+.585129100769662*ii},{
-.0874625187231375+.801322947708998*ii,
-.424393626311276-.189829556344162*ii,
-1.04426756056119+.010032611254391*ii},{
1.17293648518138+.419135675010108*ii,
-.716316717577425+.399760065697915*ii,
.20570523824176-.997858436351401*ii},{
-.87323559289995-.558018942856899*ii,
.888959981317501-.54097142448466*ii,
.0108016860631553-.590653623384341*ii},{
.44614976208454-.746970959032991*ii,
.0865526973753582+.772655502291471*ii,
-1.26212825023999-.211060561983885*ii}}}
win=dimensionPolytope(pt,wCollection(A,F))

needsPackage"Polyhedra"
vertices win
latticePoints win 
restart
loadPackage"Bertini"
bertiniPosDimSolve({jade0, jade1, jade2, jade3, jade4, jade5, jade6, jade7, jade8, jade9 })

///

-- this file is loaded by NumericalAlgebraicGeometry
-- WSet-*.m2 are loaded by NAGtypes
-- !!! perhaps merge with witness-set.m2 ???

moveSlicingVariety(MultiAffineWSet,MultiSlicingVariety) := (W,S) -> (
    E := equations W#"equations";
    A := flatten entries map slicingVariety W;
    B := flatten entries map S;
    ptsB := track(E|A, E|B, points W, 
	NumericalAlgebraicGeometry$gamma=>exp(random(0.,2*pi)*ii));
    if any(ptsB, p->status p =!= Regular) then null -- failed
    else multiAffineWSet(W#"equations",S,ptsB) 
    )

TEST /// -- an example of moving a slice for a multiAffineWSet
restart 
debug needsPackage "NAGtypes"
debug needsPackage "WitnessCollections"
errorDepth = 0
A = multiAffineSpace(CC_53,{2,2},x)
use ring A 
-- multi=homogenized parabola y-z^2=0 where y=x_(0,1) and z=x_(1,1)
F = polySystem {x_(0,1)*x_(1,2)^2-x_(1,1)^2*x_(0,2),x_(1,2)-1,x_(0,2)-1}
S = multiSlicingVariety(A, {rationalMap matrix{{x_(0,1)-1}}, rationalMap map((ring A)^0,(ring A)^1,0)})
pts = {point{{1,1,1,1}},point{{1,1,-1,1}}}
W = multiAffineWSet(F,S,pts)
dim W
degree W
S = randomSlicingVariety(A,{1,0})
W' = moveSlicingVariety(W,S)
pts' = points W'
residuals = for p in pts' list 	   
{matrix evaluate(map slicingVariety W',p), evaluate(W'#"equations",p)}
assert all(flatten residuals, r->norm r < 1e-6)
///

end
check "WitnessCollections"