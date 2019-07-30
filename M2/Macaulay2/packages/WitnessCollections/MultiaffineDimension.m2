
newPackage(
    "MultiaffineDimension",
    Version => "1.0", 
    Date => "Feb 2019", 
    Authors => {
   {Name => "Jose Israel Rodriguez",
       Email => "Jose@Math.wisc.edu",
       HomePage => "http://www.math.wisc.edu/~jose/"}
    },
    Headline => "Determines the dimension of a multiaffine variety ",
    DebuggingMode => true,
    AuxiliaryFiles => true,
    PackageImports => {"NumericalAlgebraicGeometry","Polyhedra","NAGtypes"},
    PackageExports => {"NumericalAlgebraicGeometry","Polyhedra","NAGtypes"},
  Configuration => { "RandomCoefficients"=>CC },
  CacheExampleOutput => false
)

--path=prepend("/Users/jo/Documents/GoodGit/MultiaffineDimension",path)
--loadPackage("MultiaffineDimension",Reload=>true)
--examples MultiaffineDimension

randomCC=()->random CC
randomRR=()->((-1)^(random(1,2)) *random RR)
randomZZ=()->random(1,30103)
randomValue=(kk)-> if kk===CC then randomCC() else if kk===RR then randomRR() else randomZZ() 

export { 
    "multiaffineDimension",--documented
    "restrictionCodimension","RestrictionInformation",--documented
    "sequentialDimension",
    "wellPositionedLinearSpace",
    "wellPositionedIdeal",
    "LinearSpaceType"    
        }




--###################################
-- TYPE DEFINITIONS
--###################################
--DimensionReduction=new Type of MutableHashTable
--DimensionMatroidal=new Type of MutableHashTable

--###################################
-- METHODS
--###################################
randomVector=method(Options=>{		})
randomVector(ZZ,Thing):= o->(n,R) ->apply(n,i->randomValue(R))--list of length n of randomValue
-*
--testing random vector
randomVector(3,CC)
randomVector(3,RR)
randomVector(3,ZZ)
*-

--This function returns a polyhedron. 
----INPUTS 
----F
--list of polynomials system in R=CC[x1,x2,...,xn].
---- groupVars
--list the length gens R that is used to partition the variables. if position i has the same value as position j then xi and xj are in the same variable group. 
----pt
--witness point. 
 
multiaffineDimension=method(TypicalValue=>Thing,Options=>{})
multiaffineDimension(List,List,Point) := o ->(F,G,pt)->( --(polynomial system, k variable groups, a general witness point)
    R := ring first F;
    variableGroups := apply(G, ip -> matrix{ apply( ip, j -> (gens R)_j) } );
    thePartialJacs:= apply(variableGroups, vg -> --one variable group
	sub(diff (vg,transpose matrix{F}),matrix pt)
	);
    intrinsicDimension :=  thePartialJacs/numColumns//sum - numericalRank(transpose matrix{thePartialJacs},Threshold => .0001);
    theJacs := apply(drop(subsets thePartialJacs,1),---we drop the empty set. 
	Icom -> if #Icom>0 	then matrix{Icom} else error"#Icom=0");
    subsetsIV := drop(subsets(#G),1);--all nonempty subsets of [0,1,..,#G-1]
    M := {};
    v := {};
    scan(subsetsIV,Icom->(
    	M = append(M,apply(#G,i->if member(i,Icom) then 0 else 1 ));
	DFI :=  matrix{apply(Icom,i->thePartialJacs_i)};
	v = append(v,
	    DFI//numColumns - numericalRank(transpose DFI,Threshold => .0001)
	    )));
--Inequalities.  
    M = matrix M;
--  print theVarGroups;
    v = transpose matrix {v};
  print (M,v);
--Equalities
  N:=matrix {apply(#G,i->1)};
  w:=matrix{{intrinsicDimension}};
--  print (N,w);
  P:=polyhedronFromHData(-M,-v,N,w);
  return P
    );


-*
       R = CC[x1,x2,x3,y];
       F = {x1-x2-y};
       pt = point{{2_CC,1_CC,1_CC,3_CC}};
       G={{0,1,2},{3}}
       P = multiaffineDimension(F,G,pt)
       latticePoints P

       R = CC[x0,x1,y0,y1];
       F = {x0*y0+x1*y1+1};
       pt = point{{0_CC,0_CC,1_CC,-1_CC}};--Bad point because not on the variety
       pt = point{{0_CC,1_CC,0_CC,-1_CC}};
       pt = point{{2_CC,5_CC,1_CC,-3/5_CC}};
       G={{0,1},{2,3}}
       P = multiaffineDimension(F,G,pt)
       latticePoints P

    R = CC[x1,x2,y1,y2,z1,z2,z3]
    F={x1+x2-2,
	y1+2*y2-5,
	17*y1+y2+13*z1+17*z2-23}
    pt = point{{1_CC,1_CC,1_CC,2_CC,-1_CC,1_CC,3_CC}};
    G={{0,1},{2,3},{4,5,6}}
    P = multiaffineDimension(F,G,pt)
    latticePoints P    

*-


-- Simplified versions --
-*
multiaffineDimension(List,Sequence,Point) := o ->(F,numVarsInGroup,pt)->(
    indexGroups:=toList(0..(#numVarsInGroup-1));    
    groupVars:={};
    scan(#numVarsInGroup,i->scan(numVarsInGroup_i,j->groupVars=groupVars|{i}));
--    print indexGroups;
--    print groupVars;
    multiaffineDimension(F,indexGroups,groupVars,pt)) 
--
multiaffineDimension(List,Point) := o ->(F,pt)->(
    numGens:=#gens ring first F;
    multiaffineDimension(F,toList(0..(numGens-1)),toList(0..(numGens-1)),pt)) 
*-

RestrictionInformation=new Type of MutableHashTable
LinearSpaceType=new Type of List

restrictionCodimension=method(TypicalValue=>Thing,Options=>{
	})
restrictionCodimension(Polyhedron,List) := o ->(P,permutation)->(
    if isEmpty(P) then error"P is empty." ;
    k:=ambDim P;
    m:=new RestrictionInformation from {};
    if max(permutation)=!=k-1 or k=!=#permutation then error"permutation not recognized ";
--    print 1;
    A:=   permutation/(i->apply(k,j->if i===j then 1 else 0));
    newP:=P;
    scan(k,i->(
--	    print (A_i);
--	    r:=(k-i-1);
--	    print latticePoints newP;
	    newM:= latticePoints affineImage(matrix{A_i},newP);
--    	    print newM;
	    newM= newM/entries/flatten/flatten//flatten//max;
	    if newM>0 then newM=newM-1;
--	    print newM;
	    m#i=newM;	    
	    if newM>0 and not isEmpty(newP) 
	    then (
		print 2;
		Q:=polyhedronFromHData(-matrix{A_i},-matrix {{newM}});
  	    	newP=intersection(newP,Q);
  	    	print latticePoints newP;
	    	)));
    m#"Permutation"=permutation;
    bfe:=apply(k,i->m#i);    
    m#"CodimensionRestriction"=flatten entries ((matrix A)*matrix transpose {bfe});
    LT:={};
    cr:=m#"CodimensionRestriction";
    B:=entries diagonalMatrix(apply(#cr,i->1));
    scan(#B,i->scan(cr_i,j->LT=LT|{B_i}));
    m#"RestrictionType"=new LinearSpaceType from B;
    idMatrix:= diagonalMatrix(apply(k,i->1));
    newP=affineImage(idMatrix,newP,-(matrix A)*matrix transpose {bfe});
    m#"DimensionReduction"=newP;
    return m)

restrictionCodimension(Polyhedron) := o ->(P)->restrictionCodimension(P,toList(0..(ambDim P-1)))
restrictionCodimension(Polyhedron,Nothing) := o ->(P,n)->restrictionCodimension(P,reverse toList(0..(ambDim P-1)))
    
 ------
sequentialDimension=method(TypicalValue=>Thing,Options=>{
	})
sequentialDimension(Polyhedron,List) := o ->(P,permutation)->(
    if isEmpty(P) then error"P is empty." ;
    k:=ambDim P;
    if max(permutation)=!=k-1 or k=!=#permutation then error"permutation not recognized ";
    A:=  apply(permutation,i->apply(k,j->if permutation_i>=permutation_j then 1 else 0));
--    print A;
    newP:=P;
    dimLowerBound:=1;
    linType:={};
    scan(k,i->(
	    newM:= latticePoints affineImage(matrix{A_i},newP);
	    newM= newM/entries/flatten/flatten//flatten//max;
--    	    print newM;
    	    if newM>dimLowerBound then (
		d:=(newM-dimLowerBound);
		linType=linType|apply(d,j->A_i));
	    dimLowerBound=newM));	    
    return linType)

sequentialDimension(Polyhedron) := o ->(P)->sequentialDimension(P,toList(0..(ambDim P-1)))
sequentialDimension(Polyhedron,Nothing) := o ->(P,n)->sequentialDimension(P,reverse toList(0..(ambDim P-1)))

----------------------------------------------------------------------------------------------------------------

wellPositionedLinearSpace=method(TypicalValue=>Thing,Options=>{
	})
wellPositionedLinearSpace(Polyhedron,List,List) := o ->(P,r,s)->(
--    r:=toList (0,..,(ambDim P)-1);
--    s:=toList (0,..,(ambDim P)-1);
--
    ri:=restrictionCodimension(P,r);
    LT:=new LinearSpaceType from (ri#"RestrictionType");
    Q:=ri#"DimensionReduction";
    LT=LT|sequentialDimension(Q,s);
    return LT
    )
wellPositionedLinearSpace(Polyhedron) := o ->(P)->(
    r:=s:=toList (0..(ambDim P)-1);
    wellPositionedLinearSpace(P,r,s)
    )  
wellPositionedLinearSpace(Polyhedron,Nothing) := o ->(P,n)->(
    r:=s:=reverse toList (0..(ambDim P)-1);
    wellPositionedLinearSpace(P,r,s)
    )  


wellPositionedIdeal=method(TypicalValue=>Thing,Options=>{
	})
wellPositionedIdeal(Ring,List,List,LinearSpaceType) := o ->(R,E,L,LT)->(
    indexPositions := apply(E,i->positions(L,g->i===g));
--    print indexPositions;
    vg := apply(indexPositions,ip->apply(ip,j->(gens R)_j));
--    print vg;
--    print LT;
    return apply(LT,T->(aPoly:=0; 
	    scan(#T,i->if T_i==1 
		then aPoly=(aPoly+randomCombo(R,vg_i,random(coefficientRing R))));
	    print aPoly;
	    return ideal aPoly
	    )))

wellPositionedIdeal(Ring,List,LinearSpaceType) := o ->(R,groupVars,LT)->wellPositionedIdeal(R,toList set groupVars,groupVars,LT) 
--
wellPositionedIdeal(Ring,Sequence,LinearSpaceType) := o ->(R,numVarsInGroup,LT)->(
    indexGroups:=toList(0..(#numVarsInGroup-1));    
    groupVars:={};
    scan(#numVarsInGroup,i->scan(numVarsInGroup_i,j->groupVars=groupVars|{i}));
    wellPositionedIdeal(R,indexGroups,groupVars,LT)) 
--
wellPositionedIdeal(Ring,LinearSpaceType) := o ->(R,LT)->(
    numGens:=#gens R;
    wellPositionedIdeal(R,toList(0..(numGens-1)),toList(0..(numGens-1)),LT)) 



--##########################################################################--
-- INTERNAL METHODS
--##########################################################################--
----------------------------------------
--newHyperplanes=A->for i to (numColumns A)+1 list randomVector(numRows A)
makeJac=(system,unknowns)->(--it is a list of lists of partial derivatives of a polynomial
         for i in system list for j in unknowns list  diff(j,i))

randomCombo=(R,vg,c)->(
	theLinear:=c;
	scan(vg,x->theLinear=theLinear+random(coefficientRing R)*x);
	return theLinear)


--beginDocumentation()

--load "./DOC_MD.m2";

TEST///
--load concatenate(MultiprojectiveWitnessSets#"source directory","./AEO/TST/Example1.tst.m2")
///


end
restart
loadPackage"Bertini"
makeB'InputFile(storeBM2Files,
    AffVariableGroup=>{x,y},
    B'Polynomials=>{x,y}
    )
runBertini(storeBM2Files)
readFile(storeBM2Files)
R=QQ[x,y]
bertiniPosDimSolve {x,y}


bertiniZeroDimSolve

parString=(aString)->("("|toString(aString)|")");
addSlash=(aString)->(
    if aString_-1===" " then error (aString|" cannot end with whitespace.");
    if aString_-1=!="/" then aString=aString|"/";
    return aString    )
projectDimPolytope=(P,indexGroups,downstairDim)->(
    v:=apply(#indexGroups-1,i->1)|{0};
    Q:=polyhedronFromHData(matrix{v,-v},transpose matrix{{downstairDim,-downstairDim}});
    return intersection(P,Q))
