-- Copyright 1999-2002 by Anton Leykin and Harrison Tsai

local localCohomUli
local localCohomOT
local localCohomILOaku
local computeLocalCohomOT
local localCohomRegular
local preimage
local localCohomILOTW
local Generator
------------------------------------------------------------------
-- LOCAL COHOMOLOGY
-- 
-- Caveats: does not simplify the presentations of 
--          the cohomology modules returned
-------------------------------------------------------------------
localCohom = method(Options => {Strategy => Walther, LocStrategy => null})


----------------------------------------------------------------------------------------
-- computes the local cohomology H_I(R), where I is an ideal in a polynomial ring R.
-- using Walther's algorithm
--
-- option: Strategy (sets the way localizations are computed)
----------------------------------------------------------------------------------------

localCohom(      Ideal) := HashTable => o ->    I  -> localCohom(toList (0..numgens I), I, o)
localCohom(ZZ,   Ideal) := HashTable => o -> (n,I) -> localCohom({n}, I, o)
localCohom(List, Ideal) := HashTable => o -> (l,I) -> (
     -- Promote I to the Weyl algebra if it is not already there
     if #(ring I).monoid.Options.WeylAlgebra == 0
     then I = sub(I, makeWA ring I); -- TODO: what degrees are best for the differentials?
     if (o.Strategy == Walther and o.LocStrategy === null)
     then localCohomUli (l,I)
     else (
	  R := ring I;
	  createDpairs R;
	  localCohom(l, I, R^1 / ideal R.dpairVars#1, o)
	  )
     );

localCohomUli = (l, I) -> (
     -- error checking to be added
     -- I is assumed to be an ideal in WA generated by polynomials
     f := first entries gens I;
     r := #f;
     W := ring I;
     subISets := select(subsets toList (0..r-1), s -> s =!= {});
                  
     -- Step1.
     -- Calculate J^/delta( (F_/theta)^s ) and b^/delta_(F_/theta)(s)   for all /theta
     pInfo(1, "localCohom: Computing b-functions and annihilators...");
     J := new MutableHashTable;
     bF := new MutableHashTable;
     scan(subISets, theta->(
	       Ftheta := product(theta, i->f#i);
	       J#theta = AnnFs(Ftheta);
	       bF#theta = globalBFunction(Ftheta);
	       pInfo(666, {"theta = ", theta, " F = ", Ftheta});
	       pInfo(666, {"AnnFs = ", J#theta, " bF = ", bF#theta}); 
	       ));
     
     -- Step 2.
     -- a = min integer root of all bF-s
     a := min flatten (subISets / (theta -> getIntRoots(bF#theta)));
     if a == infinity then a = 0;
     pInfo(666, {"BEST POWER = " , a});
     -- Substitute s = a in J-s
     scan(subISets, theta -> (
	       AforS := map(W, ring J#theta, vars W | matrix {{a_W}});
	       J#theta = AforS J#theta;
	       ));
     
     -- Step 3.
     -- Compute the Chech complex 
     C := new MutableList from toList ((r+1):());
     M := new MutableList from toList (r:());
     pInfo(1, "Constructing Cech complex...");
     C#0 = directSum { {} => W^1 / ideal W.dpairVars#1 };
     apply(toList(1..r), k->(
	       C#k = directSum (select(subISets, u -> #u == k) / (theta -> (
			      theta => W^1 / J#theta 
			      )));
	       
	       M#(k-1) = map (C#k, C#(k-1), (i,j) -> (
			 i0 := (indices C#(k-1))_j;
			 j0 := (indices C#k)_i;
			 if isSubset(i0, j0) 
			 then (
			      l := first toList (set j0 - set i0);
			      (-1)^(position(j0, u -> u == l)) * (f_l)^(-a)
			      )  
			 else 0_W
			 )
		    );   
	       ));     
     -- Step 4.
     -- Compute the homology of the complex
     pInfo(1, "Computing cohomology...");	  
     ret := new HashTable from (
     	  (if member(0, l) then {0 => kernel M#0} else {})
     	  | 
     	  ( delete(r,delete(0,l)) / (k -> k => homology(M#k, M#(k-1))) )  
     	  |
     	  (if member(r, l) then {r => cokernel M#(r-1)} else {})
     	  );
     ret
     );

----------------------------------------------------------------------------------------
-- computes the local cohomology H_I(D), where I is an ideal in a polynomial ring
--                                             D is a cyclic D-module
-- two versions: 
-- 1) returns cohomology in every degree  
-- 2) returns cohomology in degrees in the list passed as an argument
--
-- option:  Strategy (sets the way localizations are computed)
----------------------------------------------------------------------

localCohom(      Ideal, Module) := HashTable => o ->   (I,M) -> localCohom(toList (0..numgens I), I, M, o)
localCohom(ZZ,   Ideal, Module) := HashTable => o -> (n,I,M) -> localCohom({n}, I, M, o)
localCohom(List, Ideal, Module) := HashTable => o -> (l,I,M) -> (
     -- Promote I to the Weyl algebra if it is not already there
     if #(ring I).monoid.Options.WeylAlgebra == 0
     then (
	 D := makeWA ring I; -- TODO: what degrees are best for the differentials?
	 I = sub(I, D);
	 M = M ** D;
	 );
     pInfo (1, "localCohom: holonomicity check ...");
     if not isHolonomic M then
     error "expected a holonomic module";
     if o.Strategy == Walther then (
     	  if o.LocStrategy === null then localCohomRegular(l,I,M)
     	  else if o.LocStrategy == OaTaWa then localCohomILOTW(l,I,M)
     	  else if o.LocStrategy == Oaku then localCohomILOaku(l,I,M)
     	  )
     else if o.Strategy == OaTa then localCohomOT(l,I,M)
     else error "unknown option"
     );

----------------------------------------------------------
-- iterated localizations + localize by Oaku 
----------------------------------------------------------
localCohomILOaku = method()
localCohomILOaku(List, Ideal, Module) := (l, I, M) -> (
     -- error checking to be added
     -- I is assumed to be an ideal in WA generated by polynomials
     f := first entries gens I;
     FT := theta -> product(theta, i->f#i);
     r := #f;
     W := ring I;
     subISets := select(subsets toList (0..r-1), s -> s =!= {});
     
     L := new MutableHashTable;
     L#{} = new HashTable from {LocModule => M, Generator => 1_(ring M)}; 
     
     pInfo(1, "Constructing Cech complex...");
     C := new MutableList from toList ((r+1):());
     MM := new MutableList from toList (r:());
     C#0 = directSum { {} => M };
     -- (For this strategy only!) 
     -- keep track of powers of f_i in the localizations
     FPower := new MutableList from toList (r:0);
     scan(toList(1..r), k->(
	  dsArgs := select(subISets, u -> #u == k) / (theta -> (
		    theta' := if k == 1 then {}
		    else first select(1, subISets, u -> #u == k-1 
			 and isSubset(u,theta));
		    i := first toList (set theta - set theta');
		    papa := L#theta'#LocModule;
		    pInfo(666, {"iterated localization: ", theta', " => ", theta});
		    locPapa := computeLocalization(papa, f_i, 
			 {GeneratorPower, annFS}, 
			 new OptionTable from {Strategy => Oaku});
		    pInfo(666, {"Gen power = ", locPapa.GeneratorPower, 
			      " annFS = ", locPapa.annFS});
		    if locPapa.GeneratorPower < FPower#i 
		    then FPower#i = locPapa.GeneratorPower;
		    -- compute the locModule
		    I := locPapa.annFS;
		    subMap := map(W, ring I, vars W | matrix {{(FPower#i)_W}});
		    locIdeal := ideal subMap gens I; 
		    L#theta = new HashTable from {
			 LocModule => W^1/locIdeal, 
			 Generator => L#theta'#Generator * (f_i)^(-FPower#i)
			 };
		    theta => L#theta#LocModule
		    ));
	  C#k = directSum dsArgs;
	            	       	       
	  TempM := map (C#k, C#(k-1), 0);
	  scan(indices C#(k-1), i0->(
		    flagOK := true; 
		    scan(indices C#k, j0->(
			      if flagOK and isSubset(i0, j0)  
			      then (
				   l := first toList (set j0 - set i0);
				   gi := ((L#i0)#Generator);
				   gj := ((L#j0)#Generator);
				   if gj % gi !=0 
				   then ( 
					-- Have to recompute the previous component
					flagOK = false;
					error "Bad luck..."
					--!!! Write it sometime
					);
				   TempM = TempM 
				   + (C#k)_[j0] -- injection from j0-th component
				   * map(L#j0#LocModule, L#i0#LocModule, 
					(-1)^(position(j0, u -> u == l)) 
					* (gj//gi) -- (-1)^(...) id
					)
				   * (C#(k-1))^[i0]; -- projection onto i0-th component 
				   pInfo(666, {"multiplier: ", gj//gi}); 
				   )	   
			      )))); 
	  MM#(k-1) = TempM;
	  ));     
     -- Step 4.
     -- Compute the homology of the complex
     pInfo(1, "Computing cohomology...");	  
     ret := new HashTable from (
     	  (if member(0, l) then {0 => kernel MM#0} else {})
     	  | 
     	  ( delete(r,delete(0,l)) / (k -> k => homology(MM#k, MM#(k-1))) )  
     	  |
     	  (if member(r, l) then {r => cokernel MM#(r-1)} else {})
     	  );
     ret
     );

----------------------------------------------------------
-- iterated localizations + localize by OTW 
----------------------------------------------------------
localCohomILOTW = method()
localCohomILOTW(List, Ideal, Module) := (l, I, M) -> (
     -- error checking to be added
     -- I is assumed to be an ideal in WA generated by polynomials
     f := first entries gens I;
     FT := theta -> product(theta, i->f#i);
     r := #f;
     W := ring I;
     subISets := select(subsets toList (0..r-1), s -> s =!= {});
     
     L := new MutableHashTable;
     L#{} = new HashTable from {LocModule => M, Generator => 1_(ring M)}; 
     
     pInfo(1, "Constructing Cech complex...");
     C := new MutableList from toList ((r+1):());
     MM := new MutableList from toList (r:());
     C#0 = directSum { {} => M };
     scan(toList(1..r), k->(
          dsArgs := select(subISets, u -> #u == k) / (theta -> (
		    theta' := if k == 1 then {}
		    else first select(1, subISets, u -> #u == k-1 
			 and isSubset(u,theta));
		    i := first toList (set theta - set theta');
		    papa := L#theta'#LocModule;
		    pInfo(666, {"iterated localization: ", theta', " => ", theta});
		    locPapa := computeLocalization(papa, f_i, 
			 {LocModule, GeneratorPower}, 
			 new OptionTable from {Strategy =>OTW});
		    L#theta = new HashTable from {
			 LocModule => locPapa#LocModule, 
			 Generator => L#theta'#Generator * 
			 (f_i)^(-locPapa.GeneratorPower)
			 };
		    theta => L#theta#LocModule
		    ));
	  C#k = directSum dsArgs;
	       	       	       
	  TempM := map (C#k, C#(k-1), 0);
	  scan(indices C#(k-1), i0->(
		    scan(indices C#k, j0->(
			      if isSubset(i0, j0) 
			      then (
				   l := first toList (set j0 - set i0);
				   gi := ((L#i0)#Generator);
				   gj := ((L#j0)#Generator);
				   if gj % gi !=0 
				   then error ("Bad luck: " | toString gj | 
					" is not divisible by " | toString gi);
				   TempM = TempM 
				   + (C#k)_[j0] -- injection from j0-th component
				   * map(L#j0#LocModule, L#i0#LocModule, 
					(-1)^(position(j0, u -> u == l)) * (gj//gi) 
					-- (-1)^(...) id
					)
				   * (C#(k-1))^[i0]; -- projection onto i0-th component 
				   pInfo(666, {"multiplier: ", gj//gi}); 
				   )	   
			      )))); 
	  MM#(k-1) = TempM;
	  ));     
     -- Step 4.
     -- Compute the homology of the complex
     pInfo(1, "Computing cohomology...");	  
     ret := new HashTable from (
     	  (if member(0, l) then {0 => kernel MM#0} else {})
     	  | 
     	  ( delete(r,delete(0,l)) / (k -> k => homology(MM#k, MM#(k-1))) )  
     	  |
     	  (if member(r, l) then {r => cokernel MM#(r-1)} else {})
     	  );
     ret
     );


----------------------------------------------------------------------------------------
-- computes the local cohomology H_I(D), where I is an ideal in a polynomial ring
--                                             D is a cyclic D-module
-- caveat: not smart, does not iterate localizations. 
----------------------------------------------------------------------------------------
localCohomRegular = method()
localCohomRegular(List,Ideal,Module) := (l, I, M) -> (
     -- error checking to be added
     -- I is assumed to be an ideal in WA generated by polynomials
     f := first entries gens I;
     FT := theta -> product(theta, i->f#i);
     r := #f;
     W := ring I;
     subISets := select(subsets toList (0..r-1), s -> s =!= {});
                  
     L := new MutableHashTable;
     bF := new MutableHashTable;
     scan(subISets, theta->(
	       Ftheta := FT(theta);
	       LOC := DlocalizationAll(M, Ftheta);
	       L#theta = new HashTable from {
		    LocModule => LOC.LocModule, 
		    Generator => Power(Ftheta, LOC.GeneratorPower)
		    };	       
	       pInfo(666, {"localization: ", theta, " => ", L#theta});
	       ));
     L#{} = new HashTable from {LocModule => M, Generator => Power(1_(ring M), 0)};
     
     -- Compute the Chech complex 
     pInfo(1, "Constructing Cech complex...");
     C := new MutableList from toList ((r+1):());
     MM := new MutableList from toList (r:());
     C#0 = directSum { {} => M };
     scan(toList(1..r), k->(
	       dsArgs := select(subISets, u -> #u == k) / (theta -> 
		    theta => L#theta#LocModule
		    );
	       C#k = directSum dsArgs;
	       TempM := map (C#k, C#(k-1), 0);
	       scan(indices C#(k-1), i0->(
		    	 scan(indices C#k, j0->(
				   if isSubset(i0, j0) 
			 	   then (
					l := first toList (set j0 - set i0);
					gi := (L#i0)#Generator;
					gj := (L#j0)#Generator;
					if gj#1 > gi#1
					then error "Bad luck!"; 
					-- have to fix that: go back and recalculate 
					-- the localizations
				        TempM = TempM 
					+ (C#k)_[j0] 
					* map(L#j0#LocModule, L#i0#LocModule, 
					     (-1)^(position(j0, u -> u == l)) 
					     * 
					     (gi#0)^(gi#1-gj#1) * (f#l)^(-gj#1))
					* (C#(k-1))^[i0]; 
					)	   
				   )))); 
	       MM#(k-1) = TempM;
	       ));     
     -- Step 4.
     -- Compute the homology of the complex
     pInfo(1, "Computing cohomology...");	  
     ret := new HashTable from (
     	  (if member(0, l) then {0 => kernel MM#0} else {})
     	  | 
     	  ( delete(r,delete(0,l)) / (k -> k => homology(MM#k, MM#(k-1))) )  
     	  |
     	  (if member(r, l) then {r => cokernel MM#(r-1)} else {})
     	  );
     ret
     );
---------------------------------------------------------------------------------




-------------------------------------------------------------------------------
-- computes the preimage of a submodule M of the target of f 
-- (more precisely, M and <target f> should have the same ambient module)
------------------------------------------------------------------------------
preimage = method();
preimage (Module, Matrix) := (M, f) -> (
     T := target f;
     g := map(T/M, T);
     kernel (g*f)
     );

---------------------------------------------------------------------------
-- prunes every element of the local cohomology hashtable 
---------------------------------------------------------------------------
pruneLocalCohom = method()
pruneLocalCohom(HashTable) := HashTable => h -> (
     new HashTable from apply(keys h, i->(
     	       i => Dprune relations prune h#i
     	       ))
     )

--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
-- computes local cohomology modules using algorithm of Oaku-Takayama
-- for a holonomic D-module
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------

localCohomOT = method()
localCohomOT(Ideal, Ideal) := (I, J) -> (
     if not J.?quotient then J.quotient = (ring J)^1/J;
     localCohomOT(I, J.quotient)
     )
localCohomOT(Ideal, Module) := (I, M) -> (
     computeLocalCohomOT(I, M, 0, numgens I)
     )

localCohomOT(List, Ideal, Module) := (l, I, M) -> (
     locOut := computeLocalCohomOT(I, M, min l, max l);
     locOut = hashTable apply(keys locOut, 
	  i -> if member(i, l) then i => locOut#i);
     locOut)


computeLocalCohomOT = (I, M, n0, n1) -> (
     -- error checking to be added
     -- 1. make sure I is contained in coordinate part of Weyl alg
          
     -- preprocessing
     m := gens I;
     r := rank source gens M;
     d := numgens source m;
     W := ring M;
     createDpairs(W);
     nW := numgens W;
     n := #W.dpairVars#0;
     N := presentation M;
     -- create the auxiliary D_(n+d) ring
     t := symbol t;
     Dt := symbol Dt;
     LCW := (coefficientRing W)(monoid [(entries vars W)#0,
	  t_0 .. t_(d-1), Dt_0 .. Dt_(d-1),
	  WeylAlgebra => join(W.monoid.Options.WeylAlgebra,
	       apply(toList(0..d-1), i->(t_i=>Dt_i)) )]);
     scan(d, i -> (t_i = LCW_(t_i); Dt_i = LCW_(Dt_i)));
     nLCW := numgens LCW;
     WtoLCW := map(LCW, W, (vars LCW)_{0..nW-1});
     LCWtoW := map(W, LCW, (vars W) | matrix{toList(2*d:0_W)});
     -- weight vector for restriction to t_1 = ... = t_d = 0
     w := join( toList(n:0), toList(d:1) );
     -- create KN such that (D_{n+d}^r/KN) \cong 
     -- ( R_f[s_1..s_d]f_1^{s_1}...f_d^{s_d} \os D_n^r/N ) ??
     F := LCW^r;
     Lm := WtoLCW m;
     twistN := {};
     i := 0;
     while (i < d) do (
	  j := 0;
	  while (j < numgens F) do (
	       twistN = append( twistN, (t_i - Lm_(0,i))*(gens F)_j );
	       j = j+1; );
	  i = i+1; );
     LN1 := transpose matrix apply(twistN, i -> entries i);	  
     -- create the twistings that will be applied to N
     twistList := apply( toList(0..nLCW-1), 
	  i -> LCW_i + sum(d, j -> (LCW_i * Lm_(0,j) - 
		    Lm_(0,j) * LCW_i) * Dt_j) );
     twistMap :=  map(LCW, LCW, matrix{twistList});
     -- twist generators of N into generators of KN;
     LN2 := twistMap(WtoLCW N);
     KN := LN1 | LN2;
     KN = map(LCW^(numgens target KN), LCW^(numgens source KN), KN);
     restrictOut := computeRestriction(cokernel KN, w, n0-1, n1+1, 
	  {HomologyModules, ResToOrigRing}, hashTable{Strategy => Schreyer});
     
     -- stash the homology groups

     locOut := hashTable apply(toList(0..d), i -> (-i+d) => 
	  LCWtoW ** (restrictOut#ResToOrigRing ** restrictOut#HomologyModules#i));
     
     locOut
     )



