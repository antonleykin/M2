-- Jose: steps 1 and 2
doc ///--multiaffineDimension
  Key
    multiaffineDimension
    (multiaffineDimension,GateSystem, VariableGroup, Point)--(F,G,pt)
  Headline
    a method to determine the multiaffine (dimension polytope)
  Usage
    P = multiaffineDimension(F,G,pt)
  Inputs
    F:GateSystem
      a gate system that defines a variety
    G:List
      a list of lists of integers such that j is in the i-th list of G iff (vars F)_(0,j) is in the i-th variable group
    pt:Point
      a general point of an irreducible component of V(F)
  Outputs
    P:Polyhedron
      the multidimension of the irreducible component of V(F) containing pt
      in terms of a polymatroid polytope,
      use latticePoints P to recover the multidimension explicitly
  Description
    Text
      This method computes the multidimension by
      computing the ranks of submatrices of the Jacobian of F.
      When using exact artithmetic the method rank is used, otherwise numericalRank is used.
      Each submatrix and rank induces an inequality that defines the polytope P
      according to "A numerical toolkit for multiprojective varieties, Algorithm 2.3".
    Example
      V = matrix{declareVariable\{x1,x2,x3,y}}
      F = gateSystem(V, gateMatrix{{x1-x2-y}})
      pt = point{{2,1,1,3_CC}};
      G={{0,1,2},{3}}
      P = multiaffineDimension(F,G,pt)
      latticePoints P
      --Finite field arithmetic also works.
      V2 = matrix{declareVariable\{a1,a2,a3,b}}
      F = gateSystem(V2, gateMatrix{{a1-a2-b}})
      pt = point{{2,1,1,3_ZZ/30103}};
      G={{0,1,2},{3}}
      P = multiaffineDimension(F,G,pt)
      latticePoints P
    Text
      It is possible to go from a list of ring elements to a gate system.
    Example
      TODO
  Caveat
    We assume that the point pt is generic and that the numerical rank is computed correctly.
///

doc ///--getSequenceSC
  Key
    getSequenceSC
    (getSequenceSC,Polyhedron)--P
  Headline
    a method to compute a slicing coarsening sequence
  Usage
    getSequenceSC P
  Inputs
    P:Polyhedron
      a polyhedron representing a multiaffine dimension
  Outputs
    SCS:Sequence
      a sequence of integers and pairs of integers
  Description
    Text
      The output SCS is a sequence of integers and
      pairs of integers where an integer i in the sequence
      represents slicing in the the i-th variable group
      and a pair (i,j) represents combining the i-th and j-th variable group.
      See "A numerical toolkit for multiprojective varieties, Section 4" for more details.
    Example
      R = CC[x1,x2,x3,y];
      F = gateSystem polySystem {x1-x2-y};
      pt = point{{2,1,1,3_CC}};
      G={{0,1,2},{3}}
      P = multiaffineDimension(F,G,pt)
      getSequenceSC P
///
