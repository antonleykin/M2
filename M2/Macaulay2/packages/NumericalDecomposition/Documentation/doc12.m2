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
      a  (system need not be square)---TODO
    G:List
      a list of lists of integers such that j is in the i-th list of G iff (vars F)_(0,j) is in the i-th variable group
    pt:Point
      a general point of an irreducible component of V(F)
  Outputs
    P:Polyhedron
      the multidimension of the irreducible component of V(F) containing pt in terms of a polymatroid polytope, use latticePoints P to recover the multidimension explicitly.
  Description
    Text
      multiaffineDimension(F,G,pt)
    Example
      R = CC[x1,x2,x3,y];
      F = gateSystem polySystem {x1-x2-y};
      pt = point{{2,1,1,3_CC}};
      G={{0,1,2},{3}}
      P = multiaffineDimension(F,G,pt)
      latticePoints P
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
    Example
      R = CC[x1,x2,x3,y];
      F = gateSystem polySystem {x1-x2-y};
      pt = point{{2,1,1,3_CC}};
      G={{0,1,2},{3}}
      P = multiaffineDimension(F,G,pt)
      getSequenceSC P
  Caveat
    We assume that the point pt is generic and that the numerical rank is computed correctly.
///
