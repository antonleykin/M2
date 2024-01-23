--- status: Draft
--- author(s): Gregory G. Smith
--- notes: 

-- in Classic: (ideal, String)

document { 
     Key => ideal,
     Headline => "make an ideal",
     }
document { 
     Key => {(ideal,List), (ideal,Sequence)},
     Usage => "ideal L",
     Inputs => {
	  "L" => {"or a ", TO2("Sequence","sequence"), " of ", TO2("RingElement", "ring elements")}
	  },
     Outputs => {
	  LeftIdeal => {"which is generated by the ", TO2("List","list"), " or ",
	  TO2("Sequence","sequence"), " of ",  TO2("RingElement", "ring elements")}
	  },
     PARA {
	 "If the input (generators) is in a commutative or skew-commutative ring the type of the output is a (two-sided) ", TO "Ideal","." 
	 },
     EXAMPLE {
	  "R = ZZ/101[w,x,y,z];",
	  "ideal{x^2-w*y, x*y-w*z, x*z-y^2}",
 	  "ideal(y^2-x*z,x^2*y-z^2,x^3-y*z)",
	  "E = ZZ/2[x,y, SkewCommutative => true];",
	  "ideal(x^2,x*y)",
	  "W = QQ[x,dx, WeylAlgebra => {x => dx}];",
	  "ideal(dx*x+x*dx)",
	  "I = ideal(12,18)",
	  "mingens I"
	  },
     PARA {
	  "An empty list or sequence of generators will yield an ideal of ", TO "ZZ", ", which
	  can be promoted to another ring, if desired."
	  },
     EXAMPLE lines ///
     ideal ()
     promote(oo,R)
     ///,
     SeeAlso => {LeftIdeal, Ideal, PolynomialRing}
     }
document { 
     Key => (ideal,Matrix),
     Usage => "ideal M",
     Inputs => {
	  "M" => {"whose ", TO "entries", " are ", TO2("RingElement", "ring elements")}
	  },
     Outputs => {
	  LeftIdeal => {"which is generated by all the entries in ", TT "M"}
	  },
     PARA {
	 "If the input (generators) is in a commutative or skew-commutative ring the type of the output is a (two-sided) ", TO "Ideal","." 
	 },
     EXAMPLE {
	  "R = ZZ/7[w,x,y,z];",
	  "f = vars R",
	  "ideal f",
	  "g = matrix{{x^2-w*y, x*y-w*z, x*z-y^2},{y^2-x*z,x^2*y-z^2,x^3-y*z}}",
	  "ideal g"
	  },
     PARA{},     
     "The functions ", TO "minors", ", ", TO "pfaffians", " and ", 
     TO2("fittingIdeal", "fitting ideal"),  " are typically more useful methods of creating an ", 
     TO2("Ideal","ideal"), " from a matrix.",
     SeeAlso => {entries, flatten}
     }
document { 
     Key => {(ideal,RingElement), (ideal,Number)},
     Usage => "ideal r",
     Inputs => {
	  "r"
	  },
     Outputs => {
	  LeftIdeal => {"which is generated by the ", TO2("RingElement", "ring element"), " ", TT "r"}
	  },
      PARA {
	 "If the input (generator) is in a commutative or skew-commutative ring the type of the output is a (two-sided) ", TO "Ideal","." 
	 },
     EXAMPLE {
	  "ideal 7",
	  "R = ZZ/2[x,y];",
	  "f = x^2+y^2;",
	  "ideal f"
	  },
     PARA {
	  "The zero ideal with one generator can be made:"
	  },
     EXAMPLE "ideal 0_R",
     PARA {
	  "Alternatively, the zero ideal with no generators can be made in one of these ways:"
	  },
     EXAMPLE lines ///
     ideal id_(R^0)
     promote(ideal(), R)
     ///,
     PARA {
	  "See ", TO promote, ", ", TO (id,Module), ", and ", TO (ideal,Sequence), ".",
	  }
     }
document { 
     Key => {(ideal,Ring), (ideal,QuotientRing)},
     Headline => "returns the defining ideal",     
     Usage => "ideal R",
     Inputs => {
	  "R"
	  },
     Outputs => {
	  Ideal => {"which is the defining ideal of ", TT "R"}
	  },
     PARA {
	 "If the input (generators) is in a commutative or skew-commutative ring the type of the output is a (two-sided) ", TO "Ideal","." 
	 },
     PARA {
     	 "A ", TO2("QuotientRing","quotient ring"), " is a the quotient of its ", TO "ambient", 
     	 " ", TO2("Ring","ring"), " by its defining ideal.  Other rings have no ambient ring,
     	 and the defining ideal is its zero ideal."
	 }, 
     EXAMPLE {
	  "S = ZZ/2[x,y,z];",
	  "ideal S",
	  "R = S/(y^2-x*z,x^2*y-z^2)",
	  "ideal R",
	  "T = R/(x^3-y*z)",
	  "ideal T",
	  "ambient T",
	  "sing = singularLocus T",
	  "ideal sing",
	  "ambient sing"
	  },
     SeeAlso => {ambient, singularLocus}
     }
document { 
     Key => (ideal,MonomialIdeal),
     Headline => "converts a monomial ideal to an ideal",     
     Usage => "ideal I",
     Inputs => {
	  "I"
	  },
     Outputs => {
	  Ideal => {"which is generated by the monomials in ", TT "I"} 
	  },
     EXAMPLE {
	  "R = QQ[x,y,z];",
	  "I = monomialIdeal(x*y^2, x^2*z, y^2*z)",
     	  "ideal I"
	  },
     PARA{},
     "Most operations between ", TO2("Ideal","ideals"), " and ", TO2("MonomialIdeal","monomial ideals"), 
     " automatically perform the necessary conversions, so one rarely needs to apply the function ", 
     TT "ideal", " to a ", TO2("MonomialIdeal","monomial ideal"), ".",
     EXAMPLE {
	  "I * ideal I",
	  "I + ideal(x*y+y*z)"
	  },
     SeeAlso => {monomialIdeal, (generators,LeftIdeal)}
     }
document { 
     Key => (ideal,Module),
     Headline => "converts a module to an ideal",          
     Usage => "ideal M",
     Inputs => {
	  "M" => "which is a submodule of a free module of rank 1"
	  },
     Outputs => {
	  LeftIdeal => {"given by the generators of ", TT "M"}
	  },
     EXAMPLE {
	  "R = QQ[w,x,y,z];",
	  "f = map(R^1,R^3, matrix{{x^2-w*y, x*y-w*z, x*z-y^2}})",
	  "image f",
	  "ideal image f"
	  }, 
     SeeAlso => {(generators,Module), ambient, isSubmodule, isIdeal}
     }
