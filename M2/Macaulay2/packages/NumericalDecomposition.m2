-- -*- coding: utf-8 -*-
-- licensed under GPL v2 or any later version

newPackage (
     "NumericalDecomposition",
     Version => "1.22",
     Date => "Oct 2023",
     Headline => "numerical methods for decomposition of varieties",
     --HomePage => "",
     AuxiliaryFiles => true,
     Authors => {
	  {Name => "Anton Leykin", Email => "leykin@math.gatech.edu", 
	      HomePage => "https://people.math.gatech.edu/~aleykin3"},
	  {Name => "Tim Duff"},
	  {Name => "Jose I. Rodriguez"}
	  },
     Keywords => {"Numerical Algebraic Geometry"},
     PackageExports => {"NAGtypes",
	 "NumericalLinearAlgebra",
	 "SLPexpressions",
	 "NumericalAlgebraicGeometry",
	 },
     DebuggingMode => true
     )
-- Any symbols or functions that the user is to have access to
-- must be placed in one of the following two lists
export {"uGeneration"}
load "./NumericalDecomposition/u-generation.m2"


