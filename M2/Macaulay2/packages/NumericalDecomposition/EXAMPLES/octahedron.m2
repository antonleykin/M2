debug needsPackage "NumericalDecomposition"
declareVariable \ {x,y,z,w}
f = 1 + 2*x + 3*y^2+4*z^3+5*w^4
h = 1 + 2*x + 3*y+ 5*z + 7*w + 11*z*y + 13 * x* z + 17*x*w + 19*y*z+23*y*w +29*z*w+31*x*y*z+37*x*y*w+41*x*z*w+43*y*z*w+47*x*y*z*w
F = gateSystem(matrix{{x,y,z,w}},transpose gateMatrix{{f,h}})
G = new VariableGroup from {{0},{1},{2},{3}}

-- one point on the variety
pt = point{{4,-2/3,-1.44419, .765304}}
assert(norm evaluate(F,pt) < 0.001)

-- Dimension polytope
P = multiaffineDimension(F,G,pt)
latticePoints P
assert (
    P == convexHull transpose matrix{
    {1,1,0,0},
    {1,0,1,0},
    {1,0,0,1},
    {0,1,1,0},
    {0,1,0,1},
    {0,0,1,1}
    }
    )

-- SC-sequence
SCseq = getSequenceSC P 
assert(SCseq == ({2,3},2,{1,2},{0,1},0))


-- creates a "master" system
masterGS = makeSliceSystem(F,G,SCseq)
pp = flatten entries parameters masterGS
assert(
    gateMatrix masterGS === gateMatrix F || transpose gateMatrix{
    	{pp#0*z+pp#1*w+pp#2,
     	    pp#3*x+pp#4*y+pp#5*z+pp#6*w+pp#7
     	    }
 	}
    )

-- all above steps are executed by `witnessCurve`
blocks = {{x},{y},{z},{w}} -- TODO: we should decide between "block" and VariableGroup!!!
wc = witnessCurve(F,blocks,pt) 

-- populate a witness set 
x0 =point sub(matrix pt,CC)
populate wc
assert (length points wc == 12)

-- membership
x1 = point{{-2_CC,1,0,0}}
x2 = point{{0,0,0,0}}

-- TODO: write `membershipTest(WitnessCurve,Point)`
-- note: 
-- (1) we need a complete witness set  
-- (2) check if Point is general; if so, run a segment homotopy from Point toward the witness set  
-- (3) if nongeneric... need more care!!!
end--

restart 
needs "octahedron.m2"
    

