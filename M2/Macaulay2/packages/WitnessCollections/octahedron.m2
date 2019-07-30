restart
--needs "WitnessCollections.m2"
needsPackage "MonodromySolver"
declareVariable \ {x,y,z,w}
f = 1 + 2*x + 3*y^2+4*z^3+5*w^4
h = 1 + 2*x + 3*y+ 5*z + 7*w + 11*z*y + 13 * x* z + 17*x*w + 19*y*z+23*y*w +29*z*w+31*x*y*z+37*x*y*w+41*x*z*w+43*y*z*w+47*x*y*z*w
F = gateSystem(matrix{{x,y,z,w}},transpose gateMatrix{{f,h}})
VariableGroup = new Type of List
G = new VariableGroup from {{0},{1},{2},{3}}

pt = point{{4,-2/3,-1.44419, .765304}}

createJacobian = method()
createJacobian GateSystem := GS -> (
    F := gateMatrix GS;
    I := vars GS;
    J := diff(I,F);
    (m,n) := (numrows J,numcols J);
    JGS := gateSystem(I, flatten J);
    (m,n,JGS)
    )

evaluateJacobian = method()
evaluateJacobian (Point, ZZ, ZZ, GateSystem) := (pt, m,n,JGS) -> matrix(evaluate(JGS,matrix pt),m,n)
    
(m,n,JF)=createJacobian F
Jpt = evaluateJacobian(pt,m,n,JF)

needsPackage "Polyhedra"
-- output of step 1), input to step 2)
P = convexHull transpose matrix{
    {1,1,0,0},
    {1,0,1,0},
    {1,0,0,1},
    {0,1,1,0},
    {0,1,0,1},
    {0,0,1,1}
    }

-- output of step 2
SCseq = ({1,2},1,{0,1},{0,1})

-- output of step 3 
-*
-- outline for step 3 function
-- output: 
  -- Params: List of parameters
  -- L: GateMatrix representing the slice
makeSlices = (I,G,SCseq) := -> (
    
    
*-
sliceParam = symbol sP
sliceParamCounter = 0

newSliceParam = () -> (
    ret := sP_sliceParamCounter;
    sliceParamCounter = sliceParamCounter + 1;
    ret
    )

Params=declareVariable \ apply(8,i->newSliceParam())
L = transpose gateMatrix{
    {Params#0*y+Params#1*z+Params#2,
     Params#3*x+Params#4*y+Params#5*z+Params#6*w+Params#7
     }
 }
	
MasterGS = gateSystem(
    gateMatrix{Params},
    vars F,
    (gateMatrix F)||L
    )
x0 =point sub(matrix pt,CC)
createSeedPair(MasterGS,x0) -- bug?

-*
G=MasterGS
n := numVariables G;
m := numParameters G;
N := numFunctions G;
I := id_(CC^m);
A := random(CC^0,CC^N);
scan(m, i -> A = A || evaluate(G, point I_{i}, x0));
b := evaluate(G, point matrix 0_(CC^m), x0);
K := numericalKernel(transpose A, 1e-5) ;
offset := solve(transpose A,transpose b,ClosestFit=>true);
p0 := point(K* random(CC^(m-n), CC^1) - offset);
(p0, x0)

*-