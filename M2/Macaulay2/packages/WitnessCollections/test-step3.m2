needs "step3.m2"

-- (mostly) COPIED FROM OCTAHEDRON
declareVariable \ {x,y,z,w}
f = 1 + 2*x + 3*y^2+4*z^3+5*w^4
h = 1 + 2*x + 3*y+ 5*z + 7*w + 11*z*y + 13 * x* z + 17*x*w + 19*y*z+23*y*w +29*z*w+31*x*y*z+37*x*y*w+41*x*z*w+43*y*z*w+47*x*y*z*w
F = gateSystem(
    matrix{{x,y,z,w}},
    transpose gateMatrix{{f,h}}
    )
G = new VariableGroup from {{0},{1},{2},{3}}
pt = point{{4,-2/3,-1.44419, .765304}}
evaluate(F,pt)

(m,n,JF)=createJacobian F
Jpt = evaluateJacobian(pt,m,n,JF)

P = convexHull transpose matrix{
    {1,1,0,0},
    {1,0,1,0},
    {1,0,0,1},
    {0,1,1,0},
    {0,1,0,1},
    {0,0,1,1}
    }

SCseq =({2,3},2,{0,1},{0,1},0)


(Params, L) = makeSliceSystem(F,G,SCseq)
MasterGS = gateSystem(
    Params,
    vars F,
    (gateMatrix F) || L
    )
