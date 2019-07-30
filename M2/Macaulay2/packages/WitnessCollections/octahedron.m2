restart
needs "WitnessCollections.m2"
needsPackage "MonodromySolver"
declareVariable \ {x,y,z,w}
f = 1 + 2*x + 3*y^2+4*z^3+5*w^4
h = 1 + 2*x + 3*y+ 5*z + 7*w + 11*z*y + 13 * x* z + 17*x*w + 19*y*z+23*y*w +29*z*w+31*x*y*z+37*x*y*w+41*x*z*w+43*y*z*w+47*x*y*z*w
F = transpose gateMatrix{{f,h}}
VariableGroup = new Type of List
G = new VariableGroup from {{x},{y},{z},{w}}
I = transpose matrix G

-*
--getting the point below
R=QQ[W,T]
s = makeSLProgram(matrix G,F)
solveSystem flatten entries evaluate(s,matrix{{4,-2/3,W,T}})
*-

-*
IN: J: a gateMatrix
    I: its matrix of inputs
    pt: a point   
evalJ = (pt,I,J) -> (
    (m,n) := (numrows J,numcols J);
    JProg := makeSLProgram(I, J);
    matrix(evaluate(JProg,matrix pt),m,n)
    )    
    

J = diff(I,F)
pt = point{{4,-2/3,-1.44419, .765304}}
Jpt = evalJ(pt,I,J)
thePartialJacs = G/(g->(
        ncols := #g;
        Jpt^{g}
	)
    )


-- dim rewrite
dimensionPolytope(Point,VariableGroup,GateMatrix) := (pt,G,F) -> (
        
    