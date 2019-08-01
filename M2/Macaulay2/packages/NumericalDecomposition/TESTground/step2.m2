-- JOSE's playground
restart
needs "step1.m2"
needs "step2.m2"

declareVariable \ {x,y,z,w}
f = 1 + 2*x + 3*y^2+4*z^3+5*w^4
h = 1 + 2*x + 3*y+ 5*z + 7*w + 11*z*y + 13 * x* z + 17*x*w + 19*y*z+23*y*w +29*z*w+31*x*y*z+37*x*y*w+41*x*z*w+43*y*z*w+47*x*y*z*w
F = gateSystem(matrix{{x,y,z,w}},transpose gateMatrix{{f,h}})
G = new VariableGroup from {{0},{1},{2},{3}}
pt = point{{4,-2/3,-1.44419, .765304}}
P = multiaffineDimension(F,G,pt)
latticePoints P 
getSequenceSC (P)

declareVariable \ {d,x,y,z,w};
f = 1 + 2*x + 3*y^2+4*z^3+5*w^4;
h = 1 + 2*x + 3*y+ 5*z + 7*w + 11*z*y + 13 * x* z + 17*x*w + 19*y*z+23*y*w +29*z*w+31*x*y*z+37*x*y*w+41*x*z*w+43*y*z*w+47*x*y*z*w;
f3 = d;
F = gateSystem(matrix{{d,x,y,z,w}},transpose gateMatrix{{f,h,f3}});
G = new VariableGroup from {{0},{1},{2},{3},{4}};
pt = point{{0,4,-2/3,-1.44419, .765304}};
P = multiaffineDimension(F,G,pt);
latticePoints P ;
getSequenceSC (P)

X = {x,y,z,w,d}
declareVariable \ X
f = 1 + 2*x + 3*y^2+4*z^3+5*w^4
h = 1 + 2*x + 3*y+ 5*z + 7*w + 11*z*y + 13 * x* z + 17*x*w + 19*y*z+23*y*w +29*z*w+31*x*y*z+37*x*y*w+41*x*z*w+43*y*z*w+47*x*y*z*w
f3 = d
F = gateSystem(matrix{X},transpose gateMatrix{{f,h,f3}})
G = new VariableGroup from {{0},{1},{2},{3},{4}}
pt = point{{0,4,-2/3,-1.44419, .765304}}
P = multiaffineDimension(F,G,pt)
latticePoints P 
getSequenceSC (P)









