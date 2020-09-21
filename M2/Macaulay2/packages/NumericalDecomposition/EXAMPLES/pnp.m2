-- inhomogeneous coordinates
restart
needsPackage "NumericalDecomposition"
load "bag-of-newton.m2"
needs "../nag_supplemental.m2"
X=gateMatrix{vars{x,y,z}}
P=gateMatrix{toList vars(a_1..a_9)}
f1 = a_1*(y^2 + z^2) - a_2 *y*z - a_3
f2 = a_4*(z^2 + x^2) - a_5 *x*z - a_6
f3 = a_7*(y^2 + x^2) - a_8 *x*y - a_9
F=transpose gateMatrix{{f1,f2,f3}}
J=diff(X,F)
M=allMinors(J,3,Laplace => true); -- I get an error w/ false
D=gateSystem(P|X,F||transpose M)
elapsedTime x0 = first newtonBag(D,1)
elapsedTime Polymatroid = multiaffineDimension(D,{toList(3..11),toList(0..2)},x0)
vertices Polymatroid

-- witness set w/ "P-slicing preference"
blocks = {flatten entries P, flatten entries X}
p0 = point random(CC^0,CC^1)
Vg=squareDown(p0,x0,D)
-- witness
wc = witnessCurve(Vg, blocks, x0)
populate(wc,Verbose=>true)
--pseudo-witness
wc = witnessCurve(Vg, blocks, x0)
populate(wc,Verbose=>true, Equivalencer=>x->(m:=matrix x; point(m_{3..11})));
Tx0 = (numericalKernel(evaluateJacobian(D, x0), 1e-5));
codimD = 8 - numericalRank(Tx0^{3..11},Threshold=>1e-10)

--try different blocks
blocks = {flatten entries X, flatten entries P}
p0 = point random(CC^0,CC^1)
Vg=squareDown(p0,x0,D)
-- witness
wc = witnessCurve(Vg, blocks, x0)
populate(wc,Verbose=>true)

-- "discriminant" via FGLM
restart
FF = QQ
FFF = frac(FF[a..f])
R1=FFF[x,y,z]
R2=FFF[gens R1,MonomialOrder=>Lex]
use R1
f1 = (y^2 + z^2) - a *y*z - b
f2 = (z^2 + x^2) - c *x*z - d
f3 = (y^2 + x^2) - e *x*y - f
G = gb ideal(f1,f2,f3)
needsPackage "FGLM"
elapsedTime G2=fglm(G, R2);
COEFS = apply(
    flatten entries last coefficients first flatten entries gens G2
    ,p->sub(p, FFF))
A = denominator last COEFS
B = numerator COEFS#1
C = numerator COEFS#2
D = numerator COEFS#3
E = numerator COEFS#4
factor E
-*
            2    2           2                  2 2
o23 = (b*d*e  - b  - 2b*d - d  + 2b*f + 2d*f - f )

o23 : Expression of class Product
*-
S=QQ[aa,bb,cc,DD,ee,x]
f=aa*x^8+bb*x^6+cc*x^4+DD*x^2+ee
DISC=discriminant(f,x)
factor DISC
-*
              2   2  2  2         3  2      3  3                3       2  4      2  3            4         3                     2              2  2          2     2         4  2           2     2        2  2  2        2        2        3  3 2
o21 = (ee)(aa) (bb cc DD  - 4aa*cc DD  - 4bb DD  + 18aa*bb*cc*DD  - 27aa DD  - 4bb cc ee + 16aa*cc ee + 18bb cc*DD*ee - 80aa*bb*cc DD*ee - 6aa*bb DD ee + 144aa cc*DD ee - 27bb ee  + 144aa*bb cc*ee  - 128aa cc ee  - 192aa bb*DD*ee  + 256aa ee ) (256)
*-
-- subbing aa=>A, ..., makes each irreducible factor a square
-- can we see this numerically somehow?

-- "discriminant" via elimination (in affine chart)
restart
FF = QQ
R=FF[x,y,z,a..i,MonomialOrder=>Eliminate 3]
f1 = g*(y^2 + z^2) - a *y*z - b
f2 = h*(z^2 + x^2) - c *x*z - d
f3 = i*(y^2 + x^2) - e *x*y - f
I=ideal(f1,f2,f3)
Icrit = I+ideal det (jacobian I)^{0..2}
elapsedTime G = groebnerBasis(Icrit, Strategy => "F4");
D1 = (selectInSubring(1, G))_(0,0);
factor D1 -- square-free
degree D1


