needsPackage "MonodromySolver"
needsPackage "SLPexpressions"
needsPackage "NumericalDecomposition"
vars(h_(1,1,1)..h_(1,3,3),h_(2,1,1)..h_(2,3,3))
H1 = matrix for i from 1 to 3 list for j from 1 to 3 list h_(1,i,j)
H2 = matrix for i from 1 to 3 list for j from 1 to 3 list h_(2,i,j)
H=matrix{flatten entries H1 | flatten entries H2}
ONE = inputGate 1
ZERO = inputGate 0
I=gateMatrix for i from 1 to 3 list for j from 1 to 3 list if i==j then ONE else ZERO
det2 = M -> M_(0,0)*M_(1,1)-M_(1,0)*M_(0,1)
det3 = M -> M_(0,0)*det2(M_{1,2}^{1,2})-M_(0,1)*det2(M_{0,2}^{1,2})+M_(0,2)*det2(M_{0,1}^{1,2})
det4 = M -> M_(0,0)*det3(M_{1,2,3}^{1,2,3})-M_(0,1)*det3(M_{0,2,3}^{1,2,3})+M_(0,2)*det3(M_{0,1,3}^{1,2,3})-M_(0,3)*det3(M_{0,1,2}^{1,2,3})
H012=transpose H1*H1-I|transpose H2*H2-I;
dets = {
    det3 H012_{1,2,3} - det3 H012_{0,2,4} + det3 H012_{0,1,5},
    det3 H012_{2,3,4} - det3 H012_{1,3,5} + det3 H012_{0,4,5},
    det3 H012_{0,1,2},
    det3 H012_{3,4,5}
};

F = gateSystem(H,transpose gateMatrix{dets})

load "../nag_supplemental.m2"
d = degrees F

-*
STEP 1: get points on nonreduced component via Bezout homotopy w/ "intrinsic slicing" followed by endgame
*-
n = numVariables F
dimSlice = numFunctions F
A = random(CC^n, CC^dimSlice) 
X = transpose gateMatrix{vars {x_1,x_2,x_3,x_4}}
b = random(CC^n,CC^1)
sliceF = gateSystem(transpose X, sub(gateMatrix F, vars F, transpose(A*X+b)))
vars t
H = gateHomotopy(
    t * gateMatrix sliceF + (1-t) * (random CC)*gateMatrix for i from 0 to 3 list {X_(i,0)^(d#i) - 1},
    transpose X,
    t
    )
startSols = apply(d#0, i -> point{{exp(2*pi*ii*i/d#0),1,1,1}})
results = trackHomotopy(H, startSols)
tally(status \ results)
x = first select(1, results, x -> status x != Regular)
y = endGameCauchy(x#"H", 1.0, x)
-- y is a singular point on Vf
Jy = evaluateJacobian(sliceF,y)
yNull = numericalKernel(Jy, 1e-6)
r = numcols Jy - numcols yNull
yNull0 = yNull * random(CC^(numcols yNull), CC^1)
-*
TRY: deflating to order 1
*-
J = diff(vars sliceF, gateMatrix sliceF);
nullJ = gateMatrix{vars{a_0,a_1,a_2,a_3}}
deflatedGS = gateSystem(
    vars sliceF | nullJ,
    gateMatrix sliceF || J * transpose nullJ
    )
getSequenceSC(
    deflatedGS, 
    {{0,1,2,3}, {4,5,6,7}},
    point(matrix y | transpose yNull0)
    )
ya0 = point(matrix y | transpose yNull0)
assert(areEqual(norm evaluate(deflatedGS, ya0), 0))
numericalRank evaluateJacobian(deflatedGS, ya0)
squaredDownDeflatedGS = squareDown(point random(CC^0,CC^0), ya0, deflatedGS)

W=witnessCurve(
    squaredDownDeflatedGS, 
    flatten \ entries \ {(vars deflatedGS)_{0,1,2,3}, (vars deflatedGS)_{4,5,6,7}},
    ya0
)
end
restart
needs "homographies.m2"
populate(W,Verbose=>true) -- ????

-- sanity check w/ points on nonsingular component
x = first select(1, results, x -> status x == Regular)
numericalRank evaluateJacobian(sliceF, x) -- nonsingular
h = point(A*transpose matrix x + b)
W = witnessCurve(F, {flatten entries vars F}, h)
