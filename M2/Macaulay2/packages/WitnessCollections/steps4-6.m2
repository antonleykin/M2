needs "step3.m2"

-- override buggy code in MonodromySolver
createSeedPair (System, Point) := o -> (P, x0) -> (
    G := if instance(P, GateSystem) then P else gateSystem P.PolyMap;
    n := numVariables G;
    m := numParameters G;
    N := numFunctions G;
    I := id_(CC^m);
    A := random(CC^0,CC^N);
    scan(m, i -> A = A || evaluate(G, point I_{i}, x0));
    b := evaluate(G, point matrix 0_(CC^m), x0);
    K := numericalKernel(transpose A, 1e-5) ;
    offset := solve(transpose A,transpose b,ClosestFit=>true);
    p0 := point(K* random(CC^(numcols K), CC^1) - offset);
    (p0, x0)
    )

x0 =point sub(matrix pt,CC)
(p0, x0)=createSeedPair(MasterGS,x0)
norm evaluate(MasterGS,p0,x0)

end--

restart
needs "steps4-6.m2"

(V,npaths)=monodromySolve(MasterGS,p0,{x0},Verbose=>true,NumberOfNodes=>3)

-- STEP 5: trace test
p1 = V.BasePoint
p1Sols = points V.PartialSols
assert(#p1Sols == 12)

-- 2nd parallel slice
p2 = point{drop(coordinates p1,-1)|{random CC}}
MasterGS12 = specialize(parametricSegmentHomotopy MasterGS,transpose((matrix p1)|(matrix p2)))
p2Sols = trackHomotopy(MasterGS12, p1Sols)
p2Sols/(s->norm evaluate(MasterGS,p2,s))

-- third parallel slice
p3 = point{drop(coordinates p1,-1)|{random CC}}
MasterGS13 = specialize(parametricSegmentHomotopy MasterGS,transpose((matrix p1)|(matrix p3)))
p3Sols = trackHomotopy(MasterGS13, p1Sols)
p3Sols/(s->norm evaluate(MasterGS,p3,s))

v1=matrix{sum(p1Sols/coordinates)}
v2=matrix{sum(p2Sols/coordinates)}
v3=matrix{sum(p3Sols/coordinates)}
min first SVD matrix(v2-v1||v3-v1) -- trace test passes!

-- STEP 6 (and 7): membership test
member (Point, GateSystem, Sequence) := (x1, F, MonodromyResult) -> (
    (MasterGS, p1, p1Sols) := MonodromyResult;
    assert(instance(MasterGS,GateSystem) and instance(p1,Point) and instance(p1Sols,List));
    sliceTargParams := first createSeedPair(MasterGS, x1);
    MembershipHomotopy := specialize(parametricSegmentHomotopy MasterGS, transpose((matrix p1)|(matrix sliceTargParams)));
    memberTargets := trackHomotopy(MembershipHomotopy,p1Sols);
    any(memberTargets,x->norm(coordinates x-coordinates x1)<1e-4) -- x1 is a member
    )
    
-- example
MonodromyResult = (MasterGS, p1, p1Sols)
x1 = point{{-2_CC,1,0,0}}
member(x1,F,MonodromyResult)
x2 = point{{0,0,0,0}}
member(x2,F,MonodromyResult)