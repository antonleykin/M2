--common functions that pieces will need 
needsPackage "MonodromySolver"
needsPackage "Polyhedra"
VariableGroup = new Type of List

--Service functions
createJacobian = method()
createJacobian GateSystem := GS -> (
    F := gateMatrix GS;
    I := vars GS;
    J := diff(I,F);
    (m,n) := (numrows J,numcols J);
    JGS := gateSystem(I, flatten J);
    (m,n,JGS)
    )

evaluateJacobian (Point, ZZ, ZZ, GateSystem) := (pt, m,n,JGS) -> matrix(evaluate(JGS,matrix pt),m,n)

sliceParam = symbol sP
sliceParamCounter = 0
newSliceParam = () -> (
    ret := inputGate sliceParam_sliceParamCounter;
    sliceParamCounter = sliceParamCounter + 1;
    ret
    )

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
