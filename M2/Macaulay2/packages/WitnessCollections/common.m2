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
--pt = point{{4,-2/3,-1.44419, .765304}}
--(m,n,JF)=createJacobian F
--Jpt = evaluateJacobian(pt,m,n,JF)