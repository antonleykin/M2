needsPackage "NumericalSchubertCalculus"
setRandomSeed 2

--Problem W^6W111^2 in G(3,7)
 --a problem with 5 solutions

print("Solving problem W^6W111^2 in G(3,7)");

SchPblm = randomSchubertProblemInstance(
  {{1, 1, 1},{1, 1, 1},{1},{1},{1},{1},{1},{1}},3,7);
time S = solveSchubertProblem(SchPblm, 3,7);
assert all(S,s->checkIncidenceSolution(s, SchPblm))

 end
 ------

restart
 load"NumericalSchubertCalculus/EXA/ProblemsG37/W^6W111^2-G37.m2"

