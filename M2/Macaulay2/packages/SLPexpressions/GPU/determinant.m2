restart
debug needsPackage "SLPexpressions"
I = gateMatrix{{a,b,c,d} / declareVariable}
O = gateMatrix{{a*d-b*c}}
out = openOut "determinant.cl"
gpuCode(O,I,out)
close out
