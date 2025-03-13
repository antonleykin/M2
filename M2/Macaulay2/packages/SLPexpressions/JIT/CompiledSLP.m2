export {
    "setTryJustInTimeCompilation",
    "CompiledSLProgram", "makeCompiledSLProgram"
}
setTryJustInTimeCompilation = method()
setTryJustInTimeCompilation Boolean := v -> if v then (
    w := (run "gcc --version" == 0);
    if w then (
	print "-- SLPexpressions: Found `gcc`. Just-in-time compilation will be attempted by `makeSLProgram`.";
    	print "-- This is an experimental feature that works only for evaluation over real and complex numbers.";
	print "-- (To disable, `setTryJustInTimeCompilation false`)" 
	) else (
	print "-- SLPexpressions: Couldn't find `gcc` --- `makeSLProgram` will output `InterpretedSLProgram`.";
	);
    	TryJustInTimeCompilation = w
    ) else TryJustInTimeCompilation = false

TryJustInTimeCompilation = false

makeCompiledSLProgram = method(TypicalValue=>CompiledSLProgram)
makeCompiledSLProgram (List,List) := (inL,outL) -> (
    new CompiledSLProgram from {
	"input" => inL,
	"output" => outL,
	cache => new CacheTable 
	}
    )
makeCompiledSLProgram (GateMatrix,GateMatrix) := (inM,outM) -> makeCompiledSLProgram(flatten entries inM, flatten entries outM)
numberOfInputs CompiledSLProgram := slp -> #(slp#"input")
numberOfOutputs CompiledSLProgram := slp -> #(slp#"output")

rawSLEvaluatorK (CompiledSLProgram, Ring) := (slp, K) -> if slp.cache#?K then 
slp.cache#K else (
    typeName := (
	if K === RR_53 then "double" else 
	if K === CC_53 then "std::complex<double>" else 
    	error ("just-in-time compilation is not implemented for "| toString K) 
    	);
    fname := temporaryFileName() | "-GateSystem";
    cppName := fname | ".cpp";
    --cppName := fname | ".c";
    libName := fname | ".so";
    f := openOut cppName;
    f << "#include <complex>" << endl; 
    f << "std::complex<double> ii(0,1);" << endl;
    f << "typedef " | typeName | " C;" << endl; -- << "extern" << endl; -- the type needs to be adjusted!!!
    cCode (slp#"output", slp#"input", f);
    f << close;
    compileCommand := "g++ -shared -Wall -fPIC -Wextra -O3 -o "| libName | " " | cppName;
    --compileCommand := "gcc -shared -Wall -fPIC -Wextra -o "| libName | " " | cppName;
    print compileCommand;
    if run compileCommand > 0 then error ("error compiling a straightline program:\n"|compileCommand);      
    print get cppName;
    print libName;
    symNames := get ("!nm "|libName); 
    (a,b) := first regex("[0-9a-zA-Z_]*evaluate[0-9a-zA-Z_]*", symNames);
    print ("mangled function name: "|substring(symNames,a,b)); 
    slp.cache#K = rawCompiledSLEvaluator(libName, #(slp#"input"), #(slp#"output"),
	 raw mutableMatrix(K,0,0) -- we need to pass only the field 
	 )
    )

evaluate(CompiledSLProgram, MutableMatrix, MutableMatrix) := (slp,I,O) -> (
		--if numrows I =!= 1 or numrows O =!= 1 then error "expected matrices with 1 row";
		if numrows I * numcols I =!= #(slp#"input") then error "wrong number of inputs";
		if numrows O * numcols O =!= #(slp#"output") then error "wrong number of outputs";
		K := ring I; 
    if ring O =!= K then error "expected same Ring for input and output";
    rawSLEvaluatorEvaluate(rawSLEvaluatorK(slp,K), raw I, raw O);
    )
