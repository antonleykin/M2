-- -*- coding: utf-8 -*-
MPIPresent := Core#"private dictionary"#?"numberOfProcesses"

newPackage(
	"MPI",
    Version => "1.24.11", 
    Date => "Feb, 2024",
    	Authors => {
	     {Name => "Anton Leykin", Email => "leykin@math.gatech.edu"}
	     },
    	HomePage => "https://antonleykin.math.gatech.edu/",
    	Headline => "Message Passing Interface",
	Keywords => {"Parallel computing"},
    OptionalComponentsPresent => MPIPresent,
	AuxiliaryFiles => true
    	)

verboseLog = if debugLevel > 0 then printerr else identity    
if MPIPresent then verboseLog "success: MPI is present" else (
    verboseLog "warning: python is not present";
    verboseLog "specify --with-MPI in `configure` options and recompile M2";
    end
    )

exportFrom_Core {"sendStringWithTagMPI", "receiveStringMPI", "interruptMPI", "numberOfProcesses", "myProcessNumber"}
--exportMutable {}

-- TO DO: implement these in the kernel
export {"broadcastSend", "broadcastReceive", "sendStringMPI", "sendErrorMPI"}

export {"MPIstringTag", "MPIerrorTag"}
MPIstringTag = 1; MPIerrorTag = 2;
protect symbol MPIstringTag
protect symbol MPIerrorTag


sendStringMPI = method()
sendStringMPI (String,ZZ) := (s,p) -> sendStringWithTagMPI(s,p,MPIstringTag)

sendErrorMPI = method()
sendErrorMPI (String,ZZ) := (s,p) -> sendStringWithTagMPI(s,p,MPIerrorTag)

broadcastSend = method()
broadcastSend String := s -> for i from 1 to numberOfProcesses()-1 do sendStringMPI(s,i)

broadcastReceive = method()
installMethod(broadcastReceive, 
    () -> for i from 1 to numberOfProcesses()-1 list receiveStringMPI i)


-* 
*** PARADIGM ***

process 0 is the MANAGER (possibly interactive)
the rest of processes are WORKERS

communication: 
  * MANAGER sends WORKER a task (a string containing lines of M2 code to execute)
  * WORKER completes the task and sends the result back to MANAGER (a string containing the output of the last line) 

To build:
  * Install MPI: e.g. `brew install open-mpi`
  * `../../configure  CC=mpicc CXX=mpic++ --enable-download`  

To run: 
  * see comments in `MPI/master-worker-MPI.m2`
  
*** DONE *** 

Core functions:
* send a string to a process
* receive a string from process
* get process number (rank)
* get number of processes 


*** TO DO ***

* probe for a message (from any process)
*? "spawn" new workers
* broadcast
* write a package containing convenience parallelization routines
*? make workers prepend error messages (or anything that goes to stdin/stderr) with their ID. 
 
*** USE CASES ***

Possible use case:
  bunch of Groebner bases computed, i.e. using different monomial orders.
  run through all of them, receive from each, the initial monomials, and/or GB elements themselves.

Possible use case:
  queue of tasks: input, output, status (processor, actual status). (e.g. doRandomBinomialIdealSearch(R, 1000)), tag a task (e.g. an index).
  list of completed tasks?
  indicator function: is a worker available? (mutable list of (string sent, or null)) (or three states...)
    -- once we do receiveStringMPI(worker), we set this to null for that worked.
  freeWorker (returns null if none are free), or block+probe.
    block until one is free.

  need a loop: task manager which distributes all the tasks:
    probe workers (blocks until one is free)
      if probe tells us we have a result: we store it in the task output, set it to be free
      if one is free: grab the next task, send it to that worker.
    
E.g.: Compute random binomial ideals, find their betti numbers, and keep only the ones with regularity >= n, some fixed n.

*-

beginDocumentation()
TEST ///
master = 0
numberOfWorkers = numberOfProcesses()-1
myID = myProcessNumber()
notDone = true
if myID!=master then while true do (
    (tag,s) := receiveStringMPI master;
    << "-- " << myID << " received(tag): " << tag << endl;
    << "-- " << myID << " received(str): " << s << endl;
    r := value s;
    << "-- " << myID << " result: " << r << endl;
    sendStringMPI(toString r, master);
    ) else (
    addEndFunction(()->(for i from 1 to numberOfWorkers do sendStringMPI("exit 0",i); -*sleep 1*-));
    )

-- write the code for master below
for i from 1  to numberOfWorkers do ( 
    s = toString i | "+" | toString i;
    << "-- " << myID << " sent: " << s << endl;
    sendStringMPI(s,i);
    )
for i from 1  to numberOfWorkers do ( 
    r = receiveStringMPI i;
    << "-- " << myID << " received: " << r << endl;
    )
///
         
end

-- Here place M2 code that you find useful while developing this
-- package.  None of it will be executed when the file is loaded,
-- because loading stops when the symbol "end" is encountered.

installPackage "PackageTemplate"
installPackage("PackageTemplate", RemakeAllDocumentation=>true)
check PackageTemplate

-- Local Variables:
-- compile-command: "make -C $M2BUILDDIR/Macaulay2/packages PACKAGES=PackageTemplate pre-install"
-- End:

