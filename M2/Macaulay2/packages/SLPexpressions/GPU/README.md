# GPU SLP Project To-Do List

## Overview
Make an GIT-compiled evaluator for SLPexpressions.

## Tasks

### High Priority
- [ ] Break the procedure in `GPU/cl-SLPtoGPUs.cpp` into
  - initialization
  - registering kernels
  - execution of a kernel
  - finalization

  Compile as a DLL.
- [ ] Use `ForeignFunctions` to invoke functions from DLL.

### Medium Priority

### Low Priority

## Completed Tasks
- [x] `gpuCode` produces a "kernel" that evaluates a `GateMatrix`.
- [x] Tested `gpuCode` by compiling `GPU/cl-SLPtoGPUs.cpp` as a standalone program: 
```
g++ cl-SLPtoGPUs.cpp -framework OpenCL
```

## Notes
Platform dependent. 
[Mar 2025] Everything experiments are done on Apple M1 Pro (GPU with 16 cores).