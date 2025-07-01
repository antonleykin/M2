# RRb and RRbcell Implementation

## Summary

This document describes the implementation of RRb and RRbcell as alternative types to RR and RRcell in Macaulay2.

## What Was Added

### Core Type Definitions

1. **New Type Codes** (in `parse.d`):
   - `TCRRb := 4` - Token type code for RRb literals
   - Updated `TCstring := 5` (was 4)

2. **New Types** (in `gmp.d`):
   - `RRb := Pointer "mpfr_srcptr"` - Alternative RR type
   - `RRborNull := RRb or null` - Nullable version
   - `RRbcell := {+v:RRb}` - Cell wrapper for expressions
   - `RRbmutable := Pointer "mpfr_ptr"` - Mutable version for computations

3. **New Code Types** (in `parse.d`):
   - `realRRbCode := {+x:RRb,position:Position}` - Parse tree code
   - Added to `Code` union type

4. **New Expression Types** (in `parse.d`):
   - Added `RRbcell` to `Expr` union type

### Lexer and Parser Support

1. **Lexer** (in `lex.d`):
   - Added recognition of 'b' suffix for floating point numbers
   - Numbers like `3.14b` or `2.5p53b` are parsed as `TCRRb` type

2. **Parser** (in `parser.d`):
   - `parseRRb(s:string):RRborNull` - Parses RRb numbers (removes 'b' suffix)

3. **Converter** (in `convertr.d`):
   - Added `TCRRb` handling to convert tokens to `realRRbCode`

### Runtime Support

1. **Evaluation** (in `evaluate.d`):
   - Added `realRRbCode` evaluation support

2. **Conversion Functions** (in `util.d`):
   - `toExpr(x:RRb):Expr` - Convert RRb to expression
   - `toExpr(x:RRborNull):Expr` - Convert nullable RRb to expression

3. **Basic Operations** (in `gmp.d`):
   - `toRRb(s:string, prec:ulong):RRb` - Create RRb from string
   - `toRRb(x:double, prec:ulong):RRb` - Create RRb from double
   - `toRRb(x:RR):RRb` - Convert RR to RRb
   - `toFloat(x:RRbcell):float` - Convert to float
   - `toDouble(x:RRbcell):double` - Convert to double

4. **Debugging Support** (in `debugging.dd`):
   - Added `realRRbCode` to string conversion
   - Added `realRRbCode` to expression conversion

5. **Equality Operations** (in `actors3.d`):
   - Added RRbcell equality comparisons with ZZ, QQ, RR, RRb, RRi, CC

6. **Hash Support** (in `basic.d`):
   - Added `hash(x:RRbcell)` function

7. **Position Support** (in `common.d`):
   - Added `realRRbCode` position handling

## Current Implementation Status

### What Works
- ✅ Lexing of RRb literals (e.g., `3.14b`)
- ✅ Parsing of RRb literals 
- ✅ Basic type conversion
- ✅ Equality comparisons
- ✅ Hash functions
- ✅ Debug printing
- ✅ Expression conversion

### What's Missing
- ❌ Arithmetic operations (+, -, *, /, ^)
- ❌ Mathematical functions (sin, cos, exp, log, etc.)
- ❌ Comparison operations (<, >, <=, >=)
- ❌ String formatting and output
- ❌ Integration with the engine (e/ directory)

## Usage Example

Once fully implemented, users could write:
```macaulay2
x = 3.14b        -- Creates an RRb number
y = 2.5p53b      -- Creates an RRb with precision 53
```

## Current Limitation

The current implementation treats RRb as identical to RR at the underlying level. For RRb to be truly "alternative," additional differentiation would need to be implemented in the arithmetic and mathematical functions.

## Next Steps

To complete the implementation:

1. Add arithmetic operators for RRbcell in `actors2.dd` and `actors3.d`
2. Add mathematical functions for RRbcell in `actors3.d`
3. Add comparison operations
4. Add string formatting support in `actors4.d` and `actors5.d`
5. Consider adding engine-level support in the `e/` directory
6. Add comprehensive test cases

## Files Modified

- `M2/Macaulay2/d/parse.d` - Type definitions
- `M2/Macaulay2/d/gmp.d` - Core RRb functions
- `M2/Macaulay2/d/lex.d` - Lexer support
- `M2/Macaulay2/d/parser.d` - Parser support
- `M2/Macaulay2/d/convertr.d` - Code conversion
- `M2/Macaulay2/d/evaluate.d` - Evaluation support
- `M2/Macaulay2/d/util.d` - Utility functions
- `M2/Macaulay2/d/debugging.dd` - Debug support
- `M2/Macaulay2/d/actors3.d` - Equality operations
- `M2/Macaulay2/d/basic.d` - Hash functions
- `M2/Macaulay2/d/common.d` - Position handling