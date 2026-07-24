# ProofConcept — design conversation & rationale

This document is the "memory" of the conversation that produced the
`ProofConcept` package: the original request, how the design evolved, the
technical findings about Macaulay2 that shaped it, and the decisions made along
the way. It is a companion to `ProofConcept.m2` and is not loaded by Macaulay2.

## 1. The request

> Create a package `ProofConcept`. We'd like to have data types for
> * an implication of the type `A => B`;
> * an atomic proof: a certificate (verifying an implication);
> * a proof: a collection of implications that prove a given implication.

Three motivating examples were given:

1. **Ideal membership.** `A` = "`G` are generators of an ideal `I` in a
   polynomial ring and `f` is another polynomial"; `B` = "`f` is in `I`".
   *Certificate:* polynomials `H` such that `H·G = f`.
2. **Monomial ideal membership.** Same, but `I` is a monomial ideal.
   *Certificate:* divisibility of `f` by one of the monomial generators.
3. **Gröbner basis.** `A` = "`G` generate `I`"; `B` = "`G` is a Gröbner basis
   with respect to a monomial order". *Proof (non-atomic):* Buchberger's
   criterion.

## 2. A follow-up that reshaped the design

After the first plan was drafted, the user asked:

> I suppose your suggestion of proof steps just covers a linear chain of
> implications, correct? Can we add some simple logic (and, or, not, ...) to
> proofs?

This was the pivotal design change. The initial model had a `Proof` carrying a
flat list of `Steps` — effectively a linear chain. The revision turns a proof
into a **tree** whose internal nodes are logical inference rules and whose
leaves are certificates. Statements become **logical formulas** that can be
combined with `and`, `or`, `not`.

## 3. Technical findings from exploring Macaulay2

These determined what was idiomatic and feasible:

- **Package skeleton.** `Macaulay2/packages/PackageTemplate.m2` shows the
  expected shape: `newPackage`, `export`, body, `beginDocumentation()`,
  `document {…}`, `TEST ///…///`, `end`. New types via `T = new Type of …`;
  subtypes via `S = new Type of T`.
- **Membership certificate primitive.** `f // (gens I)` returns the coefficient
  column `H`, with the identity `f == (gens I)*(f // gens I) + (f % gens I)`
  (`m2/matrix2.m2`). Membership itself is `f % I == 0` / `isMember(f, I)`.
- **Monomial ideals.** `MonomialIdeal` is a subtype of `Ideal`, so one
  `inIdeal(RingElement, Ideal)` constructor covers both; `certify` dispatches
  separately on `Ideal` vs `MonomialIdeal`.
- **Faithful Buchberger test.** `forceGB (matrix {G})` wraps the *given*
  generators as a Gröbner basis so that `S % forceGB …` reduces using only
  those lead terms. Reducing by `ideal G` directly would always give `0` and
  prove nothing. S-polynomials are not user-exposed, so they are built by hand
  from `leadMonomial`, `leadTerm`, and an exponent-vector `lcm`.
- **Overloadable logical operators.** The interpreter (`d/actors.d`) dispatches
  `and`/`or`/`not` to installed methods when the operands are *not* Boolean
  (`binarymethod(…, andS/orS)`, `unarymethod(…, notS)`). So
  `Statement and Statement`, `Statement or Statement`, and `not Statement`
  install cleanly with the ordinary `:=` syntax (cf. `m2/integers.m2`:
  `Function and Function := …`).
- **`==>` is free.** `==>` is a right-associative binary operator
  (`d/binding.d`) with no built-in method and a LaTeX rendering of `⟹`, so it
  is used directly as the implication operator (`A ==> B`).

## 4. Decisions made (with the user)

- **Standalone, not in CI.** The user chose *not* to register the package in
  `=distributed-packages`, so it is loadable via `loadPackage "ProofConcept"`
  but is not built/checked in CI. This keeps the proof-of-concept low-risk.
- **Scope: concept + all three examples.** Core types plus `verify`, with the
  three worked cases (ideal membership, monomial-ideal membership, Buchberger),
  documented and tested — rather than a minimal skeleton or a fully extensible
  certificate framework.
- **Simple logic, not full natural deduction.** Supported rules are
  and-introduction, or-introduction, negation certificates, and the derived
  Buchberger rule. Hypothesis discharge / `or`-elimination case analysis and
  the `<==>`/`xor` connectives are intentionally left as future work.

## 5. Resulting API

Types: `Statement` (with subtypes `And`, `Or`, `Not`), `Implication`,
`Certificate`, `Proof`.

- Atomic statements: `inIdeal(f, I)`, `generates(G, I)`, `isGroebner(G)`.
- Connectives: `A and B`, `A or B`, `not A`; implication `A ==> B`
  (function form `implies(A, B)`); structural equality `A == B`.
- Certificate builders: `certify(f, I)` (ideal or monomial-ideal membership),
  `refute(f, I)` (non-membership, witnessed by the nonzero normal form).
- Proof-tree builders: `andProof(…)`, `orProof(p, S)`, `buchbergerProof(G)`.
- Operations: `verify` (independently re-checks the witness in each leaf and
  the side conditions of each rule, returning a `Boolean`), `conclusion` (the
  `Statement` a certificate/proof establishes), `proves(P, S)`.

The unifying idea: a proof is *data*, and `verify` re-derives and re-checks it
rather than trusting how it was produced.

## 6. Status

`ProofConcept.m2` was implemented, documented, and given `TEST` blocks for the
three cases and the logic layer. No compiled Macaulay2 binary was available in
the authoring environment, so the package was validated by review; the `TEST`
blocks are intended to be run with `check ProofConcept` after
`loadPackage "ProofConcept"` (or `installPackage "ProofConcept"`).
