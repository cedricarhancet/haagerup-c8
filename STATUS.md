# Scope and verification status

This file records what the computations establish and what remains open.

## Established by the recorded checks

- (F4): 6010 pivots, 131 nonzero Cartan blocks, 18368 family inputs and
  9800 additional products; exact memberships and exhaustive aggregation passed.
- All 886 zero products: 786 direct mode computations and 100 exclusions by
  Cartan energy.
- C4: seven exact commutant vectors, Virasoro matrices and primary vectors;
  280 central pivots agree with the global (F4) block.
- Zhu rewriting: 131698 overlap and 151 inclusion ambiguities checked.
- Zhu presentation modulo 8388571: dimension 762, computed from the recorded
  projected relations.

## Independent product comparisons

A second implementation of the vertex-operator products, by simultaneous
normal ordering and direct exponential coefficients, is compared with the
main implementation (iterate formula and Schur recurrence). Exact comparisons
agree for:

- 22188 products of weight-one basis states;
- 2455 stratified basis-state products in the weight window used by the proofs;
- 310 products of the six currents and the ten weight-two generators;
- 16 products of the coefficient-field basis and three Virasoro identities.

The 310 generator products have a portable replay. Their input coefficients
are read from the existing data, not reconstructed independently from the
printed definitions.

The main implementation stores lattice momenta of squared norm at most 8 and
can return an incorrect zero outside its output window. The manuscript proves
that all products used in the certificates stay inside this window; a wrapper
in checks/independent-modes rejects inputs or outputs beyond weight four.

## Open points

Agreement at two primes does not prove semisimplicity of the Zhu corner in
characteristic zero, nor that its twelve blocks survive. The computation at
the second prime is supplied but was not rerun.

The exhaustive branching and an independent categorical identification remain
open. The complete classification uses hypothesis (H_J) or (H_R), each of
which includes strong rationality, positivity of the nonvacuum lowest weights
and the braided equivalence with Z(H_3). Identifying both named constructions
uses both hypotheses.

The independent product comparisons are neither a second complete
regeneration of the (F4) coordinates nor a formal verification. Coverage of
(F4) is determined by the final complete receipt, not by intermediate counters.
