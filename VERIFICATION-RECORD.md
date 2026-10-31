# Verification provenance

This record preserves the chronological verification information moved out
of the manuscript. Version numbers refer to when checks were performed;
their retention does not claim that they were repeated in v05--v07. Statement
references agree with the v07 numbering. Use COMPUTATIONS.md to
resolve the filenames. The source code and receipts are unchanged.

## Recorded computations

The accompanying files distinguish the original computational sources from corrected verifiers and newly computed outputs. The original sources are recorded in calculs-sources/MANIFEST-sources.sha256. The following checks were rerun for version 02 and are retained with their outputs; they are not all claimed as new reruns in versions 03 or 04.

- The augmented systems defining the traces on $U_2$, their ranks and the eighteen resulting moments were checked in exact rational arithmetic. The verifier and its output are in verification-foundations/. The witness used in Proposition 5.4 is now given explicitly in its proof.
- The thirteen rejected candidates for the commutant are excluded by the type $(3,3)$ alone. The verifier check_primary33_exclusions.py computes a unit Groebner basis on each of thirty-nine rational projective charts. The same systems were checked independently with lexicographic order.
- In the pairwise closure exclusions, closure_certify_exact_intervals.py evaluates rational expressions and radicals directly in Arb at $256$ bits. It replaces the earlier conversion through unbounded decimal approximations. It produces $62$ nonzero interval witnesses. The separate verifier check_closure_coverage.py checks that they exclude $254$ of the $256$ choices. An interval containing zero is not used to exclude a candidate.
- The normalized coefficients in the identification with Jiao were checked exactly. The corrected membership computation and the generation computation give initial rank $65$ and rank $120$ after one closure step, at the residue map specified in the proof of Theorem 10.24. The corrected scripts and their outputs are in verification-certificates/.
- The modular calculation of $N'_4$ was regenerated at $p_3 = 8388571$. It gives rank $956$ in weight three and rank $6010$ for $N'_4$. The reduced matrices by Cartan weight, their pivot checks and the checks on integral reduction and integer arithmetic are supplied in verification-weight-four/. This modular rank computation is independent of the characteristic-zero membership computation.
- Identity 11.4 gives an explicit word for the affine simple current. The script check_simple_current_word.py reconstructs its exact nonzero coefficient and checks the $26$ annihilation identities used in Proposition 11.4. It does not read a weight-four pivot file.


The characteristic-zero weight-four computation was reconstructed and verified for version 03. The source manifest identifies $72$ original source or input files and the new verifiers. The complete coordinate domains contain $6010$ pivot words, all $18368$ words of the five families defining $N'_4$, and all $9800$ additional products. The final verifier checks the input labels, source identities of the pivots, integral reductions, full rank of the pivot blocks, and the exact membership identities on every rational coordinate after restriction of scalars. The $131$ Cartan blocks pass, and the two verification parts are combined by an exhaustive aggregation. A separate check verifies all $886$ zero inputs, of which $786$ are recomputed by direct modes and $100$ follow from the Cartan energy bound. The final record is checks/certificates/deliverable/COMPLETE-F4.json. Thus Proposition 10.17 establishes (F4) without an ambient Haagerup algebra or an external existence assertion. Corollary 10.22 is retained as a separate dimension proof when such an algebra is available.

The full coordinates are supplied in the two companion archives Certificats-poids-quatre-v03-partie-1.zip and Certificats-poids-quatre-v03-partie-2.zip, with common root certificats-poids-quatre-v03/. The main dossier includes the scripts, source hashes, exact block receipts, zero checks and replay instructions in checks/certificates/deliverable/. The reader rejects incomplete, truncated, duplicated or unexpected records; the final result is not inferred from the counters in a producer's log. The commutant calculation was also rerun from the regenerated central block, giving its seven exact vectors, the two Virasoro matrices and the primary vectors of Proposition 10.31; their coordinates and a second matrix verifier are supplied in checks/identification/c4/.


The targeted audit incorporated in version 04 reevaluates products by a second mode implementation. It uses simultaneous normal ordering and direct exponential coefficients, whereas the historical engine uses the iterate formula and a Schur recurrence. Exact comparisons agree for all $86^2$ pairs of weight-one basis states in output weights zero, one and two ($22188$ products), and for $2455$ stratified pairs with input and output weights at most four. It also checks $310$ products of the six currents and ten specified weight-two vectors, the $16$ products in the coefficient-field basis, and three Virasoro identities. The generator coefficients in these latter comparisons are taken from the existing setup and are not independently reconstructed from their printed definitions. These checks supplement the complete coordinate verification of (F4); they are not a second regeneration of every coordinate in that certificate. The audit also supplies a wrapper that rejects calls outside the historical engine's weight window. The scripts, portable generator coordinates and exact outputs are supplied in checks/independent-modes/; the original mode engine remains unchanged.

The new checks of version 03 are stored under checks/. They include the six joint-sector searches and the two searches for $U_2$; the dimensions of the rational projective-chart quotients in the Griess enumeration; the sixteen admissible primary lines for the case (II) candidate; thirteen unit-ideal certificates excluding the case (I) candidates; the explicit centralizer lifts and descent actions; the exact weight-three calculation; and a rational LDL enumeration of the lattice ball. The short checks of the simple-current table and of the binary lattices accompany the new analytic proofs. The new runs use Python 3.12; dependency versions and commands are recorded with each verifier, separately from the original runtime described above.

The final Zhu presentation was checked for all $131698$ overlap ambiguities and all $151$ inclusion ambiguities; the latter were not explicitly covered by the original checker. The final algebra calculation at $8388571$ was rerun from the archived relations, giving the same dimension $762$, nondegenerate trace form and twelve blocks. These checks concern the recorded presentation and do not recompute the VOA origin of each relation, prove a characteristic-zero lift, or show that every block survives in the Zhu corner. The files for the second prime and the omitted centralizer sources were recovered from the original folder, with byte counts and hashes, and are included. The second modular chain was not rerun.

Other stored computations are supplied for reproducibility without a claim that every one was independently rerun in versions 03 or 04. The numerical Heisenberg tests are auxiliary controls and are not used as proofs.


## Preparation of v05

The changes are editorial. The manuscript was compiled and inspected, its
references checked, and all retained computational files compared by
SHA-256 with the v04 inventory. Neither the complete F4 calculation nor
the modular Zhu calculations nor the independent mode comparisons were
rerun for this editorial revision.


## Preparation of v06 and v07

Version 06 restored the compact table connecting numbered statements to
computational scripts. Version 07 simplified two entries in that table and
replaced the draft data-availability wording by a release-ready reference to
the intended repository URL. The underlying sources and recorded exact/modular
results were preserved unchanged. Compilation and snapshot integrity checks
are preparation checks, not a rerun of the full F4 or Zhu computations.
