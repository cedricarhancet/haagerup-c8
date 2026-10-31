# Haagerup vertex operator algebras at central charge eight

Manuscript and reproducibility materials for the paper by Cédric Arhancet.

- [Manuscript PDF](Haagerup-c8-article-I-v07-arxiv.pdf)
- [LaTeX source](Haagerup-c8-article-I-v07-arxiv.tex)
- [Reproduction instructions](REPRODUCIBILITY.md)
- [Mathematical scope and verification status](STATUS.md)
- [Statement and coordinate index](COMPUTATIONS.md)
- [Verification provenance](VERIFICATION-RECORD.md)
- [Changes](CHANGELOG.md)

The classification is conditional on the external existence hypotheses
(H_J) or (H_R) stated in the paper. The central-charge-eight results do not
construct the proposed central-charge-two scaling limit of the Haagerup chain.

## Files

The repository contains the computational sources, exact and modular inputs,
verification scripts, recorded outputs and source inventories.

| Directory | Contents |
| --- | --- |
| calculs-sources | Scientific code and data used in the proofs |
| verification-foundations, verification-certificates, verification-weight-four | Verifiers for the foundations, the certificates and the weight-four computation |
| checks | Additional checks and their receipts |
| checks/independent-modes | Independent comparison of vertex-operator products |
| checks/certificates/deliverable | Complete verification instructions and receipts for condition (F4) |
| checks/identification/c4 | Exact commutant coordinates and matrix witnesses |
| provenance | File inventories |

The two complete weight-four coordinate archives are release assets:

1. Certificats-poids-quatre-v03-partie-1.zip
2. Certificats-poids-quatre-v03-partie-2.zip

Both are required to replay the verification of condition (F4). Their hashes
are in [ARCHIVES-F4.json](checks/certificates/deliverable/ARCHIVES-F4.json).
The release asset Haagerup-c8-v07-reproducibility.zip is a snapshot of this
repository and does not contain these two archives.

## Verification

~~~sh
python3 verify_manifest.py
~~~

This checks file integrity. Mathematical verification is described in
[REPRODUCIBILITY.md](REPRODUCIBILITY.md). Some scripts overwrite their
outputs; run them in a copy of the repository.

The manuscript contains a table relating scripts to numbered statements;
COMPUTATIONS.md gives full paths and operational details.

The verification of (F4) passed on the complete coordinate data. The
independent product comparisons agree exactly but do not regenerate every
(F4) coordinate by a second method. The characteristic-zero Zhu algebra, the
exhaustive branching and an independent categorical identification remain
open; see [STATUS.md](STATUS.md).

## Citation and reuse

Use [CITATION.cff](CITATION.cff) and specify the release or commit used.
No DOI has been assigned. Existing source attributions and notices remain in
force; no new reuse licence is granted here.
