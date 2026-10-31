# Reproducing the computations

## Inputs and environment

Use the v07 snapshot and the two unchanged v03 F4 coordinate ZIPs. Place
the coordinate ZIPs in a sibling folder named release-assets, or adapt the paths.
Verify the snapshot before making changes:

~~~sh
python3 verify_manifest.py
~~~

The recorded v03 F4 checks used Python 3.12.14. Use Python 3.12 and the pinned
exact-arithmetic packages:

~~~sh
python3.12 -m venv .venv
. .venv/bin/activate
python -m pip install -r checks/certificates/deliverable/requirements-exact.txt
~~~

On Windows, activate .venv\Scripts\Activate.ps1 in PowerShell instead.
Do not use Python's -O flag or PYTHONOPTIMIZE. Other historical computations
also use SciPy or Singular/Sage; consult their individual READMEs.

## Complete F4 replay

~~~sh
python checks/certificates/deliverable/run_f4_verification.py \
  --archive ../release-assets/Certificats-poids-quatre-v03-partie-1.zip \
    ../release-assets/Certificats-poids-quatre-v03-partie-2.zip \
  --output ../f4-replay-v07
~~~

The output directory must not exist. The launcher checks source and archive
hashes, extracts the coordinates, verifies every exact membership, checks all
zero products, and requires exhaustive aggregation. The final receipt is
../f4-replay-v07/COMPLETE-F4.json. The coordinate archives expand to about
3.15 GB, in addition to the repository and verification outputs.

Regeneration from lattice mode formulas, without archived weight-four blocks:

~~~sh
python checks/certificates/deliverable/run_f4_verification.py \
  --regenerate --output ../f4-regenerated-v07
~~~

This still uses the archived weight-three basis. Its independent regeneration
is documented by rerun_w3.py and certify_w3_validated.json in the same folder.
Full regeneration can take more than an hour.

## Other checks

Commands, inputs and limits are described in:

- checks/foundations/README-foundations-v03.txt
- checks/generators-isotropy/README.txt
- checks/identification/README-identification.txt
- checks/zhu/README-zhu-v03.txt
- the READMEs under calculs-sources.

A small independent rational-matrix check of the C4 witnesses is:

~~~sh
python checks/identification/recheck_c4_matrices.py \
  checks/identification/c4/c4-exact-witnesses.json
~~~

This checks the supplied matrices, spectra, multiplicities and independence;
it does not reconstruct VOA products. Historical outputs do not imply that
every original script was rerun during v03 or v04. See STATUS.md.

## Independent mode comparisons added to v04

See [the standalone replay instructions](checks/independent-modes/README.md).
They use Python 3.12, gmpy2 and SymPy, and compare 24643 basis-state products
and 310 actual-generator products with a second evaluator. The generator
replay also tests the guard against calls outside the audited weight window.
The original engine is retained byte for byte. These checks supplement the
full F4 certificate; they do not regenerate it by a second implementation.

The prepared v04 package was checked for file integrity, both added replay
commands were run from the packaged layout, and the manuscript was compiled.
The complete F4 and modular Zhu calculations were not rerun for v04.

## Manuscript

The LaTeX file is self-contained, with an embedded bibliography. With TeX Live
or MiKTeX and the requested packages installed:

~~~sh
pdflatex -interaction=nonstopmode -halt-on-error Haagerup-c8-article-I-v07-arxiv.tex
pdflatex -interaction=nonstopmode -halt-on-error Haagerup-c8-article-I-v07-arxiv.tex
~~~

Compilation changes generated files. Verify the original snapshot first, or
compile in a separate copy.


## Source paths and external inputs

[COMPUTATIONS.md](COMPUTATIONS.md) maps the numbered statements, coordinate
conventions and data filenames to the retained sources. The portable F4
launcher and the independent-mode replay commands are the supported entry
points described above. Some inherited exploratory scripts retain original
absolute paths or refer to external inputs; they are not all standalone
commands runnable unchanged from the repository root.

In particular, the original Riesen transport script imports defining_states.json
and definition.json from the external certificate package of Riesen. Those
external files are not redistributed here. The computed transported states
and the identification outputs are included under calculs-sources/cas-I.
See the index for their exact paths.

Version 07 preserves the underlying code and computation data byte for byte
relative to the v05 public-snapshot inventory. The v07 PDF was compiled and
package-integrity checks were performed. Full F4 and modular Zhu calculations
were not rerun for v07. The recorded v04 mode replays and earlier F4
verification remain the relevant computation receipts.
