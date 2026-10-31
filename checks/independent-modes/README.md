# Independent low-weight mode comparisons

These checks supplement, rather than replace, the complete F4 certificate.
The generator inputs are inherited from the original setup. Product evaluation
uses a second implementation based on simultaneous normal ordering.

From the repository root, with Python 3.12:

~~~sh
python -m pip install -r checks/independent-modes/requirements.txt
python checks/independent-modes/independent_modes.py \
  --engine checks/independent-modes/reference/lattice_voa.py \
  --output ../independent-modes-replay.json
python checks/independent-modes/replay_generators.py \
  --engine checks/independent-modes/reference/lattice_voa.py \
  --inputs checks/independent-modes/generator-inputs.json \
  --record checks/independent-modes/generator-modes-results.json \
  --output ../generator-modes-replay.json
~~~

The first command compares 24643 products. It also confirms the known
out-of-window counterexample of output weight nine; that mismatch is expected
and is not counted as an agreement. The second command checks the same 310
generator products with the weight guard enabled and confirms rejection of
the weight-nine call.

Re-extracting the generator inputs is optional. It needs the original
generation-forte files and their additional dependencies:

~~~sh
python checks/independent-modes/generator_modes.py \
  --source calculs-sources/generation-forte \
  --output ../generator-modes-regenerated.json
~~~

The original receipt records 16 coefficient-field products and three
Virasoro identities in addition to the 310 comparisons. The portable replay
does not claim to repeat the generator extraction.
