# Required release assets — v07

The two full weight-four coordinate archives are too large to put into the
repository tree. Upload **both**, unchanged, as assets to GitHub Release v07.
They retain `v03` in their names because the coordinate data are unchanged.

| File | Bytes | SHA-256 |
| --- | ---: | --- |
| `Certificats-poids-quatre-v03-partie-1.zip` | 420645698 | `8ea5903408d606eaaf3f98092293c2e5c7a9701bcbd7c11b645dcf3ac485ac17` |
| `Certificats-poids-quatre-v03-partie-2.zip` | 99768703 | `fa121f49c768943bbc5af1a32b22402c6bdb3cee03480f79f8ebce3dfb34e1d2` |

The verified archive list and coordinate coverage are in
[`checks/certificates/deliverable/ARCHIVES-F4.json`](checks/certificates/deliverable/ARCHIVES-F4.json).
The complete F4 replay requires **both** archives. Their recorded hashes
should be checked against the exact release uploads. Do not claim full
replayability if either asset is missing.

The accompanying reproducibility snapshot can be attached as
`Haagerup-c8-v07-reproducibility.zip`. SHA256SUMS-release-v07.txt lists the
release-file digests. No Zenodo DOI has been assigned in this preparation.
