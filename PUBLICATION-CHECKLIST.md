# Publication checklist — v07

This package prepares a GitHub release; it does not itself publish any file.

1. Create or verify the public repository:
   https://github.com/cedricarhancet/haagerup-c8
2. Upload the **contents** of the `haagerup-c8/` directory to the repository
   root. Do not publish the old `v05` PDF or LaTeX as the active manuscript.
3. Run `python3 verify_manifest.py` in the repository root.
4. Review STATUS.md and the conditional existence hypotheses in the paper.
5. Commit the snapshot and tag it `v07`; create GitHub Release `v07`.
6. Attach the following release assets:
   - `Haagerup-c8-v07-reproducibility.zip`
   - `Haagerup-c8-article-I-v07-arxiv.pdf`
   - `Haagerup-c8-article-I-v07-arxiv.tex`
   - `Certificats-poids-quatre-v03-partie-1.zip` (420645698 bytes)
   - `Certificats-poids-quatre-v03-partie-2.zip` (99768703 bytes)
   - `SHA256SUMS-release-v07.txt`
7. Check that the URLs work for an unauthenticated reader and that each
   archive's SHA-256 matches RELEASE-ASSETS.md.
8. Only once the repository **and both certificate archives** are public,
   submit the v07 manuscript to arXiv. Its Data availability paragraph uses
   present-tense availability and is **not accurate before that moment**.
9. If Zenodo assigns a DOI later, insert its *actual* DOI into the article
   and rebuild the PDF/package; do not invent a DOI.

The mathematical content and completeness claims should remain conditional
as stated in the article; package preparation is not a new independent audit.
