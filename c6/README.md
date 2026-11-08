# Haagerup symmetry in a rank-six lattice conformal field theory — suppléments publics

Suppléments de calcul de l'article **Haagerup symmetry in a rank-six lattice conformal field theory**, de [Cédric Arhancet](https://github.com/cedricarhancet). Le manuscrit de référence est la source v12, déposée sur Hexagon le 8 novembre 2026 (identifiant à ajouter à la mise en ligne). La source et le PDF sont dans [`manuscript/`](manuscript/).

Le paragraphe *Data availability* de l'article cite `lattice/verify_c6_section6_v01.py`. Ce fichier est diffusé ici à ce chemin, sans modification.

## Démarrage

Utiliser Python 3.12 ou 3.13. Dans une copie propre :

```sh
python -m pip install -r requirements.txt
python -I -X utf8 -B run_all.py --output ../../c6-run
```

Le lanceur copie les scripts dans un dossier neuf, extérieur au dépôt, les exécute, puis compare chaque sortie texte et chaque certificat JSON aux versions enregistrées ici. L'ensemble dure quelques secondes. `SHA256SUMS.txt` donne les empreintes des fichiers de ce dossier.

## Contenu

| Fichier | Contrôles en arithmétique exacte | Énoncés de la v12 |
| --- | --- | --- |
| `lattice/verify_c6_section6_v01.py` | Bases de $K$ dans $E_8^\pm$, réseau de Gram de $L$, $\det K = 39$, $q(g_K) = 17/39$, couches $14, 96, 156$, réseau $K^\flat$ à 30 racines ; réseau $N$ sans racine de déterminant 351, ses quatre droites isotropes et ses quatre surréseaux isométriques à $K$ ; action de Weyl sur les droites ; étiquettes pointées $26m+21n$, énergies minimales et tous les moments minimisants ; les huit étiquettes pertinentes ; triplets cubiques ; poids de contrainte ; objets de twist trivial ; table de génération ; réalisation de niveaux $(1,39)$ | Proposition 6.1, Proposition 6.2, Théorème 6.5, Théorème 6.8, Lemme 6.9, Théorème 6.11, Proposition 7.4, Proposition 7.6, Proposition 7.8, Proposition 7.11 et sa table, Propositions 7.12 et 7.14 |
| `lattice/verify_c6_complements_v12.py` | Certificat de masse : $B_{3,\chi} = 528$, facteur dyadique $1/56$ de Conway–Sloane IV, masse $11/2880 = 1/288 + 1/2880$, avec $\lvert O(K)\rvert = 288$ et $\lvert O(K^\flat)\rvert = 2880$ par PARI ; charge centrale $170/41$ et valeurs propres $23/41$ et $33/41$ ; valeurs $\rho^U = 4/3$ et $-1/3$, $\rho = 16/9$ et $1/9$ ; parités de cocycles ; étiquettes locales du bulk mixte | Proposition 6.2, Section 7.3, Proposition 7.8, Proposition 7.14, Proposition 7.16 |
| `independent/independent_lattice_recheck.py` | Seconde implémentation, en SymPy seul : couches, ordres $288$ et $2880$ des groupes orthogonaux par énumération directe, racines des surréseaux, masse | Proposition 6.2, Théorème 6.8 |
| `independent/independent_complements_recheck.py` | Seconde vérification des bases de $K$ dans les cas (I) et (II) | Proposition 6.1 |
| `independent/verify_weyl_isotropic_permutation.py` | Droites isotropes sur $\mathbb F_3$, permutation de Weyl, isométries $\{1,14,25,38\}$ | Lemme 6.9, Théorème 6.11 |
| `independent/verify_mixed_pointed_bulk.py` | Monodromie et étiquettes pertinentes du bulk mixte | Proposition 7.16 |

Les fichiers JSON de `independent/` sont les certificats produits par ces scripts. Le lanceur les régénère et les compare aux versions enregistrées.

## Portée

Ces scripts vérifient les calculs finis de l'article. Les arguments catégoriques et analytiques (centre plein, action de Weyl sur les algèbres vertex, Corollaire 7.15, structure unitaire, formule de masse elle-même) sont démontrés dans le texte. Ce ne sont pas des certificats informatiques. Une vérification de ces calculs ne constitue pas une expertise indépendante de ces arguments.

L'article prolonge les trois articles à charge centrale huit dont les suppléments figurent à la racine de ce dépôt.

Les licences et attributions existantes sont conservées. Le statut CC BY-SA 4.0 du manuscrit Hexagon est distinct de celui des scripts.
