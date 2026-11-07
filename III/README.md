# Article III — suppléments publics

Livraison [article-III-v10-supplements-v01](https://github.com/cedricarhancet/haagerup-c8/releases/tag/article-III-v10-supplements-v01) de **The Haagerup center at central charge eight via conformal nets**, [Cédric Arhancet](https://github.com/cedricarhancet). Référence scientifique : [Hexagon 2611.00164v1](https://hexagonmath.org/2611.00164v1), source v10. La [révision documentaire v11](manuscript/Haagerup-c8-article-3-v11.pdf) précise les archives et les liens sans changer les énoncés.

| Livraison | Contenu complet et vérifié |
| --- | --- |
| [A1/A4/C3](https://github.com/cedricarhancet/haagerup-c8/releases/download/article-III-v10-supplements-v01/Calculs-III-A1-A4-C3-v11.zip) | Identités de réseau, contrôles DJX, deux moteurs C3, tables, certificat portable, données et journaux. |
| [Certificats Zhu](https://github.com/cedricarhancet/haagerup-c8/releases/download/article-III-v10-supplements-v01/Certificats-Zhu-III-v07.zip) | Données de production, vérificateurs des deux inclusions d'idéaux, confluence, fibre spéciale, poids, branchements et lanceur de régénération. |
| [Source v11](https://github.com/cedricarhancet/haagerup-c8/releases/download/article-III-v10-supplements-v01/Haagerup-c8-article-3-v11.tex) / [PDF v11](https://github.com/cedricarhancet/haagerup-c8/releases/download/article-III-v10-supplements-v01/Haagerup-c8-article-3-v11.pdf) | Manuscrit avec références stabilisées et disponibilité des données synchronisée. |
| [Paquet arXiv v11](https://github.com/cedricarhancet/haagerup-c8/releases/download/article-III-v10-supplements-v01/Depot-arXiv-III-v11.zip) | Source v11 à la racine ; les deux archives ci-dessus, identiques octet par octet, dans `anc/`. Préparation du dépôt, sans soumission automatique. |

## Démarrage depuis une copie propre

```sh
python -I -X utf8 -B verify_manifest.py
python -m pip install -r requirements-exact.txt
python -I -X utf8 -B reproduce.py --group smoke --output ../../III-smoke
python -I -X utf8 -B reproduce.py --group C3 --output ../../III-C3
```

Le lanceur télécharge les archives et, pour C3, les sources TeX de I et III depuis des releases identifiées, puis contrôle leurs tailles, SHA-256 et manifestes internes avant utilisation. Il écrit dans un dossier neuf extérieur au dépôt. Les autres groupes sont décrits dans [REPRODUCIBILITY.md](REPRODUCIBILITY.md).

`code/` est une copie consultable des seuls scripts des deux archives ; les données binaires complètes sont dans les pièces jointes de release. Exécuter les calculs depuis les archives extraites par le lanceur, qui fournit les jeux de données attendus. `CODE-INDEX.json` donne la correspondance exacte entre copies et archives. Les nouveaux reçus sont aussi consultables dans `checks/`.

I fournit le modèle de réseau, sa structure unitaire et les fondations de génération forte et de `C_2`-cofinitude. Les grands certificats de poids quatre de I sont ceux de sa [release propre](https://github.com/cedricarhancet/haagerup-c8/releases/tag/article-I-v14-supplements-v01), avec les deux volumes conservés dans [v07](https://github.com/cedricarhancet/haagerup-c8/releases/tag/v07). Ils ne sont pas remplacés par C3 et ne sont pas dupliqués ici.

Les outils de transport et de reconnaissance de [II](https://github.com/cedricarhancet/haagerup-c8/releases/tag/article-II-v14-supplements-v03) sont établis indépendamment de `H_J`. III démontre `H_J` au théorème 6.1 par la diagonale locale et la reconstruction. A1/A4 accompagnent les identifications et les comparaisons ultérieures ; C3 accompagne la voie calculatoire secondaire du §9 et réutilise II §7.5. L'appendice A fournit les entrées algébriques de III. Une vérification des calculs finis ne constitue pas une expertise indépendante des arguments analytiques de reconstruction et de restriction du PCT.


Voir [COMPUTATIONS.md](COMPUTATIONS.md), [STATUS.md](STATUS.md), [VERSION-CONCORDANCE.md](VERSION-CONCORDANCE.md), les dépendances (`DEPENDENCIES.json`), les téléchargements (`ASSETS.json`) et le manifeste (`MANIFEST-III.json`). [CITATION.cff](CITATION.cff) et [CITATION.bib](CITATION.bib) distinguent le manuscrit de cette livraison de calculs. Les empreintes des pièces jointes sont également dans `SHA256SUMS-article-III-v10-supplements-v01.txt` sur la release. Les fichiers historiques et anciennes releases sont conservés.

Les licences et attributions existantes sont conservées. Cette livraison n'ajoute pas de licence logicielle globale aux fichiers hérités. Le statut CC BY-SA 4.0 du manuscrit Hexagon est distinct de celui de chaque composant logiciel.
