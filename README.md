# Haagerup à charge centrale huit : articles et suppléments

Suppléments publics des articles de [Cédric Arhancet](https://github.com/cedricarhancet). Chaque livraison fige ses sources, ses données, ses limites de validation et ses empreintes ; les anciennes releases restent accessibles.

| Article | Manuscrit de référence | Suppléments publics |
| --- | --- | --- |
| I — Haagerup vertex operator algebras at central charge eight | [Hexagon 2610.00733v2](https://hexagonmath.org/2610.00733v2), I v14 | [Release article-I-v14-supplements-v01](https://github.com/cedricarhancet/haagerup-c8/releases/tag/article-I-v14-supplements-v01) ; instructions ci-dessous |
| II — Haagerup symmetries of the E8 theory and conformal nets at central charge eight | [Hexagon 2611.00034v2](https://hexagonmath.org/2611.00034v2), II v14 ; [v16 documentaire](II/manuscript/Haagerup-reseaux-conformes-c8-v16-arxiv.pdf) | [Dossier II](II/README.md), [release article-II-v14-supplements-v02](https://github.com/cedricarhancet/haagerup-c8/releases/tag/article-II-v14-supplements-v02) |
| III — The Haagerup center at central charge eight via conformal nets | [Hexagon 2611.00164v1](https://hexagonmath.org/2611.00164v1), III v10 | Livraison propre à III à préparer ; les calculs de II importés dans sa voie secondaire sont fournis dans II |

I fournit notamment la construction, la structure unitaire et la C2-cofinitude utilisées parmi les fondations de III. Les outils généraux de transport et de reconnaissance de II sont établis indépendamment de H_J. III démontre H_J par sa preuve analytique ; II §7.5 fournit aussi une voie calculatoire vers une structure unitaire compatible sur la diagonale. Les calculs finis et les preuves analytiques gardent leur portée respective.

Depuis la racine du dépôt, `python verify_manifest.py` contrôle l'instantané commun. Depuis `II`, suivre son propre [README](II/README.md) ; les commandes ci-dessous restent celles de I. Pour citer la version figée de I et son manifeste d'origine, utiliser sa release et son tag propres.

## Article I



Cette livraison accompagne **Haagerup vertex operator algebras at central charge eight**, de Cédric Arhancet, dans la version [Hexagon 2610.00733v2](https://hexagonmath.org/2610.00733v2). La source diffusée est identique, octet par octet, à `Haagerup-c8-article-I-v14-arxiv.tex`.

- [Livraison figée `article-I-v14-supplements-v01`](https://github.com/cedricarhancet/haagerup-c8/releases/tag/article-I-v14-supplements-v01) ; [archive de reproduction](https://github.com/cedricarhancet/haagerup-c8/releases/download/article-I-v14-supplements-v01/Haagerup-c8-article-I-v14-supplements-v01.zip).
- [État des vérifications](STATUS.md), [énoncés et calculs](COMPUTATIONS.md), [instructions de reproduction](REPRODUCIBILITY.md).
- [Concordance des versions I–III](VERSION-CONCORDANCE.md), [empreintes des fichiers](MANIFEST-PUBLIC.json), [liste explicite des fichiers de calcul](ALLOWLIST-COMPUTATIONS.json).

## Démarrage

Utiliser Python 3.12, sans `-O` et sans `PYTHONOPTIMIZE`. Dans une copie propre :

```sh
python verify_manifest.py
python -m pip install -r checks/certificates/deliverable/requirements-exact.txt
python reproduce.py --group smoke --output ../smoke-I14
python reproduce.py --group quick --output ../quick-I14
```

`smoke` utilise seulement la bibliothèque standard. `quick` relance aussi les contrôles de fondations, de produits de champs de réseau et de la boule du réseau ; sa durée dépend de la machine. Les deux archives F4 sont nécessaires à la vérification exhaustive de poids quatre, décrite dans [REPRODUCIBILITY](REPRODUCIBILITY.md). Elles représentent environ 520 Mo compressés et 3,15 Go décompressés.

## Portée mathématique

I construit le modèle explicite `S`, sa structure unitaire (théorème 10.10), vérifie le certificat de poids quatre (proposition 10.17), puis démontre la génération forte et la `C_2`-cofinitude (théorème 10.21) et l'identification à l'algèbre de Jiao (théorème 10.24). Ces fondations sont établies dans I indépendamment de l'hypothèse `H_J`.

[III, Hexagon 2611.00164v1](https://hexagonmath.org/2611.00164v1), théorème 6.1, établit `H_J`, la complète unitarité et la comparaison VOA–réseau à partir des entrées algébriques de son annexe et des outils généraux de [II, Hexagon 2611.00034v2](https://hexagonmath.org/2611.00034v2). I apporte notamment à III l'entrée unitaire (I1), ainsi que la `C_2`-cofinitude utilisée pour ses fondations algébriques. Le certificat **(F4) de I**, défini en 10.16, est distinct de l'entrée catégorique **(F4) de III**.

Cette livraison est celle de I. Elle ne distribue ni ne prétend revérifier les certificats propres à II et à l'annexe de III. Le calcul modulaire de Zhu conservé dans I est un contrôle fini ; le passage en caractéristique zéro et les douze modules sont traités dans III.

## Archives et historique

Les deux volumes F4 sont repris exactement de [la release v07](https://github.com/cedricarhancet/haagerup-c8/releases/tag/v07). Leurs liens stables et SHA-256 sont dans [EXTERNAL-ASSETS.json](EXTERNAL-ASSETS.json). Les anciennes releases et leurs liens sont conservés. Les reçus historiques gardent leur date et leur portée ; les nouvelles exécutions figurent dans `checks/release-I-v14/`.

La sélection contient les calculs nécessaires aux énoncés publiés, leurs données, leurs témoins et leurs résultats. Aucun carnet de recherche privé n'est ajouté. Certains modules nécessaires gardent leurs anciens noms pour préserver les imports ; leur usage est précisé dans l'index.

## Citation et attribution

Citer le manuscrit et cette livraison séparément : [CITATION.cff](CITATION.cff), [CITATION.bib](CITATION.bib). [Profil de l'auteur](https://github.com/cedricarhancet).

Les licences et attributions existantes sont conservées. Le manuscrit Hexagon est sous CC BY-SA 4.0 ; cette mention ne constitue pas l'octroi d'une nouvelle licence logicielle sur tous les fichiers hérités.
