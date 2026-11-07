# Suppléments de II : symétries de Haagerup et réseaux conformes à c = 8

Cette livraison accompagne l'article de Cédric Arhancet, **Haagerup symmetries of the E8 theory and conformal nets at central charge eight**, dans la version [Hexagon 2611.00034v2](https://hexagonmath.org/2611.00034v2), source **II v14**. La v15 jointe corrige uniquement la déclaration de disponibilité des données de la copie locale v14 et y ajoute les liens de reproduction ; elle n'est pas présentée comme déjà diffusée sur Hexagon ou arXiv.

[Release figée](https://github.com/cedricarhancet/haagerup-c8/releases/tag/article-II-v14-supplements-v01) · [Archive complète de II](https://github.com/cedricarhancet/haagerup-c8/releases/download/article-II-v14-supplements-v01/Haagerup-c8-article-II-v14-supplements-v01.zip) · [État des contrôles](STATUS.md) · [Énoncés et scripts](COMPUTATIONS.md) · [Instructions de reproduction](REPRODUCIBILITY.md) · [Concordance](VERSION-CONCORDANCE.json).

## Démarrage

Depuis ce dossier `II`, avec Python 3.12, sans `-O` ni `PYTHONOPTIMIZE` :

```sh
python verify_manifest.py
python -m pip install -r requirements-exact.txt
python reproduce.py --group certificate --output ../../II-certificate
python reproduce.py --group modes --output ../../II-modes
```

Le premier contrôle utilise seulement la bibliothèque standard. Le groupe `certificate` vérifie les identités catégoriques et le certificat algébrique à partir des tenseurs fournis. Le groupe `modes` recalcule indépendamment les modes de courants du réseau, le noyau, les produits, la forme et les invariants de traces nécessaires à ce certificat. La commande `originals` rejoue les deux implémentations initiales dans une copie de travail ; voir [REPRODUCIBILITY](REPRODUCIBILITY.md).

Les sorties sont écrites hors de la livraison. Le lanceur refuse un dossier de sortie déjà utilisé et conserve les données originales. Le sous-dossier `Calculs-HU-v08` garde les noms historiques nécessaires aux imports ; son contenu public est défini par [ALLOWLIST-COMPUTATIONS.json](ALLOWLIST-COMPUTATIONS.json), et non par une copie intégrale de l'ancien ZIP.

## Correspondance mathématique

II §7.5, lemme 7.14 : caractère de dimension 496, énumération des plongements et quatre triplets commutants. Lemme 7.15 : espace de poids nul de dimension 644, noyau de dimension 16, système de traces de rang 8, 25 quadriques et quotient réduit de dimension 4. Proposition 7.16 : caractère positif et involution PCT fixant les quatre candidats. Le théorème 7.17 utilise ensuite les arguments théoriques de coset pour la complète unitarité. Le théorème général 8.6 garde sa preuve propre ; le calcul fini n'est pas à lui seul une preuve de ce théorème.

Dans [III v10, §9 et §10](https://hexagonmath.org/2611.00164v1), les mêmes calculs sont réutilisés dans la voie algébrique secondaire et le bloc C3, après fourniture des entrées (F2). La preuve analytique principale de `H_J` dans III n'utilise pas ce calcul PCT. Le bloc C3 complet de III comprend des pilotes et des contrôles de tables propres à III ; ils ne sont pas revendiqués comme distribués ou rejoués ici.

## Entrées et provenance

Les données T12 et les poids conformes viennent de I, fixé au tag [article-I-v14-supplements-v01](https://github.com/cedricarhancet/haagerup-c8/releases/tag/article-I-v14-supplements-v01), commit `813df07d2cf2a2b433d13093ccb920031fee5f9a`. La source publique I v14 est incluse dans `inputs-I` avec son SHA-256 afin de permettre la réextraction de la table. Les certificats lourds de poids quatre de I ne sont pas nécessaires au replay local de II ; la démonstration de I reste une entrée importée. Voir [DEPENDENCIES-I.json](DEPENDENCIES-I.json).

Les anciennes sorties `.out` conservées près des scripts datent du corpus du 2 novembre 2026. Les nouvelles exécutions sont dans `checks/release-II-v14`. Les anciens README et manifestes qui portaient sur d'autres arborescences ne sont pas distribués comme instructions actuelles. Une unique extraction du crochet de Lie évite de livrer un ancien script de recherche abandonné ; les deux moteurs VOA, les programmes indépendants et le certificat sont conservés.

Il s'agit de calculs exacts inspectables, avec des tests ciblés des moteurs, et non d'une vérification formelle de la théorie des VOA. Les preuves de transport de II §5, de reconnaissance de II §6 et les arguments analytiques de III restent des preuves mathématiques à lire indépendamment des certificats finis.

[Citation](CITATION.cff) · [Auteur](https://github.com/cedricarhancet). Les attributions du corpus sont conservées ; le statut CC BY-SA 4.0 du manuscrit Hexagon n'est pas une nouvelle licence logicielle sur l'ensemble des fichiers hérités.
