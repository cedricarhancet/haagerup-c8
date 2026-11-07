# Suppléments de II — livraison v03, diagnostics achevés

La v03 corrige l'encodage de deux descriptions de provenance dans VERSION-CONCORDANCE. Les calculs achevés de v02, leurs résultats et le manuscrit v16 sont conservés sans changement. Le lien de disponibilité du PDF v16 cite v02, qui reste accessible et fournit les mêmes données mathématiques. Voir RELEASE-v03.md.

Cette livraison accompagne **Haagerup symmetries of the E8 theory and conformal nets at central charge eight**, de Cédric Arhancet, dans la version [Hexagon 2611.00034v2](https://hexagonmath.org/2611.00034v2), source II v14. La v16 documentaire jointe actualise le lien de reproduction ; elle ne modifie aucun énoncé mathématique de la copie locale précédente et n'est pas présentée comme déjà diffusée sur Hexagon ou arXiv.

[Release figée](https://github.com/cedricarhancet/haagerup-c8/releases/tag/article-II-v14-supplements-v03) · [Archive autonome](https://github.com/cedricarhancet/haagerup-c8/releases/download/article-II-v14-supplements-v03/Haagerup-c8-article-II-v14-supplements-v03.zip) · [État détaillé](STATUS.md) · [Énoncés et programmes](COMPUTATIONS.md) · [Reproduction](REPRODUCIBILITY.md) · [Concordance](VERSION-CONCORDANCE.json).

La v02 termine le diagnostic de fermeture S2 interrompu dans [v01](https://github.com/cedricarhancet/haagerup-c8/releases/tag/article-II-v14-supplements-v01) et rejoue les autres diagnostics annoncés : espaces S2 des quatre candidats, PCT sur les espaces régénérés, reconstruction secondaire complète en poids deux, positivité et forme rationnelle. Quinze nouvelles commandes s'ajoutent aux 43 commandes vérifiées de v01. Le candidat 0 est contrôlé sur 1 200 paires ; les trois autres sur les 180 paires de leur parcours annoncé. Le détail des périmètres est dans STATUS.

## Démarrage

Python 3.12, sans optimisation des assertions. Depuis `II/` :

```sh
python verify_manifest.py
python verify_completed_diagnostics.py
python -m pip install -r requirements-exact.txt
python reproduce.py --group certificate --output ../../II-certificate-v02
python reproduce.py --group modes --output ../../II-modes-v02
python reproduce.py --group originals --output ../../II-originals-v02
python reproduce.py --group auxiliary --output ../../II-auxiliary-v02
python reproduce.py --group extended --output ../../II-extended-v02
```

Les sorties doivent rester hors du dépôt. Le groupe `all` désigne les quatre groupes du parcours principal ; `extended` achève les diagnostics supplémentaires. Pour maîtriser la mémoire, le lanceur `reproduce.py` exécute ce dernier avec une tâche à la fois. Le pilote direct `extended_replay.py` accepte `--workers 1` ou `--workers 2`. Les calculs lourds doivent aller jusqu'au bilan final ; les résultats sont contrôlés, au-delà du seul code de sortie.

## Entrées et rôle dans les articles

Les tables T12 et les hauteurs viennent de I, fixé au tag [article-I-v14-supplements-v01](https://github.com/cedricarhancet/haagerup-c8/releases/tag/article-I-v14-supplements-v01), commit `813df07d2cf2a2b433d13093ccb920031fee5f9a`. La source I v14 nécessaire à l'extraction est incluse dans `inputs-I`. II §7.5 fournit les calculs des lemmes 7.14–7.15 et de la proposition 7.16 ; les théorèmes de complète unitarité utilisent ensuite les arguments théoriques du manuscrit. [III v10](https://hexagonmath.org/2611.00164v1), §9 et bloc C3, importe ces calculs dans sa voie secondaire. Sa preuve analytique principale de H_J utilise la reconstruction et sa structure unitaire compatible.

La sélection de calculs est explicite dans ALLOWLIST-COMPUTATIONS. Les sources mathématiques sont inchangées depuis v01 ; les nouveaux pilotes, reçus et données régénérées sont identifiés dans SOURCE-PROVENANCE. Les `.out` historiques proches des scripts restent distincts des journaux nouveaux. Les anciennes releases et leurs liens sont conservés. Les grands certificats de I et les certificats propres à III ne sont pas dupliqués dans ce paquet.

[Citation](CITATION.cff) · [Auteur](https://github.com/cedricarhancet). Les attributions et statuts de licence antérieurs sont conservés.
