# Vérifications achevées — livraison II v14, v02

Le rejeu supplémentaire commencé le 2026-11-07T11:04:06.246919+00:00 et terminé le 2026-11-07T11:52:19.946201+00:00 a achevé ses **15 commandes**, avec codes de sortie nuls et validations des résultats. Il complète les 43 commandes de la livraison v01. Les moteurs mathématiques sont inchangés. Les reçus et journaux nouveaux sont dans `checks/release-II-v14-v02/extended` ; les reçus de `checks/release-II-v14` décrivent l'étape v01 et gardent leur portée historique.

## Fermeture en poids deux dans le premier modèle

| Candidat | Dimension S2 | Paires contrôlées | Périmètre | Échecs de fermeture |
| --- | ---: | ---: | --- | ---: |
| 0 | 120 | 1200 | Exhaustif pour ce jeu de générateurs et la base | 0 |
| 1 | 120 | 180 | Échantillonnage annoncé de la base | 0 |
| 2 | 120 | 180 | Échantillonnage annoncé de la base | 0 |
| 3 | 120 | 180 | Échantillonnage annoncé de la base | 0 |

Chaque paire est testée aux modes 1, 2 et 3. Pour le candidat 0, les dix vecteurs choisis sont confrontés aux 120 vecteurs de base : les **1 200 paires** du calcul précédemment interrompu sont entièrement parcourues, soit 3 600 contrôles de modes. Pour les autres candidats, les 180 paires correspondent au parcours annoncé de `SP.B[::7]`. Ce ne sont pas des tests de tous les modes en tous poids.

Les quatre espaces S2 et S2full sont régénérés dans des copies isolées, puis les quatre contrôles PCT sont exécutés sur ces données nouvelles. Pour chaque candidat : dimension 120, caractère U2 attendu, zéro échec de commutant, de stabilité PCT, d'automorphisme échantillonné et de positivité échantillonnée. Les JSON régénérés correspondent aux références ; les pickles R3 sont identiques en octets. Les comparaisons sont dans `ARTIFACTS.json` et les reçus de chaque candidat.

## Reconstruction dans le second modèle

`s2_space.py 0 all` reconstruit l'espace entier de poids deux. `s2_check.py 0` confirme sa dimension 120, ses multiplicités de poids et sa stabilité par le PCT dans tous les blocs. La matrice hermitienne de la tranche de dimension cinq est positive : ses cinq mineurs principaux valent exactement 50/507. Le résultat est confronté au journal historique, en ignorant seulement le temps de calcul.

`twist_check.py` retrouve les 16 choix rationnels et les 48 autres choix sur Q(sqrt(-6)), préserve les générateurs de courants et contrôle 2659 coefficients des produits échantillonnés, avec zéro échec. L'échantillonnage fixé n'est pas une preuve de fermeture en tous poids.

Le pickle de l'espace secondaire est fourni sous forme gzip dans `checks/release-II-v14-v02/generated/R3-verif/`, avec les empreintes compressée et décompressée dans `ARTIFACTS.json`. Il est généré par ce nouveau calcul et comparé à la référence historique.

## Statut de l'ancien arrêt et limites

Le lancement arrêté à 423,375 secondes dans v01 a été remplacé par une nouvelle exécution achevée. Son ancien journal partiel est conservé comme historique, avec un renvoi à cette livraison. Aucun diagnostic annoncé pour la présente passe ne reste interrompu ou en attente.

Les contrôles essentiels de v01 — certificat portable, reconstruction indépendante des tenseurs, noyaux, équations et PCT dans les deux moteurs — conservent leurs reçus. Le présent travail achève les diagnostics S2 et de forme rationnelle qui n'avaient pas été terminés. La localisation numérique ancienne des candidats n'est pas relancée : la liste est certifiée par les équations exactes et le calcul d'idéaux. Les contrôles arithmétiques bornés de §6 gardent leur portée finie. Les fondations de I, les certificats de Zhu et le bloc C3 complet de III sont des entrées ou livraisons distinctes. Ces programmes ne constituent pas une vérification formelle générale des axiomes VOA ou des preuves analytiques.

Les conditions de commutant et le PCT utilisés dans la preuve de II, proposition 7.16, sont distingués des diagnostics supplémentaires de fermeture. La distinction de rôle demeure, mais le rejeu supplémentaire est maintenant achevé.
