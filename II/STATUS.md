# Contrôles de la livraison de II — 7 novembre 2026

## Exécutions de cette passe

| Groupe | Commandes terminées | État | Reçu |
| --- | ---: | --- | --- |
| certificate | 2 | PASS | `checks/release-II-v14/certificate/RECEIPT.json` |
| modes | 1 | PASS | `checks/release-II-v14/modes/RECEIPT.json` |
| originals | 29 | PASS | `checks/release-II-v14/originals/RECEIPT.json` |
| auxiliary | 11 | PASS | `checks/release-II-v14/auxiliary/RECEIPT.json` |

Soit 43 commandes achevées, toutes avec code de sortie nul. Les résultats détaillés, les durées, stdout, stderr et les versions des dépendances sont conservés dans les reçus, séparément des anciens `.out` du corpus. Les noms de répertoires locaux de préparation sont remplacés dans les commandes archivées par des repères ; les commandes portables de reproduction figurent dans le lanceur et dans REPRODUCIBILITY.

Le vérificateur algébrique confirme le rang 8 des vingt traces, les 25 quadriques, l'égalité des idéaux sur Q et les quatre points réduits/PCT. La reconstruction indépendante confirme le noyau G de dimension 16, les 136 produits et formes et les invariants de traces des vingt espaces de multiplicité. Les deux moteurs initiaux sont rejoués ; leurs solutions et les tests PCT sont confrontés aux résultats attendus. L'énumération indépendante conserve 640 représentations orthogonales et une paire admissible.

La table a été réextraite de I v14, version figée fournie dans `inputs-I`, et les JSON de tables correspondent. Le certificat a été réexporté depuis les données R3 régénérées et comparé au JSON distribué. Les autres comparaisons de JSON et la comparaison brute des données présentes dans la copie de travail sont enregistrées dans les reçus et DATA-COMPARISON. Cette dernière table ne prétend pas que chaque donnée préexistante a été régénérée : les commandes effectivement exécutées en fixent la portée.

Sur les 34 données comparées dans DATA-COMPARISON, 24 sont identiques en octets. Les dix JSON restants ne diffèrent que par les fins de lignes LF/CRLF de Windows : leur contenu JSON et leurs octets après normalisation des fins de lignes sont identiques. Cela ne transforme pas les données préservées en données régénérées ; les commandes des reçus fixent cette distinction.

## Adaptations de la livraison

Le calcul source est repris du ZIP Calculs-HU-v08 vérifié. Le crochet `br` du fichier de recherche de signes est extrait sans changement dans `R3-verif/lie_bracket.py` ; `solve_signs.py` l'importe directement. Le programme de recherche abandonné n'est pas distribué. Ce couple a été utilisé dans la nouvelle exécution des calculs initiaux. Deux docstrings de `data/` renvoient désormais à la source I v14 et aux instructions livrées ; leur code exécutable est identique par comparaison d'arbres syntaxiques. Les moteurs et le certificat sont inchangés, avec leurs empreintes dans SOURCE-PROVENANCE.

## Diagnostic long séparé du parcours principal

Le calcul `R3/r3a_S2.py 0` est terminé. Le diagnostic exhaustif `R3/r3a_S2_full.py 0 full` a été arrêté volontairement après 423.375 secondes : la dimension 120, le caractère et les conditions de commutant sur les 120 vecteurs étaient vérifiés sans échec, mais la fermeture exhaustive n'était pas terminée. Son arrêt provoque un code de sortie non nul dans le lancement initial ; il n'est pas compté parmi les succès et le journal partiel est conservé dans `checks/release-II-v14/diagnostics`. Les treize commandes déjà terminées ont été conservées, puis les seize commandes essentielles restantes ont été exécutées dans une nouvelle copie de continuation, avec un reçu explicite.

Les autres régénérations S2 et le diagnostic de forme rationnelle ne sont pas relancés dans cette passe ; leurs données et sorties historiques sont distribuées. Les quatre exécutions de `r3b_verify.py` lisent notamment les espaces S2full fournis et en vérifient le PCT ; cela ne constitue pas leur régénération indépendante. Ces diagnostics supplémentaires sont disponibles par le groupe `extended` et ne font pas partie de `all`. Leur achèvement n'est pas utilisé dans le certificat ni dans la preuve par commutants de la proposition 7.16.

## Sorties historiques et limites

Les `.out` proches des scripts proviennent de l'archive d'entrée du 2 novembre 2026. Les anciens README v04/v07/v08, les manifestes qui portaient sur d'autres arborescences, les rapports et les tentatives abandonnées ne sont pas utilisés comme reçus de la présente passe.

Le premier modèle utilise les valeurs de localisation numérique archivées de `idem_complex.json`, puis reconstruit les candidats et vérifie exactement les équations et les produits. La recherche numérique des candidats n'est pas relancée ; leur exhaustivité résulte ici du calcul d'idéaux et du certificat exact. Des tests de VOA et de positivité utilisent des échantillons fixés ; leur réussite ne démontre pas des identités en tous poids. Les contrôles `section6/check_galois.py` et `check_m.py` sont bornés et comportent des calculs numériques ; les preuves générales restent celles du manuscrit.

Les tenseurs reconstruits indépendamment sont confrontés par leurs invariants de traces aux contraintes nécessaires, sans imposer une base de multiplicité commune. La preuve de l'existence de la diagonale, la théorie des VOA, le passage aux réseaux conformes et les hypothèses analytiques de reconstruction ne sont pas formellement vérifiés par ces programmes. Les certificats lourds de I, les certificats de Zhu et le bloc C3 complet de III ne sont pas rejoués dans cette livraison.

## Manuscrit documentaire

La v15 ne change que la déclaration de disponibilité des données de la copie locale v14. Le lanceur LaTeX local a terminé en 3 passes ; aucune référence ou citation indéfinie n'est signalée. Les pages 61–62 modifiées ont été rendues et contrôlées visuellement. La référence de calcul demeure la source II v14 de Hexagon ; les deux commandes typographiques différentes dans la copie locale sont précisées dans VERSION-CONCORDANCE.
