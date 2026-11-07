# Énoncés, entrées, programmes et résultats

Les numéros ci-dessous sont ceux de II v14, Hexagon 2611.00034v2. La v16 documentaire conserve la numérotation. Les chemins de calcul sont relatifs à `Calculs-HU-v08/`. Les commandes exactes, leurs répertoires de lancement et leur ordre sont enregistrés par `reproduce.py` ; le tableau précise leur signification.

| Énoncé ou contrôle | Fichiers / commande élémentaire | Entrées | Résultat attendu et limite |
| --- | --- | --- | --- |
| II §7.5, données T12 | `data/t12_table.py <I-v14.tex>` ; `data/t12_W1.py t12_table.json` | source publique I v14 ; tables et hauteurs du cas (II) | égalité avec les tables JSON livrées, dimensions attendues des douze espaces et U2 = 120 |
| Lemme 7.14, réduction du plongement | `data/t12_e8split.py t12_W1.json` ; `R3/t12b_optionB.py` ; `scripts-rapporteurs/F1/enum_sl2.py` | caractère de W1 et bornes d'indices de la preuve | 640 représentations orthogonales sur 13 dimensions, une paire admissible ; exclusion de l'option B |
| Lemme 7.14, réalisation concrète | `R3/sl2x4.py`, `R3/weights_W.py` ; `R3-verif/embed.py`, `R3-verif/check_W1_char.py` | racines E8+E8, cocycle, données de supports ; T12 | quatre triplets commutants, indices (2,1), (6,7), (1,2), (7,6), caractère de dimension 496 |
| Lemme 7.15, noyau | `R3/r3a_G.py` ; `R3-verif/g_space.py` | modes 0,1,2 des douze courants | espace de poids nul 644, noyau G de dimension 16 ; rang 8 des conditions de traces |
| Lemme 7.15, tenseurs | `R3/r3a_Griess.py` ; `R3-verif/griess.py` | bases de G, modes du réseau et W1 | produits, forme et actions sur les vingt espaces de multiplicité |
| Lemme 7.15, équations | `R3/r3a_groebner.py` ; `R3-verif/setup_system.py`, `solve_groebner.py` | tenseurs et hauteurs sectorielles | idéal exact, quatre points géométriques, quotient réduit ; le second système impose davantage de coefficients caractéristiques |
| Lemme 7.15, candidats | `R3/r3a_omega_exact.py` ; `R3-verif/verify_solutions.py` | sortie de localisation `R3/idem_complex.json` pour le premier moteur ; équations exactes | quatre candidats sur Q(sqrt(-6)) dans le premier modèle, vérifiés exactement ; le second modèle utilise Q(sqrt(-2),sqrt(3)). Les indices des candidats et leurs coordonnées ne sont pas identiques entre les deux jauges |
| Proposition 7.16 | `R3/r3b_chi.py`, `R3/r3b_verify.py 0..3` ; `R3-verif/phi_chi.py` | caractère positif et quatre candidats | PCT préservant les courants et fixant les candidats ; transitivité des caractères de signes dans le second modèle ; stabilité de S par l'argument théorique des commutants |
| Remarque 7.18 | `R3-verif/compare_R3.py` | sorties des deux modèles | invariants de Gram des quatre solutions ; aucune comparaison naïve des coordonnées entre les jauges |
| Certificat portable, lemme 7.15 / prop. 7.16 | `verify_omega_certificate.py omega_certificate.json` | vecteurs, tenseurs et actions fournis dans le JSON | indépendance des 16 vecteurs, rang 8, 25 quadriques, égalité des idéaux sur Q, quatre points réduits et PCT ; ne recalcule pas le noyau ni les tenseurs |
| Reconstruction indépendante du certificat | `verification-independante/independent_tensors.py <Calculs-HU>` (inclut `independent_triples.py` et `independent_kernel.py`) | certificat, `data/t12_W1.json`, formules de modes de courants | racines et cocycle, 36 modes, matrice 2876 x 644, rang 628 à deux nombres premiers avec noyau rationnel explicite de dimension 16 ; 136 produits et formes, 16 traces et 136 traces de produits pour chacun des 20 types |
| Provenance du JSON | `export_omega_certificate.py R3 <sortie.json>` | trois pickles R3 et code de l'exportateur | réexport et comparaison du contenu JSON, y compris les empreintes de provenance ; importe le premier moteur |
| Contrôles exacts de §6 | `controles-categoriels/categories-exact-check.py` ; `section6/check_moments.py`, `verify_simplification.py` | formules codées, cyclotomie d'ordre 13 et 40 | S² = I, 144 coefficients de N_b, deux sommes de Gauss, identités symboliques des moments et simplification du cas m = 2 ; la preuve de reconnaissance reste théorique |
| Contrôles arithmétiques bornés de §6 | `section6/check_galois.py`, `check_m.py` | plages d'entiers codées | tests jusqu'à n < 400 et m <= 20000 ; certaines comparaisons sont flottantes. Ne prouvent pas la borne ou le lemme pour tous les paramètres |
| Diagnostics des moteurs | `R3/test_e8e8.py` ; `R3-verif/test_voa.py` | vecteurs de petits poids, échantillons fixés | tests Virasoro, commutateurs et symétries ; portée finie, sans vérification formelle générale |
| Diagnostics S2 et forme rationnelle | `R3/r3a_S2.py`, `r3a_S2_full.py`, `r3b_verify.py` ; `R3-verif/s2_space.py 0 all`, `s2_check.py 0`, `twist_check.py` | candidats et modes des deux modèles | vérifications supplémentaires en poids deux ; les tests aléatoires et ces diagnostics ne remplacent pas la preuve des commutants de la proposition 7.16 |
| Contrôles additionnels | `scripts-rapporteurs/F1/e8facts.py`, `sectors.py`, `K1/chk.py`, `verlinde.py`, `R3P/gb_sympy.py`, `check_gb.py` | tables, matrices et tenseurs R3 archivés | validations finies auxiliaires ; R3P se lance depuis une copie de R3 où ses deux scripts sont copiés |

## Utilisation dans III

III §9 réutilise les preuves de II, lemmes 7.14–7.15 et proposition 7.16, une fois les données (F2), les niveaux et la diagonale algébrique fournis. Son bloc C3 ajoute des pilotes et une confrontation avec les données de Zhu propres à III. Le présent paquet fournit le calcul de II qui y est importé, sans revendiquer le replay des 69 sorties C3 de III. La preuve principale de H_J de III utilise la reconstruction analytique et sa structure unitaire compatible, indépendamment du calcul PCT de cette voie secondaire.


## Achèvement des diagnostics S2 dans la livraison v02

Les trois programmes R3 `r3a_S2.py`, `r3a_S2_full.py` et `r3b_verify.py` sont rejoués pour chacun des quatre candidats. Le candidat 0 passe le parcours exhaustif de 1 200 paires ; les candidats 1–3 passent les 180 paires de leur parcours annoncé. Les compteurs de commutant, de fermeture et de stabilité PCT sont vérifiés explicitement. Les quatre JSON et les espaces régénérés sont confrontés aux références de v01 ; les différences de fins de lignes n'altèrent pas leur contenu mathématique.

La reconstruction secondaire `s2_space.py 0 all` est terminée, suivie de `s2_check.py 0`. La dimension, les multiplicités, la stabilité PCT dans tous les blocs et les cinq mineurs hermitiens positifs sont contrôlés. Le diagnostic de forme rationnelle `twist_check.py` est également rejoué. Les nouveaux reçus sont dans `checks/release-II-v14-v02/extended`, et les résultats comparés dans `ARTIFACTS.json`.

La fermeture supplémentaire garde le rôle de diagnostic indiqué dans le tableau. Le fait qu'elle ait été interrompue dans v01 n'affectait pas le parcours de preuve par commutants de la proposition 7.16 ; la livraison v02 en fournit désormais un nouveau journal complet.
