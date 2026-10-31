CALCULS — CLASSIFICATION DE U_1 (HAAGERUP, c = 8), VERSION 03 — 30 octobre 2026

Environnement : Python 3.13, sympy, mpmath 1.3.0, numpy, scipy 1.18 (HiGHS, uniquement pour
u3_embedding_test.py, solution vérifiée en arithmétique exacte). Tout le reste est en arithmétique
rationnelle exacte (fractions.Fraction, sympy.Rational).
Lancer : ./run_all.sh depuis ce dossier (environ 10 minutes) ; journal dans run_all_log.txt.

MODULES (importés)
  xi_exact.py, zh3_modular_data.py   matrice fondamentale Xi (Bantay--Gannon, données d'Evans--Gannon) et S, T.
  vvmf_exact.py                      développements exacts (E2, E4, E6, Delta, J, DJ, eta^-8, 1/Delta...),
                                     formes vectorielles, permutation pi et contrôle exact de S(pi z, pi w) = S(z,w).
  obstruction.py                     principe d'obstruction : relations exactes entre coefficients d'un poids k
                                     donné (accouplement avec des formes pour la représentation duale) et
                                     éléments exacts de M_k(rho, mu).
  e6_tools.py                        poids de 27, racines de E6, constantes quartiques (3 et 18).
  e8_joint.py                        e8 explicite (Frenkel--Kac), e6 + sl3, triplets sl2, centralisateurs,
                                     traces jointes et poids joints exacts.

AXE A — PREUVE EXACTE DU THÉORÈME 7.1 (gamma = 1)
  weight2_obstruction_exact.py g M   relations exactes de poids 2 (g = 1 : 8 ; g = 0 : 9) ; COÏNCIDENT avec
                                     celles de v02 ; borne inférieure par formes exactes (rang 5 / 4) ;
                                     injectivité (Partie 3) ; d'où dim M_2 = 5 / 4. Certificats JSON.
  dump_relations.py                  listes de relations de poids 2, 4, 6 (gamma = 0, 1) + rangs exacts.
  gamma1_exclusion_exact.py          théorème gamma = 1 en arithmétique exacte (idéaux simples, centre u >= 1/9,
                                     énumération : aucune s de dimension 19).

AXE B — POIDS 4 ET 6
  weight4_formula_check.py           test numérique de la formule de Lemme 6.1 (n = 4) sur V(A2), V(sqrt6 Z), V(E8).
  weight6_formula_check.py           idem n = 4, 6 ; échec si un coefficient est perturbé.
  weight4_solve.py                   relations de poids 4 avec données explicites : Q1..Q5 en fonction de
                                     (u, T, T4, tau6, Q6) ; Q5 = (12T + 2T4 - 27u^2 - 90u + 105)/26.
  riesen_quartic_check.py            TEST : traces quartiques des espaces annoncés par Riesen (deux facteurs) : accord.

AXE C — U_1 POUR gamma = 0
  centre_gamma0.py                   s est semi-simple (Q5 = -6/13 < 0 pour h central avec p2(h) = 0).
  sl2_e6_superset.py                 sur-ensemble calculé des sl2 de e6 (27 étiquetages, indices distincts).
  e6_centralizers.py                 e6 explicite ; centralisateur de l'orbite d'indice 35 : sl2 d'indice 1 ;
                                     d'indice 9 (2A2+A1) : sl2 d'indice 3.
  e6_orbit_labels.py                 contrôle indépendant des étiquettes de Bala--Carter (dimension des orbites,
                                     indices) ; corrige l'étiquette A2+2A1 (fausse) en 2A2+A1.
  sl2_ideal_feasibility.py [--level] idéaux simples possibles : (1,0), (3,0), (9,4), (35,4).
  sl2_pair_feasibility.py            couples (méthode directe, conditions nécessaires sur d2, d3).
  sl2_pair_exact.py                  couples, décision exacte (modules de multiplicité N_a) : survivent
                                     (1,0)+(35,4) [cas I] et (3,0)+(9,4) [cas II].
  weight6_single.py [--pass2]        relations de poids 6 (6 relations) pour chaque idéal.
  weight6_mixed.py                   relations mixtes x^4y^2, x^2y^4.
  sextic_joint.py                    faisabilité jointe complète (poids <= 6) ; TEST : accord avec Riesen (d1-d3).
  weightn_single_check.py 8 10       relations de poids 8 et 10 ; identité de poids 8 vérifiée pour les 4 idéaux.
  degree8_feasibility.py             traces de degré 8 réalisables pour d2, d3.
  joint_uniqueness.py --d3           unicité des structures de R_d1, R_d2, R_d3 (une seule dans chaque cas).
  known_sectors_table.py             R_b, R_a = R_c0, R_c1, R_c2 dans les deux cas.
  u2_moments.py                      traces sur U_2 ; TEST : accord avec la décomposition de U_2 de Riesen.
  u2_embedding_test.py               U_2 : décomposition unique dans chaque cas (bornes A_2 + omega <= U_2 <= V'_2).

AXE D — VERS L'UNICITÉ
  u2_parameter_count.py              dimension de l'espace des U_2 possibles : 70 (cas I), 87 (cas II).
  cartan_conjugate.py                représentants simultanés de H' et H'' dans la Cartan standard.
  lattice_branching.py               caractères joints de V'_n (n <= 3) via le réseau E6 + A2.
  u3_moments.py                      traces sur U_3 (dim 956).
  u3_embedding_test.py               existence d'un s-module U_3 compatible (L(-1)U_2 <= U_3 <= V'_3) : oui
                                     dans les deux cas (solution HiGHS vérifiée exactement).

STATUT
  - Aucun rang numérique n'intervient plus. Les preuves reposent sur des calculs finis exacts.
  - Seule entrée non vérifiée en arithmétique exacte : la modularité des colonnes de Xi (théorème publié
    d'Evans--Gannon, confirmé numériquement à 50 chiffres dans Calculs/theorie).
  - Les tests numériques (formules de poids 4 et 6 sur des réseaux) et les comparaisons avec les données
    annoncées (Riesen) ne servent jamais d'hypothèse.
  - u3_embedding_test.py utilise un solveur flottant pour TROUVER une solution ; la solution est
    vérifiée exactement. Aucune conclusion d'infaisabilité n'est tirée de ce solveur.
