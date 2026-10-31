CALCULS — CAS (II) DES VOA DE HAAGERUP À c = 8 (VERSION 04 DE L'ARTICLE) — 30 octobre 2026

Environnement : Python 3.13, gmpy2 2.3.2, python-flint 0.9.0 (FLINT/Arb), sympy, numpy, scipy.
Tous les calculs qui servent de preuve sont exacts (rationnels, Q(sqrt13), Groebner sur Q ou Q(sqrt13)) ou
certifiés par arithmétique d'intervalles (Arb) ; ces derniers ne servent qu'à certifier qu'une quantité est NON nulle.
Les calculs modulo p (closure_modp.py) sont des vérifications exactes en caractéristique p : une égalité modulo p
n'est PAS une preuve sur le corps de nombres ; un rang modulo p minore le rang exact.

ORDRE D'EXÉCUTION (voir run_all.sh)
 1. test_lattice_voa.py          tests de l'implémentation de V(E6+A2) (dimensions 1, 86, 1370, 11586 ; Virasoro ;
                                 translation ; formule des commutateurs ; Sugawara = omega ; <omega,omega> = 4).
 2. cases.py                     réalisation de U_1 = g' + g'' dans V'_1 (cas I et II), H', H'' dans la Cartan du réseau.
 3. check_invariants.py          invariants de V'_2 (7 / 9), dim B (4 / 5), vecteurs A-primaires de poids 1 (= R_b) :
                                 contrôle indépendant des lignes R_b du tableau 1.
 4. commutant_griess.py          B = Com_{V'}(A)_2, constantes de structure, c_C (170/41, 18/5).
 5. griess_subalgebras.py        sous-algèbres de dimension 2 de B contenant e : 7 (I), 15 (II) ; charges centrales.
 6. primaries.py                 espaces Pi_tau, matrices de b_1, b_2 (b dans B) ; sorties primaries_{I,II}.pkl.
 7. isotropy_test.py             test certifié (Arb, Rump) type par type, toutes les sous-algèbres.
 8. exact_isotropy.py CAS k      analyse exacte (valeurs propres rationnelles, Groebner) ; run_exact.sh pour la boucle.
 9. exact_isotropy_s13.py CAS k  même analyse sur Q(sqrt13) (doublement rationnel) + Arb pour les facteurs non scindés.
10. isotropic_lines.py CAS k     espaces admissibles et droites isotropes (cas II k = 9 ; cas I k = 5).
11. closure_explore.py II 9      tests de fermeture par paires (flottant) : 2 choix sur 256.
12. closure_certify.py II 9      certification des rejets (forme linéaire exacte + boule Arb) : 62/72 certifiés,
                                 suffisants pour éliminer 254 choix sur 256.
13. closure_full.py II 9         fermeture complète en poids 2 des 2 choix restants (flottant).
14. closure_modp.py II 9 0 [w3]  vérification exacte modulo p : fermeture en poids 2 ; S_3 engendré et clos sous
                                 g'+g'' : rang 956 ; traces ; u_2 S_3 c U_2 et u_3 S_3 c U_1 (variables PSTART,
                                 NOC3, C3NAMES, C3U1). Choix -1 = choix rejeté (témoin).
15. weight3_explore.py II 9      exploration flottante du poids 3 (borne 932 avant clôture ; structure de dimension 956).
16. closure_explore2.py I 5      contrôle du cas (I) (flottant) : 2 choix de U_2 survivent.
17. charpolys.py, eigen_reduction.py : outils de diagnostic (polynômes caractéristiques, sous-espaces stables).

MODULES : lattice_voa.py (VOA de réseau exacte), linalg.py (algèbre linéaire flint), modp.py (réduction modulo p de
Q(sqrt13, sqrt-3) ; NB : r = sqrt(6) sqrt(-7+sqrt13) = sqrt(-3)(sqrt13 - 1)).

RÉSULTATS PRINCIPAUX
 - Cas II : sous-algèbres 0–8, 11–14 exclues (exact ou certifié) ; 9, 10 (charges 182/95, 32/19) survivent ;
   pour 9 : un seul espace admissible (dim 2) par type, 2 droites isotropes chacun, 2 choix de U_2 sur 256 après
   fermeture par paires (rejets certifiés) ; par Galois, idem pour 10.
 - Vérifié modulo p1 = 134217613 et p2 = 133217647 : dim U_2 = 120, fermeture en poids 2, dim S_3 (engendré, clos)
   = 956, neuf traces = celles de u3_moments (v03), u_2 S_3 c U_2 et u_3 S_3 c U_1 pour u = e, x et les vecteurs de
   plus haut poids des huit types primaires, y_1 S_3 c U_2 et y_2 S_3 c U_1 pour y dans U_1 (sorties
   closure_modp_II_9_w3_C3all.out [e, x, (0,4) ; processus ensuite tué faute de mémoire], _C3.out [(1,1), (3,3)],
   _C3rest.out [(1,3), (1,5), (1,7), (2,2), (2,6), U_1]).
 - Témoin négatif : avec un choix incohérent de racines carrées (p2, première exécution), la configuration n'est pas
   une solution et engendre un espace de dimension 1663 (fichier ..._branche_incoherente_temoin_negatif.out).
 - Cas I : sous-algèbres 0, 1, 3, 4 exclues ; 2, 5, 6 restent ; 5, 6 donnent (2400/1271, 70/31) (annoncé) ;
   exploration flottante : 2 choix de U_2 pour 5 (non certifié).
