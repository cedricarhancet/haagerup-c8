CALCULS — GÉNÉRATION FORTE ET C2-COFINITUDE EXACTES (v06) — 30 octobre 2026
Environnement : Python 3.13, gmpy2, python-flint 0.9 (fmpz_mat, fmpq_mat, nmod_mat), sympy, numpy, scipy.
Les modules hérités (lattice_voa.py, cases.py, closure_explore.py, griess_subalgebras.py, Kfield.py, modp.py, ...) et
les données (*.json, *.pkl) sont des copies de Calculs/seconde-realisation (v05).

PRINCIPE (exact sur K = Q(√13, √−3)) : vecteurs = combinaisons K-linéaires d'états monomiaux, calculées exactement ;
appartenance décidée bloc de poids par bloc de poids (Cartan de g' ⊕ g'') via la restriction de Weil K ≅ Q^4 : pivots
modulo p3 = 8388571 (mineur non nul ⇒ indépendance sur Q), système carré résolu exactement (flint), puis vérification
EXACTE de B c = v sur toutes les lignes en entiers (exact_engine.solve_block).

ORDRE D'EXÉCUTION (run_all.sh) ET BESOINS
 1. certify_w3.py      base exacte de S_3 (956), tests T1 (1871), T2 (2862), L(1), L(2), y(1), y(2) : 0 échec.
                       ≈ 4 min, < 1,5 Go. Écrit exact_S3_basis.pkl.
 2. count_generators.py dim A_2 = 27, dim N_3 ≥ 725 (mod p) : au plus 231 générateurs de poids 3. ≈ 1 min.
 3. certify_w4b.py     (PHASE_D non utilisé ici) pivots de N'_4 : 6010 mots (956 L(−1)S_3, 3867 y(−1)S_3, 1187
                       m(−1)U_2), écrits dans piv_w4/ (≈ 535 Mo) ; puis vecteurs tests u(0)s écrits dans tests_w4/.
                       ≈ 20 min pour les pivots, puis ≈ 4 à 37 min par vecteur u (x est le plus coûteux). < 2,5 Go.
 4. certify_w4_part.py vecteurs tests pour une partie des u, en parallèle (tests_w4_part/).
 5. solve_w4.py        résolution exacte : 9088 tests (u(0)S_3 pour les dix u, e(−1)U_2, x(−1)U_2), 0 échec. ≈ 8 min.
 6. phase_d_w4.py      phase (d) : tous les mots des cinq familles de N'_4 dans l'espace des 6010 pivots (dimension
                       exacte de N'_4). ≈ 45 min, < 3 Go.
 7. commutant_dims.py  (mod p) dim C_3 = 2, dim C_4 = 7 (blocs (0,0) de S_3, S_4). ≈ 1 min.
Ne lancer que deux calculs à la fois (mémoire effective ≈ 5 Go). Ne jamais utiliser pkill avec un motif présent dans
sa propre ligne de commande (incident du 30/10 : shell tué, sans effet sur les calculs).
Journaux d'incidents : certify_w4_run1_killed_memory.out, certify_w4_run2_pivots.out (voir Verification v06).
