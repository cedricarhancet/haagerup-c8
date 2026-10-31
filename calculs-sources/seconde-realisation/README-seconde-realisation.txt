CALCULS — SECONDE RÉALISATION (CAS II, NIVEAUX (3,13)) — VERSION 05 DE L'ARTICLE — 30 octobre 2026

Environnement : Python 3.13, gmpy2, python-flint 0.9 (fmpq_mat, nmod_mat), sympy, numpy, scipy.
Les modules lattice_voa.py, cases.py, closure_explore.py, commutant_griess.py, griess_subalgebras.py, primaries.py,
isotropic_lines.py, exact_isotropy*.py, linalg.py, modp.py et les données *.json, *.pkl sont des COPIES de
Calculs/cas-II (v04), pour que ce dossier soit autonome. zh3_modular_data.py est une copie de Calculs/theorie.

STATUTS : « exact » = arithmétique exacte dans K = Q(√13, √−3) (quadruplets de rationnels, base 1, √13, √−3, √−39 ;
restriction de Weil K ≅ Q^4 et rangs rationnels fmpq_mat) ; « mod p » = arithmétique exacte dans F_p, plongement
cohérent de K (modp.py : r = s3 (s13 − 1)) ; un rang modulo p minore le rang sur K ; une appartenance modulo p n'est
PAS une preuve sur K.

ORDRE D'EXÉCUTION (run_all.sh) ET BESOINS
 1. exact_weight2_K.py {0,1}     exact  fermeture en poids 2 (681 produits, 0 échec), rang_K U_2 = 120.  < 1 Go, ~10 s
 2. invariants_K.py 0            exact  165 couplages normalisés T^2/(N N N), tous rationnels.          < 1 Go
 3. unitary_character.py         exact  caractères positifs chi avec phi_chi(e) = −f (famille à 2 paramètres).
 4. unitarity_test.py {0,1}      exact  phi_chi(U_1) = U_1, phi_chi(U_2) = U_2.                         < 1 Go, ~10 s
 5. phi_automorphism_check.py    exact  phi_1 automorphisme (5050 produits) ; positivité poids 1–3.      ~30 s
 6. choice_conjugation.py        exact  conjugaison complexe : U_2^(0) -> U_2^(1).                       ~15 s
 7. griess_cache.py, a3_griess_invariants.py, a3_griess_automorphisms.py
                                 exact  Aut(B, e) = {1, s} ; s échange les sous-algèbres 9 et 10.       ~1 min
 8. character_coefficients.py    exact  coefficients q^0..q^6 du vecteur de caractères (gamma = 0).     ~5 s
 9. closure_w4_modp.py 0 TESTS   mod p3 = 8388571 : S_3 (956), u_1 S_3 c S_3, R_b, M_b (1250, 10630),
                                 S_4 (6010) et décomposition, tests W4.                             1–3 Go, ~1 h
10. MODE=C2 closure_w4_modp.py 0 mod p3 : identités C2 (e(−1)e, x(−1)x, e(−1)x dans N'_4).             ~1–2 Go, ~10 min
11. closure_w5_modp.py 0         mod p3 : N_5 et test (S_2)(0) S_4 c N_5 (analogue de (4.2) de Riesen). 3–5 Go, ~3 h
 Ne pas lancer plus de deux calculs à la fois (2 cœurs, 8 Go). Ne jamais utiliser pkill avec un motif présent dans
 la ligne de commande. Les caches V._MODE / V._EMODE sont vidés après chaque clé.

SORTIES PRINCIPALES : *.out (journaux), *.json (résultats). Voir Verification-theorie-Haagerup-c2-v05.txt, partie C.
