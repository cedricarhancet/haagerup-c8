CALCULS — COMMUTANT C = Com_S(A) DU CAS II (axe A, v10, 31 octobre 2026)
Lancer depuis ce dossier ; les scripts lisent les modules et données de ../generation-forte (et le point de contrôle
../seconde-realisation/ckpt_S4_0_p8388571.npz pour c1 ; la base exacte piv_w4 de S_4 pour c2).
 python3 c1_biweights.py    modulo p3 : dim C_4 = 7, spectres de L(0) pour les deux idempotents ; contrôle négatif
                            (x conjugué). Sorties : c1_biweights.json (+ c1_C4_basis_p8388571.npz, dans le scratch).
 python3 c2_exact_C4.py     EXACT sur K = Q(√13, √−3) : base exacte de C_4 (7 vecteurs), matrices exactes de L(0),
                            identités exactes, vecteurs primaires de bipoids non entiers (≈ 6 min).
                            Sorties : c2_exact_C4.json, c2_exact_C4.out (+ pickles dans le scratch v10-intermediaires).
Attention : dans les scripts, L1 / « omega_1 » désigne c = 32/19 ; dans l'article l'indice 1 désigne c = 182/95.
