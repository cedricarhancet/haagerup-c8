#!/bin/bash
# Ordre d'exécution des calculs « seconde réalisation » (cas II, sous-algèbre 9). Voir README-seconde-realisation.txt.
# Mémoire : < 1 Go pour 1–7 ; 1–3 Go pour 8–10 ; ne pas lancer plus de deux calculs à la fois (2 cœurs, 8 Go).
set -e
python3 exact_weight2_K.py 0 > exact_weight2_K_0.out          # 1. fermeture exacte en poids 2 sur K (choix 0)
python3 exact_weight2_K.py 1 > exact_weight2_K_1.out          #    idem (choix 1)
python3 invariants_K.py 0 > invariants_K_0.out                # 2. couplages normalisés (rationnels)
python3 unitary_character.py                                  # 3. caractères positifs chi (famille à 2 paramètres)
python3 unitarity_test.py 0; python3 unitarity_test.py 1      # 4. phi_chi(U_1) = U_1, phi_chi(U_2) = U_2 (exact)
python3 phi_automorphism_check.py                             # 5. phi_1 automorphisme, positivité en poids <= 3
python3 choice_conjugation.py                                 # 6. A3 : conjugaison complexe échange les deux choix
python3 griess_cache.py; python3 a3_griess_invariants.py; python3 a3_griess_automorphisms.py   # 7. A3 : Aut(B, e)
python3 character_coefficients.py > character_coefficients.out                            # 8. caractères (q^0..q^6)
PSTART=8388608 python3 closure_w4_modp.py 0 '(1, 7);(2, 6);(3, 3);(1, 5);x' > closure_w4_modp_0.out   # 9. poids 3–4, M_b
MODE=C2 PSTART=8388608 python3 closure_w4_modp.py 0 > closure_c2_modp_0.out                         # 10. test C2
PSTART=8388608 python3 closure_w5_modp.py 0 > closure_w5_modp_0.out                                 # 11. poids 5, (4.2)
