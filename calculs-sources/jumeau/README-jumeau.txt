CALCULS — JUMEAU GALOISIEN ET CENTRALISATEUR DE U_1 (v09, 31 octobre 2026)
Lancer depuis ce dossier ; les scripts lisent les modules de ../generation-forte (et ../cas-II pour j6).
 python3 j1_centralizer.py II ; python3 j2_solve.py II      centralisateur, cas II : 8 solutions, (Z/2)^3   (< 2 min)
 python3 j3_elements.py                                      les 8 éléments sur V'_1 (monomiaux)              (< 1 min)
 python3 j4_action_on_B.py                                   action sur B (cas II) : identité ou réflexion s  (< 2 min)
 python3 j5_explicit_automorphism.py 110                     relèvement explicite, 1800 contrôles, z(U_2(9,c)) = U_2(10,c), exact (≈ 5 min)
 python3 j1_centralizer.py I ; python3 j2_solve.py I ; python3 j6_caseI_action.py     cas I : (Z/2)^3, 5 ↔ 6 (≈ 3 min)
Sorties : j1_*.pkl, j2_groebner_*.pkl, j3_z_*.pkl, j4_action_on_B.json, j5_explicit_110.json, j6_*.pkl/json.
 python3 j5_explicit_automorphism.py 011                     relèvement de l'automorphisme de diagramme de A2 (tore trivial) : U_2(9,c) → U_2(10,1−c)
Note : j3_elements.py et j6_caseI_action.py transcrivent la solution de Gröbner de j2 (vérifié conforme par le
relecteur v09 ; sorties j2_groebner_*.pkl).
