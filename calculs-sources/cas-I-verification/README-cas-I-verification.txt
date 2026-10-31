CALCULS — VÉRIFICATION INDÉPENDANTE DE U^R = ⟨A, e, f, v⟩ (axe B, v10, 31 octobre 2026)
Lancer depuis ce dossier ; modules de ../cas-II ; générateurs de Riesen transportés : ../cas-I/r1_riesen_transported.pkl.
 PSTART=8388608 python3 b1_UR_closure_modp.py   p = 8 388 409 : W_2 = 120, S_3 = 956, S_4 = 6010 (≈ 12 min)
 PSTART=8000000 python3 b1_UR_closure_modp.py   p = 7 999 993 : idem
 Options : SAVEALL=dossier (tous les blocs, gros : scratch), W3P=1 (test u_(1)S_3 ⊂ S_3, u = e, f, v), NOW4=1 (arrêt
 après le poids 3), VCHOICE=conj | rand (contrôles : droite conjuguée / vecteur non isotrope).
Sorties : b1_UR_closure_p8388409.json (régénéré par la relance complète, journal b1_UR_closure_p1_rerun.out),
 b1_UR_closure_p7999993.json, *_conj.json, *_rand.json, *_W3prime.json et les journaux .out.
Incident : la première version de b1_UR_closure_p8388409.json a été écrasée par l'exécution W3P (NOW4) ; relance
complète faite (mêmes résultats que b1_UR_closure_p1.out). Blocs (0,0) de S_4 (b1_S4_00_p*.npz) : dans le scratch.
