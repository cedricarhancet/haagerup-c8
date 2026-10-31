CALCULS — CAS I : IDENTIFICATION AVEC RIESEN, EXCLUSION DE LA SOUS-ALGÈBRE 2 (v09, 31 octobre 2026)
Lancer depuis ce dossier ; les scripts lisent les modules de ../cas-II et les données de Riesen dans
../../Sources/cas-I/riesen/shared/ (definition.json, defining_states.json : données, jamais exécutées).
 python3 r1_identify_riesen.py        isométrie, w = s3 s4 s1 s3, cocycle, tore ; courants exactement identifiés (< 1 min)
 python3 r2_locate_riesen.py          e, f de Riesen dans la sous-algèbre de Griess 5 ; v sur la droite isotrope 1 (≈ 2 min)
 python3 r3_lines_and_Z.py            action des éléments de Z (j6) sur v : z_1 échange les deux droites (1,13) (≈ 3 min)
 python3 r4_sub2_explore.py           test (a) u_(2)w ∈ U_1 : vide pour (1,13) x (1,1) (documentation)
 python3 r5_sub2_testb.py 5 x x       CONTRÔLE : sous-algèbre 5, seuls les couples de droites (0,0) et (1,1) passent
 python3 r5_sub2_testb.py 2 0 2/5     sous-algèbre 2 : les 4 couples échouent (exact sympy)
 python3 r6_certify_sub2.py 2 0 2/5   même test, non-annulation CERTIFIÉE par Arb (400 bits) : sous-algèbre 2 exclue
Note : r5_sub2_testb.json est la sortie de la première version (types rationnels seulement) ; r5_testb_sub2.json et
r5_testb_sub5.json sont celles de la version générale ; r6_certify_sub2.json est le certificat.
 python3 r1b_omega_check.py           contrôle : ψ_R(ω_W) = ω, ψ_R(ω_A) = Sugawara de A, e + f = ω_C (exact)
Les droites isotropes de la sous-algèbre 2 (types (1,13), (1,1)) sont stockées dans r6_certify_sub2.json.
