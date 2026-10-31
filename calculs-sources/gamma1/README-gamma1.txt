CALCULS — VECTEUR γ = 1 (HAAGERUP c = 8), VERSION 01 — 30 octobre 2026

Environnement : Python 3.13, mpmath 1.3.0, numpy 2.5, sympy (conteneur Linux).
Tous les scripts s'exécutent depuis ce dossier ; les sorties de référence sont fournies.

Script                          Rôle                                                   Sortie
xi_exact.py                     Module : matrice fondamentale Ξ de Bantay--Gannon       (importé)
                                (récurrence exacte jusqu'à q^60, données EG (2.7)).
zh3_modular_data.py             Module : S, T de Z(H3).                                 (importé)
kac_order3.py                   Éléments d'ordre 3 de E8 (coordonnées de Kac) et        kac_order3_output.txt
                                dimensions des sous-algèbres fixes : 80, 86, 92, 134 ;
                                86 ⇔ E6 ⊕ A2 (prop. 4.1).
lattice_restriction_check.py    Caractères Θ_(i,j)/η^8 de V(E6⊕A2) = combinaisons des   lattice_restriction_check_output.txt
                                X_z^(γ) (lemme 4.3), pour γ = 0 et 1, jusqu'à q^4.
e6_rank_bound.py                Exact : aucun sous-espace de dimension 5 de la Cartan   e6_rank_bound_output.txt
                                de e6 ne rend symétriques les poids de 27 (prop. 4.4(3)).
weight2_relations.py 1|0        Espace M_2(ρ, μ) des formes de poids 2 (dim 5 / 4),     weight2_relations_gamma{1,0}_output.txt
                                relations exactes (8 / 9) ; contrôles : écart de         weight2_relations_gamma{1,0}.json
                                valeurs singulières, annulation exacte sur D(Ξ e_p),
                                modularité numérique (prop. 5.2). ~ 4 s chacun.
weight2_solve.py                Insertion des traces explicites (cor. 5.3).            weight2_solve_output.txt
gamma1_exclusion.py             Lemmes 6.2–6.3 et énumération : aucune s de dimension   gamma1_exclusion_output.txt
                                19 (théorème 6.1) ; liste des idéaux permis pour γ = 0.
gamma0_riesen_check.py          Test contre les espaces de plus bas poids ANNONCÉS par   gamma0_riesen_check_output.txt
                                Riesen (γ = 0) : accord exact pour les deux sl2.

Ordre conseillé :
  python3 kac_order3.py
  python3 lattice_restriction_check.py
  python3 e6_rank_bound.py
  python3 weight2_relations.py 1 ; python3 weight2_relations.py 0
  python3 weight2_solve.py
  python3 gamma1_exclusion.py
  python3 gamma0_riesen_check.py

Statut : weight2_relations.py contient la seule étape numérique (rang des conditions d'annulation
en τ = i et τ = e^{2πi/3}, 60 chiffres ; écart 12,0 contre 4,7e-58 pour γ = 1, 8,6 contre 3,4e-58 pour
γ = 0). Les relations sont reconnues comme rationnelles (erreur < 1e-40) puis contrôlées exactement
sur les formes explicites q d/dq (Ξ e_p). Tout le reste est en arithmétique exacte.
