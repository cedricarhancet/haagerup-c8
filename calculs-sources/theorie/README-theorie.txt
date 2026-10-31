CALCULS — HAAGERUP c = 2, VOIE THÉORIQUE, VERSION 01
29 octobre 2026

Environnement utilisé: Python 3.13, mpmath 1.3.0, numpy 2.5 (conteneur Linux). Aucune autre
dépendance. Chaque script s'exécute depuis ce dossier; les sorties de référence sont fournies.

Script                         Rôle                                              Durée   Sortie
zh3_modular_data.py            Données S, T de Z(H3) (étiquettes Evans--Gannon), < 1 s   zh3_modular_data_output.txt
                               S^2 = I, (S T_c)^3 = p S^2, unitarité,
                               dimensions quantiques, 60 chiffres.
eg_c8_verification.py          Matrice fondamentale Xi de Bantay--Gannon à       ~ 6 s   eg_c8_verification_output.txt
                               partir de Lambda et chi de EG (2.7); récurrence
                               exacte jusqu'à q^60 (entiers); test
                               Xi(-1/tau) = S Xi(tau) (erreur <= 1.4e-50);
                               vecteurs admissibles X^(gamma), gamma in {0,1}.
bg_trace_criterion.py          Critère exact de la trace [BaG07, (2.22) p. 6]:  < 1 s   bg_trace_criterion_output.txt
                               Tr Lambda = 5 = membre de droite.
e8_rank2_table.py 36           Table 1: N(K) maximal pour les 32 réseaux pairs  ~ 4 s   e8_rank2_table_output.txt,
                               binaires de déterminant <= 36; écrit aussi les           e8_rank2_table_DMAX36.txt,
                               certificats (u, v).                                       e8_rank2_certificates_DMAX36.txt
check_e8_certificates.py       Vérification indépendante et exacte des          < 1 s   check_e8_certificates_output.txt
  e8_rank2_certificates_DMAX36.txt  certificats (autre base de E8, racines
                               listées directement), et complétude de la liste
                               des formes réduites.

Commandes:
  python3 zh3_modular_data.py
  python3 eg_c8_verification.py
  python3 bg_trace_criterion.py
  python3 e8_rank2_table.py 36
  python3 check_e8_certificates.py e8_rank2_certificates_DMAX36.txt

Fichier annexe: e8_rank2_table_DMAX16.txt (première exécution avec det <= 16, cohérente avec la
table finale).

Statut: les scripts vérifient (1) la classification d'Evans--Gannon utilisée au théorème 3.1
(numérique 50 chiffres + récurrence exacte + critère exact de trace), (2) la borne N(K) >= 16
de la table 1 (certificats exacts). Ils ne vérifient ni les énoncés catégoriques de la
section 6, ni les données annoncées de Riesen et Huang utilisées à la proposition 7.2.

Script retiré: det_check.py (test de déterminant invalide, poids 84 et non 0). Il est conservé
pour traçabilité dans Notes/theorie/det_check-INVALIDE.py et n'est pas un certificat.
