# Diagnostics supplémentaires achevés

`SUMMARY.json` donne le bilan des quinze nouvelles commandes. `extended/RECEIPT.json` identifie l'environnement, les sources utilisées, le nombre de tâches et les cinq copies isolées. Chaque sous-dossier CI0, CI1, CI2, CI3 et second-engine contient son reçu, les commandes exécutées, leurs durées, codes de sortie et empreintes, ainsi que les sorties standard et d'erreur intégrales.

`ARTIFACTS.json` lie les données régénérées à leurs exécutions et aux références antérieures. Les fichiers R3 correspondants sont dans `../../Calculs-HU-v08/R3/`. Le grand espace secondaire est conservé dans `generated/R3-verif/s2_space_0_all.pkl.gz`, avec ses empreintes compressée et décompressée. `MATH-CODE-IDENTITY.json` établit que les sources mathématiques de calcul sont inchangées depuis v01.

Depuis `II/`, `python verify_completed_diagnostics.py` vérifie les reçus, compteurs et empreintes avec la bibliothèque standard. Pour refaire réellement les calculs, utiliser `python reproduce.py --group extended --output <nouveau-dossier-exterieur>`, avec les dépendances précisées dans REPRODUCIBILITY. Le contrôle des données livrées et la régénération complète sont deux opérations distinctes.
