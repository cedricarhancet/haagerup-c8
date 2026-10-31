Vérifications supplémentaires, révision v03 — 31 octobre 2026

1. check_griess_charts.py
   Recalcule les données rationnelles de Griess à partir du moteur historique,
   puis les bases de Gröbner et les dimensions des algèbres des cartes affines
   disjointes. Il vérifie séparément les points proposés et leur distinction.
   Résultat : I, dimensions (7,0,0), 7 points ; II, (11,2,1,1), 15 points.
   Les points atteignent la dimension totale du schéma fini : les listes sont
   exhaustives et les points sont simples. L'argument de comptage externe
   n'est pas nécessaire à ce contrôle.
   Vérifie aussi que les sept axes fixes engendrent H de dimension 3 et que
   det(Gram|H)=147968/2025, non nul.
   Sorties : griess-charts.json (PASS), griess-charts.log.

2. check_primary_spaces.py
   Relance exact_isotropy_s13.py pour II, candidat 9, et refuse toute exclusion
   indécidée ou tout facteur non linéaire non traité. Vérifie un seul espace
   admissible de dimension 2 pour chacun des huit types. Appelle ensuite
   isotropic_lines.py, uniquement après ce contrôle d'exhaustivité. Recalcule
   les deux poids de Virasoro exactement à partir de xi, lambda et Delta,
   puis les compare aux seize valeurs rationnelles du tableau. Vérifie les
   deux lignes distinctes, toutes les équations quadratiques et leurs normes
   bilinéaires non nulles. Ce contrôle ne repose pas sur la reconnaissance
   numérique du poids effectuée dans le script historique.
   Sorties : primary-spaces.json (PASS), primary-spaces.log.

Commandes depuis la racine du dossier complet, après installation des
versions python-flint, gmpy2, sympy indiquées dans les exigences du dossier :

python3 checks/generators-isotropy/check_griess_charts.py --source-dir calculs-sources/cas-II --output checks/generators-isotropy/griess-charts-rerun.json
python3 checks/generators-isotropy/check_primary_spaces.py --source-dir calculs-sources/cas-II --output checks/generators-isotropy/primary-spaces-rerun.json

Limite : ces scripts régénèrent les données avec les formules VOA originales.
Ils ne constituent pas une vérification indépendante du moteur de produits
VOA, ni de l'ensemble du manuscrit. Les sorties sont des résultats de relance,
pas des certificats formalisés dans un assistant de preuve.
