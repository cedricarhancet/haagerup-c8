Contrôles de la révision 03 — identification, centralisateurs et cas (I)

Dépendances : SymPy 1.14, gmpy2, python-flint, NumPy et SciPy ; mêmes versions que les contrôles v02.
Depuis la racine du dossier livré, avec calculs-sources/ contenant les sources :

python3 checks/identification/check_explicit_centralizers.py --sources calculs-sources --membership
python3 checks/identification/check_caseI_exclusions.py calculs-sources/generation-forte
python3 checks/identification/recheck_caseI_chart_certificate.py

Ces vérifications ont été exécutées pendant la révision, avec résultats PASS. Le premier script reconstruit les relèvements depuis leurs formules, sans les matrices historiques j3/j6 ; il vérifie les groupes, les courants, l'action sur B et la descente. Le deuxième produit les 13 systèmes projectifs rationnels excluant les candidats 0, 1, 3 et 4 du cas (I). Le troisième recalcule les idéaux unité en ordre lexicographique à partir du JSON, sans moteur VOA.

Les contrôles supplémentaires check-caseI-input-exactness.json et check-caseI-field.json attestent respectivement l'absence de coefficients flottants dans les données transportées et l'appartenance à K des coefficients des générateurs e, f et iv. Ces vérifications sont aussi intégrées au premier script.

inspect_centralizer_supports.py est une exploration auxiliaire qui a fourni les supports et leur indice deux dans le réseau. La preuve du manuscrit utilise directement les équations de caractères et les sous-systèmes de racines affichés.

Les caches historiques des primaires et les fonctions du moteur VOA restent des entrées. Leurs empreintes figurent dans la sortie du deuxième script. Il ne s'agit pas d'une formalisation complète du moteur.

Poids quatre du commutant : certificat régénéré (PASS, 365 secondes).

python3 checks/identification/check_commutant_c4.py --sources calculs-sources --pivots checks/identification/c4/input-pivot-block-00.pkl --output checks/identification/c4-rerun
python3 checks/identification/recheck_c4_matrices.py checks/identification/c4/c4-exact-witnesses.json

Le premier reconstruit les sept vecteurs, vérifie les 210 annulations affines et recalcule les deux matrices de Virasoro et les deux primaires. Le second vérifie indépendamment les identités entre matrices 7 × 7, les spectres, les multiplicités et le rang sept ; il ne refait pas les produits VOA.

L'entrée complète, ses étiquettes de mots et les résultats reconstruits sont dans c4/. Le fichier c4-exact-witnesses.json fournit les coefficients rationnels complets dans la base 1, sqrt(13), sqrt(-3), sqrt(-39), sans pickle ni objets dépendant d'une version de Python. L'ordre L1/L2 du code est l'inverse de celui du papier : les charges centrales sont croissantes dans le code.

Le script de régénération du bloc central est fourni parmi les contrôles de certificats. On peut aussi utiliser le bloc d'entrée figé (43 Mo) et ses étiquettes de mots. Il contient 280 mots de N4 indépendants modulo le premier choisi. Le certificat global (F4) est maintenant vérifié : 131 blocs, 6010 pivots, 18368 mots des familles et 9800 produits, sans erreur. Il établit N4′ = S4 par le théorème de génération forte. La liste des 280 pivots utilisée ici a été comparée exactement à celle du bloc correspondant du certificat global ; les listes sont identiques (reçu central-pivots-agree.json, PASS). Le présent contrôle C4 et le certificat global se complètent ; ils permettent désormais d'affirmer les dimensions et les bi-poids du commutant de S sans hypothèse (F4) laissée ouverte.
