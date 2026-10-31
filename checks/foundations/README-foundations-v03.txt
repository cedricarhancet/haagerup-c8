Vérifications des sections2–8 pour la révisionv03, 31octobre2026.

Les scripts de classification d'origine sont conservés dans calculs-sources/classification.
Les relances suivantes ont été exécutées dans une copie de ce répertoire, pour préserver les sorties originales :
  python3 joint_uniqueness.py --d3
  python3 u2_embedding_test.py
Dépendances : Python3, sympy, mpmath. Les JSON et journaux de sortie sont joints.
Résultats : six recherches sur R_d1,R_d2,R_d3, chacune avec une solution et sans plafond atteint ; deux recherches de U2, chacune avec une solution. Les sorties JSON correspondent exactement aux références v02.
Il s'agit de relances de l'algorithme fourni, pas de deux implémentations indépendantes.

Vérificateur complémentaire :
  python3 check_analytic_identities.py --sources /chemin/vers/calculs-sources/classification
Ce vérificateur lit les sources et données de classification, vérifie les constantes explicites de la nouvelle preuve de rang, les constantes quartiques, les neuf identités résiduelles de poids4, et la polarisation qui exclut trois paires de niveaux par non-intégralité. Les nouvelles preuves du texte sont analytiques ; ces tests restent des contrôles séparés.

Aucune conclusion sur les certificats de génération de poids4 n'est tirée de ces relances.
