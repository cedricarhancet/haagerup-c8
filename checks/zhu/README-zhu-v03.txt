VERIFICATIONS DE LA SECTION ZHU POUR LA V03

1. check_simple_current_table.py utilise uniquement la bibliothèque standard.
   Il vérifie les coefficients, le cocycle, les 42 termes croisés nuls,
   les normes et les bornes de poids de la nouvelle démonstration manuscrite.
   Commande : python3 check_simple_current_table.py --output result.json
   L'option --source-dir <calculs-sources/generation-forte> compare également
   les deux vecteurs du tableau avec exact_engine.setup(0). Cette comparaison
   a été exécutée et réussit (simple-current-table.json). Elle requiert les
   dépendances de l'ancien moteur exact (gmpy2, python-flint, SymPy).

2. z4_upper_algebra_modp.py est une copie inchangée du programme historique.
   Il a été relancé à partir des 19981 relations projetées conservées dans
   corner_rows_modp.pkl. Définir HAAG_PDIR comme chemin absolu de ce dossier,
   puis lancer le script avec Python, NumPy et python-flint installés.
   Le résultat est identique au JSON archivé : dimension762, forme de trace
   de rang762, centre de dimension12 et mêmes blocs/rangs des sommets.
   La génération VOA des relations et leur projection n'ont pas été relancées.

3. check_rewriting_system.py est un contrôle nouveau de la présentation
   finale archivée. Il utilise une réduction gauche/courte, différente de
   celle du programme de découverte. Il vérifie l'ordre décroissant,
   les extrémités des chemins, les 131698 recouvrements et les151 inclusions.
   Le script historique ne vérifiait pas explicitement ces inclusions.
   Elles sont indispensables puisqu'il y a5221 règles quadratiques et438
   règles cubiques. Tous les résidus sont nuls. Les chemins irréductibles
   sont239,506,5,0 aux longueurs1,2,3,4; avec12 sommets, dimension762.
   Commande : python3 check_rewriting_system.py <chemin-de-ce-dossier>
   Ce contrôle ne prouve pas à lui seul que les règles proviennent des
   identités VOA, ni que les noyaux calculés modulo p se relèvent à K.

4. Les logs et JSON du second premier4194277 ont été retrouvés dans Dropbox
   et sont fournis sous calculs-sources/zhu/p4194277/. Les invariants
   concordent, à une permutation des blocs près. Cette seconde chaîne n'a
   pas été relancée pour la présente révision.

5. Aucun de ces résultats ne prouve l'hypothèse(Z) en caractéristique zéro.
   La construction des relations exactes sur K, le modèle local et le
   contrôle de sa présentation restent nécessaires.
