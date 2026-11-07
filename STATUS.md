# État des vérifications — 7 novembre 2026

## Exécutions réalisées pendant cette livraison

| Contrôle | Résultat | Reçu |
| --- | --- | --- |
| Versions Hexagon et comparaison des sources | I v14 et III v10 identiques aux sources locales ; II v14 diffère seulement par deux commandes typographiques du résumé | `VERSION-CONCORDANCE.json` |
| Deux volumes F4 | SHA-256 conformes aux pièces jointes v07 ; tous les CRC ZIP passent, 986 fichiers de coordonnées | `EXTERNAL-ASSETS.json`, manifestes des volumes |
| Fondations | Identités analytiques, calculs rationnels modulaires et moments joints de U2 : PASS | `checks/release-I-v14/quick/RECEIPT.json` |
| Moteur de champs de réseau | 24 643 produits d'états et 310 produits de générateurs : PASS dans la fenêtre justifiée ; rejet de l'appel de poids neuf | Reçus `basis-modes.json`, `generator-modes.json` |
| Boule du réseau et courant simple | Énumération indépendante de norme au plus huit, mot explicite et annulations : PASS | Reçus `lattice-ball.json`, `simple-current-word.json` |
| Base de poids trois | Régénération : 956 vecteurs ; égalité exacte à la base archivée et contrôles T1, T2 et modes descendants : PASS | `checks/release-I-v14/certify_w3_validated.json` |
| Certificat F4 sur les coordonnées archivées | 6 010 pivots, 131 blocs, 18 368 mots de familles, 9 800 produits ; toutes les appartenances exactes et les 886 zéros : PASS | `checks/release-I-v14/f4/COMPLETE-F4.json` et journaux |
| Sorties séparées des données historiques | Commandes de démarrage : PASS dans un répertoire de sortie neuf | `checks/release-I-v14/smoke-isolation/RECEIPT.json` |

Les exécutions nouvelles sont datées et séparées dans `checks/release-I-v14/`. Les autres journaux et certificats proviennent des calculs antérieurs. L'ancien reçu F4 et son inventaire de sources sont conservés sans modification dans `provenance/`.

Le lanceur des zéros a été rendu portable sous Windows : l'adaptateur lit le compteur réel `PeakWorkingSetSize` en Kio, à la place du service Unix `resource`. Les limites de mémoire et de temps restent actives. Le vérificateur mathématique original est inchangé et son empreinte est contrôlée avant exécution ; l'adaptateur est couvert par le manifeste global et identifié dans `CHECK-METADATA.json`.

Les chemins du reçu sont sérialisés avec des séparateurs `/` pour l'agrégateur. Les étiquettes mathématiques des mots et les résultats sont conservés.

## Portée

Le replay F4 vérifie les coordonnées non nulles distribuées, les domaines exhaustifs, la provenance des pivots parmi les mots de familles et les identités exactes sur toutes les lignes. Les 886 zéros sont réévalués à partir des formules de modes ou de la borne d'énergie. **La régénération de tous les mots non nuls de poids quatre n'a pas été relancée pendant cette passe.** Le libellé historique « reconstruction » du reçu doit se lire avec cette précision. La base de poids trois, elle, a été régénérée et comparée exactement.

Les comparaisons indépendantes des produits de champs sont ciblées. Les coefficients des générateurs proviennent des données de construction. Ces contrôles ne sont ni une seconde génération exhaustive des coordonnées F4, ni une formalisation de toutes les preuves.

Les coordonnées et la relance F4 portent sur le premier choix associé à la sous-algèbre de Griess 9. Le passage aux trois autres choix est l'argument de conjugaison et de transport galoisien de la proposition 10.17 de I, et non quatre relances distinctes.

Les calculs donnent les entrées finies de la proposition 10.17. La génération forte et la `C_2`-cofinitude sont ensuite les conclusions de l'argument analytique et de l'argument sur l'algèbre de Poisson du théorème 10.21 de I ; elles ne sont pas des conséquences d'un simple rang modulo p et ne supposent pas `H_J`.

La rationalité forte, la positivité, les douze modules et le branchement sont traités dans l'annexe de III. Le théorème 6.1 de III prouve `H_J`, la complète unitarité et la comparaison VOA–réseau en utilisant les outils généraux indépendants de `H_J` de II. Cette livraison de I ne revérifie pas les calculs spécifiques à ces deux articles.

## Dépendances extérieures et archive publique

Les deux volumes F4 sont des pièces jointes stables de la release v07. Le transport originel des états de Riesen demande des fichiers de son paquet extérieur non redistribués ; ses résultats hérités sont identifiés comme tels. Certains anciens programmes supplémentaires demandent SciPy, Singular ou Sage. Aucune de ces dépendances supplémentaires n'est requise par les commandes de démarrage documentées de cette livraison.

Les SHA-256 contrôlent l'identité des fichiers ; ils ne remplacent pas l'expertise des arguments mathématiques. Les vérifications après téléchargement de la publication sont consignées séparément dans le rapport privé de livraison.
