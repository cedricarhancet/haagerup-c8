# État de validation — livraison III

Passe de publication du 7 novembre 2026. Les compteurs suivants viennent des nouveaux reçus inclus, non d'un simple examen des journaux antérieurs.

- A1/A4 : 10 commandes terminées avec succès ; 320 contrôles détaillés A4, aucun `FAUX`.
- C3 : certificat portable, contrôles catégoriques, tables I/III/Zhu et comparaison de 131 fichiers communs avec II publié réussis ; 52 fichiers propres à III distingués. Comparaison séparée avec la référence complète HU-v08 : 47 données et 59 journaux concordants.
- Zhu : données exportées contrôlées ; 131 698 chevauchements et 151 ambiguïtés d'inclusion tous résolus ; les deux inclusions d'idéaux vérifiées. Le calcul des 554 règles de complétion a abouti intégralement.
- Fibre spéciale : dimension 762, rang de la forme trace 762, centre de dimension 12 ; blocs, poids exacts et branchements concordants.
- Projection : 19 981 lignes régénérées depuis les 1 060 polynômes fournis, avec égalité exacte et contrôle modulo p.
- Depuis les modes du réseau : base abaissée et 10 premières lignes (1 155 paires) régénérées ; dix empreintes identiques à l'échantillon historique.
- Moteurs/fondations Zhu : 19 fichiers d'amont de I identiques octet par octet à la release I fixée.
- Source III v11 compilé avec références stabilisées ; six pages pertinentes ont été inspectées visuellement.

Les nouveaux reçus sont dans `checks/release-III-v10-v01/`. Ils donnent les commandes, codes de sortie, durées et empreintes des journaux réellement produits pendant cette passe. Les préfixes des dossiers locaux ont été remplacés par `[WORK]` et `[MANUSCRIPTS]`, et les fins de ligne ont été normalisées en LF. Les champs `raw_stdout_sha256`/`raw_stderr_sha256` conservent l'empreinte des journaux bruts ; `stdout_sha256`/`stderr_sha256` couvrent les journaux publics livrés.

Les journaux présents dès les archives antérieures restent historiques. En particulier, la production C3 complète, les 120 lignes de produits de modes de Zhu, le calcul des noyaux `z2x.py` et la production initiale des règles `z4x.py` ne sont pas annoncés comme régénérés pendant cette passe. Les règles et relations fournies ont fait l'objet d'un nouveau contrôle exhaustif de confluence et des deux inclusions d'idéaux.

Les bibliothèques de calcul exact et le moteur des champs de réseau restent des dépendances logicielles. Ces contrôles ne sont pas une formalisation dans un assistant de preuve. Les objets pickle fournis ne doivent être chargés qu'après vérification des archives et de leur provenance.
