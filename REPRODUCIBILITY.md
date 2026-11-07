# Reproduire la livraison de I v14

## Copie propre et environnement

Télécharger l'archive de la [release figée](https://github.com/cedricarhancet/haagerup-c8/releases/tag/article-I-v14-supplements-v01) et vérifier son SHA-256 avec le fichier `SHA256SUMS-article-I-v14-supplements-v01.txt`. Décompresser dans un chemin assez court sous Windows. Sous Python 3.12 :

Sous PowerShell, l'empreinte de l'archive se calcule avec `Get-FileHash -Algorithm SHA256 .\Haagerup-c8-article-I-v14-supplements-v01.zip` ; sous Unix, avec `sha256sum Haagerup-c8-article-I-v14-supplements-v01.zip`. Comparer la valeur au fichier de sommes joint à la release. Le fichier `.gitattributes` conserve aussi les octets lors d'un clone Git, sans conversion des fins de ligne.

```sh
python verify_manifest.py
python -m venv .venv
```

Activer `.venv/bin/activate` sous Unix, ou `.venv\Scripts\Activate.ps1` sous PowerShell, puis :

```sh
python -m pip install -r checks/certificates/deliverable/requirements-exact.txt
python reproduce.py --group smoke --output ../smoke-I14
python reproduce.py --group quick --output ../quick-I14
```

L'environnement testé est Python 3.12.14, SymPy 1.14.0, gmpy2 2.3.2, python-flint 0.9.0, NumPy 2.3.5 et mpmath 1.3.0. Chaque répertoire de sortie doit être nouveau. Ne pas utiliser `-O`, `-OO` ou `PYTHONOPTIMIZE` : des assertions font partie des vérificateurs. Les commandes sont identiques sous Windows et Unix ; adapter seulement l'activation de l'environnement et les continuations de ligne du shell.

Le manifeste doit être vérifié avant toute exécution qui pourrait créer des caches ou modifier des données. Les sorties des nouveaux contrôles sont placées hors de la copie distribuée. `reproduce.py --group smoke` est le démarrage sans paquet extérieur ; les groupes `foundations`, `modes` et `lattice` peuvent être relancés séparément.

## Certificat complet de poids quatre

Télécharger les deux volumes indiqués dans [EXTERNAL-ASSETS.json](EXTERNAL-ASSETS.json), en conservant leurs noms, dans `../release-assets`. Le téléchargement avec contrôle intégré se lance par `python fetch_assets.py --output ../release-assets` ; il demande un environnement Python avec HTTPS opérationnel. Vérifier leurs tailles et SHA-256. Exécuter cette commande sur une seule ligne si le shell ne prend pas en charge les continuations :

```sh
python checks/certificates/deliverable/run_f4_verification.py --archive ../release-assets/Certificats-poids-quatre-v03-partie-1.zip ../release-assets/Certificats-poids-quatre-v03-partie-2.zip --output ../f4-replay-I14
```

Le lanceur contrôle les empreintes, extrait exactement 986 fichiers de coordonnées, vérifie les domaines complets et les appartenances exactes sur toutes les lignes, réévalue chaque zéro et n'accepte que l'agrégation exhaustive. Le reçu final attendu est `../f4-replay-I14/COMPLETE-F4.json`, avec `status = PASS`, `pivot_words = 6010`, `cartan_blocks_checked = 131`, `family_inputs = 18368`, `product_inputs = 9800`, `independently_verified_zero_inputs = 886`, `errors = 0`.

Sous Windows, `run_zero_products_portable.py` fournit l'accès au compteur réel de mémoire maximale du processus (`PeakWorkingSetSize`, converti en Kio). Le vérificateur mathématique original des zéros est exécuté sans modification après contrôle de son empreinte ; ses limites de mémoire et de temps restent actives. L'adaptateur et le lanceur sont couverts par le manifeste global. Les 131 blocs ont été validés avant cette adaptation système, puis les zéros et l'agrégation ont été relancés séparément. Le nouveau lanceur complet n'a pas été réexécuté une seconde fois depuis zéro pendant cette passe.

Le pilote écrit les chemins des manifestes avec des séparateurs `/`, format attendu par l'agrégateur sur toutes les plateformes. Cette conversion de représentation ne change ni les étiquettes des mots, ni leurs résultats mathématiques.

Prévoir environ 3,15 Go pour les coordonnées extraites, plus les résultats et la copie de travail. Éviter de lancer plusieurs gros calculs simultanément. La durée varie selon le matériel.

## Régénérations plus fortes

Pour régénérer toutes les coordonnées de poids quatre depuis les formules du moteur, plutôt que les lire dans les archives :

```sh
python checks/certificates/deliverable/run_f4_verification.py --regenerate --output ../f4-regenerated-I14
```

Cette voie utilise encore la base S3 archivée. Sa régénération séparée se fait dans une copie du répertoire `calculs-sources/generation-forte`, en y ajoutant `checks/certificates/deliverable/rerun_w3.py`, puis en lançant `python rerun_w3.py`. Le programme compare les objets exacts avant/après et écrit `certify_w3_validated.json`. Il écrase sa copie de travail de la base, d'où la nécessité d'une copie séparée. Le reçu de cette livraison dit explicitement si cette étape a été exécutée.

La régénération complète peut durer plus d'une heure. Les contrôles ciblés du moteur ne s'y substituent pas. Les essais historiques utilisant SciPy, Singular, Sage ou les données extérieures de Riesen ne sont pas des dépendances des commandes de démarrage ci-dessus.

## Résultats et provenance

Consulter [STATUS.md](STATUS.md) et `checks/release-I-v14/` pour les exécutions nouvelles, et `provenance/` pour les reçus F4 originaux. Le nouvel inventaire F4 restreint les sources distribuées ; les sources de calcul retenues sont inchangées. Les journaux anciens présents ailleurs restent historiques et ne prouvent pas une réexécution le 7 novembre.

Le PDF de I inclus est le PDF public téléchargé depuis Hexagon ; il n'a pas été recompilé pour cette livraison. La compilation du manuscrit ne fait pas partie de la validation des calculs.
