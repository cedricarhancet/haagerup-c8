# Reproduction

Exécuter les commandes depuis `II/`, dans une copie propre de la release ou dans l'archive autonome de II. Python 3.12.14 est la version utilisée ; les versions des dépendances sont fixées dans `requirements-exact.txt`. Un environnement virtuel permet de les installer sans modifier les autres projets. Ne pas utiliser `-O`, `-OO` ou `PYTHONOPTIMIZE` : plusieurs vérifications originales emploient des assertions. Les imports locaux des scripts de `R3-verif` demandent son répertoire courant ; le lanceur le fixe explicitement et utilise `-B`, pas `-I`.

```sh
python verify_manifest.py
python -m pip install -r requirements-exact.txt
python reproduce.py --group smoke --output ../../II-smoke
python reproduce.py --group certificate --output ../../II-certificate
python reproduce.py --group modes --output ../../II-modes
python reproduce.py --group originals --output ../../II-originals
python reproduce.py --group auxiliary --output ../../II-auxiliary
```

Chaque groupe demande un dossier nouveau ou vide, extérieur à `II`. Les commandes s'appliquent aussi sous PowerShell avec le `python.exe` de l'environnement choisi. Aucun téléchargement n'a lieu pendant les calculs. Les chemins sont relatifs aux fichiers livrés ; aucune dépendance à un répertoire de l'auteur n'est nécessaire.

| Groupe | Entrées et dépendances | Opérations et portée |
| --- | --- | --- |
| `smoke` | bibliothèque standard | identités exactes de `controles-categoriels/categories-exact-check.py` ; démarrage |
| `certificate` | certificat JSON, SymPy | `smoke`, puis rang des traces, paramétrisation, quadriques, égalité d'idéaux, quatre points et tests PCT ; les tenseurs sont fournis |
| `modes` | certificat JSON, `data/t12_W1.json`, python-flint | reconstruction indépendante des racines, cocycle, triplets, caractères, noyau et tenseurs depuis les modes de courants ; aucune importation des deux moteurs initiaux |
| `originals` | source I fixée, tables, données R3 et R3-verif ; toutes les dépendances fixées | réextraction de T12 ; cœur des deux moteurs initiaux : noyaux, produits, solutions, PCT ; réexport du certificat et comparaison JSON. Les espaces S2full archivés sont lus pour les contrôles supplémentaires de `r3b_verify.py` |
| `extended` | tables et données S2 ; dépendances fixées | régénérations S2 et S2full, fermeture exhaustive du candidat 0, fermetures échantillonnées des autres, vérifications additionnelles et forme rationnelle ; diagnostic long facultatif |
| `auxiliary` | tables, données archivées, SymPy, NumPy, FLINT | identités symboliques, contrôles arithmétiques bornés, énumérations de plongements et double contrôle de Gröbner ; des contrôles de §6 restent numériques ou limités à une plage finie |
| `all` | entrées de `certificate`, `modes`, `originals`, `auxiliary` | ensemble des contrôles du parcours principal dans une même copie ; ne comprend pas `extended` et ne transforme pas les tests bornés en preuves générales |

`--seconds-per-command 1800` est la limite par défaut. Une interruption ou un dépassement produit un reçu incomplet ou en échec, pas une validation. Les calculs initiaux peuvent durer plusieurs dizaines de minutes selon la machine, en particulier les noyaux des espaces S2. Les anciens temps Linux des journaux ne prédisent pas les temps Windows.

Le lanceur copie `Calculs-HU-v08` dans `<SORTIE>/work`, enregistre stdout et stderr dans `<SORTIE>/logs`, puis écrit `RECEIPT.json` avec commandes, versions, durées et codes de sortie. Il compare les JSON régénérés des tables et du PCT, ainsi que le certificat réexporté, aux données distribuées, et vérifie que les entrées de la livraison n'ont pas changé. Les données `.pkl` initiales sont des fichiers intermédiaires ; leur concordance et les différences éventuelles d'encodage sont décrites dans STATUS. Le succès d'un ancien script qui imprime des tests n'est pas interprété seul comme une preuve : les résultats attendus et les vérifications indépendantes restent nécessaires.

## Commandes élémentaires

Depuis `II/Calculs-HU-v08` :

```sh
python -B verify_omega_certificate.py omega_certificate.json
python -B verification-independante/independent_tensors.py .
python -B controles-categoriels/categories-exact-check.py
python -B export_omega_certificate.py R3 omega_certificate-reexport.json
```

La dernière commande écrit son JSON à l'endroit choisi et réimporte le premier moteur ; elle n'est donc pas une vérification indépendante. Utiliser une copie de travail pour les programmes qui écrivent près d'eux. Depuis la racine `II`, l'extraction de la table se lance avec `data/t12_table.py` et la source de `inputs-I`, en fixant le répertoire courant à `Calculs-HU-v08/data`.

## Hypothèses importées

La réextraction confirme la correspondance des tables fournies avec les données imprimées de I v14 pour le cas (II), avec les hauteurs codées par l'extracteur. Elle ne redémontre ni les modules de I ni sa C2-cofinitude. Les 25 équations sont les contraintes nécessaires de II, lemme 7.15 ; elles ne constituent pas à elles seules une construction de l'extension diagonale. La comparaison des traces et traces de produits ne compare pas des matrices écrites dans une même base de plus haut poids.

Le corpus fini ne prouve pas les hypothèses analytiques de reconstruction de III, ni les axiomes VOA pour tous les poids. Les grands certificats F4 de I et les certificats de Zhu de III sont des livraisons séparées.
