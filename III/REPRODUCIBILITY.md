# Reproduction depuis le dépôt

Python 3.12.14 a été utilisé pendant cette passe, sans `-O`, `-OO` ni `PYTHONOPTIMIZE`. Les calculs exacts ont été exécutés avec les versions épinglées ci-dessous. A1 et les contrôles de manifestes utilisent seulement la bibliothèque standard. A4 requiert gmpy2 ; les contrôles C3 et Zhu utilisent les bibliothèques indiquées par leurs scripts. Git Bash ou un shell Bash convient aux pilotes `.sh` facultatifs. Les lanceurs Python sont utilisables directement sous Windows, Linux et macOS ; seule la configuration Windows a été rejouée pendant cette passe.


Depuis `III/`, après installation de `requirements-exact.txt` :

```sh
python -I -X utf8 -B verify_manifest.py
python -I -X utf8 -B reproduce.py --group smoke --output ../../III-smoke
python -I -X utf8 -B reproduce.py --group A --output ../../III-A1-A4
python -I -X utf8 -B reproduce.py --group C3 --output ../../III-C3
python -I -X utf8 -B reproduce.py --group zhu-verify --output ../../III-Zhu-verify
python -I -X utf8 -B reproduce.py --group zhu-projection --output ../../III-Zhu-projection
python -I -X utf8 -B reproduce.py --group zhu-modes10 --output ../../III-Zhu-modes10
```

Chaque commande exige un nouveau dossier. Les dépendances et URL sont fixées dans `ASSETS.json` ; leur validation précède l'extraction et le chargement de données. Le lanceur utilise curl s'il est disponible, sinon la bibliothèque HTTP de Python. L'accès réseau n'est nécessaire que pour ces téléchargements. Un cache local préalablement vérifié peut être fourni par `--local-assets CHEMIN` ; les mêmes empreintes restent exigées.

| Groupe | Opérations et limites |
| --- | --- |
| `smoke` | Télécharge A1/A4/C3, puis relance `A1/a1_lattices.py` avec la bibliothèque standard. |
| `A` | Relance les dix commandes A1/A4, dont les six scripts totalisant 320 contrôles détaillés. |
| `C3` | Télécharge A1/A4/C3, Zhu, II et les sources I/III fixées ; relance le certificat portable, les catégories, les tables et la comparaison des 131 fichiers communs avec II. |
| `zhu-verify` | Confluence exhaustive et inclusion directe, inclusion inverse, fibre spéciale, poids et branchements. Calcul long, notamment la confluence. |
| `zhu-confluence`, `zhu-membership`, `zhu-fibre`, `zhu-weights` | Parties séparées du profil précédent, utiles pour une expertise ciblée. |
| `zhu-projection` | Régénère les 19 981 lignes à partir des 1 060 polynômes fournis et contrôle l'égalité exacte et modulo p. |
| `zhu-modes10` | Régénère la base abaissée et les dix premières lignes de produits, soit 1 155 paires. |
| `all` | Enchaîne A, C3, Zhu verify, projection et modes10 ; ne signifie pas régénération complète de tous les produits. |
| `zhu-regenerate` | Régénère les 120 lignes, noyaux, polynômes et projection ; `--rebuild-rules` ajoute la production des règles. Cette production complète n'a pas été rejouée pendant cette publication. |

Les groupes de modes utilisent un travailleur par défaut ; `--jobs` règle leur nombre. Les calculs longs peuvent prendre plusieurs heures et produire des données volumineuses. Le reçu `RECEIPT.json` du dossier de sortie contient l'état réel, les commandes, durées, codes de sortie et empreintes de journaux. Un échec conserve les journaux et ne produit pas un état PASS.

Pour télécharger seulement les dépendances sans calcul :

```sh
python -I -X utf8 -B fetch_supplements.py --output ../../III-packages
```

Dans chacune des deux archives extraites, `verify_manifest.py` contrôle les fichiers et `verify_results.py` contrôle les reçus publics sans rejouer les calculs longs. Les README et REPRODUCIBILITY de ces archives donnent aussi les commandes élémentaires et les chaînes C3 de régénération.
