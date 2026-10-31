VERIFICATIONS COMPLEMENTAIRES — ARTICLE I, VERSION 02 — 31 octobre 2026

Cette archive distingue les sources historiques téléchargées de Dropbox
(calculs-sources/) et les vérificateurs corrigés de ce dossier. Elle ne constitue
pas une certification de l'ensemble du manuscrit. Les sources historiques
restent identiques aux octets reçus ; MANIFEST-sources.json et
MANIFEST-sources.sha256 en fixent les empreintes.

1. Exclusions des candidats de Griess, cas (II)

check_primary33_exclusions.py (révision par l'agent chargé des générateurs)
utilise seulement l'espace primaire de type (3,3), de dimension trois. Pour
chacun des treize candidats autres que 9 et 10, il écrit les conditions
quadratiques d'orthogonalité, les mineurs exprimant la condition de vecteur
propre et la condition z P(u,u)-1. Les trois cartes projectives donnent chacune
l'idéal unité sur Q. Les 39 systèmes, matrices et sorties sont dans
primary33-exclusions.json. Le calcul a été exécuté pendant cette révision.

check_primary33_certificate.py recalcule indépendamment les 39 bases de
Gröbner à partir des seuls systèmes polynomiaux enregistrés, sans le moteur du
VOA. Sortie : check-primary33-certificate.json, PASS.

2. Les 256 choix de droites pour le candidat 9

Le script historique closure_certify.py évaluait d'abord les nombres
algébriques avec sp.N avant de créer des boules Arb. Une boule contenant le
nombre décimal arrondi ne suffit pas, à elle seule, à certifier un encadrement
du nombre algébrique exact.

closure_certify_exact_intervals.py corrige ce point : il évalue récursivement
les rationnels, sommes, produits, puissances entières et racines carrées dans
Arb/Acb, et refuse les coefficients flottants ainsi que les branches de racine
non décidées. À 256 bits, la relance certifie 62 paires interdites. Les dix
autres paires proposées par l'exploration numérique restent indécises et ne
sont pas utilisées. Chaque exclusion enregistrée comporte une forme linéaire
rationnelle, le produit testé et une boule ne contenant pas zéro.

Sorties : closure-certify-exact-intervals.json et .out.
check_closure_coverage.py parcourt ensuite les 256 choix et donne un témoin
pour chacun des 254 choix exclus. Il reste exactement deux choix, complémentaires.
Sortie : check-closure-coverage.json, PASS.

3. Identification avec les générateurs de Jiao

checks/check-jiao-normalization.py, situé au niveau supérieur, vérifie les
quatorze coefficients et leur normalisation par calcul exact dans une
extension multiquadratique. i5_membership_exact.py remplace la reconnaissance
rationnelle flottante de kvec par ce calcul exact.

La relance complète de l'appartenance donne les mêmes huit valeurs booléennes
que la sortie historique : les deux générateurs normalisés appartiennent à
U2 du candidat 9, choix 0, et les six autres tests sont négatifs.
La relance de i6_generation.py donne le rang initial 65 puis le rang 120 après
une itération, modulo 8388571, avec les racines résiduelles 1649360 et 2036966.

Une seconde exécution, i6_generation_certificate.py, enregistre les 120 vecteurs
exacts et un graphe de mots indiquant leur production à partir des vecteurs
initiaux. i6_rank_certificate.json contient ce témoin. Le vérificateur
check_i6_rank_certificate.py, qui n'utilise que la bibliothèque standard,
recalcule le rang 120, contrôle les dénominateurs et l'ordre des dépendances
du graphe. Il ne remplace pas le moteur du VOA pour la vérification de chaque
produit ; la relance complète produit et calcule ces produits.

4. Borne de rang dans E6

e6_rank_bound.py a été réexécuté : 243 droites candidates, aucune donnant un
multiensemble symétrique. Sortie : e6-rank-bound-rerun.out.

5. Relancer les vérifications sélectionnées

Depuis le dossier racine du dossier livré :

  python3 -m pip install -r verification-certificates/requirements-selected.txt
  python3 verification-certificates/run_selected_checks.py --check all

Le lanceur vérifie d'abord les empreintes des sources, travaille dans une
copie temporaire et écrit les nouveaux résultats dans
verification-certificates/rerun/. Il conserve les sources et les résultats
historiques. Options --check : primary, closure, jiao, e6, all. L'option
--output permet de choisir un autre dossier de résultats.

Les calculs ont été relancés avec SymPy 1.14.0, gmpy2 2.3.2 et python-flint
0.9.0. NumPy est également nécessaire au moteur. Aucun environnement Python
installé n'est inclus dans le dossier livré.

6. Limites

Le contrôle du moteur de produits, les caches de primaires, la totalité des
calculs de classement, de poids quatre, de Zhu, de descente et d'automorphismes
ne sont pas tous certifiés par les seuls tests ci-dessus. Les contrôles de
poids quatre sont traités séparément dans verification-w4/. Les anciennes
sorties de calcul ne sont pas présentées comme des relances nouvelles.

Le dossier calculs-sources est un instantané sélectionné, non une copie
exhaustive de tous les calculs du projet. SNAPSHOT-LIMITS.json précise les
sources complémentaires non transférées ; aucune de ces sources n'est
nécessaire aux quatre relances proposées ici. Les gros blocs de poids quatre
signalés par INDEX-scratch.txt ne figurent pas dans cet instantané.

7. Courant simple : mot explicite indépendant des blocs W4

Le contrôle ajouté check_simple_current_word.py reconstruit exactement
u_(1,7)(-1)u_(2,6), vérifie son coefficient non nul et ses 26 annulations.
Il n'utilise aucun fichier piv_w4. Pour le relancer depuis la racine :

  python3 verification-certificates/check_simple_current_word.py --source-dir calculs-sources/generation-forte --output verification-certificates/simple-current-word-rerun.json

La sortie de la relance effectuée pendant la révision est conservée dans
simple-current-word.json. Ce contrôle est séparé des options du lanceur
run_selected_checks.py.
