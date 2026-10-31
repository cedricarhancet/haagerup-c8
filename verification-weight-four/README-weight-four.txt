POIDS QUATRE — CONTRÔLE DE RANG ET PORTÉE LOGIQUE — 31 octobre 2026

1. Deux vérifications à distinguer

Les anciens calculs exacts de 9088 produits et de 18194 mots ne sont pas
réexécutés ici. Les grands blocs qu'ils utilisaient n'avaient pas été déposés.
W4-audit-completude.txt et coverage-audit.json expliquent cette limite et
comptent aussi les entrées nulles omises par les anciens compteurs.

La nouvelle vérification régénère, modulo p = 8388571, le sous-espace N4'
engendré par les cinq familles de produits normaux. Elle donne une borne
inférieure de dimension en caractéristique zéro. Sous une borne supérieure
dim S4 <= 6010, cette borne remplace les anciens tests d'appartenance.

La nouvelle exécution a abouti en 525 secondes : rang 956 au poids trois,
puis rang 6010 pour N4', réparti en 131 blocs non nuls. Le vérificateur séparé
check_n4_rank_blocks.py a validé les 131 mineurs identités, de rang total 6010.
Les trois tests supplémentaires e(-1)e, x(-1)x et e(-1)x donnent un accroissement
de rang nul modulo p ; ils ne sont pas utilisés comme tests d'appartenance
en caractéristique zéro.

2. Justification du relèvement

Les coefficients sont réduits dans le corps résiduel de la place choisie de
K = Q(sqrt(13),sqrt(-3)), avec sqrt(13) -> 1649360 et sqrt(-3) -> 2036966.
Les inverses de dénominateurs sont calculés par pow(d,-1,p), qui échoue si
d est nul modulo p. Les seize coefficients des huit droites primaires ont
été comparés à leur expression exacte dans la base 1,sqrt(13),sqrt(-3),
sqrt(-39) : check_n4_inputs.json enregistre ce contrôle.

Le sous-espace de poids trois est généré par des produits de vecteurs de
U1 et U2, puis par les modes zéro des courants. Chaque opération linéaire
et chaque pivot non nul modulo p se relève dans l'anneau local. Le rang 956
obtenu au poids trois est donc une borne inférieure en caractéristique zéro.

Les générateurs de N4' sont les produits L(-1)S3, U1(-1)S3, U1(-2)U2,
U1(-3)U1 et m(-1)U2 pour les huit plus hauts poids primaires. La clôture
par les modes zéro des courants ajoute des éléments des mêmes familles
stables. Chaque vecteur réduit obtenu se relève donc dans le N4' défini
dans le manuscrit. Même une étape partielle de cette clôture ne pourrait
qu'affaiblir la borne inférieure, pas l'invalider.

Les tableaux NumPy utilisent des entiers signés de 64 bits. La dimension
ambiante maximale est dim Vprime_4 = 71052, et
71052*(8388571-1)^2 = 4999794713333434800 < 2^63-1.
Cette borne couvre les accumulations des produits matriciels. La formation
des opérateurs comporte au plus 201 termes par générateur, avant réduction
des sommes de doublons. Aucun arrondi flottant n'intervient dans les rangs.

3. Conséquence mathématique et limite

Si S est engendré par les poids un et deux d'une véritable algèbre de
Haagerup U, alors N4' est contenu dans S4, lui-même contenu dans U4, et
dim U4 = 6010 par le caractère. Une borne de rang 6010 pour N4' force
N4' = S4 = U4. Les produits U2(0)S3 et U2(-1)U2 appartiennent alors à N4'
par leur degré, sans résoudre séparément les anciens systèmes exacts.

Le même argument s'applique après l'identification avec l'algèbre de Jiao,
si son résultat d'existence et son caractère sont accordés. L'identification
utilise la clôture de poids deux et le rang de génération i6, sans la
génération forte ; il n'y a donc pas de circularité.

En revanche, le rang réduit6010 ne prouve pas à lui seul dim S4 <= 6010
pour un candidat abstrait sans cette borne supérieure. Le présent dossier
ne doit pas être cité comme une preuve inconditionnelle de cette borne.

4. Fichiers et vérificateurs

- N4_rank_run.out : journal de la nouvelle exécution.
- N4_rank_summary.json : rangs par poids de Cartan.
- N4_rank_blocks.npz : matrices réduites exportées, comprimées sans pickle.
- N4_rank_block_names.json : correspondance blocs/poids de Cartan.
- check_n4_rank_blocks.py et .json : vérification indépendante des mineurs
  identités des matrices exportées et de la somme des rangs. Ce contrôle
  porte sur l'indépendance ; l'appartenance à N4' vient du script de génération.
- closure_w4_rank_export.py : source historique du calcul de rang, avec
  export des matrices et correction du vocabulaire de réduction résiduelle.
- check_n4_inputs.py et .json : normalisations exactes et borne int64.
- strict_pickle_stream.py et reader-checks.json : lecteur refusant fichiers
  absents, troncatures, labels manquants et doublons incompatibles, destiné
  à une éventuelle reprise des anciens certificats. Ce lecteur n'est pas
  une preuve des produits nuls non fournis.

5. Relance

Depuis la racine du dossier livré :

  python3 -m pip install -r verification-weight-four/requirements-weight-four.txt
  python3 verification-weight-four/run_n4_rank.py

Le lanceur vérifie les empreintes des sources de generation-forte, utilise
une copie temporaire, fixe MODE=C2 et PSTART=8388608, puis vérifie les matrices
exportées. Les résultats sont écrits dans verification-weight-four/rerun/.
L'option --output permet de choisir un autre dossier. Les sources historiques
et les résultats enregistrés de cette révision ne sont pas écrasés.

Pour vérifier uniquement le rang des matrices enregistrées :

  python3 verification-weight-four/check_n4_rank_blocks.py

Les échecs de fichier, de format, de normalisation, de rang ou de mineur
provoquent une exception ; le vérificateur ne conclut pas par défaut à un
succès lorsque des fichiers manquent.
