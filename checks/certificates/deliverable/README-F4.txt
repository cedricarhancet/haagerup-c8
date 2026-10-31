CERTIFICATS EXACTS DE POIDS QUATRE — VERSION 03

Ce dossier accompagne l'article Haagerup-c8-article-I-v03-arxiv.tex.
Le résultat d'ensemble est exclusivement celui de COMPLETE-F4.json.
Un fichier de journal ou un compte partiel ne constitue pas ce résultat.

RÉSULTAT LIVRÉ : PASS, zéro erreur. Les 131 blocs cartaniens sont contrôlés.
Les 18368 mots de familles comprennent 18194 valeurs non nulles et 174 zéros;
les 9800 produits comprennent 9088 valeurs non nulles et 712 zéros. Les 886
zéros sont tous revérifiés séparément (786 calculs directs et 100 exclusions
par l'énergie). Les deux parties exactes ont pris respectivement 392,7 s et
390,5 s en parallèle, puis l'agrégation exhaustive a terminé sans erreur.

CONTENU MATHÉMATIQUE

K = Q(sqrt(13),sqrt(-3)). Les 6010 pivots sont des mots des cinq familles
qui définissent N4'. Les données couvrent les 18368 mots de ces familles
et les 9800 produits de F4, y compris les produits nuls. Chaque pivot est
identifié par son mot source. Le rang est contrôlé au premier 8388571, avec
sqrt(13) -> 1649360 et sqrt(-3) -> 2036966, puis les appartenances sont
vérifiées exactement sur toutes les coordonnées après restriction des
scalaires à Q. Chaque dénominateur inversé modulo le premier doit être
inversible; aucun rang flottant n'intervient.

Le producteur calcule les mots par les formules du VOA de réseau. Le
vérificateur des coordonnées établit les identités B C = X et la couverture
des domaines; il ne transforme pas à lui seul une matrice sans provenance
en un mot de VOA. Les scripts de production, les étiquettes et les sources
sont donc fournis également. Tous les produits nuls font en outre l'objet
d'un second contrôle, par calcul direct des modes ou par la borne d'énergie
Cartanienne a^2/12 + b^2/52 > 4.

FORMAT ET RÉPARTITION DES FICHIERS

Les grandes données sont dans deux archives compagnons :
Certificats-poids-quatre-v03-partie-1.zip et
Certificats-poids-quatre-v03-partie-2.zip. Leur racine commune
est certificats-poids-quatre-v03/. Elle contient f4_pivots/, f4_families/ et
les dix f4_tests_*/, avec un manifeste COMPLETE.json par domaine. Les scripts,
le reçu final, les résultats par bloc, les contrôles négatifs et les journaux
utiles figurent dans le petit dossier de vérification de l'archive principale.
La liste exacte des archives et leurs empreintes est donnée par l'inventaire
livré avec le dossier principal.

Les coordonnées sont des dictionnaires Python sérialisés avec pickle, avec
des nombres rationnels gmpy2. Le vérificateur contrôle les empreintes
indiquées dans les manifestes avant de les lire.

RELANCE À PARTIR DES COORDONNÉES FOURNIES

Installer dans un environnement Python dédié les versions indiquées dans
requirements-exact.txt. Depuis la racine du dossier principal, adapter le
chemin du présent dossier et exécuter, par exemple :

  python3 checks/certificates/deliverable/run_f4_verification.py \
    --archive Certificats-poids-quatre-v03-partie-1.zip \
      Certificats-poids-quatre-v03-partie-2.zip \
    --output /tmp/haagerup-f4-replay

Si plusieurs volumes sont fournis, les donner tous après --archive.
L'option --sources permet d'indiquer explicitement le dossier
calculs-sources/generation-forte. Le lanceur vérifie les sources, crée une
copie isolée, vérifie les empreintes des deux volumes, extrait les coordonnées,
recalcule les identités exactes,
revérifie tous les zéros, puis exige le résultat de l'agrégateur.
Le dossier --output doit être nouveau. La procédure ne modifie pas les
certificats ni les sources livrés.

RELANCE EN REPARTANT DES FORMULES DE MODES

La variante suivante n'utilise aucun bloc de poids quatre enregistré :

  python3 checks/certificates/deliverable/run_f4_verification.py \
    --regenerate --output /tmp/haagerup-f4-regeneration

Elle réutilise la base de poids trois livrée, dont la régénération exacte
indépendante est enregistrée dans certify_w3_validated.json. Le script
rerun_w3.py et le producteur historique certify_w3.py permettent aussi de
refaire cette étape. Une régénération complète peut durer plus d'une heure.

Pour refaire le contrôle du réseau par une seconde méthode, exécuter :

  python3 checks/certificates/deliverable/check_lattice_ball.py \
    --sources calculs-sources/generation-forte \
    --output /tmp/check-lattice-ball.json

Le contrôleur utilise une décomposition LDL rationnelle et des bornes
entières exactes, sans la recherche par voisins racines utilisée par le
producteur. Il retrouve les 8821 vecteurs de norme au plus 8. Le fichier
central-pivots-agree.json constate également l'égalité exacte entre les
280 pivots centraux utilisés pour le calcul du commutant C4 et les 280
pivots du bloc final de poids cartanien (0,0).

PROTOCOLE D'INTÉGRITÉ

safe_store.py ferme chaque enregistrement et relit intégralement tous les
flux avant de publier un manifeste COMPLETE. Le lecteur final impose les
noms canoniques des blocs, les empreintes, les dimensions, les étiquettes,
l'absence de doublons, l'absence de troncature et l'absence d'octets ou de
fichiers supplémentaires. Le programme refuse Python -O/PYTHONOPTIMIZE.
L'agrégateur exige la totalité des parties, des blocs et des étiquettes,
ainsi que le reçu indépendant de tous les zéros. Des injections de fautes
sur petites données temporaires contrôlent ces chemins d'échec.

JOURNAL DE LA RECONSTRUCTION

Une première écriture de cette révision a produit des flux incomplets et
quelques records tronqués, détectés par la relecture stricte; leur cause
physique n'a pas été établie. Aucun résultat de vérification globale n'a
été déclaré à partir de ces flux. La reconstruction finale conserve les
records décodés et contrôlés, recalcule exactement toutes les étiquettes
absentes, reconstruit les pivots à partir de leurs mots de familles et
réévalue leur rang. Elle est effectuée dans un répertoire temporaire isolé.
Les nombres récupérés et les offsets des troncatures sont consignés dans
les nouveaux manifestes. Les flux incomplets ne sont pas livrés comme
certificats valides. Le producteur sûr permet de refaire directement une
production complète, sans reproduire cette étape de réparation.

La réparation a recalculé 4403 mots de familles et 60 produits manquants.
Les premiers domaines e et (3,3) ont été reconstruits par repair_f4_v01.py;
les autres par repair_f4.py, qui consigne aussi une éventuelle dernière
lecture tronquée. Cette correspondance est enregistrée dans SOURCE-MANIFEST.
Les deux versions effectivement exécutées sont conservées. Leur rôle est
documentaire; --regenerate utilise le producteur sûr et ne dépend d'aucun
flux ancien incomplet.

PORTÉE DES CONTRÔLES

Il s'agit d'une preuve assistée par calcul exact, avec sources et coordonnées
reproductibles. Les empreintes assurent l'identité des fichiers; elles ne
remplacent pas les vérifications mathématiques. Les tests par injection de
fautes concernent le protocole de lecture et de couverture, et non les
identités de VOA. Les formules de modes et l'implémentation de l'arithmétique
font partie du code soumis à la lecture; aucune formalisation dans un
assistant de preuve n'est revendiquée.
