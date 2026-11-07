# Concordance des versions

| Objet | Version de référence | Livraison correspondante |
| --- | --- | --- |
| Article I | [Hexagon 2610.00733v2](https://hexagonmath.org/2610.00733v2), source I v14 | [article-I-v14-supplements-v01](https://github.com/cedricarhancet/haagerup-c8/releases/tag/article-I-v14-supplements-v01), commit `813df07d2cf2a2b433d13093ccb920031fee5f9a`. |
| Article II | [Hexagon 2611.00034v2](https://hexagonmath.org/2611.00034v2), référence scientifique II v14 ; v16 documentaire pour les liens | [article-II-v14-supplements-v03](https://github.com/cedricarhancet/haagerup-c8/releases/tag/article-II-v14-supplements-v03), commit `134ce1476c7d4ce9b336ea1d167d707870650e33`. La v03 conserve les mêmes données mathématiques que v02 et corrige deux champs de métadonnées. |
| Article III | [Hexagon 2611.00164v1](https://hexagonmath.org/2611.00164v1), source v10 confirmée par empreinte | [article-III-v10-supplements-v01](https://github.com/cedricarhancet/haagerup-c8/releases/tag/article-III-v10-supplements-v01), source/PDF v11 pour la documentation et les liens des suppléments. |
| Calculs A1/A4/C3 | Archive antérieure v10 | `Calculs-III-A1-A4-C3-v11.zip`, documentation et lanceurs révisés, nouveaux reçus ; données conservées et adaptations de portabilité identifiées dans `SOURCE-PROVENANCE.json`. |
| Certificats Zhu | Archive antérieure v06 | `Certificats-Zhu-III-v07.zip`, données mathématiques d'origine conservées, nouveaux contrôles et documentation. |
| Paquet arXiv III | Ancien paquet v10 fondé sur un source antérieur de cinq octets | `Depot-arXiv-III-v11.zip` contient exactement le source v11 et les deux archives publiées, dans `anc/`. |

La source III v10 téléchargée sur Hexagon est identique octet par octet au source local v10, SHA-256 `4d60e712e926d4aae08f87c81fe8adcb3697042b67fe58c05bd01f9f3d809e1e`. L'ancien paquet de dépôt v10 avait encore `\cal A_S\ot\cal A_S` au lieu de `\cal A_S\otimes \cal A_S` dans le résumé. Le source v11 reprend la version diffusée, actualise le §10, l'appendice A.8 et la disponibilité des données, et retire un renvoi de page devenu obsolète dans A.3, sans changer les énoncés.

La référence Hexagon ci-dessus correspond au source scientifique v10 ; la version documentaire v11 est fournie par cette release GitHub. Le paquet arXiv accompagne cette livraison en vue d'un dépôt ultérieur.

Les sources TeX de I et III nécessaires à la comparaison de tables sont des entrées externes figées par SHA-256 dans `DEPENDENCIES.json` ; le lanceur du dépôt les télécharge depuis leurs releases. Elles ne sont pas enfouies dans les archives de calcul destinées au dossier auxiliaire arXiv.
