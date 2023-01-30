#TP DIVERSITE GENETIQUE ET DIFFERENCIATION POPULATIONNELLE

Dans ce TP, nous allons explorer un jeu de données de 169 tigres génotypés à 11 locus microsatellites.
Dans un premier temps, nous allons déterminer la diversité génétique à l’intérieur de chaque population et la divergence génétique entre les 9 populations.
Puis nous allons estimer le nombre de populations le plus vraisemblable parmi l’ensemble des génotypes.

Ce jeu de données est issu de l'article suivant
https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0111207#s5


1. Calculer les indices de diversité génétique de bases:

  - Polymorphisme/fréquence allèliques
  - Diversité génetique (hétérozygotie  attendue et observée)
  - Test d'écart à l'équilibre de Hardy Weinberg
  
2. Caractériser la diversité génétique au sein des populations et la différenciation
  - diversité génétique
  - Fis et   Fst global
  - Fst par paires de populations
  
3. 
  
 
 

Nous allons effectuer ces analyses sous R. 
Installez les packages suivants:
```
install.packages("adegenet","hierfstat","pegas")
```

