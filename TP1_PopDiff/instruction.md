# TP DIVERSITE GENETIQUE ET DIFFERENCIATION POPULATIONNELLE

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
  
 
 ## Ressources et packages nécessaires

Nous allons effectuer ces analyses sous R. 
Installez les packages suivants:
```
install.packages("adegenet","hierfstat","pegas")
```
```
library("adegenet")
library("pegas")
library("hierfstat")
```

## Données

Nous allons analyser une jeu de données de microsatellites des 168 individus issus de 8 populations et génotypés à 10 microsatellites.
Nous allons importer une fichier de format 'structure' en tant qu'objet "genind"

Vue du fichier structure:

```
Pati01	Pati09	Fca304-tailed	Fca441	HDZ700	F85	Fca954	F124	Pati15	F53	Pati18
D1598	1	200	123	129	149	137	176	180	262	211	144	217
D1598	1	203	129	129	157	144	176	192	270	214	149	217
D1599	1	188	114	127	149	137	-9	171	266	211	144	225
D1599	1	188	123	129	157	139	-9	171	270	211	144	225
D1622	1	188	120	129	149	144	155	171	266	211	-9	213
D1622	1	209	123	129	157	146	164	180	266	214	-9	225
```

```
inputfile="PATH_TO_YOUR FILE"

dataset<-read.structure(inputfile)

```
follow the interactive instructions



