# TP DIVERSITE GENETIQUE ET DIFFERENCIATION POPULATIONNELLE

(source https://adegenet.r-forge.r-project.org/files/PRstats/practical-MVAintro.1.0.pdf, https://popgen.nescent.org/startMicrosatellite.html)

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
```r
install.packages("adegenet","hierfstat","pegas")
install.packages("ggplot2")
```
```r
library("adegenet")
library("pegas")
library("hierfstat")
library("ggplot2")
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

```r
inputfile="PATH_TO_YOUR FILE"

dataset<-read.structure(inputfile)

```
suivre les instructions intéractives
Vérifier que les données sont bien représentées dans l'objet genind

```r
dataset
```
```
/// GENIND OBJECT /////////

 // 169 individuals; 11 loci; 108 alleles; size: 115.7 Kb

 // Basic content
   @tab:  169 x 108 matrix of allele counts
   @loc.n.all: number of alleles per locus (range: 7-13)
   @loc.fac: locus factor for the 108 columns of @tab
   @all.names: list of allele names for each locus
   @ploidy: ploidy of each individual  (range: 2-2)
   @type:  codom
   @call: read.structure(file = input.file)

 // Optional content
   @pop: population of each individual (group size range: 1-1)

```

Explorer les différentes infos disponibles:

```r
dataset@loc.n.all
dataset@pop
```

Il manque l'info pop
 ```r
pop<-c(rep("M",14),rep("S",11),rep("P",51),rep("KPC",5),rep("K",50),rep("A",5),rep("T",11),rep("B",22))
pop(dataset)<-pop
dataset@pop


table(pop(dataset))
barplot(table(pop(dataset)), col=funky(17), las=3,
        xlab="Population", ylab="Sample size")
```

<img width="1068" alt="image" src="https://user-images.githubusercontent.com/20643860/215523567-f74f398f-d809-4174-a7ac-667260d69025.png">

## Analyses de diversité génétique

La fonction summary() permet d'obtenir les données de bases de la diversité génétique du jeu de données

```r
div<-summary(dataset)
```

```
// Number of individuals: 169
// Group sizes: 14 11 51 5 50 5 11 22
// Number of alleles per locus: 8 7 10 7 12 11 13 9 13 11 7
// Number of alleles per group: 51 42 84 41 73 41 62 49
// Percentage of missing data: 7.42 %
// Observed heterozygosity: 0.68 0.8 0.69 0.79 0.73 0.72 0.69 0.71 0.68 0.62 0.62
// Expected heterozygosity: 0.79 0.71 0.71 0.7 0.79 0.78 0.79 0.74 0.79 0.73 0.77

```

```r
table(pop(dataset))
barplot(table(pop(dataset)), col=funky(17), las=3,
        xlab="Population", ylab="Sample size")
```
 Avec le même principe essaye d'obtenir i) le nb d'allèles par locus, ii) le nb d'allèle par pop
 
 on peut tester l'ecart à HW de chaque locus
 
 ```r
 hw.test(dataset)
 ```
 
 
 ## statistiques de diversité génétique par population
 
 ```r
 dataset.hfstat <- genind2hierfstat(dataset)
basicstat <- basic.stats(dataset, diploid = TRUE, digits = 2) 
names(basicstat)
```
explorez les valeurs de cette fonction pour obtenir les données de diversité par population


un exemple

```populations<-seppop(dataset)
inbred_coef <- sapply (populations, inbreeding, res.type = "estimate") 
Fis_Bar <- sapply (inbred_coef, mean)
Fis<-as.data.frame(Fis_Bar)
set.seed(999)
c<-boot.ppfis (dat = dataset, nboot =1000, quant = c (0.0025,0.998), diploid = TRUE, dig=4)

```
faire un graphique
```
c<-c$fis.ci

ggplot(Fis, aes(x=rownames(Fis),y=Fis$Fis_Bar)) +        # ggplot2 plot with confidence intervals
  geom_point() +
  geom_errorbar(aes(ymin = c$ll, ymax = c$hl))
  
  ```
  ![image](https://user-images.githubusercontent.com/20643860/215540314-a4e38d43-75fa-44f3-a06d-7171bd0c0cc2.png)

  
 
 ## Analyses de différencation entre population
 
 
 
 
 
 
 
     
     


