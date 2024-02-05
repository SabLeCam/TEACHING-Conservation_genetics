# TP DIVERSITE GENETIQUE ET DIFFERENCIATION POPULATIONNELLE

(source https://adegenet.r-forge.r-project.org/files/PRstats/practical-MVAintro.1.0.pdf, https://popgen.nescent.org/startMicrosatellite.html)

Dans ce TP, nous allons explorer un jeu de données de 169 tigres génotypés à 11 locus microsatellites.
Dans un premier temps, nous allons déterminer la diversité génétique à l’intérieur de chaque population et la divergence génétique entre les 9 populations.
Puis nous allons estimer le nombre de populations le plus vraisemblable parmi l’ensemble des génotypes en inférant le niveau d'admixture.

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
  
3. Inférence bayésienne de la structure des populations
  
 
 ## Ressources et packages nécessaires

Nous allons effectuer ces analyses sous R. 
Installez les packages suivants:
```r
install.packages("adegenet","hierfstat","pegas")
install.packages("ggplot2","magrittr","reshape2")
install.packages(c("fields","RColorBrewer","mapplots"))

source("http://bioconductor.org/biocLite.R")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.16")
`%>%` <- magrittr::`%>%`

source("http://membres-timc.imag.fr/Olivier.Francois/POPSutilities.R")

source("http://membres-timc.imag.fr/Olivier.Francois/Conversion.R")

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("LEA")

```
```r
library("adegenet")
library("pegas")
library("hierfstat")
library("ggplot2")
library("dplyr")
library("LEA")
library("reshape2")
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
inputfile="PATH_TO_YOUR FILE" #169_tigers_structure.stru

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

```
pop1 Melghat
pop2 Satpura
pop3 Pench
pop4 Kanha-Pench corridor
pop5 Kanha
pop6 Achanakmar 
pop7 Tadoba
pop8 Bandhavgarh
```

 ```r
pop<-c(rep("M",14),rep("S",11),rep("P",51),rep("KPC",5),rep("K",50),rep("A",5),rep("T",11),rep("B",22))
pop(dataset)<-pop
dataset@pop
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
On peut récupérer plusieurs informations de ce résumé comme par exemple le nombre d'individus par pop

```r
barplot(div$n.by.pop, col=funky(17), las=3,
        xlab="Population", ylab="Sample size")
```
 Avec le même principe essaye d'obtenir i) le nb d'allèles par locus, ii) le nb d'allèle par pop
 
 on peut tester l'ecart à HW de chaque locus
 
 ```r
 hw.test(dataset)
 ```
 et teste si He et Ho sont signicativement différents (ecart à HW)
 
 ```
 bartlett.test(list(div$Hexp, div$Hobs))
 t.test(div$Hexp, div$Hobs, pair = T, var.equal = TRUE, alter = "greater")
 ```
 
 
 ## statistiques de diversité génétique par population
 
 ```r
 dataset.hfstat <- genind2hierfstat(dataset)
basicstat <- basic.stats(dataset, diploid = TRUE, digits = 2) 
names(basicstat)
```
explorez les valeurs de cette fonction pour obtenir les données de diversité par population


Estime et test le coefficient de consanguinité Fis

```
populations<-seppop(dataset)
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

  
 
 ## Analyses de différenciation entre population
 
 On compare maintenant les patrons de diversité génétiques entre les populations par rapport à la diversité globale.
 
Matrice de Fst par paire de population (estimateur du Fst de Weir de Cockerham(1984)
```r
mat.obs <- pairwise.WCfst(dataset.hfstat)
```
Représenter cette matice avec une heatmap

```r
melted_matobs <- melt(mat.obs, na.rm = TRUE)
colnames(melted_matobs)<-c("pop1","pop2","value")

ggheatmap <- ggplot2::ggplot(melted_matobs, aes(pop1, pop2, fill = value)) +
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red",  
                       midpoint = 0.15, limit = c(0,0.3), space = "Lab" ) +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))
ggheatmap
```

![image](https://user-images.githubusercontent.com/20643860/215585630-23152f8e-0602-4bd4-af86-d89fc9752bf1.png)


## Inférence bayésienne de la structure de population

```r
#input les données, fichier avec hearder différent

input.file2 = "PATH_TO_YOUR FILE" #169_tigers_structure.stru.txt"

tiger.geno<-struct2geno(input.file=input.file2, 2, FORMAT = 2,
            extra.row = 1, extra.col = 2)
```
            
On choisit le nombre de K à tester (nb de populations à inférer) et le nombre de répetition
```r
obj.snmf = snmf("tiger.geno", K=1:8, ploidy=2, entropy=T, repetition=30, alpha=100, project="new")
plot(obj.snmf, cex = 1.2, col = "lightblue", pch = 19)
```

Comment choisir le modèle avec la meilleur probabilité postérieur?

"We use SNMF’s cross-entropy criterion to infer the best estimate of K. The lower the cross-entropy, the better our model accounts for population structure. Sometimes cross-entropy continues to decline, so we might choose K where cross entropy first decreases the most."


On choisit le meilleur run pour le K considéré

```r
ce <-  cross.entropy(obj.snmf, K = 5)
ce
best_run <- which.min(ce)
best_run
```

Comment représenter les données?

```r
qmatrix = Q(obj.snmf, K = 5, run=best_run)

par(mar=c(4,4,0.5,0.5))
pops<-levels(dataset@pop)
barplot(t(qmatrix), col=RColorBrewer::brewer.pal(9,"Paired"), 
        border=NA, space=0, xlab="Individuals", 
        ylab="Admixture coefficients")
segments(x0 =14, y0 = 0, x1 = 14, y1 = 1, col = "black", lwd=3)
segments(x0 =25, y0 = 0, x1 = 25, y1 = 1, col = "black", lwd=3)
segments(x0 =76, y0 = 0, x1 = 76, y1 = 1, col = "black", lwd=3)
segments(x0 =81, y0 = 0, x1 = 81, y1 = 1, col = "black", lwd=3)
segments(x0 =131, y0 = 0, x1 = 131, y1 = 1, col = "black", lwd=3)
segments(x0 =136, y0 = 0, x1 = 136, y1 = 1, col = "black", lwd=3)
segments(x0 =147, y0 = 0, x1 = 147, y1 = 1, col = "black", lwd=3)
#Add population labels to the axis:
for (i in 1:length(pops)){
  axis(1, at=median(which(dataset@pop==pops[i])), labels=pops[i])}
  ```
  
  
  ![image](https://user-images.githubusercontent.com/20643860/215650112-ae3fb189-1511-4196-bd7e-2a6eda3db8ee.png)

  ## Compte rendu de TP:
  

1.	Une brève introduction (environ une page à 1.5 interligne ) sur le contexte problématique des petites populations isolées de tigres (dérive génétique, consanguinité) et l’importance de maintenir des corridors pour permettre l’échange d’individus.
2.	L’objectif du TP
3.	Indiquez si les populations sont en HW ou si elles ont un déficit en hétérozygotes
4.	Indiquez les Fst entre les différentes populations
5.	Quelles populations sont les plus distantes (en Fst) ?
6.	Est-ce que les populations ont de bons niveaux de  diversité génétique? 
7.	Indiquez le nombre de populations le plus vraisemblables (à l’aide de structure).
8.	Commentez brièvement ces résultats

  









 
 
 
 
 
 
 
 
 
 
 
 
     
     


