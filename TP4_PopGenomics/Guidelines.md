
# Etude de génomique des populations.

vous allez travailler sur un jeu de données issu d'un projet visant à étudier la diversité génétique de la Raie bouclée (*_Raja clavata_*) dans son aire de distribution.

Des données de RADseq ont été produites et vous travaillerez sur les données de génotypes de SNP issues de l'analyse bio-informatique de séquences RADseq.

Nous allons voir comment , au préalable de l'analyse de la diversité génétique en elle même, ces données méritent d'être explorées et filtrées afin de s'assurer de leur robustesse (taux de données manquantes, neutralité...)

Voici les packages R que nous allons utiliser dans ce TP. Les installer puis les appeler.
```r
library(vcfR)
library(pegas)
library(adegenet)
library(ggplot2)
library(plotrix)
library(dartR)
library(HardyWeinberg)
```

Ici je commence le TP avec un objet R directement car le fichier source est très volumineux (+3GB).
Voici les étapes que j'ai réalisé en amont:
- importer un fichier vcf
- le transformer en objet genlight (similaire à genind mais adapté au gros volume de données)
- vérifier que la conversion des genotypes est correcte.

<img width="1236" alt="image" src="https://user-images.githubusercontent.com/20643860/232494654-d4ff7b42-b333-4da8-83b6-ac9c19ac0698.png">


```r
#chargé l'objet snp
load(snp)
snp
```
<img width="476" alt="image" src="https://user-images.githubusercontent.com/20643860/232495027-76e20e44-3851-4ec2-9d3b-5b9bb4f6d209.png">

Explorez l'objet
```r
indNames(snp)
levels(pop(snp))

