
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
gl.report.callrate(snp)
```


<img width="635" alt="image" src="https://user-images.githubusercontent.com/20643860/232498992-e3b76569-01a0-40d1-bfa1-4f8043520900.png">

## Filtrage des données
Ici on commence l'étape primordial de filtrage des données.
On commence par filtrer les SNPs avec que des données manquantes ou ceux qui sont monomorphes.
```r
###SNP filtering####
snp <- gl.filter.allna(snp)

snp2<-gl.recalc.metrics(snp)
snp3<-gl.filter.monomorphs(snp2)
```

Filtrer les snps en fonction de la maf (minimum allele frequency) et verifier leur neutralité (écart à l'équilibre de Hardy Weinberg)
```r
filtered_maf <- gl.filter.maf(snp3, threshold=0.05, verbose=3)
gl.report.heterozygosity(filtered_maf)#get new estimation of the genetic diversity in your data
hwe<-gl.filter.hwe(filtered_maf, subset="each")
gl.report.heterozygosity(filtered_maf)
```
A chaque étape, on voit le nombre de marqueurs qui sont éliminés car ils ne remplissent pas les critères.
Une fois fini l'étape de filtre sur la "qualité" des SNPs , filtrer les données pour réduire le nombre de données manquantes. Essayer de jouer sur les valeur seuil pour trouver le meilleur compromis entre le nombre de marqueurs et le % de données manquantes.
```r
result <- gl.filter.callrate(hwe, method='loc', threshold=0.8,
                             verbose=3)
result2 <- gl.filter.callrate(result, method='ind', threshold=0.8,
                             verbose=3)
                             ```
gl.report.heterozygosity(result2)
```
Voici les estimations de diversité génétique sur notre jeu de données filtré.

<img width="647" alt="image" src="https://user-images.githubusercontent.com/20643860/232503444-ac8ba3ff-1eda-4b9b-aff0-dd8acdae3d76.png">


## Analyse de la structure de la diversité génétique

Ici on utilise la méthode de l'analyse en composantes principales (ACP)
 ```r
## perform PCA
pca <- glPca(result2, nf=4)

score1<-pca$eig[1]/ sum(pca$eig)
score2<-pca$eig[2]/ sum(pca$eig)
score3<-pca$eig[3]/ sum(pca$eig)
print(list(score1,score2,score3))

pcaeig<-data.frame(eig=pca$eig[1:20],ID=seq(1,20,by=1))
#assign('pcaeig',pcaeig,envir=.GlobalEnv) 
#print(head(pcaeig))

eig<-pcaeig$eig
ID<-pcaeig$ID

eig_plot<-ggplot(data = pcaeig,aes(y=eig ,x=ID)) +
  theme(panel.background=element_blank(),
        axis.text=element_text(family="Arial Narrow", size=14,color='grey40'),
        axis.title=element_blank(),
        panel.border=element_rect(fill=NA,colour="grey40")) +
  geom_bar(stat="identity")

print(eig_plot)
```

Représenter les résultats de l'ACP
```r
li<-pca$scores
dfpca<-data.frame(a1=li[,1],a2=li[,2],pop=as.vector(result2$pop),ind=rownames(pca$scores))
#assign('dfpca',dfpca,envir=.GlobalEnv)
#print(head(dfpca))

centroids <- aggregate(cbind(a1,a2)~pop,data=dfpca,mean)
#print(head(centroids))
gpca<-ggplot(data = dfpca, aes(x=a1,y=a2, color=pop)) +
  theme(panel.background=element_blank(),
        axis.text=element_text(family="Arial Narrow", size=14,color='grey40'),
        axis.title=element_text(family="Arial Narrow", size=14,color='grey40'),
        panel.border=element_rect(fill=NA,colour="grey40")) +
  xlab(paste0('Axis',1))+
  ylab(paste0('Axis',2))+
  geom_hline(yintercept = 0, colour = "gray65") +
  geom_vline(xintercept = 0, colour = "gray65") +
  geom_point( position= "jitter", alpha = 0.8, size = 4) +
  geom_text(data=centroids,aes(label=pop), size=5,position="jitter")

gpca
```

<img width="645" alt="image" src="https://user-images.githubusercontent.com/20643860/232504704-95ad74b0-1db3-46b9-a9ba-9bb7eeecd6a8.png">

Quels pourcentage de la variance génétique dans notre jeu de données est expliqué par l'axe 1 et 2 (composantes principales 1 et 2)?       

Pour interpréter ces resultats, il est necessaire de visualiser la localisation des sites echantillonnés

```r
length(Ech_tot$ID)
i<-seq(1,913,by=1)
ind_lib<-subset(i, (i %in% g))

geoDf<-Ech_tot[ind_lib,]
names(geoDf)<-c("ID","lat","lon","Sex","Length_cm","Disk_width_cm", "Pop")


theWorld<-ggplot2::borders("world", colour="gray88", fill="gray88")
World<-ggplot(data = geoDf, aes(x = lon, y = lat)) +
  theWorld + theme_minimal() +
  geom_point( aes(x = lon, y = lat,color=factor(geoDf$Pop)), size=2, alpha = 1) + 
  theme(axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),  
        legend.position="none")

World



Europe<-ggplot(data = geoDf, aes(x = lon, y = lat,color=factor(geoDf$Pop))) +
  theWorld + theme_minimal() +
  geom_point( size=3, alpha = 1) +  
  theme(axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank() 
  )+
  scale_color_discrete(name="Sampling zone") +
  coord_cartesian(ylim=c(35, 54), xlim = c(-30, 10))#pour voir les açores
Europe
```


<img width="638" alt="image" src="https://user-images.githubusercontent.com/20643860/232504507-c8fc51be-bd0d-47ef-b9fe-d1f95bcdbe89.png">



