# TP Définir des ESU à partir de réseaux d’haplotypes

(Source F. Dufresnes)

**Objectif**: Faire un réseau d’haplotype. Un réseau d’haplotypes permet d’explorer les liens qui existent entre les différents haplotypes présents dans une même espèce. 

Il existe différentes méthodes de construction des réseaux, nous allons explorer les réseaux d’haplotypes avec le logiciel *_PopArt_*.



## Contexte : 
Pour ce TP, nous nous inspirons de l’article de de Malgalhaes et al. 2017. 

<p float="left">
<img width="530" alt="image" src="https://user-images.githubusercontent.com/20643860/220001525-245ed199-fee1-45c6-9d3e-5d5bb0e2ade8.png">
<img width="239" alt="image" src="https://user-images.githubusercontent.com/20643860/220001561-99f6e98e-d6a5-4d13-8b17-dff30043a5af.png">
<p >



La diversité intraspécifique est corrélée au potentiel évolutif des espèces. Pourtant la diversité intraspécifique (génétique, phénotypique, etc) est souvent négligée dans les stratégies de conservation et dans la détermination des aires protégées.
Plusieurs espèces peuvent se subdiviser en ESU (Evolutionary Significant Unit) ou MU (Management. Units). La prise en compte de subdivisions intraspécifiques dans les stratégies de conservation pourrait faciliter la protection de la diversité génétique et du potentiel adaptatif des espèces.

Dans cet article, les auteurs s’intéressent à la grenouille leaf frog (*_Pithecopus ayeaye_*). Cette espèce endémique de l’écosystème *_campos rupestre_* du Brésil avait auparavant un statut d’espèce en danger critique d’extinction parce qu’elle ne se trouvait que dans deux localisations séparées et menacées par des pertes d’habitat. Toutefois la découverte de l’espèce dans de nouveaux sites fait en sorte qu’elle a perdu son statut d’espèce menacée et les efforts de conservation qui y sont associées. Les auteurs se demandent si cette décision est vraiment adaptée à la réalité de la population de Leaf frog ou si la considération de la diversité génétique, du potentiel évolutif et de la structure de population ne justifieraient pas qu’on maintienne un effort supplémentaire de conservation.

Les auteurs explorent notamment comment la population de Leaf frog peut se subdiviser en ESU. Ils concluent à la présence de 3 ESU (Pocos, Canastra et Quadrilatero). Ils cherchent ensuite à vérifier si les aires protégées actuelles permettent la protection de l’ensemble de ces ESU.

Nous allons explorer comment la diversité génétique se répartie pour deux loci de cette population. 



## Part I: faire les fichiers d’haplotypes pour PopArt
  
  *_Le fichier FASTA avec les séquences des haplotypes a été créé de la façon suivante :
Un fichier txt a été créé avec la liste des # d’accession Genbank des séquences (1 numéro par ligne). Les numéros étaient disponibles à la fin de l’article.
Dans l’outil Batch Entrez de NCBI Genbank, télécharger le fichier txt des numéros d’accession. Récupérer la liste d’accessions et l’exporter en fichier sous format FASTA.
Le fichier FASTA a été modifié au besoin pour être reconnu par le logiciel DNAsp (ex : Enlever, et ; / Y, R, W remplacés par -). Le logiciel MEGA peut être utilisé pour resauver le fichier en FASTA pour faciliter sa lecture._*
  
**Installation**
  - allez à http://www.ub.edu/dnasp/
  - cliquez sur Download DNAsp
  - installer le fichier dnasp51001.msi dans votre répertoire.
  - une fois l’installation completée, lancez le logiciel.

**Création du fichier pour le gène mitochondrial cytb**
  - Ouvrir **DNAsp**
  - Ouvrir le fichier cytb_renamed.fasta 
  - Pour le gène cytb, aller dans Data et changer les paramètres pour haploid et mitochondrial 
  - Aller dans Generate puis dans Haplotype data file. Une fenêtre va apparaitre concernant les paramètres pour générer un fichier d’haplotype. Vérifier que l’option de sortie cochée soit Nexus Haplotype data file.  
  - Cliquer sur ok.
  - Enregistrer (Ex : haplotype_cyb)
  <img width="378" alt="image" src="https://user-images.githubusercontent.com/20643860/220004835-13c1ffa2-c3bb-401e-8365-5fa84a1a5ed4.png">


  - Les résultats apparaissent dans une fenêtre . Une première partie résume votre analyse (nombre de séquences, nombre de site variables, les gaps sont-ils inclus ou non...). La partie haplotype distribution vous montre combien d’haplotypes vous avez dans vos données ainsi que leur fréquence. La diversité haplotypique est aussi donnée. Ensuite vous avez le nom des séquences qui composent chaque haplotype. 
  - Allez dans File puis Save/export data as. Choisissez le format nexus. Rentrer le nom de votre fichier avec un 2 à la fin (Ex : haplotype_cyb2).
  Vous venez de créer votre fichier d’haplotypes.

**Création du fichier pour le gène nucléaire POMC**
- Ouvrir DNAsp
- Pour ouvrir le fichier pomc_renamed.fasta, aller dans File>Open Unphase/genotype data file
- Faire Run pour reconstruire les haplotypes diploïdes.
- Aller dans Generate puis dans Haplotype data file et faire comme pour gène cytb.

**Ajout des localisations dans le fichier d’haplotypes**
- Ouvrir fichier d’haplotype (sans le 2) dans bloc note (ou autre logiciel pour lire les .txt)
- Ouvrir le fichier texte matrice_location_cytb (ou matrice_location_pomc).
- Copier le code et le coller juste avant le dernier bloc de code dans votre fichier d’haplotype (celui sans le 2).
- Enregistrer.

*_Le fichier d’haplotypes ne contient pas les données de géolocalisation. Mais les points GPS des locations sont fournies dans les suppléments de l’article et les localisations de chacun des séquences est fournie dans les suppléments (ou dans les informations de chaque séquence sur Genbank).
Une base de donnée avec les # de voucher, le site associé et les haplotypes associés pour chaque loci a été réalisé. Les informations ainsi compilées et les points GPS ont été utilisés pour ajouter une section au fichier d’haplotype qui permet d’associer à chaque localisation une fréquence de chaque haplotype à l’aide d’une matrice et d’associer chaque localisation à ses coordonnées. C’est cette partie de code que vous ajoutez au fichier._*

## Partie II : Réseaux.
  
  - Aller à : https://popart.maths.otago.ac.nz puis Download
  - Choisir votre système d’exploitation puis télécharger le fichier.
  - Installer le fichier dans votre répertoire. Cliquer sur popart.exe pour lancer l’installation (sous windows).
  - Ouvrir Popart.
  
  <img width="800" alt="image" src="https://user-images.githubusercontent.com/20643860/220002235-90be223e-7005-4c5b-8853-d57706f9bcd2.png">
  
  - Ouvrir votre fichier d’haplotypes (cyt b ou pomc)
  - Aller dans Network puis choisissez une méthode de construction (Minimum Spanning Network, epsilon=0)
  - Le réseau apparait dans la fenêtre à droite. 
  - Vous pouvez modifier les couleurs associées à chaque site en allant dans Edit > Set trait color. Choisissez un site, ok, puis choisissez sa couleur.

*_Note : Comme il y a plusieurs sites provenant de la même région, vous pourrez recolorer les haplotypes de ces sites avec la même couleur._*
  
  
  <img width="876" alt="image" src="https://user-images.githubusercontent.com/20643860/220008327-b545a361-d05c-40e2-99bb-b7ebae698fe2.png">
  
  - Explorer les résultats en coloriant d’une même couleur tous les sites d’une même ESU supposée par les auteurs (voir Table S1 dans suppléments de l’article)
  - Exporter votre réseau d’haplotype en png ou pdf en allant dans File>Export graphics


## Partie III : Réseaux en fonction du site d’échantillonnage.
  - Dans **View>Switch** to map view
  - Vous pouvez ainsi visualiser la répartition de vos haplotypes en fonction du site d’échantillonnage. Vous pouvez ajouter des légendes (échelle, rose des vents à votre carte en faisant clic droit puis Info boxes. 
  - Vous pouvez modifier les couleurs associées à chaque haplotype en cliquant sur les points de couleurs.
  - Vous pouvez colorier chaque regroupement d’haplotypes (basé sur réseau d’haplotype dans l’autre view) dans un registre d’une même couleur pour mieux visualiser comment les haplotypes les plus proches génétiquement se distribuent géographiquement.
  - Exporter votre réseau d’haplotype en png ou pdf en allant dans File>Export graphics

  <img width="1043" alt="image" src="https://user-images.githubusercontent.com/20643860/220010592-eb68eb36-2c3a-4cfb-bdae-dbb85e3644e4.png">


```
Votre rapport de TP devra répondre aux questions suivantes :
1)	Pour les deux loci, comment les groupes d’haplotypes se regroupent et quels sont les haplotypes les plus distants ?
2)	Quels sites possèdent la plus grande diversité haplotypique (considérez les deux loci) ?
3)	Est-ce que les patrons observés dans vos réseaux d’haplotypes semblent bien justifier les ESU définies par les auteurs de l’article? Justifiez.
4)	Selon vous, est-ce que les aires protégées actuelles(1) permettent de suffisamment protéger la diversité génétique de cette espèce? Justifiez. 
5)	Discutez des avantages/désavantages de l’ADN mitochondrial et de marqueurs nucléaires adaptatifs pour la caractérisation des unités évolutive significatives. Discutez de manière générale et faire des liens avec l’exemple du TP.
 ```

  *_Pour la question 4, on considère les informations suivantes (simplification des informations tirées de l’article) : -1 Le regroupement géographique de Poços n’est pas inclus dans une aire protégée. -2 Les populations regroupées dans les unités Canastra et Quadrilatero se retrouvent en partie dans une aire protégée. -3 Vous pouvez fouiller davantage ou vous en tenir à cette simplification._*

Les ESU définies par les auteurs sont les suivantes :

Canastra: Alpinopolis/Pedregulho/Sacramento/Sao Roque de Minas

Pocos: Pocos de Caldas

Quadrilatero: Lavras/Mindur/Nova Lima/Ouro Preto
  
  
  <img width="587" alt="image" src="https://user-images.githubusercontent.com/20643860/220010465-1fda5054-debb-45d1-8d99-c055404ca0a7.png">

  


