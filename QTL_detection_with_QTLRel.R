
			#--------------------------------------------------------------------------------------------------
			#   
			#		SCRIPT R : QTL DETECTION WITH QTL REL
			#
			#					Holtz Yan
			#---------------------------------------------------------------------------------------------------



#-- OBJECTIF 
	# Ce script permet de calculer le LOD de chaque marqueurs pour chaque variables phéno données.
	# Il travaille pour un chromosome donné, et appel QTLRel
	# Il sort un fichier nommé bilan_simple_marker qui servira d'entrée à l'application shiny de visualisation des QTL.

#-- FICHIER INPUT DANS L'ORDRE : 
	# Fichier de génotypage :
	# Fichier de phénotypage :
	# Carte génétique
	# Chromosome d étude
	

# -- Récupération des Arguments
#args <- commandArgs(trailingOnly = TRUE)
#fic_geno=args[1]
#fic_pheno=args[2] 
#fic_map=args[3]

fic_geno="/NAS/davem_data/EPO/SAUVEGARDE_DATA_YAN/DIC2_BYBLOS/CAPTURE_ALL_INDIV_7_2_2015/MAPPING_ON_EPO/IMPUTATION/fichier_genotypage_QTL.csv"
fic_pheno="/NAS/davem_data/EPO/SAUVEGARDE_DATA_YAN/DIC2_BYBLOS/CAPTURE_ALL_INDIV_7_2_2015/PHENOTYPAGE/bilan_pheno_Dic2_Byblos.csv"
fic_map="/NAS/davem_data/EPO/SAUVEGARDE_DATA_YAN/CARTOGRAPHIE_GENETIQUE/5_CARTE_DB/map_avec_posi_physique.txt"
chromo="1A"

# Package / libraries nécessaires :
#library(devtools)
#install_github("timflutre/rutilstimflutre") Attention, lme4 doit etre installé auparavant.
#install.packages("QTLRel")
#library(rutilstimflutre)
library(QTLRel)
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#






#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

# -------------------------------------------------
# PARTIE 1 : RECUPERATION DES 3 TABLEAUX D'ENTREE
# -------------------------------------------------


# --- Génotypage
# Le format du génotpage doit respecter les règles suivantes : 1/ Nom du fichier = fic_genotypage.csv 2/ en tete des colonnes = nom des geno 3/ données manquantes = "-" 4/ Séparateur=";" 
geno <- read.table(fic_geno, sep = ";" , header = F, na.strings = "-")
geno=as.matrix(geno)
colnames(geno)=geno[1,]
geno=as.data.frame(geno[-1 , ])
print("--- Your genotyping matrix looks correct. Dimension of the matrix are :")
print(dim(geno))

# --- Carte 
# Format de la carte : 3 colonnes : LG, nom du marqueur, position dans le LG
map <- read.table(fic_map , header=T , dec = ".", na.strings = "-" , check.names=F)
colnames(map) <- c("LG", "marqueur", "Distance","group_physique","Posi_physique")
rownames(map) <- map$marqueur
map$LG <- as.factor(map$LG)
print("--- Your genetic map looks correct. Dimension of the map are :")
print(dim(map))

# --- Phénotypage
Y=read.table(file = fic_pheno, header = TRUE, sep = ";", dec = ".", na.strings = "NA")
colnames(Y)[1]="geno"
print("--- Your Phenotyping matrix looks correct. Dimension of the matrix are :")
print(dim(Y))






#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

# -------------------------------------------------
# PARTIE 2 : VERIFICATION DES COMPATIBILITES
# -------------------------------------------------


# Vérification des données pour voir si tout va bien
verif=as.data.frame(matrix(0,6,2)) ; verif[,1]=c("number of markers in the map" , "number of markers in the genotyping matrix" , "nbr indiv dans fichier de génot" , "nbr d'indiv dans fic phénot" , "nbr geno communs carte / genotypage","nbr indiv communs genot / phenot")
verif[1,2]=nrow(map) ; verif[2,2]=ncol(geno)-1 ; verif[3,2]=nrow(geno) ; verif[4,2]=nrow(Y) ; verif[5,2]=length(colnames(geno)[colnames(geno)%in%map$marqueur==TRUE])  ; verif[6,2]=length(Y[,1][ Y[,1]%in%geno[,1]==TRUE])
print(verif)

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#











#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

# -------------------------------------------------
# PARTIE 3 : CHANGEMENT FORMAT POUR QTLREL + CREATION DE FICHIER
# -------------------------------------------------

#Geno
rownames(geno)=geno[,1] ; geno=geno[,-1] ; geno=as.matrix(geno)
geno[which(geno=="A")]<-"AA"
geno[which(geno=="B")]<-"BB"

#On va remplacer les données manquantes par un allele au hasard.
my_fun=function(x){length(x[x=="AA"])/length(!is.na(x)) }
prop=apply(geno , 2 , my_fun)
for(i in c(1:ncol(geno))){
	aa=geno[,i][is.na(geno[,i])]
	bb=rbinom(length(aa),1,prob=prop[i])
	geno[,i][is.na(geno[,i])]=c("BB","AA")[bb+1]
	}

#Phéno
rownames(Y)=Y[,1] ; Y=Y[,-1]

#Map
map=map[ , c(2,1,3)]



#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
# -------------------------------------------------
# PARTIE 3 : SELECTION PHENO + FABRICATION FICHIER CORESPONDANT
# -------------------------------------------------


# Il faut une matrice avec les memes indiv dans les matrices de phéno et de géno
y=Y[,3] ; y=as.data.frame(y) ; rownames(y)=rownames(Y) ; y=na.omit(y)
my_geno=geno[which(rownames(geno)%in%rownames(y)) , ]

#Matrice identité
I<-diag(nrow(y))

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#






#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

# -------------------------------------------------
# PARTIE 4 : PASSAGE DE QTLREL
# -------------------------------------------------

# Calcul matrice apparentement IBS
K<-genMatrix(my_geno)
# On visualise la matrice K comme ceci :
K$AA[1:10 , 1:10]
hist(K$AA[1:10 , 1:10] , breaks=20)


# le premier modele sert a calculer la variance
toto=as.matrix(y)
mod1<-estVC(y=toto,v=list(AA=K$AA,DD=NULL,HH=NULL,AD=NULL,MH=NULL,EE=I))
mod1

# Genet assoc
GWAS.mod1<-scanOne(y=toto,gdat=my_geno,vc=mod1,test="Chisq")

#Récupération des pvalues
pval.mod1<-data.frame(marqueurs=colnames(my_geno),pvalue=GWAS.mod1$p)

#Seuil de détection des LODs
seuil.mod1<-3

# Manhattan plot
pval.mod1=merge(pval.mod1 , map , by.x=1 , by.y=2 , all=T)
par(mfrow=c(4,4))
for(i in levels(pval.mod1$LG)){
	print(i)
	dat=pval.mod1[pval.mod1$LG==i , ]
	plot(-log10(dat$pvalue)~dat$Distance , ylim=c(0,5) )
	abline(h=seuil.mod1,col="red")
}


#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
