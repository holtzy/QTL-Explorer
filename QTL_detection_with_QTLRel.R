
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
args <- commandArgs(trailingOnly = TRUE)
fic_geno=args[1]
fic_pheno=args[2] 
fic_map=args[3]

#setwd("/NAS/g2pop/HOLTZ_YAN_DATA/DIC2_SILUR/QTL/PUBLI")
#fic_geno="genotypage.csv" ; fic_pheno="phenotypage.csv" ; fic_map="carte"


# Package / libraries nécessaires :
library(QTLRel)
set.seed(123)
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#










#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

# -------------------------------------------------
# PARTIE 1 : RECUPERATION DES 3 TABLEAUX D'ENTREE
# -------------------------------------------------


# --- Génotypage
# Le format du génotpage doit respecter les règles suivantes : 1/ Nom du fichier = fic_genotypage.csv 2/ en tete des colonnes = nom des geno 3/ données manquantes = "-" 4/ Séparateur=";" 
genotype<-read.table(fic_geno, sep = ";" , header = F, na.strings = "-")
genotype=as.matrix(genotype)
colnames(genotype)=genotype[1,]
genotype=as.data.frame(genotype[-1 , ])
names(genotype)[1]<-"geno"
print("--- Your genotyping matrix looks correct. Dimension of the matrix are :")
print(dim(genotype))

# --- Carte 
# Format de la carte : 3 colonnes : LG, nom du marqueur, position dans le LG
map <- read.table(fic_map , header=T , dec = ".", na.strings = "-" , check.names=F)
colnames(map) <- c("LG", "marqueur", "Distance","group_physique","Posi_physique")
rownames(map) <- map$marqueur
map$LG <- as.factor(map$LG)
print("--- Your genetic map looks correct. Dimension of the map are :")
print(dim(map))

# --- Phénotypage
BLUP<-read.table(fic_pheno, header = TRUE, sep=";")
colnames(BLUP)[1]="geno"
print("--- Your Phenotyping matrix looks correct. Dimension of the matrix are :")
print(dim(BLUP))
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#






#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

# -------------------------------------------------
# PARTIE 2 : FUSION PHENO / GENO
# -------------------------------------------------

BLUP[,1]<-as.character(BLUP[,1])
genotype[,1]<-as.character(genotype[,1])
don<-merge(BLUP,genotype, by="geno")
print("--- Nombre d'individu communs entre pheno et géno :")
dim(don)[1]

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#





#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

# ---------------------------------------------------------------------------------------------------
# PARTIE 3 : CALCUL QTL REL POUR UNE VARIABLE DONNEE --> ENCAPSULE DANS UNE FONCTION
# ---------------------------------------------------------------------------------------------------

run_my_QTLREL=function(select){

	# Sous ensemble de don pour lequel la variable select n'est pas manquante !
	don2<-don[which(!(is.na(don[,select]))),]
	
	# Et Y va etre mon vecteur de phéno pour QTLRel
	Y<-don2[,select]
	
	# Et X va etre ma matrice de génotypage avec données manquantes:
	X<-don2[,(ncol(BLUP)+1):ncol(don2)]
	
	# Je vais remplacer les trous par des A et des B au hasard. Pour choisir l'un ou l'autre je vais choisir en faisant une binomiale avec comme proba la fréquence.
	XNNA=X
	my_fun=function(x){length(x[x=="A" & !is.na(x) ])/length(!is.na(x)) }
	prop=apply(XNNA , 2 , my_fun)
	for(i in c(1:ncol(XNNA))){
 		aa=XNNA[,i][is.na(XNNA[,i])]
 		bb=rbinom(length(aa),1,prob=prop[i])
		XNNA[,i][is.na(XNNA[,i])]=c("B","A")[bb+1]
 		}
	
	# Et je remplace par les "A" et "B" par des "AA" "BB" pour QTLREL
	XNNA=as.matrix(XNNA)
	XNNA[which(XNNA=="A")]<-"AA"
	XNNA[which(XNNA=="B")]<-"BB"
	
	# Maintenant je peux calculer ma matrice de Kinship:
	K<-genMatrix(XNNA)
	
	# I = matrice identité
	I<-diag(length(Y))
	
	# le premier modele sert a calculer la variance
	mod1<-estVC(y=Y,v=list(AA=K$AA,DD=NULL,HH=NULL,AD=NULL,MH=NULL,EE=I))
	mod1
	
	# Et ma détection de QTL !
	GWAS.mod1<-scanOne(y=Y,gdat=XNNA,vc=mod1,test="Chisq")


#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#






#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

# ---------------------------------------------------------------------------------------------------
# PARTIE 4 : LA FONCTION DOIT TOUT METTRE DANS UN TABLEAU BILAN
# ---------------------------------------------------------------------------------------------------
	
	
	# Get the mean value of each allele: 
	my_function=function(x, y) {
		moy <- aggregate(y, by = list(marqueur = x), mean, na.rm = TRUE)
		resu=data.frame(moy.A=moy$x[moy$marqueur == "A"],moy.B= moy$x[moy$marqueur =="B"] , a=abs(moy[1, 2] - moy[2, 2])/2)
		return(resu)
		}
	my_mean = apply(X, MARGIN = 2, FUN = function(x,y) my_function(x, y), y=Y)
	my_mean = matrix(as.vector(unlist(my_mean)) , ncol=3 , byrow=T)	
	
	# Get the p-values and R2 of each marker
	bilan<-data.frame(marqueurs=colnames(X),pvalue=GWAS.mod1$p, r2=GWAS.mod1$v)

	# Add the mean of each allele
	bilan=cbind(bilan,my_mean)
	
	# Link with the genetic map
	bilan=merge(map,bilan,by.x=2,by.y=1,all=T)
	
	# Compute the LOD score
	bilan$LOD=-log10(bilan$pvalue)
	
	# Add variable name
	bilan$variable=colnames(BLUP)[select]
	bilan=bilan[ , c(2,1,3:ncol(bilan))]
	
	#Nom des colonnes
	colnames(bilan)=c("LG","marqueur","Distance","group_physique","Posi_physique","pvalue","R2","moy.A","moy.B","a","LOD","variable")
	return(bilan)
	
	
	# Close the function
	}

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#





#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
# -------------------------------------------------
# PARTIE 4 : CALCUL DU BILAN_SIMPLE_MARKER EN UTILISANT LA FONCTION
# -------------------------------------------------

to_check=names(BLUP)[sapply(BLUP,is.numeric)==TRUE]
print("les variables a analyser sont : ")
print(to_check)

bilan_simple_marker=data.frame()
print("variables faites : ")
for (i in to_check){
  print(i)
  a=run_my_QTLREL(which(colnames(BLUP) == i))
  bilan_simple_marker=rbind(bilan_simple_marker,a)
}

# écriture du bilan
write.csv(x=bilan_simple_marker,file="bilan_simple_marker",row.names=F)

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#




