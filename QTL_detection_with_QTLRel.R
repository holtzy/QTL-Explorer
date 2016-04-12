
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

#setwd("/NAS/davem_data/EPO/SAUVEGARDE_DATA_YAN/DIC2_SILUR/QTL/PUBLI")
#fic_geno="genotypage.csv"
#fic_pheno="phenotypage.csv"
#fic_map="carte"


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
# PARTIE 3 : CHANGEMENT FORMAT POUR QTLREL
# -------------------------------------------------

# --- Geno
rownames(geno)=geno[,1] ; geno=geno[,-1] ; geno=as.matrix(geno)
geno[which(geno=="A")]<-"AA"
geno[which(geno=="B")]<-"BB"

# On va remplacer les données manquantes par un allele au hasard.
my_fun=function(x){length(x[x=="AA"])/length(!is.na(x)) }
prop=apply(geno , 2 , my_fun)
for(i in c(1:ncol(geno))){
	aa=geno[,i][is.na(geno[,i])]
	bb=rbinom(length(aa),1,prob=prop[i])
	geno[,i][is.na(geno[,i])]=c("BB","AA")[bb+1]
	}

# --- Phéno
rownames(Y)=Y[,1] ; Y=Y[,-1]

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#







#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
# -------------------------------------------------
# PARTIE 4 : FONCTION QUI FAIT TOURNER QTLREL SUR UNE VARIABLE PHENO
# -------------------------------------------------


run_my_QTLREL=function(select){

	# ---- SELECTION PHENO + FABRICATION FICHIER CORESPONDANT
	
	# Il faut une matrice avec les memes indiv dans les matrices de phéno et de géno
	y=Y[,select]
	y=as.data.frame(y)
	rownames(y)=rownames(Y)
	y=na.omit(y)
	my_geno=geno[which(rownames(geno)%in%rownames(y)) ,  which(colnames(geno)%in%map$marqueur)]
	y=as.matrix(y)
	y=as.data.frame(y[rownames(y)%in%rownames(geno) , ])
	
	#Idem pour la carte!
	my_map=map[which(map$marqueur%in%colnames(my_geno)) , ]
	
	#Matrice identité
	I<-diag(nrow(y))
	
	
	# --- PASSAGE DE QTLREL
	
	# Calcul matrice apparentement IBS
	K<-genMatrix(my_geno)
	# On visualise la matrice K comme ceci (2=apparentement complet) :
	#K$AA[1:10 , 1:10] ; hist(K$AA[1:10 , 1:10] , breaks=20)
	
	# le premier modele sert a calculer la variance
	y=as.matrix(y)
	mod1<-estVC(y=y,v=list(AA=K$AA,DD=NULL,HH=NULL,AD=NULL,MH=NULL,EE=I))
	mod1
	
	# Genet assoc
	GWAS.mod1<-scanOne(y=y,gdat=my_geno,vc=mod1,test="Chisq")
	
	
	# Récupération des moyennes par alleles:
	my_function=function(x, y) {
		moy <- aggregate(y, by = list(marqueur = x), mean, na.rm = TRUE)
		resu=data.frame(moy.A=moy$y[moy$marqueur == "AA"],moy.B= moy$y[moy$marqueur =="BB"] , a=abs(moy[1, 2] - moy[2, 2])/2)
		return(resu)
		}
	my_mean = apply(my_geno, MARGIN = 2, FUN = function(x,y) my_function(x, y), y=y)
	my_mean = matrix(as.vector(unlist(my_mean)) , ncol=3 , byrow=T)	
	
	
	
	# --- FORMATION TABLEAU BILAN
	
	#Récupération des pvalues + r2 expliqué.
	bilan<-data.frame(marqueurs=colnames(my_geno),pvalue=GWAS.mod1$p, r2=GWAS.mod1$v)
	#Ajout des moyennes
	bilan=cbind(bilan,my_mean)
	#Ajout des données de carto
	bilan=merge(map,bilan,by.x=2,by.y=1,all=T)
	#Ajout LOD
	bilan$LOD=-log10(bilan$pvalue)
	#Ajout Variable
	bilan$variable=colnames(Y)[select]
	bilan=bilan[ , c(2,1,3:ncol(bilan))]
	
	#Nom des colonnes
	colnames(bilan)=c("LG","marqueur","Distance","group_physique","Posi_physique","pvalue","R2","moy.A","moy.B","a","LOD","variable")
	return(bilan)
	}

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#





#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
# -------------------------------------------------
# PARTIE 4 : CALCUL DU BILAN_SIMPLE_MARKER EN UTILISANT LA FONCTION
# -------------------------------------------------
to_check=names(Y)[sapply(Y,is.numeric)==TRUE]
print("les variables a analyser sont : ")
print(to_check)

bilan_simple_marker=data.frame()
print("variables faites : ")
for (i in to_check){
  print(i)
  a=run_my_QTLREL(which(colnames(Y) == i))
  bilan_simple_marker=rbind(bilan_simple_marker,a)
}

# écriture du bilan
write.csv(x=bilan_simple_marker,file="bilan_simple_marker",row.names=F)

