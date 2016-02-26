


	--- QTL Explorer ---

		A shiny application that permits to do a simple marker QTL analysis
		Yan Holtz et Alban Besnard
						yan1166@hotmail.com
						albanbesnard@hotmail.fr





1 - WHAT IS QTL EXPLORER
————————————————————————

It is an application that permits to visualize QTL with modern tool : Shiny and D3.Js (called with the plot.ly package).
To see how it looks like, have a look here :
www.

You need 3 files for QTL analysis, that are the 3 input of the application :
	- a genetic map
	- a matrix of genotyping
	- a phenotyping matrix





		
2 - INPUT FILES
———————————————

	— 1/ Genetic Map

It gives the genetic and/or the physical position of your markers (SNP / Dart / SSR…)
3 or five columns separated with tabulations. Column fields by order :

 1 - group: the chromosome or linkage group
 2 - marker: markers names
 3 - posi: position of the marker. Decimal = « . »
 4 - group_Americain: the chromosome or linkage group given by another source (optional)
 5 - posi_physique: position for this second reference (optional)
 
#Les deux dernières colonnes sont optionnelles (je crois) 
#mais une fois dans l'appli,n'essayez pas de changer en physical si vous ne les possédez pas.
 

	— 2/ Phenotyping

A file in the .csv format that gives the phenotyping features of your individuals.
One line per individual.
The first column gives the individual names. It names must be « geno »
Then you can add as many column that you need, each giving a phenotype : ex: size / precocity / color…



	— 3/ Genotyping matrix
A file in .csv format. Each line represents an individual. Each column represents a marker, and is called with the marker name.
Then alleles must be coded A, B an - for parent1, parent2 and missing data.
vous avez d'individu.



	— 4/ Example

Genetic Map :

group	marker	posi
1A	mark1	0
1A	mark2	15
1A	mark3	60.3
2A	mark4	0


Genotyping Matrix :

SNP;mark1;mark2;mark3;mark4;mark5;mark6;mark7;mark8;mark9;mark10
Cindy;A;A;A;A;A;A;A;B;B;B
Charles;A;B;A;B;B;B;B;B;B;B


Phenotyping :

geno;Taille;Poids;Age
Cindy;1.80;80;40
Charles;1.75;68;35





3 - HOW TO USE THE APP
——————————————————————

0- Perform a t-test on every marker for every phenotype with the script from_geno_pheno_map_to_bilan_simple_marker.R.
You can call it like that :
Script from_geno_pheno_map_to_bilan_simple_marker.R genotyping_file phenotyping_file  genetic_map_file

1- Create a folder for your work. In this folder create the folder : One called SHINY_APP, where you put the server.R and the ui.R files. One where
you place your data. 
Your file must be named exactly:bilan_simple_marker carte genotypage.csv phenotypage.csv

2- Open R and install the packages (available on CRAN):  plotly shiny and FactoMineR
	install.packages("plotly")
	install.packages("shiny")
	install.packages("FactoMineR")

3- Charge shiny
	library(shiny)

4- Change your working directory to the folder you created
	setwd("Home/my_appli")

5- Run the App !
	runApp("SHINY_APP_FOR_QTL_ANALYSIS")





3 - BONUS FOR EXPRESSION
——————————————————————

#Pour voir les diférentiels d'expressions selon un marqueur donné il vous faut rajouter un fichier de compte (de préférence normalisé)
#Le seul disponible pour le moment est celui sur Dic2 Silur sur un grain immature à 300°jour.

