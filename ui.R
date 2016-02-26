#####################
#
#	DEVELOPPEMENT D'UNE APPLI SHINY POUR LA VISUALISATION DES QTLS MOSAIQUES ET FUSA
#
####################





# Library
library(shiny)
library(plotly)

# Catch the existing chromosomes in the genetic map:



# Define UI for application that draws every chromosomes with diversity indexes
shinyUI(fluidPage(
  
  # Application title
  titlePanel(paste("\t\t\t\t\t\t","QTL Analysis for Mosaique and Fusariose","\n")),
  
    # ------------------------------    SIDEBAR PANNEL FOR OPTIONS -----------------------------------
  sidebarLayout(
    sidebarPanel(
      

      # INTRO 
      helpText(" 3 pages are provided : Page 1 for a first observation on every chromosomes. Page2 for a deep observation of ONE chromosome. Page 3 to study one marker in particular") ,
      br(""),

      # CHOIX DU FICHIER DE PHENOTYPAGE
      textInput("Experiment",label = "Name of the working directory: " , value="DATA_DS_FUSA"),
      submitButton("Submit"),
      selectInput("Distance",label = "Physical or genetic positions for markers ?", choices = c("cM","physical")),

      
      # CHOIX DU chromosome d'étude
	  uiOutput("choose_chromo"),
      
      # CHOIX DU SEUIL
      sliderInput("LOD_seuil", label = h3("LOD threshold printed on plots ?"), min = 2, max = 7, value = 4,step = 0.1),

	  # CHOIX D UN MARQUEUR
	  textInput("selected_marker", "marker to study ?" , value="Traes_7BS_62A8A6F7E@944"),

      # CHOIX DES CARACTERES
      uiOutput("choose_carac")
     
    ),
    
    
    
    # ------------------------------    PANEAU PRINCIPAL avec 4 onglets -----------------------------------
    mainPanel(
    
    	tabsetPanel(
 
			# --- Onglet 1 avec tous les chromosomes
 			tabPanel("All chromo",
              	# Graphique avec tous les chromosomes:
              	br(""),
				plotOutput("all_chromo" ,  height = "1800px"),
  				# ajout du tableau de vérification
      	        br(""),
				verbatimTextOutput("verif")
      	        ),
 
 
			# --- Onglet 2 : zoom sur un chromosome
 			tabPanel("zoom on chromo",
           	   	plotlyOutput("chromo_zoom" ,  height = "600px")
          	  	), 	
 				
			# --- Onglet 3 : Analyse d un marqueur
 			tabPanel("Analyse d'un marqueur",
           	   	plotlyOutput("PCA1" ,  height = "600px"),
           	   	plotOutput("PCA2" ,  height = "600px"),
				verbatimTextOutput("my_cor"),
				plotOutput("my_boxplot",  height = "600px"),
				plotOutput("my_residuals" , height = "600px")
          	  	),

			# --- Onglet 4 : Analyse de l'expression
		    tabPanel("Analyse de l'expression",
               plotlyOutput("my_expression"),
	           br(""),
               verbatimTextOutput("gene_exprime")
               )
             	   
             	   
             	   )
             )
    
  )
))
