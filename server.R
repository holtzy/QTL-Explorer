
	#####################
	#
	#	DEVELOPPEMENT D'UNE APPLI SHINY POUR LA VISUALISATION DES QTLS MOSAIQUES ET FUSA
	#
	####################



library(shiny)
library(plotly)
library(FactoMineR)

my_colors=c(rgb(242,104,34,maxColorValue = 255),rgb(254,188,18,maxColorValue = 255) ,rgb(143,199,74,maxColorValue = 255),rgb(74,132,54,maxColorValue = 255), rgb(111,145,202,maxColorValue = 255) ,rgb(236,33,39,maxColorValue = 255),rgb(165,103,40,maxColorValue = 255))



#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#



shinyServer(function(input, output) {


  
  #### Mes inputs:
  # --- chromo = liste avec tous les chromosomes choisis par l'utilisateur
  # --- to_check = liste avec les variables séléctionnés
  # --- Experiment = FUSA ou Mosaique




# --------------------------------------------------------------------------------
# 	PARTIE 1 : CHARGEMENT DU FICHIERS SELON LE REPERTOIRE INDIQUE
#--------------------------------------------------------------------------------

  # fichier de phénotypage 
  Y<- reactive({
    file_path=paste(input$Experiment,"/phenotypage.csv",sep="")
      Y=read.table(file = file_path , header = TRUE, sep = ";", dec = ".", na.strings = "NA")
      colnames(Y)[1]="geno"
      return (Y)
  })
  
  # fichier comprenant les calculs
  bilan_simple_marker <- reactive ({
    file_path=paste(input$Experiment,"/bilan_simple_marker",sep="")
    bilan_simple_marker <- read.csv(file=file_path)
    bilan_simple_marker=bilan_simple_marker[ order(bilan_simple_marker$LG,bilan_simple_marker$Distance) , ]
    return(bilan_simple_marker)
  })

  # fichier avec tous les génotypes
  geno <- reactive({
    file_path=paste(input$Experiment,"/genotypage.csv",sep="")
    geno <- read.table(file=file_path, sep = ";" , header = F, na.strings = "-")
    geno=as.matrix(geno)
    colnames(geno)=geno[1,]
    geno=as.data.frame(geno[-1 , ])
    return (geno)
  })
  
  # la carte génétique (position en cM et physical)
  map <- reactive({
    file_path=paste(input$Experiment,"/carte",sep="")
    map <- read.table(file=file_path , header=T , dec = ".", na.strings = "-" , check.names=F)
    colnames(map) <- c("LG", "marqueur", "Distance","group_Americain","Posi_physique")
    rownames(map) <- map$marqueur
    map$LG <- as.factor(map$LG)
    map=map[ order(map$LG,map$Distance) , ]
    return(map)
  })
  
  # le fichier d'expression(si il existe)
  expre <- reactive({
    file_path=paste(input$Experiment,"/expression",sep="")
    if (file.exists(file_path)) {
      expre<-read.table(file=file_path,header=T,dec=".")
    }
    
  })


  # Print de la vérification
  output$verif <- renderPrint({
    Y=Y()
    geno=geno()
    map=map()
    verif=as.data.frame(matrix(0,6,2)) ; verif[,1]=c("nbr de geno dans la map" , "nbr geno dans fichier de génot" , "nbr indiv dans fichier de génot" , "nbr d'indiv dans fic phénot" , "nbr geno communs carte / genotypage","nbr indiv communs genot / phenot")
    verif[1,2]=nrow(map) ; verif[2,2]=ncol(geno)-1 ; verif[3,2]=nrow(geno) ; verif[4,2]=nrow(Y) ; verif[5,2]=length(colnames(geno)[colnames(geno)%in%map$marqueur==TRUE])  ; verif[6,2]=length(Y[,1][ Y[,1]%in%geno[,1]==TRUE])
    verif
  })
 
 
 
 

# --------------------------------------------------------------------------------
# 	PARTIE 2 : CHOIX DYNAMIQUE DES VARIABLES + CHROMOSOME SELON LE DOSSIER SELECTIONNE
#--------------------------------------------------------------------------------

  # --- Dynamic UI for the VARIABLE to study
  output$choose_carac <- renderUI({
    # If missing input, return to avoid error later in function
    if(is.null(input$Experiment)){
      return()
    }
      
    # Sélection de toutes les variables numériques
    file_path=paste(input$Experiment,"/phenotypage.csv",sep="")
    Y=read.table(file = file_path , header = TRUE, sep = ";", dec = ".", na.strings = "NA")
    colnames(Y)[1]="geno"
    colnames <- names(Y[sapply(Y,is.numeric)==TRUE])
    
    # Create the checkboxes and select them all by default
    checkboxGroupInput("to_check", "Choose columns", choices=colnames,selected=colnames[1] )
    
  })
  


  # --- Dynamic UI for the CHROMOSOME to study
  output$choose_chromo <- renderUI({
    # If missing input, return to avoid error later in function
    if(is.null(input$Experiment)){
      return()
    }
      
    # Sélection de toutes les variables numériques
    file_path=paste(input$Experiment,"/carte",sep="")
    map <- read.table(file=file_path , header=T , dec = ".", na.strings = "-" , check.names=F)
    colnames(map) <- c("LG", "marqueur", "Distance","group_Americain","Posi_physique")
    chromo_avail=levels(map$LG)
    
    # Create the checkboxes and select them all by default
 	selectInput( "chromo", "Choose a chromosome for deep exploration?",choices =chromo_avail,selected =chromo_avail[1] )
   
  })

  
  
  
  
# --------------------------------------------------------------------------------
# 	PARTIE 3 : REALISATION DES PLOTS DE CHAQUE CHROMOSOME
#--------------------------------------------------------------------------------

  observe ({

    # GRAPHIQUE AVEC TOUS LES CHROMOSOMES
    LG=levels(bilan_simple_marker()$LG);
    output$all_chromo <- renderPlot({
      par(mfrow=c(length(LG),1))
      par(mar=c(1,3	,1,1))
      for(chromo in LG){
        debut="yes"
        num=0
        for (carac in input$to_check){
          num=num+1
          res=bilan_simple_marker()
          #On va garder les lignes de res pour lesquels on a bien un LOD (ls SNP est génotypé pour la pop)
          res=res[!is.na(res$pvalue) , ]          
          #On prend seulement les résultats correspondant à la variable en question et au chromosome en question
          res=res[res$variable== carac & res$LG == chromo, ]
          #Si l'option physique est cochée, je dois prendre les positions en pb :
          if (input$Distance == "physical"){ res$Distance <- as.numeric(res$Posi_physique) ;  res=res[res$LG==res$group_physique , ]  ; res=res[order(res$Distance),] }
          if(debut=="yes"){
            plot(res$LOD ~ res$Distance, ylim=c(0,8), xlim=c(0, max(res$Distance , na.rm=T)), xlab = "", ylab = "LOD", bty="l" , main="" , lwd=3 , type="l" , col=my_colors[num] )
            abline(h = input$LOD_seuil, col = "grey")
            text( 10 , 5 ,labels = paste(chromo," | ",nrow(res),sep="") , col="orange" , cex=2)
            abline(v=res$Distance , col="grey" , lwd=0.8)
            debut="no"
          }else{
            points(res$LOD ~ res$Distance,  lwd=3 , type="l" , col=my_colors[num] ) ; 
          }
        }
      }
    })
    
    
    
    
    
    # PLOTLY POUR LE ZOOM SUR UN CHROMOSOME :
    data=bilan_simple_marker()
    data=data[data$variable%in%as.factor(input$to_check) , ]
    data=data[data$LG==input$chromo , ]
    #On va garder les lignes de res pour lesquels on a bien un LOD (ls SNP est génotypé pour la pop)
    data=data[!is.na(data$pvalue) , ]          
    fun=function(x){return(length(x[x=="A" & !is.na(x) ]) )} ; fun2=function(x){return(length(x[x=="B" & !is.na(x) ]) )} ; fun3=function(x){return(length(x[ is.na(x) ]) )}
    tmp=apply(geno()[,-1] , 2 , fun) ; tmp2=apply(geno()[,-1] , 2 , fun2) ; tmp3=apply(geno()[,-1] , 2 , fun3)
    tmp=data.frame(SNP=names(tmp), nb_A=tmp , nb_B=tmp2 , nb_NA=tmp3 ) 	
    data=merge(data,tmp,by.x=2 , by.y=1 , all.x=T)
    if (input$Distance == "physical"){ data$Distance <- data$Posi_physique  ; data=data[data$LG==data$group_physique , ]  }	
    data=data[order(data$Distance) , ]
    my_text <- paste(data$marqueur , "\n" , "Allele A : ",round(data$moy.A,2) , " | nb :", data$nb_A , "\nAllele B : ", round(data$moy.B,2) , " | nb :", data$nb_B , "\nMissing data : ", data$nb_NA , "\nr2 : ",  round(data$R2,2) , sep="" )    
    print(head(data))
    output$chromo_zoom <- renderPlotly({ 
    	plot_ly( type="scatter", mode="markers" ) %>%  
				add_trace(data=data, y=~LOD , x=~Distance , split=~variable, hoverinfo='text+x', text=my_text) %>%
				layout( 
					xaxis = list( title="" ) , 
					annotations = list(x = 10, y = 9, text = input$chromo, showarrow = F ,font=list(color="orange") ) 
					)   
				}) 
    

    
       
    
    
    # GRAPHIQUE DE L ACP PARTIE 1
    output$PCA1 <- renderPlotly({
      
      # Préparation des données
      marker_to_study=which( colnames(geno())%in%input$selected_marker )
      don=Y()[ , c(1,which(colnames(Y())%in%input$to_check)) ]
      don=merge(don , geno()[, c(1,marker_to_study)] , by.x=1 , by.y=1 , all.x=T)
      colnames(don)[ncol(don)]="my_allele"
      don$my_allele=sub("A" , "Allele Dic2" , don$my_allele) ; 	don$my_allele=sub("B" , "Allele Silur" , don$my_allele)
      don$my_allele[is.na(don$my_allele)]="no_genotype"
      don$text_for_PCA="-"
      for (i in c(1:nrow(don)-3)){
        my_text=paste(Y()$geno[i] , "\n" , sep="")
        for( col in c(1: (ncol(don)-1) )){
          my_text=paste(my_text,colnames(don)[col]," : ",sep="")
          my_text=paste(my_text,don[i,col],"\n",sep="")
        }
        don$text_for_PCA[i]=my_text
      }
      
      # Si j ai plus de 2 variables, je fais une ACP	
      if(ncol(don)>5){
        res.PCA = PCA(don[ , c(2:(ncol(don)-2) ) ] , scale.unit=TRUE, ncp=3, graph=F) 
        AA=data.frame(ind=don[,1] , res.PCA$ind$coord)
        don=merge(don,AA , by.x=1 , by.y=1 , all.x=T)
        plot_ly(data=don, x = ~Dim.1, y = ~Dim.2, text = don$text_for_PCA ,  mode = "markers" , split=~my_allele , hoverinfo='text' , marker=list(size=rep(15,nrow(don))) )
      }
      
      # Si je n ai que 2 variables, alors je fais juste un graph y~x
      else{
        plot_ly(data=don , x=don[,2] , y=don[,3] , type="scatter" , mode="markers" , text = don$text_for_PCA , split=~my_allele , hoverinfo='text' , marker=list(size=rep(15,nrow(don))) )  %>%  layout( yaxis=list(title=colnames(don)[3]) , xaxis = list(title=colnames(don)[2])  )
      }
    })
    
    
    # GRAPHIQUE DE L ACP PARTIE 2
    output$PCA2 <- renderPlot({
      marker_to_study=which( colnames(geno())%in%input$selected_marker )
      
      # ACP et sortie plotly
      don=Y()[ , which(colnames(Y())%in%input$to_check) ]
      res.PCA = PCA(don , scale.unit=TRUE, ncp=2, graph=F) 
      AA=as.data.frame(res.PCA$ind$coord)
      plot.PCA(res.PCA, axes=c(1, 2), choix="var")
    })
    
    
    # TABLEAU DES CORRELATIONS ENTRE VARIABLES
    output$my_cor <- renderPrint({
      marker_to_study=which( colnames(geno())%in%input$selected_marker )
      don=Y()[ , which(colnames(Y())%in%input$to_check) ]
      round(cor(don , use="complete.obs"),2)
    })
    
    
    # BOXPLOT DES VARIABLES SELON L ALLELE
    output$my_boxplot <- renderPlot({
      marker_to_study=which( colnames(geno())%in%input$selected_marker )
      don=Y()[ , c(1,which(colnames(Y())%in%input$to_check)) ]
      don=merge(don , geno()[, c(1,marker_to_study)] , by.x=1 , by.y=1 , all.x=T)
      don=don[ , -1]
      colnames(don)[ncol(don)]="my_allele"
      don$my_allele=sub("A" , "Allele Dic2" , don$my_allele) ; 	don$my_allele=sub("B" , "Allele Silur" , don$my_allele)
      
      par(mfrow=c(1 , ncol(don)-1))
      for(i in c(1 : (ncol(don)-1))){
        #tracé du boxplot pour chaque phénotype
        boxplot(don[,i] ~ don[,ncol(don)] , col=my_colors[c(1,2)] , main=colnames(don)[i] )
        #calcul du pourcentage de variance expliqué 
        l<-lm(don[,i] ~ don[,ncol(don)])
        text_variance=paste("R2: ",round(summary(l)$r.squared,3))
        mtext(text=text_variance,side = 1)
      }
    })
    
    
   
# --------------------------------------------------------------------------------
# 	PARTIE 4 : GRAPHIQUE D'EXPRESSION SELON LE MARQUEUR CHOISI
#--------------------------------------------------------------------------------

    # Je commence par crééer un tableau que j appelle expression : il contiendra : nom du gène / position / resultat de l anova
    expression<- reactive({
      chromo=input$chromo
      expre=expre()
      marker=input$selected_marker
      # on prend seulement les lignes associées au chromosome choisi.
      expre <- expre [expre$chromo_BW==chromo,]
      expre <- t(expre)
      # on récupère la colonne de geno qui est associé au marker.
      allele_marker <- geno()[,c("SNP",input$selected_marker)]
      B<-merge(allele_marker,expre,by.x="SNP",by.y="row.names",all.y=T)
      # on fait l'anova et on récupère la pvalue
      row.names(B) <- B[,1]
      
      calcul_pvalue <- function(contig) {
        l<-aov(as.numeric(as.character(B[contig][[1]]))~(B[,2]))
        pvalue <- summary(l)[[1]][["Pr(>F)"]][[1]]
        return (pvalue)
      }
      calcul_nom_contig <- function(contig) {
        return (as.character(B["contig",contig]))
      }
      calcul_position <- function(contig){
        return (as.numeric(as.character(B["position_BW",contig])))
      }
      options(warn=-1)  
      pvalue <- unlist(lapply(colnames(B)[-c(1,2)],calcul_pvalue))
      nom_contig <- unlist(lapply(colnames(B)[-c(1,2)],calcul_nom_contig))
      posi <- unlist(lapply(colnames(B)[-c(1,2)],calcul_position))
      options(warn=0)
      return(data.frame(nom_contig,posi,pvalue))
    })
    
      
      # Je peux alors faire mon plot
      output$my_expression <- renderPlotly({
        expression <- expression()
        expression$logpval=-log10(expression$pvalue)
        expression=expression[expression$logpval>1.5 , ]
        plot_ly(data=expression, x=~posi, y=~logpval , text = nom_contig ,  mode = "markers" , hoverinfo="text" , marker=list(size=logpval*6 , color=logpval*2 ))
        #plot(expression$posi,-log10(expression$pvalue),xlab="distance physique",ylab="expression en -log10",main="gènes différentiellement exprimés en fonction de sa position physique")
    })
    
	 # Je fais aussi une liste des gènes exprimés    
      output$gene_exprime <- renderPrint({
        expression <- expression()
        out = expression[-log10(expression$pvalue)>input$LOD_seuil & !is.na(expression$pvalue),]  
        out = out[order(out$pvalue),]
        print (out)
      })
  })

})



