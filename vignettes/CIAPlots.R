CoinertiaPlot <- function(coin,
                  Quant = 0.9, Prop.Var = 0.9,
                  Env_Var=NULL, Env_Var2 = NULL, color=NULL, shape=NULL,
                  PtColor= "magenta",PtShape=4, PtSize=2,
                  linetype=2, LblSize=2, LabelsOpt = NULL,
                  ArrLen=0.15, ArrAngle=20){

axes <- c(1,2)
#Plot scores from each data set
XLim <- c(min(coin$mX$NorS1, coin$mY$NorS1), 
          max(coin$mX$NorS1, coin$mY$NorS1))

YLim <- c(min(coin$mX$NorS2, coin$mY$NorS2), 
          max(coin$mX$NorS2, coin$mY$NorS2))
#Get summary of distances between the two co-inertia sets
  Inert <- round(cumsum(coin$eig/sum(coin$eig)),6)
  percvar <- round((100 * coin$eig/sum(coin$eig))[axes],2)
#Find distance based on the first k axes that explain more that Prop.Var of Inertia  
  k <- min(which(Inert > Prop.Var))
  k <- min(k, dim(coin$mY)[2])#min of var explained and number of axes kept 
  Dist <- coin$mX[,1:k] - coin$mY[,1:k]
  Length <- sqrt(Dist$NorS1^2 + Dist$NorS2^2)
     if(all(substr(rownames(coin$mX), 1,1) == "X")){Char.rm <- TRUE}
     else{Char.rm <- FALSE}
  names(Length) <- rownames(coin$mX)
     if(Char.rm == TRUE){names(Length) <- substring(names(Length), first = 2)}
  MidPoint <- (coin$mX + coin$mY)/2
  Quantiles <- unique(sort(c(seq(0,1,0.25), seq(0,1, 0.1))))
  Summary <- as.data.frame(matrix(quantile(Length, probs = Quantiles), nrow = 1))
  names(Summary) <- paste("Qu_", Quantiles, sep = "")
  length(Length[Length > Summary$Qu_0.9])
#label only the points that are larger that 0.9th quantile
  QPlot <- quantile(Length, probs = Quant)
  Ind  <- which(Length > QPlot)
  Labels <- rownames(coin$mX)[Ind]
      if(Char.rm == TRUE) {Labels <- substring(Labels, first = 2)}
      if(!is.null(LabelsOpt)){Labels <- LabelsOpt[Ind]}
  df3 = data.frame(cbind(MidPoint$NorS1[Ind],  MidPoint$NorS2[Ind]), Labels)
  names(df3)[1:2] <- names(MidPoint)
#From here we change to a ggplot
  x = colnames(coin$mY)[1]
  y = colnames(coin$mY)[2]
#Create a table identifier
  Table <- c(rep(1, nrow(coin$mX)), rep(2, nrow(coin$mX)))
  df <- rbind(coin$mX, coin$mY)
  df <- cbind(df, Table)

  
  x = colnames(df)[1]
  y = colnames(df)[2]
  
  ord_map = aes_string(x = x, y = y)
  
  p <- ggplot(df, ord_map) + geom_point(na.rm = TRUE, 
                                        colour = PtColor, shape = PtShape,size = PtSize)


#Check for additional plotting variables to color arrows 
    #Save df names
  df2 <- cbind(coin$mX, coin$mY)
  names(df2) <- c("xbeg", "ybeg", "xend", "yend")
  df_Names <- names(df2)
  
  if(!is.null(Env_Var)){
    
    if(!is.null(color)){
      Var_Col <- which(names(Env_Var) == color)
      df2 <- data.frame(df2, Env_Var[,Var_Col])
      names(df2) <- c(df_Names,color)}
    
    if (!is.null(shape)) {
      Var_Shape <- which(names(Env_Var) == shape)
      df2 <- data.frame(df2, Env_Var[,Var_Shape])
      names(df2)[ncol(df2) ]<- shape
    }
    
  }#end if not null Env_Var
  
  ord_map_Arr = aes_string(x = "xbeg", y = "ybeg", xend = "xend", yend="yend",
                       color = color, shape = shape, 
                       na.rm = TRUE)
 
  p = p+ geom_segment(ord_map_Arr, data = df2, 
          arrow = arrow(angle = ArrAngle, length = unit(ArrLen, "inches")),linetype=linetype)
  if(is.null(color)){p = p+ geom_segment(ord_map_Arr, data = df2, color = PtColor,
                                         arrow = arrow(angle = ArrAngle, length = unit(ArrLen, "inches")),linetype=linetype)}
#Add percent var explained by axes on plot
  #strivar = as(c(p$label$x, p$label$y), "character")
  strivar = c(paste("Axis", axes[1], sep = ""), paste("Axis", axes[2], sep = ""))
  strivar = paste0(strivar, "   [", percvar, "%]")
  p = p + xlab(strivar[1]) + ylab(strivar[2])

#Finally add labels with
  p= p+annotate("text", x=df3$NorS1, y=df3$NorS2, label= df3$Labels, size = LblSize)

#Produce dissimilarity graph for samples data based on collection days
Dissimilarity <- data.frame(names(Length), Length)
names(Dissimilarity) <-c("Name", "Dissimilarity")
       if(!is.null(LabelsOpt)){Dissimilarity <- cbind(Dissimilarity,LabelsOpt) }
DistQuant <- rep(NA, length(Length))
              for (i in 1:(length(Summary)-1)){
                  Ind <- Length >=Summary[,i] & Length < Summary[,i+1]
                  DistQuant[Ind] <- Quantiles[i] 
                          if(i == (length(Summary)-1)) {DistQuant[Length >= Summary[,i+1]] <- Quantiles[i+1]}
              }#end for
Dissimilarity <- cbind(Dissimilarity, DistQuant)
names(Dissimilarity)[dim(Dissimilarity)[2]] <- "Quantile"

     if(!is.null(Env_Var$CollectionDays)){
         TimeDiff <- Env_Var$CollectionDays - Env_Var2$CollectionDays
         BV_Status <- Env_Var$bv
         #DistQuant <- rep(NA, length(Length))
              #for (i in 1:(length(Summary)-1)){
                  #Ind <- Length >=Summary[,i] & Length < Summary[,i+1]
                  #DistQuant[Ind] <- Quantiles[i] 
                          #if(i == (length(Summary)-1)) {DistQuant[Length >= Summary[,i+1]] <- Quantiles[i+1]}
              #}#end for

        #Dissimilarity <- data.frame(names(Length), Length, TimeDiff,BV_Status, DistQuant)
         Dissimilarity <- cbind(Dissimilarity, TimeDiff,BV_Status)
         last <- dim(Dissimilarity)[2]
         names(Dissimilarity)[(last-1):last] <- c( "TimeBtwVisits","BV")
         names(Dissimilarity)[1] <- "repeat_code"
     }#end if(!is.null(Env_Var$CollectionDays)
Dissimilarity <- Dissimilarity[order(Dissimilarity$Dissimilarity), ]

return(list(Summary = Summary,  Dissimilarity = Dissimilarity, p=p))
}#End of Finction CoinertiaPlot

##############################################################################
PlotDissimilarity <- function(Dissimilarity,Title="", label = "point",
                          plot_colors =c("red",  "blue", "green")){
  #browser()
A <- abs(Dissimilarity$TimeBtwVisits[Dissimilarity$Quantile < 0.5])
Ind <- Dissimilarity$Quantile >= 0.5 & Dissimilarity$Quantile <0.9 
B <- abs(Dissimilarity$TimeBtwVisits[Ind])
C <- abs(Dissimilarity$TimeBtwVisits[Dissimilarity$Quantile >= 0.9])

summary(A)
summary(B)
summary(C)

DissimilTime <- c(A,B,C)

DissimilQuantile <- c(rep("Less_Than_Median", length(A)),
          rep("Greater_Than_Median", length(B)),
          rep("Greater_Than_Q90", length(C)))


A <- as.character(Dissimilarity$BV[Dissimilarity$Quantile < 0.5])
B <- as.character(Dissimilarity$BV[Ind])
C <- as.character(Dissimilarity$BV[Dissimilarity$Quantile >= 0.9])

BV_Status <- c(A,B,C)

A <- as.character(Dissimilarity$repeat_code[Dissimilarity$Quantile < 0.5])
B <- as.character(Dissimilarity$repeat_code[Ind])
C <- as.character(Dissimilarity$repeat_code[Dissimilarity$Quantile >= 0.9])

repeat_code <- c(A,B,C)

df <- data.frame(repeat_code, DissimilTime, DissimilQuantile, BV_Status)
df$BV_Status <- as.factor(df$BV_Status)
df$DissimilQuantile <- factor(df$DissimilQuantile)
#diamonds$cut <- factor(diamonds$cut, levels = rev(levels(diamonds$cut)))
ds <- ddply(df, .(DissimilQuantile), summarise, mean = mean(DissimilTime), sd = sd(DissimilTime))
myColors <- brewer.pal(length(levels(df$BV_Status)),"Set1")
names(myColors) <- levels(df$BV_Status)
colScale <- scale_colour_manual(name = "BV_History",values = myColors)

Groups <- as.character(df$DissimilQuantile)
Levels <- unique(df$DissimilQuantile)
Lab1 <- which(df$DissimilQuantile == Levels[1])
Lab2 <- which(DissimilQuantile == Levels[2])
Lab3 <- which(df$DissimilQuantile == Levels[3])
Groups[Lab1] <- "G1"
Groups[Lab2] <- "G2"
Groups[Lab3] <- "G3"

df <- data.frame(df, Groups)
df$Groups <- as.factor(df$Groups)

Groups <- c("G2", "G3", "G1")
ds <- data.frame(ds, Groups)
df$BV_Status <- factor(df$BV_Status, levels = c("Yes", "No", "NS"))
if(label == "text"){
Plot <- ggplot(df,aes(x =Groups,  y=DissimilTime, color = BV_Status)) + geom_text(aes(x = DissimilQuantile, y = DissimilTime, label=repeat_code, color = BV_Status), size = 2, data=df, parse = T)}

else{Plot <- ggplot(df,aes(x =Groups,  y=DissimilTime, color = BV_Status)) + geom_point(stat="identity") + geom_point(data = ds, aes(y = mean), colour = 'black', size = 3)}

Levels <- unique(DissimilQuantile)
Lab1 <- which(DissimilQuantile == Levels[1])
Lab2 <- which(DissimilQuantile == Levels[2])
Lab3 <- which(DissimilQuantile == Levels[3])
Labels <- rep(0, length(DissimilQuantile))
Labels[Lab1] <- "Q<50"
Labels[Lab2] <- "50 < Q < 90"
Labels[Lab3] <- "Q>90"
#browser()
#browser()
Cols <- c("Yes" = plot_colors[1], "No" = plot_colors[2],"NS" = plot_colors[3])
Brks <- c("Yes", "No", "NS")
Labs <- c("Yes", "No", "NS")
    
Plot <- Plot + ggtitle(Title) + scale_x_discrete(name = "Quantile") + scale_y_continuous(name = "Time Between Visits") +scale_colour_manual(name = "BV History", values = Cols, breaks = Brks, labels = Labs)
#ggsave(file = paste(path, name, ".pdf", sep = ""), Plot)

return(Plot)
}


################################################################
PlotCoinVars <- function(coin, tab1 = "Table1", tab2 = "Table2", 
                         Labels1 = NULL, Labels2 = NULL, label = TRUE,
                         hjust = 0, vjust = -1.5, PtSize = 2, LblSize = 2){

  if(is.null(Labels1)){Labels1 <- 1:nrow(coin$co)} #tab1 labels
  if(is.null(Labels2)){Labels2 <- 1:nrow(coin$li)} #tab2 labels
  
  #extract scores for each table
  x = colnames(coin$co)[1]
  y = colnames(coin$co)[2]
  
  #first table data
  df1 <- data.frame(coin$co$Comp1, coin$co$Comp2, Labels1, rep(tab1, nrow(coin$co)))
  rownames(df1) <- rownames(coin$co)
  names(df1) <- c(x,y,  "Labels", "Table")
  
  #second table data
  df2 <- data.frame(coin$li$Axis1, coin$li$Axis2, Labels2, rep(tab2, nrow(coin$li)))
  rownames(df2) <- rownames(coin$li)
  names(df2) <- c(x,y, "Labels", "Table")
  
  #conbine two tables for plotting
  df <- rbind(df1, df2)
  ord_map = aes_string(x = x, y = y, color = "Table", shape = "Table")
  
  CW_X <- ggplot(df, ord_map) + geom_point( size = PtSize) + 
          scale_color_manual(values = c("red", "blue")) + theme_bw()+
          xlab("") + ylab("") +theme(legend.title=element_blank())
  
  if(label == TRUE) { 
      #Last thing: Fix order of labels
      lbl_map = aes_string(x = x, y = y, label = "Labels")  
      CW_X <- CW_X + geom_text(data = df, mapping = lbl_map, size = LblSize, vjust = vjust, hjust = hjust)
  }
  return(CW_X)
}
#####################################################################



###################################################################
#this plotting function corresponds to the same plot presented in co-inertia papers
#but it contains improved graphics
###################################################################
PlotCW <- function(coin, name, path, color = "red",
            Title1 = "Canonical Weights for the First Visit",
            Title2 = "Canonical Weights for the Second Visit",
            Labels1 = NULL, Labels2 = NULL, scale = FALSE,
            PtShape=2, PtSize=2,
            linetype=1, linesize=0.4,
            LblSize=2, 
            ArrLen=0.15, ArrAngle=20,
            TitleSize = 8){
    
  library(grid)
  library(gridExtra)
  
  if(is.null(Labels1)){Labels1 <- 1:nrow(coin$co)}
  if(is.null(Labels2)){Labels2 <- 1:nrow(coin$li)}
  
  x = colnames(coin$co)[1]
  y = colnames(coin$co)[2]

  ord_map = aes_string(x = x, y = y)
  df <- data.frame(coin$co$Comp1, coin$co$Comp2, 
          rep(0, length(coin$co$Comp1)), rep(0, length(coin$co$Comp2)),
          Labels1)
  rownames(df) <- rownames(coin$co)
  names(df) <- c(x,y, "x_0", "y_0", "Labels1")
  arrows_map =  aes_string(x = "x_0", y = "y_0",
                xend = x, yend = y)
    #browser()
CW_X <- ggplot(coin$co, ord_map) + geom_point(colour = color, 
                                              shape = PtShape,
                                              size = PtSize)
CW_X <- CW_X + scale_x_continuous(name="") + scale_y_continuous(name="")
CW_X <- CW_X + geom_segment(data = df, mapping = arrows_map ,
            size = linesize, 
            colour = color,
            arrow = arrow(angle = ArrAngle, 
                          length = unit(ArrLen, "inches")),
                          linetype=linetype)

#Last thing: Fix order of labels
lbl_map = aes_string(x = x, y = y, label = "Labels1")  
CW_X <- CW_X + geom_text(data = df, mapping = lbl_map, size = LblSize)
CW_X <- CW_X + ggtitle(Title1) +theme(plot.title = element_text(size = TitleSize, colour = "black"),
                 panel.background = element_rect(fill = "white"),
                 panel.grid.major = element_line(colour = "grey90"))


#Canonical Weight for the second data set
x = colnames(coin$li)[1]
y = colnames(coin$li)[2]
ord_map = aes_string(x = x, y = y)
df <- data.frame(coin$li$Axis1, coin$li$Axis2, 
                 rep(0, length(coin$li$Axis1)), rep(0, length(coin$li$Axis2)),
                 Labels2)
rownames(df) <- rownames(coin$li)
names(df) <- c(x,y, "x_0", "y_0", "Labels2")
arrows_map =  aes_string(x = "x_0", y = "y_0",
                         xend = x, yend = y)

CW_Y <- ggplot(coin$li, ord_map) + geom_point(colour = color,
                                              shape = PtShape,
                                              size = PtSize)
CW_Y <- CW_Y + scale_x_continuous(name="") + scale_y_continuous(name="")
CW_Y <- CW_Y + geom_segment(data = df, mapping = arrows_map, 
                            size = linesize, linetype = linetype, colour = color,
                            arrow = arrow(angle = ArrAngle, 
                                          length = unit(ArrLen, "inches")))
lbl_map = aes_string(x = x, y = y, label = "Labels2")  
CW_Y <- CW_Y + geom_text(data = df, mapping = lbl_map, size = LblSize)
CW_Y <- CW_Y + ggtitle(Title2) +theme(plot.title = element_text(size = TitleSize, colour = "black"),
                                      panel.background = element_rect(fill = "white"),
                                      panel.grid.major = element_line(colour = "grey90"))
  
#browser()

if(scale == TRUE){
  
  xmin <- min(coin$co$Comp1, coin$li$Axis1)
  xmax <- max(coin$co$Comp1, coin$li$Axis1)
  
  ymin<- min(coin$co$Comp2, coin$li$Axis2)
  ymax<- max(coin$co$Comp2, coin$li$Axis2)
  
  CW_X <- CW_X + scale_y_continuous(name="", limits=c(ymin,ymax)) + 
         scale_x_continuous(name="", limits=c(xmin,xmax))
  
  CW_Y <- CW_Y + scale_y_continuous(name="", limits=c(ymin,ymax)) + 
    scale_x_continuous(name="", limits=c(xmin,xmax))
}

pushViewport(viewport(layout = grid.layout(1, 2)))
print(CW_X, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(CW_Y, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))

name2 <- paste(name, "_CW", sep = "")
pdf(paste(path, name2, ".pdf", sep = ""))
grid.arrange(CW_X, CW_Y, ncol=2,nrow=1)
dev.off()
 return(list(CW_1 = CW_X, CW_2 = CW_Y))
}
#####################################################################
