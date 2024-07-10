library("ggplot2")
library("patchwork")
library("gridExtra")
library("grid")
library("DescTools")
library("rjson")
setwd("C:/Users/thele/OneDrive/Documents/Anaconda/envs/IT_TF/src/")
source("AUC.R")

setwd("C:/Users/thele/OneDrive/Documents/Anaconda/envs/IT_TF/src/Output/Loss/")

varvec = as.character(c(0,0.25,0.5,0.75,1))
title = "Loss of information in various encoding strategies after Generations of Mutations"
jsonParameters = fromJSON(file = "C:/Users/thele/OneDrive/Documents/Anaconda/envs/IT_TF/src/settings.json" )
MutationRate = 0.01
subtitle = "95% confidence intervals using Standard Error of the Mean are visualized through the shaded area"
xlab = "# Generations of Mutations"
ylab = "Information (Bits)"
filename = "/Users/thele/OneDrive/Documents/Anaconda/envs/IT_TF/src/Output/Figures/Paper/SuppFigure_1.png"
fileDir = "./Paper/SuppFigure_1/"
lossdf = list()
graphs = list()
maxInformation = 0

#Adds all 4 strategies to lossdf list as a dataframe and updates maxInformation as largest amount of information
for (var in varvec){
  popType = as.character(var)
  lossdf[[popType]] = read.csv(paste(fileDir, "population_100_128_2_1bit_",var, "IC.csv", sep = ""))
  
  if (lossdf[[popType]]$TI[1] > maxInformation){
    maxInformation = lossdf[[popType]]$TI[1]
  }
  lossdf[[popType]]$Method = rep(paste(popType,"IC",sep = ""), nrow(lossdf[[popType]]))
  lossdf[[popType]]$TIsemMAX = lossdf[[popType]]$TI + 2*lossdf[[popType]]$SEM_TI
  lossdf[[popType]]$TIsemMIN = lossdf[[popType]]$TI - 2*lossdf[[popType]]$SEM_TI
}

hybridDataFrame = do.call(rbind, lossdf) 
hybridDataFrame$Method = factor(hybridDataFrame$Method, levels = c("1IC", "0.75IC", "0.5IC", "0.25IC", "0IC"))

#Iterates through lossdf and adds ggplots for each df to graphs
var1 <- ggplot(hybridDataFrame) + geom_line(aes(x = X..Generation.after.mutation, y = TI, color = Method)) + scale_color_manual(values = c("red", "darkred","purple","blue","darkblue")) + labs(colour = "Information Encoding Strategy") + ylim(0,maxInformation) 
var1 = var1 + geom_ribbon(aes(x = X..Generation.after.mutation, ymax = TIsemMAX, ymin = TIsemMIN, fill = Method), alpha = 0.25, show.legend = FALSE) + scale_fill_manual(values = c("red", "darkred","purple","blue","darkblue"))
var1 = var1 + xlab("# of Generations of Mutation") + ylab("Information (bits)")
finalGraph = var1 + plot_annotation(title = title, subtitle = subtitle) + theme_bw()

#Finds asymptote
asymptotevec <- list()
ypos = 0.95
for (var in varvec){
  popType = as.character(var)
  asymptotevec[[popType]] = AUC_InformationCurve(lossdf[[popType]]$X..Generation.after.mutation, lossdf[[popType]]$TI)
  
  #Computes asymptote AUC as #num generations * min(Value)
  aucvalue = asymptotevec[[popType]]$AUCval - (max(lossdf[[popType]]$X..Generation.after.mutation) * asymptotevec[[popType]]$TIasym)
  finalGraph = finalGraph + annotate("label", x = 250, y = ypos, label = paste(popType, "IC curve AUC:", round(aucvalue, digits = 2), sep = ""))
  ypos = ypos - 0.08
}


ggsave(filename, scale = 1.25)