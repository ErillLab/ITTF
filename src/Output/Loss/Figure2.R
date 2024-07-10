library("ggplot2")
library("patchwork")
library("gridExtra")
library("grid")
library("DescTools")
library("rjson")
setwd("C:/Users/thele/OneDrive/Documents/Anaconda/envs/IT_TF/src/")
source("AUC.R")

setwd("C:/Users/thele/OneDrive/Documents/Anaconda/envs/IT_TF/src/Output/Loss/")

varvec = as.character(c(0.5,1,1.5))
MIvec = vapply(varvec, function(x) paste(x, "MI", sep = ""), FUN.VALUE = character(length = 1L))
ICvec = vapply(varvec, function(x) paste(x, "IC", sep = ""), FUN.VALUE = character(length = 1L))
togethervec = c(MIvec,ICvec)
title = "Information loss in MI-encoding population and IC-encoding population after Generations of Mutations"
jsonParameters = fromJSON(file = "C:/Users/thele/OneDrive/Documents/Anaconda/envs/IT_TF/src/settings.json" )
MutationRate = 0.01
subtitle = "95% confidence intervals using Standard Error of the Mean are visualized through the shaded area"
xlab = "# Generations of Mutations"
ylab = "Information (Bits)"
filename = "/Users/thele/OneDrive/Documents/Anaconda/envs/IT_TF/src/Output/Figures/Paper/Figure1.png"

lossdf = list()
graphs = list()
maxInformation = 0

#Adds all 4 strategies to lossdf list as a dataframe and updates maxInformation as largest amount of information
for (pop in togethervec){
  #lossdf[[as.character(var)]] = read.csv(paste("./ConstrictOneModeInformation/population_100_100_5_0.5_", var, "IC_MImode.csv", sep = ""))
  lossdf[[pop]] = read.csv(paste("./Paper/Figure_2/population_100_128_2_1bit_",pop, ".csv", sep = ""))
  colnames(lossdf[[pop]]) <- c("Generation", "IC", "MI", "TI", "SEM_IC", "SEM_MI")
  lossdf[[pop]][["Method"]] = rep(pop, nrow(lossdf[[pop]]))
  lossdf[[pop]][["Pane"]] = rep(substring(pop,first = 1, last = nchar(sample)-2), nrow(lossdf[[pop]]))
  #colnames(lossdf[[as.character(var)]]) <- c("Generation", "KLD")
  
  if (lossdf[[pop]]$TI[1] > maxInformation){
    maxInformation = lossdf[[pop]]$TI[1]
  }
  lossdf[[pop]]$ICsemMAX = lossdf[[pop]]$IC + 2*lossdf[[pop]]$SEM_IC
  lossdf[[pop]]$ICsemMIN = lossdf[[pop]]$IC - 2*lossdf[[pop]]$SEM_IC
  
  lossdf[[pop]]$MIsemMAX = lossdf[[pop]]$MI + 2*lossdf[[pop]]$SEM_MI
  lossdf[[pop]]$MIsemMIN = lossdf[[pop]]$MI - 2*lossdf[[pop]]$SEM_MI
}


hybridDataFrame = do.call(rbind, lossdf)

#Iterates through lossdf and adds ggplots for each df to graphs
var1 <- ggplot(hybridDataFrame) + facet_wrap(~Pane) + geom_line(aes(x = Generation, y = IC, colour = "IC")) + geom_line(aes(x = Generation, y = MI, colour = "MI")) + scale_color_manual(values = c("red", "blue")) + 
 labs(colour = "Information Encoding Strategy") + ylim(0,maxInformation)
var1 = var1 + geom_ribbon(aes(x = Generation, ymax = ICsemMAX, ymin = ICsemMIN), alpha = 0.25, fill = "red") 
var1 = var1 + geom_ribbon(aes(x = Generation, ymax = MIsemMAX, ymin = MIsemMIN), alpha = 0.25, fill = "blue")
var1 = var1 + xlab("# of Generations of Mutation") + ylab("Information (bits)")
finalGraph = var1 + plot_annotation(title = title, subtitle = subtitle) + theme_bw()

#Finds asymptote
asymptotevec <- list()
asymptotevec[["1IC"]] = AUC_InformationCurve(hybridDataFrame$Generation, hybridDataFrame$IC)
asymptotevec[["1MI"]] = AUC_InformationCurve(hybridDataFrame$Generation, hybridDataFrame$MI)

#Computes asymptote AUC as #num generations * min(Value)
oneICAUC = asymptotevec[["1IC"]]$AUCval - (max(lossdf[["1IC"]]$Generation) * asymptotevec[["1IC"]]$TIasym)
oneMIAUC = asymptotevec[["1MI"]]$AUCval - (max(lossdf[["1MI"]]$Generation) * asymptotevec[["1MI"]]$TIasym)

finalGraph = finalGraph + annotate("label", x = 250, y = 0.95, label = paste("IC curve AUC:", round(oneICAUC, digits = 2)))
finalGraph + annotate("label", x = 250, y = 0.86, label = paste("MI curve AUC:", round(oneMIAUC, digits = 2)))

ggsave(filename, scale = 1.25)