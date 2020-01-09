library(tidyverse)
library(ggplot2)
library(plotly)
library(plyr)
library(data.table)
library(readr)
library(dplyr)
library(ggfortify)
library(factoextra)
library(mdatools)



#Grab the file names of the 20 spectra for 785nm measured Calcium Nitrate
myfiles <- list.files(path = "C:/Users/jays/Desktop/Converted_Spectra/785/Ca_Nitrate", 
                      pattern = "*.txt", full.names = TRUE)
myfiles

#Instantiate an empty data frame to jam the raw data into
ca_nitrate_785 <- data.frame()

#Loop through the files to build one dataframe
for (lambda in c(1:20)) {
  location <- myfiles[lambda]
  temp <- fread(file = location, 
                sep = ";", fill = TRUE, header = TRUE) %>% select(`Dark Subtracted`)
  
  temp <- transpose(temp)
  
  ca_nitrate_785[lambda, 1:582] <- temp
}

#Instantiate an empty data frame to jam the raw data into
mg_nitrate_785 <- data.frame()

#Repeat above steps for the 20 spectra of 785nm measured Magnesium Nitrate
myfiles <- list.files(path = "C:/Users/jays/Desktop/Converted_Spectra/785/Mg_Nitrate", 
                      pattern = "*.txt", full.names = TRUE)

#Loop through the files
for (lambda in c(1:20)) {
  location <- myfiles[lambda]
  temp <- fread(file = location, 
                sep = ";", fill = TRUE, header = TRUE) %>% select(`Dark Subtracted`)
  
  temp <- transpose(temp)
  
  mg_nitrate_785[lambda, 1:582] <- temp
}

ca_calibration <- ca_nitrate_785[5:20, ]
mg_calibration <- mg_nitrate_785[5:20, ]

#Slap on an identification column
ca_nitrate_785 <- ca_nitrate_785 %>% mutate(compound = "ca")
mg_nitrate_785 <- mg_nitrate_785 %>% mutate(compound = "mg")

#Create a master dataframe of all 40 spectra
all_785_spectra <- rbind(ca_nitrate_785, mg_nitrate_785)

#Remove outlier
all_785_spectra <- all_785_spectra[-c(24), ]

ca_model <- simca(ca_calibration, 'calcium nitrate', ncomp = 3, cv = 1, center = TRUE, scale = TRUE)

plot(ca_model)
summary(ca_model)

plotme_ca <- data.frame(ca_model$calres$scores)

caplotly <- plot_ly(data = plotme_ca, x = ~Comp.1, y = ~Comp.2, z = ~Comp.3) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'Comp 1'),
                      yaxis = list(title = 'Comp 2'),
                      zaxis = list(title = 'Comp 3')))

caplotly

mg_model <- simca(mg_calibration, 'magnesium nitrate', ncomp = 3, cv = 1, center = TRUE, scale = TRUE)

plot(mg_model)
summary(mg_model)

plotme_mg <- data.frame(mg_model$calres$scores)

mgplotly <- plot_ly(data = plotme_mg, x = ~Comp.1, y = ~Comp.2, z = ~Comp.3) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'Comp 1'),
                      yaxis = list(title = 'Comp 2'),
                      zaxis = list(title = 'Comp 3')))

mgplotly

nitrates <- simcam(models = list(ca_model, mg_model), info = 'Nitrate SIMCA')
summary(nitrates)
plot(nitrates)
plotCooman(nitrates)
plotResiduals(nitrates)
plotDiscriminationPower(nitrates)
#plotModellingPower(nitrates)
#plotModelDistance(nitrates)