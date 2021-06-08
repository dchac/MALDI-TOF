###install all listed packages before running script, libraries only need to be loaded at the beginning of any R session
###-----------------------------------------------------------------------
#library(devtools)
#install_github("vqv/ggbiplot")

library(ggbiplot)
library(vegan)
library(ecodist)
library(labdsv)
library(ape)
library(ade4)
library(smacof)
library(tidyverse)
library(cluster)
library(factoextra)
library(dendextend)
library(ggplot2)
library(rgl)
library(pca3d)
library(signal)
library(rlist)
library(rwt)
library(MALDIquant)
library(MALDIquantForeign)
library(pvclust)
library(FactoMineR)


### read in spectra data files
###---------------------------------------------------------------------



#Set working directory where all PMF files are located.
#File naming convention is as follows "species_condition_biological.replicate_technical.replicate"

## This command makes a list of sample names from files. 

filenames <- list.files( pattern="*.txt", full.names=FALSE)
samples   <- tools::file_path_sans_ext(basename(filenames)) 



### convert and concatenate data to class MassSpectra from maldiquant package
###---------------------------------------------------------------------

data <- c()
for (i in 1:length(filenames))
{
  spectrum    <- import( filenames[i] )
  data[[i]] <- spectrum
  
}


### collapse sample names by technical or biological replicate number for spectra averaging
###---------------------------------------------------------------------

technical_reps <- c()

for (i in 1:length(samples))
{
  temp <- unlist(strsplit(samples[i],"_"))
  temp <- paste(temp[1:3], collapse = "_") #change to indices representing the data you want
  technical_reps[i] = temp#eg. [1:3] to average over technical replicates and [1:2] to average over biological replicates
}


### create data.frame with relevant sample data for plotting 
###---------------------------------------------------------------------

meta_data <- data.frame(species = c(), 
                        cond = c()#, 
                        #bioRep = c()  #comment out if averaging over biological replicates
                        )  

tech_reps <- unique(technical_reps)

for (i in 1:length(tech_reps))
{
  temp <- unlist(strsplit(tech_reps[i],"_"))
  temp <- data.frame(species = temp[1], 
                     cond = temp[2], 
                     bioRep = temp[3] #comment out if averaging over biological replicates
                    )
  meta_data = rbind(meta_data, temp)
}


### spectral analysis on PMF data 
###---------------------------------------------------------------------

data <- lapply(data, function (ll) ll[[1]])
names(data) <- technical_reps

any(sapply(data, isEmpty))

spectra <- trim(data, range = c(3000, 15000))
spectra <- transformIntensity(spectra, method = "sqrt")
spectra <- smoothIntensity(spectra, method = "SavitzkyGolay", halfWindowSize = 25)
#spectra <- removeBaseline(spectra, method = "TopHat", halfWindowSize = 5)
spectra <- removeBaseline(spectra, method = "SNIP", iterations = 100)
spectra <- calibrateIntensity(spectra, method = "TIC")
spectra <- alignSpectra(spectra, SNR = 3, tolerance = 0.001, warpingMethod = "lowess")
spectra <- averageMassSpectra(spectra, labels = technical_reps, method = "mean")
peaks   <- detectPeaks(spectra, method = "MAD", SNR = 3)
peaks   <- binPeaks(peaks, tolerance = 0.001)
names(peaks) <- unique(technical_reps)
#names(peaks) <- samples


### convert data class MassPeaks to class matrix for PCA/dendrogram analysis and plotting 
###---------------------------------------------------------------------

dataMatrix <- intensityMatrix(peaks, spectra)


### calculate PCA
###---------------------------------------------------------------------

pca.data <- prcomp(dataMatrix, center = TRUE, scale. = TRUE)

### plot PCA with 95% CI ellispe on clusters
###---------------------------------------------------------------------

summ = summary(pca.data)
pc1  = round(summ$importance[2,1]*100, digits = 1)
pc2  = round(summ$importance[2,2]*100, digits = 1)

g <- ggbiplot(pca.data, obs.scale = 1, choices = 1:2, var.scale = 1, 
              groups = meta_data$cond, # change to meta_data$... based on desired cluster classification (eg cond, species)
              ellipse = TRUE, 
              var.axes = FALSE, 
              ellipse.prob = 0.95,
              shape = 3)
g <- g + scale_color_discrete(name = '')
g <- g + scale_shape_discrete(name = '')
g <- g + geom_point(aes(colour = meta_data$cond, shape = meta_data$cond), size = 4)
g <- g + theme_bw(base_size = 24)
g <- g + theme(legend.title = element_blank()) + labs(title = "", x = paste("PC1 (",pc1,"% explained var.)", sep = ""), y = paste("PC2 (",pc2,"% explained var.)", sep = ""))
g <- g + theme(aspect.ratio = 1, legend.position = c(0.15,0.9), legend.background = element_rect(linetype = "solid", color = "black"), legend.key.height = unit(1, "cm")) 
print(g)



### plot PCA with 95% CI ellispe on cluster means
###---------------------------------------------------------------------

#meta_data$cond <- as.character(meta_data$cond)
#meta_data$cond[meta_data$cond == "DH5a"] <- "DH5-alpha"
#meta_data$cond[meta_data$cond == "EPEC"] <- "Enteropathogenic E. coli"

summ = summary(pca.data)
pc1  = round(summ$importance[2,1]*100, digits = 1)
pc2  = round(summ$importance[2,2]*100, digits = 1)
g <- fviz_pca_ind(pca.data, 
                  col.ind = meta_data$cond, # change to meta_data$... based on desired cluster classification (eg cond, species)
                  addEllipses = TRUE, 
                  geom.ind = "point", pointsize = 4,
                  ellipse.type = "confidence",
                  ellipse.level = 0.95,
                  ggtheme = theme_bw(base_size = 24))

g <- g + theme(legend.title = element_blank()) + labs(title = "", x = paste("PC1 (",pc1,"% explained var.)", sep = ""), y = paste("PC2 (",pc2,"% explained var.)", sep = ""))

print(g)



### calculate and plot dendrogram
###---------------------------------------------------------------------

data.dis <- dist(dataMatrix, 'euclidian')
labels(data.dis) <- names(peaks)

dend <- data.dis %>%
  hclust(method = "ward.D") %>%
  as.dendrogram %>%
  set("labels_cex", 0.5) %>%
  set("labels_col", k = 3) %>%
  set("branches_k_color", k = 3) %>%
  plot
#change k=# to indicate number of colored branches

statMat <- t(dataMatrix)
result <- pvclust(scale(na.omit(statMat)), method.dist = "euclidian", method.hclust = "ward.D", nboot = 10000)
plot(result, print.pv = TRUE, print.num = FALSE)
pvrect(result, alpha = 0.95)

#trtlab <- c(37, 26)[factor(substr(labels(dend),1,1))]
#trtlab
dend <- result$hclust %>%
  as.dendrogram %>%
  set("labels_cex", 0.5) %>%
  set("labels_col", k = 2) %>%
  set("branches_k_color", k = 3) %>%
  set("leaves_pch", 19) %>%
  set("leaves_col", meta_data$cond) %>%
  set("branches_lwd", 3) %>%
  plot (horiz = TRUE)


pvrect(result, alpha = 0.95)


