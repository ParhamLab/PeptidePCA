# GDrive is where my Google drive is mapped in my home directory
# work is my subfolder on Google drive 
# hugo is the folder in work that I shared with you
setwd("~/GDrive/Illing") 
library("ggplot2", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
library("FactoMineR", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
library("PeptidePCA", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
library(RWebLogo)
library(stats)
source('pep.pca.plot.R')
source('read.ligands.R')
colors = c( "#E69F00", "#0072B2", "#CC79A7", "#009E73")

# I put the data in pepties directory in hugo
for (l in 9:11) {
  
  lmers = paste(toString(l), "mers", sep="")
  
  # for example "peptides/9mers"
  directory = paste("peptides/", lmers, sep="")
  print(directory)
  ligands.B57 = read.ligands(directory)
  features.B57 = conv.features.list(ligands.B57, convmat.24, l)
  features.B57.pca = pep.pca(features.B57)
  
  pep.pca.plot(pca=features.B57.pca, type="density", colors=colors)
  ggsave(paste("pdf/", lmers, "_pca_density.pdf", sep=""))
  
  pdf(paste("pdf/", lmers, "_pca_vars.pdf", sep=""))
  plot(features.B57.pca,choix="var",select="contrib 10")
  dev.off()
  
  pca = features.B57.pca
  df = data.frame(pca$ind$coord)
  sizes = list()
  pca_blank = pca
  levels(pca_blank$HLA) <- c(levels(pca_blank$HLA), 'all_alleles')
  pca_blank$HLA[] <- 'all_alleles'
  
  for (d in 2) {
    for (k in 2:4) {
      # run clustering on subset only
      res = kmeans(df[, 1:d], k, nstart=10, iter.max=100)
      sequences = pca$sequence
      assignments = res$cluster
      pep.pca.plot(pca=features.B57.pca, type="density", colors=colors)
      for (hla in unique(pca$HLA)) {
        cluster_size_file = paste('clusters/cluster_sizes_', lmers,'_',  toString(k), '_clusters_usingPC1-', toString(d),
                                  '_for_', hla, '.txt', sep="")
        pca_copy = pca
        # get hla subset
        others = which(pca$HLA != hla)
        cat(c(hla, '\n'), file=cluster_size_file)
        for (c in 1:k) {
          subset = which(assignments==c & pca$HLA==hla)
          hla_cluster_name = paste(hla, '_cluster_', c, sep="")
          levels(pca_copy$HLA)<-c(levels(pca_copy$HLA), hla_cluster_name)
          cluster_name = paste('cluster_', toString(c), '_of_', toString(k), '_usingPC1-', toString(d),
                               '_for_', hla, sep="")  
          cluster_sequences = sequences[subset]
          seq_hla_df = data.frame(seq=as.character(cluster_sequences))
          write.csv(seq_hla_df, file=paste('clusters/', cluster_name, sep=""))
          weblogo(as.character(cluster_sequences), size='medium', open=FALSE,
                  file.out=paste('pdf/logos/', cluster_name, '.pdf', sep=""))
          pca_copy$HLA[subset] <- hla_cluster_name
          cat(c('cluster:', toString(c), toString(length(subset)), ' '),
              file=cluster_size_file, append=TRUE)
        }
        levels(pca_copy$HLA) <- c(levels(pca_copy$HLA), 'other_alleles')
        pca_copy$HLA[others] <- 'other_alleles'
        pep.pca.plot(pca=pca_copy, type="density")
        ggsave(paste("pdf/densities/", toString(k), '_clusters_usingPC1-', toString(d),
                                       '_for_', hla, "_pca_density.pdf", sep=""))

        pep.pca.plot(pca=pca_copy, type='both')
        ggsave(paste("pdf/points/", toString(k), '_clusters_usingPC1-', toString(d),
                     '_for_', hla, "_pca_points.pdf", sep=""))
        
        cat('\n', file=cluster_size_file, append=TRUE)
      }
    }
  }
}
  
  