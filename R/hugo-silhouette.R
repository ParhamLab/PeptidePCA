# GDrive is where my Google drive is mapped in my home directory
# work is my subfolder on Google drive 
# Illing is the folder you shared with me
setwd("~/GDrive/Illing") 
library("ggplot2", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
library("FactoMineR", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
library("PeptidePCA", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
library(RWebLogo)
library(stats)
library(cluster)
source('pep.pca.plot.R')
source('read.ligands.R')
colors = c( "#E69F00", "#0072B2", "#CC79A7", "#009E73")
drawlogos = FALSE
# I put the data in pepties directory in hugo
pep_silhouette = function (df, sequences, filename_prefix) {
  for (d in 2:5) {
    D = daisy(df[,1:d])
    sis = c()
    for (k in 2:10) {
      res = kmeans(df[,1:d], k, nstart=10, iter.max=100)
      filename = paste(filename_prefix, '_', toString(k), '_using_PC1-', toString(d), '.pdf', sep="")
      pdf(filename)
      si <- silhouette(res$cluster, D)
      plot(si, col=1:k, border=NA)
      sis <- append(sis, summary(si)$si.summary["Mean"])
      dev.off()
    }
    pdf(paste(filename_prefix, '_summary_using_PC1-',toString(d), '.pdf', sep=""))
    plot(x=2:10, y=sis)
    title(main=paste('Silhouettes for different k in k-means (PC1-', toString(d), ')', sep=""), 
          xlab='number of clusters (k)', ylab='Silhouette')
    dev.off()
  }
}
for (l in 9:11) {
  
  lmers = paste(toString(l), "mers", sep="")
  # for example "peptides/9mers"
  directory = paste("peptides/", lmers, sep="")
  print(directory)
  ligands.B57 = read.ligands(directory)
  features.B57 = conv.features.list(ligands.B57, convmat.24, l)
  features.B57.pca = pep.pca(features.B57)
  pca = features.B57.pca
  df = data.frame(pca$ind$coord)
  for (hla in unique(pca$HLA)) {
    # get hla subset
    subset = which(pca$HLA == hla)
    # extract hla subset data
    df_subset = df[subset, ]
    # collect hla subset sequences
    sequences = pca$sequence[subset]        
    filename_prefix = paste('pdf/silhouette/', lmers, '_', hla, '_silhouette', sep='')
    pep_silhouette(df_subset, sequences, filename_prefix)
  }
  filename_prefix = paste('pdf/silhouette/', lmers, '_all','_silhouette', sep='')
  pep_silhouette(df, pca$sequence, filename_prefix)
}

