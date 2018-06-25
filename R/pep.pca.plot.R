pep.pca.plot = function (pca,
                         dim1 = 1,
                         dim2 = 2,
                         type = "density",
                         colors = c(
                           "#000000",
                           "#E69F00",
                           "#56B4E9",
                           "#009E73",
                           "#F0E442",
                           "#0072B2",
                           "#D55E00",
                           "#CC79A7"
                         ))
{
  PC.labels = paste0("PC.", 1:ncol(pca$ind$coord))
  pca.df = data.frame(pca$ind$coord)
  colnames(pca.df) = PC.labels
  pca.df$Peptide = as.character(pca$sequence)
  pca.df$HLA = pca$HLA
  
  subset = which(pca$HLA != "other_alleles")
  pca.subset.df = data.frame(pca$ind$coord[subset,])
  colnames(pca.subset.df) = PC.labels
  pca.subset.df$Peptide = as.character(pca$sequence[subset])
  pca.subset.df$HLA = pca$HLA[subset]
  
  
  palette = scale_color_manual(values = colors)
  pca.plot = ggplot()
  pca.settings = ggplot() + scale_alpha_continuous(limits = c(0, 1e-05)) +
    labs(
      title = paste0(
        "PCA for Peptides by HLA (",
        length(levels(pca$HLA)),
        " HLA, ",
        length(PC.labels),
        " features)"
      ),
      x = paste0("PC", dim1, ": ", round(pca$eig[, 2][dim1], digits = 4), "%"),
      y = paste0("PC", dim2, ": ", round(pca$eig[, 2][dim2], digits = 4), "%")
    ) +
    theme(
      title = element_text(
        size = 16,
        face = "bold",
        angle = 0
      ),
      axis.text = element_text(
        size = 14,
        face = "bold",
        angle = 0
      ),
      axis.line = element_blank(),
      panel.grid = element_blank(),
      panel.background = element_blank()
    ) +
    palette + guides(
      color = guide_legend(
        label.position = "bottom",
        title.theme = element_text(
          size = 18,
          face = "italic",
          angle = 0
        ),
        label.theme = element_text(
          size = 12,
          face = "bold",
          angle = 0
        ),
        override.aes = list(size = 3)
      )
    )
  if (type == "points") {
    set.seed(129)
    pca.df = pca.df[sample(1:nrow(pca.df), nrow(pca.df)),]
    pca.plot = pca.settings + geom_point(
      data = pca.df,
      mapping = aes(x = pca.df[, dim1], y = pca.df[, dim2],
                    color = HLA),
      size = 4,
      alpha = 0.5
    )
  }
  else if (type == "density") {
    pca.plot = pca.settings + geom_density2d(
      data = pca.df,
      mapping = aes(x = pca.df[, dim1], y = pca.df[, dim2],
                    color = HLA),
      size = 1,
      alpha = 0.8
    )
  }
  else if (type == "both") {
    print("background")
    pca.plot = pca.settings + geom_density2d(
      data = pca.df,
      mapping = aes(x = pca.df[, dim1], y = pca.df[, dim2], color='HLA'),
      size = 1,
      alpha = 0.8
    )
    print("points")
    pca.plot = pca.plot + 
      geom_point(
        data = pca.subset.df,
        mapping = aes(x = pca.subset.df[, dim1], y = pca.subset.df[, dim2],
                      color = HLA),
        size = 1,
        alpha = 0.5
      )
  }
  return(pca.plot)
}