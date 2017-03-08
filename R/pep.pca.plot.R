#' Plot LAPP PCA
#'
#' \code{pep.pca.plot()} takes in a modified PCA object returned by pep.pca() and plots it.
#' Uses ggplot2.
#'
#' @param pca modified PCA object returned by pep.pca(). Can take up to 8 different HLA alleles.
#'
#' @param type Describes the type of plot to make.
#' Either 'points' (default) for scatter plo,t or 'density' for a 2D KDE density graph.
#'
#' @param dim1 Index of PC to plot along the x-axis. Defaults to 1
#'
#' @param dim2 Index of PC to plot along the y-axis. Defaults to 2
#'
#' @param type How the PCA should be plotted. Either 'points' for scatterplot or
#' 'density' for an automated 2D KDE plot (default).
#'
#' @param colors character that describes how each of the clusters should be colored.
#' Defaults to the colorblind palette proposed by Okabe and Ito (2002) http://jfly.iam.u-tokyo.ac.jp/color/
#'
#' @export

pep.pca.plot = function(pca, dim1 = 1, dim2 = 2, type = "density",
                        colors = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")) {

  PC.labels = paste0("PC.", 1:ncol(pca$ind$coord))

  pca.df = data.frame(pca$ind$coord)
  colnames(pca.df) = PC.labels
  pca.df$Peptide = as.character(pca$sequence)
  pca.df$HLA = pca$HLA

  palette = scale_color_manual(values = colors)

  # pca plot

  pca.plot= ggplot()

  pca.settings = ggplot() +
    scale_alpha_continuous(limits = c(0, 1e-05)) +
    labs(title = paste0("PCA for Peptides by HLA (", length(levels(pca$HLA)), " HLA, ", length(PC.labels), " features)"),
         x = paste0("PC", dim1, ": ", round(pca$eig$per[dim1], digits = 4), "%"),
         y = paste0("PC", dim2, ": ", round(pca$eig$per[dim2], digits = 4), "%")) +
    theme(title = element_text(size = 16, face = "bold", angle = 0),
          axis.text = element_text(size = 14, face = "bold", angle = 0),
          axis.line= element_blank(),
          panel.grid = element_blank(), panel.background = element_blank()) +
    palette +
    guides(color = guide_legend(label.position = "bottom",
                                title.theme = element_text(size = 18, face = "italic", angle = 0),
                                label.theme = element_text(size = 12, face = "bold", angle = 0),
                                override.aes = list(size = 3)))

  if (type == "points") {

    # scramble rows of PC table so that plotting them looks nicer
    set.seed(129)
    pca.df = pca.df[sample(1:nrow(pca.df), nrow(pca.df)), ]

    pca.plot= pca.settings +
      geom_point(data = pca.df, mapping = aes(x = pca.df[, dim1], y = pca.df[, dim2], color = HLA),
                 size = 4, alpha = 0.5)

  } else if (type == "density") {

    pca.plot= pca.settings +
      geom_density2d(data = pca.df, mapping = aes(x = pca.df[, dim1], y = pca.df[, dim2], color = HLA),
                     size = 1, alpha = 0.8)

  }
  return(pca.plot)
}

