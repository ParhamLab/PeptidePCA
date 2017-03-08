#' Perform PCA on HLA features
#'
#' \code{pep.pca()} takes in a data.frame from conv.features.list() and performs PCA.
#' The results are a modified form of results from PCA() from the FactoMineR package.
#'
#' @param features data.frame object returned by make.features.list().
#'
#' @export

pep.pca = function(features) {
  features.data = subset(features, select = -c(sequence, HLA))

  # drop any constant features
  drop.index = apply(features.data, 2, var, na.rm = TRUE) == 0
  features.drop = features.data[, !drop.index]

  # manually scale features
  features.scale = scale(features.drop)

  # perform PCA with FactoMineR
  pca = PCA(features.scale, scale.unit = FALSE, ncp = ncol(features.scale), graph = FALSE)

  # create loading matrix
  pca.loading = sweep(pca$var$coord, 2, sqrt(pca$eig[1:ncol(pca$var$coord), 1]), FUN = "/")

  # add additional items to output
  pca$sequence = features$sequence
  pca$HLA = features$HLA
  pca$drop.inds = drop.index
  pca$centering = attr(features.scale, "scaled:center")
  pca$scaling = attr(features.scale, "scaled:scale")
  pca$loading = pca.loading

  return(pca)
}
