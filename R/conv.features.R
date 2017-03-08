#' Create converted feature matrix using a list of ligands and a conversion table
#'
#' \code{conv.features()} takes in a list of ligands from one HLA and returns
#' a converted feature matrix of the HLA's ligands using a conversion matrix.
#'
#' @param ligands List of HLA ligand sequences from one HLA allele.
#'
#' @param conv.table Feature conversion table created by make.conv.table().
#'
#' @param n Length of peptides in ligands to convert into features.
#'
#' @export

conv.features = function(ligands, conv.table, n) {
  n.feat = ncol(conv.table)

  nmers = ligands[which(sapply(ligands, nchar) == n)]

  n.seq = length(nmers)

  # blank matrix of properties to preallocate space
  features = as.data.frame(matrix(nrow = n.seq, ncol = n.feat * n))

  # for every character in the ligand...
  for (char in 1:n) {
    # for every AA feature...
    for (feat in 1:n.feat) {
      # blank vector of new pos by prop feature to preacllocate space
      new.feature = vector(mode = "numeric", length = n.seq)
      # for every peptide...
      for (seq in 1:n.seq) {
        AA.index = match(substr(nmers[seq], char, char), rownames(conv.table))
        new.feature[seq] = conv.table[AA.index, feat]
      }
      # add new.feature to features
      feat.ind = (char - 1) * n.feat + feat
      features[, feat.ind] = new.feature
      colnames(features)[feat.ind] = paste0("AA", char, ".", colnames(conv.table)[feat])
    }
  }
  # add rownames to pep matrix
  rownames(features) = rownames(ligands)

  return(features)
}
