#' Create converted feature matrix using a list of list of ligands and a conversion matrix
#' Uses conv.features()
#'
#' \code{conv.features.list()} takes in a file directory and returns a list of lists of each HLA's ligands
#'
#' @param ligands.list List of list of HLA allele ligand sequences created by read.ligands().
#'
#' @param conv.table Feature conversion matrix created by read.AA.table().
#'
#' @param n Length of peptides in ligands to convert into features.
#'
#' @export

conv.features.list = function(ligands.list, conv.table, n) {

  features.total = data.frame()

  for (hla in 1:length(ligands.list)) {
    ligands = ligands.list[[hla]]
    nmers = ligands[which(sapply(ligands, nchar) == n)]
    features.allele = conv.features(nmers, conv.table, n)  # helper fxn
    features.named = cbind(data.frame(sequence = nmers, HLA = rep(names(ligands.list)[hla], length(nmers))), features.allele)
    features.total = rbind(features.total, features.named)
    print(paste(names(ligands.list[hla]), " feature conversion completed"))
  }
  return(features.total)
}
