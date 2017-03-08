#' Read in all ligand files in a folder
#'
#' \code{read.ligands()} takes in a file directory and returns a list of lists of each HLA's ligands
#'
#' @param lig.dir file directory with one or more files of txt files of peptide sequences.
#'      Each file is named after the HLA allele, and one peptide sequence is present per line.
#'
#' @param lig.names character list of HLA names for the ligand files. If not supplied, the names
#' of the .txt files will be used.
#'
#' @export

read.ligands = function(lig.dir, lig.names= NA) {
  lig.files = list.files(lig.dir, full.names = TRUE)
  ligands = lapply(lig.files, scan, what = " ", sep = "\n")

  # erase duplicated ligands
  for (i in 1:length(ligands)) {
    ligands[[i]] = unique(ligands[[i]])
  }

  # erase ligands with inappropriate characters
  forbidden= c("B", "J", "O", "U", "X", "Z")
  valid.chars = LETTERS[!LETTERS %in% forbidden]
  for (i in 1:length(ligands)) {
    ligands[[i]] = ligands[[i]][which(
      sapply(ligands[[i]], FUN = function(x)
        sum(unlist(strsplit(x, "")) %in% valid.chars) == nchar(x)))]
  }

  # name the ligands
  if(is.na(lig.names)){
    # drop last 4 characters ('.txt') to create names
    names(ligands) = gsub(".{4}$", "", list.files(lig.dir))

  } else {
    names(ligands) = lig.names
  }

  return(ligands)
}
