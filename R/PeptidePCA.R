#' PeptidePCA
#'
#' PeptidePCA contains a small suite of functions that can be used to recreate PCA figures from Hilton et al. (2017).
#' HLA peptide sequences can be transformed into new features using an amino acid conversion matrix and then plotted
#' using PCA. The packages ggplot2 and FactoMineR are required.
#'
#' @name LAPP
#' @docType package
#' @import ggplot2 FactoMineR
NULL

#' Amino acid conversion matrix with 24 properties
#'
#' Conversion matrix used to transform n-mer ligand sequences into 24 * n position-by-property features.
#' 4 properties are physicochemical: molecular weight, surface area, hydrophobicity index, and isoelectric point.
#' 20 properties describe identity of the amino acid at a position.
#'
#' @docType data
#' @keywords datasets
#' @name convmat.24
#' @usage data(convmat.24)
#' @format A data frame with 20 rows (amino acids) and 24 columns (properties)
#' @references \url{to do}
NULL

#' Amino acid conversion matrix with 4 properties
#'
#' Conversion matrix used to transform n-mer ligand sequences into 4 * n position-by-property features.
#' 4 properties are physicochemical: molecular weight, surface area, hydrophobicity index, and isoelectric point.
#'
#' @docType data
#' @keywords datasets
#' @name convmat.4
#' @usage data(convmat.4)
#' @format A data frame with 20 rows (amino acids) and 4 columns (properties)
NULL

#' Ligand sequences from 4 HLA-A alleles
#'
#' List of list of sequences from HLA-A*01:01, HLA-A*02:01, HLA-A*03:01, and HLA-A*11:01.
#' All sequences from obtained from the publically shared sequences from IEDB on November 2016.
#'
#' @docType data
#' @keywords datasets
#' @name ligands.4A
#' @usage data(ligands.4A)
#' @format A list of 4 character lists
NULL

