#' Read in amino acid table
#'
#' \code{make.conv.table()} takes in a csv file to create a matrix of amino acid feature to be used in conversion
#'
#' @param AA.file csv file to be converted into the amino acid conversion matrix. See given example files for the format of these input files.
#'
#' @export

make.conv.table = function(conv.table.file) {
  # first column should be AA symbols
  conv.table = read.csv(conv.table.file, row.names = 1, check.names= FALSE)
  return(conv.table)
}
