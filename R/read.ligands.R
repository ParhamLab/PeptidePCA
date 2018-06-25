read.ligands = function(lig.dir, lig.names= NA) {
  lig.files = list.files(lig.dir, full.names = TRUE)
  print(lig.files)
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