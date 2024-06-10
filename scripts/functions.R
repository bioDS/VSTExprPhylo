rename_categories <- function(data, expected_char_num) {
  categories=sort(unique(as.vector(unlist(data))))
  
  char_num=length(categories)-1
  
  if (char_num < expected_char_num-1) {
    print(paste("Renaming chategories as only the following categories are present:", toString(categories), sep=" "))
    categories = sort(categories)
    replace=0:char_num
    replace = replace[seq_along(categories)]
    data = phyloRNA::replace(data, categories, replace)
  } else {
    print("The number of categories as expected, nothing to be done.")
  }
  
  data
}

remove_zeroexrp_cells_genes <- function(data, sample) {
    n_cells <- ncol(data)
    data <- data[colSums(data)>0]
    data <- data[rowSums(data)>0,]
    
    if (ncol(data) < n_cells) {
        print(paste("WARNING: Some spots were removed in sample",
        sample, "with", as.character(ncol(data)), "remaning spots.", sep=" "))
    }
    
    print(paste("The number of genes in sample", sample, "is", as.character(nrow(data)), "after removing zero-expressed genes.", sep = " "))
    
    data
    
}
