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
