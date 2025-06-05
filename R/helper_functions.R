#' check_id_names
#'
#' This function checks that identifier names in abundance data are actually present in geneset data
#'
#' @param abundance_data data object of class 'abundance_data'
#' @param geneset_data data object of class 'geneset_data'
#' @return True/false value if id names are correct
#'
check_id_names = function(abundance_data, geneset_data){
  if(!inherits(abundance_data, "abundance_data")) stop("abundance_data object must be of class 'abundance_data'")
  if(!inherits(geneset_data, "geneset_data")) stop("geneset_data object must be of class 'geneset_data'")

  row_names = rownames(abundance_data)
  geneset_names = unique(as.vector(geneset_data$matrix))
  overlap = intersect(row_names, geneset_names)

  if(all(row_names %in% geneset_names)){message("all identifiers in 'abundance_data' are present in 'geneset_data'")}
  else{message(paste("only", length(overlap), "identifiers in 'abundance_data' are present in 'geneset_data'", sep = " "))}

}


#' compare_id_names
#'
#' This function compares identifier names between two abundance data objects
#'
#' @param abundance_data1 first data object of class 'abundance_data'
#' @param abundance_data2 second data object of class 'abundance_data'
#' @return check number of identifiers in common
#'
#' @noRd
compare_id_names = function(abundance_data1, abundance_data2){
  if(!inherits(abundance_data1, "abundance_data") & !inherits(abundance_data2, "abundance_data")) stop("abundance_data1 and abundance_data2 objects must be of class 'abundance_data'")

  row_names1 = rownames(abundance_data1)
  row_names2 = rownames(abundance_data2)

  overlap = intersect(row_names1, row_names2)

  if(length(overlap) == 0){message("'abundance_data1' and 'abundance_data2' have zero identifiers in common")}
  else{message(paste("only", length(overlap), "identifiers in common between 'abundance_data1' and 'abundance_data2'", sep = " "))}
}


#' compare_column_names
#'
#' This function compares column names between two abundance data objects
#'
#' @param abundance_data1 first data object of class 'abundance_data'
#' @param abundance_data2 second data object of class 'abundance_data'
#' @return reports number of columns in common
#'
compare_column_names = function(abundance_data1, abundance_data2){
  if(!inherits(abundance_data1, "abundance_data") & !inherits(abundance_data2, "abundance_data")) stop("abundance_data1 and abundance_data2 objects must be of class 'abundance_data'")

  col_names1 = colnames(abundance_data1)
  col_names2 = colnames(abundance_data2)

  overlap = intersect(col_names1, col_names2)

  if(length(overlap) == 0){message("'abundance_data1' and 'abundance_data2' have zero column names in common")}
  else{message(paste("only", length(overlap), "column names in common between 'abundance_data1' and 'abundance_data2'", sep = " "))}
}


