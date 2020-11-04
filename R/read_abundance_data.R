#' read_abundance_data
#'
#' This function reads in abundance data
#'
#' @param filename the directory path to a file containing the abundance data
#' @param data_type a character string, either 'protdata' or 'phosphoprotdata' which describes the type of abundance data
#' @param infocol_name a character string describing the name of the column extra identifier info in abundance data of type 'phosphoprotdata'
#'
#'
#' @export
#'
read_abundance_data = function(filename, data_type, infocol_name = NULL){

  #check that data_type is one of the appropriate types
  if(!is.element(data_type, c("protdata", "phospho_protdata"))) stop("data_type must be one of 'protdata', 'phospho_protdata'")
  if(data_type == "phospho_protdata" & is.null(infocol_name)) stop("'infocol_name must be provided with 'phospho_protdata' objects")

  df = read.table(filename, sep="\t", row.names=1, header=1, stringsAsFactors=F)

  #assign "data_type" attribute to df data object
  if(data_type == "phospho_protdata" & !is.null(infocol_name)){
    if(infocol_name %in% colnames(df)){attr(df, "infocol_name") = infocol_name}
    else stop("infocol_name is not one of the column names in the data input through 'filename'")

    attr(df, "data_type") = data_type

  }

  class(df) = c("abundance_data", "data.frame")

  return(df)

}


