# missing scripts
detag <- function(thing, splitchar="_") {
  if (class(thing)=='character' && length(thing)==1) {
    return(strsplit(thing, splitchar)[[1]][2])
  }
  else if (class(thing) %in% c('character', 'list')) {
    return(sapply(thing, function (n) strsplit(n, splitchar)[[1]][2]))
  }
}

extract_pname <- function(thing) {
  if (class(thing)=='character' && length(thing)==1) {
    return(strsplit(thing, "-")[[1]][1])
  }
  else if (class(thing) %in% c('character', 'list')) {
    return(sapply(thing, function (n) strsplit(n, "-")[[1]][1]))
  }
}