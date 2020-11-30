#' @export
`%nin%` <- Negate("%in%")


mismatch_df = function (x, y, on = NULL){
  if (is.null(on)) {
    on <- dplyr::intersect(names(x), names(y))
    message("Matching on: ", paste(on, collapse = ", "))
  }
  keys <- plyr::join.keys(x, y, on)
  x[keys$x %nin% keys$y, , drop = FALSE]
}

#' @importFrom utils download.file
download_data = function(path, is.rds = F){
  ftp_base = "ftp://ftp.ebi.ac.uk/pub/contrib/petsalaki/Wadie_et_al"
  if (!is.rds){
    tmp = tempfile()
    data = download.file(file.path(ftp_base, path), destfile = tmp)
    return(tmp)
  }
  else{
    data = readRDS(url(file.path(ftp_base, path), "rb"))
    return(data)
  }
}
