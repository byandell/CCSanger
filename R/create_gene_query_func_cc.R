#' @export
create_gene_query_func_cc <- function(dbfile) {
  function(chr, start, stop) {
    # This only gets features with names, which are genes. 
    CCSanger::get_mgi_features(chr, start, stop, with_name = FALSE,
                               sql_file = dbfile)
  }
}