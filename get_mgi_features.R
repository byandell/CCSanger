#' Convert gzipped file of MGI genes to SQLite database
#'
#' Read file, filter to MGI identified genes and copy to database.
#' If database does not exist, create it.
#' The gzipped gene file is 24.8Mb,or 368.2Mb uncompressed.
#' The R object with only type=="gene" and known Name is 10.5Mb.
#' Storing in SQLite, file is 1.5Mb.
#'
#' @param gene_tbl character string of tbl to be used in database
#' @param sql_file character string name of SQLite file
#' @param gz_file character string path to MGI gene file
#'
#' @return tbl collected from SQLite database (invisible)
#'
#' @examples
#' dontrun(mgigenes_gz2sql("mgi_gene"))
#'
#' @export
mgigenes_gz2sql <- function(gene_tbl = "mgi_gene",
                            sql_file = "mgi_db.sqlite",
                            gz_file = file.path("ftp://ftp.jax.org",
                                                "SNPtools",
                                                "genes",
                                                "MGI.20130703.sorted.txt.gz")) {
  ## Create or use existing SQLite database.
  my_db <- src_sqlite(sql_file,
                      create = !file.exists(sql_file))
  ## Read gzipped file from MGI. Names not included.
  mgi_gene <- read_tsv(gz_file, comment="#", col_names=FALSE)
  names(mgi_gene) <- c("seqid", "source", "type", "start", "stop",
                       "score", "strand", "phase", "ID", "Name",
                       "Parent", "Dbxref", "mgiName", "bioType")
  invisible(copy_to(my_db,
                   mgi_gene,
                   name = gene_tbl,
                   temporary = FALSE,
                   indexes = list(names(mgi_gene))))
}

#' Pull MGI gene tbl from SQLite database
#'
#' Get subset of MGI gene table from SQLite using chromosome interval.
#'
#' @param chr chromosome number
#' @param start_bp starting position in bp
#' @param end_bp ending position in bp
#' @param with_name only pull named features if TRUE (default)
#' @param gene_tbl character string of tbl to be used in database
#' @param sql_file character string name of SQLite file (default "mgi_db.sqlite")
#'
#' @return tbl from SQLite database
#'
#' @examples
#' dontrun(get_mgi_features(9, 104*1e6, 109*1e6))
#'
#' @export
get_mgi_features <- function(chr, start_bp, end_bp,
                          with_name = TRUE,
                          gene_tbl = "mgi_gene",
                          sql_file = "mgi_db.sqlite") {
  check_interval(chr, start_bp, end_bp)
  start_bp <- convert_bp(start_bp)
  end_bp <- convert_bp(end_bp)
  out <- tbl(src_sqlite(sql_file), gene_tbl) %>%
    filter(seqid == chr,
           start >= start_bp,
           stop <= end_bp)
  if(with_name)
    out <- out %>% filter(!is.na(Name))

  ## Collect from database and add class.
  out <- collect(out)
  class(out) <- c("feature_tbl", class(out))
  out
}
