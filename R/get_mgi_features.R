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
#' @importFrom dplyr copy_to src_sqlite
#' @importFrom readr read_tsv
mgigenes_gz2sql <- function(gene_tbl = "mgi_gene",
                            sql_file = "mgi_db.sqlite",
                            gz_file = file.path("ftp://ftp.jax.org",
                                                "SNPtools",
                                                "genes",
                                                "MGI.20130703.sorted.txt.gz")) {
  ## Create or use existing SQLite database.
  my_db <- dplyr::src_sqlite(sql_file,
                      create = !file.exists(sql_file))
  ## Read gzipped file from MGI. Names not included.
  mgi_gene <- readr::read_tsv(gz_file, comment="#", col_names=FALSE)
  names(mgi_gene) <- c("chr", "source", "type", "start", "stop",
                       "score", "strand", "phase", "ID", "Name",
                       "Parent", "Dbxref", "mgiName", "bioType")
  invisible(dplyr::copy_to(my_db,
                   mgi_gene,
                   name = gene_tbl,
                   temporary = FALSE,
                   indexes = list(names(mgi_gene))))
}

#' Pull MGI gene tbl from SQLite database
#'
#' Get subset of MGI gene table from SQLite using chromosome interval.
#' This is being replaced by \code{query_variants}. See package qtl2db.
#'
#' @param chr_id chromosome number
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
#' @importFrom dplyr collect filter rename src_sqlite
#' 
get_mgi_features <- function(chr_id, start_bp, end_bp,
                          with_name = TRUE,
                          gene_tbl = "mgi_gene",
                          sql_file = "mgi_db.sqlite") {
  check_interval(chr, start_bp, end_bp)
  start_bp <- convert_bp(start_bp)
  end_bp <- convert_bp(end_bp)
  
  # kludge for now; make sure first column has name "chr"
  out <- tbl(dplyr::src_sqlite(sql_file), gene_tbl)
  chk <- dplyr::collect(head(out, 1))
  if(names(chk)[1] == "seqid")
    out <- dplyr::rename(out, chr = seqid)
  
  out <- dplyr::filter(out,
                       chr == chr_id,
                       start >= start_bp,
                       stop <= end_bp)
  if(with_name)
    out <- dplyr::filter(out, !is.na(Name))

  ## Collect from database and add class.
  out <- dplyr::collect(out)
  class(out) <- c("feature_tbl", class(out))
  out
}
