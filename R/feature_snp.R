#' Match features with SNPs
#'
#' Find features that overlap with SNPs
#'
#' @param snp_tbl tbl of SNPs from \code{assoc.map}
#' @param feature_tbl tbl of feature information from \code{\link{get_mgi_features}}
#' @param extend extend region for SNPs in bp (default 5000)
#'
#' @return tbl of features covering SNPs
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @examples
#' \dontrun{get_feature_snp(snp_tbl, feature_tbl)}
#'
#' @export
#' @importFrom dplyr arrange bind_rows distinct filter group_by inner_join 
#' mutate select summarize ungroup
#'
get_feature_snp <- function(snp_tbl, feature_tbl, extend=5000) {
  snp_tbl$pos <- convert_bp(snp_tbl$pos)

  feature_tbl <- dplyr::filter(feature_tbl, 
                               type %in% c("exon","gene","pseudogene","pseudogenic_exon"))
  if(!nrow(feature_tbl))
    return(NULL)

  ## Pull features that overlap with SNPs.
  ## Group features by Dbxref for each gene and get bp range.
  genes <- dplyr::ungroup(
    dplyr::summarize(
      dplyr::group_by(feature_tbl, Dbxref), 
      min_bp = min(start) - extend,
      max_bp = max(stop) + extend))
  if(!nrow(genes))
    return(NULL)
  
  ## For each gene, get rows of snp_tbl in bp range.
  tmpfn <- function(x,snp_tbl) {
    keep <- (x[1] <= snp_tbl$pos) & (x[2] >= snp_tbl$pos)
    snp_tbl[keep,]
  }
  tmp <- apply(as.matrix(genes[,2:3]), 1, tmpfn, snp_tbl)
  names(tmp) <- genes$Dbxref
  tmp <- dplyr::bind_rows(tmp, .id = "Dbxref")
  out <- dplyr::distinct(
    dplyr::arrange(
      dplyr::select(
        dplyr::mutate(
          dplyr::inner_join(feature_tbl, tmp, by = "Dbxref"),
          SNP = snp_id), 
        -snp_id),
      Name, type), 
    SNP, start, stop, strand, type, .keep_all=TRUE)

  if(!nrow(out))
    return(NULL)

  class(out) <- c("feature_snp", class(out))
  out
}

#' Summary of features with SNP information
#'
#' @param object tbl of feature information from \code{\link{get_feature_snp}}
#' @param ... additional parameters ignored
#'
#' @return tbl of feature summaries by type
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @method summary feature_snp
#' @rdname feature_snp
#' @export
#' @importFrom dplyr group_by summarize ungroup
#' 
summary.feature_snp <- function(object, ...) {
  dplyr::ungroup(
    dplyr::summarize(
      dplyr::group_by(object, type), 
      count = n(),
      minbp = min(start),
      maxbp = max(stop)))
}

