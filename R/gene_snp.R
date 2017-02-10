#' Match genes with SNPs
#'
#' Find features that overlap with SNPs
#'
#' @param snp_tbl tbl of SNPs from \code{assoc.map}
#' @param feature_tbl tbl of feature information from \code{\link{get_mgi_features}}
#' @param feature_snp tbl of feature information from \code{\link{get_feature_snp}}
#'
#' @return tbl of genes covering SNPs
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords hplot
#'
#' @examples
#' \dontrun{get_gene_snp(snp_tbl, feature_tbl)}
#'
#' @export
get_gene_snp <- function(snp_tbl, feature_tbl,
                         feature_snp =
                           get_feature_snp(snp_tbl, feature_tbl, 0)) {
  out <- feature_snp %>%
    filter(type=="gene" & !is.na(Name)) %>%
    distinct(SNP, start, stop, strand, .keep_all=TRUE) %>%
    mutate(gene=Name) %>%
    select(SNP,pos,lod,gene,start,stop,strand) %>%
    mutate(in_gene=(start<=pos & pos<=stop))
  
  class(out) <- c("gene_snp", class(out))
  out
  
}
#' Summary of genes overlapping SNPs
#'
#' @param object tbl of feature information from \code{\link{get_feature_snp}}
#' @param ... additional parameters ignored
#'
#' @return tbl of feature summaries by type
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @method summary gene_snp
#' @rdname gene_snp
#' @export
summary.gene_snp <- function(object, ...) {
  object %>%
    group_by(gene, start, stop) %>%
    summarize(snp_count=n(),
              min_lod=min(lod),
              max_lod=max(lod))
}
