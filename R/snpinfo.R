#' Get SNP and InDel information in window around peak
#'
#' Get SNP and InDel information from SQLite databases within \code{window_Mbp} of \code{peak_Mbp} on \code{chri_id}
#'
#' @param chr_id chromosome identifier
#' @param peak_Mbp position in Mbp of peak
#' @param window_Mbp half-width of \code{window} around \code{peak_Mbp}
#' @param datapath path to Derived Data
#' @param info_type type of info: one of c("snps","indels")
#'
#' @return table of class \code{snpinfo} of top_snps near maximum lod
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @examples
#' \dontrun{get_snpinfo(chr_id, peak_Mbp, window_Mbp)}
#'
#' @export
#' @importFrom dplyr collect filter src_sqlite tbl
get_snpinfo <- function(chr_id, peak_Mbp, window_Mbp,
                        datapath, info_type = c("snps","indels")) {
  start_val <- peak_Mbp - window_Mbp
  end_val <- peak_Mbp + window_Mbp
  ## identify table in SQLite, and filter to chr and bp interval.
  info_type <- match.arg(info_type)
  file.sqlite <- paste0("ccfounder", info_type, ".sqlite")
  my_db <- dplyr::src_sqlite(file.path(datapath, file.sqlite))
  snpinfo <- dplyr::filter(
    dplyr::tbl(my_db, info_type),
    chr == chr_id,
    pos_Mbp >= start_val,
    pos_Mbp <= end_val)
  ## Collect from database.
  snpinfo <- dplyr::collect(snpinfo, n=1e+6)
  class(snpinfo) <- c("snpinfo", class(snpinfo))
  snpinfo
}
#' Summary of snpinfo object
#'
#' @param object of class \code{near_snps}
#' @param ... other arguments not used
#'
#' @return table summary
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @examples
#' \dontrun{summary(object)}
#'
#' @method summary snpinfo
#' @rdname snpinfo
#' @export
#' @importFrom dplyr arrange desc group_by mutate summarize ungroup
summary.snpinfo <- function(object, ...) {
  dplyr::arrange(
    dplyr::mutate(
      dplyr::ungroup(
        dplyr::summarize(
          dplyr::group_by(object, sdp),
          pct = round(100 * n() / nrow(snpinfo), 2),
          min_Mbp = min(pos_Mbp),
          max_Mbp = max(pos_Mbp))),
      pattern = sdp_to_pattern(sdp)),
    dplyr::desc(pct))
}

#' Get SNP genotype probabilities in window around peak
#'
#' Get SNP information from SQLite database within \code{window_Mbp} of \code{peak_Mbp} on \code{chri_id}
#'
#' @param chr_id chromosome identifier
#' @param peak_Mbp position in Mbp of peak
#' @param window_Mbp half-width of \code{window} around \code{peak_Mbp}
#' @param phename names of phenotypes
#' @param probs_obj object of class \code{\link[qtl2geno]{calc_genoprob}} for \code{chr_id}
#' @param datapath path to Derived Data
#'
#' @return table of class \code{near_snps} of top_snps near maximum lod
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @examples
#' \dontrun{get_snpprobs(chr_id, peak_Mbp, window_Mbp, scan_obj, probs_obj, datapath)}
#'
#' @export
#' @importFrom dplyr bind_rows mutate
#' @importFrom qtl2scan genoprob_to_snpprob
get_snpprobs <- function(chr_id=NULL, peak_Mbp=NULL, window_Mbp=NULL,
                         phename, probs_obj, datapath) {
  if(is.null(chr_id) | is.null(peak_Mbp) | is.null(window_Mbp))
    return(NULL)

  if(window_Mbp == 0) {
    window_Mbp <- 3
    cat(file=stderr(),
        "\nNo window_Mbp provided -- set to 3\n")
  }
  if(peak_Mbp == 0) {
    cat(file=stderr(),
        "\nNo peak_Mbp provided -- set to midpoint\n")
    peak_Mbp <- mean(range(probs_obj$map[[1]]))
  }
  snpinfo <- dplyr::mutate(
    get_snpinfo(chr_id, peak_Mbp, window_Mbp, datapath),
    svs_type = "SNP")
  indelinfo <- dplyr::select(
    dplyr::mutate(
      get_snpinfo(chr_id, peak_Mbp, window_Mbp, datapath, info_type = "indels"),
      svs_type = "Indel"),
    -allele)
  svsinfo <- get_svs8(chr_id, peak_Mbp, window_Mbp, datapath)
  snpinfo <- dplyr::bind_rows(snpinfo, indelinfo, svsinfo)
  ## Need names pos and snp for genoprob_to_snpprob.
  snpinfo <- dplyr::mutate(snpinfo,
                           pos = pos_Mbp,
                           snp = snp_id,
                           svs_type = factor(svs_type))
  qtl2scan::genoprob_to_snpprob(probs_obj, as.data.frame(snpinfo))
}
