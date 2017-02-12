#' Get SVS information in window around peak
#'
#' Get SVS information from RDS database within \code{window_Mbp} of \code{peak_Mbp} on \code{chri_id}
#'
#' @param chr_id chromosome identifier
#' @param peak_Mbp position in Mbp of peak
#' @param window_Mbp half-width of \code{window} around \code{peak_Mbp}
#' @param datapath path to Derived Data
#'
#' @return table of class \code{snpinfo} of top_snps near maximum
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @examples
#' \dontrun{get_snpinfo(chr_id, peak_Mbp, window_Mbp)}
#'
#' @export
#' @importFrom dplyr filter mutate select tbl_df
get_svs8 <- function(chr_id, peak_Mbp, window_Mbp, datapath) {
  start_val <- convert_bp(peak_Mbp - window_Mbp)
  end_val <- convert_bp(peak_Mbp + window_Mbp)

  svs8_len <- dplyr::filter(
    dplyr::tbl_df(
      readRDS(file.path(datapath, "svs8_len.rds"))),
    chrom == chr_id,
    start >= start_val,
    end <= end_val)

  ## Extract CC Founder names
  cc_founders <- names(svs8_len)[-(1:4)]

  svs8_len <-
    dplyr::mutate(svs8_len,
                  snp_id = paste0(chrom, ":", start, "-", end, ":", type),
                  chr = chrom,
                  pos_Mbp = convert_bp(start, FALSE))
  svs8_len$alleles <- paste(svs8_len$type,
                            apply(svs8_len[,cc_founders], 1,
                              function(x) paste(unique(x[x>0]), collapse=",")),
                            sep = ":")
  tmp2 <- 2^(0:7)
  svs8_len$sdp <- apply(svs8_len[,cc_founders], 1,
        function(x) sum(tmp2*(x>0)))
  svs8_len <- dplyr::select(
    dplyr::rename(svs8_len, svs_type = type),
    snp_id, chr, pos_Mbp, alleles, sdp, svs_type)
  class(svs8_len) <- c("snpinfo", class(svs8_len))
  svs8_len
}
