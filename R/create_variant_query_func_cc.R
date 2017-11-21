#' @export
create_variant_query_func_cc <- function(dbpath) {
  function(chr, start, end) {
    peak_Mbp <- (start + end) / 2
    window_Mbp <- peak_Mbp - start
    snpinfo <- 
      dplyr::mutate(
        CCSanger::get_snpinfo(
          chr, peak_Mbp, window_Mbp, dbpath),
        svs_type = "snp")
    indelinfo <- 
      dplyr::select(
        dplyr::mutate(
          CCSanger::get_snpinfo(
            chr, peak_Mbp, window_Mbp, dbpath, 
            info_type = "indels"),
          svs_type = "indel"),
        -allele)
    svsinfo <- 
      dplyr::mutate(
        CCSanger::get_svs8(
          chr, peak_Mbp, window_Mbp, dbpath),
        chr = as.character(chr),
        snp_id = ifelse("rs" == stringr::str_sub(snp_id, 1, 2),
                        snp_id, 
                        paste0(
                          "SV_", 
                          stringr::str_replace_all(
                            stringr::str_replace(
                              snp_id, ":[A-Z].*",""),
                            ":|-", "_"))))
    # "SV_1_34002149_34002151" "1:34002149-34002151:INS"
    dplyr::arrange(
      dplyr::rename(
        dplyr::bind_rows(
          snpinfo, indelinfo, svsinfo),
        pos = pos_Mbp,
        type = svs_type,
        consequence = csq),
      pos)
  }
}
