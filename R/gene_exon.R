#' Get exons for set of genes
#'
#' Match up exon start,stop,strand with genes.
#'
#' @param top_snps_tbl tbl from \code{\link[qtl2pattern]{top_snps_all}} or \code{\link[qtl2scan]{top_snps}}
#' @param sql_filename path to \code{\link{get_mgi_features}}
#'
#' @return tbl of exon and gene features
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @examples
#' \dontrun{get_gene_exon_snp(top_snps_tbl, sql_filename)}
#'
#' @export
#' @rdname gene_exon
#' @importFrom dplyr arrange desc distinct mutate select
get_gene_exon_snp <- function(top_snps_tbl, sql_filename) {
  ## Only need distinct snp_id.
  top_snps_tbl <- dplyr::arrange(
    dplyr::select(
      dplyr::distinct(top_snps_tbl, snp_id, .keep_all=TRUE),
      -pheno),
    pos_Mbp)

  if(is.null(top_snps_tbl))
    return(NULL)
  if(!nrow(top_snps_tbl))
    return(NULL)

  chr_id <- as.character(unique(top_snps_tbl$chr))
  if(length(chr_id) != 1)
    stop("need exactly 1 chromosome in top_snps_tbl")
  range_Mbp <- range(top_snps_tbl$pos_Mbp)
  feature_tbl <- get_mgi_features(chr_id,
                                  range_Mbp[1], range_Mbp[2],
                                  with_name = FALSE,
                                  sql_file = sql_filename)
  gene_snp <- get_gene_snp(
    dplyr::select(
      dplyr::mutate(top_snps_tbl, pos = pos_Mbp), 
      snp_id,pos,lod),
    feature_tbl)
  get_gene_exon(feature_tbl, gene_snp)
}

#' Get exons for set of genes
#'
#' Match up exon start,stop,strand with genes.
#'
#' @param feature_tbl tbl of features from \code{\link{get_mgi_feature}}
#' @param gene_snp tbl of genes with SNPs IDs from \code{\link{match_feature_snp}}
#'
#' @return tbl of exon and gene features
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @examples
#' \dontrun{get_gene_exon(feature_snp)}
#'
#' @export
#' @rdname gene_exon
#' @importFrom dplyr bind_rows distinct filter
get_gene_exon <- function(feature_tbl, gene_snp) {
  exons <- list()
  ## Need to get unique genes -- duplication with SNPs.
  gene_snp <- dplyr::distinct(gene_snp, gene, .keep_all=TRUE)

  ## Use gene name as Name from feature_tbl.
  for(exoni in seq_len(nrow(gene_snp))) {
    genei <- gene_snp$gene[exoni]
    ## get genei and exons spanning genei
    exons[[genei]] <- dplyr::filter(
      dplyr::filter(feature_tbl,
                    type %in% c("exon","gene"),
                    start >= gene_snp$start[exoni],
                    stop <= gene_snp$stop[exoni]),
      (!is.na(Name) & Name==genei) | type=="exon")
    strandi <- gene_snp$strand[exoni]
    if(strandi != ".")
      exons[[genei]] <- dplyr::filter(exons[[genei]], strand==strandi)
  }
  out <- dplyr::distinct(
    dplyr::bind_rows(exons, .id="gene"),
    start, stop, strand, .keep_all=TRUE)
  class(out) <- c("gene_exon", class(out))

  out
}

#' Plot of exons for a gene with SNPs
#'
#' Uses \code{\link{gene_plot}} to plot genes, exons, mRNA with SNPs.
#'
#' @param exon_tbl tbl of feature information from \code{\link{get_mgi_features}}
#' @param top_snps_tbl tbl of SNP data from \code{\link[qtl2scan]{top_snps}}
#' @param plot_now plot now if TRUE
#' @param ... arguments passed along to \code{\link{gene_plot}}
#'
#' @return list of ggplots (see \code{\link{gene_plot}})
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords hplot
#'
#' @method plot gene_exon
#' @rdname gene_exon
#' @export
#' @importFrom dplyr group_by summarize ungroup
#' @importFrom ggplot2 ggtitle
plot.gene_exon <- function(gene_exon, top_snps_tbl=NULL, plot_now=TRUE,
                           genes = unique(gene_exon$gene), ...) {
  p <- list()
  ## Reduce to max lod per SNP
  if(!is.null(top_snps_tbl)) {
    top_snps_tbl <- dplyr::ungroup(
      dplyr::summarize(
        dplyr::group_by(top_snps_tbl, snp_id,pos_Mbp),
        lod=max(lod)))
  }
  for(genei in genes) {
    p[[genei]] <-
      plot.feature_tbl(
        dplyr::mutate(
          dplyr::filter(gene_exon, gene==genei),
          gene = NA),
        top_snps_tbl = top_snps_tbl,
        str_rect="",
        ...) +
      ggplot2::ggtitle(genei)
  }
  if(plot_now & length(p)) {
    for(genei in genes)
      print(p[[genei]])
  }
  invisible(p)
}

#' Summary of exons for a gene with SNPs
#'
#' Returns table of gene and its exons.
#'
#' @param gene_exon tbl of feature information from \code{\link{get_gene_exon}}
#' @param gene_name name of gene as character string
#' @param top_snps_tbl table of top SNPs in region from \code{\link[qtl2scan]{top_snps}}
#' @param extra_bp extra region beyond gene for SNPs (in bp)
#'
#' @return tbl of summary
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @method summary gene_exon
#' @rdname gene_exon
#' @export
#' @importFrom dplyr arrange desc distinct filter group_by mutate select summarize ungroup
summary.gene_exon <- function(gene_exon, gene_name=NULL,
                              top_snps_tbl = NULL,
                              extra_bp = 5000) {
  ## Want to add columns for each phenotype
  ## with number of SNPs within extra_bp of each gene.
  ## How to do this cleverly?
  if(!is.null(gene_name)) {
    out <- dplyr::distinct(
      dplyr::mutate(
        dplyr::arrange(
          dplyr::select(
            dplyr::filter(gene_exon, gene==gene_name),
            gene, source, type, start, stop, strand),
          dplyr::desc(type)),
        start = convert_bp(start, FALSE),
        stop = convert_bp(stop, FALSE)),
      start, stop, strand, .keep_all=TRUE)
  } else {
    out <- dplyr::ungroup(
      dplyr::summarize(
        dplyr::group_by(
          dplyr::distinct(
            dplyr::filter(gene_exon, type != "gene"),
            start, stop, strand, .keep_all=TRUE),
          gene),
        exons = n(),
        min.len = min(stop-start),
        max.len = max(stop-start),
        sum.len = sum(stop-start),
        min_Mbp = convert_bp(min(start), FALSE),
        max_Mbp = convert_bp(max(stop), FALSE),
        strand = strand[1]))
    if(!is.null(top_snps_tbl)) {
      ## Goal: add columns to out for each pheno in top_snps_tbl.
      ## Column should have number of SNPs within extra_bp of gene.
      top_snps_tbl <- dplyr::select(top_snps_tbl, pheno, pos_Mbp, lod)
      pheno_names <- sort(unique(top_snps_tbl$phen))
      outlim <- out[,c("min_Mbp","max_Mbp")]
      extra_bp <- convert_bp(extra_bp, FALSE)
      outlim[,1] <- outlim[,1] - extra_bp
      outlim[,2] <- outlim[,2] + extra_bp
      for(pheno_val in pheno_names) {
        out[[pheno_val]] <-
          apply(outlim, 1,
                function(x,y) {
                  in_region <- y$pos_Mbp >= x[1] & y$pos_Mbp <= x[2]
                  if(any(in_region))
                    max(y$lod[in_region])
                  else
                    NA
                },
                dplyr::filter(top_snps_tbl, pheno == pheno_val))
      }
    }
  }
  out
}
#' @rdname gene_exon
#' @export
subset.gene_exon <- function(x, gene_val, ...) {
  x <- dplyr::filter(x, gene %in% gene_val)
  class(x) <- c("gene_exon", class(x))
  x
}
#' @rdname gene_exon
#' @export
#' @importFrom ggplot2 ggtitle
plot_gene_exon <- function(gene_exon_tbl, top_snps_tbl, gene_name, pheno) {
  if(nrow(top_snps_tbl)) {
    p <- plot(gene_exon_tbl, top_snps_tbl,
              FALSE, genes = gene_name)
    p[[1]] + ggplot2::ggtitle(paste(gene_name, "SNPs for", pheno))
  } else
    plot_null()
}
