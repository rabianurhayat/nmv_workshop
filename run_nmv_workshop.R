#!/usr/bin/env Rscript
# Minimal ONT methylation demo (NanoMethViz + dmrseq), mm10 chr7 subset.
# Runs end-to-end: download -> build NMR -> light DMR around Peg3/Kcnq1ot1/Cdkn1c -> save figs & PPTX.

set.seed(1234)

# ---- helpers: install/load ----
pkg_ok <- function(p) suppressWarnings(require(p, character.only=TRUE, quietly=TRUE))
inst <- function(pkgs){
  if (!pkg_ok("BiocManager")) install.packages("BiocManager")
  bioc <- c("NanoMethViz","dmrseq","plyranges","bsseq","BiocParallel","GenomicRanges","AnnotationHub")
  for (p in pkgs){
    if (!pkg_ok(p)){
      if (p %in% bioc) BiocManager::install(p, ask=FALSE, update=FALSE) else install.packages(p)
      library(p, character.only=TRUE)
    }
  }
}
inst(c("tidyverse","NanoMethViz","dmrseq","plyranges","bsseq","BiocParallel","GenomicRanges","officer"))

# ---- folders ----
dir.create("input", showWarnings = FALSE)
dir.create("data",  showWarnings = FALSE)
dir.create("figs",  showWarnings = FALSE)

# ---- data: zenodo modBAMs if missing ----
if (!length(list.files("input", pattern="bam$"))){
  message("Downloading demo input (Zenodo)…")
  options(timeout=600)
  url <- "https://zenodo.org/records/12747551/files/input.tar.gz?download=1"
  download.file(url, "input.tar.gz", mode="wb", method="libcurl")
  utils::untar("input.tar.gz", exdir = ".")
  file.remove("input.tar.gz")
}
bam_files <- dir("input", pattern="bam$", full.names=TRUE)
samples   <- readr::read_tsv("input/sample_anno.tsv", show_col_types = FALSE)

# ---- exon anno (mm10 chr7) ----
exon_anno <- NanoMethViz::get_exons_mm10() |> dplyr::filter(chr == "chr7")

# ---- build NMR (via modbam_to_tabix once) ----
if (!file.exists("data/methy.tsv.bgz")) {
  message("Converting modBAM -> tabix (one-time)…")
  mbr <- NanoMethViz::ModBamResult(
    methy   = NanoMethViz::ModBamFiles(paths = bam_files, samples = samples$sample),
    samples = samples, exons = exon_anno, mod_code = "m"
  )
  NanoMethViz::modbam_to_tabix(mbr, "data/methy.tsv.bgz")
}
nmr <- NanoMethViz::NanoMethResult(methy="data/methy.tsv.bgz", samples=samples, exons=exon_anno)

# ---- quick figure helper ----
save_gene_plot <- function(obj, gene, regions=NULL, w=10, h=6, dpi=300){
  p <- NanoMethViz::plot_gene(obj, gene, anno_regions = regions)
  ggplot2::ggsave(file.path("figs", paste0(gene, ".png")), p, width=w, height=h, dpi=dpi)
  ggplot2::ggsave(file.path("figs", paste0(gene, ".pdf")), p, width=w, height=h)
  invisible(p)
}

# ---- baseline peg3 plot ----
save_gene_plot(nmr, "Peg3")

# ---- BSseq (cached) ----
if (!file.exists("data/bss.rds")) {
  message("Converting NMR -> BSseq (one-time)…")
  bss <- NanoMethViz::methy_to_bsseq(nmr)
  saveRDS(bss, "data/bss.rds")
} else bss <- readRDS("data/bss.rds")

# keep chr7 only + soft coverage filter
bss <- bss[as.character(GenomeInfoDb::seqnames(rowRanges(bss)))=="chr7", ]
cov <- bsseq::getCoverage(bss)
pat <- grepl("pat", colnames(cov)); mat <- grepl("mat", colnames(cov))
keep <- rowSums(cov[,pat] > 0) >= 2 & rowSums(cov[,mat] > 0) >= 2
bss <- bss[keep, ]
pData(bss)$condition <- NanoMethViz::samples(nmr)$group

# ---- target promoters (TSS ±10 kb) ----
gene_anno_gr <- plyranges::as_granges(dplyr::rename(NanoMethViz::exons_to_genes(NanoMethViz::exons(nmr)), seqnames="chr"))
targets <- c("Peg3","Kcnq1ot1","Cdkn1c")
proms <- gene_anno_gr[gene_anno_gr$symbol %in% targets] |>
         plyranges::anchor_5p() |> plyranges::mutate(width=1) |> plyranges::stretch(10000)

keep_rows <- GenomicRanges::countOverlaps(rowRanges(bss), proms, ignore.strand=TRUE) > 0
bss_small <- bss[keep_rows, , drop=FALSE]

# ---- dmrseq (serial, minimal permutations for compatibility) ----
BiocParallel::register(BiocParallel::SerialParam())

perm_name <- intersect(c("maxPerms","nPerms","B"), names(formals(dmrseq::dmrseq)))[1]
dmr_args <- list(
  bss_small, testCovariate="condition",
  cutoff=0.03, minNumRegion=5,
  BPPARAM=BiocParallel::SerialParam()
)
if (!is.na(perm_name)) dmr_args[[perm_name]] <- 1L
regions <- do.call(dmrseq::dmrseq, dmr_args)
saveRDS(regions, "data/regions_promoters.rds")

# ---- annotate, save tables, make plots ----
ov  <- NanoMethViz::join_overlap_intersect(regions, proms)
sig <- tibble::as_tibble(ov) |> dplyr::rename(chr="seqnames") |> dplyr::arrange(qval)
readr::write_tsv(sig, "figs/DMRs_TSS10k.tsv")

options("NanoMethViz.highlight_col"="green")
save_gene_plot(nmr, "Peg3",     sig)
save_gene_plot(nmr, "Kcnq1ot1", sig)
save_gene_plot(nmr, "Cdkn1c",   sig)

# ---- PPTX with slides ----
library(officer); library(magrittr)
doc <- read_pptx()
for (g in targets){
  f <- file.path("figs", paste0(g, ".png"))
  if (file.exists(f)) {
    doc <- doc %>%
      add_slide(layout="Title and Content", master="Office Theme") %>%
      ph_with(value=paste("Methylation:", g), location=ph_location_type(type="title")) %>%
      ph_with(value=external_img(f, width=9, height=5), location=ph_location_type(type="body"))
  }
}
print(doc, target="figs/methylation_plots.pptx")

# ---- session info ----
utils::capture.output(sessionInfo(), file = "figs/sessionInfo.txt")
message("Done. See ./figs for outputs.")
