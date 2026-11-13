#!/usr/bin/env Rscript
#
# attach_bcr_to_seurat.R - add per-cell BCR columns to an existing GEX Seurat object.
#
# The R-side twin of SCRATCH-TCRmerge, for cohorts small enough to live in a Seurat
# object (R's dgCMatrix caps at ~2.1e9 non-zeros; for anything bigger use TCRmerge on the
# .h5ad, receptor = BCR). BCR data is per-cell METADATA, so this adds bcr_* columns to
# @meta.data; no matrix is touched and no cell is added or removed.
#
# Barcode reconciliation follows TCRmerge: every candidate key is scored against the
# Seurat cell names and the best match rate wins; library-type words (GEX, BCR, VDJ ...)
# are dropped from the sample part of both sides as extra candidates; a sample is never
# dropped from the key, because 10x barcodes collide across libraries.
#
# If the object already carries TCR columns (tcr_has_tcr_match), cells matched by both
# get vdj_dual_receptor = TRUE: in one GEM almost always a T+B doublet. Flagged only.
#
# Usage: attach_bcr_to_seurat.R --seurat gex.rds --cells post_qc_cells.tsv
#                               [--prefix bcr_] [--library-tokens GEX,TCR,VDJ,BCR,ADT]
#                               [--min-match-rate 0.01] [--out bcr_gex_seurat.rds]

suppressPackageStartupMessages({
    library(Seurat)
    library(data.table)
})

.args <- commandArgs(trailingOnly = TRUE)
.get <- function(flag, default = NA_character_) {
    i <- match(flag, .args)
    if (!is.na(i) && length(.args) > i) .args[i + 1] else default
}
seu_path <- .get("--seurat"); cells_path <- .get("--cells")
prefix   <- .get("--prefix", "bcr_")
tokens   <- .get("--library-tokens", "GEX,TCR,VDJ,BCR,ADT")
min_rate <- as.numeric(.get("--min-match-rate", "0.01"))
out      <- .get("--out", "bcr_gex_seurat.rds")
msg <- function(...) message("[attach_bcr] ", ...)

seu   <- readRDS(seu_path)
cells <- as.data.frame(fread(cells_path, sep = "\t", na.strings = ""))
obs   <- colnames(seu)
msg(sprintf("Seurat: %s cells; BCR table: %s cells", format(length(obs), big.mark = ","),
            format(nrow(cells), big.mark = ",")))

tok <- toupper(trimws(strsplit(tokens, ",")[[1]])); tok <- tok[nzchar(tok)]
norm_sample <- function(s) vapply(strsplit(as.character(s), "[-_]"), function(p) {
    k <- p[!toupper(p) %in% tok]; if (length(k)) paste(k, collapse = "-") else paste(p, collapse = "-")
}, character(1))
# <sample>_<barcode...> -> <normalised sample>_<barcode...>; names without "_" are left alone.
norm_name <- function(n) {
    has <- grepl("_", n)
    s <- sub("_[^_]*$", "", n); t <- sub("^.*_", "", n)
    ifelse(has, paste0(norm_sample(s), "_", t), n)
}

bc <- cells$barcode; sm <- cells$sample
cand <- list(
    "cell_id"          = cells$cell_id,
    "cell_id-1"        = paste0(cells$cell_id, "-1"),
    "sample_barcode-1" = paste0(sm, "_", bc, "-1"),
    "sample_barcode"   = paste0(sm, "_", bc),
    "sample-barcode"   = paste0(sm, "-", bc)
)
obs_norm <- character(0)
if (length(tok)) {
    nn <- norm_name(obs)
    dup <- nn %in% nn[duplicated(nn)]                     # ambiguous after normalising: excluded
    obs_norm <- setNames(obs[!dup], nn[!dup])
    cand[["libnorm:sample_barcode-1"]] <- paste0(norm_sample(sm), "_", bc, "-1")
    cand[["libnorm:sample_barcode"]]   <- paste0(norm_sample(sm), "_", bc)
}
resolve <- function(name, keys) if (startsWith(name, "libnorm:")) unname(obs_norm[keys]) else ifelse(keys %in% obs, keys, NA)
rep <- do.call(rbind, lapply(names(cand), function(n) {
    hit <- resolve(n, cand[[n]])
    data.frame(candidate = n, matched = sum(!is.na(hit)), total = length(hit),
               exact = !startsWith(n, "libnorm:"))
}))
rep$match_rate <- round(rep$matched / pmax(rep$total, 1), 6)
rep <- rep[order(-rep$matched, -rep$exact), ]
fwrite(rep, "bcr_barcode_join_report.tsv", sep = "\t")
print(rep, row.names = FALSE)

best <- rep$candidate[1]
if (rep$matched[1] == 0 || rep$match_rate[1] < min_rate)
    stop(sprintf("best candidate '%s' matched %d/%d BCR cells (< %.2f%%). Seurat cells look like %s; BCR cell_ids like %s",
                 best, rep$matched[1], rep$total[1], 100 * min_rate,
                 paste(head(obs, 2), collapse = ", "), paste(head(cells$cell_id, 2), collapse = ", ")))
cells$.obs <- resolve(best, cand[[best]])
cells <- cells[!is.na(cells$.obs) & !duplicated(cells$.obs), ]
msg(sprintf("using '%s': %s BCR cells matched a GEX cell", best, format(nrow(cells), big.mark = ",")))

md <- cells; rownames(md) <- md$.obs; md$.obs <- NULL
names(md) <- paste0(prefix, names(md))
clash <- intersect(names(md), colnames(seu@meta.data))
if (length(clash)) stop(length(clash), " '", prefix, "' columns already in @meta.data (e.g. ",
                        paste(head(clash, 3), collapse = ", "), ") - object already has BCR data")
seu <- AddMetaData(seu, md)
flag <- paste0(prefix, "has_bcr_match")
seu[[flag]] <- colnames(seu) %in% rownames(md)

tcr_flags <- grep("has_tcr_match$", colnames(seu@meta.data), value = TRUE)
if (length(tcr_flags)) {
    tcr <- Reduce(`|`, lapply(tcr_flags, function(c) as.logical(seu@meta.data[[c]]) %in% TRUE))
    seu$vdj_dual_receptor <- tcr & seu@meta.data[[flag]]
    msg(sprintf("%d cells carry both TCR and BCR -> vdj_dual_receptor (likely doublets; flagged only)",
                sum(seu$vdj_dual_receptor)))
}
msg(sprintf("%s of %s GEX cells carry BCR data (%.2f%%)", format(sum(seu@meta.data[[flag]]), big.mark = ","),
            format(ncol(seu), big.mark = ","), 100 * mean(seu@meta.data[[flag]])))
saveRDS(seu, out)
msg("wrote ", out)
