#!/usr/bin/env Rscript
#
# merge_bcr_object.R - per-cell BCR objects at the pre- and post-QC stage.
#
# The BCR counterpart of TCRtoolkit's bin/merge_vdj_object.R, writing the same bundle so
# SCRATCH-TCRmerge (receptor = BCR) and R users treat both receptors alike:
#
#   <stage>_cells.tsv         one row per cell. THE joinable table: TCRmerge adds it to
#                             the GEX .h5ad's .obs under a bcr_ prefix.
#   <stage>_summary.tsv       one row per sample.
#   <stage>_combineBCR.rds    scRepertoire::combineBCR() list - hand it to
#                             combineExpression() if you work in Seurat directly.
#   <stage>_seurat.rds        Seurat object, one cell per barcode, BCR fields in
#                             @meta.data, backed by a 2-row all-zero PLACEHOLDER assay:
#                             VDJ libraries have no genes. Do not normalise/cluster it.
#                             Its purpose is the metadata; see attach_bcr_to_seurat.R or
#                             SCRATCH-TCRmerge for joining onto real expression.
#
# pre_qc  = every IG contig Cell Ranger called (Cell Ranger annotations)
# post_qc = cells passing BCR QC, with IMGT re-annotation, SCOPer clones, and - when
#           those steps ran - SHM and lineage-tree membership
#
# Cell key: cell_id = <sample>_<barcode without -1>, exactly TCRtoolkit's convention, so
# TCRmerge's barcode reconciliation and library-word normalisation (…-GEX vs …-BCR)
# apply unchanged.
#
# Usage:
#   merge_bcr_object.R --contigs-pre all_contigs_before_qc.tsv
#                      --contigs-post all_contigs_after_qc.tsv
#                      --airr a_germ-pass.tsv,b_germ-pass.tsv
#                      [--lineage a_lineage_summary.tsv,...] [--unmutated-cutoff 0.01]
#                      [--samplesheet samplesheet.valid.csv]
#
# Sample metadata: every samplesheet column that is not a file path or a pipeline key
# (timepoint, tissue, batch, condition, patient_id, ...) is copied onto each cell of
# that sample, so the BCR tables can be split by timepoint on their own, before any
# GEX merge. Columns the cell tables already have are never overwritten.

suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
})

.args <- commandArgs(trailingOnly = TRUE)
.get <- function(flag, default = NA_character_) {
    i <- match(flag, .args)
    if (!is.na(i) && length(.args) > i) .args[i + 1] else default
}
split_list <- function(x) {
    if (is.na(x) || !nzchar(x)) return(character(0))
    p <- trimws(strsplit(x, ",", fixed = TRUE)[[1]])
    p[nzchar(p) & file.exists(p) & file.size(p) > 0]
}
opt <- list(
    pre      = .get("--contigs-pre"),
    post     = .get("--contigs-post"),
    airr     = split_list(.get("--airr", "")),
    lineage  = split_list(.get("--lineage", "")),
    unmut    = as.numeric(.get("--unmutated-cutoff", "0.01")),
    sheet    = .get("--samplesheet", ""),
    project  = .get("--project", "bcr")
)
msg <- function(...) message("[merge_bcr] ", ...)

read_tsv <- function(p) as.data.frame(fread(p, sep = "\t", colClasses = "character", na.strings = c("", "NA")))
num <- function(x) suppressWarnings(as.numeric(x))
gene <- function(x) {                                   # "IGHV3-23*01,IGHV3-23D*01" -> "IGHV3-23"
    x <- sub(",.*$", "", as.character(x)); sub("\\*.*$", "", x)
}
isotype_class <- function(c_call) {
    g <- toupper(gene(c_call))
    dplyr::case_when(
        is.na(g) | g == ""      ~ NA_character_,
        grepl("^IGHM", g)       ~ "IgM",
        grepl("^IGHD", g)       ~ "IgD",
        grepl("^IGHG", g)       ~ "IgG",
        grepl("^IGHA", g)       ~ "IgA",
        grepl("^IGHE", g)       ~ "IgE",
        TRUE                    ~ NA_character_)
}
strip_junction <- function(j) ifelse(is.na(j) | nchar(j) < 3, NA_character_, substr(j, 2, nchar(j) - 1))
size_category <- function(n) cut(n, breaks = c(0, 1, 5, 20, 100, Inf),
                                 labels = c("Single", "Small", "Medium", "Large", "Hyperexpanded"))
top_by_umi <- function(df, umi_col) df[order(-num(df[[umi_col]])), , drop = FALSE] %>%
    distinct(cell_id, .keep_all = TRUE)

# ── PRE-QC cells (Cell Ranger calls) ──────────────────────────────────────────
pre <- read_tsv(opt$pre)
msg(sprintf("pre-QC: %s IG contigs, %s cells", format(nrow(pre), big.mark = ","),
            format(n_distinct(pre$cell_id), big.mark = ",")))
lt <- function(df) df$chain %in% c("IGK", "IGL")
h <- top_by_umi(pre[pre$chain == "IGH", ], "umis")
l <- top_by_umi(pre[lt(pre), ], "umis")
pre_cells <- pre %>% group_by(cell_id) %>%
    summarise(sample = first(sample), subject = first(subject),
              barcode = sub("-\\d+$", "", first(barcode)),
              n_contigs = n(), n_heavy = sum(chain == "IGH"), n_light = sum(chain %in% c("IGK", "IGL")),
              passed_bcr_qc = any(tolower(passed_bcr_qc) == "true"),
              .groups = "drop") %>%
    left_join(transmute(h, cell_id, cdr3_heavy = cdr3, cdr3_heavy_nt = cdr3_nt, ighv = gene(v_gene),
                        ighd = gene(d_gene), ighj = gene(j_gene), isotype = gene(c_gene),
                        umis_heavy = num(umis)), by = "cell_id") %>%
    left_join(transmute(l, cell_id, light_locus = chain, cdr3_light = cdr3, cdr3_light_nt = cdr3_nt,
                        iglv = gene(v_gene), iglj = gene(j_gene), umis_light = num(umis)), by = "cell_id") %>%
    mutate(isotype_class = isotype_class(isotype),
           multi_heavy = n_heavy > 1, multi_light = n_light > 1,
           has_bcr = TRUE, paired_bcr = !is.na(cdr3_heavy) & !is.na(cdr3_light),
           CTaa = paste0(ifelse(is.na(cdr3_heavy), "NA", paste0("H:", cdr3_heavy)), "|",
                         ifelse(is.na(cdr3_light), "NA", paste0("L:", cdr3_light)))) %>%
    as.data.frame()

# ── POST-QC cells (IMGT + clones + SHM + lineage) ────────────────────────────
airr <- bind_rows(lapply(opt$airr, read_tsv))
if (!nrow(airr)) stop("no post-QC AIRR rows - did every sample fail BCR QC?")
umi <- if ("umi_count" %in% names(airr)) "umi_count" else "consensus_count"
ah <- top_by_umi(airr[airr$locus == "IGH", ], umi)
al <- top_by_umi(airr[airr$locus %in% c("IGK", "IGL"), ], umi)
has <- function(col) col %in% names(airr)
pick <- function(df, col) if (col %in% names(df)) df[[col]] else rep(NA_character_, nrow(df))

post_cells <- airr %>% group_by(cell_id) %>%
    summarise(sample = first(sample), subject = first(subject),
              barcode = sub("-\\d+$", "", first(barcode)),
              n_contigs = n(), clone_id = first(na.omit(clone_id)),
              clone_subgroup = if ("clone_subgroup" %in% names(airr)) first(na.omit(clone_subgroup)) else NA,
              .groups = "drop") %>%
    left_join(data.frame(cell_id = ah$cell_id,
                         cdr3_heavy = strip_junction(ah$junction_aa), junction_heavy_nt = ah$junction,
                         ighv = gene(ah$v_call), ighd = gene(pick(ah, "d_call")), ighj = gene(ah$j_call),
                         ighv_allele = sub(",.*$", "", ah$v_call),
                         isotype = gene(pick(ah, "c_call")), umis_heavy = num(ah[[umi]]),
                         mu_freq_heavy = num(pick(ah, "mu_freq")), mu_count_heavy = num(pick(ah, "mu_count")),
                         stringsAsFactors = FALSE), by = "cell_id") %>%
    left_join(data.frame(cell_id = al$cell_id, light_locus = al$locus,
                         cdr3_light = strip_junction(al$junction_aa), junction_light_nt = al$junction,
                         iglv = gene(al$v_call), iglj = gene(al$j_call), umis_light = num(al[[umi]]),
                         mu_freq_light = num(pick(al, "mu_freq")), stringsAsFactors = FALSE),
              by = "cell_id") %>%
    mutate(isotype_class = isotype_class(isotype),
           has_bcr = TRUE, paired_bcr = !is.na(cdr3_heavy) & !is.na(cdr3_light),
           CTaa = paste0(ifelse(is.na(cdr3_heavy), "NA", paste0("H:", cdr3_heavy)), "|",
                         ifelse(is.na(cdr3_light), "NA", paste0("L:", cdr3_light))),
           CTgene = paste0(ighv, ".", ighj, ".", isotype, "|", iglv, ".", iglj)) %>%
    as.data.frame()

# Pre-QC chain-multiplicity flags: QC keeps one heavy and one light per cell, so post-QC
# counts would read FALSE by construction. Carry the before-QC truth across.
post_cells <- post_cells %>%
    left_join(transmute(pre_cells, cell_id, multi_heavy_pre_qc = multi_heavy,
                        multi_light_pre_qc = multi_light), by = "cell_id")

# Clone sizes: cells per clone in the subject, and in this sample; samples the clone spans.
cl <- post_cells %>% filter(!is.na(clone_id))
by_subject <- cl %>% count(clone_id, name = "clone_size")
by_sample  <- cl %>% count(clone_id, sample, name = "clone_size_sample")
n_samples  <- cl %>% distinct(clone_id, sample) %>% count(clone_id, name = "clone_n_samples")
post_cells <- post_cells %>% left_join(by_subject, by = "clone_id") %>%
    left_join(by_sample, by = c("clone_id", "sample")) %>% left_join(n_samples, by = "clone_id") %>%
    mutate(clone_size_category = as.character(size_category(clone_size)))

if (has("mu_freq")) {
    post_cells$shm_status <- ifelse(is.na(post_cells$mu_freq_heavy), NA_character_,
                                    ifelse(post_cells$mu_freq_heavy < opt$unmut, "unmutated", "mutated"))
} else {
    post_cells$mu_freq_heavy <- NULL; post_cells$mu_count_heavy <- NULL; post_cells$mu_freq_light <- NULL
}
if (length(opt$lineage)) {
    lin <- bind_rows(lapply(opt$lineage, read_tsv))
    post_cells$in_lineage_tree <- !is.na(post_cells$clone_id) & post_cells$clone_id %in% lin$clone_id
}

msg(sprintf("post-QC: %s cells, %s paired, %s clones (%s expanded)",
            format(nrow(post_cells), big.mark = ","), format(sum(post_cells$paired_bcr), big.mark = ","),
            format(n_distinct(cl$clone_id), big.mark = ","), format(sum(by_subject$clone_size > 1), big.mark = ",")))

# ── per-sample summaries ──────────────────────────────────────────────────────
summarise_stage <- function(cells) {
    # Column checks happen out here: inside summarise(), `cells` is the n() column below.
    has <- names(cells)
    cells %>% group_by(sample, subject) %>%
        summarise(cells = n(), paired = sum(paired_bcr),
                  pct_paired = round(100 * sum(paired_bcr) / n(), 2),
                  unique_CTaa = n_distinct(CTaa),
                  unique_clones = if ("clone_id" %in% has) n_distinct(na.omit(clone_id)) else NA_integer_,
                  cells_in_expanded_clones = if ("clone_size" %in% has) sum(clone_size > 1, na.rm = TRUE) else NA_integer_,
                  median_mu_freq_heavy = if ("mu_freq_heavy" %in% has) round(median(mu_freq_heavy, na.rm = TRUE), 5) else NA_real_,
                  IgM = sum(isotype_class == "IgM", na.rm = TRUE), IgD = sum(isotype_class == "IgD", na.rm = TRUE),
                  IgG = sum(isotype_class == "IgG", na.rm = TRUE), IgA = sum(isotype_class == "IgA", na.rm = TRUE),
                  IgE = sum(isotype_class == "IgE", na.rm = TRUE), .groups = "drop") %>%
        as.data.frame()
}

# ── writers ───────────────────────────────────────────────────────────────────
write_stage <- function(stage, cells, contigs) {
    fwrite(cells, paste0(stage, "_cells.tsv"), sep = "\t", na = "")
    fwrite(summarise_stage(cells), paste0(stage, "_summary.tsv"), sep = "\t", na = "")

    if (requireNamespace("scRepertoire", quietly = TRUE)) {
        # combineBCR parses Cell Ranger's own columns positionally - feed it the raw frame.
        cr_cols <- intersect(c("barcode", "is_cell", "contig_id", "high_confidence", "length", "chain",
                               "v_gene", "d_gene", "j_gene", "c_gene", "full_length", "productive",
                               "cdr3", "cdr3_nt", "reads", "umis", "raw_clonotype_id", "raw_consensus_id"),
                             names(contigs))
        lst <- split(contigs[, cr_cols, drop = FALSE], contigs$sample)
        lst <- lapply(lst, function(d) { d$reads <- num(d$reads); d$umis <- num(d$umis); d })
        fun <- get("combineBCR", asNamespace("scRepertoire")); fml <- names(formals(fun))
        args <- list(lst, samples = names(lst))
        # Clones are SCOPer's (in *_cells.tsv); skip scRepertoire's own related-clone search.
        if ("call.related.clones" %in% fml) args$call.related.clones <- FALSE
        if ("removeNA" %in% fml)    args$removeNA    <- FALSE
        if ("removeMulti" %in% fml) args$removeMulti <- FALSE
        if ("filterMulti" %in% fml) args$filterMulti <- FALSE
        obj <- tryCatch(do.call(fun, args), error = function(e) {
            msg("WARN combineBCR() failed for ", stage, ": ", conditionMessage(e)); NULL })
        if (!is.null(obj)) saveRDS(obj, paste0(stage, "_combineBCR.rds"))
    } else msg("WARN scRepertoire unavailable - no ", stage, "_combineBCR.rds")

    if (requireNamespace("Seurat", quietly = TRUE)) {
        tryCatch({
            # Two rows, not one: Seurat 5's Assay5 rejects a single-row layer
            # ("Layers must be two-dimensional objects"). Same as TCRtoolkit.
            m <- Matrix::Matrix(0, nrow = 2, ncol = nrow(cells), sparse = TRUE,
                                dimnames = list(c("PLACEHOLDER-no-GEX-1", "PLACEHOLDER-no-GEX-2"), cells$cell_id))
            md <- cells; rownames(md) <- md$cell_id
            seu <- Seurat::CreateSeuratObject(counts = m, meta.data = md, project = opt$project,
                                              assay = "BCR", min.cells = 0, min.features = 0)
            seu@misc$provenance <- list(stage = stage, created = as.character(Sys.time()),
                note = paste("VDJ-only: no gene expression. The counts assay is a zero placeholder.",
                             "Do not normalize, run PCA/UMAP or cluster on this object."))
            saveRDS(seu, paste0(stage, "_seurat.rds"))
        }, error = function(e) msg("WARN Seurat object failed for ", stage, ": ", conditionMessage(e)))
    }
    msg(sprintf("wrote %s_* (%s cells)", stage, format(nrow(cells), big.mark = ",")))
}

# ── per-sample metadata from the samplesheet ──────────────────────────────────
sheet_meta <- NULL
if (!is.na(opt$sheet) && nzchar(opt$sheet) && file.exists(opt$sheet) && file.size(opt$sheet) > 0) {
    sh <- as.data.frame(fread(opt$sheet, sep = ",", colClasses = "character", na.strings = c("", "NA")))
    skip <- c("subject", "cr_id", "modality", "fastq_1", "fastq_2", "contig_annotations", "contig_fasta")
    keep <- setdiff(names(sh), skip)
    if ("sample" %in% keep && length(keep) > 1) {
        sh <- sh[, keep, drop = FALSE]
        # FASTQ sheets have one row per lane; one row per sample, first non-empty value.
        sheet_meta <- sh %>% group_by(sample) %>%
            summarise(across(everything(), ~ first(na.omit(.x))), .groups = "drop") %>%
            as.data.frame()
        msg("sample metadata from samplesheet: ", paste(setdiff(names(sheet_meta), "sample"), collapse = ", "))
    }
}
add_sample_meta <- function(cells) {
    if (is.null(sheet_meta)) return(cells)
    cols <- c("sample", setdiff(names(sheet_meta), names(cells)))
    if (length(cols) < 2) return(cells)
    out <- left_join(cells, sheet_meta[, cols, drop = FALSE], by = "sample")
    # metadata sits after the identity columns, before the receptor fields
    id <- intersect(c("cell_id", "sample", "subject", "barcode"), names(out))
    out[, c(id, setdiff(cols, "sample"), setdiff(names(out), c(id, cols))), drop = FALSE]
}
pre_cells  <- add_sample_meta(pre_cells)
post_cells <- add_sample_meta(post_cells)

post_contigs <- read_tsv(opt$post)
write_stage("pre_qc", pre_cells, pre)
write_stage("post_qc", post_cells, post_contigs)
msg("done")
