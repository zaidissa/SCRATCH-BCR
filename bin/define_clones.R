#!/usr/bin/env Rscript
#
# define_clones.R - B-cell clonal grouping for ONE subject (all of its samples pooled).
#
# Clones are called per subject, not per sample: the same B-cell lineage recurs across a
# patient's timepoints and tissues, and that sharing is usually the question. Different
# subjects never share clones by construction, so they are processed independently.
#
#   1. threshold  shazam::distToNearest on heavy-chain junctions (single-cell aware), then
#                 shazam::findThreshold. 'auto' tries the GMM fit, then the density
#                 method, then --fallback-threshold - small subjects often have too few
#                 clonally related cells for a bimodal distance distribution to exist.
#   2. clones     scoper::hierarchicalClones, single-linkage on length-normalised
#                 Hamming distance within V/J/junction-length groups, heavy chain only.
#                 Light chains do not split clones here (scoper >= 1.6 deprecates
#                 split_light); dowser::resolveLightChains then assigns light-chain
#                 subgroups (clone_subgroup) within each clone, the current Immcantation
#                 recipe.
#
# Threshold sanity: a fitted threshold above --warn-above (0.25) tends to merge unrelated
# cells into one clone; below --warn-below (0.05) it tends to split real clones. Neither
# stops the run - the subject's _threshold.tsv gets a threshold_note, the log a WARNING,
# and the distance plot shows where the cut fell. Re-run with --clonal_threshold <num>.
#
# Usage: define_clones.R --subject S --threshold auto|<num> --fallback-threshold 0.15
#                        --method hierarchical|identical --nproc N
#                        [--warn-above 0.25] [--warn-below 0.05] <airr_qc.tsv> [...]

suppressPackageStartupMessages({
    library(airr)
    library(dplyr)
    library(shazam)
    library(scoper)
})

.args <- commandArgs(trailingOnly = TRUE)
.get <- function(flag, default = NA_character_) {
    i <- match(flag, .args)
    if (!is.na(i) && length(.args) > i) .args[i + 1] else default
}
flags   <- c("--subject", "--threshold", "--fallback-threshold", "--method", "--nproc",
             "--warn-above", "--warn-below")
inputs  <- .args[!(.args %in% flags) & !(seq_along(.args) - 1) %in% match(flags, .args)]
subject <- .get("--subject", "subject")
thr_arg <- .get("--threshold", "auto")
fallback <- as.numeric(.get("--fallback-threshold", "0.15"))
warn_hi <- as.numeric(.get("--warn-above", "0.25"))
warn_lo <- as.numeric(.get("--warn-below", "0.05"))
method  <- .get("--method", "hierarchical")
nproc   <- as.integer(.get("--nproc", "1"))
msg <- function(...) message("[define_clones] ", subject, ": ", ...)

# airr::read_rearrangement() types only AIRR-standard columns. Change-O's --extended
# columns (v/d/j_germline_length, *_start/*_end) and our own come back as character, and
# dowser/shazam then fail on arithmetic ("invalid 'type' (character) of argument").
fix_types <- function(db) {
    num <- grep("(_start|_end|_length|_count|^mu_freq|^mu_count)", names(db), value = TRUE)
    for (c in num) if (is.character(db[[c]])) db[[c]] <- suppressWarnings(as.numeric(db[[c]]))
    db
}

db <- bind_rows(lapply(inputs, function(f) {
    if (file.size(f) == 0) return(NULL)
    fix_types(read_rearrangement(f))
}))
msg(sprintf("%d contigs, %d cells from %d sample file(s)",
            nrow(db), n_distinct(db$cell_id), length(inputs)))

thr_row <- data.frame(subject = subject, n_cells = n_distinct(db$cell_id),
                      threshold = NA_real_, threshold_method = NA_character_,
                      threshold_note = NA_character_)

if (nrow(db) == 0 || !any(db$locus == "IGH")) {
    msg("no heavy chains - writing an empty clone table")
    db$clone_id <- character(0)[seq_len(nrow(db))]
    write_rearrangement(db, paste0(subject, "_clone-pass.tsv"))
    write.table(thr_row, paste0(subject, "_threshold.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    quit(save = "no")
}

# ── 1. threshold ──────────────────────────────────────────────────────────────
threshold <- suppressWarnings(as.numeric(thr_arg))
if (method == "identical") {
    threshold <- 0; thr_row$threshold_method <- "identical"
} else if (!is.na(threshold)) {
    thr_row$threshold_method <- "user"
} else {
    dist <- tryCatch(
        distToNearest(db, sequenceColumn = "junction", vCallColumn = "v_call",
                      jCallColumn = "j_call", model = "ham", normalize = "len",
                      cellIdColumn = "cell_id", locusColumn = "locus",
                      onlyHeavy = TRUE, VJthenLen = FALSE, nproc = nproc),
        error = function(e) { msg("distToNearest failed: ", conditionMessage(e)); NULL })
    d <- if (!is.null(dist)) dist$dist_nearest[!is.na(dist$dist_nearest)] else numeric(0)
    fit_ok <- function(t) !is.null(t) && length(t@threshold) && is.finite(t@threshold) &&
        t@threshold > 0 && t@threshold < 0.5
    for (m in c("gmm", "density")) {
        if (length(d) < 20) break                      # too few neighbours to fit anything
        t <- tryCatch(findThreshold(d, method = m), error = function(e) NULL)
        if (fit_ok(t)) { threshold <- t@threshold; thr_row$threshold_method <- m; break }
    }
    if (is.na(threshold)) {
        threshold <- fallback; thr_row$threshold_method <- "fallback"
        msg(sprintf("no usable threshold fit (%d nearest-neighbour distances); using fallback %.3f",
                    length(d), fallback))
    }
    if (!is.null(dist) && length(d) && requireNamespace("ggplot2", quietly = TRUE)) {
        p <- ggplot2::ggplot(data.frame(d = d), ggplot2::aes(d)) +
            ggplot2::geom_histogram(binwidth = 0.02, fill = "grey60", colour = "white") +
            ggplot2::geom_vline(xintercept = threshold, colour = "firebrick", linetype = 2) +
            ggplot2::labs(title = paste0(subject, ": heavy-chain distance to nearest"),
                          subtitle = sprintf("threshold %.3f (%s)", threshold, thr_row$threshold_method),
                          x = "Hamming distance (length-normalised)", y = "cells") +
            ggplot2::theme_bw()
        ggplot2::ggsave(paste0(subject, "_dist_to_nearest.png"), p, width = 6, height = 4, dpi = 150)
    }
}
thr_row$threshold <- threshold
msg(sprintf("threshold %.3f (%s)", threshold, thr_row$threshold_method))
# Only fitted thresholds are second-guessed: a user value is deliberate, and 'identical'
# is 0 by definition.
if (thr_row$threshold_method %in% c("gmm", "density", "fallback")) {
    note <- if (threshold > warn_hi) {
        sprintf("high: %.3f > %.2f may merge unrelated cells into one clone", threshold, warn_hi)
    } else if (threshold < warn_lo) {
        sprintf("low: %.3f < %.2f may split real clones", threshold, warn_lo)
    } else NA_character_
    if (!is.na(note)) {
        thr_row$threshold_note <- note
        msg("WARNING threshold ", note, " - check ", subject, "_dist_to_nearest.png and ",
            "consider --clonal_threshold <num>")
    }
}

# ── 2. clones ─────────────────────────────────────────────────────────────────
clones <- if (method == "identical") {
    identicalClones(db, method = "nt", junction = "junction", v_call = "v_call",
                    j_call = "j_call", cell_id = "cell_id", locus = "locus",
                    only_heavy = TRUE, split_light = FALSE, summarize_clones = FALSE, nproc = nproc)
} else {
    hierarchicalClones(db, threshold = threshold, method = "nt", linkage = "single",
                       normalize = "len", junction = "junction", v_call = "v_call",
                       j_call = "j_call", cell_id = "cell_id", locus = "locus",
                       only_heavy = TRUE, split_light = FALSE, summarize_clones = FALSE,
                       nproc = nproc)
}
clones <- as.data.frame(clones)
# Light-chain subgroups within each heavy-chain clone (cells whose light chains disagree
# land in different subgroups). Best effort: older dowser lacks it.
if (requireNamespace("dowser", quietly = TRUE) &&
    "resolveLightChains" %in% getNamespaceExports("dowser") && any(clones$locus != "IGH")) {
    clones <- tryCatch(as.data.frame(dowser::resolveLightChains(clones, nproc = nproc)),
                       error = function(e) { msg("resolveLightChains skipped: ", conditionMessage(e)); clones })
}
# Subject-scoped ids: scoper numbers clones 1..n per call, so two subjects would both
# have a clone "1". Prefixing makes clone_id unique across the whole cohort.
clones$clone_id <- ifelse(is.na(clones$clone_id), NA_character_,
                          paste0(subject, "_", clones$clone_id))

n_clones <- n_distinct(clones$clone_id, na.rm = TRUE)
sizes <- clones %>% filter(!is.na(clone_id)) %>% distinct(clone_id, cell_id) %>% count(clone_id)
thr_row$n_clones <- n_clones
thr_row$n_expanded_clones <- sum(sizes$n > 1)
thr_row$largest_clone_cells <- if (nrow(sizes)) max(sizes$n) else 0L
msg(sprintf("%d clones, %d expanded (>1 cell), largest %d cells",
            n_clones, thr_row$n_expanded_clones, thr_row$largest_clone_cells))

write_rearrangement(clones, paste0(subject, "_clone-pass.tsv"))
write.table(thr_row, paste0(subject, "_threshold.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
