#!/usr/bin/env Rscript
#
# germline_shm.R - clonal germline reconstruction and somatic hypermutation, ONE subject.
#
#   germlines  dowser::createGermlines builds one inferred germline per CLONE (not per
#              sequence) from the IMGT references, with D-region masking. Lineage trees
#              root on this germline, so it runs whenever --run-shm or --run-lineage is on.
#   SHM        shazam::observedMutations against germline_alignment_d_mask, over the
#              whole V(D)J (regionDefinition = NULL), as a frequency (mu_freq) and a
#              count (mu_count). Also split into CDR/FWR replacement vs silent frequency
#              over the IMGT V segment, which is what selection analyses read.
#
# A clone whose germline cannot be built (unknown allele, trimmed V start) keeps its rows
# with NA germline and NA mutation values rather than being dropped.
#
# Usage: germline_shm.R --subject S --input clone-pass.tsv --germline-dir DIR
#                       --species human --run-shm true --nproc N

suppressPackageStartupMessages({
    library(airr)
    library(dplyr)
    library(dowser)
    library(shazam)
})

.args <- commandArgs(trailingOnly = TRUE)
.get <- function(flag, default = NA_character_) {
    i <- match(flag, .args)
    if (!is.na(i) && length(.args) > i) .args[i + 1] else default
}
subject <- .get("--subject", "subject")
input   <- .get("--input")
gdir    <- .get("--germline-dir", "/usr/local/share/germlines/imgt")
species <- .get("--species", "human")
run_shm <- tolower(.get("--run-shm", "true")) %in% c("true", "1", "yes")
nproc   <- as.integer(.get("--nproc", "1"))
msg <- function(...) message("[germline_shm] ", subject, ": ", ...)

# airr::read_rearrangement() types only AIRR-standard columns. Change-O's --extended
# columns (v/d/j_germline_length, *_start/*_end) and our own come back as character, and
# dowser/shazam then fail on arithmetic ("invalid 'type' (character) of argument").
fix_types <- function(db) {
    num <- grep("(_start|_end|_length|_count|^mu_freq|^mu_count)", names(db), value = TRUE)
    for (c in num) if (is.character(db[[c]])) db[[c]] <- suppressWarnings(as.numeric(db[[c]]))
    db
}
out <- paste0(subject, "_germ-pass.tsv")

db <- if (file.size(input) > 0) fix_types(read_rearrangement(input)) else data.frame()
if (!nrow(db) || all(is.na(db$clone_id))) {
    msg("no clones - passing the table through")
    write_rearrangement(db, out); quit(save = "no")
}

ref_dir <- file.path(gdir, species, "vdj")
if (!dir.exists(ref_dir)) stop("IMGT germline directory not found: ", ref_dir,
                               " (set --imgt_germlines / --species)")
refs <- readIMGT(dir = ref_dir, quiet = TRUE)

cloned   <- db[!is.na(db$clone_id), , drop = FALSE]
uncloned <- db[ is.na(db$clone_id), , drop = FALSE]
# One locus at a time: createGermlines stops on the first locus it cannot handle (e.g. an
# empty reference set), which would otherwise cost every other locus its germlines too.
germ <- bind_rows(lapply(split(cloned, cloned$locus), function(d) {
    tryCatch(createGermlines(d, references = refs, locus = "locus", clone = "clone_id",
                             nproc = nproc, verbose = 0),
             error = function(e) {
                 msg(sprintf("createGermlines failed for %s (%d sequences): %s",
                             d$locus[1], nrow(d), conditionMessage(e))); NULL })
}))
if (!nrow(germ)) germ <- cloned[0, , drop = FALSE]
# createGermlines drops rows it cannot reconstruct; put them back with NA germlines.
lost <- cloned[!cloned$sequence_id %in% germ$sequence_id, , drop = FALSE]
if (nrow(lost)) msg(sprintf("%d sequences in %d clones had no reconstructable germline",
                            nrow(lost), n_distinct(lost$clone_id)))
db <- bind_rows(germ, lost, uncloned)

if (run_shm) {
    ok <- !is.na(db$germline_alignment_d_mask) & nzchar(db$germline_alignment_d_mask)
    mut <- db[ok, , drop = FALSE]
    if (nrow(mut)) {
        mut <- observedMutations(mut, sequenceColumn = "sequence_alignment",
                                 germlineColumn = "germline_alignment_d_mask",
                                 regionDefinition = NULL, frequency = TRUE,
                                 combine = TRUE, nproc = nproc)            # -> mu_freq
        mut <- observedMutations(mut, sequenceColumn = "sequence_alignment",
                                 germlineColumn = "germline_alignment_d_mask",
                                 regionDefinition = NULL, frequency = FALSE,
                                 combine = TRUE, nproc = nproc)            # -> mu_count
        mut <- tryCatch(
            observedMutations(mut, sequenceColumn = "sequence_alignment",
                              germlineColumn = "germline_alignment_d_mask",
                              regionDefinition = IMGT_V, frequency = TRUE,
                              combine = FALSE, nproc = nproc),             # mu_freq_cdr_r ...
            error = function(e) { msg("CDR/FWR split skipped: ", conditionMessage(e)); mut })
    }
    db <- bind_rows(mut, db[!ok, , drop = FALSE])
    h <- db$locus == "IGH" & !is.na(db$mu_freq)
    msg(sprintf("SHM on %d/%d sequences; median heavy-chain mu_freq %.4f",
                sum(!is.na(db$mu_freq)), nrow(db),
                if (any(h)) median(db$mu_freq[h]) else NA_real_))
}

write_rearrangement(db, out)
