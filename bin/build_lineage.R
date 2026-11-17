#!/usr/bin/env Rscript
#
# build_lineage.R - B-cell lineage trees for ONE subject with dowser.
#
# Only clones with >= --min-seqs distinct sequences get a tree (a 1-2 sequence "tree" has
# no topology to read). --max-clones caps the work at the N largest clones (0 = all),
# because runtime grows with both clone count and clone size.
#
# Build methods, fastest first:
#   pratchet  maximum parsimony (phangorn). Minutes. Default.
#   pml       maximum likelihood (phangorn). Several times slower.
#   igphyml   IgPhyML, a B-cell-specific SHM model. Most accurate, hours on big subjects.
#
# Chains: 'auto' builds heavy+light (HL) trees when most cells are paired and falls back
# to heavy-only (H) if that fails. HL trees resolve more of the topology, because the
# light chain carries independent mutations.
#
# Outputs: <subject>_trees.rds (dowser airrClone tibble incl. phylo objects),
#          <subject>_lineage_summary.tsv (one row per tree), <subject>_top_trees.pdf
#
# Usage: build_lineage.R --subject S --input germ-pass.tsv --build pratchet --min-seqs 3
#                        --max-clones 0 --chain auto --plot-top 10 --nproc N
#                        [--igphyml /path/to/igphyml]

suppressPackageStartupMessages({
    library(airr)
    library(dplyr)
    library(dowser)
})

.args <- commandArgs(trailingOnly = TRUE)
.get <- function(flag, default = NA_character_) {
    i <- match(flag, .args)
    if (!is.na(i) && length(.args) > i) .args[i + 1] else default
}
subject  <- .get("--subject", "subject")
input    <- .get("--input")
build    <- .get("--build", "pratchet")
min_seqs <- as.integer(.get("--min-seqs", "3"))
max_cl   <- as.integer(.get("--max-clones", "0"))
chain    <- toupper(.get("--chain", "auto"))
plot_top <- as.integer(.get("--plot-top", "10"))
igphyml  <- .get("--igphyml", "/usr/local/share/igphyml/src/igphyml")
nproc    <- as.integer(.get("--nproc", "1"))
msg <- function(...) message("[lineage] ", subject, ": ", ...)

# airr::read_rearrangement() types only AIRR-standard columns. Change-O's --extended
# columns (v/d/j_germline_length, *_start/*_end) and our own come back as character, and
# dowser/shazam then fail on arithmetic ("invalid 'type' (character) of argument").
fix_types <- function(db) {
    num <- grep("(_start|_end|_length|_count|^mu_freq|^mu_count)", names(db), value = TRUE)
    for (c in num) if (is.character(db[[c]])) db[[c]] <- suppressWarnings(as.numeric(db[[c]]))
    db
}

summary_file <- paste0(subject, "_lineage_summary.tsv")
empty_summary <- function(note) {
    msg(note)
    write.table(data.frame(subject = character(0), clone_id = character(0), tree_id = character(0), n_seqs = integer(0),
                           n_cells = integer(0), tree_tips = integer(0), tree_length = numeric(0),
                           chain = character(0), build = character(0)),
                summary_file, sep = "\t", quote = FALSE, row.names = FALSE)
    quit(save = "no")
}

db <- if (file.size(input) > 0) fix_types(read_rearrangement(input)) else data.frame()
if (!nrow(db) || !"germline_alignment_d_mask" %in% names(db)) empty_summary("no germlines - no trees")
db <- db[!is.na(db$clone_id) & !is.na(db$germline_alignment_d_mask), , drop = FALSE]

# Clones big enough to be worth a tree, largest first, capped.
eligible <- db %>% filter(locus == "IGH") %>% count(clone_id, name = "n") %>%
    filter(n >= min_seqs) %>% arrange(desc(n))
if (max_cl > 0) eligible <- head(eligible, max_cl)
if (!nrow(eligible)) empty_summary(sprintf("no clone has >= %d heavy chains - no trees", min_seqs))
db <- db[db$clone_id %in% eligible$clone_id, , drop = FALSE]
msg(sprintf("%d clones eligible (>= %d sequences%s)", nrow(eligible), min_seqs,
            if (max_cl > 0) paste0(", capped at ", max_cl) else ""))

traits <- intersect(c("sample", "c_call"), names(db))
# dowser's single-cell recipe: heavy+light trees need light-chain subgroups
# (resolveLightChains, normally already run in DEFINE_CLONES) and split_light/no collapse.
fmt <- function(ch) {
    if (ch == "HL") {
        if (!"clone_subgroup" %in% names(db)) db <- resolveLightChains(db, nproc = nproc)
        formatClones(db, chain = "HL", traits = traits, minseq = min_seqs, cell = "cell_id",
                     locus = "locus", heavy = "IGH", split_light = TRUE, collapse = FALSE,
                     nproc = nproc)
    } else {
        formatClones(db[db$locus == "IGH", , drop = FALSE], chain = "H", traits = traits,
                     minseq = min_seqs, collapse = TRUE, nproc = nproc)
    }
}
if (chain == "AUTO") {
    paired <- db %>% group_by(cell_id) %>% summarise(p = any(locus == "IGH") && any(locus != "IGH"))
    chain <- if (mean(paired$p) >= 0.5) "HL" else "H"
}
clones <- tryCatch(fmt(chain), error = function(e) {
    if (chain == "HL") { msg("HL formatting failed (", conditionMessage(e), "); using heavy only"); NULL }
    else stop(e)
})
if (is.null(clones)) { chain <- "H"; clones <- fmt("H") }
if (!nrow(clones)) empty_summary("formatClones returned no clones")

t0 <- Sys.time()
trees <- if (build == "igphyml") {
    getTrees(clones, build = "igphyml", exec = igphyml, nproc = nproc)
} else {
    getTrees(clones, build = build, nproc = nproc)
}
msg(sprintf("%d %s trees built with %s in %.1f min", nrow(trees), chain, build,
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
saveRDS(trees, paste0(subject, "_trees.rds"))

# HL trees are built per light-chain SUBGROUP, and dowser names them <clone_id>_<subgroup>
# (P1_2 -> P1_2_1). Map every tree back to its parent clone so clone_id joins to cells;
# keep the tree's own id in tree_id.
parent <- setNames(as.character(db$clone_id), as.character(db$clone_id))
if ("clone_subgroup" %in% names(db)) {
    sg <- unique(db[!is.na(db$clone_subgroup), c("clone_id", "clone_subgroup")])
    parent <- c(parent, setNames(as.character(sg$clone_id), paste0(sg$clone_id, "_", sg$clone_subgroup)))
}
tree_id <- as.character(trees$clone_id)
cells_per_clone <- db %>% distinct(clone_id, cell_id) %>% count(clone_id, name = "n_cells")
summ <- data.frame(
    subject   = subject,
    clone_id  = ifelse(tree_id %in% names(parent), parent[tree_id], tree_id),
    tree_id   = tree_id,
    n_seqs    = trees$seqs,
    tree_tips = vapply(trees$trees, function(t) length(t$tip.label), integer(1)),
    tree_length = vapply(trees$trees, function(t) sum(t$edge.length), numeric(1)),
    chain = chain, build = build, stringsAsFactors = FALSE
) %>% left_join(cells_per_clone, by = "clone_id")
write.table(summ, summary_file, sep = "\t", quote = FALSE, row.names = FALSE)

if (plot_top > 0 && requireNamespace("ggtree", quietly = TRUE)) {
    top <- head(trees[order(-trees$seqs), ], plot_top)
    tip <- if ("c_call" %in% traits) "c_call" else NULL
    plots <- tryCatch(plotTrees(top, tips = tip, tipsize = 2), error = function(e) NULL)
    if (!is.null(plots)) {
        pdf(paste0(subject, "_top_trees.pdf"), width = 7, height = 6)
        for (i in seq_along(plots)) print(plots[[i]] + ggplot2::ggtitle(top$clone_id[i]))
        invisible(dev.off())
    }
}
