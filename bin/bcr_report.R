#!/usr/bin/env Rscript
#
# bcr_report.R - one self-contained HTML report for a SCRATCH-BCR run.
#
# Sections: sample status (why each sample is in, flagged or out), cell QC, isotypes,
# clonal expansion, clones shared across samples/time points, SHM, IGHV usage, the
# clonal threshold, lineage trees. Every figure is drawn twice (light and dark page) and
# embedded as PNG, with its numbers in a table beside it, so the file opens anywhere,
# offline, with nothing else.
#
# Usage: bcr_report.R --cells post_qc_cells.tsv --qc bcr_qc_summary.tsv
#            [--status sample_status.tsv] [--thresholds a_threshold.tsv,...]
#            [--dist-plots a_dist_to_nearest.png,...] [--lineage a_lineage_summary.tsv,...]
#            [--trees a_trees.rds,...] [--project NAME] [--version V] [--out bcr_report.html]

suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(base64enc)
})

.args <- commandArgs(trailingOnly = TRUE)
.get <- function(flag, default = "") {
    i <- match(flag, .args)
    if (!is.na(i) && length(.args) > i) .args[i + 1] else default
}
files <- function(flag) {
    f <- strsplit(.get(flag), ",", fixed = TRUE)[[1]]
    f[nzchar(f) & file.exists(f) & !grepl("^NO_FILE", basename(f))]
}
opt <- list(cells = .get("--cells"), qc = .get("--qc"), status = files("--status"),
            thresholds = files("--thresholds"), dist = files("--dist-plots"),
            lineage = files("--lineage"), trees = files("--trees"),
            project = .get("--project", "bcr"), version = .get("--version", ""),
            out = .get("--out", "bcr_report.html"))

read_tsv <- function(f) {
    if (!length(f) || !file.exists(f) || file.size(f) == 0) return(data.frame())
    read.delim(f, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE,
               na.strings = c("", "NA"), quote = "")
}
bind_files <- function(fs) bind_rows(lapply(fs, function(f) mutate(read_tsv(f), across(everything(), as.character))))
num <- function(x) suppressWarnings(as.numeric(gsub("[,%]", "", as.character(x))))
has <- function(d, ...) all(c(...) %in% names(d))

cells  <- read_tsv(opt$cells)
qc     <- read_tsv(opt$qc)
status <- if (length(opt$status)) read_tsv(opt$status[1]) else data.frame()
thr    <- bind_files(opt$thresholds)
for (c in c("stage", "reason", "flags", "detail", "cells_passed", "pct_paired", "cr_estimated_number_of_cells",
            "cr_valid_barcodes", "cr_reads_mapped_to_any_v_d_j_gene", "cr_fraction_reads_in_cells"))
    if (nrow(status) && !c %in% names(status)) status[[c]] <- NA
lin    <- bind_files(opt$lineage)

tidytext_reorder <- function(x, by, within) {
    k <- paste(x, within, sep = "___")
    factor(k, levels = unique(k[order(within, by)]))
}

# ── colours: the reference data-viz palette, light and dark steps ───────────────
MODES <- list(
    light = list(surface = "#fcfcfb", ink = "#0b0b0b", ink2 = "#52514e", muted = "#898781",
                 grid = "#e1e0d9", axis = "#c3c2b7",
                 cat = c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"),
                 other = "#c3c2b7"),
    dark  = list(surface = "#1a1a19", ink = "#f0efec", ink2 = "#c3c2b7", muted = "#898781",
                 grid = "#2c2c2a", axis = "#383835",
                 cat = c("#3987e5", "#d95926", "#199e70", "#c98500", "#d55181", "#008300", "#9085e9", "#e66767"),
                 other = "#55544f"))
SEQ <- c("#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b")
ISO <- c("IgM", "IgD", "IgG", "IgA", "IgE")
iso_cols <- function(m) c(setNames(MODES[[m]]$cat[1:5], ISO), Unknown = MODES[[m]]$other)

theme_mode <- function(m) {
    k <- MODES[[m]]
    theme_minimal(base_size = 11, base_family = "sans") +
        theme(plot.background = element_rect(fill = k$surface, colour = NA),
              panel.background = element_rect(fill = k$surface, colour = NA),
              panel.grid.major = element_line(colour = k$grid, linewidth = 0.3),
              panel.grid.minor = element_blank(),
              axis.line.x = element_line(colour = k$axis, linewidth = 0.4),
              text = element_text(colour = k$ink2), axis.text = element_text(colour = k$ink2),
              strip.text = element_text(colour = k$ink, face = "bold", hjust = 0),
              legend.position = "top", legend.justification = "left",
              legend.title = element_blank(), legend.key.size = unit(10, "pt"),
              plot.margin = margin(6, 12, 6, 6))
}

# ── samples: order and short labels ─────────────────────────────────────────────
# Within a patient, order by time point when it is numeric; label "patient · tp".
sample_info <- function() {
    src <- if (nrow(status)) status else if (nrow(qc)) qc else cells
    s <- distinct(data.frame(sample = src$sample,
                             subject = if ("subject" %in% names(src)) src$subject else NA_character_,
                             stringsAsFactors = FALSE))
    if (has(cells, "sample", "subject")) {
        cs <- distinct(cells[, c("sample", "subject")])
        hit <- match(s$sample, cs$sample)
        s$subject[!is.na(hit)] <- cs$subject[hit[!is.na(hit)]]
    }
    s <- s[!duplicated(s$sample), ]
    tp <- if (has(cells, "sample", "timepoint")) distinct(cells[, c("sample", "timepoint")]) else NULL
    if (!is.null(tp)) s <- left_join(s, tp[!duplicated(tp$sample), ], by = "sample") else s$timepoint <- NA
    s$tp_num <- num(gsub("^[^0-9.-]*", "", s$timepoint))
    s <- s[order(s$subject, is.na(s$tp_num), s$tp_num, s$sample), ]
    short <- ifelse(!is.na(s$timepoint), paste0(s$subject, " tp", s$timepoint), s$sample)   # ASCII: PNG text in the container is not UTF-8
    s$label <- if (anyDuplicated(short)) s$sample else short
    s
}
S <- sample_info()
lab <- function(x) factor(S$label[match(x, S$sample)], levels = rev(S$label))
labs_fwd <- function(x) factor(S$label[match(x, S$sample)], levels = S$label)

# ── HTML helpers ────────────────────────────────────────────────────────────────
esc <- function(x) { x <- as.character(x); x[is.na(x)] <- ""
    x <- gsub("&", "&amp;", x, fixed = TRUE); x <- gsub("<", "&lt;", x, fixed = TRUE)
    gsub(">", "&gt;", x, fixed = TRUE) }
b64 <- function(f, type = "image/png") paste0("data:", type, ";base64,", base64encode(f))

fig <- function(make, w = 8, h = 4, alt = "", caption = "") {
    uri <- vapply(c("light", "dark"), function(m) {
        f <- tempfile(fileext = ".png")
        p <- tryCatch(make(m), error = function(e) { message("figure skipped: ", conditionMessage(e)); NULL })
        if (is.null(p)) return("")
        png(f, width = w * 144, height = h * 144, res = 144, bg = MODES[[m]]$surface)
        print(p); invisible(dev.off()); b64(f)
    }, character(1))
    if (!nzchar(uri["light"])) return("")
    sprintf('<figure><picture><source srcset="%s" media="(prefers-color-scheme: dark)"><img src="%s" alt="%s" width="%d" height="%d"></picture>%s</figure>',
            uri["dark"], uri["light"], esc(alt), as.integer(round(w * 100)), as.integer(round(h * 100)),
            if (nzchar(caption)) paste0("<figcaption>", caption, "</figcaption>") else "")
}

table_html <- function(d, cls = "", num_cols = NULL) {
    if (!nrow(d)) return('<p class="empty">Nothing to show.</p>')
    num_cols <- num_cols %||% names(d)[vapply(d, is.numeric, logical(1))]
    th <- paste0("<th", ifelse(names(d) %in% num_cols, ' class="n"', ""), ">", esc(names(d)), "</th>", collapse = "")
    rows <- apply(d, 1, function(r) paste0("<td", ifelse(names(d) %in% num_cols, ' class="n"', ""), ">",
                                          ifelse(grepl("^<span", r), r, esc(r)), "</td>", collapse = ""))
    sprintf('<div class="tw"><table class="%s"><thead><tr>%s</tr></thead><tbody>%s</tbody></table></div>',
            cls, th, paste0("<tr>", rows, "</tr>", collapse = ""))
}
`%||%` <- function(a, b) if (is.null(a)) b else a
numbers <- function(d) paste0("<details><summary>Show the numbers</summary>", table_html(d), "</details>")
section <- function(id, title, lead, ...) {
    body <- paste0(..., collapse = "")
    sprintf('<section id="%s"><h2>%s</h2>%s%s</section>', id, title,
            if (nzchar(lead)) paste0('<p class="lead">', lead, "</p>") else "", body)
}
badge <- function(s) {
    cls <- c(ok = "good", ok_low_quality = "warn", failed = "bad")[s]
    icon <- c(ok = "&#10003;", ok_low_quality = "!", failed = "&#10005;")[s]
    txt <- c(ok = "ok", ok_low_quality = "low quality", failed = "failed")[s]
    ifelse(is.na(cls), esc(s), sprintf('<span class="b %s"><i>%s</i>%s</span>', cls, icon, txt))
}
pct <- function(a, b) ifelse(b > 0, round(100 * a / b, 1), NA)
fmt <- function(x) format(x, big.mark = ",", trim = TRUE)
nsamp <- nrow(S)
hgt <- function(per = 0.32, base = 1.4, n = nsamp) max(2.4, base + per * n)

# ── 1. overview tiles ───────────────────────────────────────────────────────────
n_cells <- nrow(cells)
n_pair  <- if (has(cells, "paired_bcr")) sum(as.logical(cells$paired_bcr), na.rm = TRUE) else NA
n_clone <- if (has(cells, "clone_id")) n_distinct(na.omit(cells$clone_id)) else NA
n_exp   <- if (has(cells, "clone_id", "clone_size")) n_distinct(cells$clone_id[!is.na(cells$clone_size) & cells$clone_size >= 2]) else NA
st_tab  <- if (nrow(status)) table(factor(status$status, c("ok", "ok_low_quality", "failed"))) else NULL
tile <- function(v, l, sub = "") sprintf('<div class="tile"><div class="v">%s</div><div class="l">%s</div>%s</div>',
                                         v, l, if (nzchar(sub)) paste0('<div class="s">', sub, "</div>") else "")
tiles <- paste0('<div class="tiles">',
    tile(fmt(nsamp), "samples", if (!is.null(st_tab)) sprintf("%d ok · %d flagged · %d failed",
         st_tab[["ok"]], st_tab[["ok_low_quality"]], st_tab[["failed"]]) else ""),
    tile(fmt(n_distinct(S$subject)), "patients"),
    tile(fmt(n_cells), "B cells after QC", if (!is.na(n_pair)) sprintf("%s%% heavy + light paired", round(100 * n_pair / max(n_cells, 1))) else ""),
    tile(fmt(n_clone), "clones", if (!is.na(n_exp)) sprintf("%s expanded (2+ cells)", fmt(n_exp)) else ""),
    tile(fmt(nrow(lin)), "lineage trees"),
    "</div>")

# ── 2. sample status ────────────────────────────────────────────────────────────
sec_status <- if (nrow(status)) {
    d <- status %>% mutate(label = S$label[match(sample, S$sample)]) %>%
        transmute(Sample = sprintf('<span class="sl">%s</span><small>%s</small>', esc(label), esc(sample)),
                  id = sample, Status = badge(status),
                  `Stage / reason` = ifelse(status == "failed", paste0(stage, ": ", reason), ""),
                  Flags = gsub(",", ", ", ifelse(is.na(flags), "", flags)),
                  `Cells (Cell Ranger)` = num(cr_estimated_number_of_cells),
                  `Cells passed QC` = num(cells_passed), `% paired` = num(pct_paired),
                  `Valid barcodes` = cr_valid_barcodes, `Reads on V(D)J` = cr_reads_mapped_to_any_v_d_j_gene,
                  `Reads in cells` = cr_fraction_reads_in_cells, Detail = detail)
    d <- d[order(match(d$id, S$sample)), ]; d$id <- NULL
    d <- d[, vapply(d, function(x) any(!is.na(x) & nzchar(as.character(x))), logical(1)), drop = FALSE]
    section("status", "Sample status",
        "Every sample in the samplesheet: whether it reached the results, and why not. Flags mark weak samples that were <b>kept</b>; results from them rest on fewer or noisier cells.",
        table_html(d, "status", c("Cells (Cell Ranger)", "Cells passed QC", "% paired")),
        '<p class="note">Flag thresholds: few_cells &lt; 50 cells · low_valid_barcodes &lt; 75% · low_vdj_mapping &lt; 50% of reads · low_reads_in_cells &lt; 50% · low_read_depth &lt; 5,000 read pairs per cell · low_pairing &lt; 30% of cells paired. Failure reasons are read from Cell Ranger\'s error message; its full log is in <code>cellranger/&lt;sample&gt;/cellranger.log</code>.</p>')
} else ""

# ── 3. cell QC ──────────────────────────────────────────────────────────────────
sec_qc <- if (nrow(qc) && has(qc, "cells_total", "cells_passed")) {
    q <- qc %>% mutate(cells_total = num(cells_total), cells_passed = num(cells_passed),
                       failed = cells_total - cells_passed)
    long <- rbind(data.frame(sample = q$sample, what = "passed QC", n = q$cells_passed),
                  data.frame(sample = q$sample, what = "removed", n = q$failed))
    long$what <- factor(long$what, c("removed", "passed QC"))
    f1 <- fig(function(m) {
        k <- MODES[[m]]
        ggplot(long, aes(n, lab(sample), fill = what)) +
            geom_col(width = 0.68, colour = k$surface, linewidth = 0.4) +
            scale_fill_manual(values = c("passed QC" = k$cat[1], removed = k$other), breaks = c("passed QC", "removed")) +
            scale_x_continuous(expand = expansion(c(0, 0.05)), labels = scales::comma) +
            labs(x = "cells (Cell Ranger calls)", y = NULL) + theme_mode(m)
    }, 8, hgt(), "Cells per sample that passed or were removed by BCR QC")
    rc <- grep("^contigs_failed_", names(q), value = TRUE)
    f2 <- ""; reasons <- data.frame()
    if (length(rc)) {
        reasons <- do.call(rbind, lapply(rc, function(c) data.frame(sample = q$sample,
                     reason = sub("contigs_failed_", "", c), n = num(q[[c]]))))
        reasons <- reasons[!is.na(reasons$n) & reasons$n > 0, ]
        if (nrow(reasons)) {
            ord <- names(sort(tapply(reasons$n, reasons$reason, sum), decreasing = TRUE))
            reasons$reason <- factor(reasons$reason, ord)
            f2 <- fig(function(m) {
                k <- MODES[[m]]
                ggplot(reasons, aes(n, lab(sample), fill = reason)) +
                    geom_col(width = 0.68, colour = k$surface, linewidth = 0.4) +
                    scale_fill_manual(values = setNames(c(k$cat, rep(k$other, 8))[seq_along(ord)], ord)) +
                    scale_x_continuous(expand = expansion(c(0, 0.05)), labels = scales::comma) +
                    labs(x = "contigs removed", y = NULL) + theme_mode(m)
            }, 8, hgt(), "Contigs removed per sample by reason")
        }
    }
    tab <- q %>% transmute(Sample = S$label[match(sample, S$sample)], `Cells` = cells_total,
                           `Passed` = cells_passed, `% passed` = pct(cells_passed, cells_total),
                           `Paired H+L` = num(if ("cells_paired_heavy_light" %in% names(q)) cells_paired_heavy_light else NA),
                           `% paired` = num(if ("pct_paired" %in% names(q)) pct_paired else NA))
    section("qc", "Cell QC",
        "A cell is kept when it has one productive, annotated heavy chain (and its best light chain). Reasons are counted per contig, in the order the rules are applied.",
        f1, f2, if (nzchar(f2)) '<p class="note">no_heavy: light chain only · multi_heavy: two or more heavy chains (likely doublet) · non_productive · not_annotated (no IgBLAST call) · low_umi · extra_light: a cell\'s light chains beyond the best one.</p>' else "",
        numbers(tab))
} else ""

# ── 4. isotypes ─────────────────────────────────────────────────────────────────
sec_iso <- if (n_cells && has(cells, "isotype_class", "sample")) {
    iso <- cells %>% mutate(isotype_class = ifelse(isotype_class %in% ISO, isotype_class, "Unknown")) %>%
        count(sample, isotype_class) %>% group_by(sample) %>% mutate(pct = 100 * n / sum(n)) %>% ungroup() %>%
        mutate(isotype_class = factor(isotype_class, c(ISO, "Unknown")))
    lv <- levels(droplevels(iso$isotype_class))
    f <- fig(function(m) {
        ggplot(iso, aes(pct, lab(sample), fill = isotype_class)) +
            geom_col(width = 0.68, colour = MODES[[m]]$surface, linewidth = 0.4) +
            scale_fill_manual(values = iso_cols(m)[lv], breaks = lv) +
            scale_x_continuous(expand = expansion(c(0, 0.01)), labels = function(x) paste0(x, "%")) +
            labs(x = "% of B cells", y = NULL) + theme_mode(m)
    }, 8, hgt(), "Isotype composition of B cells per sample")
    tab <- iso %>% select(-pct) %>% tidyr::pivot_wider(names_from = isotype_class, values_from = n, values_fill = 0) %>%
        mutate(Sample = S$label[match(sample, S$sample)]) %>% select(Sample, any_of(c(ISO, "Unknown")))
    section("isotypes", "Isotypes",
        "Heavy-chain constant region per cell. IgG and IgA dominate in antigen-experienced, class-switched B cells and plasma cells; IgM/IgD mark naive or early responses.",
        f, numbers(tab))
} else ""

# ── 5. clonal expansion ─────────────────────────────────────────────────────────
sec_clones <- if (n_cells && has(cells, "clone_id", "clone_size", "sample")) {
    cs <- cells %>% filter(!is.na(clone_id)) %>%
        mutate(size = num(if ("clone_size_sample" %in% names(cells)) clone_size_sample else clone_size),
               grp = cut(size, c(0, 1, 2, 5, 20, Inf), c("1 cell", "2 cells", "3-5", "6-20", ">20"))) %>%
        count(sample, grp) %>% group_by(sample) %>% mutate(pct = 100 * n / sum(n)) %>% ungroup()
    gl <- levels(cs$grp)
    f <- fig(function(m) {
        ggplot(cs, aes(pct, lab(sample), fill = grp)) +
            geom_col(width = 0.68, colour = MODES[[m]]$surface, linewidth = 0.4) +
            scale_fill_manual(values = setNames(SEQ[c(2, 3, 4, 5, 6)], gl), drop = FALSE) +
            scale_x_continuous(expand = expansion(c(0, 0.01)), labels = function(x) paste0(x, "%")) +
            labs(x = "% of B cells, by the size of their clone in that sample", y = NULL) + theme_mode(m)
    }, 8, hgt(), "Share of B cells in clones of each size per sample")
    top <- cells %>% filter(!is.na(clone_id)) %>% group_by(subject, clone_id) %>%
        summarise(Cells = n(), Samples = paste(unique(S$label[match(sample, S$sample)][order(match(sample, S$sample))]), collapse = ", "),
                  Isotypes = paste(names(sort(table(isotype_class), decreasing = TRUE)), collapse = "/"),
                  IGHV = names(sort(table(ighv), decreasing = TRUE))[1] %||% NA,
                  `CDR3 (heavy)` = names(sort(table(cdr3_heavy), decreasing = TRUE))[1] %||% NA,
                  `Median SHM %` = if ("mu_freq_heavy" %in% names(cells)) round(100 * median(num(mu_freq_heavy), na.rm = TRUE), 1) else NA,
                  `In tree` = if ("in_lineage_tree" %in% names(cells)) ifelse(any(as.logical(in_lineage_tree)), "yes", "") else "",
                  .groups = "drop") %>%
        filter(Cells >= 2) %>% arrange(subject, desc(Cells)) %>% group_by(subject) %>% slice_head(n = 15) %>% ungroup() %>%
        rename(Patient = subject, Clone = clone_id)
    section("clones", "Clonal expansion",
        "Clones are called per patient (same V and J genes, same junction length, similar junction). An expanded clone, many cells from one ancestor, points to an active antigen-driven response.",
        f, numbers(cs %>% transmute(Sample = S$label[match(sample, S$sample)], `Clone size` = grp, Cells = n, `%` = round(pct, 1))),
        "<h3>Largest clones</h3>", table_html(top, "", c("Cells", "Median SHM %")))
} else ""

# ── 6. clones shared across samples / time points ───────────────────────────────
sec_shared <- if (n_cells && has(cells, "clone_id", "sample", "subject")) {
    parts <- lapply(split(cells[!is.na(cells$clone_id), ], cells$subject[!is.na(cells$clone_id)]), function(d) {
        p <- d$subject[1]; samples <- S$sample[S$sample %in% d$sample]
        if (length(samples) < 2) return(sprintf("<h3>%s</h3><p class=\"empty\">One sample only: nothing to share.</p>", esc(p)))
        cc <- count(d, clone_id, sample)
        multi <- cc %>% group_by(clone_id) %>% filter(n_distinct(sample) >= 2) %>% ungroup()
        # Pairwise: clones seen in both samples.
        sets <- split(cc$clone_id, cc$sample)
        pw <- expand.grid(a = samples, b = samples, stringsAsFactors = FALSE)
        pw$n <- mapply(function(a, b) length(intersect(sets[[a]], sets[[b]])), pw$a, pw$b)
        pw$diag <- pw$a == pw$b
        f_pw <- fig(function(m) {
            k <- MODES[[m]]
            ggplot(pw, aes(labs_fwd(b), lab(a))) +
                geom_tile(aes(fill = ifelse(diag, NA, n)), colour = k$surface, linewidth = 1) +
                geom_text(aes(label = ifelse(diag, paste0(n, "\nclones"), n),
                              colour = ifelse(!diag & n > max(pw$n[!pw$diag], 1) * 0.55, "w", "i")), size = 3.2, lineheight = 0.9) +
                scale_fill_gradient(low = SEQ[1], high = SEQ[6], limits = c(0, max(pw$n[!pw$diag], 1)),
                                    na.value = k$grid, guide = "none") +
                scale_colour_manual(values = c(w = "#ffffff", i = k$ink), guide = "none") +
                labs(x = NULL, y = NULL) + theme_mode(m) +
                theme(panel.grid = element_blank(), axis.line.x = element_blank(),
                      axis.text.x = element_text(angle = 30, hjust = 1))
        }, min(8, 2.4 + 0.9 * length(samples)), min(7, 1.6 + 0.75 * length(samples)),
            sprintf("Clones shared between each pair of samples for patient %s", p),
            "Diagonal: clones in that sample. Off-diagonal: clones found in both samples.")
        f_tr <- ""
        if (nrow(multi)) {
            keep <- multi %>% group_by(clone_id) %>% summarise(t = sum(n)) %>% arrange(desc(t)) %>% slice_head(n = 30)
            tr <- cc %>% filter(clone_id %in% keep$clone_id) %>%
                tidyr::complete(clone_id, sample = samples, fill = list(n = 0)) %>%
                mutate(clone_id = factor(clone_id, rev(keep$clone_id)))
            f_tr <- fig(function(m) {
                k <- MODES[[m]]
                ggplot(tr, aes(labs_fwd(sample), clone_id)) +
                    geom_tile(aes(fill = ifelse(n > 0, n, NA)), colour = k$surface, linewidth = 0.8) +
                    geom_text(aes(label = ifelse(n > 0, n, ""), colour = ifelse(n > max(tr$n) * 0.55, "w", "i")), size = 2.8) +
                    scale_fill_gradient(low = SEQ[1], high = SEQ[6], limits = c(0, max(tr$n)),
                                        na.value = k$grid, name = "cells") +
                    scale_colour_manual(values = c(w = "#ffffff", i = k$ink), guide = "none") +
                    labs(x = NULL, y = NULL) + theme_mode(m) +
                    theme(panel.grid = element_blank(), axis.line.x = element_blank(), legend.position = "right",
                          legend.title = element_text(), axis.text.x = element_text(angle = 30, hjust = 1))
            }, min(9, 3 + 0.9 * length(samples)), max(2.4, 1.2 + 0.24 * nrow(keep)),
                sprintf("Cells per sample for the clones of patient %s found in two or more samples", p),
                sprintf("%d clone%s found in 2+ samples%s; cells per sample. Read left to right to follow a clone over time.",
                        n_distinct(multi$clone_id), if (n_distinct(multi$clone_id) == 1) "" else "s",
                        if (n_distinct(multi$clone_id) > 30) " (largest 30 shown)" else ""))
        } else f_tr <- '<p class="empty">No clone was found in more than one sample.</p>'
        paste0("<h3>", esc(p), "</h3>", f_pw, f_tr)
    })
    section("shared", "Clones across samples and time points",
        "The same clone in several samples of one patient means a B-cell response that persists, or moves between sites. Counts are cells per sample.",
        paste(parts, collapse = ""))
} else ""

# ── 7. SHM ──────────────────────────────────────────────────────────────────────
sec_shm <- if (n_cells && has(cells, "mu_freq_heavy", "isotype_class")) {
    sh <- cells %>% mutate(mu = 100 * num(mu_freq_heavy),
                           isotype_class = factor(ifelse(isotype_class %in% ISO, isotype_class, "Unknown"), c(ISO, "Unknown"))) %>%
        filter(!is.na(mu))
    if (nrow(sh)) {
        lv <- levels(droplevels(sh$isotype_class))
        f1 <- fig(function(m) {
            k <- MODES[[m]]
            ggplot(sh, aes(mu, lab(sample))) +
                geom_boxplot(outlier.shape = NA, width = 0.55, colour = k$ink2, fill = NA, linewidth = 0.4) +
                geom_jitter(aes(colour = isotype_class), height = 0.18, width = 0, size = 1.3, alpha = 0.75) +
                scale_colour_manual(values = iso_cols(m)[lv], breaks = lv) +
                labs(x = "heavy-chain V-region mutation (% of positions)", y = NULL) + theme_mode(m)
        }, 8, hgt(0.42), "Somatic hypermutation per cell by sample, coloured by isotype")
        f2 <- fig(function(m) {
            k <- MODES[[m]]
            ggplot(sh, aes(mu, factor(isotype_class, rev(lv)), colour = isotype_class)) +
                geom_boxplot(outlier.shape = NA, width = 0.55, fill = NA, linewidth = 0.5) +
                geom_jitter(height = 0.18, width = 0, size = 1.2, alpha = 0.6) +
                scale_colour_manual(values = iso_cols(m)[lv], guide = "none") +
                labs(x = "heavy-chain V-region mutation (% of positions)", y = NULL) + theme_mode(m)
        }, 8, max(2.2, 1.2 + 0.4 * length(lv)), "Somatic hypermutation by isotype")
        tab <- sh %>% group_by(Sample = S$label[match(sample, S$sample)], Isotype = isotype_class) %>%
            summarise(Cells = n(), `Median %` = round(median(mu), 2), `IQR low` = round(quantile(mu, .25), 2),
                      `IQR high` = round(quantile(mu, .75), 2), .groups = "drop")
        section("shm", "Somatic hypermutation",
            "Mutations in the heavy-chain V gene against the inferred germline. Higher SHM means more rounds of selection in germinal centres; class-switched cells usually carry more than IgM/IgD cells.",
            f1, f2, numbers(tab))
    } else ""
} else ""

# ── 8. IGHV usage ───────────────────────────────────────────────────────────────
sec_v <- if (n_cells && has(cells, "ighv", "subject")) {
    v <- cells %>% filter(!is.na(ighv)) %>% count(subject, ighv) %>% group_by(subject) %>%
        mutate(pct = 100 * n / sum(n)) %>% arrange(desc(n), .by_group = TRUE) %>% slice_head(n = 15) %>% ungroup()
    np <- n_distinct(v$subject)
    f <- fig(function(m) {
        k <- MODES[[m]]
        ggplot(v, aes(pct, tidytext_reorder(ighv, pct, subject))) +
            geom_col(width = 0.68, fill = k$cat[1]) +
            facet_wrap(~subject, scales = "free_y", ncol = min(np, 3)) +
            scale_y_discrete(labels = function(x) sub("___.*$", "", x)) +
            scale_x_continuous(expand = expansion(c(0, 0.05)), labels = function(x) paste0(x, "%")) +
            labs(x = "% of B cells", y = NULL) + theme_mode(m)
    }, min(10, 3.4 * min(np, 3)), 1.2 + 0.24 * 15 * ceiling(np / 3), "Most used IGHV genes per patient")
    section("ighv", "IGHV gene usage", "The 15 most used heavy-chain V genes per patient. One gene dominating often reflects a single large clone.",
            f, numbers(v %>% transmute(Patient = subject, IGHV = ighv, Cells = n, `%` = round(pct, 1))))
} else ""

# ── 9. clonal threshold ─────────────────────────────────────────────────────────
sec_thr <- if (nrow(thr) || length(opt$dist)) {
    imgs <- paste0(vapply(opt$dist, function(f) sprintf('<figure class="ext"><img src="%s" alt="Distance-to-nearest histogram for %s"><figcaption>%s</figcaption></figure>',
                                                         b64(f), esc(sub("_dist_to_nearest.png$", "", basename(f))),
                                                         esc(sub("_dist_to_nearest.png$", "", basename(f)))), character(1)), collapse = "")
    t2 <- if (nrow(thr)) thr %>% mutate(across(any_of(c("n_cells", "n_clones", "n_expanded_clones", "largest_clone_cells")), as.integer),
                                        threshold = num(threshold)) %>% rename_with(~ gsub("_", " ", .x)) else data.frame()
    section("threshold", "Clonal threshold",
        "Junction distance below which two sequences join the same clone. <i>density</i> / <i>gmm</i>: fitted from the dip between related and unrelated sequences in the patient's data. <i>fallback</i>: too few cells or no clear dip, so the default was used; check the histogram and set <code>--clonal_threshold</code> if needed.",
        table_html(t2), imgs)
} else ""

# ── 10. lineage trees ───────────────────────────────────────────────────────────
sec_lin <- if (nrow(lin) || length(opt$trees)) {
    tr_figs <- ""
    if (length(opt$trees) && requireNamespace("ggtree", quietly = TRUE) && requireNamespace("dowser", quietly = TRUE)) {
        tr_figs <- paste(vapply(opt$trees, function(f) {
            trees <- tryCatch(readRDS(f), error = function(e) NULL)
            if (is.null(trees) || !nrow(trees)) return("")
            top <- head(trees[order(-trees$seqs), ], 6)
            p_id <- sub("_trees.rds$", "", basename(f))
            tips <- if ("sample" %in% names(top$data[[1]]@data)) "sample" else NULL
            smp <- if (!is.null(tips)) unique(unlist(lapply(top$data, function(x) x@data$sample))) else character(0)
            smp <- S$sample[S$sample %in% smp]
            fig(function(m) {
                k <- MODES[[m]]
                cols <- setNames(rep(k$cat, length.out = length(smp)), smp)
                pl <- dowser::plotTrees(top, tips = tips, tipsize = 2.2)
                pl <- lapply(seq_along(pl), function(i) pl[[i]] + ggtitle(sprintf("%s - %d seqs", top$clone_id[i], top$seqs[i])) +
                    theme(plot.background = element_rect(fill = k$surface, colour = NA),
                          panel.background = element_rect(fill = k$surface, colour = NA),
                          text = element_text(colour = k$ink2), plot.title = element_text(colour = k$ink, size = 10),
                          legend.position = if (i == 1) "top" else "none", legend.title = element_blank(),
                          legend.text = element_text(colour = k$ink2, size = 7)))
                if (length(smp)) pl <- lapply(pl, function(g) suppressMessages(
                    g + scale_colour_manual(values = cols, breaks = smp, labels = S$label[match(smp, S$sample)], drop = FALSE)))
                patchwork::wrap_plots(pl, ncol = min(3, length(pl)))
            }, 10, 3.4 * ceiling(length(seq_len(nrow(top))) / 3) + 0.5,
               sprintf("Lineage trees of the largest clones of patient %s", p_id),
               sprintf("%s: the %d largest trees. Tips are cells, coloured by sample; the root is the inferred germline; branch length is mutations per site.",
                       esc(p_id), nrow(top)))
        }, character(1)), collapse = "")
    }
    lt <- if (nrow(lin)) lin %>% mutate(across(any_of(c("n_seqs", "n_cells", "tree_tips")), as.integer),
                                        tree_length = round(num(tree_length), 4)) %>%
        arrange(subject, desc(n_seqs)) else data.frame()
    section("lineage", "Lineage trees",
        "Phylogenies of the larger clones (3+ heavy chains), showing how a clone diversified by somatic hypermutation from its germline.",
        tr_figs, if (nrow(lt)) paste0("<h3>All trees</h3>", table_html(lt, "", c("n_seqs", "n_cells", "tree_tips", "tree_length"))) else "")
} else ""

# ── page ────────────────────────────────────────────────────────────────────────
nav <- c(status = "Sample status", qc = "Cell QC", isotypes = "Isotypes", clones = "Clonal expansion",
         shared = "Shared clones", shm = "SHM", ighv = "IGHV usage", threshold = "Clonal threshold", lineage = "Lineage")
secs <- list(status = sec_status, qc = sec_qc, isotypes = sec_iso, clones = sec_clones, shared = sec_shared,
             shm = sec_shm, ighv = sec_v, threshold = sec_thr, lineage = sec_lin)
present <- names(secs)[nzchar(unlist(secs))]
css <- '
:root{color-scheme:light;--page:#f9f9f7;--surface:#fcfcfb;--ink:#0b0b0b;--ink2:#52514e;--muted:#898781;--line:#e1e0d9;--ring:rgba(11,11,11,.10);--accent:#2a78d6;--good:#0ca30c;--warn:#fab219;--bad:#d03b3b;--goodt:#006300;--warnt:#7a5300;--badt:#b02a2a}
@media (prefers-color-scheme:dark){:root:not([data-theme="light"]){color-scheme:dark;--page:#0d0d0d;--surface:#1a1a19;--ink:#f0efec;--ink2:#c3c2b7;--muted:#898781;--line:#2c2c2a;--ring:rgba(255,255,255,.10);--accent:#3987e5;--goodt:#0ca30c;--warnt:#fab219;--badt:#e66767}}
:root[data-theme="dark"]{color-scheme:dark;--page:#0d0d0d;--surface:#1a1a19;--ink:#f0efec;--ink2:#c3c2b7;--muted:#898781;--line:#2c2c2a;--ring:rgba(255,255,255,.10);--accent:#3987e5;--goodt:#0ca30c;--warnt:#fab219;--badt:#e66767}
*{box-sizing:border-box}body{margin:0;background:var(--page);color:var(--ink);font:15px/1.5 system-ui,-apple-system,"Segoe UI",sans-serif}
main{max-width:1080px;margin:0 auto;padding:24px 16px 64px}header h1{font-size:26px;margin:0 0 4px}header p{margin:0;color:var(--ink2)}
nav{display:flex;flex-wrap:wrap;gap:6px;margin:18px 0 8px}nav a{font-size:13px;color:var(--ink2);text-decoration:none;padding:4px 10px;border:1px solid var(--line);border-radius:999px;background:var(--surface)}nav a:hover{color:var(--accent);border-color:var(--accent)}
.tiles{display:grid;grid-template-columns:repeat(auto-fit,minmax(170px,1fr));gap:12px;margin:18px 0}
.tile{background:var(--surface);border:1px solid var(--ring);border-radius:10px;padding:14px 16px}.tile .v{font-size:28px;font-weight:600}.tile .l{color:var(--ink2);font-size:13px}.tile .s{color:var(--muted);font-size:12px;margin-top:2px}
section{background:var(--surface);border:1px solid var(--ring);border-radius:12px;padding:20px 22px;margin:18px 0}
h2{font-size:19px;margin:0 0 6px}h3{font-size:15px;margin:22px 0 8px}.lead{color:var(--ink2);margin:0 0 14px;max-width:76ch}
.note{color:var(--muted);font-size:12.5px;max-width:90ch}.empty{color:var(--muted);font-style:italic}
figure{margin:12px 0}figure img{max-width:100%;height:auto;display:block}figcaption{color:var(--muted);font-size:12.5px;margin-top:4px}
figure.ext{display:inline-block;max-width:480px;margin-right:12px;vertical-align:top}figure.ext img{border-radius:6px;background:#fff}
.tw{overflow-x:auto;margin:6px 0}table{border-collapse:collapse;font-size:13px;width:100%}th,td{text-align:left;padding:6px 10px;border-bottom:1px solid var(--line);vertical-align:top}
th{color:var(--ink2);font-weight:600;white-space:nowrap}td.n,th.n{text-align:right;font-variant-numeric:tabular-nums}tbody tr:hover{background:color-mix(in srgb,var(--accent) 6%,transparent)}
table.status td:first-child{min-width:190px}td small{display:block;color:var(--muted);font-size:11.5px;max-width:240px;overflow-wrap:anywhere}table.status td:last-child{color:var(--muted);font-size:12px;max-width:340px}
details{margin:6px 0 4px}summary{cursor:pointer;color:var(--accent);font-size:13px}
.b{display:inline-flex;align-items:center;gap:5px;font-size:12px;font-weight:600;white-space:nowrap}.b i{font-style:normal;display:inline-grid;place-items:center;width:16px;height:16px;border-radius:50%;color:#fff;font-size:10px}
.b.good{color:var(--goodt)}.b.good i{background:var(--good)}.b.warn{color:var(--warnt)}.b.warn i{background:var(--warn);color:#0b0b0b}.b.bad{color:var(--badt)}.b.bad i{background:var(--bad)}
code{font-size:12.5px;background:color-mix(in srgb,var(--ink) 7%,transparent);padding:1px 5px;border-radius:4px}footer{color:var(--muted);font-size:12px;margin-top:28px}
'
html <- paste0('<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">',
    '<title>SCRATCH-BCR report</title><style>', css, '</style></head><body><main><header><h1>SCRATCH-BCR report</h1><p>',
    esc(opt$project), " · ", format(Sys.time(), "%Y-%m-%d %H:%M"), if (nzchar(opt$version)) paste0(" · SCRATCH-BCR ", esc(opt$version)) else "",
    '</p></header><nav>', paste0('<a href="#', present, '">', nav[present], "</a>", collapse = ""), "</nav>", tiles,
    paste(unlist(secs[present]), collapse = ""),
    '<footer>Per-cell table: <code>bridge/merged_bcr_object/post_qc_cells.tsv</code> · per-sample status: <code>BCR_QC/tables/sample_status.tsv</code>. Generated by bin/bcr_report.R.</footer></main></body></html>')
writeLines(html, opt$out)
message("[report] wrote ", opt$out, " (", length(present), " sections, ", round(file.size(opt$out) / 1024), " KB)")
