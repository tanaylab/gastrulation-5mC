# Locations of the data and external resources used by the notebooks and scripts.
# Edit these to point at your copies. Notebooks source this file in their first cell
# ("paths.R" from the repo root, "../paths.R" from prep/ and E9.5_multiome.Fig6/).

# Original analysis directory: holds data/, figs/ and revision_e95/data/.
# Notebooks setwd() here, so relative ./data/... and ./scripts/... paths resolve against it.
embdnmt_dir <- "/net/mraid20/ifs/wisdom/tanay_lab/tgdata/users/ronisto/proj/DNMT_final/embdnmt"

# WT gastrulation atlas (spatwt, iteration 10.5.5)
spatwt_dir <- "/net/mraid20/export/tgdata/users/talso/eexplorer_ng/scrna_db/spatwt/iterative/iteration-10.5.5"

# misha genome database (mm10 trackdb)
mm10_trackdb <- "/net/mraid20/export/tgdata/db/tgdb/mm10/trackdb"

# 10x multiome / ATAC pipeline outputs
tg_atac_dir <- "/net/mraid20/export/tgdata/db/tgdb/tg_atac"
multiome_dir <- "/net/mraid20/export/tgdata/users/ronisto/proj/multiome_all"

# Motif models (motif_reg project) and development copies of tanaylab R packages
motif_reg_dir <- "/home/aviezerl/proj/motif_reg"
aviezerl_src_dir <- "/home/aviezerl/src"

# Repository root: the nearest directory, from the working directory up, that holds paths.R.
# Code under scripts/ and E9.5_multiome.Fig6/ is sourced from here.
repo_dir <- local({
    d <- normalizePath(getwd())
    while (!file.exists(file.path(d, "paths.R")) && dirname(d) != d) d <- dirname(d)
    d
})

# Everything the code writes (figures, RDS caches, tables) goes under out_dir.
out_dir <- file.path(repo_dir, "git_test")

# Map a path the code writes to (relative to the working directory, or absolute) onto out_dir,
# keeping its layout relative to embdnmt_dir (or repo_dir), and create its parent directory.
# Absolute paths outside both go to out_dir/external/<original path>; paths already under out_dir are kept.
out_path <- function(p, create = TRUE) {
    roots <- normalizePath(c(embdnmt_dir, repo_dir))
    od <- normalizePath(out_dir, mustWork = FALSE)
    resolve <- function(x) {
        x <- path.expand(x)
        if (!startsWith(x, "/")) x <- file.path(getwd(), x)
        # normalize the deepest existing directory, so different spellings of the same place match
        d <- x
        rest <- character(0)
        while (!dir.exists(d) && dirname(d) != d) {
            rest <- c(basename(d), rest)
            d <- dirname(d)
        }
        do.call(file.path, as.list(c(normalizePath(d), rest[rest != "."])))
    }
    res <- vapply(as.character(p), function(x) {
        x <- resolve(x)
        if (x == od || startsWith(x, paste0(od, "/"))) return(x)
        root <- roots[x == roots | startsWith(x, paste0(roots, "/"))][1]
        rel <- if (is.na(root)) file.path("external", sub("^/", "", x)) else if (x == root) "" else substring(x, nchar(root) + 2)
        f <- if (rel == "") out_dir else file.path(out_dir, rel)
        if (create) dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
        f
    }, character(1), USE.NAMES = FALSE)
    res
}

# Cached intermediates. Notebooks read intermediates (e.g. data/output/mod.RDS, chisq_screen.RDS) from the
# original data directory, as they were when the figures were made. Code in prep/ and the figure notebooks can
# regenerate them, but writes go to git_test/. With use_regenerated_inputs = TRUE, in_path() reads the
# regenerated git_test/ copy of a file when one exists, and the original otherwise.
use_regenerated_inputs <- FALSE
in_path <- function(p) {
    if (!isTRUE(use_regenerated_inputs)) return(p)
    q <- out_path(p, create = FALSE)
    ifelse(file.exists(q), q, p)
}

# Figures. show_fig() always draws a figure in the notebook; it also writes it (through out_path())
# when save_figure is TRUE. Figure notebooks set save_figure in their first cell. `file` may name
# several files (e.g. a .png and a .pdf of the same figure); the figure is drawn once and written to each.
# show = FALSE only saves (used where a loop makes too many figures to draw, e.g. the Fig 4 cluster profiles).
#   show_fig(file, plot, ...)          a ggplot or grid object; saved with ggsave(out_path(file), plot, ...)
#   show_fig(file, draw = { ... }, ...) base-graphics code; to save, it runs again on a pdf()/png()
#                                        device opened with out_path(file) and ...
#   show_fig(file, draw = f(...))      f saves through its own `filename` argument (e.g. pheatmap): drawn with
#                                        filename = NA, saved with filename = out_path(file)
# In Jupyter, the inline figure is a PNG with the proportions and layout of the saved file (same size in
# inches, so fonts and spacing match), at a resolution that keeps its longest side <= fig_inline_px pixels.
# Outside Jupyter it is drawn on the current device.
save_figure <- FALSE
fig_inline_px <- 300   # longest side of inline figures, in pixels
fig_inline_dpi <- 48   # resolution of inline figures that are small enough
# plots drawn directly (not through show_fig) use the repr.plot.width/height a cell sets, at this resolution
if (isTRUE(getOption("jupyter.in_kernel"))) options(repr.plot.res = fig_inline_dpi)
show_fig <- function(file, plot = NULL, draw, show = TRUE, ...) {
    env <- parent.frame()
    dots <- list(...)
    # Inline display: render a PNG at the saved size, then hand it to the notebook.
    display <- getOption("show_fig.display")
    if (!is.function(display) && isTRUE(getOption("jupyter.in_kernel")))
        display <- function(f) IRdisplay::display_png(file = f)
    inline <- show && is.function(display)
    show_png <- function(dims, render) {
        dpi <- min(fig_inline_dpi, fig_inline_px / max(dims))
        tmp <- tempfile(fileext = ".png")
        grDevices::png(tmp, width = dims[1], height = dims[2], units = "in", res = dpi)
        res <- tryCatch(render(), finally = grDevices::dev.off())
        display(tmp)
        unlink(tmp)
        res
    }
    if (!missing(draw)) {
        expr <- substitute(draw)
        fn <- if (is.call(expr)) tryCatch(eval(expr[[1]], env), error = function(e) NULL)
        if (is.function(fn) && "filename" %in% names(formals(fn))) {
            res <- NULL
            if (inline) {
                # the size f() chooses for its file: write it once to a scratch pdf and read its page size
                tmp <- tempfile(fileext = ".pdf")
                expr$filename <- tmp
                eval(expr, env)
                dims <- .pdf_size(tmp)
                unlink(tmp)
                expr$filename <- NA
                if ("silent" %in% names(formals(fn))) expr$silent <- TRUE
                res <- show_png(dims, function() {
                    r <- eval(expr, env)
                    if (!is.null(r$gtable)) { grid::grid.newpage(); grid::grid.draw(r$gtable) }
                    r
                })
                expr$silent <- NULL
            } else if (show) {
                expr$filename <- NA
                res <- eval(expr, env)
            }
            if (isTRUE(save_figure)) for (f in file) {
                expr$filename <- out_path(f)
                eval(expr, env)
            }
            return(invisible(res))
        }
        res <- if (inline) show_png(.dev_size(file[1], dots), function() eval(expr, env)) else if (show) eval(expr, env)
        if (isTRUE(save_figure)) for (f in file) {
            dev <- switch(tolower(tools::file_ext(f)), pdf = grDevices::pdf, png = grDevices::png,
                stop("show_fig: unsupported file type: ", f))
            dev(out_path(f), ...)
            tryCatch(eval(expr, env), finally = grDevices::dev.off())
        }
        return(invisible(res))
    }
    draw_plot <- function() { if (inherits(plot, "ggplot")) print(plot) else { grid::grid.newpage(); grid::grid.draw(plot) } }
    if (inline) show_png(.ggsave_size(dots), draw_plot) else if (show) draw_plot()
    if (isTRUE(save_figure)) for (f in file) ggplot2::ggsave(out_path(f), plot, ...)
    invisible(plot)
}

# Saved sizes, in inches (width, height). Arguments are matched as the saving functions match them
# (partial names such as w = / h = included).
.ggsave_size <- function(dots) {
    a <- do.call(function(width = NA, height = NA, units = "in", dpi = 300, scale = 1, ...)
        list(width = width, height = height, units = units, dpi = dpi, scale = scale), dots)
    to_in <- switch(a$units, `in` = 1, cm = 1 / 2.54, mm = 1 / 25.4, px = 1 / a$dpi)
    d <- c(a$width, a$height) * to_in
    d[is.na(d)] <- 7   # ggsave's default: the size of a new device
    d * a$scale
}
.dev_size <- function(f, dots) {
    if (tolower(tools::file_ext(f)) == "pdf")
        return(do.call(function(width = 7, height = 7, ...) c(width, height), dots))
    a <- do.call(function(width = 480, height = 480, units = "px", pointsize = 12, bg = "white", res = NA, ...)
        list(width = width, height = height, units = units, res = res), dots)
    res <- if (is.na(a$res)) 72 else a$res
    c(a$width, a$height) * switch(a$units, `in` = 1, cm = 1 / 2.54, mm = 1 / 25.4, px = 1 / res)
}
.pdf_size <- function(f) {   # page size of a pdf written by R's pdf() device
    b <- readBin(f, "raw", file.size(f))
    b[b == as.raw(0)] <- as.raw(32)
    m <- regmatches(rawToChar(b), regexpr("/MediaBox *\\[[^]]*\\]", rawToChar(b), useBytes = TRUE))
    box <- as.numeric(strsplit(trimws(gsub("[^0-9. ]", " ", m)), " +")[[1]])
    box[3:4] / 72
}

# Parent project of embdnmt: embryo metadata, the metacell database (scdb/) and the transcriptional
# time inference (my_data/; used by prep/transcriptional_time.ipynb and Fig S2D in ExtData_RNA.ipynb)
dnmt_final_dir <- "/net/mraid20/ifs/wisdom/tanay_lab/tgdata/users/ronisto/proj/DNMT_final"

# mm10 chromosome sizes for wigToBigWig (genome browser tracks, Fig4.ipynb)
mm10_chrom_sizes <- file.path(embdnmt_dir, "hub_p2p_rnp/mm10/mm10.chrom.sizes")
