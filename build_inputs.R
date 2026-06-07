#!/usr/bin/env Rscript

# build_inputs.R — extraction half of LineagePhenotyping, ported from
# CompareDivTime.pl + ComparePositions.pl.
#
# Consumes a freeze directory (per-embryo dats/*.csv copied flat, as produced
# by `embryodb phenotyping freeze`, or assembled by hand) and produces the
# tables run_pipeline.R reads:
#   <name>DivTimeNorm.tsv, <name>CCLengthNorm.tsv,
#   <name>CCLengthMinTerminal.tsv, <name>positions.txt
# plus parity extras (CCLength, DivTime, CCLinNorm, STATS, CellsVsTime,
# radialpositions) that the Perl scripts also emitted.
#
# Usage:
#   Rscript build_inputs.R <freeze_dir> <name> [output_dir] [list_file]
#
# Defaults:
#   output_dir = <freeze_dir>/<name>
#   list_file  = <freeze_dir>/<name>.list   (one series_name per line; order
#                preserved — it sets the column order of the *Norm tables)
#
# Reference data shipped in-repo: SupplementalTable2_DivisionTimes.txt is NOT
# needed here (CompareDivTime.pl reads it but never uses %birth/%division for
# the emitted tables); the Sulston series itself is the regression baseline.
#
# Missing TIME handling: if TIME<series>.csv is absent, frame index is
# converted to minutes via `minutes_per_timepoint` (read from a
# `minutes_per_timepoint:` line in <freeze_dir>/configs/<name>.yaml, or the
# env var BUILD_INPUTS_MPT). Series with neither TIME nor a fallback are
# skipped, matching CompareDivTime.pl.

SULSTON <- "20081128_sulston"

HARDCODED_PARENTS <- c(
  AB = "P0", P1 = "P0", EMS = "P1", P2 = "P1",
  E = "EMS", MS = "EMS", C = "P2", P3 = "P2",
  D = "P3", P4 = "P3", Z2 = "P4", Z3 = "P4"
)

FOUNDERS <- c(
  "ABala", "ABalp", "ABara", "ABarp", "ABpla", "ABplp",
  "ABpra", "ABprp", "E", "MS", "C", "D"
)

# --- helpers ---------------------------------------------------------------

# Perl int(): truncate toward zero. SigFigs truncates to `figs` decimals.
sig_figs <- function(value, figs) {
  scale <- 10^figs
  trunc(value * scale) / scale
}

# Stringify a number the way Perl does for these truncated values: integers
# print without a decimal, fractions print their shortest exact form.
fmt_num <- function(x) {
  # format="fg" + digits=15 matches Perl's %.15g rounding; formatC right-pads
  # to the digit width, so strip it to get Perl's bare-number stringification.
  trimws(formatC(x, format = "fg", digits = 15, drop0trailing = TRUE))
}

# Byte/ASCII order, matching Perl's default sort.
asort <- function(x) sort(x, method = "radix")

read_minutes_per_timepoint <- function(freeze_dir, name) {
  env <- Sys.getenv("BUILD_INPUTS_MPT", "")
  if (nzchar(env)) return(as.numeric(env))
  cfg <- file.path(freeze_dir, "configs", paste0(name, ".yaml"))
  if (file.exists(cfg)) {
    for (ln in readLines(cfg, warn = FALSE)) {
      m <- regmatches(ln, regexec("^\\s*minutes_per_timepoint:\\s*([0-9.]+)", ln))[[1]]
      if (length(m) == 2) return(as.numeric(m[2]))
    }
  }
  NA_real_
}

# --- argument parsing ------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript build_inputs.R <freeze_dir> <name> [output_dir] [list_file]")
}
freeze_dir <- args[1]
name <- args[2]
output_dir <- if (length(args) >= 3) args[3] else file.path(freeze_dir, name)
list_file <- if (length(args) >= 4) args[4] else file.path(freeze_dir, paste0(name, ".list"))

if (!file.exists(list_file)) {
  stop(sprintf("list file not found: %s", list_file))
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

mpt <- read_minutes_per_timepoint(freeze_dir, name)

# Series in list order (drop blanks). Sulston is prepended for the division
# analysis; positions use the list order sorted (see below).
list_series <- trimws(readLines(list_file, warn = FALSE))
list_series <- list_series[nzchar(list_series)]

dats_path <- function(series, kind) {
  # kind: "CD", "ACD", "TIME" -> "<kind><series>.csv"; "AuxInfo" -> "<series>AuxInfo.csv"
  if (kind == "AuxInfo") {
    file.path(freeze_dir, paste0(series, "AuxInfo.csv"))
  } else {
    file.path(freeze_dir, paste0(kind, series, ".csv"))
  }
}

# ===========================================================================
# Part 1 — CompareDivTime.pl port
# ===========================================================================

series_dt <- c(SULSTON, list_series)

# MIN[[series]] / MAX[[series]] are named numeric vectors keyed by cell.
MIN <- list()
MAX <- list()
# parents: named character; seeded with hardcoded founders, extended per cell.
parents <- HARDCODED_PARENTS
# CellsVsTime accumulation (parity extra)
cvt_lines <- c("Series\tMinutes\tTPs\tnCells\tAB\tE\tMS\tC\tD")

for (series in series_dt) {
  cd <- dats_path(series, "CD")
  timefile <- dats_path(series, "TIME")
  if (!file.exists(cd)) {
    message(sprintf("Couldn't find %s - SKIPPING %s", cd, series))
    next
  }
  is_sulston <- identical(series, SULSTON)
  have_time <- file.exists(timefile)
  if (!have_time && !is_sulston && is.na(mpt)) {
    message(sprintf("No TIME for %s and no minutes_per_timepoint - SKIPPING", series))
    next
  }

  # timelookup: frame index -> minutes
  timelookup <- NULL
  if (have_time) {
    tl <- readLines(timefile, warn = FALSE)
    tp <- integer(0); mins <- numeric(0)
    for (ln in tl) {
      f <- strsplit(ln, "\t", fixed = TRUE)[[1]]
      if (length(f) < 3) next
      this_tp <- as.integer(f[2])
      seconds <- as.numeric(f[3])
      tp <- c(tp, this_tp)
      mins <- c(mins, trunc(seconds / 6) / 10)
    }
    timelookup <- setNames(mins, as.character(tp))
  }

  cd_lines <- readLines(cd, warn = FALSE)
  cellmin <- numeric(0); cellmax <- numeric(0)
  # per-time lineage counts for CellsVsTime
  cnt_total <- list(); cnt_AB <- list(); cnt_E <- list()
  cnt_MS <- list(); cnt_C <- list(); cnt_D <- list()
  bump <- function(lst, key) {
    k <- as.character(key)
    lst[[k]] <- if (is.null(lst[[k]])) 1L else lst[[k]] + 1L
    lst
  }

  for (ln in cd_lines) {
    f <- strsplit(ln, ",", fixed = TRUE)[[1]]
    if (length(f) < 3) next
    cell <- f[2]
    if (cell == "cell") next
    time <- f[3]
    timenum <- suppressWarnings(as.numeric(time))
    if (!is_sulston) {
      key <- as.character(time)
      if (have_time) {
        if (!is.na(timelookup[key]) && key %in% names(timelookup)) {
          timenum <- unname(timelookup[key])
        }
      } else {
        # minutes_per_timepoint fallback (frame index -> minutes)
        timenum <- timenum * mpt
      }
    }

    # lineage counts vs time
    if (grepl("^AB", cell)) cnt_AB <- bump(cnt_AB, timenum)
    else if (grepl("^E", cell) && cell != "EMS") cnt_E <- bump(cnt_E, timenum)
    else if (grepl("^MS", cell)) cnt_MS <- bump(cnt_MS, timenum)
    else if (grepl("^C", cell)) cnt_C <- bump(cnt_C, timenum)
    else if (grepl("^D", cell)) cnt_D <- bump(cnt_D, timenum)
    cnt_total <- bump(cnt_total, timenum)

    # standard-naming parent for novel cells: chop last char (Perl: set
    # $parents{$cell} unless already defined).
    if (!cell %in% names(parents)) {
      parents[cell] <- substr(cell, 1, nchar(cell) - 1)
    }

    # min/max observed time
    cur_min <- cellmin[cell]
    if (is.na(cur_min) || cur_min > timenum) cellmin[cell] <- timenum
    cur_max <- cellmax[cell]
    if (is.na(cur_max) || cur_max < timenum) cellmax[cell] <- timenum
  }

  MIN[[series]] <- cellmin
  MAX[[series]] <- cellmax

  # CellsVsTime output rows (numeric time order)
  times <- asort_num <- sort(as.numeric(names(cnt_total)))
  tp_i <- 1L
  for (t in times) {
    if (!grepl("[0-9]", as.character(t))) next
    k <- as.character(t)
    g <- function(lst) { v <- lst[[k]]; if (is.null(v)) 0L else v }
    cvt_lines <- c(cvt_lines, paste(
      series, fmt_num(t), tp_i, g(cnt_total), g(cnt_AB), g(cnt_E),
      g(cnt_MS), g(cnt_C), g(cnt_D), sep = "\t"
    ))
    tp_i <- tp_i + 1L
  }
}

# progeny: for cell in sorted(parents keys): progeny[parents[cell]] = cell
progeny <- character(0)
for (cell in asort(names(parents))) {
  progeny[parents[[cell]]] <- cell
}

# Cell universe = every observed cell plus every cell that appears as a
# parent. The parent values pull in never-observed ancestors (notably the
# root P0), matching CompareDivTime.pl which iterates keys+values of %parents.
all_cells <- asort(unique(c(
  unlist(lapply(MIN, names), use.names = FALSE),
  names(parents),
  unname(parents)
)))

has_in <- function(map, series, cell) {
  v <- map[[series]]
  !is.null(v) && cell %in% names(v) && !is.na(v[cell])
}
has_parent <- function(cell, series) {
  if (!cell %in% names(parents)) return(FALSE)
  has_in(MIN, series, parents[[cell]])
}
has_progeny <- function(cell, series) {
  if (!cell %in% names(progeny)) return(FALSE)
  has_in(MIN, series, progeny[[cell]])
}
getv <- function(map, series, cell) {
  v <- map[[series]][cell]
  unname(v)
}

# Least-squares fit matching Statistics::Descriptive (plain-double sums,
# left-to-right) so truncated outputs match byte-for-byte.
lsq <- function(xs, ys) {
  n <- length(ys)
  if (n < 2) return(c(q = 0, m = 0))
  sx <- 0; sxx <- 0; sxy <- 0; sy <- 0
  for (i in seq_len(n)) {
    sy <- sy + ys[i]
  }
  for (i in seq_len(n)) {
    x <- xs[i]; y <- ys[i]
    sxx <- sxx + x * x
    sxy <- sxy + x * y
    sx <- sx + x
  }
  denom <- n * sxx - sx * sx
  if (abs(denom) <= 0) return(c(q = 0, m = 0))
  m <- (n * sxy - sx * sy) / denom
  q <- (sxx * sy - sx * sxy) / denom
  c(q = q, m = m)
}

# slope/offset/lineageSlopes per series
slopes <- setNames(numeric(length(series_dt)), series_dt)
offsets <- setNames(numeric(length(series_dt)), series_dt)
lineage_slopes <- list()  # lineage_slopes[[series]][[founder]]

for (series in series_dt) {
  series_times <- numeric(0)
  sulston_times <- numeric(0)
  lin_bins <- list(); sul_bins <- list()
  for (cell in all_cells) {
    if (has_parent(cell, series) && has_progeny(cell, series) &&
        has_progeny(cell, SULSTON)) {
      series_times <- c(series_times, getv(MAX, series, cell))
      sulston_times <- c(sulston_times, getv(MAX, SULSTON, cell))
      for (f in FOUNDERS) {
        if (grepl(f, cell, fixed = TRUE)) {
          lin_bins[[f]] <- c(lin_bins[[f]], getv(MAX, series, cell))
          sul_bins[[f]] <- c(sul_bins[[f]], getv(MAX, SULSTON, cell))
        }
      }
    }
  }
  if (length(series_times) == 0) next
  fit <- lsq(sulston_times, series_times)
  slopes[series] <- fit["m"]
  offsets[series] <- fit["q"]
  ls <- list()
  for (f in FOUNDERS) {
    if (is.null(lin_bins[[f]]) || length(lin_bins[[f]]) < 2) {
      ls[[f]] <- fit["m"]
    } else {
      lf <- lsq(sul_bins[[f]], lin_bins[[f]])
      ls[[f]] <- if (lf["m"] == 0) fit["m"] else lf["m"]
    }
  }
  lineage_slopes[[series]] <- ls
}

# --- emit division tables --------------------------------------------------

tsv_header <- paste0("Cell\t", paste(series_dt, collapse = "\t"))

cc_lines <- tsv_header
divtime_lines <- tsv_header
ccnorm_lines <- tsv_header
divtimenorm_lines <- tsv_header
ccterminal_lines <- tsv_header
cclinnorm_lines <- tsv_header

NA_STR <- "NA"

for (cell in all_cells) {
  cc_row <- cell; dt_row <- cell
  ccn_row <- cell; dtn_row <- cell; cct_row <- cell; ccl_row <- cell
  for (series in series_dt) {
    hp <- has_parent(cell, series)
    hpr <- has_progeny(cell, series)
    mx <- getv(MAX, series, cell); mn <- getv(MIN, series, cell)
    # Perl reads an undefined MAX as 0 in numeric context; this only bites
    # never-observed ancestors that still have observed progeny (e.g. the
    # zygote P0, whose division time is 0). CC length keeps NA (needs MIN).
    mxd <- if (is.na(mx)) 0 else mx
    cc_raw <- if (is.na(mx)) NA else (mx - mn)
    slope <- slopes[series]
    cc <- if (!is.na(slope) && slope != 0 && !is.na(cc_raw)) cc_raw / slope else 0

    # CCLength.tsv (raw cc, terminal -> NA)
    if (!hpr || !hp) {
      cc_row <- paste0(cc_row, "\t", NA_STR)
    } else {
      cc_row <- paste0(cc_row, "\t", fmt_num(sig_figs(cc_raw, 1)))
    }
    # DivTime.tsv (raw max, terminal -> NA)
    if (!hpr) {
      dt_row <- paste0(dt_row, "\t", NA_STR)
    } else {
      dt_row <- paste0(dt_row, "\t", fmt_num(sig_figs(mxd, 1)))
    }

    # CCLengthNorm + CCLengthMinTerminal + CCLinNorm
    if (!hp || !hpr) {
      ccn_row <- paste0(ccn_row, "\t", NA_STR)
      if (!is.na(cc) && cc > 0) {
        cct_row <- paste0(cct_row, "\t", fmt_num(sig_figs(cc, 1)))
      } else {
        cct_row <- paste0(cct_row, "\t", NA_STR)
      }
      ccl_row <- paste0(ccl_row, "\t", NA_STR)
    } else {
      ccn_row <- paste0(ccn_row, "\t", fmt_num(sig_figs(cc, 1)))
      cct_row <- paste0(cct_row, "\t", NA_STR)
      printed <- FALSE
      ls <- lineage_slopes[[series]]
      for (f in FOUNDERS) {
        if (grepl(f, cell, fixed = TRUE)) {
          printed <- TRUE
          lslope <- if (!is.null(ls)) ls[[f]] else NA
          if (is.na(lslope) || lslope <= 0) lslope <- 1
          ccl_row <- paste0(ccl_row, "\t", fmt_num(sig_figs((mx - mn) / lslope, 1)))
        }
      }
      if (!printed) {
        ccl_row <- paste0(ccl_row, "\t", fmt_num(sig_figs(cc, 1)))
      }
    }

    # DivTimeNorm
    thisslope <- if (is.na(slope) || slope == 0) 1 else slope
    if (!hpr) {
      dtn_row <- paste0(dtn_row, "\t", NA_STR)
    } else {
      dtn_row <- paste0(dtn_row, "\t",
                        fmt_num(sig_figs((mxd - offsets[series]) / thisslope, 1)))
    }
  }
  cc_lines <- c(cc_lines, cc_row)
  divtime_lines <- c(divtime_lines, dt_row)
  ccnorm_lines <- c(ccnorm_lines, ccn_row)
  divtimenorm_lines <- c(divtimenorm_lines, dtn_row)
  ccterminal_lines <- c(ccterminal_lines, cct_row)
  cclinnorm_lines <- c(cclinnorm_lines, ccl_row)
}

# STATS.tsv
stats_lines <- paste0("series\toffset\tslope\t", paste(FOUNDERS, collapse = "\t"))
for (series in series_dt) {
  ls <- lineage_slopes[[series]]
  founder_vals <- vapply(FOUNDERS, function(f) {
    v <- if (!is.null(ls)) ls[[f]] else slopes[series]
    fmt_num(sig_figs(v, 3))
  }, character(1))
  stats_lines <- c(stats_lines, paste(
    series, fmt_num(sig_figs(offsets[series], 1)),
    fmt_num(sig_figs(slopes[series], 3)),
    paste(founder_vals, collapse = "\t"), sep = "\t"
  ))
}

w <- function(lines, suffix) {
  writeLines(lines, file.path(output_dir, paste0(name, suffix)))
}
w(cvt_lines, "CellsVsTime.tsv")
w(cc_lines, "CCLength.tsv")
w(divtime_lines, "DivTime.tsv")
w(ccnorm_lines, "CCLengthNorm.tsv")
w(divtimenorm_lines, "DivTimeNorm.tsv")
w(ccterminal_lines, "CCLengthMinTerminal.tsv")
w(cclinnorm_lines, "CCLinNorm.tsv")
w(stats_lines, "STATS.tsv")

# ===========================================================================
# Part 2 — ComparePositions.pl port
# ===========================================================================

# positions use the list series, sorted byte order (no sulston).
series_pos <- asort(list_series)

# X/Y/Z[[cell]][[series]] keyed by time; centroids per series.
PX <- new.env(); PY <- new.env(); PZ <- new.env()
xavg <- setNames(numeric(0), character(0))
yavg <- xavg; zavg <- xavg
pos_series_seen <- character(0)
# Per-cell observed times, accumulated as chunks and split once at the end to
# avoid O(k^2) vector growth.
cell_chunks <- list(); time_chunks <- list()

store <- function(envv, cell, series, time, val) {
  key <- paste0(cell, "\r", series, "\r", time)
  assign(key, val, envir = envv)
}
fetch <- function(envv, cell, series, time) {
  key <- paste0(cell, "\r", series, "\r", time)
  if (exists(key, envir = envv, inherits = FALSE)) get(key, envir = envv) else NA
}

for (series in series_pos) {
  acd <- dats_path(series, "ACD")
  if (!file.exists(acd)) {
    message(sprintf("No ACD for %s - skipping positions", series))
    next
  }
  pos_series_seen <- c(pos_series_seen, series)
  acd_lines <- readLines(acd, warn = FALSE)
  n <- length(acd_lines)
  cellv <- character(n); timev <- numeric(n); k <- 0L
  xsum <- 0; ysum <- 0; zsum <- 0; avgcount <- 0
  for (ln in acd_lines) {
    f <- strsplit(ln, ",", fixed = TRUE)[[1]]
    if (length(f) < 11) next
    cell <- f[2]; t <- f[3]
    x <- suppressWarnings(as.numeric(f[10]))
    y <- suppressWarnings(as.numeric(f[11]))
    z <- suppressWarnings(as.numeric(f[9]))
    if (is.na(x) || x == 0) next
    # ACD already in microns; xcor=ycor=zcor=1, no rotation.
    store(PX, cell, series, t, x)
    store(PY, cell, series, t, y)
    store(PZ, cell, series, t, z)
    k <- k + 1L
    cellv[k] <- cell; timev[k] <- as.numeric(t)
    xsum <- xsum + x; ysum <- ysum + y; zsum <- zsum + z
    avgcount <- avgcount + 1
  }
  if (k > 0L) {
    cell_chunks[[series]] <- cellv[seq_len(k)]
    time_chunks[[series]] <- timev[seq_len(k)]
  }
  if (avgcount > 0) {
    xavg[series] <- xsum / avgcount
    yavg[series] <- ysum / avgcount
    zavg[series] <- zsum / avgcount
  } else {
    message(sprintf("NO DATA FOR %s", series))
  }
}

# cells_times[[cell]] = numeric vector of all observed times for that cell.
cells_times <- split(unlist(time_chunks, use.names = FALSE),
                     unlist(cell_chunks, use.names = FALSE))

pos_header <- "cell\ttime"
radial_header <- "cell\ttime"
for (series in pos_series_seen) {
  pos_header <- paste0(pos_header, "\t", paste0("X_", series), "\t",
                       paste0("Y_", series), "\t", paste0("Z_", series))
  radial_header <- paste0(radial_header, "\t", paste0("X_", series), "\t",
                          paste0("R_", series), "\t", paste0("THETA_", series))
}

PI <- 4 * atan2(1, 1)
all_pos_cells <- asort(names(cells_times))
# Accumulate per-cell row blocks in a list, collapse once, to avoid the
# O(n^2) cost of growing a vector with c() across ~100k+ rows.
pos_acc <- vector("list", length(all_pos_cells))
radial_acc <- vector("list", length(all_pos_cells))
for (ci in seq_along(all_pos_cells)) {
  cell <- all_pos_cells[ci]
  times <- sort(unique(cells_times[[cell]]))
  prow_vec <- character(length(times))
  rrow_vec <- character(length(times))
  for (ti in seq_along(times)) {
    time <- times[ti]
    tstr <- fmt_num(time)
    prow <- paste0(cell, "\t", tstr)
    rrow <- paste0(cell, "\t", tstr)
    for (series in pos_series_seen) {
      xv <- fetch(PX, cell, series, tstr)
      if (!is.na(xv)) {
        xnew <- sig_figs(xv - xavg[series], 1)
        ynew <- sig_figs(fetch(PY, cell, series, tstr) - yavg[series], 1)
        znew <- sig_figs(fetch(PZ, cell, series, tstr) - zavg[series], 1)
        r <- sqrt(ynew^2 + znew^2)
        if (ynew > 0) {
          theta <- atan(znew / ynew)
        } else if (ynew < 0) {
          theta <- if (znew >= 0) atan(znew / ynew) + PI else atan(znew / ynew) - PI
        } else {
          theta <- if (znew > 0) PI / 2 else if (znew < 0) -PI / 2 else 0
        }
        prow <- paste0(prow, "\t", fmt_num(unname(xnew)), "\t",
                       fmt_num(unname(ynew)), "\t", fmt_num(unname(znew)))
        rrow <- paste0(rrow, "\t", fmt_num(unname(xnew)), "\t",
                       fmt_num(sig_figs(r, 1)), "\t", fmt_num(sig_figs(theta, 3)))
      } else {
        prow <- paste0(prow, "\t\t\t")
        rrow <- paste0(rrow, "\t\t\t")
      }
    }
    prow_vec[ti] <- prow
    rrow_vec[ti] <- rrow
  }
  pos_acc[[ci]] <- prow_vec
  radial_acc[[ci]] <- rrow_vec
}
pos_lines <- c(pos_header, unlist(pos_acc, use.names = FALSE))
radial_lines <- c(radial_header, unlist(radial_acc, use.names = FALSE))

writeLines(pos_lines, file.path(output_dir, paste0(name, "positions.txt")))
writeLines(radial_lines, file.path(output_dir, paste0(name, "radialpositions.txt")))

message(sprintf("build_inputs.R: wrote outputs to %s", output_dir))
