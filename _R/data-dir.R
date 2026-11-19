#'
#' Folder for cluster outputs needed by local analyses (owner's rule,
#' 2026-10-30): `/Volumes/lan_images4/tany_time` if that drive is
#' connected, otherwise `/Volumes/LAN_4TB_SSD/tany_time`. When both are
#' connected they're kept in sync, so either one has everything.
#'
#' Usage: `source("_R/data-dir.R")`, then `tany_dir("kraken")`, etc.
#'

tany_dir <- function(...) {
    roots <- c("/Volumes/lan_images4/tany_time",
               "/Volumes/LAN_4TB_SSD/tany_time")
    roots <- roots[dir.exists(roots)]
    if (length(roots) == 0) {
        stop("Connect drive lan_images4 or LAN_4TB_SSD (needs a tany_time ",
             "folder)")
    }
    file.path(roots[[1]], ...)
}

#' The first 2026-11 BayPass rerun (18 pools, with SN-1-95) was moved to
#' `archive-18pools/` on both drives (owner, 2026-11-19). Scripts that
#' analyse it set `ARCHIVE_18POOLS <- TRUE` before sourcing this file, so
#' their `baypass-*` folders resolve there.
if (exists("ARCHIVE_18POOLS") && isTRUE(ARCHIVE_18POOLS)) {
    tany_dir_current <- tany_dir
    tany_dir <- function(...) {
        a <- c(...)
        if (length(a) > 0 && grepl("^baypass-", a[[1]])) {
            tany_dir_current("archive-18pools", ...)
        } else {
            tany_dir_current(...)
        }
    }
}
