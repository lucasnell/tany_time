#'
#' Fetch mitogenomes from NCBI for the other species in the competitive-
#' mapping reference (`_bash/00-references/02-hologenome`) whose assemblies
#' lack one (`_bash/tests/other-mito-check.sh`). Run from the repo root.
#'
#' Chosen 2026-10-28 from nuccore records with `mitochondrion[filter]` and
#' 10-30 kb (search results in `_data/other-mitogenomes.csv`):
#' prefer a record from the same BioProject as the nuclear assembly, then
#' RefSeq (`NC_`), then the longest.
#'   Cripar  Chironomus riparius: OU815687.1 (BioProject 849090, as the
#'           assembly GCA_917627325; another record, LS991863.1, is from
#'           an older project)
#'   Saterr  Smittia aterrima: NC_069312.1 (RefSeq copy of ON838254.1; the
#'           assembly GCA_033063855.1 has no BioProject mitogenome)
#' The C. tentans (2023) and P. akamusi assemblies already contain theirs
#' (OY754306.1, JAFFZY010000103.1); `02-hologenome` relabels those.
#'
#' Output: `_bash/references/other-mitogenomes.fasta`, sequences named
#' `mt_<code>__<accession>`; `_data/other-mitogenomes.csv`.
#'

library(tidyverse)
library(rentrez)

species <- c(Cripar = "Chironomus riparius", Saterr = "Smittia aterrima")
chosen <- c(Cripar = "OU815687.1", Saterr = "NC_069312.1")

hits <- imap(species, \(sp, code) {
    term <- paste0("\"", sp, "\"[Organism] AND mitochondrion[filter] AND ",
                   "10000:30000[SLEN]")
    s <- entrez_search(db = "nuccore", term = term, retmax = 50)
    Sys.sleep(0.4)
    if (s$count == 0) return(NULL)
    su <- entrez_summary(db = "nuccore", id = s$ids)
    if (!is.null(su$uid)) su <- list(su)
    map(su, \(x) tibble(code = code, accession = x$accessionversion,
                        length = x$slen, bioproject = x$projectid,
                        title = x$title)) |>
        list_rbind()
}) |>
    list_rbind() |>
    mutate(chosen = accession %in% chosen)

write_csv(hits, "_data/other-mitogenomes.csv")
print(hits)
stopifnot(all(chosen %in% hits$accession))

fa <- imap_chr(chosen, \(acc, code) {
    x <- entrez_fetch(db = "nuccore", id = acc, rettype = "fasta")
    Sys.sleep(0.4)
    sq <- str_remove_all(str_remove(x, "^>[^\n]*\n"), "\\s")
    paste0(">mt_", code, "__", acc, "\n",
           paste(str_extract_all(sq, ".{1,80}")[[1]], collapse = "\n"))
})
writeLines(fa, "_bash/references/other-mitogenomes.fasta")
