
#'
#' Fetch reference sequences from NCBI for the mitogenome assembly
#' (`_bash/00-references/01-mitogenome`) and the barcode
#' contamination index (`00-references/02-hologenome`). Run from the repo root.
#'
#' Outputs (in `_bash/references/`, transferred to cluster jobs):
#'   tgraci-mito-bait.fasta  all T. gracilentus mitochondrial records, used
#'                           to recruit mitochondrial nanopore reads
#'   tgraci-coi.fasta        T. gracilentus COI barcodes, to validate the
#'                           assembled mitogenome
#'   coi-barcodes.fasta      one COI barcode per abundant trap taxon
#'                           (Einarsson et al. 2026, Table 1), named
#'                           `bc_<code>__<accession>`
#' and `_data/coi-barcode-selection.csv` (chosen records and candidates).
#'
#' Barcode choice per taxon: COI records 600-700 bp; prefer specimens from
#' Iceland (country or geo_loc_name); then the longest; then the lowest
#' accession.
#'
#' Procladius islandicus has no sequences in NCBI (checked 2026-10-12), so
#' Procladius is represented at genus level (`bc_Procla`): one barcode
#' (same rule) for each Procladius species with >= 4 COI barcodes
#' (600-700 bp) from Norway, Finland, or Sweden.
#'

library(tidyverse)
library(rentrez)

out_dir <- "_bash/references"
dir.create(out_dir, showWarnings = FALSE)

taxa <- tribble(
    ~code,     ~species,
    "Mlindr",  "Micropsectra lindrothi",
    "Matrof",  "Micropsectra atrofasciata",
    "Csylve",  "Cricotopus sylvestris",
    "Cislan",  "Chironomus islandicus",
    "Eminor",  "Eukiefferiella minor",
    "Svitta",  "Simulium vittatum")

coi_terms <- paste0("(COI[Gene] OR COX1[Gene] OR CO1[Gene] OR ",
                    "\"cytochrome oxidase subunit 1\"[Title] OR ",
                    "\"cytochrome c oxidase subunit I\"[Title])")

#' Parse GenBank flat-file text into a tibble (accession, length, country,
#' sequence) for each record.
parse_gb <- function(txt) {
    recs <- str_split(txt, "\n//\n?")[[1]]
    recs <- recs[str_detect(recs, "^LOCUS|\nLOCUS")]
    map(recs, \(r) {
        acc <- str_match(r, "\nVERSION\\s+(\\S+)")[, 2]
        country <- str_match(r, "/(?:country|geo_loc_name)=\"([^\"]+)\"")[, 2]
        sq <- str_match(r, "ORIGIN([\\s\\S]*)$")[, 2] |>
            str_remove_all("[^acgtnACGTNrykmswbdhvRYKMSWBDHV]") |>
            toupper()
        tibble(accession = acc, length = nchar(sq), country = country,
               sequence = sq)
    }) |>
        list_rbind()
}

fetch_gb <- function(term, retmax = 300) {
    ids <- entrez_search(db = "nuccore", term = term, retmax = retmax)$ids
    if (length(ids) == 0) return(tibble())
    # Fetch in batches to keep requests small:
    split(ids, ceiling(seq_along(ids) / 100)) |>
        map(\(x) {
            Sys.sleep(0.4)
            parse_gb(entrez_fetch(db = "nuccore", id = x, rettype = "gb",
                                  retmode = "text"))
        }) |>
        list_rbind()
}

write_fasta <- function(names, seqs, file) {
    lines <- map2(names, seqs, \(n, s) c(paste0(">", n),
                                         str_sub(s, seq(1, nchar(s), 80),
                                                 pmin(seq(1, nchar(s), 80) + 79,
                                                      nchar(s)))))
    write_lines(unlist(lines), file)
}


# =============================================================================*
# T. gracilentus mitochondrial bait and COI validation ----
# =============================================================================*

tg_mito <- fetch_gb("\"Tanytarsus gracilentus\"[Organism] AND mitochondrion[Filter]")
stopifnot(nrow(tg_mito) > 0)
write_fasta(tg_mito$accession, tg_mito$sequence,
            file.path(out_dir, "tgraci-mito-bait.fasta"))

tg_coi <- tg_mito |> filter(accession %in% c("MZ659759.1", "MN679642.1"))
stopifnot(nrow(tg_coi) == 2)
write_fasta(tg_coi$accession, tg_coi$sequence,
            file.path(out_dir, "tgraci-coi.fasta"))


# =============================================================================*
# COI barcodes for trap taxa ----
# =============================================================================*

cands <- taxa |>
    mutate(recs = map(species, \(s) {
        fetch_gb(sprintf("\"%s\"[Organism] AND %s", s, coi_terms))
    })) |>
    unnest(recs) |>
    filter(length >= 600, length <= 700) |>
    mutate(iceland = !is.na(country) & str_detect(country, "Iceland"))

chosen <- cands |>
    arrange(code, desc(iceland), desc(length), accession) |>
    group_by(code) |>
    mutate(n_candidates = n(), n_iceland = sum(iceland)) |>
    slice(1) |>
    ungroup()

#' Procladius at genus level (see top). Searched by country so the number
#' of records stays manageable.
procla <- c("Norway", "Finland", "Sweden") |>
    map(\(cn) {
        fetch_gb(sprintf(paste("Procladius[Organism] AND %s AND",
                               "600:700[SLEN] AND %s[All Fields]"),
                         coi_terms, cn), retmax = 500)
    }) |>
    list_rbind() |>
    distinct(accession, .keep_all = TRUE) |>
    filter(length >= 600, length <= 700,
           str_detect(country, "Norway|Finland|Sweden"))
procla_org <- entrez_summary(db = "nuccore",
                             id = procla$accession)
procla$species <- map_chr(procla_org, \(x) x$organism)
procla_chosen <- procla |>
    group_by(species) |>
    filter(n() >= 4) |>
    arrange(desc(length), accession) |>
    mutate(n_candidates = n(), n_iceland = 0L, iceland = FALSE) |>
    slice(1) |>
    ungroup() |>
    mutate(code = "Procla")
chosen <- bind_rows(chosen, procla_chosen)

missing <- setdiff(taxa$code, chosen$code)
if (length(missing) > 0) warning("No barcode for: ", paste(missing, collapse = ", "))

write_fasta(paste0("bc_", chosen$code, "__", chosen$accession),
            chosen$sequence, file.path(out_dir, "coi-barcodes.fasta"))

chosen |>
    select(code, species, accession, length, country, iceland,
           n_candidates, n_iceland) |>
    write_csv("_data/coi-barcode-selection.csv")

print(select(chosen, code, species, accession, length, country, n_candidates,
             n_iceland))
cat(sprintf("T. gracilentus mitochondrial bait: %i records, %i bp total\n",
            nrow(tg_mito), sum(tg_mito$length)))
