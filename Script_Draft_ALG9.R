#working directory
setwd("~/Desktop/Genomenon")

#csv imports
evi <- read.csv("evidence-7067-2026-02-04.csv", header = TRUE)
qa <- read.csv("qa-7067-2026-02-04.csv", header = TRUE)

#install packages
install.packages("ggplot2")
install.packages("dplyr")
install.packages("readr")
install.packages("forcats")
# Install packages and grab the reference sequence
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("GenomicFeatures", "GenomicRanges"))
BiocManager::install("txdbmaker")
BiocManager::install("GenomeInfoDbData")
install.packages("RMariaDB")

library(readr)
library(ggplot2)
library(stringr)
library(forcats)
library(GenomicFeatures)
library(GenomicRanges)
library(GenomeInfoDbData)
library(txdbmaker)
library(jsonlite)
library(RMariaDB)
library(dplyr)

#resolve package conflicts
conflicted::conflicts_prefer(dplyr::select)
conflicted::conflicts_prefer(dplyr::slice)
conflicted::conflicts_prefer(dplyr::first)
conflicted::conflicts_prefer(dplyr::lag)


# just to make sure that the files were read in correctly. "glimpse" function is kind of similar to print.
glimpse(evi)
glimpse(qa)

# concatenate the evidence csv file and the qa csv file to add in the variant_type (ie, missense, nonsense, frameshift, etc) and the cdna if available.
# --- Canonical transcript we care about for cDNA strings ---
target_tx <- "NM_024740.2"

# --- Build a QA lookup that picks the canonical RefSeq cDNA (pipe-delimited) ---
qa_cdna_best <- qa %>%
  dplyr::select(gene, variant_desc, variant_type, cdna_variant) %>%
  mutate(cdna_variant = na_if(cdna_variant, "")) %>%
  tidyr::separate_rows(cdna_variant, sep = "\\|") %>%
  mutate(
    cdna_variant = str_trim(cdna_variant),
    tx = str_extract(cdna_variant, "^NM_\\d+\\.\\d+"),
    hgvs_c = str_replace(cdna_variant, "^NM_\\d+\\.\\d+:", "")  # drop transcript prefix
  ) %>%
  group_by(gene, variant_desc) %>%
  arrange(desc(tx == target_tx), tx, hgvs_c) %>%   # canonical first; deterministic fallback
  slice(1) %>%
  ungroup() %>%
  select(gene, variant_desc, variant_type, hgvs_c)

# --- Join evidence + QA (adds variant_type and canonical hgvs_c) ---
evi2 <- evi %>%
  left_join(
    qa_cdna_best,
    by = c("gene_symbol" = "gene", "mutation" = "variant_desc")
  )

# clean the dataframe so we just have the columns and info we want. 
evi_clean <- evi2 %>%
  select(
    gene_symbol,
    mutation,
    variant_type,
    classification,
    hgvs_c,
    all_cats,
    evidence_type,
    unweighted,
    number_of_probands
  )

# Creates a dataframe with PS3 data. The weird symbols around the PS3 are a regex incase in the future we have a modification to ps3 we may like to pull, like "PS3_M" which is a downgrade to ps3 moderate.
evi_ps3 <- evi_clean %>%
  mutate(
    ps3_applied_row =
      !is.na(all_cats) & str_detect(all_cats, "(^|\\|)ps3(\\||$|_)")
  )

# make a dataframe that contains information about which variants have weighted probands (without losing variants that may have no probands at all)
evi_weighted <- evi_ps3 %>%
  mutate(
    is_weighted_proband =
      evidence_type == "unrelated_proband" & unweighted == FALSE
  )

# Creates dataframe without duplicates. It should look like:
# One row per variant, with:
# variant identity (gene_symbol, mutation)
# variant_type
# classification
# hgvs_c (cdna)
# weighted_probands (numeric)
# any_weighted_proband (TRUE/FALSE)
# ps3_applied (TRUE/FALSE)
# n_evidence_rows (quick qc check)

variant_summary <- evi_weighted %>% 
  group_by(gene_symbol, mutation) %>% # treats all evidence rows for the same variant (as there is multiple) as belonging together.
  summarise( # make a summary for each variant
    variant_type = first(variant_type), # carry on the variant type from our og file
    classification = first(classification), # carry on the classification from our og file
    
    hgvs_c = first(hgvs_c), # adds in cdna description if available
    weighted_probands = sum(number_of_probands[is_weighted_proband], na.rm = TRUE), # adds up all the weighted probands
    any_weighted_proband = any(is_weighted_proband, na.rm = TRUE), # is there any weighted proband evidence for this variant?
    
    ps3_applied = any(ps3_applied_row, na.rm = TRUE), # from our all_cats step. basically asks if there was ps3 applied for the variant.
    
    n_evidence_rows = n(), # how many evidence rows were for the variant. just a quick qc check
    .groups = "drop"
  )

# gets reference genome info from ucsc
txdb_refgene_hg19 <- txdbmaker::makeTxDbFromUCSC(
  genome = "hg19",
  tablename = "refGene"
)

target_tx <- "NM_024740.2" # This is the canonical refseq transcript we care about for our test gene (ALG9).
target_tx_novers <- sub("\\.[0-9]+$", "", target_tx) # just in case the .4 version isnt in the ucsc table, this will at least grab the correct id (NM_024740)

# Get all transcripts so we can confirm the name is present. Basically, gives us every refseq transcript in hg19.
tx_all <- transcripts(txdb_refgene_hg19) 

# Find the canonical transcript (with or without version). Basically, from all the transcripts available, find and keep the one that matches my canonical transcript.
tx_canonical <- tx_all[tx_all$tx_name %in% c(target_tx, target_tx_novers)]

# if the canonical transcript isn't available, pls tell us
if (length(tx_canonical) == 0) {
  stop("Couldn't find ", target_tx, " (or ", target_tx_novers, ") in this TxDb. Check tx_name field.")
}

# Prefer the versioned accession if both exist
if (length(tx_canonical) > 1) {
  v <- tx_canonical[tx_canonical$tx_name == target_tx]
  tx_canonical <- if (length(v) > 0) v else tx_canonical[1]
}

tx_name <- tx_canonical$tx_name[1] # save the name of the canonical transcript
tx_name # prints the name of the transcript just so were sure it grabbed the right one.

# Lets get the exon and intron boundaries
ex_by_tx <- exonsBy(txdb_refgene_hg19, by = "tx", use.names = TRUE) # the "exonsBy() function grabs exon coordinates
in_by_tx <- intronsByTranscript(txdb_refgene_hg19, use.names = TRUE) # the intronsByTranscript () function computes introns as the genomic gaps between exons. basically, infer the coordinates between the exons.

exons_canonical <- ex_by_tx[[tx_name]] # this gets the exons for our specified canonical transcript, which we made above
introns_canonical <- in_by_tx[[tx_name]] # this gets the introns for our specified canonical transcript, which we made above

stopifnot(!is.null(exons_canonical)) # basically if it couldn't find exons, this will stop the code and let us know if something is wrong

exons_canonical # print/show me the exon structure
introns_canonical # print/show me the intron structure


# begin making the diagram, first by making a function to pull in our information about our exons and introns.

make_gene_diagram <- function(exons_gr, introns_gr, intron_cap = 400, exon_cap = Inf) {
  
  ex <- as.data.frame(exons_gr) %>%
    arrange(start, end) %>%
    transmute(type = "exon", start = start, end = end, exon_rank = exon_rank)
  
  intr <- as.data.frame(introns_gr) %>%
    arrange(start, end) %>%
    transmute(type = "intron", start = start, end = end)
  
  # Interleave exon1, intron1, exon2, intron2, ..., exonN
  pieces <- vector("list", length(ex) + length(intr))
  idx <- 1
  for (i in seq_len(nrow(ex))) {
    pieces[[idx]] <- ex[i, ]; idx <- idx + 1
    if (i <= nrow(intr)) { pieces[[idx]] <- intr[i, ]; idx <- idx + 1 }
  }
  pieces <- bind_rows(pieces)
  
  pieces <- pieces %>%
    mutate(
      genomic_len = end - start + 1,
      display_len = case_when(
        type == "intron" ~ pmin(genomic_len, intron_cap),
        type == "exon"   ~ pmin(genomic_len, exon_cap),
        TRUE ~ genomic_len
      ),
      display_start = cumsum(lag(display_len, default = 0)) + 1,
      display_end = display_start + display_len - 1
    )
  
  pieces
}

gene_pieces <- make_gene_diagram(exons_canonical, introns_canonical,intron_cap = 400, exon_cap = 800)
# pulls in our exons and introns from earlier, caps intron sizes so that way big old introns dont take up the entire diagram.
gene_pieces

# ---- 0) Confirm transcript naming matches TxDb list names ----
ex_by_tx <- exonsBy(txdb_refgene_hg19, by = "tx", use.names = TRUE)
cds_by_tx <- cdsBy(txdb_refgene_hg19, by = "tx", use.names = TRUE)

# If you set tx_name <- "NM_001791.4" earlier, but the TxDb stores "NM_001791",
# pick the matching key automatically:
if (!tx_name %in% names(ex_by_tx)) {
  tx_base <- sub("\\.\\d+$", "", tx_name)
  if (tx_base %in% names(ex_by_tx)) {
    message("tx_name not found in TxDb names; switching to: ", tx_base)
    tx_name <- tx_base
  } else {
    stop("tx_name not found in TxDb exonsBy() names: ", tx_name)
  }
}


# Refresh canonical exons/introns from the adjusted tx_name
exons_canonical <- ex_by_tx[[tx_name]]
introns_canonical <- intronsByTranscript(txdb_refgene_hg19, use.names = TRUE)[[tx_name]]

# ---- 1) Build exon spliced-coordinate table ----
exons_ord <- as.data.frame(exons_canonical) %>%
  arrange(exon_rank) %>%
  mutate(
    exon_len = end - start + 1,
    exon_tx_start = cumsum(lag(exon_len, default = 0)) + 1,
    exon_tx_end   = exon_tx_start + exon_len - 1
  )

tx_strand <- as.character(unique(strand(exons_canonical)))
if (length(tx_strand) != 1) tx_strand <- "+"

tx_to_genomic <- function(tx_pos, exons_ord, strand = "+") {
  hit <- exons_ord %>% filter(exon_tx_start <= tx_pos, exon_tx_end >= tx_pos)
  if (nrow(hit) == 0) return(NA_integer_)
  if (strand == "+") as.integer(hit$start[1] + (tx_pos - hit$exon_tx_start[1])) else
    as.integer(hit$end[1] - (tx_pos - hit$exon_tx_start[1]))
}

genomic_to_tx <- function(gpos, exons_ord, strand = "+") {
  hit <- exons_ord %>% filter(start <= gpos, end >= gpos)
  if (nrow(hit) == 0) return(NA_integer_)
  if (strand == "+") as.integer(hit$exon_tx_start[1] + (gpos - hit$start[1])) else
    as.integer(hit$exon_tx_start[1] + (hit$end[1] - gpos))
}

# ---- 2) Compute CDS start/end in transcript coordinates ----
if (!tx_name %in% names(cds_by_tx)) stop("tx_name not found in cdsBy() names: ", tx_name)
cds_gr <- cds_by_tx[[tx_name]]
stopifnot(!is.null(cds_gr), length(cds_gr) > 0)

cds_df <- as.data.frame(cds_gr) %>% arrange(start, end)
cds_tx_start <- genomic_to_tx(min(cds_df$start), exons_ord, strand = tx_strand)
cds_tx_end   <- genomic_to_tx(max(cds_df$end),   exons_ord, strand = tx_strand)
stopifnot(!is.na(cds_tx_start), !is.na(cds_tx_end))

# ---- 3) HGVS c. parser (base R; robust) ----
parse_hgvs_cdna <- function(hgvs) {
  if (is.na(hgvs)) return(NULL)
  hgvs <- trimws(hgvs)
  if (substr(hgvs, 1, 2) != "c.") return(NULL)
  
  s <- substr(hgvs, 3, nchar(hgvs))
  token <- regmatches(s, regexpr("^[*0-9\\-\\+]+", s))
  if (length(token) == 0 || is.na(token) || token == "" || token == "*" || token == "-") return(NULL)
  
  offset <- 0L
  off_match <- regmatches(token, regexpr("[+-]-?[0-9]+$", token))
  if (!is.na(off_match) && length(off_match) > 0 && off_match != "") {
    offset <- suppressWarnings(as.integer(off_match))
    token <- sub("[+-]-?[0-9]+$", "", token)
  }
  
  if (substr(token, 1, 1) == "*") {
    base <- suppressWarnings(as.integer(sub("^\\*", "", token)))
    if (is.na(base)) return(NULL)
    return(list(kind = "c_star", base = base, offset = offset))
  }
  if (substr(token, 1, 1) == "-") {
    base <- suppressWarnings(as.integer(token))
    if (is.na(base)) return(NULL)
    return(list(kind = "c_neg", base = base, offset = offset))
  }
  
  base <- suppressWarnings(as.integer(token))
  if (is.na(base)) return(NULL)
  list(kind = "c_pos", base = base, offset = offset)
}

cdna_to_genomic_safe <- function(hgvs_c, exons_ord, cds_tx_start, cds_tx_end, strand = "+") {
  tryCatch({
    parsed <- parse_hgvs_cdna(hgvs_c)
    if (is.null(parsed)) return(NA_integer_)
    
    tx_pos <- dplyr::case_when(
      parsed$kind == "c_pos"  ~ cds_tx_start + (parsed$base - 1L),
      parsed$kind == "c_neg"  ~ cds_tx_start + parsed$base,
      parsed$kind == "c_star" ~ cds_tx_end + parsed$base,
      TRUE ~ NA_integer_
    )
    
    gpos_exonic <- tx_to_genomic(tx_pos, exons_ord, strand = strand)
    if (is.na(gpos_exonic)) return(NA_integer_)
    
    if (strand == "+") as.integer(gpos_exonic + parsed$offset) else as.integer(gpos_exonic - parsed$offset)
  }, error = function(e) NA_integer_)
}

genomic_to_display <- function(gpos, gene_pieces) {
  if (is.na(gpos)) return(NA_real_)
  row <- gene_pieces %>% filter(start <= gpos, end >= gpos)
  if (nrow(row) == 0) return(NA_real_)
  frac <- (gpos - row$start[1]) / (row$genomic_len[1])
  row$display_start[1] + frac * row$display_len[1]
}

# ---- 4) Unit test (should NOT be NA for a coding HGVS) ----
message("Unit test gpos for c.476C>T: ",
        cdna_to_genomic_safe("c.476C>T", exons_ord, cds_tx_start, cds_tx_end, strand = tx_strand))

# ---- 5) Map variant_summary ----
variant_plot_df <- variant_summary %>%
  mutate(
    hgvs_for_map = dplyr::coalesce(hgvs_c, if_else(substr(mutation, 1, 2) == "c.", mutation, NA_character_)),
    gpos = vapply(
      hgvs_for_map,
      function(x) cdna_to_genomic_safe(x, exons_ord, cds_tx_start, cds_tx_end, strand = tx_strand),
      integer(1)
    ),
    display_x = vapply(gpos, genomic_to_display, numeric(1), gene_pieces = gene_pieces),
    mapped = !is.na(display_x)
  )

print(variant_plot_df %>% count(mapped))

# fix variants that may not be mapping as they were missing cdnas
aa_from_p <- function(p_hgvs) {
  if (is.na(p_hgvs)) return(NA_integer_)
  p_hgvs <- trimws(p_hgvs)
  if (substr(p_hgvs, 1, 2) != "p.") return(NA_integer_)
  suppressWarnings(as.integer(stringr::str_extract(p_hgvs, "\\d+")))
}

variant_plot_df <- variant_summary %>%
  mutate(
    aa_pos = vapply(mutation, aa_from_p, integer(1)),
    
    # Prefer canonical cDNA from QA; fallback to mutation if it's already c.; else if p. and we have AA, approximate cDNA base
    hgvs_for_map = dplyr::coalesce(
      hgvs_c,
      if_else(substr(mutation, 1, 2) == "c.", mutation, NA_character_),
      if_else(!is.na(aa_pos), paste0("c.", (aa_pos - 1L) * 3L + 1L, "A>G"), NA_character_)
      # NOTE: the "A>G" is a dummy change; we only need the POSITION for mapping
    ),
    
    gpos = vapply(
      hgvs_for_map,
      function(x) cdna_to_genomic_safe(x, exons_ord, cds_tx_start, cds_tx_end, strand = tx_strand),
      integer(1)
    ),
    display_x = vapply(gpos, genomic_to_display, numeric(1), gene_pieces = gene_pieces),
    mapped = !is.na(display_x)
  )
left_bucket_x <- min(gene_pieces$display_start) - 200

variant_plot_df <- variant_plot_df %>%
  mutate(
    is_upstream_unmapped = !mapped & !is.na(hgvs_for_map) & startsWith(hgvs_for_map, "c.-"),
    display_x2 = if_else(is_upstream_unmapped, left_bucket_x, display_x),
    mapped2 = mapped | is_upstream_unmapped
  )

# how aggressively to bin x positions (smaller = more stacking groups)
bin_width <- 10

variant_plot_df2 <- variant_plot_df %>%
  filter(mapped2) %>%
  mutate(
    x_plot = display_x2,
    x_bin = round(x_plot / bin_width) * bin_width
  ) %>%
  group_by(x_bin) %>%
  arrange(classification, ps3_applied, mutation) %>%
  mutate(
    stack_i = row_number(),
    y_plot = 0.35 + 0.08 * (stack_i - 1) # baseline + step per stacked variant
  ) %>%
  ungroup()

# Grab the list of variant types (ie, missense, intronic synonymous) so we can use them in our plot. 
# this facets the graph by variant type (ie, makes 4 separate graphs for 4 variant types)
# if we have a ton of variant types though this can look messy. so im keeping this code in case we want to go that route.
# variant_plot_df2 <- variant_plot_df2 %>%
#  mutate(
#    variant_type_facet = fct_infreq(variant_type)
#  )


# ---- 6) Plot ----
exon_boxes   <- gene_pieces %>% filter(type == "exon")
intron_lines <- gene_pieces %>% filter(type == "intron")

# Make a stable mapping from variant_type -> shape codes
vt_levels <- sort(unique(variant_plot_df2$variant_type))

shape_pool <- c(21, 22, 24, 23, 25, 8, 4, 3, 1, 2)  # extend if you ever need more. the numbers are just different shapes (circle, triangle, etc)
vt_shape_map <- setNames(shape_pool[seq_along(vt_levels)], vt_levels)


# Make sure ps3_applied is TRUE/FALSE (not "TRUE"/"FALSE" strings)
variant_plot_df2 <- variant_plot_df2 %>%
  mutate(ps3_applied = as.logical(ps3_applied))

p <- ggplot(data = variant_plot_df2) +
  geom_segment(
    data = intron_lines,
    aes(x = display_start, xend = display_end, y = 0, yend = 0),
    linewidth = 0.6
  ) +
  geom_rect(
    data = exon_boxes,
    aes(xmin = display_start, xmax = display_end, ymin = -0.25, ymax = 0.25),
    linewidth = 0.6, fill = "grey50", color = "grey20"
  ) +
  geom_segment(
    aes(x = x_plot, xend = x_plot, y = 0.18, yend = y_plot, color = classification),
    linewidth = 0.3, alpha = 0.7
  ) +
  geom_point(
    aes(
      x = x_plot,
      y = y_plot,
      color = classification,
      shape = variant_type,
      fill  = ps3_applied
    ),
    size = 2.8,
    stroke = 0.9,
    alpha = 0.95
  ) +
  geom_text(
    data = variant_plot_df2 %>% dplyr::filter(is_upstream_unmapped),
    aes(x = x_plot + 30, y = y_plot, label = hgvs_for_map),
    size = 3, hjust = 0
  ) +
  scale_shape_manual(values = vt_shape_map) +
  scale_fill_manual(
    values = c(`TRUE` = "black", `FALSE` = "white"),
    name = "PS3 applied"
  ) +
  guides(
    fill = guide_legend(
      override.aes = list(
        shape = 21,
        color = "black",
        stroke = 0.9,
        size = 3
      )
    )
  ) +
  scale_y_continuous(NULL, breaks = NULL) +
  labs(
    x = "5′ → 3′ (schematic)",
    title = paste0(unique(variant_plot_df2$gene_symbol), " variants (", tx_name, ", hg19)"),
    subtitle = "Color = classification; Shape = variant type; Fill = PS3 applied"
  ) +
  theme_classic()

print(p)

# Make plot reliant on p. nomenclature instead of cdna
variant_plot_protein <- variant_plot_df2 %>%
  mutate(
    aa_pos = case_when(
      str_detect(mutation, "^p\\.") ~ as.integer(
        str_extract(mutation, "(?<=\\D)\\d+(?=\\D)")
      ),
      TRUE ~ NA_integer_
    ),
    is_protein_mapped = !is.na(aa_pos)
  )

variant_plot_protein <- variant_plot_protein %>%
  filter(is_protein_mapped) %>%
  arrange(aa_pos, mutation) %>%
  group_by(aa_pos) %>%
  mutate(stack_i = row_number()) %>%
  ungroup() %>%
  mutate(
    x_plot_aa = aa_pos,
    y_plot_aa = 0.35 + 0.08 * (stack_i - 1)
  )

p_protein <- ggplot(data = variant_plot_protein) +
  # protein baseline
  geom_segment(
    aes(
      x = min(x_plot_aa),
      xend = max(x_plot_aa),
      y = 0,
      yend = 0
    ),
    linewidth = 0.6
  ) +
  # lollipop stems
  geom_segment(
    aes(
      x = x_plot_aa,
      xend = x_plot_aa,
      y = 0,
      yend = y_plot_aa,
      color = classification
    ),
    linewidth = 0.3,
    alpha = 0.7
  ) +
  # lollipop heads
  geom_point(
    aes(
      x = x_plot_aa,
      y = y_plot_aa,
      color = classification,
      shape = variant_type,
      fill = ps3_applied
    ),
    size = 2.8,
    stroke = 0.9,
    alpha = 0.95
  ) +
  scale_shape_manual(values = vt_shape_map) +
  scale_fill_manual(
    values = c(`TRUE` = "black", `FALSE` = "white"),
    name = "PS3 applied"
  ) +
  guides(
    fill = guide_legend(
      override.aes = list(
        shape = 21,
        color = "black",
        stroke = 0.9,
        size = 3
      )
    )
  ) +
  scale_y_continuous(NULL, breaks = NULL) +
  labs(
    x = "Protein position (amino acid)",
    title = paste0(unique(variant_plot_protein$gene_symbol),
                   " variants (", tx_name, ")"),
    subtitle = "Protein-space lollipop: Color = classification; Shape = variant type; Fill = PS3 applied"
  ) +
  theme_classic()

print(p_protein)

# Load UniProt JSON (from your uploaded file)
uniprot <- jsonlite::fromJSON("alg9_domains.json")

domains_df <- tibble(
  feature_type = uniprot$features$type,
  aa_start = as.integer(uniprot$features$location$start$value),
  aa_end   = as.integer(uniprot$features$location$end$value),
  ligand   = uniprot$features$ligand$name
) %>%
  filter(!is.na(aa_start), !is.na(aa_end)) %>%
  mutate(
    aa_min = pmin(aa_start, aa_end),
    aa_max = pmax(aa_start, aa_end),
    label = ifelse(
      feature_type,
      paste0(feature_type, " (", ligand, ")")
    )
  )

# attempting to build a plot that goes by cdna position and aa position.

cds_by_tx <- cdsBy(txdb_refgene_hg19, by = "tx", use.names = TRUE)
cds_canonical <- cds_by_tx[[tx_name]]
stopifnot(!is.null(cds_canonical))

cds_df <- as.data.frame(cds_canonical) %>%
  arrange(start, end) %>%
  mutate(seg_len = end - start + 1) %>%
  mutate(
    cds_start_bp = cumsum(lag(seg_len, default = 0)) + 1,
    cds_end_bp   = cds_start_bp + seg_len - 1,
    aa_start     = floor((cds_start_bp - 1) / 3) + 1,
    aa_end       = floor((cds_end_bp - 1) / 3) + 1
  )

protein_lane_y0 <- 0 # offsets the protein variants from the noncoding variants
noncoding_lane_y0 <- -0.2

variant_plot_protein2 <- variant_plot_df2 %>%
  mutate(
    aa_pos = case_when(
      str_detect(mutation, "^p\\.") ~ as.integer(str_extract(mutation, "(?<=\\D)\\d+(?=\\D)")),
      TRUE ~ NA_integer_
    ),
    is_protein_mapped = !is.na(aa_pos),
    
    # For protein mapped variants, x is AA position
    x_aa = aa_pos,
    
    # label for noncoding lane (use hgvs_for_map if you have it; else mutation)
    nc_label = if_else(!is_protein_mapped, coalesce(hgvs_for_map, hgvs_c, mutation), NA_character_)
  )


# Protein lane stacking (stack variants that share same AA)
protein_df <- variant_plot_protein2 %>%
  filter(is_protein_mapped) %>%
  arrange(x_aa, mutation) %>%
  group_by(x_aa) %>%
  mutate(stack_i = row_number()) %>%
  ungroup() %>%
  mutate(
    y_base = protein_lane_y0,
    y_plot_aa = y_base + 0.35 + 0.08 * (stack_i - 1)  # same style as your gene plot
  )

# Noncoding lane stacking (no AA coordinate -> bucket them along the left with small offsets)
# Here we set x_aa to (min_aa - 5) so they appear left of AA1, and stack downward.
min_aa <- suppressWarnings(min(protein_df$x_aa, na.rm = TRUE))
if (!is.finite(min_aa)) min_aa <- 1

noncoding_df <- variant_plot_protein2 %>%
  filter(!is_protein_mapped) %>%
  arrange(variant_type, mutation) %>%
  mutate(
    x_aa = min_aa - 5L  # a “noncoding column” left of the protein
  ) %>%
  group_by(x_aa) %>%
  mutate(stack_i = row_number()) %>%
  ungroup() %>%
  mutate(
    y_base = noncoding_lane_y0,
    y_plot_aa = y_base - 0.05 - 0.08 * (stack_i - 1)  # stack downward from the NC baseline
  )

variant_plot_protein_plot <- bind_rows(protein_df, noncoding_df)

# compute baseline ranges once (avoids ggplot warnings)
protein_xmin <- min(protein_df$x_aa, na.rm = TRUE)
protein_xmax <- max(protein_df$x_aa, na.rm = TRUE)

noncoding_xmin <- min(noncoding_df$x_aa, na.rm = TRUE)
noncoding_xmax <- max(noncoding_df$x_aa, na.rm = TRUE)

# Label CDS exons using transcript exon rank (ie numbers exons correctly due to there being non-coding exons in some genes)
cds_df_labeled <- cds_df %>%
  dplyr::arrange(aa_start, aa_end) %>%   # stable ordering for display
  dplyr::mutate(
    exon_mid = (aa_start + aa_end) / 2,
    exon_label = exon_rank               # or paste0("Ex", exon_rank) if you prefer
  )

# Attempting to set up so that the fill of the shape matches the colors of the lollipops
classification_levels <- sort(unique(variant_plot_df2$classification))

classification_color_map <- setNames(
  scales::hue_pal()(length(classification_levels)),
  classification_levels
)

variant_plot_protein_plot <- variant_plot_protein_plot %>%
  dplyr::mutate(
    fill_color = ifelse(ps3_applied, as.character(classification), NA)
  )

# ---- Prep for plot aesthetics ----

# Make sure these are logical
variant_plot_protein_plot <- variant_plot_protein_plot %>%
  dplyr::mutate(
    any_weighted_proband = as.logical(any_weighted_proband),
    ps3_applied = as.logical(ps3_applied)
  )

# Legend-only data. Updating to show ps3 = filled and dotted = weighted
ps3_legend_df <- data.frame(
  ps3_key = factor(c("PS3 applied", "PS3 not applied"),
                   levels = c("PS3 applied", "PS3 not applied")),
  x = -Inf, y = -Inf
)

weighted_legend_df <- data.frame(
  weight_key = factor(c("Weighted proband", "Unweighted proband"),
                      levels = c("Weighted proband", "Unweighted proband")),
  x = -Inf, xend = -Inf,
  y = -Inf, yend = -Inf
)

# Directional intron chevrons between CDS exons (not to scale). this is so it looks more accurate and exons aren't squished together.
# Assumes you have cds_df_labeled with aa_start/aa_end and exon_rank
intron_chevrons_df <- cds_df_labeled %>%
  dplyr::arrange(aa_start, aa_end) %>%
  dplyr::mutate(
    chevron_x = aa_end + 0.6   # place slightly to the right of each exon end
  ) %>%
  dplyr::slice(-dplyr::n()) %>%  # no chevron after last exon
  dplyr::mutate(
    intron_key = factor("Intron (not to scale)", levels = "Intron (not to scale)"),
    chevron_label = "\u203A\u203A"  # '››' (directional chevrons)
  )

# Legend-only dummy row (kept off-plot)
intron_legend_df <- data.frame(
  intron_key = factor("Intron (not to scale)", levels = "Intron (not to scale)"),
  x = -Inf, y = -Inf,
  chevron_label = "\u203A\u203A"
)

intron_legend_df <- data.frame(
  intron_key = factor("Intron (not to scale)", levels = "Intron (not to scale)"),
  x = -Inf, y = -Inf
)


# finally make the updated plot

p_protein <- ggplot() +
  
  # --- legend-only: PS3 applied (filled vs hollow) ---
  # We use alpha instead of fill to avoid conflicting with your classification fill scale.
  geom_point(
    data = ps3_legend_df,
    aes(x = x, y = y, alpha = ps3_key),
    inherit.aes = FALSE,
    shape = 21, size = 3,
    color = "black", fill = "black"
  ) +
  
  # --- legend-only: weighted proband (dotted vs solid) ---
  geom_segment(
    data = weighted_legend_df,
    aes(x = x, xend = xend, y = y, yend = yend, linetype = weight_key),
    inherit.aes = FALSE,
    linewidth = 0.8,
    color = "black"
  ) +
  
  # --- legend-only: intron (not to scale) ---
  # Uses linetype (NOT color) so it won't interfere with classification colors.
  geom_segment(
    data = data.frame(
      line_key = factor("Intron (\u203A\u203A not to scale)", levels = "Intron (\u203A\u203A not to scale)"),
      x = -Inf, xend = -Inf, y = -Inf, yend = -Inf
    ),
    aes(x = x, xend = xend, y = y, yend = yend, linetype = line_key),
    inherit.aes = FALSE,
    linewidth = 0.8,
    color = "black"
  ) +
  
  # --- protein baseline ---
  annotate(
    "segment",
    x = protein_xmin, xend = protein_xmax,
    y = protein_lane_y0, yend = protein_lane_y0,
    linewidth = 0.6
  ) +
  # --- noncoding baseline ---
  annotate(
    "segment",
    x = noncoding_xmin, xend = noncoding_xmax,
    y = noncoding_lane_y0, yend = noncoding_lane_y0,
    linewidth = 0.6
  ) +
  # label the noncoding lane
  annotate(
    "text",
    x = noncoding_xmin, y = noncoding_lane_y0 - 0.10,
    label = "Non-coding",
    size = 3, vjust = 1, hjust = 0
  ) +
  
  # --- UniProt feature track ---
  geom_rect(
    data = domains_df,
    aes(xmin = aa_min, xmax = aa_max, ymin = -0.22, ymax = -0.14),
    inherit.aes = FALSE,
    alpha = 0.9
  ) +
  geom_text(
    data = domains_df,
    aes(x = (aa_min + aa_max) / 2, y = -0.24, label = label),
    inherit.aes = FALSE,
    size = 2.8, vjust = 1
  ) +
  
  # --- CDS exon track (coding-only) ---
  geom_rect(
    data = cds_df_labeled,
    aes(xmin = aa_start, xmax = aa_end, ymin = -0.12, ymax = 0.12),
    inherit.aes = FALSE,
    fill = "grey60", color = "grey30", linewidth = 0.4
  ) +
  
  # --- intron chevrons between CDS exons (directional; not to scale) ---
  # Make them more apparent (bigger, bold, black, slightly above baseline)
  geom_text(
    data = intron_chevrons_df,
    aes(x = chevron_x, y = 0.02, label = chevron_label),
    inherit.aes = FALSE,
    size = 6,
    fontface = "bold",
    color = "black"
  ) +
  
  # exon ranks (true transcript exon numbering)
  geom_text(
    data = cds_df_labeled,
    aes(x = exon_mid, y = 0.14, label = exon_rank),
    inherit.aes = FALSE,
    size = 2.8
  ) +
  
  # --- stems (dotted if weighted proband) ---
  geom_segment(
    data = variant_plot_protein_plot,
    aes(
      x = x_aa, xend = x_aa,
      y = y_base, yend = y_plot_aa,
      color = classification,
      linetype = any_weighted_proband
    ),
    linewidth = 0.3,
    alpha = 0.7
  ) +
  
  # --- heads (fill matches classification IF ps3_applied) ---
  geom_point(
    data = variant_plot_protein_plot,
    aes(
      x = x_aa,
      y = y_plot_aa,
      color = classification,
      shape = variant_type,
      fill  = fill_color
    ),
    size = 2.8,
    stroke = 0.9,
    alpha = 0.95
  ) +
  
  # --- optional noncoding labels ---
  geom_text(
    data = noncoding_df,
    aes(x = x_aa + 1, y = y_plot_aa, label = nc_label),
    inherit.aes = FALSE,
    size = 2.8, hjust = 0
  ) +
  
  # --- scales ---
  scale_shape_manual(values = vt_shape_map, name = "variant_type") +
  scale_color_manual(values = classification_color_map, name = "classification") +
  scale_fill_manual(values = classification_color_map, na.value = "white", guide = "none") +
  
  # PS3 legend (alpha-only so it doesn't touch your fill scale)
  scale_alpha_manual(
    name = "PS3",
    values = c("PS3 applied" = 1, "PS3 not applied" = 0),
    breaks = c("PS3 applied", "PS3 not applied"),
    guide = guide_legend(override.aes = list(shape = 21, fill = "black", color = "black"))
  ) +
  
  # linetype scale: supports dummy legend labels + TRUE/FALSE from your data
  scale_linetype_manual(
    name = "Line key",
    values = c(
      "Weighted proband" = "dotted",
      "Unweighted proband" = "solid",
      "Intron (\u203A\u203A not to scale)" = "twodash",
      "TRUE" = "dotted",
      "FALSE" = "solid"
    ),
    breaks = c("Weighted proband", "Unweighted proband", "Intron (\u203A\u203A not to scale)")
  ) +
  
  scale_y_continuous(NULL, breaks = NULL) +
  labs(
    x = "Protein position (AA)",
    title = paste0(unique(variant_plot_df2$gene_symbol), " variants (", tx_name, ")"),
    subtitle = "Protein lane: AA-mapped variants. Lower lane: non-coding variants."
  ) +
  theme_classic() +
  theme(
    legend.position = "right",
    legend.box = "vertical"
  )

print(p_protein)