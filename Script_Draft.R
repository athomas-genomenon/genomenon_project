#working directory
setwd("~/Desktop/Genomenon")

#csv imports
evidence <- read.csv("evidence-7293-2026-01-08.csv", header = TRUE, sep = ",")
qa <- read.csv("qa-7293-2026-01-08.csv", header = TRUE, sep = ",")

#install packages
install.packages("ggplot2")
install.packages("dplyr")
library(ggplot2)
library(dplyr)
library(stringr)

#removing columns we don't want
cleaned_evidence = subset(evidence, select = c("gene_symbol", "mutation", "classification", "all_cats"))

#concatenated files
concat_evidence <- evidence %>%
  left_join(
    qa %>% select(gene, variant_desc, variant_type),
    by = c("gene_symbol" = "gene", "mutation" = "variant_desc")
  )

#cleaned to include just the columns we want
cleaned_csvs <- concat_evidence %>%
  select(
    gene_symbol,
    mutation,
    variant_type,
    classification,
    all_cats,
    evidence_type,
    unweighted,
    number_of_probands
  )

#ps3 applied or not
evidence_ps3 <- cleaned_csvs %>%
  mutate(
    ps3_applied_row =
      !is.na(all_cats) & str_detect(all_cats, "(^|\\|)ps3(\\||$|_)")
  )

#keeps weighted probands
evidence_weighted <- evidence_ps3 %>%
  mutate(
    is_weighted_proband =
      evidence_type == "unrelated_proband" & unweighted == FALSE
  )

# Creates dataframe without duplicates. It should look like:
# One row per variant, with:
# variant identity (gene_symbol, mutation)
# variant_type
# classification
# weighted_probands (numeric)
# any_weighted_proband (TRUE/FALSE)

variant_summary <- evidence_weighted %>% 
  group_by(gene_symbol, mutation) %>% # treats all evidence rows for the same variant (as there is multiple) as belonging together.
  summarise( # make a summary for each variant
    variant_type = first(variant_type), # carry on the variant type from our og file
    classification = first(classification), # carry on the classification from our og file
    
    weighted_probands = sum(number_of_probands[is_weighted_proband], na.rm = TRUE), # adds up all the weighted probands
    any_weighted_proband = any(is_weighted_proband, na.rm = TRUE), # is there any weighted proband evidence for this variant?
    
    ps3_applied = any(ps3_applied_row, na.rm = TRUE), 
    
    n_evidence_rows = n(),
    .groups = "drop"
  )

#notes for building the lollipop plot itself:
#what is y axis? number of variants?
#shapes for variant classification
#colors for variant type
#colored in shape for if it meets ps3/not colored in
#5 to 3 on x-axis (with exons/introns labeled)

lollipop <- ggplot(variant_summary, aes(x = mutation, y = weighted_probands))
            
#labs(lollipop, title =  "weighted proband variant plot", x = "variants", y = "number of weighted probands")




#how to input json data?