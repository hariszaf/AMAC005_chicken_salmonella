library(dplyr)
library(stringr)  # str_replace_all()
library(tibble)  # row_names_to_columns

# Load metadata (treatments)
metadata = readxl::read_excel("microbetag_analysis/sample_metadata.xlsx", sheet = "sample_metadata")
transposed_df <- as.data.frame(t(metadata))
colnames(transposed_df) <- transposed_df[1, ]  # Set first row as column names
metadata_df <- transposed_df[-1, ]
metadata_df <- tibble::rownames_to_column(metadata_df, var="sample")

# Load metabolites data
metabolites = readxl::read_excel("microbetag_analysis/sample_metadata.xlsx", sheet = "mets_of_interest")
transposed_df <- as.data.frame(t(metabolites))
colnames(transposed_df) <- transposed_df[1, ]  # Set first row as column names
metabolites_df <- transposed_df[-1, ]
metabolites_df <-tibble::rownames_to_column(metabolites_df, var="individual")
metabolites_df <- metabolites_df %>%
  mutate(individual = str_replace_all(individual, "M",""))

# Merge metadata with metabolites in consideration
all_metadata <- right_join(metadata_df, metabolites_df, by="individual")
all_metadata_t <- as.data.frame(t(all_metadata))

# Save to file
write.table(all_metadata_t, "microbetag_analysis/microbetag_input/all_metadata_mets_transp.tsv",
            sep = "\t", quote=FALSE, row.names = TRUE, col.names = FALSE)

# ---------------------------------------------------------------------------------

# Split per day
day_7_metadata <- all_metadata %>% filter(day==7)    # n = 15
day_14_metadata <- all_metadata %>% filter(day==14)  # n = 21
day_21_metadata <- all_metadata %>% filter(day==21)  # n = 25
day_28_metadata <- all_metadata %>% filter(day==28)  # n = 25
day_35_metadata <- all_metadata %>% filter(day==35)  # n = 25
# Total number of samples = 111

# Load genome counts
genome_count_table <- read.table("microbetag_analysis/microbetag_input/salm_multio_abd_table_gc.tsv", sep = "\t", header = TRUE)

day_7_abd <- genome_count_table %>%
  select(any_of(c("genome", day_7_metadata$sample, "classification")))
write.table(day_7_abd, "microbetag_analysis/microbetag_input/per_day/abd_day_7.tsv", sep = "\t", col.names = TRUE, row.names = FALSE)

day_14_abd <- genome_count_table %>%
  select(any_of(c("genome", day_14_metadata$sample, "classification")))
write.table(day_14_abd, "microbetag_analysis/microbetag_input/per_day/abd_day_14.tsv", sep = "\t", col.names = TRUE, row.names = FALSE)

day_21_abd <- genome_count_table %>%
  select(any_of(c("genome", day_21_metadata$sample, "classification")))
write.table(day_21_abd, "microbetag_analysis/microbetag_input/per_day/abd_day_21.tsv", sep = "\t", col.names = TRUE, row.names = FALSE)

day_28_abd <- genome_count_table %>%
  select(any_of(c("genome", day_28_metadata$sample, "classification")))
write.table(day_28_abd, "microbetag_analysis/microbetag_input/per_day/abd_day_28.tsv", sep = "\t", col.names = TRUE, row.names = FALSE)

day_35_abd <- genome_count_table %>%
  select(any_of(c("genome", day_35_metadata$sample, "classification")))
write.table(day_35_abd, "microbetag_analysis/microbetag_input/per_day/abd_day_35.tsv", sep = "\t", col.names = TRUE, row.names = FALSE)

# Remember FlashWeave only considers .tsv extension for metadata files
write.table(t(select(day_7_metadata,  -c(sample, individual))), "microbetag_analysis/microbetag_input/per_day/all_metadata_day_7.tsv", sep = "\t", quote=FALSE, row.names = TRUE, col.names = FALSE)
write.table(t(select(day_14_metadata, -c(sample, individual))), "microbetag_analysis/microbetag_input/per_day/all_metadata_day_14.tsv", sep = "\t", quote=FALSE, row.names = TRUE, col.names = FALSE)
write.table(t(select(day_21_metadata, -c(sample, individual))), "microbetag_analysis/microbetag_input/per_day/all_metadata_day_21.tsv", sep = "\t", quote=FALSE, row.names = TRUE, col.names = FALSE)
write.table(t(select(day_28_metadata, -c(sample, individual))), "microbetag_analysis/microbetag_input/per_day/all_metadata_day_28.tsv", sep = "\t", quote=FALSE, row.names = TRUE, col.names = FALSE)
write.table(t(select(day_35_metadata, -c(sample, individual))), "microbetag_analysis/microbetag_input/per_day/all_metadata_day_35.tsv", sep = "\t", quote=FALSE, row.names = TRUE, col.names = FALSE)


# ---------------------------------------------------------------------------------

# Split per treatment
tg_1_metadata <- all_metadata %>% filter(treatment=="TG1")  # n = 20
tg_2_metadata <- all_metadata %>% filter(treatment=="TG2")  # n = 23
tg_3_metadata <- all_metadata %>% filter(treatment=="TG3")  # n = 24
tg_4_metadata <- all_metadata %>% filter(treatment=="TG4")  # n = 22
tg_5_metadata <- all_metadata %>% filter(treatment=="TG5")  # n = 22
# Total number of samples = 111

treatment_groups <- list(
  "TG1" = tg_1_metadata,
  "TG2" = tg_2_metadata,
  "TG3" = tg_3_metadata,
  "TG4" = tg_4_metadata,
  "TG5" = tg_5_metadata
)

for (tg_name in names(treatment_groups)) {
  tg_abd <- genome_count_table %>%
    select(any_of(c("genome", treatment_groups[[tg_name]]$sample, "classification")))

  output_path <- paste0("microbetag_analysis/microbetag_input/per_treatment/abd_", tg_name, ".tsv")
  write.table(tg_abd, output_path, sep = "\t", col.names = TRUE, row.names = FALSE)

  output_path <- paste0("microbetag_analysis/microbetag_input/per_treatment/all_metadata_", tg_name, ".tsv")
  write.table(t(select(treatment_groups[[tg_name]], -c(sample, individual))), output_path, sep = "\t", quote=FALSE, row.names = TRUE, col.names = FALSE)
}


# Treatment of missing values
# FlashWeave currently does not support missing data, please remove all samples with missing entries (both in OTU and meta data tables) prior to running FlashWeave.

# I had to remove samples with null values since flashweave cannot handle them
# TG1: D300477
# TG0: D300530 and D300531
# also, we had to make sure samples appear in the same order in the abundance table and the metadata file (as they did already)


# PoultryStar taxa
# Bifidobacterium animalis, Enterococcus faecium, Lactobacillus salivarius, Lactobacillus reuteri and Pediococcus acidilactici
# GEXTRA:bin_000001,  GEXTRA:bin_000002, GPB:bin_000025, GEXTRA:bin_000004,  GEXTRA:bin_000006
