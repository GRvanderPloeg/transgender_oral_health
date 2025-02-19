library(tidyverse)
library(biomformat)

# Primary dataset
df = read_biom("./Data/GOH TRANS Roel/GOH TRANS Roel/Miseq/UNOISE_all_incl_pos_ctrl/GoHT_all_zotu_table_sorted.biom")

# Count data extraction
counts = biom_data(df) %>% as.matrix() %>% t() %>% as_tibble()

# Taxonomy extraction
taxa = cbind(rownames(df), observation_metadata(df)) %>% as_tibble()
colnames(taxa) = c("zOTU", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# Sample metadata extraction
mapping = read.csv("./Data/GOH TRANS Roel/GOH TRANS Roel/Miseq/mapping_file_GOHTRANS_20240116.txt", sep="\t") %>% as_tibble()
sampleNames = biom_data(df) %>% as.matrix() %>% colnames()
sampleInfo = sampleNames %>% as.tibble() %>% left_join(mapping, by=c("value"="X.SampleID"))

# Extra metadata that I'm currently not using
temp = read.csv("./Data/GOH-TRANS_csv_export_20240205114955/GOH-TRANS_export_20240205.csv", sep=";") %>% as_tibble()

# Remove controls
counts = counts[33:nrow(counts),]
sampleInfo = sampleInfo[33:nrow(sampleInfo),]

# Repair subject names
sampleInfo$subject = str_split_fixed(sampleInfo$value, "\\.", 3)[,2]

# Repair timepoints
sampleInfo = sampleInfo %>% mutate(newTimepoint = 0)
sampleInfo[sampleInfo$Timepoint == "Baseline_V1", "newTimepoint"] = 0
sampleInfo[sampleInfo$Timepoint == "3_months_V2", "newTimepoint"] = 3
sampleInfo[sampleInfo$Timepoint == "6_months_V3", "newTimepoint"] = 6
sampleInfo[sampleInfo$Timepoint == "9_months_V3", "newTimepoint"] = 9
sampleInfo[sampleInfo$Timepoint == "1_year_V4", "newTimepoint"] = 12

# Add hormone levels to tongue data
ylabs = c("Total_protein_ug_ml", "S_IgA_ug_ml", "MUC_5_B_ug_ml", "Amylase_U_ml", "Lysozyme_U_ml", "Chitinase_AU_ml", "Estradiol_pmol_ml", "LH_U_L", "SHBG_nmol_L", "Testosterone_nmol_L", "Free_testosterone_Vermeulen_pmol_L")
hormoneInfo = sampleInfo %>% filter(Niche == "Saliva") %>% select(-Niche) %>% select(c(subject, newTimepoint, ylabs))
sampleInfo = sampleInfo %>% select(-all_of(ylabs)) %>% left_join(hormoneInfo)

# Export full datasets
write.table(counts, "./Data/counts_fixed.csv", col.names=FALSE, row.names=FALSE)
write.table(sampleInfo, "./Data/sampleInfo_fixed.csv")
write.table(taxa, "./Data/taxonomy_fixed.csv")

# Split datasets
tongue = counts[sampleInfo$Niche == "Tongue",]
tongueSampleMeta = sampleInfo[sampleInfo$Niche == "Tongue",]
saliva = counts[sampleInfo$Niche == "Saliva",]
salivaSampleMeta = sampleInfo[sampleInfo$Niche == "Saliva",]

# Feature selection not performed to allow for CLR at a later point

# Sample selection based on total counts
tongue_threshold = 5000
saliva_threshold = 5000

tongueSampleSelection = rowSums(tongue) >= tongue_threshold
salivaSampleSelection = rowSums(saliva) >= saliva_threshold

# Do sample selection immediately as it will not affect CLR
tongue = tongue[tongueSampleSelection,]
tongueSampleMeta = tongueSampleMeta[tongueSampleSelection,]

saliva = saliva[salivaSampleSelection,]
salivaSampleMeta = salivaSampleMeta[salivaSampleSelection,]

# Save split data
write.table(tongue, "./Data/tongueCounts.csv", sep=",", row.names=FALSE, col.names=FALSE)
write.table(tongueSampleMeta, "./Data/tongueSampleMeta.csv", sep=",", row.names=FALSE, col.names=FALSE)
write.table(saliva, "./Data/salivaCounts.csv", sep=",", row.names=FALSE, col.names=FALSE)
write.table(salivaSampleMeta, "./Data/salivaSampleMeta.csv", sep=",", row.names=FALSE, col.names=FALSE)

