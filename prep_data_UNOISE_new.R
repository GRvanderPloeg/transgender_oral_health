library(tidyverse)
library(biomformat)

# UNOISE tongue
df2 = read_biom("./Data/GoHTrans_rerun_Bernd/UNOISE_TP_ctrls/GoHT_TP_zotu_table_sorted.biom")

# UNOISE saliva
df3 = read_biom("./Data/GoHTrans_rerun_Bernd/UNOISE_US_ctrls/GoHT_US_zotu_table_sorted.biom")

# Count data extraction
countsTongue = biom_data(df2) %>% as.matrix() %>% t() %>% as_tibble()
countsSaliva = biom_data(df3) %>% as.matrix() %>% t() %>% as_tibble()

# Taxonomy extraction
taxaTongue = cbind(rownames(df2), observation_metadata(df2)) %>% as_tibble()
taxaSaliva = cbind(rownames(df3), observation_metadata(df3)) %>% as_tibble()
colnames(taxaTongue) = c("zOTU", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
colnames(taxaSaliva) = c("zOTU", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# Sample metadata extraction
mapping = read.csv("./Data/GOH TRANS Roel/GOH TRANS Roel/Miseq/mapping_file_GOHTRANS_20240116.txt", sep="\t") %>% as_tibble()
sampleNamesTongue = biom_data(df2) %>% as.matrix() %>% colnames()
sampleNamesSaliva = biom_data(df3) %>% as.matrix() %>% colnames()

sampleInfoTongue = sampleNamesTongue %>% as.tibble() %>% left_join(mapping, by=c("value"="X.SampleID"))
sampleInfoSaliva = sampleNamesSaliva %>% as.tibble() %>% left_join(mapping, by=c("value"="X.SampleID"))

# Extra metadata that I'm currently not using
temp = read.csv("./Data/GOH-TRANS_csv_export_20240205114955/GOH-TRANS_export_20240205.csv", sep=";") %>% as_tibble()

# Remove controls
countsTongue = countsTongue[17:nrow(countsTongue),]
countsSaliva = countsSaliva[9:nrow(countsSaliva),]

sampleInfoTongue = sampleInfoTongue[17:nrow(sampleInfoTongue),]
sampleInfoSaliva = sampleInfoSaliva[9:nrow(sampleInfoSaliva),]

# Repair subject names
sampleInfoTongue$subject = str_split_fixed(sampleInfoTongue$value, "\\.", 3)[,2]
sampleInfoSaliva$subject = str_split_fixed(sampleInfoSaliva$value, "\\.", 3)[,2]

# Repair timepoints
sampleInfoTongue = sampleInfoTongue %>% mutate(newTimepoint = 0)
sampleInfoTongue[sampleInfoTongue$Timepoint == "Baseline_V1", "newTimepoint"] = 0
sampleInfoTongue[sampleInfoTongue$Timepoint == "3_months_V2", "newTimepoint"] = 3
sampleInfoTongue[sampleInfoTongue$Timepoint == "6_months_V3", "newTimepoint"] = 6
sampleInfoTongue[sampleInfoTongue$Timepoint == "9_months_V3", "newTimepoint"] = 9
sampleInfoTongue[sampleInfoTongue$Timepoint == "1_year_V4", "newTimepoint"] = 12

sampleInfoSaliva = sampleInfoSaliva %>% mutate(newTimepoint = 0)
sampleInfoSaliva[sampleInfoSaliva$Timepoint == "Baseline_V1", "newTimepoint"] = 0
sampleInfoSaliva[sampleInfoSaliva$Timepoint == "3_months_V2", "newTimepoint"] = 3
sampleInfoSaliva[sampleInfoSaliva$Timepoint == "6_months_V3", "newTimepoint"] = 6
sampleInfoSaliva[sampleInfoSaliva$Timepoint == "9_months_V3", "newTimepoint"] = 9
sampleInfoSaliva[sampleInfoSaliva$Timepoint == "1_year_V4", "newTimepoint"] = 12

# Add hormone levels to tongue data
ylabs = c("Total_protein_ug_ml", "S_IgA_ug_ml", "MUC_5_B_ug_ml", "Amylase_U_ml", "Lysozyme_U_ml", "Chitinase_AU_ml", "Estradiol_pmol_ml", "LH_U_L", "SHBG_nmol_L", "Testosterone_nmol_L", "Free_testosterone_Vermeulen_pmol_L")
hormoneInfo = sampleInfoSaliva %>% select(c(subject, newTimepoint, ylabs))
sampleInfoTongue = sampleInfoTongue %>% select(-all_of(ylabs)) %>% left_join(hormoneInfo) %>% select(all_of(colnames(sampleInfoSaliva)))

# Export full datasets
write.table(countsTongue, "./Data/20240503_UNOISE_new/tongueCounts_fixed.csv", col.names=FALSE, row.names=FALSE)
write.table(countsSaliva, "./Data/20240503_UNOISE_new/salivaCounts_fixed.csv", col.names=FALSE, row.names=FALSE)
write.table(sampleInfoTongue, "./Data/20240503_UNOISE_new/sampleInfoTongue_fixed.csv")
write.table(sampleInfoSaliva, "./Data/20240503_UNOISE_new/sampleInfoSaliva_fixed.csv")
write.table(taxaTongue, "./Data/20240503_UNOISE_new/taxonomyTongue_fixed.csv")
write.table(taxaSaliva, "./Data/20240503_UNOISE_new/taxonomySaliva_fixed.csv")

# Feature selection not performed to allow for CLR at a later point

# Sample selection based on total counts
tongue_threshold = 5000
saliva_threshold = 5000

tongueSampleSelection_new = rowSums(countsTongue) >= tongue_threshold
salivaSampleSelection_new = rowSums(countsSaliva) >= saliva_threshold

# Do sample selection immediately as it will not affect CLR
countsTongue_subset = countsTongue[tongueSampleSelection_new,]
sampleInfoTongue_subset = sampleInfoTongue[tongueSampleSelection_new,]

countsSaliva_subset = countsSaliva[salivaSampleSelection_new,]
sampleInfoSaliva_subset = sampleInfoSaliva[salivaSampleSelection_new,]

# Save split data
write.table(countsTongue_subset, "./Data/20240503_UNOISE_new/tongueCounts.csv", sep=",", row.names=FALSE, col.names=FALSE)
write.table(sampleInfoTongue_subset, "./Data/20240503_UNOISE_new/tongueSampleMeta.csv", sep=",", row.names=FALSE, col.names=FALSE)
write.table(countsSaliva_subset, "./Data/20240503_UNOISE_new/salivaCounts.csv", sep=",", row.names=FALSE, col.names=FALSE)
write.table(sampleInfoSaliva_subset, "./Data/20240503_UNOISE_new/salivaSampleMeta.csv", sep=",", row.names=FALSE, col.names=FALSE)

