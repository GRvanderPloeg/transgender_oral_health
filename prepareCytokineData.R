library(tidyverse)

df = read.csv("./Data/results_multiplex_sorted_2.csv") %>% as_tibble()

# Repair timepoints
df = df %>% mutate(newTimepoint = 0)
df[df$Timepoint == "Baseline", "newTimepoint"] = 0
df[df$Timepoint == "3months", "newTimepoint"] = 3
df[df$Timepoint == "1year", "newTimepoint"] = 12

# Separate into numeric data and sample metadata
sampleMeta = df %>% select(Participat_code, GenderID, newTimepoint)
numericData = df %>% select(-X.Sample_ID, -Participat_code, -GenderID, -Timepoint, -newTimepoint, -Grouped_cytokines, -Ratio_IL_1b_IL_1RA)
features = colnames(numericData)

# Combine sampleMeta with hormone levels for easy access in MATLAB
hormones = read.csv("./Data/20240503_UNOISE_new/salivaSampleMeta.csv", header=FALSE) %>% as_tibble() %>% select(V23, V24, V20, V17)
colnames(hormones) = c("Participat_code", "newTimepoint", "Testosterone", "Estradiol")
sampleMeta = sampleMeta %>% left_join(hormones)

# Save
write.table(sampleMeta, "./Data/20241209_cytokines_sampleMeta.csv", row.names=FALSE, col.names=FALSE)
write.table(features, "./Data/20241209_cytokines_featureMeta.csv", row.names=FALSE, col.names=FALSE)
write.table(numericData, "./Data/20241209_cytokines.csv", row.names=FALSE, col.names=FALSE)
