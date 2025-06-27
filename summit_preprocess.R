setwd("~/Desktop/DNABERT_2/AR_bert2")
library(tidyverse)
bed <- read.table("VCap_DMSO_6h_AR_narrow_summits_t2t.bed", sep = "\t", header = F)
head(bed)
bed <- bed %>% select(V1, V2, V3, V5)
colnames(bed) <- c("original_seq_name", "start","end", "score")
ref <- read_tsv("t2t_sequence_report.tsv")
head(ref)
ref <- ref %>% select(`RefSeq seq accession`, `UCSC style name`)
head(ref)
colnames(ref) <- c("original_seq_name", "new_seq_name")
df <- merge(bed, ref, by="original_seq_name", all.x = TRUE)
head(df)
df <- df %>% select(new_seq_name,start, end, score)
head(df)
df$start <- df$start - 75
df$end <- df$end + 75
head(df)
dim(df) # 63476     4
summary(df$score)
df <- df %>% arrange(-score)
df <- df[1:5000,]
dim(df) #5000    4
write.table(df, file = "AR_peaks_VCap.bed", sep = "\t", row.names = F,
            col.names = F, quote = F)
# calculate GC contents
tmp <- read.table("AR_pos_peaks_gc_content.txt", sep = "\t", header = F)
head(tmp)
mean(tmp$V6) #0.4218848
tmp_2 <- read.table("random_non_AR_overlapping_regions_gc.txt", sep = "\t", header = F)
dim(tmp_2) # 9998    15
head(tmp_2)
tmp_3 <- tmp_2 %>% filter(V8 >= 0.37) %>% filter(V8 <= 0.47)
dim(tmp_3) # 3584    15
head(tmp_3)
tmp_3 <- tmp_3 %>% select(V1, V2, V3)
write.table(tmp_3, file = "random_non_AR_overlapping_regions_filtered_gc.bed", sep = "\t",
            col.names = F, row.names = F, quote = F)

# extract names and sequences from a fasta file
# awk '/^>/ {if (seq) print name"\t"seq; name=$0; seq=""; next} {seq=seq $0} END {print name"\t"seq}' your_fasta_file.fasta > output_sequences.tsv
pos <- read.table("AR_pos_peaks_sequences.tsv", sep = "\t", header = F)
head(pos)
View(pos)
neg <- read.table("random_non_AR_overlapping_regions_filtered_gc.tsv", sep = "\t", header = F)
View(neg)
pos$V3 <- 1
neg$V3 <- 0
total <- rbind(pos,neg)
dim(total) # 8584     3
View(total)

set.seed(123) 
shuffled_df <- total[sample(nrow(total)), ]
View(shuffled_df)
write.table(shuffled_df, file = "dataset.tsv", sep = "\t", row.names = F, col.names = F, quote = F)

shuffled_df <- shuffled_df %>% select(V2, V3)
colnames(shuffled_df) <- c("sequence", "label")
n <- nrow(shuffled_df)
train_size <- floor(0.7 * n)
val_size <- floor(0.15 * n)
train_df <- shuffled_df[1:train_size, ]
val_df <- shuffled_df[(train_size + 1):(train_size + val_size), ]
test_df <- shuffled_df[(train_size + val_size + 1):n, ]
head(train_df)
write.csv(train_df, "train.csv", row.names = FALSE)
write.csv(val_df, "dev.csv", row.names = FALSE)
write.csv(test_df, "test.csv", row.names = FALSE)

# Read in TR sequences file
test <- read.table("TRs_w_padding_sequences.tsv", sep = "\t", header = F)
head(test)
test <- test %>% select(V2)
colnames(test) <- "sequence"
head(test)
dim(test)
sub <- test[5001:12372,]
write.csv(test, "TRs_inquiry.csv", row.names = FALSE)
write.csv(sub, "test.csv", row.names = FALSE)

# Compile TRs_inquiry results csv files
r1 <- read.csv("prediction_results/TRs_inquiry_results_5000.csv")
head(r1)
r2 <- read.csv("prediction_results/TRs_inquiry_results_rest.csv")
r <- rbind(r1,r2)
dim(r) #12372     2
sum(r$prediction == 0) #10101
sum(r$prediction == 1) # 2271




