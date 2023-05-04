series_matrix_data <- read.table("htseq_counts_SRR_labeled.txt",sep="\t",header=TRUE,nrows=57011,quote="\"")
numerical_data <- series_matrix_data[,2:13]

head(numerical_data)
data_normalized <- scale(numerical_data)
head(data_normalized)
corr_matrix <- cor(data_normalized)

library(ggcorrplot)
library(ggplot2)
library(GGally)

ggcorr(corr_matrix, label=TRUE, hjust = 0.75, size = 2)
