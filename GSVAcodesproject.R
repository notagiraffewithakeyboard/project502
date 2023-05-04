series_matrix_data <- read.table("htseq_counts_SRR_labeled.txt", sep = "\t", skip = 0, header = TRUE, nrows = 57010, quote = "\"", row.names = 1)

################# loading

hallmark_gene_sets <- msigdbr::msigdbr(
  species = "mouse", # Can change this to what species you need
  category = "H" # Only hallmark gene sets
)

hallmarks_list <- split(
  hallmark_gene_sets$gene_symbol, # The genes we want split into pathways
  hallmark_gene_sets$gs_name # The pathways made as the higher levels of the list
)


library(tibble)
library(GSVA)
library(Biobase)
library(GSEABase)
library(GSVAdata)
library(dplyr)
library(biomaRt)
library(AnnotationDbi)
library(RColorBrewer)
library(dplyr)
library(org.Mm.eg.db)
library(pheatmap)


##alt usefulsets,not used 
data(c2BroadSets)
class(c2BroadSets)
c2BroadSets
data_matrix <- as.matrix(series_matrix_data)
data2 <- as.data.frame(c2BroadSets)


#########


mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
gene_info <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"), filters = "ensembl_gene_id", values = ensembl_ids, mart = mart)


###fix number of rows in original data###

test <- tolower(rownames(series_matrix_data))
series_matrix_data <- unique(series_matrix_data)

series_matrix_data <- series_matrix_data[complete.cases(series_matrix_data), ]
series_matrix_data <- series_matrix_data[!duplicated(rownames(series_matrix_data)), ]

series_matrix_data <- distinct(series_matrix_data, .keep_all = TRUE)


##add IDS

ensembl_ids <- rownames(series_matrix_data)
ensembl <- useMart("ensembl")
ensembl <- useDataset("mmusculus_gene_ensembl", mart = ensembl)
ensembl_mapping <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"), 
                         filters = "ensembl_gene_id", 
                         values = ensembl_ids, 
                         mart = ensembl)

print(ensembl_mapping)

org <- org.Mm.eg.db


idx <- match(ensembl_ids, ensembl_mapping$ensembl_gene_id)
gene_names <- ensembl_mapping$external_gene_name[idx]


gene_names


# Create unique IDs for duplicate gene names

gene_freq <- table(gene_names)
dup_ids <- rep(NA, length(gene_names))
dup_names <- names(gene_freq[gene_freq > 1])

for (i in dup_names) {
  idx <- which(gene_names == i)
  dup_ids[idx] <- seq_along(idx)
}

dup_ids[is.na(dup_ids)] <- 1

rownames(series_matrix_data) <- ifelse(dup_ids == 1, gene_names, paste(gene_names, dup_ids, sep="_"))


####at last####

# Convert the data frame to a matrix
expr_mat <- as.matrix(series_matrix_data)

# Unused gene set, here for GMT alternative
gene_sets <- getGmt("GOBP_MUSCLE_STRUCTURE_DEVELOPMENT.v2023.1.Mm.GMT")
gene_sets_list <- as.list(gene_sets)

gsva_result <- gsva(expr_mat, gene_sets_list)


#data used

gsva_results1 <- gsva(
  expr_mat,
  hallmarks_list,
  method = "gsva",
  # Appropriate for our vst transformed data
  kcdf = "Gaussian",
  # Minimum gene set size
  min.sz = 15,
  # Maximum gene set size
  max.sz = 500,
  # Compute Gaussian-distributed scores
  mx.diff = TRUE,
  # Don't print out the progress bar
  verbose = FALSE
)

pheatmap(t(gsva_results1), scale="column", col = terrain.colors(256))