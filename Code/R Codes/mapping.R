if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("org.Hs.eg.db")

library(dplyr)

library(org.Hs.eg.db)

# __________________________________________________
# Mapping expression data

mrna_tpm_liu <- read.table("Liu/processed_exp_data_liu.tsv", header = TRUE, sep = "\t")
mrna_tpm_ravi <- read.table("Ravi/processed_exp_data_ravi.tsv", header = TRUE, sep="\t")

# Get gene names
gene_names_liu <- colnames(mrna_tpm_liu)
gene_names_ravi <- colnames(mrna_tpm_ravi)


entrez_ids_liu <- mapIds(org.Hs.eg.db, keys = c(gene_names_liu), column = "ENTREZID", keytype="SYMBOL")
entrez_ids_ravi <- mapIds(org.Hs.eg.db, keys = c(gene_names_ravi), column = "ENTREZID", keytype="SYMBOL")


mapping_df_liu <- data.frame(gene_name = gene_names_liu, entrez_id = entrez_ids_liu)

mapping_df_ravi <- data.frame(gene_name = gene_names_ravi, entrez_id = entrez_ids_ravi)

write.table(mapping_df_liu, file = "mapping_liu_exp.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(mapping_df_ravi, file = "mapping_ravi_exp.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# ----------------------

mapped_merged <- read.table("Merged/mapped_merged_exp.tsv", header = TRUE, sep = "\t", check.names = FALSE)
gene_names_merged <- c(colnames(mapped_merged))
merged_gene_symbols <- mapIds(org.Hs.eg.db, keys = gene_names_merged, column = "SYMBOL", keytype = "ENTREZID")


merged_df = data.frame(entrez_id = gene_names_merged, gene_name = merged_gene_symbols)
write.table(merged_df, file = "mapping_from_entrez_to_symbol_for_merged_data_exp.tsv", sep = "\t", quote = FALSE, row.names = FALSE)


# ____________________________________________________________________

# Merging SNP data
snp_liu <- read.table("Liu/processed_snp_data_liu.tsv", header = TRUE, sep = "\t")
snp_ravi <- read.table("Ravi/processed_snp_data_ravi.tsv", header = TRUE, sep="\t")

snp_names_liu <- colnames(snp_liu)
snp_names_liu_replaced <- gsub("_SNP", "", snp_names_liu)
snp_names_ravi <- colnames(snp_ravi)
snp_names_ravi_replaced <- gsub("_SNP", "", snp_names_ravi)


entrez_ids_snp_liu <- mapIds(org.Hs.eg.db, keys = c(snp_names_liu_replaced), column = "ENTREZID", keytype="SYMBOL")
entrez_ids_snp_ravi <- mapIds(org.Hs.eg.db, keys = c(snp_names_ravi_replaced), column = "ENTREZID", keytype="SYMBOL")


mapping_df_snp_liu <- data.frame(gene_name = snp_names_liu, entrez_id = entrez_ids_snp_liu)

mapping_df_snp_ravi <- data.frame(gene_name = snp_names_ravi, entrez_id = entrez_ids_snp_ravi)

write.table(mapping_df_snp_liu, file = "mapping_liu_snp.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(mapping_df_snp_ravi, file = "mapping_ravi_snp.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
# --------------------------------------------------

mapped_merged_snp <- read.table("Merged/mapped_merged_snp.tsv", header = TRUE, sep = "\t", check.names = FALSE)
snp_names_merged <- c(colnames(mapped_merged_snp))
merged_snp_symbols <- mapIds(org.Hs.eg.db, keys = snp_names_merged, column = "SYMBOL", keytype = "ENTREZID")


merged_snp_df = data.frame(entrez_id = snp_names_merged, gene_name = merged_snp_symbols)
write.table(merged_snp_df, file = "mapping_from_entrez_to_symbol_for_merged_data_snp.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# ____________________________________________________________________
# Merging CNA data

cna_liu <- read.table("Liu/processed_cna_data_liu.tsv", header = TRUE, sep = "\t")
cna_ravi <- read.table("Ravi/processed_cna_data_ravi.tsv", header = TRUE, sep="\t")

cna_names_liu <- colnames(cna_liu)
cna_names_liu_replaced <- gsub("_CN", "", cna_names_liu)
cna_names_ravi <- colnames(cna_r)
cna_names_ravi_replaced <- gsub("_CN", "", cna_names_ravi)


entrez_ids_cna_liu <- mapIds(org.Hs.eg.db, keys = c(cna_names_liu_replaced), column = "ENTREZID", keytype="SYMBOL")
entrez_ids_cna_ravi <- mapIds(org.Hs.eg.db, keys = c(cna_names_ravi_replaced), column = "ENTREZID", keytype="SYMBOL")


mapping_df_cna_liu <- data.frame(gene_name = cna_names_liu, entrez_id = entrez_ids_cna_liu)

mapping_df_cna_ravi <- data.frame(gene_name = cna_names_ravi, entrez_id = entrez_ids_cna_ravi)

write.table(mapping_df_cna_liu, file = "mapping_liu_cna.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(mapping_df_cna_ravi, file = "mapping_ravi_cna.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# --------------------------------------------------

mapped_merged_cna <- read.table("Merged/mapped_merged_cna.tsv", header = TRUE, sep = "\t", check.names = FALSE)
cna_names_merged <- c(colnames(mapped_merged_cna))
merged_cna_symbols <- mapIds(org.Hs.eg.db, keys = cna_names_merged, column = "SYMBOL", keytype = "ENTREZID")


merged_cna_df = data.frame(entrez_id = cna_names_merged, gene_name = merged_cna_symbols)
write.table(merged_cna_df, file = "mapping_from_entrez_to_symbol_for_merged_data_cna.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# _________________________________________________________________________

# Adjusting gene names for the validation dataset
validation_mrna = read.table("Validation/data_mrna_seq_rpkm.txt", header = TRUE, sep = "\t")
gene_ids <- as.character(validation_mrna$Entrez_Gene_Id)
gene_symbols <- mapIds(org.Hs.eg.db, keys = c(gene_ids), column = "SYMBOL", keytype="ENTREZID")

merged_df_validation = data.frame(entrez_id = gene_ids, gene_name = gene_symbols)
write.table(merged_df_validation, file = "mapping_from_entrez_to_symbol_for_validation_data.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
