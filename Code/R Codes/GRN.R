if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GENIE3")

install.packages('doRNG')
install.packages('doParallel')
install.packages('tidyverse')
library(doParallel)
library(doRNG)
# Load the necessary library
library(GENIE3)


####################
tf <- read.table('TF_names_v_1.01.txt', header = FALSE, sep = "\t")
tf <- tf[[1]]
target <- read.table('genes_exp_response_associated_union_merged.tsv', header = TRUE, sep = "\t")
target <- gsub('-','.',target[[1]])
expression_data_responders <- read.table('merged_res.tsv', header = TRUE, sep = "\t", row.names = 1)
expression_data_nonresponders <- read.table('merged_nonres.tsv', header = TRUE, sep = "\t", row.names = 1)
genes <- colnames(expression_data_responders)
gene_tf <- intersect(tf, genes)
# non_responders = t(subset(expression_data, ICB.Response == 0, select = -c(dataset, ICB.Response)))
# responders = t(subset(expression_data, ICB.Response == 1, select = -c(dataset, ICB.Response)))

#expression_1 is the expression file for the responders
# u should have expression_2 as well ready and prepared, such that when u run for the responders u just change the file name
# it should have genes as rows and samples as columns (opposite to the usual)
# this is the only way it would be used with GENIE3.
# So... you'd do your pre-processing stuff first (in python or whatever) and then pass a pre-processed file and run the code

expression_data_responders = t(expression_data_responders)
expression_data_nonresponders = t(expression_data_nonresponders)

# Run GENIE3 on the matrix (this is the step where you'll carefully adjust the parameters :D)

non_responders_gene_network <- GENIE3(expression_data_nonresponders, regulators = gene_tf, targets = target, K = 'all', nCores = 16)

responders_gene_network <- GENIE3(expression_data_responders, regulators = gene_tf, targets = target, K = 'all', nCores = 16)

write.table(non_responders_gene_network, file = 'non_responders_gene_network.tsv', sep = "\t", row.names = TRUE, col.names = TRUE)
write.table(responders_gene_network, file = 'responders_gene_network.tsv', sep = "\t", row.names = TRUE, col.names = TRUE)

non_responders_linkList <- getLinkList(non_responders_gene_network, threshold = mean(non_responders_gene_network)+3*sd(non_responders_gene_network))#, threshold = mean(gene_network)+10*sd(gene_network) )
responders_linkList <- getLinkList(responders_gene_network, threshold = mean(responders_gene_network)+3*sd(responders_gene_network))

write.table(non_responders_linkList, file = 'nonresponders_linkList.tsv', row.names = FALSE, col.names = TRUE, sep='\t')
write.table(responders_linkList, file = 'responders_linkList.tsv', row.names = FALSE, col.names = TRUE, sep='\t')
