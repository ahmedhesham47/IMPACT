#import libraries

library(dplyr)
library(ggfortify)
library(limma)

# functions to convert tpm to log and vice versa
logtoTPM <- function(x) {return(log2(x+1))}
TPMtolog <- function(x) {return((2^x) - 1)}

# reading the expression file
dfl <- read.csv('merged_exp.tsv', sep = '\t')
dfl$dataset <- ifelse(grepl("Sample", dfl$Sample.Identifier), "Liu", "Ravi")
#extracting the batch info
batch1 <- dfl[['dataset']]

# transposing the dataframe to suit the input of removebatcheffect function
df0 <- t(subset(dfl, select = -c(dataset, Sample.Identifier)))

# log2 transformation of the tpm data
df1 <- as.data.frame(df0) %>% mutate_if(is.numeric, logtoTPM)

# creating the default design matrix with number of columns equal to the expression matrix
design <- matrix(1,ncol(df1),1)

# Remove the batch effect using removeBatchEffect with the design matrix
df2 <- removeBatchEffect(df1, batch = batch1, design = design)

# returning the dataframe back to tpm 
df3 <- t(as.data.frame(df2))
df4 <- as.data.frame(df3) %>% mutate_if(is.numeric, TPMtolog)
df5 <-  cbind(df4, subset(dfl, select = c(dataset, Sample.Identifier)))


# pca plot for BEFORE batch correction
pca_res <- prcomp( subset(dfl, select = -c(dataset, Sample.Identifier) ), scale. = TRUE)
autoplot(pca_res, data = dfl, colour = 'dataset')

# pca plot for AFTER batch correction
pca_resf <- prcomp(df4, scale. = TRUE)
autoplot(pca_resf, data = df5, colour = 'dataset')

write.table(df5, file = "batch_corrected_final.tsv", sep = "\t", row.names = FALSE, col.names = TRUE)




