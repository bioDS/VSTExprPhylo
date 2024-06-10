library("phyloRNA")
library(Seurat)
library(stringr)

##Add any scripts or files that are needed for this script in the repo and also in the description of analyse.sh.

#enter root directory here
root_dir="../VSTExprPhylo"


# enter correct paths to the phyloRNAanalysis scripts
import::from("../phyloRNAanalysis/src/expr.r", process_expression)
import::from("../phyloRNAanalysis/src/utils.r", read_table, write_table, table2fasta)

import::from("../../scripts/functions.R", rename_categories)
import::from("../../scripts/functions.R", remove_zeroexrp_cells_genes)

setwd(file.path(root_dir, "patient2/Tumour2"))

#enter path to aggregated data here
path_to_tumour2_data ="../outs/filtered_feature_bc_matrix"

data_t2=Read10X(data.dir=path_to_tumour2_data)

data_t2 <- as.data.frame(data_t2)

outdir = file.path("preprocess", "expr")

hdi = c(0.6, 0.9)

prefix="all"

if (!file.exists(outdir)) {
  dir.create(outdir)
}

result = list(
  intervals = file.path(outdir, paste(prefix, "intervals", "txt", sep=".")),
  discretized = file.path(outdir, paste(prefix, "discretized", "txt", sep="."))
)


data = process_expression(
    data_t2,
    hdi = hdi,
    minGene=0,
    minUMI=0,
    trim = FALSE, normalize = FALSE,
    intervals = result$intervals
)

write_table(data, result$discretized)
data = read_table(result$discretized)


## normalized and discretized data

#copied from phyloRNA but without log1p
expr_normalize <- function(data, scale_factor) {
  data[data == 0] = NA
  totals = colSums(data, na.rm=TRUE)
  data = t(t(data)/totals) * scale_factor
  data
}

data_t2_norm <- expr_normalize(data_t2, 10000)

outdir_norm = file.path("preprocess", "expr_norm")

if (!file.exists(outdir_norm)) {
  dir.create(outdir_norm)
}

result_norm = list(
  intervals = file.path(outdir_norm, paste(prefix, "intervals", "txt", sep=".")),
  discretized = file.path(outdir_norm, paste(prefix, "discretized", "txt", sep="."))
)

data_norm_disc = process_expression(
  data_t2_norm, hdi = hdi,
  minGene = 0, minUMI = 0,
  trim = FALSE, normalize = FALSE,
  intervals = result_norm$intervals
)

write_table(data_norm_disc, result_norm$discretized)


## find HVG

data_t2_nz <- data_t2[colSums(data_t2)>0]
data_t2_nz <- data_t2_nz[rowSums(data_t2_nz)>0,]

data_t2_nz_na <- data_t2_nz
data_t2_nz_na[data_t2_nz_na==0] <- NA
gene_var <- apply(data_t2_nz_na, 1, function(x) var(x, na.rm=T))

hvg2_names <- names(which(gene_var>2)) # 198 in total

write(hvg2_names, "../../HVGs/p2_tumour2_HVGs_var_less_2.txt")



## Testset sample barcodes:

testset_new_barcodes <- scan("testset_new_barcodes.txt", what=character())
testset_barcodes <- str_replace(testset_new_barcodes, "-mg1|-mg2|-ms", "-1")





## Testset all

test_data_all <- data[, testset_barcodes]

names(test_data_all) <- testset_new_barcodes

gene_names <- row.names(test_data_all)
test_data_all[test_data_all=="-"] <- 0
test_data_all <- as.data.frame(lapply(test_data_all, as.numeric))
names(test_data_all) <- str_replace(names(test_data_all), "\\.", "-")
row.names(test_data_all) <- gene_names

test_data_all <- remove_zeroexrp_cells_genes(test_data_all, "P2 testset all" )

test_data_all[test_data_all==1] <- 2

print("Checking whether categories need to be shifted in P2 testset all sample")
test_data_all <- rename_categories(test_data_all, 6)

write_fasta(tab2seq(test_data_all, margin=2), file="fasta/testset_zero.fasta")


## Testset HVG

test_data_hvg2 <- test_data_all[row.names(test_data_all) %in% hvg2_names, ]

test_data_hvg2 <- remove_zeroexrp_cells_genes(test_data_hvg2, "P2 testset HVG" )

print("Checking whether categories need to be shifted in P2 testset HVG sample")
test_data_hvg2 <- rename_categories(test_data_hvg2, 6)

write_fasta(tab2seq(test_data_hvg2, margin=2), file="fasta/testset_zero_hvg.fasta")


#Testset continuous HVG

data_test_cont_hvg2 <- data_t2[row.names(data_t2) %in% hvg2_names,testset_barcodes]
names(data_test_cont_hvg2) <- testset_new_barcodes

data_test_cont_hvg2  <- remove_zeroexrp_cells_genes(data_test_cont_hvg2 , "P2 testset continuous HVG" )


if (!file.exists("continuous")) {
  dir.create("continuous")
}

#create keys:
cat(colnames(data_test_cont_hvg2), file="continuous/keys_testset_hvg.txt", sep=" ")

# create taxonset:
cat(paste("<taxon id=\"", colnames(data_test_cont_hvg2), "\" spec=\"Taxon\"/>", sep=""), file="continuous/taxonset_testset_hvg.txt", sep="\n")

# create counts:
cat(unlist(data_test_cont_hvg2), file="continuous/counts_testset_hvg.txt", sep = " ")


## Testset normalized HVG

data_norm_disc = read_table(result_norm$discretized)

data_testset_norm <- data_norm_disc[testset_barcodes]
names(data_testset_norm) <- testset_new_barcodes

gene_names <- row.names(data_testset_norm)
data_testset_norm[data_testset_norm=="-"] <- 0
data_testset_norm <- as.data.frame(lapply(data_testset_norm, as.numeric))
names(data_testset_norm) <- str_replace(names(data_testset_norm), "\\.", "-")
row.names(data_testset_norm) <- gene_names

data_testset_norm_hvg <- data_testset_norm[row.names(data_testset_norm) %in% hvg2_names, ]
data_testset_norm_hvg <- remove_zeroexrp_cells_genes(data_testset_norm_hvg, "P2 testset normalised HVG" )

print("Checking whether categories need to be shifted in P2 testset normalised HVG sample")
data_testset_norm_hvg <- rename_categories(data_testset_norm_hvg, 6)

write_fasta(tab2seq(data_testset_norm_hvg, margin=2), file="fasta/testset_zero_norm_hvg.fasta")
