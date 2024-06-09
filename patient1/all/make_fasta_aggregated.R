library("phyloRNA")
library("Seurat")
library("stringr")

#enter root directory here
root_dir="../VSTExprPhylo"

#enter the correct path to the phyloRNAanalysis
import::from("../phyloRNAanalysis/src/expr.r", preprocess_expression, process_expression)
import::from("../phyloRNAanalysis/src/utils.r", read_table, write_table, table2fasta)
import::from("../phyloRNAanalysis/src/filter.r", density_filtering)

import::from("../../scripts/functions.R", rename_categories)
import::from("../../scripts/functions.R", remove_zeroexrp_cells_genes)

#enter path to aggregated data here
path_to_aggr_data = "../outs/filtered_feature_bc_matrix"

data_all=Read10X(data.dir=path_to_aggr_data)

data_all <- as.data.frame(data_all)

data_p1_all <- data_all[ grepl("-5", names(data_all), fixed = TRUE) | grepl("-6", names(data_all), fixed = TRUE) | grepl("-7", names(data_all), fixed = TRUE) | grepl("-8", names(data_all), fixed = TRUE)]

outdir = file.path("preprocess", "expr_acrosspt")
hdi = c(0.6, 0.9)

prefix="all"

if (!file.exists(outdir)) {
  dir.create(outdir)
}

result = list(
  intervals = file.path(outdir, paste(prefix, "intervals", "txt", sep=".")),
  discretized = file.path(outdir, paste(prefix, "discretized", "txt", sep="."))
)

data_acrosspt = process_expression(
  data_p1_all, hdi = hdi,
  minGene = 0, minUMI = 0,
  trim = FALSE, normalize = FALSE,
  intervals = result$intervals
)

write_table(data_acrosspt, result$discretized)
data_acrosspt = read_table("preprocess/expr_acrosspt/all.discretized.txt")


## normalized and discretized data

#copied from phyloRNA but without log1p
expr_normalize <- function(data, scale_factor) {
  data[data == 0] = NA
  totals = colSums(data, na.rm=TRUE)
  data = t(t(data)/totals) * scale_factor
  data
}

data_p1_norm <- expr_normalize(data_p1_all, 10000)

outdir_norm = file.path("preprocess", "expr_acrosspt_norm")

if (!file.exists(outdir_norm)) {
  dir.create(outdir_norm)
}

result_norm = list(
  intervals = file.path(outdir_norm, paste(prefix, "intervals", "txt", sep=".")),
  discretized = file.path(outdir_norm, paste(prefix, "discretized", "txt", sep="."))
)

data_norm_disc = process_expression(
  data_p1_norm, hdi = hdi,
  minGene = 0, minUMI = 0,
  trim = FALSE, normalize = FALSE,
  intervals = result_norm$intervals
)

write_table(data_norm_disc, result_norm$discretized)
data_norm_disc = read_table(result_norm$discretized)


## find HVG

data_p1_nz <- data_p1_all[colSums(data_p1_all)>0]
data_p1_nz <- data_p1_nz[rowSums(data_p1_nz)>0,]

data_p1_nz_na <- data_p1_nz
data_p1_nz_na[data_p1_nz_na==0] <- NA
gene_var <- apply(data_p1_nz_na, 1, function(x) var(x, na.rm=T))

hvg_names <- names(which(gene_var>10))

print(paste("The number of highly variable genes for patient 1 is:", as.character(length(hvg_names)), sep=" "))

write(hvg_names, file.path(root_dir, "HVGs", "p1_HVGs.txt"))

## Similar quality sample
# The following commented piece of code uses highest quality normal spots and total counts
# to randomly sample similar quality spots that are then written in
# stc_barcodes.txt and stc_new_barcodes.txt files

# data_p1_normal1 <- data_all[ grepl("-5", names(data_all), fixed = TRUE)]
# data_p1_normal2 <- data_all[ grepl("-6", names(data_all), fixed = TRUE)]
# data_p1_tumour1 <- data_all[ grepl("-7", names(data_all), fixed = TRUE)]
# data_p1_tumour2 <- data_all[ grepl("-8", names(data_all), fixed = TRUE)]
#
# # adjust path here
# annotation_file_path=""
# annotation <- read.csv(annotation_file_path, col.names=c("Barcode", "Annotation"))
#
# barcodes_mg <- annotation$Barcode[annotation$Annotation=="Malignant Glands"]
# barcodes_bce <- annotation$Barcode[annotation$Annotation=="Benign Colon Epithelium"]
# 
# barcodes_normal1_b <- intersect(barcodes_bce, names(data_p1_normal1))
# barcodes_normal2_b <- intersect(barcodes_bce, names(data_p1_normal2))
# barcodes_tumour1_b <- intersect(barcodes_bce, names(data_p1_tumour1))
# barcodes_tumour1_m <- intersect(barcodes_mg, names(data_p1_tumour1))
# barcodes_tumour2_b <- intersect(barcodes_bce, names(data_p1_tumour2))
# barcodes_tumour2_m <- intersect(barcodes_mg, names(data_p1_tumour2))
# 
# 
# total_count_p1_normal1_b <- colSums(data_p1_normal1[barcodes_normal1_b])
# total_count_p1_normal2_b <- colSums(data_p1_normal2[barcodes_normal2_b])
# total_count_p1_tumour1_b <- colSums(data_p1_tumour1[barcodes_tumour1_b])
# total_count_p1_tumour1_m <- colSums(data_p1_tumour1[barcodes_tumour1_m])
# total_count_p1_tumour2_b <- colSums(data_p1_tumour2[barcodes_tumour2_b])
# total_count_p1_tumour2_m <- colSums(data_p1_tumour2[barcodes_tumour2_m])
# 
# 
# upper=quantile(density(total_count_p1_normal2_b), probs=0.95)
# lower=1500
# 
# sample_stc <- function(x, upper, lower, num) {
#   x_sample_from <- sort(x[x < upper & x > lower])
#   if (num < length(x_sample_from)) {
#       names <- names(x_sample_from[sample(1:length(x_sample_from), size=num)])
#   } else {
#       names <- names(x_sample_from)
#   }
#   names
#  
# }
# 
# stc_n1_b_barcodes <- sample_stc(total_count_p1_normal1_b, upper, lower, 5)
# stc_n2_b_barcodes <- sample_stc(total_count_p1_normal2_b, upper, lower, 5)
# stc_t1_b_barcodes <- sample_stc(total_count_p1_tumour1_b, upper, lower, 5)
# stc_t2_b_barcodes <- sample_stc(total_count_p1_tumour2_b, upper, lower, 5)
# stc_t1_m_barcodes <- sample_stc(total_count_p1_tumour1_m, upper, lower, 5)
# stc_t2_m_barcodes <- sample_stc(total_count_p1_tumour2_m, upper, lower, 5)
# 
# 
# stc_barcodes <- c(stc_n1_b_barcodes, stc_n2_b_barcodes, stc_t1_b_barcodes, stc_t2_b_barcodes, stc_t1_m_barcodes, stc_t2_m_barcodes)
# write(stc_barcodes, "stc_barcodes.txt")
# 
stc_barcodes <- scan("stc_barcodes.txt", what=character())
 
# stc_new_barcodes <- c(str_replace(stc_n1_b_barcodes, "-5", "-normal1_b"), str_replace(stc_n2_b_barcodes, "-6", "-normal2_b"),
#                              str_replace(stc_t1_b_barcodes, "-7", "-tumour1_b"), str_replace(stc_t2_b_barcodes, "-8", "-tumour2_b"),
#                              str_replace(stc_t1_m_barcodes, "-7", "-tumour1_m"), str_replace(stc_t2_m_barcodes, "-8", "-tumour2_m"))
# write(stc_new_barcodes, "stc_new_barcodes.txt")

stc_new_barcodes <- scan("stc_new_barcodes.txt", what=character())



## Similar all

data_stc <- data_acrosspt[stc_barcodes]
names(data_stc) <- stc_new_barcodes

gene_names <- row.names(data_stc)
data_stc[data_stc=="-"] <- 0
data_stc <- as.data.frame(lapply(data_stc, as.numeric))
names(data_stc) <- str_replace(names(data_stc), "\\.", "-")
row.names(data_stc) <- gene_names


data_stc <- remove_zeroexrp_cells_genes(data_stc, "P1 stc all" )

print("Checking whether categories need to be shifted in P1 stc all sample")
data_stc <- rename_categories(data_stc, 6)

write_fasta(tab2seq(data_stc, margin=2), file="fasta/stc_zero.fasta")



## Similar HVG

data_stc_hvg <- data_stc[row.names(data_stc) %in% hvg_names, ]
data_stc_hvg <- remove_zeroexrp_cells_genes(data_stc_hvg, "P1 stc HVG" )

print("Checking whether categories need to be shifted in P1 stc HVG sample")
data_stc_hvg <- rename_categories(data_stc_hvg, 6)

write_fasta(tab2seq(data_stc_hvg, margin=2), file="fasta/stc_zero_hvg.fasta")



## Similar continuous HVG

data_stc_cont_hvg <- data_p1_all[hvg_names,stc_barcodes]
names(data_stc_cont_hvg ) <- stc_new_barcodes

data_stc_cont_hvg  <- remove_zeroexrp_cells_genes(data_stc_cont_hvg , "P1 stc continuous HVG" )


#create keys:
cat(colnames(data_stc_cont_hvg), file="continuous/keys_stc_hvg.txt", sep=" ")

# create taxonset:
cat(paste("<taxon id=\"", colnames(data_stc_cont_hvg), "\" spec=\"Taxon\"/>", sep=""), file="continuous/taxonset_stc_hvg.txt", sep="\n")

# create counts:
cat(unlist(data_stc_cont_hvg), file="continuous/counts_stc_hvg.txt", sep = " ")



## Similar normalized HVG

data_acrosspt_norm = read_table("preprocess/expr_acrosspt_norm/all.discretized.txt")

data_stc_norm <- data_acrosspt_norm[stc_barcodes]
names(data_stc_norm) <- stc_new_barcodes

gene_names <- row.names(data_stc_norm)
data_stc_norm[data_stc_norm=="-"] <- 0
data_stc_norm <- as.data.frame(lapply(data_stc_norm, as.numeric))
names(data_stc_norm) <- str_replace(names(data_stc_norm), "\\.", "-")
row.names(data_stc_norm) <- gene_names

data_stc_norm_hvg <- data_stc_norm[row.names(data_stc_norm) %in% hvg_names, ]
data_stc_norm_hvg <- remove_zeroexrp_cells_genes(data_stc_norm_hvg, "P1 stc normalised HVG" )

print("Checking whether categories need to be shifted in P1 stc normalised HVG sample")
data_stc_norm_hvg <- rename_categories(data_stc_norm_hvg, 6)

write_fasta(tab2seq(data_stc_norm_hvg, margin=2), file="fasta/stc_zero_norm_hvg.fasta")
