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

setwd(file.path(root_dir, "patient2/all"))

#enter path to aggregated data here
path_to_aggr_data = "/outs/filtered_feature_bc_matrix"

data_all=Read10X(data.dir=path_to_aggr_data)

data_all <- as.data.frame(data_all)

data_419_all <- data_all[ grepl("-1", names(data_all), fixed = TRUE) | grepl("-2", names(data_all), fixed = TRUE) | grepl("-3", names(data_all), fixed = TRUE) | grepl("-4", names(data_all), fixed = TRUE)]

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
  data_419_all, hdi = hdi,
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

data_419_norm <- expr_normalize(data_419_all, 10000)

outdir_norm = file.path("preprocess", "expr_acrosspt_norm")

if (!file.exists(outdir_norm)) {
  dir.create(outdir_norm)
}

result_norm = list(
  intervals = file.path(outdir_norm, paste(prefix, "intervals", "txt", sep=".")),
  discretized = file.path(outdir_norm, paste(prefix, "discretized", "txt", sep="."))
)

data_norm_disc = process_expression(
  data_419_norm, hdi = hdi,
  minGene = 0, minUMI = 0,
  trim = FALSE, normalize = FALSE,
  intervals = result_norm$intervals
)

write_table(data_norm_disc, result_norm$discretized)
data_norm_disc = read_table(result_norm$discretized)


## find HVG

data_419_nz <- data_419_all[colSums(data_419_all)>0]
data_419_nz <- data_419_nz[rowSums(data_419_nz)>0,]

data_419_nz_na <- data_419_nz
data_419_nz_na[data_419_nz_na==0] <- NA
gene_var <- apply(data_419_nz_na, 1, function(x) var(x, na.rm=T))

hvg_names <- names(which(gene_var>10)) # a total of 115 high variance genes

print(paste("The number of highly variable genes for patient 2 is:", as.character(length(hvg_names)), sep=" "))

write(hvg_names, file.path(root_dir, "HVGs", "p2_HVGs.txt"))


## Selected sample
# The following commented piece of code converts manually selected spots with annotation from selected_annotation.csv file to 
# selected_barcodes.txt and selected_new_barcodes.txt files  

# annotation_selected <- read.csv("selected_annotation.csv", col.names=c("Barcode", "Annotation"))
# 
# barcodes_normal <- annotation_selected$Barcode[annotation_selected$Annotation=="Epithelium"]
# barcodes_tumour1_m <- annotation_selected$Barcode[annotation_selected$Annotation=="Malignant_gland_t1"]
# barcodes_tumour1_b <- annotation_selected$Barcode[annotation_selected$Annotation=="Epithelium_t1"]
# barcodes_tumour2_m <- annotation_selected$Barcode[annotation_selected$Annotation=="Malignant_gland_t2"]
# barcodes_tumour2_b <- annotation_selected$Barcode[annotation_selected$Annotation=="Epithelium_t2"]
# barcodes_ovary_m <- annotation_selected$Barcode[annotation_selected$Annotation=="Malignant_gland_ovary"]
# 
# selected_barcodes <- c(barcodes_normal, barcodes_tumour1_m, barcodes_tumour1_b, barcodes_tumour2_m, barcodes_tumour2_b, barcodes_ovary_m)
# write(selected_barcodes, "selected_barcodes.txt")

selected_barcodes <- scan ("selected_barcodes.txt", what=character())
# 
# selected_new_barcodes <- c(str_replace(barcodes_normal, "-3", "-normal_b"), str_replace(barcodes_tumour1_m, "-1", "-tumour1_m"), str_replace(barcodes_tumour1_b, "-1", "-tumour1_b"), str_replace(barcodes_tumour2_m, "-2", "-tumour2_m"), str_replace(barcodes_tumour2_b, "-2", "-tumour2_b"), str_replace(barcodes_ovary_m, "-4", "-ovary_m"))
# write(selected_new_barcodes, "selected_new_barcodes.txt")

selected_new_barcodes <- scan ("selected_new_barcodes.txt", what=character())



## Selected all

data_selected <- data_acrosspt[selected_barcodes]
names(data_selected) <- selected_new_barcodes

gene_names <- row.names(data_selected)
data_selected[data_selected=="-"] <- 0
data_selected <- as.data.frame(lapply(data_selected, as.numeric))
names(data_selected) <- str_replace(names(data_selected), "\\.", "-")
row.names(data_selected) <- gene_names


data_selected <- remove_zeroexrp_cells_genes(data_selected, "P2 selected all" )

print("Checking whether categories need to be shifted in P2 selected all sample")
data_selected <- rename_categories(data_selected, 6)

write_fasta(tab2seq(data_selected, margin=2), file="fasta/selected_zero.fasta")



## Selected HVG 

data_selected_hvg <- data_selected[row.names(data_selected) %in% hvg_names, ]
data_selected_hvg <- remove_zeroexrp_cells_genes(data_selected_hvg, "P2 selected HVG" )

print("Checking whether categories need to be shifted in P2 selected HVG sample")
data_selected_hvg <- rename_categories(data_selected_hvg, 6)

write_fasta(tab2seq(data_selected_hvg, margin=2), file="fasta/selected_zero_hvg.fasta")



## Selected continuous HVG

data_selected_cont_hvg <- data_419_all[hvg_names,selected_barcodes]
names(data_selected_cont_hvg ) <- selected_new_barcodes

data_selected_cont_hvg  <- remove_zeroexrp_cells_genes(data_selected_cont_hvg , "P2 selected continuous HVG" )

if (!file.exists("continuous")) {
  dir.create("continuous")
}

#create keys:
cat(colnames(data_selected_cont_hvg), file="continuous/keys_selected_hvg.txt", sep=" ")

# create taxonset:
cat(paste("<taxon id=\"", colnames(data_selected_cont_hvg), "\" spec=\"Taxon\"/>", sep=""), file="continuous/taxonset_selected_hvg.txt", sep="\n")

# create counts:
cat(unlist(data_selected_cont_hvg), file="continuous/counts_selected_hvg.txt", sep = " ")



## Selected normalized HVG

data_acrosspt_norm = read_table("preprocess/expr_acrosspt_norm/all.discretized.txt")

data_selected_norm <- data_acrosspt_norm[selected_barcodes]
names(data_selected_norm) <- selected_new_barcodes

gene_names <- row.names(data_selected_norm)
data_selected_norm[data_selected_norm=="-"] <- 0
data_selected_norm <- as.data.frame(lapply(data_selected_norm, as.numeric))
names(data_selected_norm) <- str_replace(names(data_selected_norm), "\\.", "-")
row.names(data_selected_norm) <- gene_names

data_selected_norm_hvg <- data_selected_norm[row.names(data_selected_norm) %in% hvg_names, ]
data_selected_norm_hvg <- remove_zeroexrp_cells_genes(data_selected_norm_hvg, "P2 selected normalised HVG" )

print("Checking whether categories need to be shifted in P2 selected normalised HVG sample")
data_selected_norm_hvg <- rename_categories(data_selected_norm_hvg, 6)

write_fasta(tab2seq(data_selected_norm_hvg, margin=2), file="fasta/selected_zero_norm_hvg.fasta")




## Highest quality sample
# The following commented piece of code uses pathologist annotation and total counts to select highest quality spots
# that are then written in t5q_barcodes.txt and t5q_new_barcodes.txt

#data_419_normal <- data_all[ grepl("-3", names(data_all), fixed = TRUE)]
#data_419_tumour1 <- data_all[ grepl("-1", names(data_all), fixed = TRUE)]
#data_419_tumour2 <- data_all[ grepl("-2", names(data_all), fixed = TRUE)]

#data_419_ovary <- data_all[ grepl("-4", names(data_all), fixed = TRUE)]

#enter path to the annotation
#annotation <- read.csv(".../Pathologist_guided_annotations.csv", col.names=c("Barcode", "Annotation"))

#barcodes_mg <- annotation$Barcode[annotation$Annotation=="Malignant Glands"]
#barcodes_bce <- annotation$Barcode[annotation$Annotation=="Benign Colon Epithelium"]

#barcodes_tumour1_b <- intersect(barcodes_bce, names(data_419_tumour1))
#barcodes_tumour1_m <- intersect(barcodes_mg, names(data_419_tumour1))
#barcodes_tumour2_b <- intersect(barcodes_bce, names(data_419_tumour2))
#barcodes_tumour2_m <- intersect(barcodes_mg, names(data_419_tumour2))
#barcodes_normal_b <- intersect(barcodes_bce, names(data_419_normal))
#barcodes_ovary_m <- intersect(barcodes_mg, names(data_419_ovary))

#total_count_419_normal_b <- colSums(data_419_normal[barcodes_normal_b])
#total_count_419_tumour1_b <- colSums(data_419_tumour1[barcodes_tumour1_b])
#total_count_419_tumour1_m <- colSums(data_419_tumour1[barcodes_tumour1_m])
#total_count_419_tumour2_b <- colSums(data_419_tumour2[barcodes_tumour2_b])
#total_count_419_tumour2_m <- colSums(data_419_tumour2[barcodes_tumour2_m])
#total_count_419_ovary_m <- colSums(data_419_ovary[barcodes_ovary_m])

#total_count_419_tumour_m <- c(total_count_419_tumour1_m,total_count_419_tumour2_m)
#total_count_419_tumour_b <- c(total_count_419_tumour1_b,total_count_419_tumour2_b)

#normal_419_b_exclude <- names(total_count_419_normal_b[total_count_419_normal_b>quantile(density(total_count_419_normal_b), probs=0.99)])
#tumour_419_b_exclude <- names(total_count_419_tumour_b[total_count_419_tumour_b>quantile(density(total_count_419_tumour_b), probs=0.99)])
#tumour_419_m_exclude <- names(total_count_419_tumour_m[total_count_419_tumour_m>quantile(density(total_count_419_tumour_m), probs=0.99)])
#ovary_419_m_exclude <- names(total_count_419_ovary_m[total_count_419_ovary_m>quantile(density(total_count_419_ovary_m), probs=0.99)])

#total_count_419_normal_b <- total_count_419_normal_b[!names(total_count_419_normal_b) %in% normal_419_b_exclude]
#total_count_419_tumour1_b <- total_count_419_tumour1_b[!names(total_count_419_tumour1_b) %in% tumour_419_b_exclude]
#total_count_419_tumour2_b <- total_count_419_tumour2_b[!names(total_count_419_tumour2_b) %in% tumour_419_b_exclude]
#total_count_419_tumour1_m <- total_count_419_tumour1_m[!names(total_count_419_tumour1_m) %in% tumour_419_m_exclude]
#total_count_419_tumour2_m <- total_count_419_tumour2_m[!names(total_count_419_tumour2_m) %in% tumour_419_m_exclude]
#total_count_419_ovary_m <- total_count_419_ovary_m[!names(total_count_419_ovary_m) %in% ovary_419_m_exclude]

#t5q_normal_b_barcodes <- names(sort(total_count_419_normal_b, decreasing=T)[1:5])
#t5q_tumour1_b_barcodes <- names(sort(total_count_419_tumour1_b, decreasing=T)[1:5])
#t5q_tumour2_b_barcodes <- names(sort(total_count_419_tumour2_b, decreasing=T)[1:5])
#t5q_tumour1_m_barcodes <- names(sort(total_count_419_tumour1_m, decreasing=T)[1:5])
#t5q_tumour2_m_barcodes <- names(sort(total_count_419_tumour2_m, decreasing=T)[1:5])
#t5q_ovary_m_barcodes <- names(sort(total_count_419_ovary_m, decreasing=T)[1:5])

#t5q_barcodes <- c(t5q_normal_b_barcodes, t5q_tumour1_b_barcodes, t5q_tumour2_b_barcodes, t5q_tumour1_m_barcodes,
#                  t5q_tumour2_m_barcodes, t5q_ovary_m_barcodes)

#write(t5q_barcodes, "t5q_barcodes.txt")
t5q_barcodes <- scan("t5q_barcodes.txt", what=character())

#t5q_new_barcodes <- c(str_replace(t5q_normal_b_barcodes, "-3", "-normal_b"),
#                      str_replace(t5q_tumour1_b_barcodes, "-1", "-tumour1_b"),
#                      str_replace(t5q_tumour2_b_barcodes, "-2", "-tumour2_b"),
#                      str_replace(t5q_tumour1_m_barcodes, "-1", "-tumour1_m"),
#                      str_replace(t5q_tumour2_m_barcodes, "-2", "-tumour2_m"),
#                      str_replace(t5q_ovary_m_barcodes, "-4", "-ovary_m"))

#write(t5q_new_barcodes, "t5q_new_barcodes.txt")
t5q_new_barcodes <- scan("t5q_new_barcodes.txt", what=character())
                      

## Highest All

data_t5q <- data_acrosspt[t5q_barcodes]
names(data_t5q) <- t5q_new_barcodes

gene_names <- row.names(data_t5q)
data_t5q[data_t5q=="-"] <- 0
data_t5q <- as.data.frame(lapply(data_t5q, as.numeric))
names(data_t5q) <- str_replace(names(data_t5q), "\\.", "-")
row.names(data_t5q) <- gene_names


data_t5q <- remove_zeroexrp_cells_genes(data_t5q, "P2 t5q all" )

print("Checking whether categories need to be shifted in P2 t5q all sample")
data_t5q <- rename_categories(data_t5q, 6)

write_fasta(tab2seq(data_t5q, margin=2), file="fasta/t5q_zero.fasta")



## Highest HVG

data_t5q_hvg <- data_t5q[row.names(data_t5q) %in% hvg_names, ]
data_t5q_hvg <- remove_zeroexrp_cells_genes(data_t5q_hvg, "P2 t5q HVG" )

print("Checking whether categories need to be shifted in P2 t5q HVG sample")
data_t5q_hvg <- rename_categories(data_t5q_hvg, 6)

write_fasta(tab2seq(data_t5q_hvg, margin=2), file="fasta/t5q_zero_hvg.fasta")



## Highest continuous HVG

data_t5q_cont_hvg <- data_419_all[hvg_names,t5q_barcodes]
names(data_t5q_cont_hvg ) <- t5q_new_barcodes

data_t5q_cont_hvg  <- remove_zeroexrp_cells_genes(data_t5q_cont_hvg , "P2 t5q continuous HVG" )

#create keys:
cat(colnames(data_t5q_cont_hvg), file="continuous/keys_t5q_hvg.txt", sep=" ")

# create taxonset:
cat(paste("<taxon id=\"", colnames(data_t5q_cont_hvg), "\" spec=\"Taxon\"/>", sep=""), file="continuous/taxonset_t5q_hvg.txt", sep="\n")

# create counts:
cat(unlist(data_t5q_cont_hvg), file="continuous/counts_t5q_hvg.txt", sep = " ")



## Highest normalized HVG

data_acrosspt_norm = read_table("preprocess/expr_acrosspt_norm/all.discretized.txt")

data_t5q_norm <- data_acrosspt_norm[t5q_barcodes]
names(data_t5q_norm) <- t5q_new_barcodes

gene_names <- row.names(data_t5q_norm)
data_t5q_norm[data_t5q_norm=="-"] <- 0
data_t5q_norm <- as.data.frame(lapply(data_t5q_norm, as.numeric))
names(data_t5q_norm) <- str_replace(names(data_t5q_norm), "\\.", "-")
row.names(data_t5q_norm) <- gene_names

data_t5q_norm_hvg <- data_t5q_norm[row.names(data_t5q_norm) %in% hvg_names, ]
data_t5q_norm_hvg <- remove_zeroexrp_cells_genes(data_t5q_norm_hvg, "P2 t5q normalised HVG" )

print("Checking whether categories need to be shifted in P2 t5q normalised HVG sample")
data_t5q_norm_hvg <- rename_categories(data_t5q_norm_hvg, 6)

write_fasta(tab2seq(data_t5q_norm_hvg, margin=2), file="fasta/t5q_zero_norm_hvg.fasta")





## Similar quality sample
# The following commented piece of code uses highest quality normal spots and total counts
# to randomly sample similar quality spots that are then written in
# stc_barcodes.txt and stc_new_barcodes.txt files 

# upper=max(total_count_419_normal_b[top5quality_normal_b_barcodes])+3
# lower=min(total_count_419_normal_b[top5quality_normal_b_barcodes])-3
# 
# sample_stc <- function(x, upper, lower) {
#   x_sample_from <- sort(x[x < upper & x > lower])
#   names(x_sample_from[c(1, sample(2:(length(x_sample_from)-1), size=3), length(x_sample_from))])
# }
# 
# stc_t1_b_barcodes <- sample_stc(total_count_419_tumour1_b, upper, lower)
# stc_t1_m_barcodes <- sample_stc(total_count_419_tumour1_m, upper, lower)
# stc_t2_b_barcodes <- sample_stc(total_count_419_tumour2_b, upper, lower)
# stc_t2_m_barcodes <- sample_stc(total_count_419_tumour2_m, upper, lower)
# stc_o_m_barcodes <- sample_stc(total_count_419_ovary_m, upper, lower)
# 
# stc_barcodes <- c(t5q_normal_b_barcodes, stc_t1_b_barcodes, stc_t1_m_barcodes, stc_t2_b_barcodes, stc_t2_m_barcodes, stc_o_m_barcodes) 
# write(stc_barcodes, "stc_barcodes.txt")

stc_barcodes <- scan("stc_barcodes.txt", what=character())
# 
# stc_new_barcodes <- c(str_replace(t5q_normal_b_barcodes, "-3", "-normal_b"), str_replace(stc_t1_b_barcodes, "-1", "-tumour1_b"), 
# str_replace(stc_t1_m_barcodes, "-1", "-tumour1_m"), str_replace(stc_t2_b_barcodes, "-2", "-tumour2_b"), 
# str_replace(stc_t2_m_barcodes, "-2", "-tumour2_m"), str_replace(stc_o_m_barcodes, "-4", "-ovary_m")) 
# 
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


data_stc <- remove_zeroexrp_cells_genes(data_stc, "P2 stc all" )

print("Checking whether categories need to be shifted in P2 stc all sample")
data_stc <- rename_categories(data_stc, 6)

write_fasta(tab2seq(data_stc, margin=2), file="fasta/stc_zero.fasta")



## Similar HVG

data_stc_hvg <- data_stc[row.names(data_stc) %in% hvg_names, ]
data_stc_hvg <- remove_zeroexrp_cells_genes(data_stc_hvg, "P2 stc HVG" )

print("Checking whether categories need to be shifted in P2 stc HVG sample")
data_stc_hvg <- rename_categories(data_stc_hvg, 6)

write_fasta(tab2seq(data_stc_hvg, margin=2), file="fasta/stc_zero_hvg.fasta")



## Similar continuous HVG

data_stc_cont_hvg <- data_419_all[hvg_names,stc_barcodes]
names(data_stc_cont_hvg ) <- stc_new_barcodes

data_stc_cont_hvg  <- remove_zeroexrp_cells_genes(data_stc_cont_hvg , "P2 stc continuous HVG" )


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
data_stc_norm_hvg <- remove_zeroexrp_cells_genes(data_stc_norm_hvg, "P2 stc normalised HVG" )

print("Checking whether categories need to be shifted in P2 stc normalised HVG sample")
data_stc_norm_hvg <- rename_categories(data_stc_norm_hvg, 6)

write_fasta(tab2seq(data_stc_norm_hvg, margin=2), file="fasta/stc_zero_norm_hvg.fasta")
