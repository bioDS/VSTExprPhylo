# Creating fasta files for beast analyses
# This script requiers
# 1. file tumour2.h5 to be in path_to_h5s folder
# 2. R
# 3. PhyloRNA, Seurat, stringr package to be installed for R
# 4. PhyloRNAanalysis repo https://github.com/bioDS/phyloRNAanalysis + adjust paths to expr.r, utils.r and filter.r in make_fasta.R script
# 5. Adjust all paths in make_fasta_testset.R

mkdir fasta
mkdir preprocess

Rscript make_fasta_testset.R > out

# after running this script based on the output, insert frequencies and number of genes manually to make_xml.sh and run it. It will creat BEAST directory with all beast analyses
# make_xml.sh

# Then beast analyses were run manually. After checking for convergence, MCC trees were calculated using the following command
# path_to_tree_annotator=""
# java -Xmx30G -cp $path_to_tree_annotator -burnin 10 -heights keep $basename".trees" "MCCkeep_"$basename".tree"
