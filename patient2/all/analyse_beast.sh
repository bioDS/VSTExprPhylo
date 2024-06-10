# Creating fasta files for beast analyses
# This script requiers
# 1. R
# 2. PhyloRNA, Seurat, stringr packages to be installed for R
# 3. PhyloRNAanalysis repo https://github.com/bioDS/phyloRNAanalysis + adjust paths to expr.r, utils.r and filter.r in make_fasta.R script
# 4. Adjust all paths in make_fasta_aggregated.R

mkdir preprocess
mkdir fasta

Rscript make_fasta_aggregated.R  > out

# after running this script based on the output, insert frequencies and number of genes manually to make_xml.sh and run it. It will creat BEAST directory with all beast analyses
# make_xml.sh

# Then beast analyses were run manually. After checking for convergence, MCC trees were calculated using the following command
# path_to_tree_annotator=""
# java -Xmx30G -cp $path_to_tree_annotator -burnin 10 -heights keep $basename".trees" "MCCkeep_"$basename".tree"
