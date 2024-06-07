#P2 50% density analysis
# This script requiers
# 1. files tumour1.h5, tumour2.h5, normal.h5, ovary.h5 to be in path_to_h5s folder
# 2. R
# 3. PhyloRNA package to be installed for R
# 4. PhyloRNAanalysis repo https://github.com/bioDS/phyloRNAanalysis + adjust paths to expr.r, utils.r and filter.r in make_fasta.R script
# 5. iqtree
# 6. adjust path_to_h5s here

#insert path to the folder with .h5 files:
path_to_h5s=""

Rscript make_fasta.R  $path_to_h5s"tumour1.h5" $path_to_h5s"tumour2.h5" $path_to_h5s"normal.h5" $path_to_h5s"ovary.h5"

mkdir IQ-TREE
mkdir IQ-TREE/filtered05

cp fasta/expr05.fasta IQ-TREE/filtered05

cd IQ-TREE/filtered05

iqtree2 -s expr05.fasta -st MORPH -m ORDERED+ASC
