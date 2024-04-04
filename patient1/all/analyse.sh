#P1 50% density analysis
# This script requiers
# 1. files tumour1.h5, tumour2.h5, normal1.h5, normal2.h5 to be in path_to_h5s folder
# 2. R
# 3. PhyloRNA package to be installed for R
# 4. PhyloRNAanalysis repo https://github.com/bioDS/phyloRNAanalysis + adjust paths to expr.r, utils.r and filter.r in make_fasta.R script
# 5. iqtree
# 6. adjust path_to_h5s here

#insert path to the folder with .h5 files:
path_to_h5s=""

Rscript make_fasta.R  $path_to_h5s"tumour1.h5" $path_to_h5s"tumour2.h5" $path_to_h5s"normal1.h5" $path_to_h5s"normal2.h5"

mkdir iqtree
mkdir iqtree/filtered05

cp fasta/expr05.fasta iqtree/filtered05

cd iqtree/filtered05

iqtree -s expr05.fasta -st MORPH -m ORDERED+ASC
