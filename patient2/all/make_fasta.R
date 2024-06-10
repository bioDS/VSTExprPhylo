library("phyloRNA")

# enter correct paths to the phyloRNAanalysis scripts
import::from("../phyloRNAanalysis/src/expr.r", preprocess_expression)
import::from("../phyloRNAanalysis/src/utils.r", read_table, write_table, table2fasta)
import::from("../phyloRNAanalysis/src/filter.r", density_filtering)

import::from("../../scripts/functions.R", rename_categories)

# enter correct path
setwd("../VSTExprPhylo/patient2/all")

args = commandArgs(trailingOnly=TRUE)

if (length(args)<1) {
  stop("There has to be at least one argument: the path to one or more .h5 file.n", call.=FALSE)
}

hdi = c(0.6, 0.9)

expr_preprocessed = preprocess_expression(
    h5 = args,
    hdi = hdi,
    minGene=0,
    minUMI=0,
    outdir = file.path("preprocess", "expr"),
    prefix = "all"
)

filtdir = "filtered"
prefix = "expr"
mkdir(filtdir)

data = read_table(expr_preprocessed$discretized)
density=0.5
filter = density_filtering(data, density=density, empty="-", outdir=filtdir, prefix=prefix)

filtered=read_table(filter)

filtered[filtered=="-"] <- "0"

filtered <- rename_categories(filtered, 6)
write_table(filtered, filter)

expr_fasta = table2fasta(filter, outdir="fasta")
