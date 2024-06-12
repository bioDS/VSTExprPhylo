library(ape)
library(stringr)
library(phytools)
library(ggplot2)
library(gridExtra) 
library(ggtree)
library(ips)
library(treeio)
library(Seurat)

process_tip_labels <- function(tree) {
  types <- str_split_i(tree$tip.label, "-", 2) 
  types <- str_replace_all(types, c("normal_b" = "Benign Normal", "normal1_b" = "Benign Normal", "normal2_b" = "Benign Normal",
                                    "tumour1_b" = "Benign Tumour", "tumour2_b" = "Benign Tumour", 
                                    "tumour1_m" = "Malignant Tumour", "tumour2_m" = "Malignant Tumour", 
                                    "ovary_m" = "Malignant Metastasis"))
  names(types) <- tree$tip.label
  types
}

count_cancer_clades <- function(tree, cancer=T) {
  types<-process_tip_labels(tree)
  n=tree$Nnode
  tumour_clades <- c()
  descendant_tumour_clades <- c()
  for (i in 1:(2*n+1)) {
    descendants <- getDescendants(tree, i)
    tips<-tree$tip.label[descendants[descendants <= n+1]]
    current_types<- unique(types[tips])
    if (cancer) {
      if ( !"Benign Normal" %in% current_types & !"Benign Tumour" %in% current_types) {
        tumour_clades <- c(tumour_clades,i)
        descendant_tumour_clades <-c(descendant_tumour_clades, descendants[descendants!=i])
      }  
    } else {
      if ( !"Benign Normal" %in% current_types & !"Benign Tumour" %in% current_types & !"Malignant Tumour" %in% current_types) {
        tumour_clades <- c(tumour_clades,i)
        descendant_tumour_clades <-c(descendant_tumour_clades, descendants[descendants!=i])
      }
    }
    
  }
  tumour_clades <- tumour_clades[!tumour_clades %in% descendant_tumour_clades]
  
  length(tumour_clades)
}

calculate_cc_mc_for_all <- function(tree_dists) {
  clades_df <- data.frame(
    Name = character(), 
    meanCC = numeric(), 
    meanMC = numeric()
  )
  
  for (i in 1:length(tree_dists)) {
    tree_dist_nobn <- tree_dists[[i]][round(length(tree_dists[[i]])/10):length(tree_dists[[i]])]  
    cc <- lapply(tree_dist_nobn, count_cancer_clades)
    mc <- lapply(tree_dist_nobn, function(x) count_cancer_clades(x, cancer=F))
    new_row <- data.frame(Name=names(tree_dists)[i], meanCC=mean(unlist(cc)), meanMC=mean(unlist(mc)))
    clades_df <- rbind(clades_df, new_row)
  }
  clades_df
}

external_to_total_edge_length <- function(trees) {
  ratio =c()
  for (tree in trees) {
    ntaxa=length(tree$tip.label)
    ratio <- c(ratio, sum(tree$edge.length[tree$edge[,2]<ntaxa+1])/sum(tree$edge.length))
  }
  ratio
}


phylogenetic_to_tc_plot <- function(tc, tree, legend_move_v=2, xlab=TRUE,  ylab=TRUE, analysis="", title="", plot_only=TRUE) {
  
  dist_tc <- dist(tc, method = "euclidean", diag = TRUE, upper = TRUE)
  
  distances= cophenetic.phylo(tree)
  distances <- distances[labels(dist_tc), labels(dist_tc)]
  
  dist = as.dist(distances)
  
  df <- data.frame(Phylo = as.vector(dist), TC = as.vector(dist_tc))
  
  y_label <- if (ylab) {
    paste(analysis, "Phylogenetic distance", sep="\n")  
  } else {
    "" 
  }
  
  x_label <- if(xlab) {
    "Total counts difference"
  } else {
    ""
  }
  
  p <- ggplot(df, aes(x = TC, y = Phylo)) +
    geom_point(size=0.5) +
    geom_smooth(method = "lm", se = FALSE, color = "red3") +
    labs(x = x_label, y = y_label, title = title)+       
    theme(title = element_text(face = "bold", size = 14), 
          axis.text = element_text(size = 12),    
          axis.title = element_text(size = 14))
  
  if (plot_only) {
    return(p)
  }
  
  correlation_coef <- cor(df$TC, df$Phylo, method = "pearson")
  
  p_with_label<- p+ geom_text(x = max(df$TC), y = min(df$Phylo), label = paste("R = ", round(correlation_coef, 2)), 
                              hjust = 1, vjust = 0, size = 6, color = "red3")

  print(p_with_label)
  
}


average_phy_dist_bwpairs <- function(tree, barcodes, P1=F, norm=T){
  
  dist <- cophenetic.phylo(tree)
  
  tumour1_m <-grep("tumour1_m", barcodes, value = TRUE)
  tumour2_m <-grep("tumour2_m", barcodes, value = TRUE)
  tumour1_b <-grep("tumour1_b", barcodes, value = TRUE)
  tumour2_b <-grep("tumour2_b", barcodes, value = TRUE)
  if (P1) {
    normal1 <- grep("normal1", barcodes, value = TRUE)
    normal2 <- grep("normal2", barcodes, value = TRUE)
  }

  
  if (length(tumour1_m)<2 | length(tumour2_m)<2 | length(tumour1_b)<2 | length(tumour2_b)<2) {
    print("WARNING: one of the tumour samples is empty or only have one cell")
  }
  if (P1) {
    if (length(normal1)<2 | length(normal2)<2) {
      print("WARNING: one of the normal samples is empty or only have one cell")
    }
  }
  
  dist_t1_m <- dist[tumour1_m, tumour1_m]
  dist_t2_m <- dist[tumour2_m, tumour2_m]
  dist_t_m <- dist[tumour1_m, tumour2_m]
  
  dist_t1_b <- dist[tumour1_b, tumour1_b]
  dist_t2_b <- dist[tumour2_b, tumour2_b]
  dist_t_b <- dist[tumour1_b, tumour2_b]
  
  between <- c(mean(dist_t_m), mean(dist_t_b))
  within <- c(mean(c(dist_t1_m[lower.tri(dist_t1_m)], dist_t2_m[lower.tri(dist_t2_m)])),mean(c(dist_t1_b[lower.tri(dist_t1_b)], dist_t2_b[lower.tri(dist_t2_b)])))
  
  if (P1) {
    dist_n1 <- dist[normal1, normal1]
    dist_n2 <- dist[normal2, normal2]
    dist_n <- dist[normal1, normal2]
    
    between <- c(between,mean(dist_n))
    within <- c(within, mean(c(dist_n1[lower.tri(dist_n1)], dist_n2[lower.tri(dist_n2)])))
    
  } else {
    between <- c(between,NA)
    within <- c(within,NA)
  }
  
  if (norm) {
    tree_height <- max(node.depth.edgelength(tree))
    between <- between/tree_height/2
    within <- within/tree_height/2
  }
  
  list(between, within)
  
}

average_phy_dist_bwpairs_ml <- function(tree, P1=F, norm=T){
  
  dist <- cophenetic.phylo(tree)
  barcodes <- tree$tip.label
  
  if (P1) {
    tumour1 <- grep("tumour1", barcodes, value = TRUE)
    tumour2 <- grep("tumour2", barcodes, value = TRUE)
    normal1 <- grep("normal1", barcodes, value = TRUE)
    normal2 <- grep("normal2", barcodes, value = TRUE)
  } else {
    tumour1 <- grep("-1", barcodes, value = TRUE)
    tumour2 <- grep("-2", barcodes, value = TRUE)
  }
  
  if (length(tumour1)<2 | length(tumour2)<2) {
    print("WARNING: one of the tumour samples is empty or only have one cell")
  }
  if (P1) {
    if (length(normal1)<2 | length(normal2)<2) {
      print("WARNING: one of the normal samples is empty or only have one cell")
    }
  }
  
  dist_t1 <- dist[tumour1, tumour1]
  dist_t2 <- dist[tumour2, tumour2]
  dist_t <- dist[tumour1, tumour2]
  
  between <- c(mean(dist_t))
  within <- c(mean(c(dist_t1[lower.tri(dist_t1)], dist_t2[lower.tri(dist_t2)])))
  
  if (P1) {
    dist_n1 <- dist[normal1, normal1]
    dist_n2 <- dist[normal2, normal2]
    dist_n <- dist[normal1, normal2]
    
    between <- c(between,mean(dist_n))
    within <- c(within, mean(c(dist_n1[lower.tri(dist_n1)], dist_n2[lower.tri(dist_n2)])))
    
  } else {
    between <- c(between,NA)
    within <- c(within,NA)
  }
  
  if (norm) {
    max_dist <- max(dist)
    between <- between/max_dist
    within <- within/max_dist
  }
  
  list(between, within)
  
}


plot_tree_tc_single_slice <- function(tree, tc, title="", tc_bar_width= 0.05,
                                      pdf_file=NULL, log=T, offset_second_bar=0.5, support=T) {
  
  sample_data <- data.frame(ID=names(tc), TotalCount=unname(tc)) 
  sample_data$Type <- rep("Malignant stroma", nrow(sample_data))
  sample_data$Type[grepl("-mg1", sample_data$ID, fixed=TRUE)] <- "Malignant region 1"
  sample_data$Type[grepl("-mg2", sample_data$ID, fixed=TRUE)] <- "Malignant region 2"
  
  p <- ggtree(tree, size=1.5) %<+% sample_data+
    aes(color=I(Type))+
    scale_color_manual(name = "Type", breaks = c("Malignant region 1", "Malignant region 2", "Malignant stroma"),
                       labels = c("Malignant region 1", "Malignant region 2", "Malignant stroma"),
                       values = c("red", "darkred", "orange"))+       
    theme(legend.title = element_text(face = "bold", size = 14),   
          legend.text=element_text(face = "bold", size = 12),  
          legend.position = "bottom",    
          legend.box = "horisontal",
          legend.margin = margin())
  
  if (support) {
    p <- p+geom_nodelab(aes(x=branch, label=round(posterior, 2)), hjust=-.5, size=4, color="black")
  }
  
  if (log) {
    
    log_tc_df <- as.data.frame(log(sample_data$TotalCount))
    row.names(log_tc_df) <- sample_data$ID
    
    p3 <- gheatmap(p, log_tc_df,
                   width = tc_bar_width,
                   colnames = FALSE)+ scale_fill_continuous(name = "Log total UMI",
                                                            low = "lightgreen", high = "darkblue",
                                                            na.value = "white")+
      theme(legend.position = "bottom",
            legend.title = element_text(size = 14),
            legend.text = element_text(size = 12),
            legend.box = "vertical", legend.margin = margin())
    
  } else {
    tc_df <- as.data.frame(sample_data$TotalCount)
    row.names(tc_df) <- sample_data$ID
    
    p3 <- gheatmap(p, tc_df, width = tc_bar_width, colnames = FALSE) + 
      scale_fill_continuous(name = "Total UMI", low = "yellow", high = "red", na.value = "white")+
      theme(legend.position = "bottom",
            legend.title = element_text(size = 14),
            legend.text = element_text(size = 12),
            legend.box = "vertical", legend.margin = margin())+
      ggtitle(title)
  }
  
  if (is.null(pdf_file)) {
    p3 
  } else {
    pdf(pdf_file)
    print(p3) 
    dev.off()
  }
}


plot_tree_tc <- function(tree, barcodes, barcodes_new_names, title="", tc_bar_width= 0.05, 
                         offset_second_bar=0.5, log=TRUE, pdf_file=NULL, custom_breaks=NULL) {
  tc <- total_counts[barcodes,]
  names(tc) <- barcodes_new_names
  
  sample_data <- data.frame(ID=barcodes_new_names, TotalCount=unname(tc)) 
  sample_data$Type <- rep("Benign Normal", nrow(sample_data))
  
  sample_data$Type[grepl("tumour1_m", sample_data$ID, fixed=TRUE)] <- "Malignant Tumour"
  sample_data$Type[grepl("tumour2_m", sample_data$ID, fixed=TRUE)] <- "Malignant Tumour"
  sample_data$Type[grepl("tumour1_b", sample_data$ID, fixed=TRUE)] <- "Benign Tumour"
  sample_data$Type[grepl("tumour2_b", sample_data$ID, fixed=TRUE)] <- "Benign Tumour"
  sample_data$Type[grepl("ovary_m", sample_data$ID, fixed=TRUE)] <- "Malignant Metastasis"    

  
  p <- ggtree(tree, size=1.5) %<+% sample_data+
    aes(color=I(Type))+
    scale_color_manual(name = "Type", breaks = c("Benign Normal", "Benign Tumour", "Malignant Tumour", "Malignant Metastasis"),
                       labels = c("Benign\nNormal", "Benign\nTumour", "Malignant\nTumour", "Malignant\nMetastasis"),
                       values = c("blue", "cyan", "red", "darkred"))+ 
    geom_nodelab(aes(x=branch, label=round(posterior, 2)), hjust=-.5, size=4, color="black")+       
    theme(legend.title = element_text(size = 14),   
          legend.text=element_text(size = 12),  
          legend.position = "bottom",    
          legend.box = "horizontal",
          legend.margin = margin())
  
  if (is.null(custom_breaks)) {
    custom_breaks <- c(ceiling(min(sample_data$TotalCount)/ 1000)*1000, 
                       round(min(sample_data$TotalCount)/1000 + (max(sample_data$TotalCount)-min(sample_data$TotalCount))/2000)*1000,
                       round(max(sample_data$TotalCount)/1000)*1000) 
  }
  
  
  
  if (log) {
    
    log_tc_df <- as.data.frame(log(sample_data$TotalCount))
    row.names(log_tc_df) <- sample_data$ID
    
    p3 <- gheatmap(p, log_tc_df,
                   width = tc_bar_width,
                   colnames = FALSE)+ scale_fill_continuous(name = "Log total UMI",
                                                            low = "lightgreen", high = "darkblue",
                                                            na.value = "white")+
      theme(legend.position = "bottom",
            legend.title = element_text(size = 14),
            legend.text = element_text(size = 12),
            legend.box = "vertical", legend.margin = margin())
  } else {
    tc_df <- as.data.frame(sample_data$TotalCount)
    row.names(tc_df) <- sample_data$ID
    
    p3 <- gheatmap(p, tc_df, width = tc_bar_width, colnames = FALSE) + 
      scale_fill_continuous(name = "Total UMI", low = "lightgreen", high = "darkblue", na.value = "white",
                            breaks = custom_breaks)+
      theme(legend.position = "bottom",
            legend.title = element_text(size = 14),
            legend.text = element_text(size = 12),
            legend.box = "vertical", legend.margin = margin())+
      ggtitle(title)
  }
  
  if (is.null(pdf_file)) {
    p3 
  } else {
    pdf(pdf_file)
    print(p3) 
    dev.off()
  }
  
}

# adjust path to the pateient 2 tumour 2 data and total_counts.tsv matrix which was calculated for aggregated data
path_to_tumour2_data = "/outs/filtered_feature_bc_matrix"
tc_file_path="../total_counts.tsv"

analyses <- c("zero", "zero_hvg", "zero_norm_hvg")
p2_analyses <- c(paste("BDStrictOrdinal", c(paste("selected", analyses, sep="_"), 
                       paste("t5q", analyses, sep="_"), 
                       paste("stc", analyses, sep="_")), sep="_"),
                 paste("BDStrictContinuous", c("selected", "t5q", "stc"), "hvg", sep="_"))

p2_names <-  c(paste("selected", analyses, sep="_"), 
               paste("t5q", analyses, sep="_"), 
               paste("stc", analyses, sep="_"), 
               paste(c("selected", "t5q", "stc"), "cont_hvg", sep="_"))

names(p2_analyses) <- paste("p2", p2_names, sep="_")

p2_analyses <- paste("../patient2/all/BEAST", p2_analyses, p2_names, sep="/")
p2_analyses <- paste(p2_analyses, "trees", sep=".")
names(p2_analyses) <- p2_names

p1_analyses <- paste("../patient1/all/BEAST/", c(rep("BDStrictOrdinal", 3), "BDStrictContinuous"), "_stc_", c(analyses, "hvg"), "/stc_", c(analyses, "cont_hvg"), ".trees", sep="")
          
names(p1_analyses) <- paste("p1_stc", c(analyses, "cont_hvg"), sep="_")

p1_trees <- lapply(p1_analyses, read.nexus)
names(p1_trees) <- names(p1_analyses)

p2_trees <- lapply(p2_analyses, read.nexus)
names(p2_trees) <- names(p2_analyses)

clades_p1 <- calculate_cc_mc_for_all(p1_trees)
clades_p2 <- calculate_cc_mc_for_all(p2_trees)
 
clades_df <- rbind(clades_p2, clades_p1)

write.csv(clades_df, "clades.csv", quote=F)
#clades_df <- read.csv("clades.csv", row.names=1)

All <- clades_df[! grepl("hvg", clades_df$Name), ]
HVG <- clades_df[grepl("zero_hvg", clades_df$Name), ]
NormHVG <- clades_df[grepl("zero_norm_hvg", clades_df$Name), ]
ContHVG <- clades_df[grepl("cont_hvg", clades_df$Name), ]

df_cc <- data.frame(All=All$meanCC, HVG=HVG$meanCC, NormHVG=NormHVG$meanCC, ContHVG=ContHVG$meanCC,
                    Allm=All$meanMC, HVGm=HVG$meanMC, NormHVGm=NormHVG$meanMC, ContHVGm=ContHVG$meanMC,
                    Analysis=c("P2 selected", "P2 highest", "P2 similar", "P1 similar"))

df_cc$Analysis <- factor(df_cc$Analysis, levels = c("P2 selected", "P2 highest", "P2 similar", "P1 similar"))

cc_plot4 <- ggplot(df_cc, aes(x = Analysis)) +
  geom_point(aes(y = HVG,
                 color = "HVG", shape = "HVG"), size = 4) +
  geom_point(aes(y = All,
                 color = "All", shape = "All"), size = 4) +
  geom_point(aes(y = NormHVG,
                 color = "NormHVG", shape = "NormHVG"), size = 4) +
  geom_point(aes(y = ContHVG,
                   color = "ContHVG", shape = "ContHVG"), size = 4) +
  labs(x = "Analysis", y = "Mean number of clades") +
  scale_color_manual(values = c("#00BFC4", "#F8766D", "#B0C4DE", "#FFD700"), breaks=c("HVG", "All", "NormHVG", "ContHVG"), labels=c("HVGs", "All", "HVGs norm", "HVGs cont"), name = "Genes") +
  scale_shape_manual(values = c(17, 19, 18, 16), breaks=c("HVG", "All", "NormHVG", "ContHVG"),labels=c("HVGs", "All", "HVGs norm", "HVGs cont"), name = "Genes") +
  theme(legend.text = element_text(size = 14), 
        legend.title = element_text(size = 16),
        axis.text = element_text(size = 12),    
        axis.title = element_text(size = 14),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(size = 16)) + ggtitle("a) Malignant clades")

mc_plot4 <- ggplot(df_cc, aes(x = Analysis)) +
  geom_point(aes(y = HVGm,
                 color = "HVGm", shape = "HVGm"), size = 4) +
  geom_point(aes(y = Allm,
                 color = "Allm", shape = "Allm"), size = 4) +
  geom_point(aes(y = NormHVGm,
                 color = "NormHVGm", shape = "NormHVGm"), size = 4) +
  geom_point(aes(y = ContHVGm,
                 color = "ContHVGm", shape = "ContHVGm"), size = 4) +
  labs(x = "Analysis", y = "Mean number of clades") +
  scale_color_manual(values = c("#00BFC4", "#F8766D", "#B0C4DE", "#FFD700"), breaks=c("HVGm", "Allm", "NormHVGm", "ContHVGm"), labels=c("HVGs", "All", "HVGs norm", "HVGs cont"), name = "Genes") +
  scale_shape_manual(values = c(17, 19, 18, 16), breaks=c("HVGm", "Allm", "NormHVGm", "ContHVGm"),labels=c("HVGs", "All", "HVGs norm", "HVGs cont"), name = "Genes") +
  theme(legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14),
        axis.text = element_text(size = 12),    
        axis.title = element_text(size = 14),
        axis.text.x = element_text(angle = 90, hjust = 0.5),
        plot.title = element_text(size = 16))+ ggtitle("b) Metastatic clades")


legend <- cowplot::get_legend(cc_plot4 + theme(legend.position = "right"))


## read in MCC trees
p2_tree_filenames <- c(paste("BDStrictOrdinal", c(paste("testset", analyses, sep="_"), 
                                            paste("selected", analyses, sep="_"), 
                                            paste("t5q", analyses, sep="_"), 
                                            paste("stc", analyses, sep="_")), sep="_"),
                 paste("BDStrictContinuous", c("testset", "selected", "t5q", "stc"), "hvg", sep="_"))

p2_names_all <-  c(paste("testset", analyses, sep="_"), 
               paste("selected", analyses, sep="_"), 
               paste("t5q", analyses, sep="_"), 
               paste("stc", analyses, sep="_"), 
               paste(c("testset", "selected", "t5q", "stc"), "cont_hvg", sep="_"))


prefixes<- c(rep("../patient2/tumour2/BEAST", 3), rep("../patient2/all/BEAST", 13))
prefixes[13] = "../patient2/tumour2/BEAST"

p2_tree_filenames <- paste(prefixes, p2_tree_filenames, paste("MCCkeep", p2_names_all, sep="_"), sep="/")
p2_tree_filenames<- paste(p2_tree_filenames, "tree", sep=".")
names(p2_tree_filenames) <- paste("p2", p2_names_all, sep="_")

p1_tree_filenames <- paste("../patient1/all/BEAST/",
                           c(rep("BDStrictOrdinal", 3), "BDStrictContinuous"), "_stc_", c(analyses, "hvg"), "/MCCkeep_stc_", c(analyses, "cont_hvg"), ".tree", sep="")

names(p1_tree_filenames) <- paste("p1_stc", c(analyses, "cont_hvg"), sep="_")


p1_mcctrees <- lapply(p1_tree_filenames, read.nexus)
names(p1_mcctrees) <- names(p1_tree_filenames)


p2_mcctrees <- lapply(p2_tree_filenames, read.nexus)
names(p2_mcctrees) <- names(p2_tree_filenames)


Moravec_na_hvg_tree <- read.nexus("../Trees/Moravec_58cells_hvg_na.tree")
Moravec_na_tree <- read.nexus("../Trees/Moravec_58cells_na.tree")



hvg_trees <- append(c(p2_mcctrees[grepl("zero_hvg", names(p2_mcctrees))], p1_mcctrees[grepl("zero_hvg", names(p1_mcctrees))]), list(Moravec_na_hvg_tree))
names(hvg_trees)[[length(hvg_trees)]] <- "Moravec_na_hvg_tree"

all_trees <- append(c(p2_mcctrees[!grepl("hvg", names(p2_mcctrees))], p1_mcctrees[!grepl("hvg", names(p1_mcctrees))]), list(Moravec_na_tree))
names(all_trees)[[length(all_trees)]] <- "Moravec_na_tree"

hgv_length_ratio = external_to_total_edge_length(hvg_trees)
length_ratio =external_to_total_edge_length(all_trees)
df_ratio <- data.frame(All=length_ratio, HVG=hgv_length_ratio, Analysis=c("P2 testset", "P2 selected", "P2 highest", "P2 similar", "P1 similar", "Moravec et al.\n2021"))

analysis_order <- c("P2 testset", "P2 selected", "P2 highest", "P2 similar", "P1 similar", "Moravec et al.\n2021")
df_ratio$Analysis <- factor(df_ratio$Analysis, levels = analysis_order)

length_ratio_plot <- ggplot(df_ratio, aes(x = Analysis)) +
  geom_point(aes(y = HVG, color = "HVG", shape = "HVG"), size = 4) +
  geom_point(aes(y = All, color = "All", shape = "All"), size = 4) +
  labs(x = "Analysis", y = "External to total branch length ratio") +
  labs(color = "Genes", shape = "Genes")+
  theme(legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14),
        axis.text = element_text(size = 12),    
        axis.title = element_text(size = 14),
        axis.text.x = element_text(angle = 90, hjust = 0.5),
        plot.title = element_text(size = 16))+ ggtitle("c) Branch lengths")


clade_branches_plot <-grid.arrange(
  arrangeGrob(
    cc_plot4 + theme(legend.position = "none"),
    legend,
    mc_plot4 + theme(legend.position = "none"),
    length_ratio_plot + theme(legend.position = "none"),
    ncol = 2,
    widths = c(1, 1),
    heights = c(1, 1.3)))

ggsave("clades_branches.pdf", plot = clade_branches_plot, width=6, height=7)

total_counts <- read.table(tc_file_path)

data_t2=Read10X(data.dir=path_to_tumour2_data)

data_t2 <- as.data.frame(data_t2)

t2_total_counts <- colSums(data_t2)

## Testset trees and tc

testset_barcodes_new_names <- p2_mcctrees[["p2_testset_zero_hvg"]]$tip.label
testset_barcodes <- str_replace(testset_barcodes_new_names, "-ms|-mg1|-mg2", "-2")
testset_tc <- t2_total_counts[str_replace(testset_barcodes, "-2", "-1")]
names(testset_tc) <- testset_barcodes_new_names

## P2 selected trees and tc

selected_barcodes <- scan("../patient2/all/selected_barcodes.txt", what=character())
selected_barcodes_new_names <- scan("../patient2/all/selected_new_barcodes.txt", what=character())
selected_tc <- total_counts[selected_barcodes,]
names(selected_tc) <- selected_barcodes_new_names

## P2 t5q trees and tc

t5q_barcodes <- scan("../patient2/all/t5q_barcodes.txt", what=character())
t5q_barcodes_new_names <- scan("../patient2/all/t5q_new_barcodes.txt", what=character())
t5q_tc <- total_counts[t5q_barcodes,]
names(t5q_tc) <- t5q_barcodes_new_names

## P2 stc trees and tc

stc_barcodes <- scan("../patient2/all/stc_barcodes.txt", what=character())
stc_barcodes_new_names <- scan("../patient2/all/stc_new_barcodes.txt", what=character())
stc_tc <- total_counts[stc_barcodes,]
names(stc_tc) <- stc_barcodes_new_names

## P1 stc
 
p1_stc_barcodes <-scan("../patient1/all/stc_barcodes.txt", what=character())
p1_stc_barcodes_new_names <-scan("../patient1/all/stc_new_barcodes.txt", what=character())
p1_stc_tc <- total_counts[p1_stc_barcodes,]
names(p1_stc_tc) <- p1_stc_barcodes_new_names


## phylogenetic vs total counts against distance plots:

plotR11 <- phylogenetic_to_tc_plot(testset_tc, p2_mcctrees[["p2_testset_zero"]], 3, FALSE, TRUE, analysis="All", title="P2 testset", plot_only=F)
plotR12 <- phylogenetic_to_tc_plot(testset_tc, p2_mcctrees[["p2_testset_zero_hvg"]], 3, TRUE, TRUE, analysis="HVGs", "", plot_only=F)

plotR21 <- phylogenetic_to_tc_plot(selected_tc, p2_mcctrees[["p2_selected_zero"]], 3, FALSE, FALSE, analysis="All", "P2 manually selected", plot_only=F)
plotR22 <- phylogenetic_to_tc_plot(selected_tc, p2_mcctrees[["p2_selected_zero_hvg"]], 3, TRUE, FALSE, analysis="HVGs", "", plot_only=F)

plotR31 <-phylogenetic_to_tc_plot(t5q_tc, p2_mcctrees[["p2_t5q_zero"]], 3, FALSE, FALSE, analysis="All", "P2 highest quality", plot_only=F)
plotR32 <-phylogenetic_to_tc_plot(t5q_tc, p2_mcctrees[["p2_t5q_zero_hvg"]], 3, TRUE, FALSE, analysis="HVGs", "", plot_only=F)

plotR41 <-phylogenetic_to_tc_plot(stc_tc, p2_mcctrees[["p2_stc_zero"]], 3, FALSE, FALSE, analysis="All", "P2 similar quality", plot_only=F)
plotR42 <-phylogenetic_to_tc_plot(stc_tc, p2_mcctrees[["p2_stc_zero_hvg"]], 3, TRUE, FALSE, analysis="HVGs", "", plot_only=F)

plotR51 <-phylogenetic_to_tc_plot(p1_stc_tc, p1_mcctrees[["p1_stc_zero"]], 3, FALSE, FALSE, analysis="All", "P1 similar quality", plot_only=F)
plotR52 <-phylogenetic_to_tc_plot(p1_stc_tc, p1_mcctrees[["p1_stc_zero_hvg"]], 3, TRUE, FALSE, analysis="HVGs", "", plot_only=F)

phylo_tc_arrangedR<- grid.arrange(plotR11, plotR21, plotR31, plotR41, plotR51, 
                                  plotR12, plotR22, plotR32, plotR42, plotR52, 
                                  ncol=5)

ggsave("phylo_tc_R.pdf", phylo_tc_arrangedR, width = 16, height = 8)



## find average distance between the pairs of serial sections

p1_stc_dist <- average_phy_dist_bwpairs(p1_mcctrees[["p1_stc_zero_hvg"]], p1_stc_barcodes_new_names, P1=T, norm=T)
stc_dist <- average_phy_dist_bwpairs(p2_mcctrees[["p2_stc_zero_hvg"]], stc_barcodes_new_names, P1=F, norm=T)
selected_dist <- average_phy_dist_bwpairs(p2_mcctrees[["p2_selected_zero_hvg"]], selected_barcodes_new_names, P1=F, norm=T)
t5q_dist <- average_phy_dist_bwpairs(p2_mcctrees[["p2_t5q_zero_hvg"]], t5q_barcodes_new_names, P1=F, norm=T)


p1_ml_tree <- read.tree("../Trees/p1_ml.treefile")
p2_ml_tree <- read.tree("../Trees/p2_ml.treefile")

p1_ml_dist <- average_phy_dist_bwpairs_ml(p1_ml_tree, P1=T, norm=T)
p2_ml_dist <- average_phy_dist_bwpairs_ml(p2_ml_tree, P1=F, norm=T)

dist_across_sec_df <- data.frame(BetweenSections=c(p1_stc_dist[[1]], stc_dist[[1]], selected_dist[[1]], t5q_dist[[1]]), 
                                 AveWithinSection=c(p1_stc_dist[[2]], stc_dist[[2]], selected_dist[[2]], t5q_dist[[2]]), 
                                 Tissue=c("Tumour_m", "Tumour_b", "Normal", "Tumour_m", "Tumour_b", "Normal", "Tumour_m", "Tumour_b", "Normal", "Tumour_m", "Tumour_b", "Normal"),
                                 Sample = c(rep("P1 similar",3), rep("P2 similar",3), rep("P2 selected",3), rep("P2 highest",3)),
                                 Patient=c(rep("P1",3), rep("P2",9)))

dist_across_sec_df_ml <- data.frame(BetweenSections=c(p1_ml_dist[[1]], p2_ml_dist[[1]]), 
                                    AveWithinSection=c(p1_ml_dist[[2]], p2_ml_dist[[2]]), 
                                    Tissue=c("Tumour", "Normal", "Tumour","Normal"),
                                    Sample = c(rep("P1 50% density",2), rep("P2 50% density",2)),
                                    Patient=c(rep("P1",2), rep("P2",2)))

dist_across_sec_df <- rbind(dist_across_sec_df, dist_across_sec_df_ml)

## exclude tumour_b P1 analysis

dist_across_sec_df <- dist_across_sec_df[dist_across_sec_df$Sample!="P1 similar" | dist_across_sec_df$Tissue != "Tumour_b",]

plot_sections <- ggplot(dist_across_sec_df, aes(y = BetweenSections, x = AveWithinSection, color = Sample, shape=Tissue))+geom_point(size=3)+ 
  geom_abline(intercept = 0, slope = 1, color = "blue") +
  labs(y = "Average pairwise distances between two sections",x = "Average pairwise distances within same sections",  title = "Phylogenetic distance within and between sections")


ggsave("sections_separation.pdf", plot_sections, width=6, height=5)


## plot trees:

p1_beast_trees <- lapply(p1_tree_filenames, read.beast)
names(p1_beast_trees) <- names(p1_tree_filenames)

p2_beast_trees <- lapply(p2_tree_filenames, read.beast)
names(p2_beast_trees) <- names(p2_tree_filenames)

if (!file.exists("tree_pics")) {
  dir.create("tree_pics")
}

plot_tree_tc(p2_beast_trees[["p2_selected_zero"]], selected_barcodes, selected_barcodes_new_names, log=T,
             pdf_file="tree_pics/selected_zero_tree_support.pdf")
plot_tree_tc(p2_beast_trees[["p2_selected_zero_hvg"]], selected_barcodes, selected_barcodes_new_names, log=T,
             pdf_file="tree_pics/selected_zero_hvg_tree_support.pdf")
plot_tree_tc(p2_beast_trees[["p2_selected_cont_hvg"]], selected_barcodes, selected_barcodes_new_names, 
             pdf_file="tree_pics/selected_cont_hvg_tree_support.pdf")
plot_tree_tc(p2_beast_trees[["p2_selected_zero_norm_hvg"]], selected_barcodes, selected_barcodes_new_names, log=T,
             offset_second_bar=40, tc_bar_width= 0.03,
             pdf_file="tree_pics/selected_zero_norm_hvg_tree_support.pdf")

plot_tree_tc(p2_beast_trees[["p2_stc_cont_hvg"]], stc_barcodes, stc_barcodes_new_names, log=F, custom_breaks=c(300, 350, 400),
             pdf_file="tree_pics/stc_cont_hvg_tree_support.pdf")
plot_tree_tc(p2_beast_trees[["p2_stc_zero_hvg"]], stc_barcodes, stc_barcodes_new_names, log=F, custom_breaks=c(300, 350, 400),
             pdf_file="tree_pics/stc_zero_hvg_tree_support.pdf")
plot_tree_tc(p2_beast_trees[["p2_stc_zero"]], stc_barcodes, stc_barcodes_new_names, log=F, custom_breaks=c(300,350,400),
             pdf_file="tree_pics/stc_zero_tree_support.pdf")
plot_tree_tc(p2_beast_trees[["p2_stc_zero_norm_hvg"]], stc_barcodes, stc_barcodes_new_names, custom_breaks=c(300,350,400),
             log=F,
             pdf_file="tree_pics/stc_zero_norm_hvg_tree_support.pdf")

plot_tree_tc(p2_beast_trees[["p2_t5q_zero"]], t5q_barcodes, t5q_barcodes_new_names, custom_breaks=c(500,4000,8000), log=F,
             pdf_file="tree_pics/t5q_zero_tree_support.pdf")
plot_tree_tc(p2_beast_trees[["p2_t5q_zero_hvg"]], t5q_barcodes, t5q_barcodes_new_names, custom_breaks=c(500,4000,8000), log=F,
             offset_second_bar=0.8, pdf_file="tree_pics/t5q_zero_hvg_tree_support.pdf")
plot_tree_tc(p2_beast_trees[["p2_t5q_cont_hvg"]], t5q_barcodes, t5q_barcodes_new_names, 
             offset_second_bar=10, pdf_file="tree_pics/t5q_cont_hvg_tree_support.pdf")

plot_tree_tc(p2_beast_trees[["p2_t5q_zero_norm_hvg"]], t5q_barcodes, t5q_barcodes_new_names, custom_breaks=c(500,4000,8000),
             log=F,
             pdf_file="tree_pics/t5q_zero_norm_hvg_tree_support.pdf")

plot_tree_tc(p1_beast_trees[["p1_stc_zero_hvg"]], p1_stc_barcodes, p1_stc_barcodes_new_names, custom_breaks=c(1600,1800,2000),
             pdf_file="tree_pics/p1_stc_zero_hvg_tree_support.pdf", log=F)
plot_tree_tc(p1_beast_trees[["p1_stc_zero"]], p1_stc_barcodes, p1_stc_barcodes_new_names, custom_breaks=c(1600,1800,2000),
             pdf_file="tree_pics/p1_stc_zero_tree_support.pdf", log=F)
plot_tree_tc(p1_beast_trees[["p1_stc_cont_hvg"]], p1_stc_barcodes, p1_stc_barcodes_new_names, custom_breaks=c(1600,1800,2000),
             log=F, pdf="tree_pics/p1_stc_cont_hvg_tree_support.pdf")
plot_tree_tc(p1_beast_trees[["p1_stc_zero_norm_hvg"]], p1_stc_barcodes, p1_stc_barcodes_new_names, custom_breaks=c(1600,1800,2000),
             log=F, pdf="tree_pics/p1_stc_zero_norm_hvg_tree_support.pdf")

plot_tree_tc_single_slice(p2_beast_trees[["p2_testset_zero_hvg"]], testset_tc, offset_second_bar=0.3, 
                          pdf="tree_pics/testset_zero_hvg_tree_support.pdf")
plot_tree_tc_single_slice(p2_beast_trees[["p2_testset_zero"]], testset_tc, offset_second_bar=0.3, 
                          pdf="tree_pics/testset_zero_tree_support.pdf")
plot_tree_tc_single_slice(p2_beast_trees[["p2_testset_cont_hvg"]], testset_tc, offset_second_bar=0.3, 
                          pdf="tree_pics/testset_cont_hvg_tree_support.pdf")
plot_tree_tc_single_slice(p2_beast_trees[["p2_testset_zero_norm_hvg"]], testset_tc, offset_second_bar=0.3, 
                          pdf="tree_pics/testset_zero_norm_hvg_tree_support.pdf")

