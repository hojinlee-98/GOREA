
### human tutorial ####

library(dplyr)
library(plyr)
library(fgsea)
library(tibble)
library(ggplot2)
library(GO.db)
library(GOSemSim)
library(WriteXLS)
library(colorRamp2)
library(simplifyEnrichment)
library(ComplexHeatmap)
library(org.Hs.eg.db) # human 

setwd("/Users/hojin/Dropbox/project/GOREA/20250401/human/")
source("/Users/hojin/Dropbox/project/GOREA/20250401/human/20250401_gorea_function_human_hj.R")

# 1. Setting environment for analysis ----

localdir <- "/Users/hojin/Dropbox/project/GOREA/20250401/" # this directory has to be parents directory of GeneOntology directory
gorea_enviromnet(localdir)

# 2. example ---

## 2.1 make test data ----
test <- sample(GOID_TERM$GOID, 500, replace = F) 
input_df <- data.frame(GOID = test)
input_df$NES <- sample(seq(0.1, 4, 0.01), replace = T, size = nrow(input_df))

head(input_df) 

## 2.2 outlier plot (additional step) ----
# before starting clustering steps, you can check a plot for the number of small clusters depending on cutoff (the cutoff is the value that you can assign in the gorea function as a parameter)
w <- gorea_sim_mat(input = input_df, godata_GO = godata_GO)
gorea_outlier_plot(w = w)

## 2.3 GOREA main function ----
# input_df; input dataframe from tools such as fgsea and enrichGO.
# k_val; when increasing this value, the number of cluster is increased.
# cutoff; if you want to remove large amount of small clusters, decrease this cutoff. but, according to simplifyenrichment, 0.85 is a default value. 
# top_ancestor_annotation; for general description, top panel for large GOBP terms is implemented.
# top_ancestor; clustered large GO terms will be split according to this number. and the high ranked large GO terms based on the number of input GOBP terms as child terms for each cluster can be illustrated in the top panel. 
res <- gorea(input = input_df,
             k_val = 10, # considering your total number of GOBP terms.
             godata_GO = godata_GO,
             cutoff = 0.85, # default (you can change this value, according to the results from gorea_outlier_plot function)
             outlier_detect = T,
             min_cluster = 3,
             representative_term_level_cutoff = 1, GO_explain = 3,
             score = "NES", # "NES" or "Overlap_freq"
             filename1 = "testfile1_hj.xlsx",
             filename2 = "testfile2_hj.xlsx",
             heatmap_filename = "testplot_hj.png",
             plot = T,
             heatmap_width = 40, heatmap_height = 30,
             ancestor_annotation = T,
             right_annotation_font_size = 10,
             cluster_font_size = 4,
             top_ancestor_annotation = T,
             top_ancestor = 3,
             color = c("gold"))

### mouse tutorial ####

library(dplyr)
library(plyr)
library(fgsea)
library(tibble)
library(ggplot2)
library(GO.db)
library(GOSemSim)
library(WriteXLS)
library(colorRamp2)
library(simplifyEnrichment)
library(ComplexHeatmap)
library(org.Mm.eg.db)

setwd("/Users/hojin/Dropbox/project/GOREA/20250401/mouse/")
source("/Users/hojin/Dropbox/project/GOREA/20250401/mouse/20250401_gorea_function_mouse_hj.R")

# 1. Setting environment for analysis ----

localdir <- "/Users/hojin/Dropbox/project/GOREA/20250401/" # this directory has to be parents directory of GeneOntology directory
gorea_enviromnet(localdir)

# 2. example ---

## 2.1 make test data ----
test <- sample(GOID_TERM$GOID, 500, replace = F) 
input_df <- data.frame(GOID = test)
input_df$NES <- sample(seq(0.1, 4, 0.01), replace = T, size = nrow(input_df))

head(input_df) 

## 2.2 outlier plot (additional step) ----
# before starting clustering steps, you can check a plot for the number of small clusters depending on cutoff (the cutoff is the value that you can assign in the gorea function as a parameter)
w <- gorea_sim_mat(input = input_df, godata_GO = godata_GO)
gorea_outlier_plot(w = w)


## 2.3 GOREA main function ----
# input_df; input dataframe from tools such as fgsea and enrichGO.
# k_val; when increasing this value, the number of cluster is increased.
# cutoff; if you want to remove large amount of small clusters, decrease this cutoff. but, according to simplifyenrichment, 0.85 is a default value. 
# top_ancestor_annotation; for general description, top panel for large GOBP terms is implemented.
# top_ancestor; clustered large GO terms will be split according to this number. and the high ranked large GO terms based on the number of input GOBP terms as child terms for each cluster can be illustrated in the top panel. 
res <- gorea(input = input_df,
             k_val = 10, # considering your total number of GOBP terms.
             godata_GO = godata_GO,
             cutoff = 0.85, # default (you can change this value, according to the results from gorea_outlier_plot function)
             outlier_detect = T,
             min_cluster = 3,
             representative_term_level_cutoff = 1, GO_explain = 3,
             score = "NES", # "NES" or "Overlap_freq"
             filename1 = "testfile1_hj.xlsx",
             filename2 = "testfile2_hj.xlsx",
             heatmap_filename = "testplot_hj.png",
             plot = T,
             heatmap_width = 40, heatmap_height = 30,
             ancestor_annotation = T,
             right_annotation_font_size = 10,
             cluster_font_size = 4,
             top_ancestor_annotation = T,
             top_ancestor = 3,
             color = c("gold"))

### mouse tutorial v7.5.1 ####

library(dplyr)
library(plyr)
library(fgsea)
library(tibble)
library(ggplot2)
library(GO.db)
library(GOSemSim)
library(WriteXLS)
library(colorRamp2)
library(simplifyEnrichment)
library(ComplexHeatmap)
library(org.Mm.eg.db)

setwd("/Users/hojin/Dropbox/project/GOREA/20250401/mouse_7.5.1/")
source("/Users/hojin/Dropbox/project/GOREA/20250401/mouse_7.5.1/20250401_gorea_function_mouse_v7.5.1_hj.R")

# 1. Setting environment for analysis ----

localdir <- "/Users/hojin/Dropbox/project/GOREA/20250401/" # this directory has to be parents directory of GeneOntology directory
gorea_enviromnet(localdir)

# 2. example ---

## 2.1 make test data ----
test <- sample(GOID_TERM$GOID, 500, replace = F) 
input_df <- data.frame(GOID = test)
input_df$NES <- sample(seq(0.1, 4, 0.01), replace = T, size = nrow(input_df))

head(input_df) 

## 2.2 outlier plot (additional step) ----
# before starting clustering steps, you can check a plot for the number of small clusters depending on cutoff (the cutoff is the value that you can assign in the gorea function as a parameter)
w <- gorea_sim_mat(input = input_df, godata_GO = godata_GO)
gorea_outlier_plot(w = w)


## 2.3 GOREA main function ----
# input_df; input dataframe from tools such as fgsea and enrichGO.
# k_val; when increasing this value, the number of cluster is increased.
# cutoff; if you want to remove large amount of small clusters, decrease this cutoff. but, according to simplifyenrichment, 0.85 is a default value. 
# top_ancestor_annotation; for general description, top panel for large GOBP terms is implemented.
# top_ancestor; clustered large GO terms will be split according to this number. and the high ranked large GO terms based on the number of input GOBP terms as child terms for each cluster can be illustrated in the top panel. 
res <- gorea(input = input_df,
             k_val = 10, # considering your total number of GOBP terms.
             godata_GO = godata_GO,
             cutoff = 0.85, # default (you can change this value, according to the results from gorea_outlier_plot function)
             outlier_detect = T,
             min_cluster = 3,
             representative_term_level_cutoff = 1, GO_explain = 3,
             score = "NES", # "NES" or "Overlap_freq"
             filename1 = "testfile1_hj.xlsx",
             filename2 = "testfile2_hj.xlsx",
             heatmap_filename = "testplot_hj.png",
             plot = T,
             heatmap_width = 40, heatmap_height = 30,
             ancestor_annotation = T,
             right_annotation_font_size = 10,
             cluster_font_size = 4,
             top_ancestor_annotation = T,
             top_ancestor = 3,
             color = c("gold"))
