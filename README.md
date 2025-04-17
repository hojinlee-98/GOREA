# GOREA

## This tool is for summarizing GOBP and extracting meaningful biological context.

In GeneOntology directory, processed databases using GOBP originated from MSigDB (https://www.gsea-msigdb.org/gsea/msigdb) were contained.

To re-construct the database for specific version of MSigDB, utilize the following simple code. 

```bash
python ./scripts/20250124_go_term_id_mapping_hj.py [json file from MSigdb] [output file name]
```
- `20250124_go_term_id_mapping_hj.py` script is included in scripts directory under GeneOntology.

### Example
```bash
python ./scripts/20250124_go_term_id_mapping_hj.py ./v2024.1.Hs/GOBP/c5.go.bp.v2024.1.Hs.json ./v2024.1.Hs/GOBP/20250124_c5.go.bp.v2024.1.Hs.term.id_mapping_hj.txt
```

## Key algorithm for GOREA
<img width="830" alt="image" src="https://github.com/user-attachments/assets/17632e57-f23a-41e7-b8ee-aa371045a855" />

To conduct GOREA, a dataframe containing significant GOBP terms with proportion of overlapping genes relative to the total number of genes in each GOBP term or NES must be assigned as input data. First, using the input data, clustering step is conducted, wherein combined method that we devised to apply the positive aspects of binary cut and hierarchical clustering is utilized. To define representative terms, information about ancestor terms and levels of GOBP terms was used (see Supplementary Data). Specifically, the following steps were applied iteratively within each cluster: (1) Identify the common ancestor term at the highest level that covers a subset of the input GOBP terms. (2) For any remaining GOBP terms not explained by the representative term from step (1), the process is repeated on the remaining terms.


## The result of GOREA
<img width="913" alt="image" src="https://github.com/user-attachments/assets/8a36e81f-62d6-4c69-884d-c987c90faaf9" />

The clustering results are displayed as a heatmap, with the representative terms shown on the right. Clusters, as shown between the heatmap and the representative terms, are ordered in descending order based on the average proportion of overlapping genes or the absolute value of NES. To make more general observations about the resulting GOBP terms, we created a panel above the heatmap. First, we defined broad GOBP terms (see Supplementary Data). Broad GOBP terms including input GOBP terms as their child term are clustered. For each cluster, the broad GOBP terms that contain the highest number of significant GOBP terms as their child term are selected and displayed. On the right side of the broad GOBP terms’ panel, the percentage of GOBP terms that each broad GOBP terms encompasses as child term are indicated.


## Run GOREA
The following code is used to perform GOREA analysis.

```R

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
```

## Real example for GOREA analysis
<img width="1093" alt="image" src="https://github.com/user-attachments/assets/48ac6b32-1cc1-4c67-b243-6392b434f26f" />


