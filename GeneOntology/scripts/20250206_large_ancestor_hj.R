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

# 0. database load ----
gorea_enviromnet <- function(directory = NULL) {
  
  GOGMT <- gmtPathways(paste0(directory, "/GeneOntology/v2024.1.Hs/GOBP/c5.go.bp.v2024.1.Hs.symbols.gmt"))
  GOID_TERM <- read.table(paste0(directory, "/GeneOntology/v2024.1.Hs/GOBP/20250124_c5.go.bp.v2024.1.Hs.term.id_mapping_hj.txt"),
                          header = T,
                          sep = "\t") # this file is made from json file
  GO.ANCESTOR <<- readRDS(paste0(directory, "/GeneOntology/v2024.1.Hs/GOBP/20250124_c5.go.bp.v2024.1.Hs_ANCESTOR_df_hj.rds"))
  GO.LEVEL <<- readRDS(paste0(directory, "/GeneOntology/v2024.1.Hs/GOBP/20250124_c5.go.bp.v2024.1.Hs_LEVEL_df_hj.rds"))
  
  GO.ANCESTOR.LEVEL <- merge(GO.LEVEL, GO.ANCESTOR, by.x = "GOID", by.y = "ANCESTOR_GOID", all.y = T)
  colnames(GO.ANCESTOR.LEVEL) <- c("ANCESTOR_GOID", "ANCESTOR_LEVEL", "TARGET_GOID", "ANCESTOR_GOTERM")
  GO.ANCESTOR.LEVEL <<- GO.ANCESTOR.LEVEL %>% dplyr::filter(ANCESTOR_LEVEL != "all")
  GO.ANCESTOR.LEVEL$ANCESTOR_LEVEL %>% is.na() %>% which() # there are no NA data. 
  
  GO.ANCESTOR.SUBSET <<- GOID_TERM$GOID[which(GOID_TERM$GOID %in% GO.ANCESTOR$TARGET_GOID)] # this gobp is used for GSEA.
  GOID_TERM <<- GOID_TERM %>% dplyr::filter(GOID %in% GO.ANCESTOR.SUBSET)
  GOGMT <<- GOGMT[GOID_TERM$GOTERM] # utlize this object (7597 GOBPs) on fgsea 
  
  #if (!("hsGO" %in% ls(envir=.GlobalEnv))) {
  #  set.seed(1234); hsGO <<- godata('org.Hs.eg.db', ont="BP", computeIC = T) # make IC score using GOsemsim package 
  #}
}

gorea_enviromnet("/Users/hojin/Dropbox/project/enrichmate/20250206/")


# 1. large category ----

msigDB_GO_ancestor <- GO.ANCESTOR.LEVEL %>% dplyr::filter(TARGET_GOID %in% GOID_TERM$GOID)

level1_2_count_df <- msigDB_GO_ancestor %>%
  dplyr::filter(ANCESTOR_LEVEL %in% c(1,2)) %>%
  dplyr::group_by(ANCESTOR_GOID, ANCESTOR_LEVEL, ANCESTOR_GOTERM) %>%
  dplyr::summarise(n=n()) %>%
  dplyr::arrange(desc(n))

# histogram
level1_2_count_df%>% 
  ggplot(aes(x = n)) + 
  geom_histogram(fill = "red3", color = "black") +
  theme_bw() +
  theme(aspect.ratio = 1, panel.border = element_rect(linewidth = 1), axis.text = element_text(colour = "black"))

# 2. outlier test ----
library(outliers)

outlier_test_df <- data.frame()
level1_2_count_df %>% head()

res_tmp <- grubbs.test(level1_2_count_df$n, type = 10) # outlier test for a largest value.
pval_tmp <- res_tmp$p.value
go_val_tmp <- level1_2_count_df$n[1]
df_tmp <- data.frame(Ancestor_count = go_val_tmp, Pval = pval_tmp)
outlier_test_df <- rbind(outlier_test_df, df_tmp)
test_values <- level1_2_count_df$n

for (i in 1:113) { # Only 1~114 values can be tested, because the test allow at least 7 values and different value for each other.
  test_values <- test_values[-1]
  go_val_tmp <- test_values[1]
  res_tmp <- grubbs.test(test_values, type = 10) 
  pval_tmp <- res_tmp$p.value
  df_tmp <- data.frame(Ancestor_count = go_val_tmp, Pval = pval_tmp)
  outlier_test_df <- rbind(outlier_test_df, df_tmp)
}

# multiple testing correction 
outlier_test_df$Padj <- p.adjust(outlier_test_df$Pval, method = "bonferroni")
outlier_test_df <- outlier_test_df %>% dplyr::mutate(outlier = case_when(Padj < 0.05 ~ T,
                                                                         T ~ F))

outlier_test_df %>% head()

# histogram
level1_2_count_df <- level1_2_count_df %>% dplyr::mutate(Outlier = case_when(n %in% c(5005, 3335, 3092) ~ T,
                                                                             T ~ F))
pdf("/Users/hojin/Dropbox/project/enrichmate/20250206/20250206_ancestor_outlier_grubbs_test_hist_hj.pdf", width = 4.5, height = 4.5)
level1_2_count_df%>% 
  ggplot(aes(x = n)) + 
  xlab("Number of times considered an ancestor") +
  ylab("Count") +
  geom_histogram(aes(fill = Outlier), color = "black") +
  theme_bw() +
  theme(aspect.ratio = 1, panel.border = element_rect(linewidth = 1), axis.text = element_text(colour = "black")) +
  scale_fill_manual(values = c("darkgrey", "red3"))
dev.off()

# 3. large GO term ----
msigDB_GO_ancestor_large <- msigDB_GO_ancestor %>%
  dplyr::filter(!ANCESTOR_GOTERM %in% c("cellular process", "biological regulation", "regulation of biological process")) %>%
  dplyr::filter(ANCESTOR_LEVEL %in% c(1,2)) %>%
  dplyr::group_by(ANCESTOR_GOID, ANCESTOR_LEVEL, ANCESTOR_GOTERM) %>%
  dplyr::summarise(n=n()) %>% dplyr::arrange(desc(n))


saveRDS(msigDB_GO_ancestor_large, "/Users/hojin/Dropbox/project/enrichmate/20250206/GeneOntology/v2024.1.Hs/GOBP/20250206_msigDB_GO_ancestor_large_hj.rds")


### mouse ####
# 0. database load ----
gorea_enviromnet <- function(directory = NULL) {
  
  GOGMT <- gmtPathways(paste0(directory, "/GeneOntology/v2024.1.Mm/GOBP/m5.go.bp.v2024.1.Mm.symbols.gmt"))
  GOID_TERM <- read.table(paste0(directory, "/GeneOntology/v2024.1.Mm/GOBP/20250124_m5.go.bp.v2024.1.Mm.term.id_mapping_hj.txt"),
                          header = T,
                          sep = "\t") # this file is made from json file
  GO.ANCESTOR <<- readRDS(paste0(directory, "/GeneOntology/v2024.1.Mm/GOBP/20250124_m5.go.bp.v2024.1.Mm_ANCESTOR_df_hj.rds"))
  GO.LEVEL <<- readRDS(paste0(directory, "/GeneOntology/v2024.1.Mm/GOBP/20250124_m5.go.bp.v2024.1.Mm_LEVEL_df_hj.rds"))
  
  GO.ANCESTOR.LEVEL <- merge(GO.LEVEL, GO.ANCESTOR, by.x = "GOID", by.y = "ANCESTOR_GOID", all.y = T)
  colnames(GO.ANCESTOR.LEVEL) <- c("ANCESTOR_GOID", "ANCESTOR_LEVEL", "TARGET_GOID", "ANCESTOR_GOTERM")
  GO.ANCESTOR.LEVEL <<- GO.ANCESTOR.LEVEL %>% dplyr::filter(ANCESTOR_LEVEL != "all")
  GO.ANCESTOR.LEVEL$ANCESTOR_LEVEL %>% is.na() %>% which() # there are no NA data. 
  
  GO.ANCESTOR.SUBSET <<- GOID_TERM$GOID[which(GOID_TERM$GOID %in% GO.ANCESTOR$TARGET_GOID)] # this gobp is used for GSEA.
  GOID_TERM <<- GOID_TERM %>% dplyr::filter(GOID %in% GO.ANCESTOR.SUBSET)
  GOGMT <<- GOGMT[GOID_TERM$GOTERM] # utlize this object (7597 GOBPs) on fgsea 
  
  #if (!("hsGO" %in% ls(envir=.GlobalEnv))) {
  #  set.seed(1234); hsGO <<- godata('org.Hs.eg.db', ont="BP", computeIC = T) # make IC score using GOsemsim package 
  #}
}

gorea_enviromnet("/Users/hojin/Dropbox/project/enrichmate/20250206/")


# 1. large category ----

msigDB_GO_ancestor <- GO.ANCESTOR.LEVEL %>% dplyr::filter(TARGET_GOID %in% GOID_TERM$GOID)

level1_2_count_df <- msigDB_GO_ancestor %>%
  dplyr::filter(ANCESTOR_LEVEL %in% c(1,2)) %>%
  dplyr::group_by(ANCESTOR_GOID, ANCESTOR_LEVEL, ANCESTOR_GOTERM) %>%
  dplyr::summarise(n=n()) %>%
  dplyr::arrange(desc(n))

# histogram
level1_2_count_df%>% 
  ggplot(aes(x = n)) + 
  geom_histogram(fill = "red3", color = "black") +
  theme_bw() +
  theme(aspect.ratio = 1, panel.border = element_rect(linewidth = 1), axis.text = element_text(colour = "black"))

# 2. outlier test ----
library(outliers)

outlier_test_df <- data.frame()
level1_2_count_df %>% head()

res_tmp <- grubbs.test(level1_2_count_df$n, type = 10) # outlier test for a largest value.
pval_tmp <- res_tmp$p.value
go_val_tmp <- level1_2_count_df$n[1]
df_tmp <- data.frame(Ancestor_count = go_val_tmp, Pval = pval_tmp)
outlier_test_df <- rbind(outlier_test_df, df_tmp)
test_values <- level1_2_count_df$n

for (i in 1:90) { # Only 1~114 values can be tested, because the test allow at least 7 values and different value for each other.
  test_values <- test_values[-1]
  go_val_tmp <- test_values[1]
  res_tmp <- grubbs.test(test_values, type = 10) 
  pval_tmp <- res_tmp$p.value
  df_tmp <- data.frame(Ancestor_count = go_val_tmp, Pval = pval_tmp)
  outlier_test_df <- rbind(outlier_test_df, df_tmp)
}

# multiple testing correction 
outlier_test_df$Padj <- p.adjust(outlier_test_df$Pval, method = "bonferroni")
outlier_test_df <- outlier_test_df %>% dplyr::mutate(outlier = case_when(Padj < 0.05 ~ T,
                                                                         T ~ F))

outlier_test_df %>% head()

# histogram
level1_2_count_df <- level1_2_count_df %>% dplyr::mutate(Outlier = case_when(n %in% c(5099) ~ T,
                                                                             T ~ F))
pdf("/Users/hojin/Dropbox/project/enrichmate/20250206/20250206_ancestor_outlier_grubbs_test_hist_mouse_hj.pdf", width = 4.5, height = 4.5)
level1_2_count_df%>% 
  ggplot(aes(x = n)) + 
  xlab("Number of times considered an ancestor") +
  ylab("Count") +
  geom_histogram(aes(fill = Outlier), color = "black") +
  theme_bw() +
  theme(aspect.ratio = 1, panel.border = element_rect(linewidth = 1), axis.text = element_text(colour = "black")) +
  scale_fill_manual(values = c("darkgrey", "red3"))
dev.off()

# 3. large GO term ----
msigDB_GO_ancestor_large <- msigDB_GO_ancestor %>%
  dplyr::filter(!ANCESTOR_GOTERM %in% c("cellular process")) %>%
  dplyr::filter(ANCESTOR_LEVEL %in% c(1,2)) %>%
  dplyr::group_by(ANCESTOR_GOID, ANCESTOR_LEVEL, ANCESTOR_GOTERM) %>%
  dplyr::summarise(n=n()) %>% dplyr::arrange(desc(n))


saveRDS(msigDB_GO_ancestor_large, "/Users/hojin/Dropbox/project/enrichmate/20250206/GeneOntology/v2024.1.Mm/GOBP/20250206_msigDB_GO_ancestor_large_hj.rds")
