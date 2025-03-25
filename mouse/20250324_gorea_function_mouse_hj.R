
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
  
  GO.LARGE <<- readRDS(paste0(directory, "/GeneOntology/v2024.1.Mm/GOBP/20250206_msigDB_GO_ancestor_large_hj.rds"))
  
  if (!("godata_GO" %in% ls(envir=.GlobalEnv))) {
    set.seed(1234); godata_GO <<- godata('org.Mm.eg.db', ont="BP", computeIC = T) # make IC score using GOsemsim package 
  }
}

gorea_sim_mat <- function(input = NULL,
                          sim_method = "Rel",
                          godata_GO = NULL) {
  
  goids <- input$GOID
  set.seed(1234); w <- GOSemSim::mgoSim(goids, goids, semData = godata_GO, measure = sim_method, combine = NULL) # w is similarity matrix. 
  na_goid <- names(which(is.na(w[1,]))) # select NA goid.
  goids <- goids[which(!goids %in% na_goid)]
  set.seed(1234); w <- GOSemSim::mgoSim(goids, goids, semData = godata_GO, measure = sim_method, combine = NULL) # w is similarity
  
  return(w)
}

gorea_outlier_plot <- function(w = NULL,
                               min_cluster = 3) {
  
  cutoff <- seq(0.6, 0.9, by = 0.01)
  cutoff = cutoff[cutoff >= 0.5 & cutoff <= 1]
  s1 = s2 = s3 = NULL
  
  for (i in seq_along(cutoff)) {
    set.seed(1234)
    cl = binary_cut(w, cutoff = cutoff[i])
    s1[i] = difference_score(w, cl)
    tb = table(cl)
    s2[i] = length(tb)
    s3[i] = sum(tb < min_cluster)
    
  }
  
  df1 <- data.frame(cutoff, s2)
  colnames(df1) <- c("method", "value")
  df2 <- data.frame(cutoff, s3)
  colnames(df2) <- c("method", "value")
  df1$type <- "All sizes"
  df2$type <- paste("Size", "<", min_cluster)
  df <- rbind(df1, df2)
  
  p <- df %>% ggplot(aes(x = method, y = value, color = type)) +
    geom_point() +
    theme_classic() +
    xlab("Cutoff") +
    ylab("# of clusters") +
    ggtitle("Outlier detection") +
    scale_color_discrete("Type of clusters") +
    theme(aspect.ratio = 0.5,
          plot.title = element_text(hjust = 0.5, size = 15)) + 
    scale_color_manual(values = c("darkgrey", "red3"))
  
  return(p) # select outlier cutoff using this plot.
  
}


gorea_clustering <- function(cutoff = 0.84,
                             w = NULL,
                             k_val = NULL,
                             input = input_df,
                             outlier_detect = T,
                             min_cluster = 3) {
  
  if (outlier_detect == T) {
    
    # binary clustering
    set.seed(1234); cl <- binary_cut(mat = w, cutoff = cutoff, partial = F)
    bc_res <- data.frame(GOID = colnames(w), cluster = cl) 
    outlier_clusters <- bc_res %>% dplyr::group_by(cluster) %>% dplyr::summarise(n = n()) %>% dplyr::filter(n < min_cluster) %>% dplyr::pull(cluster)
    
    set.seed(1234); pseudoname_outliers <- sample(10000:20000, size = length(outlier_clusters), replace = F)
    pseudoname_outliers_df <- data.frame(pseudoname = pseudoname_outliers, cluster = outlier_clusters) 
    pseudoname_outliers_df <- merge(pseudoname_outliers_df, bc_res, by = "cluster", all.x = T)
    
  } else {
    
    pseudoname_outliers_df <- NULL
    
  }
  
  # hclust ward.D
  set.seed(1234); dmat <- stats::as.dist(1 - w) # distance matrix.
  set.seed(1234); hc <- stats::hclust(dmat, method = "ward.D") # h.clust using ward.D.
  
  hc_res_list <- list()
  set.seed(1234); hclust_wardD <- stats::cutree(hc, k = k_val)
  hc_res_df <- data.frame(GOID = names(hclust_wardD),
                          cluster = hclust_wardD,
                          row.names = NULL)
  
  hc_res_df <- hc_res_df[hc$order, ] # ordering
  rownames(hc_res_df) <- NULL
  hc_res_df$key <- rownames(hc_res_df)
  hc_res_df <- hc_res_df %>% dplyr::arrange(as.numeric(key))
  
  hc_cluster_df <- merge(hc_res_df, input, by = "GOID", all.x = T)
  hc_cluster_df <- hc_cluster_df %>% dplyr::arrange(as.numeric(key))
  hc_cluster_df <- hc_cluster_df %>% dplyr::select(-c("key"))
  
  final_res_df <- merge(hc_res_df, GOID_TERM, by = "GOID", all.x = T)
  final_res_df <- merge(final_res_df, GO.LEVEL, by = "GOID", all.x = T) # get LEVEL information
  final_res_df <- final_res_df %>% dplyr::arrange(cluster) # arrange cluster
  
  my_list <- list("pseudoname_outliers_df" = pseudoname_outliers_df,
                  "hc_cluster_df" = hc_cluster_df,
                  "final_res_df" = final_res_df)
  
  return(my_list)
  
}

gorea_representative_terms <- function(gorea_clustering_res = NULL,
                                       k_val = NULL,
                                       outlier_detect = T,
                                       representative_term_level_cutoff = 1,
                                       GO_explain = 2
) {
  
  pseudoname_outliers_df <- gorea_clustering_res$pseudoname_outliers_df
  final_res_df <- gorea_clustering_res$final_res_df
  
  GO.ANCESTOR.LEVEL_subset <- GO.ANCESTOR.LEVEL %>% dplyr::filter(ANCESTOR_LEVEL >= representative_term_level_cutoff)
  
  if (is.null(pseudoname_outliers_df)) {
    
    outlier_detect = F
    
  } else {
    
    if (nrow(pseudoname_outliers_df) == 0) {
      
      outlier_detect = F
      
    }
  }
  
  if (outlier_detect == T) {
    final_res_df <- final_res_df %>% dplyr::filter(!(GOID %in% pseudoname_outliers_df$GOID))
  }
  
  df <- data.frame()
  ans_df <- data.frame()
  first_rep_c_goid_df <- data.frame()
  
  for (c_tmp in sort(unique(final_res_df$cluster))) {
    final_res_df_c_tmp <- final_res_df %>% dplyr::filter(cluster == c_tmp)
    ancestor_df_tmp <- GO.ANCESTOR.LEVEL_subset %>% dplyr::filter(TARGET_GOID %in% final_res_df_c_tmp$GOID)
    ancestor_df_tmp$cluster <- c_tmp
    ans_df <- rbind(ans_df, ancestor_df_tmp)
    
    seed_c_tmp <- final_res_df_c_tmp
    remained_goid <- final_res_df_c_tmp$GOID
    
    
    while (TRUE) {
      
      explained_go <- floor(nrow(seed_c_tmp)/GO_explain)
      ancestor_sub_df_tmp <- GO.ANCESTOR.LEVEL_subset %>% dplyr::filter(TARGET_GOID %in% seed_c_tmp$GOID)
      
      ancestor_sub_df_filtered <- ancestor_sub_df_tmp %>%
        dplyr::group_by(ANCESTOR_GOID, ANCESTOR_LEVEL) %>%
        dplyr::summarise(n=n()) %>%
        dplyr::arrange(desc(n)) %>%
        dplyr::filter(n > explained_go)
      
      if (nrow(ancestor_sub_df_filtered) != 0) {
        
        max_ancestor_tmp <- ancestor_sub_df_filtered %>%
          dplyr::ungroup(ANCESTOR_GOID, ANCESTOR_LEVEL) %>%
          dplyr::filter(ANCESTOR_LEVEL == max(ANCESTOR_LEVEL)) %>%
          dplyr::filter(n == max(n))
        
        first_rep_c_goid <- ancestor_sub_df_tmp %>%
          dplyr::filter(ANCESTOR_GOID %in% max_ancestor_tmp$ANCESTOR_GOID) %>%
          dplyr::pull(TARGET_GOID) %>% unique()
        
        # get the common ancestor terms and child term 
        first_rep_c_goid_df_tmp <- ancestor_sub_df_tmp %>% dplyr::filter(TARGET_GOID %in% first_rep_c_goid &
                                                                           ANCESTOR_GOID %in% max_ancestor_tmp$ANCESTOR_GOID)
        first_rep_c_goid_df_tmp$cluster <- c_tmp
        first_rep_c_goid_df <- rbind(first_rep_c_goid_df, first_rep_c_goid_df_tmp) # ancestor and child dataframe
        
        remained_goid <- seed_c_tmp %>%
          dplyr::filter(!GOID %in% first_rep_c_goid) %>%
          dplyr::pull(GOID)
        
        seed_c_tmp <- seed_c_tmp %>% dplyr::filter(GOID %in% remained_goid)
        
        df_tmp <- data.frame(cluster = rep(c_tmp, nrow(max_ancestor_tmp)),
                             max_level = rep(unique(max_ancestor_tmp$ANCESTOR_LEVEL), nrow(max_ancestor_tmp)),
                             n = rep(unique(max_ancestor_tmp$n), nrow(max_ancestor_tmp)),
                             total_gobp = rep(nrow(final_res_df_c_tmp), nrow(max_ancestor_tmp)),
                             n_total_gobp = rep(paste(c(unique(max_ancestor_tmp$n), nrow(final_res_df_c_tmp)), collapse = " out of "), nrow(max_ancestor_tmp)),
                             ANCESTOR_GOID = max_ancestor_tmp$ANCESTOR_GOID)
        
        df <- rbind(df, df_tmp)
        
        if (length(remained_goid) == 1) {
          df_tmp <- data.frame(cluster = rep(c_tmp, nrow(seed_c_tmp)),
                               max_level = rep(".", nrow(seed_c_tmp)),
                               n = rep(1, nrow(seed_c_tmp)),
                               total_gobp = rep(nrow(final_res_df_c_tmp), nrow(seed_c_tmp)),
                               n_total_gobp = rep(paste(c(1, nrow(final_res_df_c_tmp)), collapse = " out of "), nrow(seed_c_tmp)),
                               ANCESTOR_GOID = seed_c_tmp$GOID)
          df <- rbind(df, df_tmp)
          
          break
        } else if (length(remained_goid) == 0) {
          break
        }
        
      } else {
        
        df_tmp <- data.frame(cluster = rep(c_tmp, nrow(seed_c_tmp)),
                             max_level = rep(".", nrow(seed_c_tmp)),
                             n = rep(1, nrow(seed_c_tmp)),
                             total_gobp = rep(nrow(final_res_df_c_tmp), nrow(seed_c_tmp)),
                             n_total_gobp = rep(paste(c(1, nrow(final_res_df_c_tmp)), collapse = " out of "), nrow(seed_c_tmp)),
                             ANCESTOR_GOID = seed_c_tmp$GOID)
        
        df <- rbind(df, df_tmp)
        
        break
        
      }
      
    }
    
  }
  
  rep_df <- df
  
  return(rep_df)
  
}

gorea_getorder <- function(w = NULL,
                           gorea_clustering_res = NULL,
                           input = NULL,
                           score = NULL,
                           outlier_detect = T) {
  
  pseudoname_outliers_df <- gorea_clustering_res$pseudoname_outliers_df
  hc_cluster_df <- gorea_clustering_res$hc_cluster_df
  
  if (is.null(pseudoname_outliers_df)) {
    
    outlier_detect = F
    
  } else {
    
    if (nrow(pseudoname_outliers_df) == 0) {
      
      outlier_detect = F
      
    }
  }
  
  # if you do not want to assign outliers, change outlier_detect
  if (outlier_detect == T) {
    
    pseudoname_outliers_df <- merge(pseudoname_outliers_df, input, by = "GOID", all.x = T)
    
    pseudoname_outliers_df <- pseudoname_outliers_df %>% dplyr::select(-c("cluster")) %>% dplyr::rename("cluster" = pseudoname) 
    
    # subset w matrix to remove outliers 
    outliers_rm_goids <- rownames(w)[!rownames(w) %in% pseudoname_outliers_df$GOID]
    w <- w[outliers_rm_goids, outliers_rm_goids]
    
    # final order
    hc_cluster_df <- hc_cluster_df %>% dplyr::filter(!(GOID %in% pseudoname_outliers_df$GOID))
    hc_cluster_df$order <- as.numeric(rownames(hc_cluster_df))
    bc_cluster_df <- pseudoname_outliers_df %>% dplyr::arrange(cluster)
    bc_cluster_tb <- table(pseudoname_outliers_df$cluster)
    bc_cluster_numb <- length(bc_cluster_tb)
    bc_cluster_orig_name <- as.numeric(names(bc_cluster_tb))
    max_cl <- max(hc_cluster_df$cluster)
    cl_mold <- data.frame(new_cluster = seq(max_cl+1, max_cl+bc_cluster_numb),
                          cluster = bc_cluster_orig_name)
    bc_cluster_df <- merge(pseudoname_outliers_df, cl_mold, by = "cluster", all.x = T) %>%
      dplyr::select(-c("cluster")) %>%
      dplyr::rename("cluster" = new_cluster)
    bc_cluster_df$order <- max(as.numeric(rownames(hc_cluster_df))) + as.numeric(rownames(bc_cluster_df))
    final_order <- hc_cluster_df
    rownames(final_order) <- final_order$GOID
    score_sort <- final_order %>% dplyr::group_by(cluster) %>% dplyr::summarise(score_mean = mean(!!sym(score))) %>% dplyr::arrange(desc(score_mean))
    final_order <- merge(final_order, score_sort, by = "cluster", all.x = T)
    final_order$cluster <- factor(x = final_order$cluster, levels = score_sort$cluster)
    final_order <- final_order %>% dplyr::arrange(cluster, order)
    
    # orig order
    orig_order <- data.frame(GOID = colnames(w), orig_order = seq(1, nrow(w)), row.names = colnames(w))
    orig_order <- orig_order %>% dplyr::filter(!GOID %in% pseudoname_outliers_df$GOID)
    orig_order <- orig_order[final_order$GOID, ]
    
    mylist <- list("bc_cluster_df" = bc_cluster_df, "final_order" = final_order, "orig_order" = orig_order, "w" = w)
    
  } else if (outlier_detect == F) {
    
    # final order
    hc_cluster_df$order <- as.numeric(rownames(hc_cluster_df))
    max_cl <- max(hc_cluster_df$cluster)
    final_order <- hc_cluster_df
    rownames(final_order) <- final_order$GOID
    score_sort <- final_order %>% dplyr::group_by(cluster) %>% dplyr::summarise(score_mean = mean(!!sym(score))) %>% dplyr::arrange(desc(score_mean))
    final_order <- merge(final_order, score_sort, by = "cluster", all.x = T)
    final_order$cluster <- factor(x = final_order$cluster, levels = score_sort$cluster)
    final_order <- final_order %>% dplyr::arrange(cluster, order)
    
    # orig order
    orig_order <- data.frame(GOID = colnames(w), orig_order = seq(1, nrow(w)), row.names = colnames(w))
    orig_order <- orig_order[final_order$GOID, ]
    
    mylist <- list("final_order" = final_order, "orig_order" = orig_order, "w" = w)
    
  }
  
  return(mylist)
  
}


gorea_summary <- function(gorea_getorder_res = NULL,
                          score = NULL,
                          filename1 = "total_GOBP.xlsx",
                          filename2 = "representative_term.xlsx",
                          rep_df = NULL
) {
  
  final_order <- gorea_getorder_res$final_order
  orig_order <- gorea_getorder_res$orig_order
  
  # target GOID and ancestor GOID
  df <- merge(final_order, GOID_TERM, by = "GOID", all.x = T) %>% dplyr::arrange(cluster, order)
  c_orig <- df$cluster[(!duplicated(df$cluster))]
  c_df <- data.frame(clusters_orig = c_orig, clusters_anno = seq(1, length(c_orig)))
  df <- merge(df, c_df, by.x = "cluster", by.y = "clusters_orig", all.x = T)
  df$GOTERM <- gsub(x = df$GOTERM, pattern = "GOBP_", replacement = "") %>%
    gsub(x = ., pattern = "_", replacement = " ") %>% tolower()
  
  WriteXLS(df, ExcelFileName = filename1) # this dataframe includes GO terms and cluster information.
  
  # target
  df2 <- merge(df, rep_df, by.x = "cluster", by.y = "cluster", all.x = T)
  vals <- c("clusters_anno", "score_mean", "max_level", "ANCESTOR_GOID", "n_total_gobp", "n", "total_gobp", score)
  df2 <- df2 %>% dplyr::select(all_of(vals))
  df2 <- df2[!duplicated(df2[c("clusters_anno", "score_mean", "max_level", "ANCESTOR_GOID", "n_total_gobp")]),]
  df2 <- df2 %>% dplyr::filter(!is.na(ANCESTOR_GOID)) # filtering out outliers
  
  df3 <- data.frame()
  for (i in seq(1:nrow(df2))) {
    df2_tmp <- df2[i,]
    
    GOTERM_tmp <- GO.ANCESTOR.LEVEL %>%
      dplyr::filter(ANCESTOR_GOID == df2_tmp$ANCESTOR_GOID) %>%
      dplyr::slice(1) %>%
      dplyr::pull(ANCESTOR_GOTERM)
    df2_tmp$GOBP_ANCESTOR_GOTERM <- GOTERM_tmp
    df3 <- rbind(df3, df2_tmp)
  }
  
  df3 <- df3 %>% dplyr::arrange(clusters_anno, desc(n))
  
  # final results
  WriteXLS(df3, ExcelFileName = filename2) # this dataframe is the final results. put the table with final heatmap.
  
  mylist <- list(SummaryGO_1 = df, SummaryGO_2 = df3)
  
  return(mylist)
  
}



gorea_getHeatmap <- function(k_val = NULL,
                             gorea_getorder_res = NULL,
                             gorea_summary_res = NULL,
                             outlier_detect = T,
                             plot = T,
                             heatmap_filename = "test.png",
                             heatmap_width = 40,
                             heatmap_height = 30,
                             ancestor_annotation = T,
                             right_annotation_font_size = 10,
                             right_annotation_fontface = "bold",
                             cluster_font_size = 2,
                             top_ancestor_annotation = T,
                             top_ancestor = 3,
                             color = c("gold", "darkgrey")) {
  
  final_order <- gorea_getorder_res$final_order
  orig_order <- gorea_getorder_res$orig_order
  w <- gorea_getorder_res$w
  
  SummaryGO_1 <- gorea_summary_res$SummaryGO_1
  SummaryGO_2 <- gorea_summary_res$SummaryGO_2
  
  # frequency for clusters
  cl <- c()
  freq_cl <- c()
  
  for (anode in unique(final_order$cluster)) {
    
    cl_tmp <- as.numeric(table(final_order$cluster)[as.character(anode)])
    cl <- append(cl, cl_tmp)
    term_sum <- length(final_order$cluster)
    freq_cl_tmp <- cl_tmp/term_sum
    freq_cl <- append(freq_cl, freq_cl_tmp)
    
  }
  
  col_fun <- colorRamp2(seq(0, quantile(w[w > 0], 0.975), length = length(c("white","red"))), c("white","red"))
  
  if (plot == T) {
    
    png(heatmap_filename, width = heatmap_width, height = heatmap_height, units = "cm", res = 1200)
    
  }
  
  if (ancestor_annotation == T) {
    
    rep_anno_plot <- SummaryGO_2 %>%
      dplyr::arrange(clusters_anno, desc(n)) %>%
      dplyr::group_by(clusters_anno) %>%
      dplyr::filter(n == max(n))
    
    rep_anno_plot2 <- data.frame()
    
    for (c_tmp in sort(unique(rep_anno_plot$clusters_anno))) {
      
      rep_anno_plot_tmp <- rep_anno_plot %>%
        dplyr::filter(clusters_anno == c_tmp)
      
      if (rep_anno_plot_tmp[1,]$n != 1) {
        
        plot_anno_tmp <- paste(rep_anno_plot_tmp$GOBP_ANCESTOR_GOTERM, collapse = "\n")
        rep_anno_plot_tmp$plot_anno <- plot_anno_tmp 
        
        rep_anno_plot2 <- rbind(rep_anno_plot2, rep_anno_plot_tmp)
        
      } else if (rep_anno_plot_tmp[1,]$n == 1) {
        
        if (nrow(rep_anno_plot_tmp) < 4) {
          
          plot_anno_tmp <- paste(rep_anno_plot_tmp$GOBP_ANCESTOR_GOTERM, collapse = "\n")
          rep_anno_plot_tmp$plot_anno <- plot_anno_tmp 
          
          rep_anno_plot2 <- rbind(rep_anno_plot2, rep_anno_plot_tmp)
          
        } else {
          
          rep_anno_plot_tmp$plot_anno <- "" 
          
          rep_anno_plot2 <- rbind(rep_anno_plot2, rep_anno_plot_tmp)
          
        }
        
      }
      
    }
    
    rep_anno_plot2 <- rep_anno_plot2[!duplicated(rep_anno_plot2[c("clusters_anno", "n")]),]
    
    # if some clusters were removed by detecting outliers, k_val has to be changed
    k_val <- SummaryGO_1$clusters_anno %>% max()
    
    textbox_list <- list()
    for (i in 1:k_val) {
      name_tmp <- as.character(i)
      goids_tmp <- SummaryGO_1 %>% dplyr::filter(clusters_anno == i) %>% dplyr::pull(GOID) %>% unique()
      textbox_list[[name_tmp]] <- orig_order %>% dplyr::filter(GOID %in% goids_tmp) %>% dplyr::pull(orig_order)
    }
    
    anno_list <- list()
    for (i in 1:k_val) {
      name_tmp <- as.character(i)
      anno_text <- rep_anno_plot2 %>% dplyr::filter(clusters_anno == i) %>% dplyr::pull(plot_anno) %>% unique()
      df_tmp <- data.frame(anno_text, col = "black", fontsize = right_annotation_font_size, fontface = right_annotation_fontface)
      anno_list[[name_tmp]] <- df_tmp
    }
    
    p1 <- Heatmap(mat = w, col = col_fun, name = "Similarity",
                  row_order = orig_order$orig_order,
                  column_order = orig_order$orig_order,
                  row_title = NULL,
                  border_gp = gpar(col = "black", lwd = 2),
                  show_column_dend = F,
                  show_row_dend = F,
                  show_row_names = F,
                  show_column_names = F,
                  row_dend_width = unit(2, "cm"),
                  width = unit(15, "cm"),
                  height = unit(15, "cm"), 
                  #left_annotation = rowAnnotation(ggplot1 = anno_empty(height = unit(15, "cm"), width = unit(1, "cm")),
                  #                                show_annotation_name = FALSE, gp = gpar(col = "white", lwd = 2)),
                  right_annotation = rowAnnotation(ggplot1 = anno_empty(height = unit(15, "cm"), width = unit(0.7, "cm")),
                                                   show_annotation_name = FALSE, gp = gpar(col = "white", lwd = 2),
                                                   textbox = anno_textbox(textbox_list, anno_list),
                                                   gap = unit(1, "mm"))
    )
    
  } else if (ancestor_annotation == F) {
    
    p1 <- Heatmap(mat = w, col = col_fun, name = "Similarity",
                  row_order = orig_order$orig_order,
                  column_order = orig_order$orig_order,
                  row_title = NULL,
                  border_gp = gpar(col = "black", lwd = 2),
                  show_column_dend = F,
                  show_row_dend = F,
                  show_row_names = F,
                  show_column_names = F,
                  row_dend_width = unit(2, "cm"),
                  width = unit(15, "cm"),
                  height = unit(15, "cm"),
                  left_annotation = rowAnnotation(ggplot1 = anno_empty(height = unit(15, "cm"), width = unit(0.7, "cm")),
                                                  show_annotation_name = FALSE, gp = gpar(col = "white")),
                  right_annotation = rowAnnotation(ggplot2 = anno_empty(border = FALSE, height = unit(15, "cm"), width = unit(1, "cm")),
                                                   show_annotation_name = FALSE, gp = gpar(col = "white")
                  ))
    
  }
  
  if (top_ancestor_annotation == T) {
    
    res <- GO.ANCESTOR.LEVEL %>% dplyr::filter(TARGET_GOID %in% gorea_summary_res$SummaryGO_1$GOID) %>%
      dplyr::filter(ANCESTOR_GOID %in% GO.LARGE$ANCESTOR_GOID) %>%
      dplyr::group_by(ANCESTOR_GOTERM, ANCESTOR_LEVEL, ANCESTOR_GOID) %>% dplyr::summarise(n=n()) %>% dplyr::arrange(desc(n))
    
    res_filtered <- res %>% head(n = top_ancestor)
    GO.ANCESTOR.LEVEL_subset <- GO.ANCESTOR.LEVEL %>% dplyr::filter(ANCESTOR_LEVEL %in% c(1,2))
    
    df1 <- final_order
    df2 <- merge(final_order, GO.ANCESTOR.LEVEL_subset, by.x = "GOID", by.y = "TARGET_GOID", all.x = T)
    
    large_go <- c()
    
    for (i in 1:nrow(res_filtered)) {
      mycolname <- res_filtered[i,]$ANCESTOR_GOTERM
      child_terms <- df2 %>% dplyr::filter(ANCESTOR_GOID %in% res_filtered[i,]$ANCESTOR_GOID) %>% dplyr::pull(GOID)
      df1 <- df1 %>% dplyr::mutate(!!sym(mycolname) :=  case_when(GOID %in% child_terms ~ "child",
                                                                  T ~ "X"))
      n_child <- df1 %>% dplyr::filter(!!sym(mycolname) == "child") %>% nrow()
      pt_child <- round((n_child/nrow(df1)) * 100, 2)
      mycolname2 <- paste0(mycolname, " (", pt_child, "%)")
      df1 <- df1 %>% dplyr::rename(!!sym(mycolname2) := !!sym(mycolname))
      
      large_go <- c(large_go, mycolname2)
      
    }
    
    rownames(df1) <- df1$GOID
    df1 <- df1[rownames(w), ]
    
    df1 <- df1 %>% dplyr::select(all_of(large_go))
    df1 <- df1[,rev(colnames(df1))]
    names(color) <- c("child", "X")
    
    collist <- list()
    for (col_tmp in colnames(df1)) {
      collist[[col_tmp]] <- color
    }
    
    ht1 <- HeatmapAnnotation(df = df1,
                             gap = unit(1, "mm"),
                             col = collist, border = F,
                             annotation_name_gp = gpar(fontsize = 10)
    )
    
  } else if (top_ancestor_annotation == F) {
    
    ht1 <- NULL
    
  }
  
  if (top_ancestor_annotation == T) {
    
    print(ht1 %v% p1)
    
    for (col_tmp in colnames(df1)) {
      
      decorate_annotation(col_tmp, {
        grid.rect(gp = gpar(lwd = 2, col = "black", fill = NA)) 
      })
      
    }
    
  } else if (top_ancestor_annotation == F) {
    
    print(p1)
    
  }
  
  # make lines 
  freq_cl_cumsum <- cumsum(freq_cl)
  freq_cl_cumsum <- freq_cl_cumsum[-length(freq_cl_cumsum)]
  
  gap_filler <- 0.0001
  
  decorate_heatmap_body("Similarity", {
    for (freq_cl in freq_cl_cumsum) {
      grid.lines(c(freq_cl + gap_filler, freq_cl + gap_filler), c(0, 1), gp = gpar(lty = 1, lwd = 1.5))
    }
  })
  
  freq_cl_cumsum_rev <- cumsum(rev(freq_cl))
  freq_cl_cumsum_rev <- freq_cl_cumsum_rev[-length(freq_cl_cumsum_rev)]
  
  decorate_heatmap_body("Similarity", {
    for (freq_cl in freq_cl_cumsum_rev) {
      grid.lines(c(1, 0), c(freq_cl + gap_filler, freq_cl + gap_filler), gp = gpar(lty = 1, lwd = 1.5))
    }
  })
  
  # make left cluster box
  coor_df <- data.frame(x1 = 0, x2 = 5,
                        y1 = c(0,freq_cl_cumsum_rev), y2 = c(freq_cl_cumsum_rev, 1))
  
  total_clusters <- final_order %>% dplyr::group_by(cluster) %>% dplyr::summarise(n=n()) %>% dplyr::select(cluster) %>% table() %>% names()
  
  if (outlier_detect == T) {
    p <- ggplot() +
      geom_rect(data = coor_df,
                mapping = aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2, fill = y1),
                color =  "black",
                alpha = 1) +
      scale_fill_gradient(low = "ivory", high = "green4") +
      scale_x_continuous(expand = c(0,0)) +
      scale_y_continuous(expand = c(0,0)) +
      theme_void() +
      theme(panel.background = element_blank(),
            panel.grid.major= element_blank(),
            panel.grid.minor = element_blank(),
            axis.line = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank(),
            axis.text = element_blank(),
            legend.position = "none",
            plot.background = element_blank()) +
      labs(x=NULL, y=NULL) +
      geom_text(data = coor_df,
                aes(x=(x1+x2)/2, y=y1+((y2-y1)/2), label = rev(c(paste0(sort(as.numeric(total_clusters)))))), size = cluster_font_size, fontface = "bold")
    
  } else if (outlier_detect == F) {
    
    p <- ggplot() +
      geom_rect(data = coor_df,
                mapping = aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2, fill = y1),
                color =  "black",
                alpha = 1) +
      scale_fill_gradient(low = "ivory", high = "green4") +
      scale_x_continuous(expand = c(0,0)) +
      scale_y_continuous(expand = c(0,0)) +
      theme_void() +
      theme(panel.background = element_blank(),
            panel.grid.major= element_blank(),
            panel.grid.minor = element_blank(),
            axis.line = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank(),
            axis.text = element_blank(),
            legend.position = "none",
            plot.background = element_blank()) +
      labs(x=NULL, y=NULL) +
      geom_text(data = coor_df,
                aes(x=(x1+x2)/2, y=y1+((y2-y1)/2), label = rev(c(paste0(sort(as.numeric(total_clusters)))))), size = cluster_font_size, fontface = "bold")
    
  }
  
  # merge plot
  decorate_annotation("ggplot1", {
    vp = current.viewport()$name
    print(p, vp = vp)
  })
  
  decorate_annotation("ggplot1", {
    grid.rect(gp = gpar(lwd = 2, col = "black", fill = NA)) 
  })
  
  if (plot == T) {
    dev.off()
  }
  
}

gorea <- function(input = input_df, k_val = NULL, godata_GO = NULL,
                  cutoff = 0.84,
                  outlier_detect = T,
                  min_cluster = 3,
                  representative_term_level_cutoff = 1, GO_explain = 3,
                  score = NULL, filename1 = NULL, filename2 = NULL, heatmap_filename = NULL,
                  plot = T,
                  heatmap_width = 40, heatmap_height = 30,
                  ancestor_annotation = T,
                  right_annotation_font_size = 10,
                  cluster_font_size = 4,
                  top_ancestor_annotation = T,
                  top_ancestor = 3,
                  color = c("gold", "darkgrey") ) {
  
  w <- gorea_sim_mat(input = input,
                     sim_method = "Rel",
                     godata_GO = godata_GO)
  
  
  gorea_clustering_res <- gorea_clustering(cutoff = cutoff,
                                           w = w,
                                           k_val = k_val,
                                           input = input_df,
                                           outlier_detect = outlier_detect,
                                           min_cluster = min_cluster)
  
  gorea_representative_terms_res <- gorea_representative_terms(gorea_clustering_res = gorea_clustering_res,
                                                               k_val = k_val,
                                                               outlier_detect = outlier_detect,
                                                               representative_term_level_cutoff = representative_term_level_cutoff,
                                                               GO_explain = GO_explain)
  
  gorea_getorder_res <- gorea_getorder(w = w,
                                       gorea_clustering_res = gorea_clustering_res,
                                       input = input_df,
                                       score = score,
                                       outlier_detect = outlier_detect)
  
  gorea_summary_res <- gorea_summary(gorea_getorder_res = gorea_getorder_res,
                                     score = score,
                                     filename1 = filename1,
                                     filename2 = filename2,
                                     rep_df = gorea_representative_terms_res)
  
  
  gorea_getHeatmap_res <- gorea_getHeatmap(k_val = k_val,
                                           gorea_getorder_res = gorea_getorder_res,
                                           gorea_summary_res = gorea_summary_res,
                                           outlier_detect = outlier_detect,
                                           plot = plot,
                                           heatmap_filename = heatmap_filename,
                                           heatmap_width = heatmap_width,
                                           heatmap_height = heatmap_height,
                                           ancestor_annotation = ancestor_annotation,
                                           right_annotation_font_size = right_annotation_font_size,
                                           right_annotation_fontface = "bold",
                                           cluster_font_size = cluster_font_size,
                                           top_ancestor_annotation = top_ancestor_annotation,
                                           top_ancestor = top_ancestor,
                                           color = color)
  
  final_list <- list(w = w,
                     gorea_clustering_res = gorea_clustering_res,
                     gorea_representative_terms_res = gorea_representative_terms_res,
                     gorea_getorder_res = gorea_getorder_res,
                     gorea_summary_res = gorea_summary_res)
  
  return(final_list)
}