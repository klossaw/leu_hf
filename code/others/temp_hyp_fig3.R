pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes",
          "SingleCellExperiment", "RColorBrewer", "vroom", "jhtools", "glue", "jhuanglabHyperion",
          "openxlsx", "ggsci", "patchwork", "cowplot", "tidyverse", "dplyr", "rstatix")
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}
project <- "hyperion"
dataset <- "qzhang"
species <- "human"
workdir <- glue::glue("~/projects/{project}/analysis/{dataset}/{species}/figures/manuscripts/fig3")
setwd(workdir)


sce <- readr::read_rds("/cluster/home/jhuang/projects/hyperion/analysis/qzhang/human/steinbock/rds/combined/all_anno.rds")
meta_clu <- readr::read_rds("/cluster/home/yjliu_jh/share/meta_clu3.rds")
col_new <- colData(sce) %>% as_tibble()
col_meta <- unique(meta_clu$meta_color) %>% `names<-`(unique(meta_clu$meta_short))
metadata <- metadata(sce)
metadata[["differentiation_degree"]]$diff_degree <- factor(metadata[["differentiation_degree"]]$diff_degree, levels = c("low", "middle", "high"))
metadata[["stages_tme"]]$stage <- factor(metadata[["stages_tme"]]$stage)
levels(metadata[["stages_tme"]]$stage) = c("RPC", "BRPC_LAPC", "MPC")


keycols_list <- c(list("stage"), list("treatment_type"), list("os_group_24"))
compare_groups <- c("stages_tme", "neoadj_vs_direct_surgery", "os_analysis")
names(keycols_list) <- compare_groups
meta_cluster_names <- rev(c("Macro_M1_M2", "Immune_enriched", "Macro_M1_like", "Immune_mixed",
                            "Stroma_mCAF", "Stroma_Macro",
                            "Stroma_CAF", "Tumor_boundary", "Bulk_tumor"))
total_group <- col_new %>% dplyr::select(c(sample_id, sample_tiff_id, patient_id, community_name, stype2, meta_cluster)) %>%
  na.omit() |> dplyr::filter(meta_cluster %in% meta_cluster_names)
#---------------------------------------- main figures fig3a----------------------------------------
plist <- list()
# total meta cluster ----------------------------------------
total_group_pro <- total_group %>% group_by(meta_cluster, sample_tiff_id) %>% summarise(pro = n()/nrow(total_group))
total_group_pro$meta_cluster <- fct_relevel(total_group_pro$meta_cluster, meta_cluster_names)
p1.1 <- ggplot(total_group_pro, aes(x = pro, y = meta_cluster, fill = meta_cluster)) +
  geom_col(width = 0.8) +
  theme_classic() +
  scale_fill_manual(values = col_meta,
                    labels = vars(meta_cluster)) +
  theme(strip.placement  = "outside",
        panel.spacing    = unit(3, "points"),
        strip.background = element_blank(),
        strip.text       = element_text(face = "bold", size = 5),
        axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10),
        legend.position="none")
ggsave("community_total_metacluster.pdf", p1, width = 3, height = 4)
plist[["total"]] <- p1.1


# -------- all total for alteration? and for new p1  --------
 all_total <- total_group_pro %>% group_by(meta_cluster) %>% summarise(allpro = sum(pro))


##   cell count proportion? cluster numbers proportion? or alteration? 


# 1 ===== count ======
plist <- list()
plist[["total"]] <- p1.1

# compare groups ----------------------------------------
for(compare_group in compare_groups){
  cli::cli_h1(compare_group)
  groups <- metadata[[compare_group]]
  dt_groups <- total_group %>% dplyr::filter(patient_id %in% groups$patient_id)
  if (compare_group == "chemo_outcome_before") {
    dt_groups <- dt_groups %>% dplyr::filter(stype2 %in% c("before_chemo","puncture_pdac"))
  } else if (compare_group %notin% c("paired_after_neoadj", "mpc_primary_metastasis", "tumor_para")) {
    dt_groups <- dt_groups %>% dplyr::filter(stype2 %in% c("after_chemo","tumor"))
  }
  dat_plot <- dt_groups %>%
    group_by(meta_cluster, sample_tiff_id) |>
    summarise(nc = n()) |>
    group_by(sample_tiff_id) |>
    dplyr::mutate(nt = sum(nc)) |>
    dplyr::mutate(pro = nc/nt) |> ungroup()
  dat_plot$meta_cluster <- factor(dat_plot$meta_cluster, levels = meta_cluster_names)
  dat_plot <- dat_plot |> tidyr::complete(meta_cluster, sample_tiff_id, fill = list(pro = 0)) |>
    dplyr::mutate(pro = sd_censor(pro, range = c(-3, 3)))
  metainfo <- col_new[, c("sample_tiff_id", "patient_id")] |> distinct()
  dat_plot <- left_join(dat_plot, metainfo, by = "sample_tiff_id") |>
    left_join(groups, by = "patient_id") %>% na.omit() |> distinct()
  # dat_plot <- left_join(dat_plot, groups, by = "patient_id") %>% na.omit() |> distinct()
  #dat_plot <- dat_plot |> dplyr::filter(pro < 0.1)
  gp_key <- keycols_list[[compare_group]]
  cli::cli_alert_info(gp_key)
  dat_plot[[gp_key]] <- factor(dat_plot[[gp_key]])
  dat_plot$meta_cluster <- fct_relevel(dat_plot$meta_cluster, meta_cluster_names)
  p <- ggboxplot(dat_plot, x = "meta_cluster", y = "pro", fill = gp_key, outlier.shape = NA,
                 palette = pal_nejm("default")(3), xlab = NULL) + theme_classic() +
    theme(strip.placement  = "outside",
          panel.spacing    = unit(3, "points"),
          strip.background = element_blank(),
          strip.text       = element_text(face = "bold", size = 5),
          axis.text.x = element_text(size = 10), axis.text.y = element_blank(),
          legend.position="right") +
    labs(x= NULL, y = NULL)
  exp1 <- expr(pro ~ !!ensym(gp_key))
  stat_test <- dat_plot |>
    group_by(meta_cluster) %>% rstatix::t_test(eval(exp1), p.adjust.method = "none")
  stat_test$p.adj.signif <- ifelse(stat_test$p >= 0.05, "ns", 
                                   ifelse(stat_test$p >= 0.01, "*", 
                                          ifelse(stat_test$p >= 0.001, "**", 
                                                 ifelse(stat_test$p >= 0.05, "***", "****"))))
  stat_test <- stat_test %>%
    add_xy_position(x = "meta_cluster", dodge = 0.8)
  p1 <- p +
    stat_pvalue_manual(
      stat_test, tip.length = 0.01, hide.ns = T, label = "p.adj.signif",
      coord.flip = TRUE
    ) + coord_flip()
  ggsave(glue("metacluster_{compare_group}.pdf"), p1, width = 6, height = 6)
  plist[[compare_group]] <- p1
}



pc <- (plist[[1]]|plist[[2]]|plist[[3]]|plist[[4]]) + plot_layout(guides = 'collect', widths = c(2, 2.5, 2.5, 2.5))
ggsave("fig3a1x.pdf", pc, width = 16, height = 6)





# 1 ===== count/mcs ======
plist <- list()
plist[["total"]] <- p1.1

# compare groups ----------------------------------------
for(compare_group in compare_groups){
  cli::cli_h1(compare_group)
  groups <- metadata[[compare_group]]
  dt_groups <- total_group %>% dplyr::filter(patient_id %in% groups$patient_id)
  if (compare_group == "chemo_outcome_before") {
    dt_groups <- dt_groups %>% dplyr::filter(stype2 %in% c("before_chemo","puncture_pdac"))
  } else if (compare_group %notin% c("paired_after_neoadj", "mpc_primary_metastasis", "tumor_para")) {
    dt_groups <- dt_groups %>% dplyr::filter(stype2 %in% c("after_chemo","tumor"))
  }
  dt_groups <- unique(dt_groups)
  dat_plot <- dt_groups %>%
    group_by(meta_cluster, sample_id) |>
    summarise(nc = n()) |>
    group_by(sample_id) |>
    dplyr::mutate(nt = sum(nc)) |>
    dplyr::mutate(pro = nc/nt) |> ungroup()
  
  #dat_plot <- dt_groups %>%
   # group_by(meta_cluster, community_name, sample_id) |>
   # summarise(nc = n()) |>
   # group_by(meta_cluster, sample_id) |>
   # dplyr::mutate(nt = sum(nc)) |>
  #  dplyr::mutate(pro = nc/nt) |> ungroup()
  dat_plot$meta_cluster <- factor(dat_plot$meta_cluster, levels = meta_cluster_names)
  dat_plot <- dat_plot |> tidyr::complete(meta_cluster, sample_id, fill = list(pro = 0)) |>
    dplyr::mutate(pro = sd_censor(pro, range = c(-3, 3)))
  metainfo <- col_new[, c("sample_id", "patient_id")] |> distinct()
  dat_plot <- left_join(dat_plot, metainfo, by = "sample_id") |>
    left_join(groups, by = "patient_id") %>% na.omit() |> distinct()
  # dat_plot <- left_join(dat_plot, groups, by = "patient_id") %>% na.omit() |> distinct()
  #dat_plot <- dat_plot |> dplyr::filter(pro < 0.1)
  gp_key <- keycols_list[[compare_group]]
  cli::cli_alert_info(gp_key)
  dat_plot[[gp_key]] <- factor(dat_plot[[gp_key]])
  dat_plot$meta_cluster <- fct_relevel(dat_plot$meta_cluster, meta_cluster_names)
  p <- ggboxplot(dat_plot, x = "meta_cluster", y = "pro", fill = gp_key, outlier.shape = NA,
                 palette = pal_nejm("default")(3), xlab = NULL) + theme_classic() +
    theme(strip.placement  = "outside",
          panel.spacing    = unit(3, "points"),
          strip.background = element_blank(),
          strip.text       = element_text(face = "bold", size = 5),
          axis.text.x = element_text(size = 10), axis.text.y = element_blank(),
          legend.position="right") +
    labs(x= NULL, y = NULL)
  exp1 <- expr(pro ~ !!ensym(gp_key))
  stat_test <- dat_plot |>
    group_by(meta_cluster) %>% rstatix::t_test(eval(exp1), p.adjust.method = "none")
  stat_test <- stat_test %>%
    add_xy_position(x = "meta_cluster", dodge = 0.8)
  p1 <- p +
    stat_pvalue_manual(
      stat_test, tip.length = 0.01, hide.ns = T,
      coord.flip = TRUE
    ) + coord_flip()
  ggsave(glue("metacluster_{compare_group}.pdf"), p1, width = 6, height = 6)
  plist[[compare_group]] <- p1
}



pc <- (plist[[1]]|plist[[2]]|plist[[3]]|plist[[4]]) + plot_layout(guides = 'collect', widths = c(2, 2.5, 2.5, 2.5))
ggsave("fig3a3.pdf", pc, width = 16, height = 6)





temp <- left_join(total_group, all_total)
dat_plot <- dt_groups %>%
  group_by(meta_cluster, sample_id, allpro) |>
  summarise(nc = n()) |>
  group_by(sample_id) |>
  dplyr::mutate(nt = sum(nc)) |>
  dplyr::mutate(pro = nc/nt) |> ungroup()




# 1 ===== alteration proportion ======
plist <- list()
plist[["total"]] <- p1.1

# compare groups ----------------------------------------
for(compare_group in compare_groups){
  cli::cli_h1(compare_group)
  groups <- metadata[[compare_group]]
  dt_groups <- total_group %>% dplyr::filter(patient_id %in% groups$patient_id)
  if (compare_group == "chemo_outcome_before") {
    dt_groups <- dt_groups %>% dplyr::filter(stype2 %in% c("before_chemo","puncture_pdac"))
  } else if (compare_group %notin% c("paired_after_neoadj", "mpc_primary_metastasis", "tumor_para")) {
    dt_groups <- dt_groups %>% dplyr::filter(stype2 %in% c("after_chemo","tumor"))
  }
  dt_groups <- left_join(dt_groups, all_total)
  dat_plot <- dt_groups %>%
    group_by(meta_cluster, sample_id, allpro) |>
    summarise(nc = n()) |>
    group_by(sample_id) |>
    dplyr::mutate(nt = sum(nc)) |>
    dplyr::mutate(prox = nc/nt) |> 
    dplyr::mutate(pro = ((prox - allpro) / allpro)) |> ungroup()
  
  dat_plot$meta_cluster <- factor(dat_plot$meta_cluster, levels = meta_cluster_names)
  dat_plot <- dat_plot |> tidyr::complete(meta_cluster, sample_id, fill = list(pro = 0)) |>
     dplyr::mutate(pro = sd_censor_max(pro, range = c(-3, 3)))
  metainfo <- col_new[, c("sample_id", "patient_id")] |> distinct()
  dat_plot <- left_join(dat_plot, metainfo, by = "sample_id") |>
    left_join(groups, by = "patient_id") %>% na.omit() |> distinct()
  # dat_plot <- left_join(dat_plot, groups, by = "patient_id") %>% na.omit() |> distinct()
  #dat_plot <- dat_plot |> dplyr::filter(pro < 0.1)
  gp_key <- keycols_list[[compare_group]]
  cli::cli_alert_info(gp_key)
  dat_plot[[gp_key]] <- factor(dat_plot[[gp_key]])
  dat_plot$meta_cluster <- fct_relevel(dat_plot$meta_cluster, meta_cluster_names)
  p <- ggboxplot(dat_plot, x = "meta_cluster", y = "pro", fill = gp_key, outlier.shape = NA,
                 palette = pal_nejm("default")(3), xlab = NULL) + theme_classic() +
    theme(strip.placement  = "outside",
          panel.spacing    = unit(3, "points"),
          strip.background = element_blank(),
          strip.text       = element_text(face = "bold", size = 5),
          axis.text.x = element_text(size = 10), axis.text.y = element_blank(),
          legend.position="right") +
    labs(x= NULL, y = NULL)
  exp1 <- expr(pro ~ !!ensym(gp_key))
  stat_test <- dat_plot |>
    group_by(meta_cluster) %>% rstatix::t_test(eval(exp1), p.adjust.method = "none")
  stat_test <- stat_test %>%
    add_xy_position(x = "meta_cluster", dodge = 0.8)
  p1 <- p +
    stat_pvalue_manual(
      stat_test, tip.length = 0.01, hide.ns = T,
      coord.flip = TRUE
    ) + coord_flip()
  ggsave(glue("metacluster_{compare_group}.pdf"), p1, width = 6, height = 6)
  plist[[compare_group]] <- p1
}



pc <- (plist[[1]]|plist[[2]]|plist[[3]]|plist[[4]]) + plot_layout(guides = 'collect', widths = c(2, 2.5, 2.5, 2.5))
ggsave("fig3a5.pdf", pc, width = 16, height = 6)















# compare groups -x ----------------------------------------
plist <- list()
plist[["total"]] <- p1.1

for(compare_group in compare_groups){
  cli::cli_h1(compare_group)
  groups <- metadata[[compare_group]]
  dt_groups <- total_group %>% dplyr::filter(patient_id %in% groups$patient_id)
  if (compare_group == "chemo_outcome_before") {
    dt_groups <- dt_groups %>% dplyr::filter(stype2 %in% c("before_chemo","puncture_pdac"))
  } else if (compare_group %notin% c("paired_after_neoadj", "mpc_primary_metastasis", "tumor_para")) {
    dt_groups <- dt_groups %>% dplyr::filter(stype2 %in% c("after_chemo","tumor"))
  }
  dat_plot <- dt_groups %>%
    dplyr::mutate(tt = n()) |>
    group_by(meta_cluster, sample_id, tt) |>
    summarise(nc = n()) |>
    group_by(meta_cluster) |>
    dplyr::mutate(nt = sum(nc)) |>
    dplyr::mutate(pro = nc) |> ungroup()
  dat_plot$meta_cluster <- factor(dat_plot$meta_cluster, levels = meta_cluster_names)
  dat_plot <- dat_plot |> tidyr::complete(meta_cluster, sample_id, fill = list(pro = 0)) |>
    dplyr::mutate(pro = sd_censor(pro, range = c(-3, 3)))
  metainfo <- col_new[, c("sample_id", "patient_id")] |> distinct()
  dat_plot <- left_join(dat_plot, metainfo, by = "sample_id") |>
    left_join(groups, by = "patient_id") %>% na.omit() |> distinct()
  # dat_plot <- left_join(dat_plot, groups, by = "patient_id") %>% na.omit() |> distinct()
  #dat_plot <- dat_plot |> dplyr::filter(pro < 0.1)
  gp_key <- keycols_list[[compare_group]]
  cli::cli_alert_info(gp_key)
  dat_plot[[gp_key]] <- factor(dat_plot[[gp_key]])
  dat_plot$meta_cluster <- fct_relevel(dat_plot$meta_cluster, meta_cluster_names)
  p <- ggboxplot(dat_plot, x = "meta_cluster", y = "pro", fill = gp_key, outlier.shape = NA,
                 palette = pal_nejm("default")(3), xlab = NULL) + theme_classic() +
    theme(strip.placement  = "outside",
          panel.spacing    = unit(3, "points"),
          strip.background = element_blank(),
          strip.text       = element_text(face = "bold", size = 5),
          axis.text.x = element_text(size = 10), axis.text.y = element_blank(),
          legend.position="right") +
    labs(x= NULL, y = NULL)
  exp1 <- expr(pro ~ !!ensym(gp_key))
  stat_test <- dat_plot |>
    group_by(meta_cluster) %>% rstatix::t_test(eval(exp1), p.adjust.method = "none")
  stat_test <- stat_test %>%
    add_xy_position(x = "meta_cluster", dodge = 0.8)
  p1 <- p +
    stat_pvalue_manual(
      stat_test, tip.length = 0.01, hide.ns = T,
      coord.flip = TRUE
    ) + coord_flip()
  ggsave(glue("metacluster_{compare_group}.pdf"), p1, width = 6, height = 6)
  plist[[compare_group]] <- p1
}

pc <- (plist[[1]]|plist[[2]]|plist[[3]]|plist[[4]]) + plot_layout(guides = 'collect', widths = c(2, 2.5, 2.5, 2.5))
ggsave("fig3a7.pdf", pc, width = 16, height = 6)



