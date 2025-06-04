if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")

library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)


results_df <- read.csv("../results/merged_results/krogan_lab_hpidb_results_deduplicated.csv")
results_df_high_iptm <- results_df[results_df$iptm >= 0.5, ]
human_proteins <- results_df_high_iptm$human_protein

ids_df <- bitr(human_proteins, fromType = "UNIPROT",
                         toType = "SYMBOL",
                         OrgDb = org.Hs.eg.db)
human_genes <- ids_df$SYMBOL

ego <- enrichGO(gene = results_df_high_iptm$human_protein,
                OrgDb = org.Hs.eg.db,
                keyType = "UNIPROT",
                ont = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05, # previously 0.01. 0.05 is the default
                qvalueCutoff  = 0.2) # previously 0.05. 0.2 is the default
ego <- setReadable(ego, OrgDb = org.Hs.eg.db)
#ego_summary = summary(ego)
#barplot(ego, showCategory=20) 
ego_simplified <- simplify(ego, cutoff=0.7, by="p.adjust", select_fun=min)
ego_simpl_summary = summary(ego_simplified)
plt <- barplot(ego_simplified, fontsize=12, showCategory=10)
ggsave('go_enrichment_simplified_high_iptm.svg',
       plot=plt, device='svg', width=6, height=5,
       path='../results/plots_apms_hpidb/')
