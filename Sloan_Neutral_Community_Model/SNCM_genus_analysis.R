################################################################################
## SNCM Genus-Level Analysis (Per Niche × Treatment)
## Refined for GitHub repository
################################################################################

# -------------------------------
# Setup
# -------------------------------
# Packages
library(phyloseq)
library(MicEco)
library(dplyr)
library(ggplot2)

# Paths (relative)
ps_file <- "data/phyloseq_object.rds"
meta_file <- "data/metadata.csv"
out_dir <- "out"

# Create output folder if missing
if(!dir.exists(out_dir)) dir.create(out_dir)

# -------------------------------
# Load data
# -------------------------------
ps <- readRDS(ps_file)
meta <- read.csv(meta_file, row.names = 1)
sample_data(ps) <- sample_data(meta)

# Clean phyloseq: remove mitochondria & chloroplast
ps_clean <- subset_taxa(ps,
                        !(Family %in% c("Mitochondria")) &
                        !(Order %in% c("Chloroplast")) &
                        !(Class %in% c("Chloroplast")))

# Aggregate to Genus level
ps_genus <- tax_glom(ps_clean, taxrank = "Genus")

# -------------------------------
# Split by Treatment
# -------------------------------
meta_df <- data.frame(sample_data(ps_genus))
treatments <- unique(meta_df$Treatment)
ps_treat_list <- lapply(treatments, function(tr) {
  prune_samples(rownames(subset(meta_df, Treatment == tr)), ps_genus)
})
names(ps_treat_list) <- treatments

# -------------------------------
# SNCM function
# -------------------------------
run_sncm <- function(ps_sub) {
  otu <- as.data.frame(otu_table(ps_sub))
  if(taxa_are_rows(ps_sub)) otu <- t(otu)
  otu <- as.data.frame(lapply(otu, as.numeric))
  
  # Filter low-prevalence taxa
  otu <- otu[, colSums(otu > 0) >= 3, drop = FALSE]
  if(ncol(otu) < 5) return(NULL)
  
  fit <- tryCatch(neutral.fit(otu), error=function(e) NULL)
  if(is.null(fit)) return(NULL)
  
  # Manual pseudo-R²
  obs <- fit[[2]]$freq
  pred <- fit[[2]]$freq.pred
  fit$model$manual_R2 <- 1 - sum((obs - pred)^2)/sum((obs - mean(obs))^2)
  
  # Classify taxa
  fit[[2]]$fit_class <- ifelse(fit[[2]]$freq > fit[[2]]$Upper, "Above",
                               ifelse(fit[[2]]$freq < fit[[2]]$Lower, "Below", "Neutral"))
  return(fit)
}

# -------------------------------
# Run SNCM for each treatment
# -------------------------------
sncm_list <- lapply(ps_treat_list, run_sncm)

# Combined dataset
fit_combined <- run_sncm(ps_genus)
sncm_list[["Combined"]] <- fit_combined

# -------------------------------
# Combine taxa-level SNCM results
# -------------------------------
sncm_taxa_all <- lapply(names(sncm_list), function(nm) {
  fit <- sncm_list[[nm]]
  if(is.null(fit)) return(NULL)
  df <- fit[[2]]
  df$Genus <- rownames(df)
  df$Treatment <- nm
  tax_df <- as.data.frame(tax_table(ps_genus))
  tax_df$Genus <- rownames(tax_df)
  merge(df, tax_df, by="Genus", all.x=TRUE)
})
sncm_taxa_df <- do.call(rbind, sncm_taxa_all)
sncm_taxa_df <- na.omit(sncm_taxa_df)
write.csv(sncm_taxa_df, file.path(out_dir, "SNCM_AllNiches_Genus_fitclass.csv"), row.names=FALSE)

# -------------------------------
# m + R² summary
# -------------------------------
sncm_m_values <- lapply(names(sncm_list), function(nm) {
  fit <- sncm_list[[nm]]
  if(is.null(fit)) return(data.frame(Treatment=nm, m=NA, R2=NA))
  data.frame(Treatment=nm, m=fit[[1]][1,"m"], R2=fit$model$manual_R2)
})
m_summary <- do.call(rbind, sncm_m_values)
write.csv(m_summary, file.path(out_dir, "SNCM_AllNiches_Genus_m_R2_summary.csv"), row.names=FALSE)

# -------------------------------
# Faceted genus-level plots per treatment
# -------------------------------
for(tr in unique(sncm_taxa_df$Treatment)){
  df_plot <- sncm_taxa_df %>% filter(Treatment == tr)
  p <- ggplot(df_plot, aes(x=p, y=freq)) +
    geom_point(aes(color=fit_class), alpha=0.7, size=2) +
    geom_line(aes(y=freq.pred), color="black") +
    geom_line(aes(y=Upper), linetype="dashed", color="black") +
    geom_line(aes(y=Lower), linetype="dashed", color="black") +
    scale_x_log10() +
    scale_y_continuous(limits=c(0,1)) +
    scale_color_manual(values=c("Above"="#E64B35","Neutral"="#4DBBD5","Below"="#00A087")) +
    theme_minimal(base_size=14) +
    labs(title=paste("SNCM (Genus) -", tr),
         x="Mean Relative Abundance (log10 scale)",
         y="Occurrence Frequency",
         color="Fit Category")
  
  ggsave(file.path(out_dir, paste0("SNCM_", tr, "_Genus_faceted.png")), p, width=8, height=6)
}

# -------------------------------
# 8-Panel combined plot with bold m + R²
# -------------------------------
sncm_taxa_df <- sncm_taxa_df %>%
  mutate(Panel = Treatment)

panel_stats <- m_summary %>%
  mutate(Panel = Treatment) %>%
  filter(Panel %in% unique(sncm_taxa_df$Panel)) %>%
  left_join(
    sncm_taxa_df %>%
      group_by(Panel) %>% summarize(xpos = min(p, na.rm=TRUE)*0.5),
    by="Panel"
  )

p_8panel <- ggplot(sncm_taxa_df, aes(x=p, y=freq)) +
  geom_point(aes(color=fit_class), alpha=0.7, size=2) +
  geom_line(aes(y=freq.pred), color="black") +
  geom_line(aes(y=Upper), linetype="dashed", color="black") +
  geom_line(aes(y=Lower), linetype="dashed", color="black") +
  scale_x_log10() +
  scale_y_continuous(limits=c(0,1)) +
  scale_color_manual(values=c("Above"="#E64B35","Neutral"="#4DBBD5","Below"="#00A087")) +
  facet_wrap(~Panel, ncol=2) +
  geom_text(
    data=panel_stats,
    aes(x=xpos, y=1.05, label=paste0("m=", round(m,3), ", R²=", round(R2,2))),
    inherit.aes=FALSE, hjust=0, vjust=0, fontface="bold", size=4
  ) +
  theme_minimal(base_size=14) +
  labs(title="SNCM (Genus) per Treatment",
       x="Mean Relative Abundance (log10 scale)",
       y="Occurrence Frequency",
       color="Fit Category")

ggsave(file.path(out_dir, "SNCM_AllNiches_8Panel_mR2_bold.png"), p_8panel, width=12, height=16)

# -------------------------------
# Save session info for reproducibility
# -------------------------------
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))

################################################################################
# END OF SCRIPT
################################################################################
