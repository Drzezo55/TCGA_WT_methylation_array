
# journey through specific DMRs
GATA3_cpgs <- ann450k$Name[grepl("GATA3", ann450k$UCSC_RefGene_Name)]
GATA3_cpgs
length(GATA3_cpgs)
beta_subset <- beta_mat[row.names(beta_mat) %in% GATA3_cpgs, ]
annotation_col <- data.frame(Survival = metadata$ten_year_survival)
rownames(annotation_col) <- metadata$barcode
ordered_samples <- metadata$barcode[order(metadata$ten_year_survival)]
beta_subset_ordered <- beta_subset[, ordered_samples]
annotation_col <- annotation_col[ordered_samples, , drop = FALSE]
# Add gap where group changes
group_vector <- annotation_col$Survival
gap_index <- which(diff(as.numeric(group_vector)) != 0)

pheatmap(beta_subset_ordered,
         annotation_col = annotation_col,
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         gaps_col = gap_index,
         show_rownames = FALSE, show_colnames = FALSE,
         main = "Methylation of  CpGs")
#

#
top_annot <- HeatmapAnnotation(
  group = annotation_col$Survival,
  col = list(group = c("alive" = "#1b9e77", "dead" = "#d95f02")),
  annotation_name_side = "left",
  annotation_legend_param = list(title = "Survival"),
  show_annotation_name = TRUE
)
column_order <- order(annotation_col$Survival)

Heatmap(beta_subset,
        name = "Beta",
        top_annotation = top_annot,
        cluster_columns = FALSE,
        cluster_rows = TRUE,
        show_column_names = FALSE,
        show_row_names = FALSE,
        row_title = "survival", column_order = column_order,
        column_split = annotation_col$Survival,  # separates blocks
        column_title = NULL)
#

cpgs <-  GATA3_cpgs[GATA3_cpgs %in% row.names(beta_mat)]
cpgs <- as.vector(cpgs)
sig_cpgs <- row.names(DMPs)[DMPs$adj.P.Val < 0.05]  # or P.Value if unadjusted
cpgs <- intersect(GATA3_cpgs, sig_cpgs)
cpgs <- cpgs[cpgs %in% rownames(beta_mat)]  # ensure they exist in beta_mat
cpgs <- as.vector(cpgs)

for (cpg in cpgs) {
  group1 <- beta_mat[cpg, metadata$ten_year_survival == "dead"]
  group2 <- beta_mat[cpg, metadata$ten_year_survival == "alive"]
  p <- t.test(group1, group2)$p.value

  boxplot(beta_mat[cpg, ] ~ metadata$ten_year_survival,
          main = paste("CpG:", cpg, "p =", signif(p, 3)),
          ylab = "Beta value",
          xlab = "Survival Group",
          col = c("blue", "orange"))
}



#

cpg <- cpgs[2]
df <- data.frame(beta = beta_mat[row.names(beta_mat) %in% cpg , ], survival = metadata$ten_year_survival)
ggplot(df, aes(x = survival, y = beta, fill = survival)) +
  geom_boxplot() + theme_minimal() +
  labs(title = paste("Methylation of", cpg))
#

# Boxplot with p-value
p1 <- ggplot(df, aes(x = survival, y = beta, color = survival)) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.7) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.4, color = "black") +
  stat_compare_means(method = "t.test", label = "p.format") +
  labs(title = paste("Jitter Plot of Methylation of", cpg)) +
  theme_minimal()


# Violin plot with p-value
p2 <- ggplot(df, aes(x = survival, y = beta, fill = survival)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +  # Optional: box inside violin
  stat_compare_means(method = "t.test", label = "p.format") +
  theme_minimal() +
  labs(title = paste("Violin Plot of Methylation of", cpg))

# Combine both plots
grid.arrange(p1, p2, ncol = 2)


cpg <- cpgs[1]  # First CpG
df <- data.frame(beta = beta_mat[cpg, ], survival = metadata$ten_year_survival)

p <- ggplot(df, aes(x = survival, y = beta, fill = survival)) +
  geom_boxplot() + theme_minimal() +
  labs(title = paste("Methylation of", cpg))

ggplotly(p)



#
DMPs$logFC <- DMPs$logFC  # already calculated
DMPs$adj.P.Val <- DMPs$adj.P.Val  # already present

ggplot(DMPs, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), col = "red") +
  geom_vline(xintercept = c(-1, 1), col = "blue") +
  theme_minimal() +
  labs(title = "Volcano Plot of DMPs", x = "Log2 Fold Change", y = "-Log10 FDR")
#


# Create a summary data.frame
summary_df <- lapply(cpgs, function(cpg) {
  df <- data.frame(
    beta = beta_mat[cpg, ],
    survival = metadata$ten_year_survival
  )

  # Calculate group means
  group_means <- df %>%
    group_by(survival) %>%
    summarize(mean_beta = mean(beta), .groups = "drop")

  mean_alive <- group_means$mean_beta[group_means$survival == "alive"]
  mean_dead  <- group_means$mean_beta[group_means$survival == "dead"]

  # Wilcoxon test
  p_val <- wilcox.test(beta ~ survival, data = df)$p.value

  # Determine which group has higher methylation
  if (is.na(mean_alive) || is.na(mean_dead)) {
    higher_group <- NA
    diff <- NA
  } else {
    higher_group <- ifelse(mean_alive > mean_dead, "alive", "dead")
    diff <- abs(mean_alive - mean_dead)
  }

  data.frame(
    CpG = cpg,
    p_value = round(p_val, 4),
    mean_alive = round(mean_alive, 4),
    mean_dead = round(mean_dead, 4),
    higher_group = higher_group,
    diff = round(diff, 4)
  )
}) %>% bind_rows()

# View the summary
head(summary_df)

# Assume you already have:
# - beta_mat (matrix of beta values: rows = CpGs, columns = samples)
# - metadata with barcode and ten_year_survival

# Step 1: Subset metadata to match beta_mat columns
metadata <- metadata[match(colnames(beta_mat), metadata$barcode), ]

# Step 2: Calculate mean methylation per group
group_means <- t(apply(beta_mat, 1, function(x) {
  tapply(x, metadata$ten_year_survival, mean, na.rm = TRUE)
}))

# Step 3: Compute difference and which group is higher
diff_df <- as.data.frame(group_means)
diff_df$CpG <- rownames(diff_df)
diff_df$diff <- diff_df$alive - diff_df$dead
diff_df$higher_group <- ifelse(diff_df$diff > 0, "alive", "dead")

# Optional: Select top differentially methylated CpGs
diff_df_top <- diff_df[order(abs(diff_df$diff), decreasing = TRUE), ][1:20, ]

# Step 4: Visualize

ggplot(diff_df_top, aes(x = reorder(CpG, diff), y = diff, fill = higher_group)) +
  geom_col() +
  coord_flip() +
  labs(title = "Top Differentially Methylated CpGs",
       x = "CpG site",
       y = "Alive - Dead Methylation Difference") +
  scale_fill_manual(values = c("alive" = "#1f77b4", "dead" = "#d62728")) +
  theme_minimal()
# Filter statistically significant CpGs
sig_DMPs <- DMPs[DMPs$adj.P.Val < 0.05, ]
sig_DMPs <- DMPs[DMPs$adj.P.Val < 0.05, ]

# Identify direction of methylation difference
sig_DMPs$higher_group <- ifelse(sig_DMPs$logFC > 0, "alive", "dead")

# Count how many CpGs are more methylated in each group
ggplot(sig_DMPs, aes(x = higher_group, fill = higher_group)) +
  geom_bar() +
  labs(title = "Direction of Differential Methylation",
       x = "Group with Higher Methylation",
       y = "Number of Significant CpGs") +
  scale_fill_manual(values = c("alive" = "#1f77b4", "dead" = "#d62728")) +
  theme_minimal()
#



# Filter top significant CpGs (e.g., top 20 by adjusted p-value)
top_cpgs <- sig_DMPs %>%
  arrange(adj.P.Val) %>%
  slice(1:40) %>%
  pull(Name)

# Subset beta matrix
beta_subset <- beta_mat[top_cpgs, ]
beta_long <- melt(beta_subset)
colnames(beta_long) <- c("CpG", "Sample", "Beta")

# Add survival status
beta_long <- beta_long %>%
  left_join(metadata[, c("barcode", "ten_year_survival")], 
            by = c("Sample" = "barcode"))

# Convert survival to factor
beta_long$ten_year_survival <- factor(beta_long$ten_year_survival, levels = c("alive", "dead"))

# Plot using ggplot
ggplot(beta_long, aes(x = ten_year_survival, y = Beta, fill = ten_year_survival)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  facet_wrap(~ CpG, scales = "free_y", ncol = 4) +
  theme_minimal() +
  labs(title = "Violin Plot of Top 20 Differentially Methylated CpGs",
       x = "Survival Group",
       y = "Beta Value") +
  scale_fill_manual(values = c("alive" = "#1f77b4", "dead" = "#d62728")) +
  theme(strip.text = element_text(size = 8))
