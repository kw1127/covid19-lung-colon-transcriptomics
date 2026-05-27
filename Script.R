# =========================================================================================================
# Transcriptomic profile of fatal COVID-19 cases in lung and colon
# Wu et al. "Transcriptional and proteomic insights into the host response in fatal COVID-19 cases"
#
# Workflow:
#   1. Load and prepare counts/metadata
#   2. Pre-fit QC and model diagnostics
#   3. Pooled, tissue-adjusted differential expression (effect common to both organs)
#   4. Functional enrichment of the pooled signature (ORA + GSEA)
#   5. Tissue-resolved analysis (stratified DE, interaction test, per-tissue enrichment)
#   6. Assemble and export all figures as PNGs 
# =========================================================================================================

# ---- Libraries ----
library(dplyr)
library(tibble)
library(readr)
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(ashr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GO.db)
library(pheatmap)
library(RColorBrewer)
library(EnhancedVolcano)
library(enrichplot)
library(patchwork)
library(ggplotify)

# Shared theme: applied to every panel for a consistent look across figures
panel_theme <- theme(
  plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
  plot.subtitle = element_text(size = 9),
  plot.tag = element_text(size = 14, face = "bold"))

theme_set(theme_classic(base_size = 11))


# =========================================================================================================
# 1. Load and prepare data
# =========================================================================================================

raw_counts <- read_tsv("~/Practice_R/Data_bulk/lung_colon_covid19_counts.tsv")
raw_metadata <- read_tsv("~/Practice_R/Data_bulk/lung_colon_covid19_metadata.tsv")

# Counts matrix: gene IDs as row names, drop the gene-name column for now
counts <- raw_counts %>%
  dplyr::select(-"Gene Name") %>%
  column_to_rownames(var = "Gene ID")

# Metadata: keep disease status and tissue, with reference levels first
metadata <- raw_metadata %>%
  column_to_rownames(var = "Run") %>%
  dplyr::select(Disease = 'Sample Characteristic[disease]',
                Location = 'Sample Characteristic[organism part]')

metadata$Disease <- factor(metadata$Disease, levels = c("normal", "COVID-19"))
metadata$Location <- factor(metadata$Location, levels = c("lung", "colon"))

# Counts columns and metadata rows must align before building the dataset
stopifnot(all(colnames(counts) == rownames(metadata)))

# Pooled model: tissue as a covariate so the disease term estimates the
# disease effect COMMON to both organs (the shared signature)
dds <- DESeqDataSetFromMatrix(counts, metadata, design = ~ Location + Disease)

# Keep genes with at least 10 counts in at least 3 samples
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]


# =========================================================================================================
# 2. Pre-fit QC and model diagnostics
# =========================================================================================================

# Variance-stabilising transform for QC visualisation only (blind to design)
vsd <- vst(dds, blind = TRUE)

# PCA: do samples separate by disease status and tissue?
p_pca <- plotPCA(vsd, intgroup = c("Disease", "Location")) +
  labs(color = "Group") +
  ggtitle("PCA by disease status and tissue") +
  theme_classic() +
  panel_theme

# Sample-to-sample distances: a second view of overall structure
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$Disease, vsd$Location, sep = "-")
colnames(sampleDistMatrix) <- NULL

colours <- colorRampPalette(rev(brewer.pal(9, "YlOrRd")))(255)

ph_dist <- pheatmap(sampleDistMatrix,
                    clustering_distance_rows = sampleDists,
                    clustering_distance_cols = sampleDists,
                    col = colours,
                    main = "Sample-sample distance heatmap",
                    fontsize_row = 7,
                    show_colnames = FALSE,
                    silent = TRUE)
p_dist <- as.ggplot(ph_dist)

# Library size and size factors: check normalisation is well behaved
dds <- estimateSizeFactors(dds)
libsize <- colSums(counts(dds))

libdf <- tibble(sample = names(libsize), lib = libsize / 1e6) %>%
  arrange(lib) %>%
  mutate(sample = factor(sample, levels = sample))

p_libsize <- ggplot(libdf, aes(sample, lib)) +
  geom_col(fill = "grey70") +
  geom_hline(yintercept = median(libdf$lib), colour = "red", linetype = 2) +
  annotate("text", x = 1, y = median(libdf$lib),
           label = "median", hjust = 0, vjust = -0.5, size = 3, colour = "red") +
  labs(y = "Total counts (millions)", x = NULL, title = "Library size per sample") +
  theme_classic() + panel_theme +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 6))

sfdf <- tibble(sf = sizeFactors(dds), lib = libsize)
p_sizefactor <- ggplot(sfdf, aes(sf, lib)) +
  geom_point(shape = 1) +
  labs(x = "Size factor", y = "Library size", title = "Size factor vs library size") +
  theme_classic() + panel_theme

# Fit the model and check dispersion estimates
dds <- DESeq(dds)
p_disp <- as.ggplot(function() {
  par(mar = c(5, 5, 4, 2))
  plotDispEsts(dds, main = "Dispersion estimates", cex = 0.45, legend = FALSE)
  legend("bottomright", legend = c("gene-est", "fitted", "final"),
         col = c("black", "red", "dodgerblue"), pch = 16, cex = 0.7, bty = "n")})

normalised_counts <- counts(dds, normalized = TRUE)


# =========================================================================================================
# 3. Pooled differential expression: COVID-19 vs normal (tissue-adjusted)
# =========================================================================================================

res <- results(dds, contrast = c("Disease", "COVID-19", "normal"), alpha = 0.05)
summary(res)

# Shrink fold-changes for ranking and visualisation
res_shrunk <- lfcShrink(dds, contrast = c("Disease", "COVID-19", "normal"), type = "ashr")

p_ma <- as.ggplot(function() {
  plotMA(res_shrunk, ylim = c(-5, 5),
         main = "MA: COVID-19 vs normal (ashr-shrunk)")})

# Annotate with gene symbols, dropping genes with no padj or no symbol
gene_map <- raw_counts %>%
  dplyr::select(gene_id = "Gene ID", gene_name = "Gene Name") %>%
  distinct(gene_id, .keep_all = TRUE)

res_df <- as.data.frame(res_shrunk)
res_df$symbol <- gene_map$gene_name[match(rownames(res_df), gene_map$gene_id)]
res_df <- res_df[!is.na(res_df$padj), ]
res_df <- res_df[!is.na(res_df$symbol) & res_df$symbol != "", ]

symbols <- res_df$symbol[match(rownames(res_shrunk), rownames(res_df))]

# Significant DEGs, split by direction
sig <- subset(res_df, padj < 0.05 & abs(log2FoldChange) > 1)
sig_up <- subset(sig, log2FoldChange > 0)
sig_down <- subset(sig, log2FoldChange < 0)

# Volcano: label the top 10 genes in each direction
top_genes <- c(head(sig_up[order(sig_up$padj), "symbol"], 10),
               head(sig_down[order(sig_down$padj), "symbol"], 10))

p_volcano <- EnhancedVolcano(res_shrunk,
                             lab = symbols, selectLab = top_genes,
                             x = "log2FoldChange", y = "padj",
                             xlab = "log2 (fold change)",
                             ylab = "-log10 (p.adjust)",
                             pCutoff = 0.05, FCcutoff = 1,
                             labSize = 3, pointSize = 1.2,
                             drawConnectors = TRUE,
                             widthConnectors = 0.4,
                             colConnectors = "grey50",
                             max.overlaps = Inf,
                             min.segment.length = 0,
                             boxedLabels = FALSE,
                             legendPosition = "bottom",
                             title = "COVID-19 vs normal",
                             subtitle = "Top up- and down-regulated genes (tissue-adjusted)")

# Heatmap of the top 30 DEGs (15 per direction)
top_deg <- rbind(head(sig_up[order(sig_up$padj), ], 15),
                 head(sig_down[order(sig_down$padj), ], 15))

mat <- assay(vsd)[rownames(top_deg), ]
rownames(mat) <- top_deg$symbol

annot <- as.data.frame(colData(vsd)[, c("Disease", "Location")])
ann_colours <- list(
  Disease = c("normal" = "green", "COVID-19" = "orange"),
  Location = c("lung" = "blue", "colon" = "purple"))

ph_deg <- pheatmap(mat,
                   scale = "row",
                   annotation_col = annot,
                   annotation_colors = ann_colours,
                   show_colnames = FALSE,
                   clustering_distance_cols = "correlation",
                   color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(255),
                   main = "Top 30 DEGs (COVID-19 vs normal)",
                   silent = TRUE,
                   fontsize_row = 7)

p_deg <- as.ggplot(ph_deg)


# =========================================================================================================
# 4. Functional enrichment of the pooled signature
# =========================================================================================================

# ---- Over-representation analysis (GO BP), split by direction ----
universe_entrez <- bitr(res_df$symbol, "SYMBOL", "ENTREZID", org.Hs.eg.db)$ENTREZID
up_entrez <- bitr(sig_up$symbol, "SYMBOL", "ENTREZID", org.Hs.eg.db)$ENTREZID
down_entrez <- bitr(sig_down$symbol, "SYMBOL", "ENTREZID", org.Hs.eg.db)$ENTREZID

ego_up <- enrichGO(up_entrez, universe = universe_entrez,
                   OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
                   ont = "BP", pAdjustMethod = "BH",
                   qvalueCutoff = 0.05, readable = TRUE)
ego_down <- enrichGO(down_entrez, universe = universe_entrez,
                     OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
                     ont = "BP", pAdjustMethod = "BH",
                     qvalueCutoff = 0.05, readable = TRUE)

# Collapse redundant terms before plotting
ego_up <- simplify(ego_up, cutoff = 0.7, by = "p.adjust", select_fun = min)
ego_down <- simplify(ego_down, cutoff = 0.7, by = "p.adjust", select_fun = min)

p_ego_up <- dotplot(ego_up, showCategory = 15) +
  ggtitle("GO BP - upregulated in COVID-19") + panel_theme

p_ego_down <- dotplot(ego_down, showCategory = 15) +
  ggtitle("GO BP - downregulated in COVID-19") + panel_theme

# ---- GSEA (GO BP) on the full ranked gene list ----
# Uses the whole ranking, so it can recover coordinated but sub-threshold programmes
sym2entrez <- bitr(res_df$symbol, "SYMBOL", "ENTREZID", org.Hs.eg.db)
res_df$entrez <- sym2entrez$ENTREZID[match(res_df$symbol, sym2entrez$SYMBOL)]

ranks_df <- res_df[!is.na(res_df$entrez), ]
ranks_df <- ranks_df[!duplicated(ranks_df$entrez), ]

ranks <- ranks_df$log2FoldChange
names(ranks) <- ranks_df$entrez
ranks <- sort(ranks, decreasing = TRUE)

gse <- gseGO(ranks, OrgDb = org.Hs.eg.db, ont = "BP",
             keyType = "ENTREZID", pvalueCutoff = 0.05)

p_gse <- dotplot(gse, showCategory = 8, split = ".sign", label_format = 30) +
  facet_grid(. ~ .sign) +
  ggtitle("GSEA: GO BP (COVID-19 vs normal)") +
  panel_theme

p_ridge <- ridgeplot(gse, showCategory = 10) +
  labs(x = expression(log[2]~fold~change)) +
  ggtitle("GSEA fold-change distributions") +
  panel_theme

# Running-enrichment plot for the top gene set (collagen fibril organisation)
fig_gsea_single <- as.ggplot(
  gseaplot2(gse, geneSetID = 1, title = "", color = "firebrick", base_size = 12)) +
  plot_annotation(title = "Figure 6. GSEA enrichment: collagen fibril organization")


# =========================================================================================================
# 5. Tissue-resolved analysis: COVID-19 effects in lung and colon separately
# =========================================================================================================

# ---------------------------------------------------------------------------------------------------------
# 5a. Stratified per-tissue differential expression
#   The pooled model gives the effect COMMON to both organs; here we estimate the
#   effect WITHIN each organ, re-estimating dispersions and size factors per tissue.
# ---------------------------------------------------------------------------------------------------------

# Sanity check: stratified DE needs adequate replicates per group
print(table(metadata$Location, metadata$Disease))

run_tissue_de <- function(tissue) {
  samp <- rownames(metadata)[metadata$Location == tissue]
  md <- metadata[samp, , drop = FALSE]
  md$Disease <- droplevels(md$Disease)
  
  dds_t <- DESeqDataSetFromMatrix(counts[, samp], md, design = ~ Disease)
  dds_t <- dds_t[rowSums(counts(dds_t) >= 10) >= 3, ]   # re-filter within tissue
  dds_t <- DESeq(dds_t)
  
  out <- as.data.frame(
    lfcShrink(dds_t, contrast = c("Disease", "COVID-19", "normal"), type = "ashr"))
  out$gene_id <- rownames(out)
  out$symbol <- gene_map$gene_name[match(out$gene_id, gene_map$gene_id)]
  out
}

res_lung <- run_tissue_de("lung")
res_colon <- run_tissue_de("colon")

sig_lung <- subset(res_lung, !is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 1)
sig_colon <- subset(res_colon, !is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 1)

# Per-tissue volcanoes (same styling as the pooled volcano)
make_volcano <- function(res_t, title) {
  up <- subset(res_t, !is.na(padj) & padj < 0.05 & log2FoldChange > 1)
  dn <- subset(res_t, !is.na(padj) & padj < 0.05 & log2FoldChange < -1)
  lab <- c(head(up[order(up$padj), "symbol"], 8),
           head(dn[order(dn$padj), "symbol"], 8))
  EnhancedVolcano(res_t, lab = res_t$symbol, selectLab = lab,
                  x = "log2FoldChange", y = "padj",
                  xlab = "log2 (fold change)",
                  ylab = "-log10 (p.adjust)",
                  pCutoff = 0.05, FCcutoff = 1,
                  labSize = 3, pointSize = 1.0,
                  drawConnectors = TRUE, widthConnectors = 0.4,
                  colConnectors = "grey50", max.overlaps = Inf,
                  min.segment.length = 0, legendPosition = "bottom",
                  title = title, subtitle = NULL)
}

p_volcano_lung <- make_volcano(res_lung, "COVID-19 vs normal - lung")
p_volcano_colon <- make_volcano(res_colon, "COVID-19 vs normal - colon")

# ---------------------------------------------------------------------------------------------------------
# 5b. Formal tissue-dependence test + four-way gene classification
# ---------------------------------------------------------------------------------------------------------

dds_int <- DESeqDataSetFromMatrix(counts, metadata, design = ~ Location * Disease)
dds_int <- dds_int[rowSums(counts(dds_int) >= 10) >= 3, ]
dds_int <- DESeq(dds_int)

# Interaction coefficient name is version/level dependent, so fetch it rather than hard-code
int_name <- grep("^Location.*Disease", resultsNames(dds_int), value = TRUE)
message("Interaction coefficient: ", int_name)   

res_int_df <- as.data.frame(results(dds_int, name = int_name, alpha = 0.05))
res_int_df$gene_id <- rownames(res_int_df)

# Merge per-tissue fold-changes and classify each gene
cmp <- merge(
  res_lung[, c("gene_id", "symbol", "log2FoldChange", "padj")],
  res_colon[, c("gene_id", "log2FoldChange", "padj")],
  by = "gene_id", suffixes = c("_lung", "_colon"),
  all = TRUE)  

cmp <- cmp[!(is.na(cmp$padj_lung) & is.na(cmp$padj_colon)), ]

cmp$padj_int <- res_int_df$padj[match(cmp$gene_id, res_int_df$gene_id)]
cmp$min_padj <- pmin(cmp$padj_lung, cmp$padj_colon, na.rm = TRUE)

lfc_cut <- 1
p_cut <- 0.05
cmp <- cmp %>%
  mutate(
    sig_lung = !is.na(padj_lung) & padj_lung < p_cut & abs(log2FoldChange_lung) > lfc_cut,
    sig_colon = !is.na(padj_colon) & padj_colon < p_cut & abs(log2FoldChange_colon) > lfc_cut,
    int_sig = !is.na(padj_int) & padj_int < p_cut,
    class = case_when(
      sig_lung & sig_colon &
        sign(log2FoldChange_lung) != sign(log2FoldChange_colon) ~ "Discordant",
      sig_lung & sig_colon ~ "Shared (concordant)",
      sig_lung & !sig_colon ~ "Lung-specific",
      sig_colon & !sig_lung ~ "Colon-specific",
      TRUE ~ "Not significant"))

print(table(cmp$class))
# Of the tissue-specific calls, how many are backed by a significant interaction term?
print(with(subset(cmp, class %in% c("Lung-specific", "Colon-specific")), table(class, int_sig)))

# Lung-vs-colon scatter: concordant genes track the diagonal, tissue-specific genes sit on an axis
class_cols <- c("Shared (concordant)" = "#1b9e77", "Lung-specific" = "#377eb8",
                "Colon-specific" = "#984ea3", "Discordant" = "#e41a1c",
                "Not significant" = "grey80")

# Label the most significant genes per class, plus a few fibrosis markers of interest
fib_genes <- c("SPP1", "PLOD2", "SERPINE2", "COL1A1", "COL1A2", "COL3A1")
lab_df <- cmp %>%
  filter(class != "Not significant") %>%
  group_by(class) %>%
  slice_min(min_padj, n = 6, with_ties = FALSE) %>%
  ungroup() %>%
  bind_rows(filter(cmp, symbol %in% fib_genes)) %>%
  distinct(gene_id, .keep_all = TRUE)

p_class <- ggplot(cmp, aes(log2FoldChange_lung, log2FoldChange_colon)) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey60") +
  geom_point(aes(colour = class), size = 1, alpha = 0.7) +
  ggrepel::geom_text_repel(data = lab_df, aes(label = symbol),
                           size = 2.8, max.overlaps = Inf, min.segment.length = 0,
                           segment.colour = "grey60") +
  scale_colour_manual(values = class_cols, name = NULL) +
  labs(x = expression(log[2]~FC~~lung), y = expression(log[2]~FC~~colon),
       title = "COVID-19 effect: lung vs colon",
       subtitle = "Concordant genes track the diagonal; tissue-specific genes sit along an axis") +
  theme_classic() + panel_theme + theme(legend.position = "right")

# ---------------------------------------------------------------------------------------------------------
# 5c. Per-tissue enrichment: where does each programme live?
# ---------------------------------------------------------------------------------------------------------

to_entrez <- function(s) na.omit(bitr(s, "SYMBOL", "ENTREZID", org.Hs.eg.db)$ENTREZID)

# GO BP across all four tissue x direction groups in one comparison
gene_groups <- list(
  Lung_up = to_entrez(subset(sig_lung, log2FoldChange > 0)$symbol),
  Lung_down = to_entrez(subset(sig_lung, log2FoldChange < 0)$symbol),
  Colon_up = to_entrez(subset(sig_colon, log2FoldChange > 0)$symbol),
  Colon_down = to_entrez(subset(sig_colon, log2FoldChange < 0)$symbol))

cc <- compareCluster(gene_groups, fun = "enrichGO",
                     universe = universe_entrez, OrgDb = org.Hs.eg.db,
                     keyType = "ENTREZID", ont = "BP",
                     pAdjustMethod = "BH", qvalueCutoff = 0.05)
cc <- simplify(cc, cutoff = 0.7, by = "p.adjust", select_fun = min)

# 2. check the group sizes INSIDE cc — this is what the plot will print
as.data.frame(cc) |> dplyr::count(Cluster)

# 3. rebuild the plot object
p_cc <- dotplot(cc, showCategory = 4, label_format = 30) +
  ggtitle("Figure 9. GO BP enrichment by tissue and direction") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        axis.text.y = element_text(size = 10)) +
  panel_theme

# Enrichment-map networks of the up-regulated programme per tissue.
# A connected web indicates a coherent programme; isolated/absent nodes a diffuse signal.
run_ego <- function(sig_t, direction = c("up", "down")) {
  direction <- match.arg(direction)
  genes <- if (direction == "up")
    subset(sig_t, log2FoldChange > 0)$symbol
  else
    subset(sig_t, log2FoldChange < 0)$symbol
  
  ent <- na.omit(bitr(genes, "SYMBOL", "ENTREZID", org.Hs.eg.db)$ENTREZID)
  if (length(ent) < 5) return(NULL)             # too few genes to test
  
  ego <- enrichGO(ent, universe = universe_entrez, OrgDb = org.Hs.eg.db,
                  keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH",
                  qvalueCutoff = 0.05, readable = TRUE)
  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) return(NULL)
  simplify(ego, cutoff = 0.7, by = "p.adjust", select_fun = min)
}

ego_lung <- run_ego(sig_lung, "up")
ego_colon <- run_ego(sig_colon, "up")

# Report how many terms each tissue returned before plotting
cat("Lung up-regulated enriched terms:",
    if (is.null(ego_lung)) 0 else nrow(as.data.frame(ego_lung)), "\n")
cat("Colon up-regulated enriched terms:",
    if (is.null(ego_colon)) 0 else nrow(as.data.frame(ego_colon)), "\n")

make_emap <- function(ego, title) {
  n <- if (is.null(ego)) 0 else nrow(as.data.frame(ego))
  if (n < 2) {   # emapplot needs at least two terms; sparsity is itself a result
    return(as.ggplot(grid::textGrob(
      sprintf("%s\n(%d enriched term%s - no network)", title, n, ifelse(n == 1, "", "s")),
      gp = grid::gpar(fontsize = 10))))
  }
  ego@result$Description <- stringr::str_wrap(ego@result$Description, width = 18)
  old <- ggplot2::update_geom_defaults("text", list(size = 2.3))   # shrink node labels
  on.exit(ggplot2::update_geom_defaults("text", old), add = TRUE)  # restore afterwards
  emapplot(pairwise_termsim(ego), showCategory = 8) +
    ggtitle(title) + panel_theme
}

p_emap_lung <- make_emap(ego_lung, "Lung - up-regulated GO BP")
p_emap_colon <- make_emap(ego_colon, "Colon - up-regulated GO BP")


# =========================================================================================================
# 6. Assemble and export figures as PNGs
# =========================================================================================================

# Figure 1: PCA + sample-distance heatmap
fig_qc1 <- (p_pca | p_dist) +
  plot_layout(widths = c(1, 2)) +
  plot_annotation(
    tag_levels = "A",
    title = "Figure 1. Quality control and model diagnostics") &
  theme(plot.tag = element_text(size = 12, face = "bold"))

ggsave("fig_1_qc.png", fig_qc1, width = 15, height = 7, dpi = 300, bg = "white")

# Figure 2: library size, size factors, dispersion, MA
fig_qc2 <- wrap_plots(p_libsize, p_sizefactor, p_disp, p_ma, ncol = 2) +
  plot_annotation(
    tag_levels = "A",
    title = "Figure 2. Further quality control") &
  theme(plot.tag = element_text(size = 14, face = "bold"))

ggsave("fig_2_qc.png", fig_qc2, width = 12, height = 10, dpi = 300, bg = "white")

# Figure 3: pooled volcano + top-DEG heatmap
fig_de <- (p_volcano | p_deg) +
  plot_layout(widths = c(1, 1.3)) +
  plot_annotation(
    tag_levels = "A",
    title = "Figure 3. Differential expression: COVID-19 vs normal") &
  theme(plot.tag = element_text(size = 14, face = "bold"))

ggsave("fig_3_de.png", fig_de, width = 14, height = 8, dpi = 300, bg = "white")

# Figure 4: GO BP over-representation by direction
fig_ora <- (p_ego_up / p_ego_down) +
  plot_layout(guides = "collect", heights = c(15, 3)) +
  plot_annotation(
    tag_levels = "A",
    title = "Figure 4. GO BP over-representation by direction") &
  theme(plot.tag = element_text(size = 14, face = "bold"),
        legend.position = "right")

ggsave("fig_4_ora.png", fig_ora, width = 10, height = 12, dpi = 300, bg = "white")

# Figure 5: GSEA dot plot + ridge plot
fig_gsea <- (p_gse | p_ridge) +
  plot_layout(widths = c(1, 1.1)) +
  plot_annotation(
    tag_levels = "A",
    title = "Figure 5. GSEA of GO BP terms") &
  theme(plot.tag = element_text(size = 14, face = "bold"))

ggsave("fig_5_gsea.png", fig_gsea, width = 15, height = 9, dpi = 300, bg = "white")

# Figure 6: GSEA running-enrichment for the top gene set
ggsave("fig_6_gsea_single.png", fig_gsea_single, width = 15, height = 11, dpi = 300, bg = "white")

# Figure 7: per-tissue volcanoes
fig_tissue_de <- (p_volcano_lung | p_volcano_colon) +
  plot_layout(guides = "collect") +
  plot_annotation(
    tag_levels = "A",
    title = "Figure 7. Differential expression within each tissue") &
  theme(plot.tag = element_text(size = 14, face = "bold"), legend.position = "bottom")

ggsave("fig_7_tissue_de.png", fig_tissue_de, width = 14, height = 8, dpi = 300, bg = "white")

# Figure 8: lung-vs-colon fold-change scatter
ggsave("fig_8_lung_vs_colon.png", p_class, width = 9, height = 7, dpi = 300, bg = "white")

# Figure 9: GO BP enrichment across the four tissue x direction groups
ggsave("fig_9_tissue_ora.png", p_cc, width = 11, height = 13, dpi = 300, bg = "white")

# Figure 10: enrichment-map networks per tissue
fig_emap <- (p_emap_lung | p_emap_colon) +
  plot_layout(guides = "collect") +
  plot_annotation(
    tag_levels = "A",
    title = "Figure 10. Enrichment-map networks of up-regulated GO BP terms by tissue") &
  theme(plot.tag = element_text(size = 14, face = "bold"), legend.position = "right")

ggsave("fig_10_emap.png", fig_emap, width = 14, height = 7, dpi = 300, bg = "white")