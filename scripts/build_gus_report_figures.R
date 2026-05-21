#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(forcats)
  library(tibble)
  library(scales)
  library(patchwork)
  library(ggalluvial)
  library(svglite)
})

args <- commandArgs(trailingOnly = TRUE)

arg_get <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (is.na(idx) || idx == length(args)) return(default)
  args[[idx + 1]]
}

report_dir <- arg_get("--report-dir", "results/summary/pident95_cov90/report")
tables_dir <- file.path(report_dir, "tables")
figures_dir <- file.path(report_dir, "figures")
abundance_dir <- file.path(dirname(report_dir), "abundance_all_subject")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

standard_clade_order <- c("FMN", "L1", "L2", "NL", "NC", "hybrid", "m1,2", "mL1", "mL2")
clade_colors <- c(
  FMN = "#1b9e77",
  L1 = "#d95f02",
  L2 = "#7570b3",
  NL = "#e7298a",
  NC = "#66a61e",
  hybrid = "#e6ab02",
  `m1,2` = "#a6761d",
  mL1 = "#1f78b4",
  mL2 = "#b2df8a"
)

theme_report <- function() {
  theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 15),
      plot.subtitle = element_text(size = 11, color = "#52606d"),
      axis.title = element_text(face = "bold", size = 12),
      legend.title = element_text(face = "bold"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "right",
      plot.margin = margin(10, 12, 10, 12)
    )
}

save_svg <- function(plot, path, width, height) {
  ggsave(
    filename = path,
    plot = plot,
    width = width,
    height = height,
    device = svglite::svglite,
    bg = "white",
    limitsize = FALSE
  )
}

read_abundance_matrix <- function(path) {
  df <- read_tsv(path, show_col_types = FALSE)
  if (ncol(df) < 2) stop("Unexpected matrix format: ", path)
  rn <- df[[1]]
  df[[1]] <- NULL
  mat <- as.matrix(df)
  rownames(mat) <- rn
  storage.mode(mat) <- "double"
  mat
}

read_metadata <- function(path) {
  read_tsv(path, show_col_types = FALSE) %>%
    mutate(
      clade_std = recode(
        clade,
        mL1mL2 = "m1,2",
        mL2Bacteroides = "mL2",
        mL2Parabacteroides = "mL2",
        .default = clade
      ),
      taxon_label = if_else(
        !is.na(genus) & genus != "",
        genus,
        if_else(!is.na(phylum) & phylum != "", phylum, "unclear")
      )
    )
}

palette_subjects <- function(subjects) {
  setNames(viridisLite::turbo(length(subjects)), subjects)
}

collapse_taxonomy <- function(clade_tax) {
  tax_order <- clade_tax %>%
    group_by(taxon) %>%
    summarise(total_weight = sum(weight), .groups = "drop") %>%
    arrange(desc(total_weight)) %>%
    pull(taxon)

  keep <- tax_order[tax_order != "Other"]
  keep <- unique(c(keep[keep != "unclear"], "unclear"))
  keep <- head(keep, 10)
  keep <- unique(c(keep, "unclear"))

  clade_tax %>%
    mutate(
      taxon_plot = if_else(taxon %in% keep, taxon, "Other"),
      taxon_plot = factor(taxon_plot, levels = c(setdiff(keep, "unclear"), if ("unclear" %in% keep) "unclear" else NULL, "Other"))
    )
}

total_of_genome <- read_abundance_matrix(file.path(abundance_dir, "total_gus_abundance_of_genome_percent_all_subject.tsv"))
candidate_of_genome <- read_abundance_matrix(file.path(abundance_dir, "candidate_gus_abundance_of_genome_percent_all_subject.tsv"))
total_within <- read_abundance_matrix(file.path(abundance_dir, "total_gus_abundance_within_gus_db_percent_all_subject.tsv"))
candidate_within <- read_abundance_matrix(file.path(abundance_dir, "candidate_gus_abundance_within_gus_db_percent_all_subject.tsv"))

total_meta <- read_metadata(file.path("results/reference/gus_279_plus_fmn_metadata.taxonomy.tsv"))

subject_order <- colnames(total_of_genome)[order(colSums(total_of_genome), decreasing = TRUE)]
subject_colors <- palette_subjects(subject_order)

fmt_pct <- function(x, digits = 3) formatC(x, format = "f", digits = digits)
fmt_count <- function(x) formatC(x, format = "d", big.mark = ",")

subject_burden <- tibble(
  subject = subject_order,
  total_gus_of_genome_percent = colSums(total_of_genome[, subject_order, drop = FALSE]),
  candidate_gus_of_genome_percent = colSums(candidate_of_genome[, subject_order, drop = FALSE]),
  total_gus_detected_count = colSums(total_within[, subject_order, drop = FALSE] > 0),
  candidate_gus_detected_count = colSums(candidate_within[, subject_order, drop = FALSE] > 0)
) %>%
  mutate(subject = factor(subject, levels = subject_order))

prevalence_count_plot <- function(mat, title, fill) {
  prev <- rowSums(mat > 0)
  plot_df <- tibble(prevalence = prev) %>%
    count(prevalence, name = "n") %>%
    complete(prevalence = 1:18, fill = list(n = 0))

  ggplot(plot_df, aes(x = prevalence, y = n)) +
    geom_col(width = 0.82, fill = fill, color = "white", linewidth = 0.25) +
    scale_x_continuous(breaks = 1:18, limits = c(0.5, 18.5)) +
    scale_y_continuous(breaks = pretty_breaks()) +
    labs(
      title = title,
      x = "Detection prevalence (subjects out of 18)",
      y = "GUS genes (count)"
    ) +
    theme_report() +
    theme(plot.title = element_text(size = 13))
}

figure1_total <- ggplot(
  subject_burden %>%
    select(subject, total_gus_of_genome_percent, candidate_gus_of_genome_percent) %>%
    pivot_longer(-subject, names_to = "dataset", values_to = "percent") %>%
    mutate(dataset = recode(dataset,
      total_gus_of_genome_percent = "Total GUS",
      candidate_gus_of_genome_percent = "Candidate GUS"
    )),
  aes(x = subject, y = percent, fill = dataset)
) +
  geom_col(position = position_dodge(width = 0.72), width = 0.65, color = "white", linewidth = 0.2) +
  coord_flip() +
  scale_fill_manual(values = c("Total GUS" = "#2f6fdf", "Candidate GUS" = "#f28e2b")) +
  scale_y_continuous(labels = label_number(accuracy = 0.001)) +
  labs(
    title = "Subject-level GUS abundance",
    x = "Subject",
    y = "Genome percent (%)",
    fill = NULL
  ) +
  theme_report() +
  theme(legend.position = "top", panel.grid.major.y = element_blank())

figure1_counts <- ggplot(
  subject_burden %>%
    select(subject, total_gus_detected_count, candidate_gus_detected_count) %>%
    pivot_longer(-subject, names_to = "dataset", values_to = "count") %>%
    mutate(dataset = recode(dataset,
      total_gus_detected_count = "Total GUS",
      candidate_gus_detected_count = "Candidate GUS"
    )),
  aes(x = subject, y = count, fill = dataset)
) +
  geom_col(position = position_dodge(width = 0.72), width = 0.65, color = "white", linewidth = 0.2) +
  coord_flip() +
  scale_fill_manual(values = c("Total GUS" = "#2f6fdf", "Candidate GUS" = "#f28e2b")) +
  scale_y_continuous(breaks = pretty_breaks()) +
  labs(
    title = "Detected GUS count per subject",
    x = "Subject",
    y = "Detected genes (count)",
    fill = NULL
  ) +
  theme_report() +
  theme(legend.position = "top", panel.grid.major.y = element_blank())

figure1_summary <- sprintf(
  "Across 18 subjects, total GUS load ranges from %s%% to %s%% of the genome and candidate GUS ranges from %s%% to %s%%; detected counts range from %s to %s total GUS genes and %s to %s candidate genes.",
  fmt_pct(min(subject_burden$total_gus_of_genome_percent)),
  fmt_pct(max(subject_burden$total_gus_of_genome_percent)),
  fmt_pct(min(subject_burden$candidate_gus_of_genome_percent)),
  fmt_pct(max(subject_burden$candidate_gus_of_genome_percent)),
  fmt_count(min(subject_burden$total_gus_detected_count)),
  fmt_count(max(subject_burden$total_gus_detected_count)),
  fmt_count(min(subject_burden$candidate_gus_detected_count)),
  fmt_count(max(subject_burden$candidate_gus_detected_count))
)

figure1 <- figure1_total / figure1_counts + plot_annotation(tag_levels = "A")
save_svg(figure1, file.path(figures_dir, "figure1_subject_burden.svg"), 12.5, 11.0)

total_prev_plot <- prevalence_count_plot(total_of_genome, "Total GUS prevalence", "#2f6fdf")
candidate_prev_plot <- prevalence_count_plot(candidate_of_genome, "Candidate GUS prevalence", "#f28e2b")

total_core_genes <- sum(rowSums(total_of_genome > 0) >= 9)
candidate_core_genes <- sum(rowSums(candidate_of_genome > 0) >= 9)
figure2_summary <- sprintf(
  "Total GUS prevalence is sparse: %s of %s genes are core (detected in at least 9/18 subjects) and none are universal; the candidate set contains %s core genes.",
  fmt_count(total_core_genes),
  fmt_count(nrow(total_of_genome)),
  fmt_count(candidate_core_genes)
)

core_ids <- names(which(rowSums(total_of_genome > 0) >= 9))
core_meta <- total_meta %>%
  filter(reference_id %in% core_ids) %>%
  arrange(desc(match(reference_id, core_ids))) %>%
  mutate(reference_id = factor(reference_id, levels = core_ids))

core_mat <- total_of_genome[core_ids, subject_order, drop = FALSE]
heat_long <- as.data.frame(as.table(core_mat), stringsAsFactors = FALSE) %>%
  rename(reference_id = Var1, subject = Var2, tpm = Freq) %>%
  mutate(
    reference_id = factor(reference_id, levels = core_ids),
    subject = factor(subject, levels = subject_order),
    log_tpm = log10(tpm + 0.001)
  ) %>%
  left_join(core_meta %>% select(reference_id, clade_std), by = "reference_id")

strip_df <- core_meta %>%
  transmute(reference_id, clade_std = factor(clade_std, levels = standard_clade_order)) %>%
  arrange(match(reference_id, core_ids))

heatmap_plot <- ggplot(heat_long, aes(x = subject, y = reference_id, fill = log_tpm)) +
  geom_tile(color = "white", linewidth = 0.15) +
  scale_fill_gradientn(
    colours = viridisLite::mako(100),
    name = "Genome-wide TPM\nlog10(TPM + 0.001)"
  ) +
  scale_x_discrete(position = "top", expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(
    title = "Core total GUS heatmap",
    x = "Subject",
    y = "GUS gene"
  ) +
  theme_report() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 0, vjust = 0),
    axis.text.y = element_text(size = 7),
    legend.position = "right"
  )

clade_strip_plot <- ggplot(strip_df, aes(x = 1, y = reference_id, fill = clade_std)) +
  geom_tile(color = "white", linewidth = 0.15) +
  scale_fill_manual(values = clade_colors, drop = TRUE, name = "Clade") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  theme_void(base_size = 12) +
  theme(
    legend.position = "bottom",
    plot.margin = margin(6, 2, 6, 8)
  )

figure2_top <- total_prev_plot + candidate_prev_plot + plot_layout(ncol = 2, guides = "collect")
figure2_bottom <- clade_strip_plot + heatmap_plot + plot_layout(widths = c(0.08, 0.92), guides = "collect")
figure2 <- figure2_top / figure2_bottom + plot_layout(heights = c(1.0, 2.5)) + plot_annotation(tag_levels = "A")
save_svg(figure2, file.path(figures_dir, "figure2_prevalence_and_heatmap.svg"), 15.0, 12.5)

clade_subject_long <- as.data.frame(total_within[, subject_order, drop = FALSE]) %>%
  mutate(reference_id = rownames(total_within)) %>%
  left_join(total_meta %>% select(reference_id, clade_std), by = "reference_id") %>%
  pivot_longer(cols = all_of(subject_order), names_to = "subject", values_to = "percent") %>%
  mutate(
    subject = factor(subject, levels = subject_order),
    clade_std = factor(clade_std, levels = standard_clade_order)
  ) %>%
  group_by(subject, clade_std) %>%
  summarise(percent = sum(percent), .groups = "drop") %>%
  filter(percent > 0)

clade_point_plot <- ggplot(clade_subject_long, aes(x = clade_std, y = percent)) +
  geom_boxplot(fill = "#e9eef5", color = "#4c566a", outlier.shape = NA, width = 0.62) +
  geom_point(aes(color = subject), position = position_jitter(width = 0.16, height = 0), size = 2.3, alpha = 0.85) +
  scale_color_manual(values = subject_colors) +
  scale_y_continuous(labels = label_number(accuracy = 0.001)) +
  labs(
    title = "Subject-level clade abundance",
    x = "Clade",
    y = "Abundance (% of GUS database)",
    color = "Subject"
  ) +
  theme_report() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 35, hjust = 1)
  ) +
  guides(color = guide_legend(nrow = 3, byrow = TRUE, override.aes = list(size = 3, alpha = 1)))

save_svg(clade_point_plot, file.path(figures_dir, "figure3_clade_abundance_boxplot.svg"), 12.8, 8.2)

clade_stack_plot <- clade_subject_long %>%
  ggplot(aes(x = subject, y = percent, fill = clade_std)) +
  geom_col(width = 0.72, color = "white", linewidth = 0.2) +
  scale_fill_manual(values = clade_colors, drop = TRUE, name = "Clade") +
  scale_y_continuous(labels = label_number(accuracy = 0.001)) +
  labs(
    title = "Subject-by-subject clade composition",
    x = "Subject",
    y = "Abundance (% of GUS database)"
  ) +
  theme_report() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  ) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))

save_svg(clade_stack_plot, file.path(figures_dir, "figure4_clade_stacked_bars.svg"), 14.8, 7.8)

tax_res_plot <- read_tsv(file.path(tables_dir, "taxonomy_resolution_total.tsv"), show_col_types = FALSE) %>%
  mutate(category = factor(category, levels = c("genus", "phylum", "unclear"))) %>%
  ggplot(aes(x = category, y = count, fill = category)) +
  geom_col(width = 0.65, color = "white", linewidth = 0.2) +
  geom_text(aes(label = sprintf("%d (%.1f%%)", count, percent)), vjust = -0.35, size = 4) +
  scale_fill_manual(values = c(genus = "#2f6fdf", phylum = "#3fa34d", unclear = "#e15759")) +
  scale_y_continuous(breaks = pretty_breaks()) +
  labs(
    title = "Taxonomy resolution among detected total GUS",
    x = "Taxonomic resolution",
    y = "Detected GUS genes (count)",
    fill = NULL
  ) +
  theme_report() +
  theme(legend.position = "none")

save_svg(tax_res_plot, file.path(figures_dir, "figure5_taxonomy_resolution.svg"), 8.4, 6.2)

clade_tax <- read_tsv(file.path(tables_dir, "clade_taxonomy_abundance.tsv"), show_col_types = FALSE) %>%
  mutate(clade = factor(clade, levels = standard_clade_order)) %>%
  filter(weight > 0)

clade_tax_plot_df <- collapse_taxonomy(clade_tax)
clade_tax_agg <- clade_tax_plot_df %>%
  group_by(clade, taxon_plot) %>%
  summarise(pct = sum(pct), .groups = "drop")
tax_levels <- levels(clade_tax_plot_df$taxon_plot)
tax_colors <- setNames(viridisLite::turbo(length(tax_levels)), tax_levels)
if ("Other" %in% names(tax_colors)) tax_colors["Other"] <- "#B3B3B3"
if ("unclear" %in% names(tax_colors)) tax_colors["unclear"] <- "#666666"

clade_tax_plot <- clade_tax_agg %>%
  ggplot(aes(x = clade, y = pct, fill = taxon_plot)) +
  geom_col(width = 0.74, color = "white", linewidth = 0.2) +
  geom_text(
    data = clade_tax_agg %>% filter(pct >= 15),
    aes(label = taxon_plot),
    position = position_stack(vjust = 0.5),
    color = "white",
    size = 3.1,
    fontface = "bold",
    show.legend = FALSE
  ) +
  scale_fill_manual(values = tax_colors, drop = FALSE, name = "Taxonomy label") +
  scale_y_continuous(labels = label_number(accuracy = 0.1)) +
  coord_cartesian(ylim = c(0, 100)) +
  labs(
    title = "Abundance-weighted taxonomy composition by clade",
    x = "Clade",
    y = "Abundance within clade (%)",
    fill = "Taxonomy label"
  ) +
  theme_report() +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    legend.position = "bottom"
  ) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))

save_svg(clade_tax_plot, file.path(figures_dir, "figure6_clade_taxonomy_composition.svg"), 14.6, 8.6)

sankey_df <- read_tsv(file.path(tables_dir, "sankey_links.tsv"), show_col_types = FALSE) %>%
  mutate(
    left = factor(left, levels = standard_clade_order),
    right = factor(right, levels = tax_levels)
  ) %>%
  filter(weight > 0)

sankey_plot <- ggplot(sankey_df, aes(axis1 = left, axis2 = right, y = weight)) +
  geom_alluvium(aes(fill = left), alpha = 0.65, width = 0.15, knot.pos = 0.35, color = NA) +
  geom_stratum(width = 0.15, fill = "grey92", color = "grey50", linewidth = 0.25) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3, color = "#1f2933") +
  scale_x_discrete(limits = c("Clade", "Taxonomy label"), expand = c(0.12, 0.06)) +
  scale_fill_manual(values = clade_colors, drop = TRUE, name = "Clade") +
  scale_y_continuous(labels = label_number(accuracy = 0.1)) +
  labs(
    title = "Clade-to-taxonomy alluvial summary",
    x = "Feature layer",
    y = "Mean abundance (% of GUS database)",
    fill = "Clade"
  ) +
  theme_report() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(face = "bold")
  )

save_svg(sankey_plot, file.path(figures_dir, "figure7_sankey.svg"), 14.8, 8.4)
