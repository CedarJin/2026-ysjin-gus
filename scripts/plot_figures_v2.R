#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(forcats)
  library(scales)
  library(RColorBrewer)
  library(readxl)
})

args <- commandArgs(trailingOnly = TRUE)
arg_get <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (is.na(idx) || idx == length(args)) return(default)
  args[[idx + 1]]
}

report_dir <- arg_get("--report-dir", "results/summary/pident95_cov90/report")
tables_dir <- file.path(report_dir, "tables")
abundance_dir <- file.path(dirname(report_dir), "abundance_all_subject")
figures_dir <- file.path(report_dir, "figures_v2")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

standard_clade_order <- c("FMN", "L1", "L2", "NL", "NC", "hybrid", "mL1,2", "mL1", "mL2")

active_display <- c(
  WP_022501306 = "Eubacterium",
  `6MVG_A` = "R. gnavus",
  WP_012740861 = "E. eligens",
  `NP_416134.1` = "E. coli",
  `WP_003467686.1` = "C. perfringens",
  WP_000966715 = "S. agalactiae",
  WP_217761416 = "P. copri",
  WP_007841259 = "B. dorei",
  `6NCY_A` = "F. saccharivorans",
  WP_270534885 = "B. mass",
  `3CMG_A` = "B. fragilis"
)

active_colors <- c(
  "Eubacterium" = "#66C2A5",
  "R. gnavus" = "#1B9E77",
  "E. eligens" = "#FC8D62",
  "E. coli" = "#D95F02",
  "C. perfringens" = "#E6A05C",
  "S. agalactiae" = "#A63603",
  "P. copri" = "#7570B3",
  "B. dorei" = "#E7298A",
  "F. saccharivorans" = "#E6AB02",
  "B. mass" = "#A6CEE3",
  "B. fragilis" = "#1F78B4"
)

sra <- read_csv("reference_from_paper/SraRunTable.csv", show_col_types = FALSE)
run_to_sample <- setNames(
  sprintf("S%02d", as.integer(sra$`Sample Name`)),
  sra$Run
)

theme_v2 <- function(base_size = 16) {
  theme_classic(base_size = base_size) +
    theme(
      axis.title = element_text(face = "plain"),
      axis.text = element_text(color = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      legend.title = element_text(face = "plain"),
      legend.background = element_blank(),
      legend.key = element_blank(),
      plot.margin = margin(8, 12, 8, 8)
    )
}

save_v2 <- function(plot, stem, width, height, dpi = 300) {
  pdf_path <- file.path(figures_dir, paste0(stem, ".pdf"))
  png_path <- file.path(figures_dir, paste0(stem, ".png"))
  ggsave(pdf_path, plot, width = width, height = height, device = cairo_pdf, bg = "white")
  ggsave(png_path, plot, width = width, height = height, dpi = dpi, bg = "white")
}

# ---- Figure 3: clade abundance boxplot (patched total GUS = Pollet 279 + FMN 28) ----
total_meta <- read_tsv(
  "results/reference/gus_279_plus_fmn_metadata.tsv",
  show_col_types = FALSE
) %>%
  mutate(
    clade_std = recode(
      clade,
      mL1mL2 = "mL1,2",
      mL2Bacteroides = "mL2",
      mL2Parabacteroides = "mL2",
      .default = clade
    )
  )

total_within <- read_tsv(
  file.path(abundance_dir, "total_gus_abundance_within_gus_db_percent_all_subject.tsv"),
  show_col_types = FALSE
)

clade_long <- total_within %>%
  left_join(
    total_meta %>% transmute(sseqid = reference_id, clade = clade_std),
    by = "sseqid"
  ) %>%
  pivot_longer(-c(sseqid, clade), names_to = "subject", values_to = "percent") %>%
  group_by(subject, clade) %>%
  summarise(percent = sum(percent), .groups = "drop") %>%
  complete(subject, clade = standard_clade_order, fill = list(percent = 0)) %>%
  mutate(sample = recode(subject, !!!run_to_sample)) %>%
  filter(sample != "S02") %>%
  mutate(
    clade = factor(clade, levels = standard_clade_order),
    sample = factor(sample, levels = sort(unique(sample)))
  )

subject_levels <- levels(clade_long$sample)
subject_colors <- setNames(
  c(
    brewer.pal(8, "Set1"),
    brewer.pal(8, "Set2"),
    brewer.pal(8, "Dark2")
  )[seq_along(subject_levels)],
  subject_levels
)

figure3 <- ggplot(clade_long, aes(x = clade, y = percent)) +
  geom_point(
    aes(color = sample),
    position = position_jitter(width = 0.18, height = 0, seed = 1),
    size = 2.2,
    alpha = 0.9
  ) +
  geom_boxplot(
    fill = "transparent",
    color = "black",
    outlier.shape = NA,
    width = 0.62,
    linewidth = 0.5
  ) +
  scale_color_manual(values = subject_colors, name = "Subject") +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.05))) +
  labs(
    x = "Clade",
    y = "Abundance (% of GUS database)"
  ) +
  theme_v2() +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))

save_v2(figure3, "figure3_clade_abundance_boxplot", width = 11.2, height = 7.2)

if (!is.null(arg_get("--only")) && arg_get("--only") == "figure3") {
  message("Wrote figures to ", normalizePath(figures_dir))
  quit(save = "no")
}

# ---- Figure 3 variant: Pollet 279 + Simpson Loop_Class == FMN (n=28) ----
# Does not overwrite the original figure3 files.
simpson_summary <- read_xlsx(
  "reference_from_paper/Supplemental_Table_2.xlsx",
  sheet = "GUS_Summary"
)
fmn28_ids <- paste0(
  simpson_summary$Gene_ID[simpson_summary$Loop_Class == "FMN"],
  "-FMN"
)
pollet_meta <- read_tsv("results/reference/gus_279_metadata.tsv", show_col_types = FALSE) %>%
  mutate(
    clade_std = recode(
      clade,
      mL1mL2 = "mL1,2",
      mL2Bacteroides = "mL2",
      mL2Parabacteroides = "mL2",
      .default = clade
    )
  )
keep_ids <- unique(c(pollet_meta$reference_id, fmn28_ids))
clade_lookup <- bind_rows(
  pollet_meta %>% transmute(sseqid = reference_id, clade = clade_std),
  tibble(sseqid = fmn28_ids, clade = "FMN")
)

# ---- Figure 3 variant: Pollet 279 + Zhang candidate 18 + Simpson FMN 28 ----
candidate_meta <- read_tsv(
  "results/reference/candidate_gus_metadata.tsv",
  show_col_types = FALSE
) %>%
  mutate(
    clade_std = recode(
      clade,
      mL1mL2 = "mL1,2",
      mL2Bacteroides = "mL2",
      mL2Parabacteroides = "mL2",
      .default = clade
    )
  )
clade_lookup_plus <- bind_rows(
  clade_lookup,
  candidate_meta %>% transmute(sseqid = reference_id, clade = clade_std)
)
plus_clade_order <- c("FMN", "L1", "L2", "NL", "NC", "hybrid", "mL1,2", "mL1", "mL2")

total_genome <- read_tsv(
  file.path(abundance_dir, "total_gus_abundance_gus_genome_all_subject.tsv"),
  show_col_types = FALSE
)
candidate_genome <- read_tsv(
  file.path(abundance_dir, "candidate_gus_abundance_gus_genome_all_subject.tsv"),
  show_col_types = FALSE
)
missing_cand <- setdiff(candidate_meta$reference_id, candidate_genome$sseqid)
if (length(missing_cand) > 0) {
  zeros <- matrix(
    0,
    nrow = length(missing_cand),
    ncol = ncol(candidate_genome) - 1,
    dimnames = list(NULL, names(candidate_genome)[-1])
  )
  candidate_genome <- bind_rows(
    candidate_genome,
    bind_cols(tibble(sseqid = missing_cand), as.data.frame(zeros))
  )
}

combined_genome <- bind_rows(
  total_genome %>% filter(sseqid %in% keep_ids),
  candidate_genome
)
message(
  "Pollet279+candidate18+FMN28: ",
  n_distinct(combined_genome$sseqid), " abundance rows; ",
  nrow(candidate_meta), " candidate IDs"
)

clade_long_plus <- combined_genome %>%
  left_join(clade_lookup_plus, by = "sseqid") %>%
  pivot_longer(-c(sseqid, clade), names_to = "subject", values_to = "tpm") %>%
  group_by(subject, clade) %>%
  summarise(tpm = sum(tpm), .groups = "drop") %>%
  group_by(subject) %>%
  mutate(percent = tpm / sum(tpm) * 100) %>%
  ungroup() %>%
  complete(subject, clade = plus_clade_order, fill = list(tpm = 0, percent = 0)) %>%
  mutate(sample = recode(subject, !!!run_to_sample)) %>%
  filter(sample != "S02") %>%
  mutate(
    clade = factor(clade, levels = plus_clade_order),
    sample = factor(sample, levels = sort(unique(sample)))
  )

figure3_plus <- ggplot(clade_long_plus, aes(x = clade, y = percent)) +
  geom_point(
    aes(color = sample),
    position = position_jitter(width = 0.18, height = 0, seed = 1),
    size = 2.2,
    alpha = 0.9
  ) +
  geom_boxplot(
    fill = "transparent",
    color = "black",
    outlier.shape = NA,
    width = 0.62,
    linewidth = 0.5
  ) +
  scale_color_manual(values = subject_colors, name = "Subject") +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.05))) +
  labs(
    x = "Clade",
    y = "Abundance (% of GUS database)"
  ) +
  theme_v2() +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))

save_v2(
  figure3_plus,
  "figure3_clade_abundance_boxplot_pollet279_fmn28",
  width = 11.2,
  height = 7.2
)

# ---- Figure 4b: active GUS stacked composition ----
sample_order <- c(
  "S11", "S19", "S18", "S09", "S07", "S15", "S05", "S03",
  "S20", "S16", "S13", "S04", "S12", "S08", "S14", "S10", "S17", "S21"
)

active_wide <- read_tsv(
  file.path(abundance_dir, "active_gus_abundance_within_gus_db_percent_all_subject.tsv"),
  show_col_types = FALSE
)

missing_ids <- setdiff(names(active_display), active_wide$sseqid)
if (length(missing_ids) > 0) {
  zeros <- matrix(
    0,
    nrow = length(missing_ids),
    ncol = ncol(active_wide) - 1,
    dimnames = list(NULL, names(active_wide)[-1])
  )
  active_wide <- bind_rows(
    active_wide,
    bind_cols(tibble(sseqid = missing_ids), as.data.frame(zeros))
  )
}

active_long <- active_wide %>%
  mutate(protein = recode(sseqid, !!!active_display)) %>%
  filter(protein %in% names(active_colors)) %>%
  pivot_longer(-c(sseqid, protein), names_to = "subject", values_to = "percent") %>%
  mutate(sample = recode(subject, !!!run_to_sample)) %>%
  filter(sample != "S02") %>%
  mutate(
    sample = factor(sample, levels = sample_order[sample_order %in% unique(sample)]),
    protein = factor(protein, levels = names(active_colors))
  )

figure4b <- ggplot(active_long, aes(x = sample, y = percent, fill = protein)) +
  geom_col(width = 0.78, color = NA) +
  scale_fill_manual(values = active_colors, drop = FALSE, name = "Active GUS") +
  scale_y_continuous(
    breaks = c(0, 25, 50, 75, 100),
    expand = c(0, 0)
  ) +
  coord_cartesian(ylim = c(0, 100), clip = "off") +
  labs(
    x = "Subject",
    y = "Abundance (% of active GUS set)"
  ) +
  theme_v2() +
  theme(legend.position = "right") +
  guides(fill = guide_legend(title = "Active GUS"))

save_v2(figure4b, "figure4b_active_gus_composition", width = 12.4, height = 7.0)

message("Wrote figures to ", normalizePath(figures_dir))
