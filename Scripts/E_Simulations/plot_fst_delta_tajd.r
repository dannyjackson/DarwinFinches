#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(purrr)
})

# ------------------------------------------------------------
# Settings
# ------------------------------------------------------------
root_dir <- "output"
outfile  <- "fst_vs_delta_tajd_grid.png"

keep_decline <- NULL  # e.g., c(1.00, 0.95, 0.90, 0.85)
keep_sel     <- NULL  # e.g., c(0.00, 0.25, 0.50, 0.75, 1.00)

# ------------------------------------------------------------
# Find all summary_stats.tsv files
# Expected layout: output/<decline_rate>/selection_<sel>/summary_stats.tsv
# ------------------------------------------------------------
files <- list.files(
  path = root_dir,
  pattern = "^summary_stats\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)

if (length(files) == 0) stop("No summary_stats.tsv files found under: ", root_dir)

# ------------------------------------------------------------
# Parse decline_rate and selection from path
# ------------------------------------------------------------
parse_from_path <- function(path) {
  rel <- gsub(paste0("^", root_dir, "/"), "", path)

  decline_chr <- str_split(rel, "/", simplify = TRUE)[1]
  sel_chr <- str_match(rel, "/selection_([0-9]+\\.[0-9]+|[0-9]+)")[, 2]

  tibble(
    file = path,
    decline_rate = suppressWarnings(as.numeric(decline_chr)),
    sel_s        = suppressWarnings(as.numeric(sel_chr))
  )
}

meta <- map_dfr(files, parse_from_path) %>%
  filter(!is.na(decline_rate), !is.na(sel_s))

if (nrow(meta) == 0) {
  stop(
    "Found summary_stats.tsv files, but none matched expected path structure:\n",
    "  output/<decline_rate>/selection_<sel>/summary_stats.tsv"
  )
}

if (!is.null(keep_decline)) meta <- meta %>% filter(decline_rate %in% keep_decline)
if (!is.null(keep_sel))     meta <- meta %>% filter(sel_s %in% keep_sel)
if (nrow(meta) == 0) stop("After filtering, no files remain to plot.")

# ------------------------------------------------------------
# Reader with forced types (prevents rep type mismatch)
# ------------------------------------------------------------
# If your TSV has extra columns sometimes, they will be guessed via .default
forced_types <- cols(
  rep           = col_character(), # <-- key fix
  sel_s         = col_double(),
  decline_rate  = col_double(),
  model         = col_character(),
  seed          = col_character(),
  num_sites     = col_double(),
  num_mutations = col_double(),
  pi_t1         = col_double(),
  pi_t2         = col_double(),
  dxy           = col_double(),
  fst_hudson    = col_double(),
  tajd_t1       = col_double(),
  tajd_t2       = col_double(),
  delta_tajd    = col_double(),
  .default      = col_guess()
)

read_one <- function(file, decline_rate_path, sel_s_path) {
  df <- read_tsv(file, col_types = forced_types, show_col_types = FALSE)

  probs <- problems(df)
  if (nrow(probs) > 0) {
    message("\nParsing issues in: ", file)
    print(probs)
  }

  need <- c("fst_hudson", "delta_tajd")
  missing <- setdiff(need, names(df))
  if (length(missing) > 0) {
    stop("Missing required columns in ", file, ": ", paste(missing, collapse = ", "))
  }

  df %>%
    mutate(
      file = file,
      decline_rate = decline_rate_path,  # overwrite from directory name
      sel_s = sel_s_path                 # overwrite from directory name
    )
}

df_all <- pmap_dfr(meta %>% rename(decline_rate_path = decline_rate, sel_s_path = sel_s),
                   read_one)

# Order facets: selection left->right (ascending), decline top->bottom (descending)
sel_levels <- sort(unique(df_all$sel_s), decreasing = FALSE)
dec_levels <- sort(unique(df_all$decline_rate), decreasing = TRUE)

df_all <- df_all %>%
  mutate(
    sel_s = factor(sel_s, levels = sel_levels),
    decline_rate = factor(decline_rate, levels = dec_levels)
  )

# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------
p <- ggplot(df_all, aes(x = fst_hudson, y = delta_tajd)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
  geom_point(size = 0.25, alpha = 0.8, color = "blue") +
  facet_grid(
    rows = vars(decline_rate),
    cols = vars(sel_s),
    labeller = labeller(
      decline_rate = function(x) paste0("decline ", x),
      sel_s = function(x) paste0("s=", x)
    )
  ) +
  theme_classic(base_size = 14) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 11)
  ) +
  labs(
    x = "FST (Hudson)",
    y = expression(Delta * " Tajima's D"),
    title = "FST vs ΔTajima’s D across selection and population decline"
  )

ncol_facets <- length(sel_levels)
nrow_facets <- length(dec_levels)

ggsave(
  filename = outfile,
  plot = p,
  width = max(6, 3 * ncol_facets),
  height = max(5, 2.6 * nrow_facets),
  dpi = 300
)

message("\nRead ", nrow(meta), " file(s). Saved plot to: ", outfile)

# ------------------------------------------------------------
# Combined violin figure:
# rows = decline_rate
# cols = metric (FST left, ΔTajD right)
# within each panel: violins by selection (colored by selection)
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

# Long format with a metric column
df_long <- df_all %>%
  pivot_longer(
    cols = c(fst_hudson, delta_tajd),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = factor(
      metric,
      levels = c("fst_hudson", "delta_tajd"),
      labels = c("FST (Hudson)", "Δ Tajima's D")
    )
  )

p_violin_combo <- ggplot(df_long, aes(x = sel_s, y = value, fill = sel_s)) +
  geom_violin(scale = "width", trim = TRUE) +
  stat_summary(fun = median, geom = "crossbar", width = 0.6, fatten = 0) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  facet_wrap(
    vars(decline_rate, metric),
    ncol = 2,
    scales = "free_y",
    labeller = labeller(
      decline_rate = function(x) paste0("decline ", x)
    )
  ) +
  theme_classic(base_size = 14) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 11),
    legend.position = "none"
  ) +
  labs(
    x = "Selection coefficient (s)",
    y = NULL,
    title = "Distributions of FST and ΔTajima’s D by selection and population decline"
  )


# Save: width depends on #selection (violins), height depends on #decline (rows)
n_sel  <- length(levels(df_all$sel_s))
n_decl <- length(levels(df_all$decline_rate))

out_violin <- sub("\\.png$", ".violin_bySel_twoMetrics.png", outfile)

ggsave(
  filename = out_violin,
  plot = p_violin_combo,
  width = max(8, 1.2 * n_sel * 2),   # *2 because two metric columns
  height = max(5, 2.2 * n_decl),
  dpi = 300
)

message("Saved combined violin plot to: ", out_violin)


# ------------------------------------------------------------
# NEW: Combined scatter by decline_rate
# - For each decline_rate: overlay all sel_s values together
# - Color points by selection coefficient
# - Facet rows = decline_rate (stacked vertically)
# ------------------------------------------------------------

out_scatter_byDecline <- sub("\\.png$", ".scatter_byDecline_coloredBySel.png", outfile)

p_scatter_byDecline <- ggplot(df_all, aes(x = fst_hudson, y = delta_tajd)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +

  # 1) draw non-neutral first (bottom)
  geom_point(
    data = dplyr::filter(df_all, as.numeric(as.character(sel_s)) != 0),
    aes(color = sel_s),
    size = 0.25, alpha = 0.8
  ) +

  # 2) draw neutral last (top)
  geom_point(
    data = dplyr::filter(df_all, as.numeric(as.character(sel_s)) == 0),
    aes(color = sel_s),
    size = 0.25, alpha = 0.8
  ) +

  facet_grid(rows = vars(decline_rate), labeller = labeller(
    decline_rate = function(x) paste0("decline ", x)
  )) +
  theme_classic(base_size = 14) +
  theme(strip.background = element_blank(), strip.text = element_text(size = 11)) +
  labs(
    x = "FST (Hudson)",
    y = expression(Delta * " Tajima's D"),
    color = "Selection (s)",
    title = "FST vs ΔTajima’s D (all selections overlaid), faceted by population decline"
  )

ggsave(
  filename = out_scatter_byDecline,
  plot = p_scatter_byDecline,
  width = max(6, 7),                 # single column; keep reasonable width
  height = max(5, 2.6 * nrow_facets),# scales with # decline rows (you already computed nrow_facets)
  dpi = 300
)

message("Saved combined scatter (colored by selection) to: ", out_scatter_byDecline)
