suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(stringr)
})

# -----------------------------
# Settings
# -----------------------------
pops <- c("cra_pre", "cra_post")
metric_col <- "neg_log_pvalues_one_tailed"

# If your chromo column is like "NC_044572.1" instead of "1", this still works,
# but axis labels will be those strings.
# -----------------------------
# Helpers
# -----------------------------

read_raisD <- function(pop) {
  f <- paste0(pop, ".raisd.50.Ztransformed.csv")

  # Whitespace-delimited file
  df <- readr::read_table2(f)

  metric_col <- "neg_log_pvalues_one_tailed"

  # sanity check with a helpful error
  missing <- setdiff(c("chromo", "position", metric_col), names(df))
  if (length(missing) > 0) {
    stop(
      "Missing expected columns in ", f, ": ",
      paste(missing, collapse = ", "),
      "\nColumns found: ", paste(names(df), collapse = ", ")
    )
  }

  df %>%
    dplyr::transmute(
      pop = pop,
      chromo = as.character(chromo),
      position = as.numeric(position),
      value = as.numeric(.data[[metric_col]])
    ) %>%
    dplyr::filter(!is.na(position), !is.na(value))
}



make_cumpos <- function(df) {
  # Order chromosomes in a sensible way:
  # - if they’re purely numeric, order numerically
  # - otherwise order lexicographically
  chrom_levels <- df %>%
    distinct(chromo) %>%
    mutate(
      chrom_is_num = str_detect(chromo, "^[0-9]+$"),
      chrom_num = ifelse(chrom_is_num, as.numeric(chromo), NA_real_)
    ) %>%
    arrange(desc(chrom_is_num), chrom_num, chromo) %>%
    pull(chromo)

  df2 <- df %>%
    mutate(chromo = factor(chromo, levels = chrom_levels)) %>%
    group_by(chromo) %>%
    summarise(chr_len = max(position, na.rm = TRUE), .groups = "drop") %>%
    arrange(chromo) %>%
    mutate(chr_start = lag(cumsum(chr_len), default = 0)) %>%
    select(chromo, chr_len, chr_start)

  df %>%
    mutate(chromo = factor(as.character(chromo), levels = levels(df2$chromo))) %>%
    left_join(df2, by = "chromo") %>%
    mutate(cum_pos = position + chr_start)
}

# -----------------------------
# Load + prep
# -----------------------------
df_pre  <- read_raisD("cra_pre")
df_post <- read_raisD("cra_post")

df_all <- bind_rows(df_pre, df_post)

# Build a single coordinate system so pre/post align perfectly
df_all_cum <- make_cumpos(df_all)

# Axis centers for chromosome labels
axis_df <- df_all_cum %>%
  distinct(chromo, chr_start, chr_len) %>%
  mutate(center = chr_start + chr_len / 2) %>%
  arrange(chromo)

# Facet data:
#   - pre only
#   - post only
#   - overlay contains BOTH, with an extra grouping column "overlay_pop"
df_facets <- bind_rows(
  df_all_cum %>%
    filter(pop == "cra_pre") %>%
    mutate(panel = "cra_pre", overlay_pop = "cra_pre"),
  df_all_cum %>%
    filter(pop == "cra_post") %>%
    mutate(panel = "cra_post", overlay_pop = "cra_post"),
  df_all_cum %>%
    mutate(panel = "cra_pre_over_cra_post", overlay_pop = pop)
) %>%
  mutate(
    panel = factor(panel, levels = c("cra_pre", "cra_post", "cra_pre_over_cra_post"))
  )

df_facets <- df_facets %>%
  mutate(
    overlay_pop = factor(overlay_pop, levels = c("cra_post", "cra_pre"))
  )

# Alternating chromosome shading (optional but helps readability)
shade_df <- axis_df %>%
  mutate(idx = row_number()) %>%
  filter(idx %% 2 == 0) %>%
  transmute(
    xmin = chr_start,
    xmax = chr_start + chr_len,
    ymin = -Inf,
    ymax = Inf
  )

# -----------------------------
# Plot
# -----------------------------
p <- ggplot(df_facets, aes(x = cum_pos, y = value)) +
  # alternating chr shading
  geom_rect(
    data = shade_df,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    alpha = 0.08
  ) +
  # points; in overlay panel, color by which pop
  geom_point(
    aes(color = overlay_pop),
    size = 0.5,
    alpha = 0.7
  ) +
  facet_wrap(~ panel, ncol = 1, scales = "free_y") +
  scale_x_continuous(
    breaks = axis_df$center,
    labels = as.character(axis_df$chromo),
    expand = c(0.005, 0.005)
  ) +
  labs(
    x = "Chromosome",
    y = metric_col,
    color = NULL,
    title = paste0("Manhattan: ", metric_col, " (RAiSD windows)"),
    subtitle = "Top: cra_pre | Middle: cra_post | Bottom: cra_pre overlaid on cra_post"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 0, vjust = 0.5),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )

ggsave("manhattan_neg_log_pvalues_one_tailed_3facets.png", p, width = 14, height = 10, dpi = 300)
ggsave("manhattan_neg_log_pvalues_one_tailed_3facets.pdf", p, width = 14, height = 10)
