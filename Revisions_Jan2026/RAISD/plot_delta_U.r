library(ggplot2)
library(dplyr)
library(readr)
library(stringr)

# If df got defined as a function earlier, remove it
if (exists("df")) rm(df)

# Read data
df <- read_tsv(
  "cra_delta_mu_50kb.tsv",
  col_names = c("chr", "start", "end", "mu_pre", "mu_post", "delta_mu"),
  show_col_types = FALSE
)

# Parse chr ordering + window midpoint
df <- df %>%
  mutate(
    chr = as.character(chr),
    chr_num = case_when(
      str_detect(chr, "^\\d+$")    ~ as.numeric(chr),
      str_detect(chr, "^(\\d+)A$") ~ as.numeric(str_match(chr, "^(\\d+)A$")[,2]) + 0.5,
      chr %in% c("Z", "chrZ")      ~ 100,
      chr %in% c("W", "chrW")      ~ 101,
      chr %in% c("MT","chrM","M")  ~ 102,
      TRUE                         ~ 999
    ),
    mid = (start + end) / 2
  ) %>%
  arrange(chr_num, chr)

# Freeze chromosome order
chr_levels <- df %>%
  distinct(chr, chr_num) %>%
  arrange(chr_num, chr) %>%
  pull(chr)

df <- df %>%
  mutate(chr = factor(chr, levels = chr_levels))

# Compute cumulative x positions in numeric chr order
chr_sizes <- df %>%
  group_by(chr) %>%
  summarise(
    chr_len = max(end),
    chr_num = first(chr_num),
    .groups = "drop"
  ) %>%
  arrange(chr_num, chr) %>%
  mutate(cum_start = lag(cumsum(chr_len), default = 0))

df <- df %>%
  left_join(chr_sizes %>% select(chr, chr_len, cum_start), by = "chr") %>%
  mutate(cum_pos = mid + cum_start)

axis_df <- chr_sizes %>%
  mutate(center = cum_start + chr_len / 2)

# Manhattan plot
p <- ggplot(df, aes(x = cum_pos, y = delta_mu, color = chr)) +
  geom_point(size = 0.6, alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  scale_x_continuous(breaks = axis_df$center, labels = axis_df$chr) +
  scale_color_manual(values = rep(c("grey30", "grey65"), length.out = nlevels(df$chr))) +
  labs(
    x = "Chromosome",
    y = expression(Delta * mu),
    title = expression("Manhattan plot of " * Delta * mu * " (post \u2212 pre)")
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(p)

ggsave("cra_delta_mu_50kb.manhattan.png", p, width = 12, height = 4, dpi = 300)
