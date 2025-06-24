Make manhattan plots of TENM3 genes

### TENM3
#### NC_044574.1:32971132-34255509

# TENM3.fst.tsv

# CRA
# Load required packages, installing if necessary
required_packages <- c("qqman", "hexbin", "readr", "ggrepel", "ggplot2", "dplyr", "RColorBrewer", "data.table")
installed_packages <- rownames(installed.packages())

cat("Checking required packages...\n")
for (pkg in required_packages) {
  if (!(pkg %in% installed_packages)) {
    install.packages(pkg, repos = "http://cran.us.r-project.org")
  }
  library(pkg, character.only = TRUE)
}

color1 <- "#61ABAC"
color2 <- "#61ABAC"
input <- "TENM3.fst.tsv"
metric <- "fst"
pop_name<- "cra"

# Define parameters
cat("Reading in file...\n")
# Read file
data <- fread(input, sep = "\t", na.strings = c("", "NA"), data.table = TRUE)

# Prepare data for plotting
cat("Preparing data for plotting...\n")
data$chr <- factor(data$chr, levels = c(1, "1A", 2:4, "4A", 5:29, "Z"))

plot_data <- data %>%
  group_by(chr) %>%
  summarise(chr_len = max(midPos)) %>%
  mutate(tot = cumsum(chr_len) - chr_len) %>%
  select(-chr_len) %>%
  left_join(data, by = "chr") %>%
  arrange(chr, midPos) %>%
  mutate(BPcum = midPos + tot)

axisdf <- plot_data %>%
  group_by(chr) %>%
  summarize(center = mean(BPcum))


# Plot
cat("Generating plot...\n")
(294682,294682)(33849547,33849547)(33849547,33849548)   4       33849547        2       0.656030 # CpG, between exon 4,5

plot_data$fst[plot_data$fst < 0] <- 0


plot_data$color_flag <- ifelse(round(plot_data$fst, 6) == 0.656030, "highlight", as.character(plot_data$chr))
highlight_color <- "red"
chr_colors <- rep(c(color1, color2), length.out = length(unique(plot_data$chr)))
names(chr_colors) <- unique(plot_data$chr)


# Add the highlight color
color_values <- c(chr_colors, highlight = highlight_color)


ggplot(plot_data, aes(x = BPcum, y = !!sym(metric))) +
  # Background points (non-highlighted)
  geom_point(
    data = subset(plot_data, color_flag != "highlight"),
    aes(color = color_flag),
    alpha = 0.2, size = 8
  ) +
  # Foreground points (highlighted)
  geom_point(
    data = subset(plot_data, color_flag == "highlight"),
    color = "red",  # explicitly set highlight color
    alpha = 1, size = 8
  ) +
  scale_color_manual(values = color_values) +
  scale_x_continuous(
    breaks = range(plot_data$BPcum, na.rm = TRUE)
  ) +
  scale_y_continuous(
    breaks = range(plot_data[[metric]], na.rm = TRUE),
    expand = c(0, 0)
  ) +
  labs(x = NULL, y = NULL) +  # <--- this line removes axis labels
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )


ggsave(filename = file.path("TENM3.fst.cra.png"), 
       width = 10, height = 5, units = "in")

cat("Script completed successfully!\n")


# FOR

# Load required packages, installing if necessary
required_packages <- c("qqman", "hexbin", "readr", "ggrepel", "ggplot2", "dplyr", "RColorBrewer", "data.table")
installed_packages <- rownames(installed.packages())

cat("Checking required packages...\n")
for (pkg in required_packages) {
  if (!(pkg %in% installed_packages)) {
    install.packages(pkg, repos = "http://cran.us.r-project.org")
  }
  library(pkg, character.only = TRUE)
}

color1 <- "#E68580"
color2 <- "#E68580"
input <- "TENM3.fst.tsv"
metric <- "fst"

# Define parameters
cat("Reading in file...\n")
# Read file
data <- fread(input, sep = "\t", na.strings = c("", "NA"), data.table = TRUE)

# Prepare data for plotting
cat("Preparing data for plotting...\n")
data$chr <- factor(data$chr, levels = c(1, "1A", 2:4, "4A", 5:29, "Z"))

plot_data <- data %>%
  group_by(chr) %>%
  summarise(chr_len = max(midPos)) %>%
  mutate(tot = cumsum(chr_len) - chr_len) %>%
  select(-chr_len) %>%
  left_join(data, by = "chr") %>%
  arrange(chr, midPos) %>%
  mutate(BPcum = midPos + tot)

axisdf <- plot_data %>%
  group_by(chr) %>%
  summarize(center = mean(BPcum))


# Plot
cat("Generating plot...\n")



(289895,289895)(33149057,33149057)(33149057,33149058)   4       33149057        2       0.347509 # CpG, between exon 12 and 13
(296736,296736)(34069815,34069815)(34069815,34069816)   4       34069815        2       0.412365 # CpG site, between exon 3 and 4
(293493,293493)(33652266,33652266)(33652266,33652267)   4       33652266        2       0.416604 # CpG site, between exon 6 and 7
(293514,293514)(33653784,33653784)(33653784,33653785)   4       33653784        2       0.429384 # CpG site, between exon 6 and 7


plot_data$fst[plot_data$fst < 0] <- 0

highlight_vals <- c(0.347509, 0.412365, 0.416604, 0.429384)
plot_data$color_flag <- ifelse(round(plot_data$fst, 6) %in% highlight_vals, "highlight", as.character(plot_data$chr))


highlight_color <- "red"
chr_colors <- rep(c(color1, color2), length.out = length(unique(plot_data$chr)))
names(chr_colors) <- unique(plot_data$chr)

# Add the highlight color
color_values <- c(chr_colors, highlight = highlight_color)

ggplot(plot_data, aes(x = BPcum, y = !!sym(metric))) +
  # Background points (non-highlighted)
  geom_point(
    data = subset(plot_data, color_flag != "highlight"),
    aes(color = color_flag),
    alpha = 0.2, size = 8
  ) +
  # Foreground points (highlighted)
  geom_point(
    data = subset(plot_data, color_flag == "highlight"),
    color = "red",  # explicitly set highlight color
    alpha = 1, size = 8
  ) +
  scale_color_manual(values = color_values) +
  scale_x_continuous(
    breaks = range(plot_data$BPcum, na.rm = TRUE)
  ) +
  scale_y_continuous(
    breaks = range(plot_data[[metric]], na.rm = TRUE),
    expand = c(0, 0)
  ) +
  labs(x = NULL, y = NULL) +  # <--- this line removes axis labels
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )


ggsave(filename = file.path("TENM3.fst.for.png"), 
       width = 10, height = 5, units = "in")

cat("Script completed successfully!\n")

# PAR

# Load required packages, installing if necessary
required_packages <- c("qqman", "hexbin", "readr", "ggrepel", "ggplot2", "dplyr", "RColorBrewer", "data.table")
installed_packages <- rownames(installed.packages())

cat("Checking required packages...\n")
for (pkg in required_packages) {
  if (!(pkg %in% installed_packages)) {
    install.packages(pkg, repos = "http://cran.us.r-project.org")
  }
  library(pkg, character.only = TRUE)
}

color1 <- "#ACC772"
color2 <- "#ACC772"
input <- "TENM3.fst.tsv"
metric <- "fst"

# Define parameters
cat("Reading in file...\n")
# Read file
data <- fread(input, sep = "\t", na.strings = c("", "NA"), data.table = TRUE)

# Prepare data for plotting
cat("Preparing data for plotting...\n")
data$chr <- factor(data$chr, levels = c(1, "1A", 2:4, "4A", 5:29, "Z"))

plot_data <- data %>%
  group_by(chr) %>%
  summarise(chr_len = max(midPos)) %>%
  mutate(tot = cumsum(chr_len) - chr_len) %>%
  select(-chr_len) %>%
  left_join(data, by = "chr") %>%
  arrange(chr, midPos) %>%
  mutate(BPcum = midPos + tot)

axisdf <- plot_data %>%
  group_by(chr) %>%
  summarize(center = mean(BPcum))


# Plot
cat("Generating plot...\n")


(298112,298112)(34248653,34248653)(34248653,34248654)   4       34248653        2       0.363803 # CpG, between exon 1 and 2
(298012,298012)(34235390,34235390)(34235390,34235391)   4       34235390        2       0.408015 # CpG, between exon 1 and 2
(297941,297941)(34225827,34225827)(34225827,34225828)   4        v        2       0.485842 # CpG, between exon 1 and 2

highlight_color <- "red"
chr_colors <- rep(c(color1, color2), length.out = length(unique(plot_data$chr)))
names(chr_colors) <- unique(plot_data$chr)


highlight_vals <- c(0.363803, 0.408015, 0.485842)
plot_data$color_flag <- ifelse(round(plot_data$fst, 6) %in% highlight_vals, "highlight", as.character(plot_data$chr))


# Add the highlight color
color_values <- c(chr_colors, highlight = highlight_color)


ggplot(plot_data, aes(x = BPcum, y = !!sym(metric))) +
  # Background points (non-highlighted)
  geom_point(
    data = subset(plot_data, color_flag != "highlight"),
    aes(color = color_flag),
    alpha = 0.2, size = 8
  ) +
  # Foreground points (highlighted)
  geom_point(
    data = subset(plot_data, color_flag == "highlight"),
    color = "red",  # explicitly set highlight color
    alpha = 1, size = 8
  ) +
  scale_color_manual(values = color_values) +
  scale_x_continuous(
    breaks = range(plot_data$BPcum, na.rm = TRUE)
  ) +
  scale_y_continuous(
    breaks = range(plot_data[[metric]], na.rm = TRUE),
    expand = c(0, 0)
  ) +
  labs(x = NULL, y = NULL) +  # <--- this line removes axis labels
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

ggsave(filename = file.path("TENM3.fst.par.png"), 
       width = 10, height = 5, units = "in")

cat("Script completed successfully!\n")