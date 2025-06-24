Make manhattan plots of NUDT4 genes

### NUDT4
#### NC_044586.1:28599179-28626256	

# NUDT4.fst.tsv

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
input <- "NUDT4.fst.tsv"
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
(234395,234395)(28603313,28603313)(28603313,28603314)   1A      28603313        2       0.244754 # CpG, between exon 4 and 5

plot_data$fst[plot_data$fst < 0] <- 0


plot_data$color_flag <- ifelse(round(plot_data$fst, 6) == 0.244754, "highlight", as.character(plot_data$chr))
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


ggsave(filename = file.path("NUDT4.fst.cra.png"), 
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
input <- "NUDT4.fst.tsv"
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

(234581,234581)(28621591,28621591)(28621591,28621592)   1A      28621591        2       0.205953 # CpG, between exon 1 and 2


plot_data$fst[plot_data$fst < 0] <- 0

highlight_vals <- c(0.205953)
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


ggsave(filename = file.path("NUDT4.fst.for.png"), 
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
input <- "NUDT4.fst.tsv"
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


(234496,234496)(28600154,28600154)(28600154,28600155)   1A      28600154        2       0.241479 # 5`UTR, CpG
# (234506,234506)(28601802,28601802)(28601802,28601803)   1A      28601802        2       0.246576 # 5`UTR, not CpG
(234568,234568)(28614781,28614781)(28614781,28614782)   1A      28614781        2       0.249894 # CpG, between exon 1 and 2
(234505,234505)(28601647,28601647)(28601647,28601648)   1A      28601647        2       0.319059 # CpG within 5'UTR


highlight_color <- "red"
chr_colors <- rep(c(color1, color2), length.out = length(unique(plot_data$chr)))
names(chr_colors) <- unique(plot_data$chr)


highlight_vals <- c(0.241479, 0.249894, 0.319059)
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

ggsave(filename = file.path("NUDT4.fst.par.png"), 
       width = 10, height = 5, units = "in")

cat("Script completed successfully!\n")