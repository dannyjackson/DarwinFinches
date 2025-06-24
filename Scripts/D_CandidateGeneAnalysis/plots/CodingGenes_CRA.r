Make manhattan plots of coding genes in CRA

### HPSE
#### NC_044574.1:59804605-59815556	

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
input <- "HPSE.fst.tsv"
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
# second highest snp, 59815456, is in the 3' noncoding but transcribed region of exon 12 (3'UTR)

plot_data$fst[plot_data$fst < 0] <- 0


plot_data$color_flag <- ifelse(round(plot_data$fst, 6) == 0.327240, "highlight", as.character(plot_data$chr))
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


ggsave(filename = file.path("HPSE.fst.cra.png"), 
       width = 10, height = 5, units = "in")

cat("Script completed successfully!\n")


## PODXL
## NC_044586.1:67365108-67411859


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
input <- "PODXL.fst.tsv"
metric <- "fst"
pop_name<- "cra"

# Define parameters
cat("Reading in file...\n")
# Read file
data <- fread(input, sep = "\t", na.strings = c("", "NA"), data.table = TRUE)
colnames(data) <- c("region", "chr", "midPos", "nSNPS", "fst")
# Prepare data for plotting
cat("Preparing data for plotting...\n")
data$chr <- factor(data$chr, levels = c(1, "1A", 2:4, "4A", 5:29, "Z"))

data = data[-1,]


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
# (577585,577585)(67381351,67381351)(67381351,67381352)   1A      67381351        2       0.448900 # within exon 2, near 67381354

plot_data$fst[plot_data$fst < 0] <- 0


plot_data$color_flag <- ifelse(round(plot_data$fst, 6) == 0.448900, "highlight", as.character(plot_data$chr))
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


ggsave(filename = file.path("PODXL.fst.cra.png"), 
       width = 10, height = 5, units = "in")

cat("Script completed successfully!\n")

