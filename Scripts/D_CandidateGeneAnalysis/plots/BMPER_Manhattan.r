Make manhattan plots of BMPER genes

### BMPER
#### NC_044572.1:46478775-46625955

# BMPER.fst.tsv

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
input <- "BMPER.fst.tsv"
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
(400163,400163)(46511099,46511099)(46511099,46511100)   2       46511099        2       0.246331 # CpG site!!! Noncoding, between exon 3 and 4

plot_data$fst[plot_data$fst < 0] <- 0


plot_data$color_flag <- ifelse(round(plot_data$fst, 6) == 0.246331, "highlight", as.character(plot_data$chr))
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


ggsave(filename = file.path("BMPER.fst.cra.png"), 
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
input <- "BMPER.fst.tsv"
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


(400810,400810)(46580828,46580828)(46580828,46580829)   2       46580828        2       0.303085 # CpG site between exon 13 and 14 (just after 13)
(400849,400849)(46584598,46584598)(46584598,46584599)   2       46584598        2       0.327521 # CpG site between exon 9 and 10
(400813,400813)(46580935,46580935)(46580935,46580936)   2       46580935        2       0.331082 # CpG site between exon 13 and 14 (just after 13)
(400924,400924)(46602867,46602867)(46602867,46602868)   2       46602867        2       0.376037 # CpG site between exon 14 and 15
(400783,400783)(46576940,46576940)(46576940,46576941)   2       46576940        2       0.402065 # CpG site between exon 12 and 13

plot_data$fst[plot_data$fst < 0] <- 0

highlight_vals <- c(0.303085, 0.327521, 0.331082, 0.376037, 0.402065)
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


ggsave(filename = file.path("BMPER.fst.for.png"), 
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
input <- "BMPER.fst.tsv"
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


(400836,400836)(46578649,46578649)(46578649,46578650)   2       46578649        2       0.245141 # CpG site between exon 12 and 13
(400310,400310)(46513495,46513495)(46513495,46513496)   2       46513495        2       0.257733 # CpG site between exon 3 and 4, very close to 4 (46514490, 995 sites)
(400794,400794)(46574902,46574902)(46574902,46574903)   2       46574902        2       0.258011 # CpG, adjacent to noncoding nonCpG site two rows down, between exon 9 and 10

highlight_color <- "red"
chr_colors <- rep(c(color1, color2), length.out = length(unique(plot_data$chr)))
names(chr_colors) <- unique(plot_data$chr)


highlight_vals <- c(0.245141, 0.257733, 0.258011)
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

ggsave(filename = file.path("BMPER.fst.par.png"), 
       width = 10, height = 5, units = "in")

cat("Script completed successfully!\n")