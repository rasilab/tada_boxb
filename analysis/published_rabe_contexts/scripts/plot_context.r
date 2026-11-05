suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 6)
annotations_file <- args[1]
source_file <- args[2]
panel_file <- args[3]
png_file <- args[4]
pdf_file <- args[5]
summary_files <- args[-(1:5)]
font_family <- "Helvetica"
font_size <- 7
panel_width <- 1.0
annotation_width <- 1.5
figure_height <- 1.9
png_dpi <- 96
bases <- c("A", "C", "G", "T")
rna_labels <- c("A", "C", "G", "U")
scale_step <- 10
fill_low <- "#F7FBFF"
fill_high <- "#6BAED6"
fill_missing <- "grey80"
number_color <- "#111111"
legend_title <- "Mean % edited"

annotations <- read_csv(annotations_file, show_col_types = FALSE)
data <- bind_rows(lapply(summary_files, read_csv, show_col_types = FALSE))
stopifnot(!anyDuplicated(annotations$sample_name),
          setequal(data$sample_name, annotations$sample_name),
          all(table(data$sample_name) == length(bases)^2),
          !anyDuplicated(data[c("sample_name", "upstream_base", "downstream_base")]),
          all(data$upstream_base %in% bases), all(data$downstream_base %in% bases),
          !anyNA(data$n_sites), !anyNA(data$sum_editing_fraction),
          all(data$n_sites >= 0), all(data$sum_editing_fraction >= 0),
          all(data$sum_editing_fraction <= data$n_sites))
data <- left_join(data, annotations, by = "sample_name")
condition_order <- annotations %>% distinct(condition, condition_order) %>% arrange(condition_order)
source <- data %>%
  group_by(condition, upstream_base, downstream_base) %>%
  summarise(n_samples = n(), n_samples_with_sites = sum(n_sites > 0),
            sum_editing_fraction = sum(sum_editing_fraction), n_sites = sum(n_sites),
            .groups = "drop") %>%
  mutate(mean_percent = if_else(n_sites > 0, 100 * sum_editing_fraction / n_sites, NA_real_))
write_csv(source, source_file, na = "")
dir.create(dirname(panel_file), showWarnings = FALSE, recursive = TRUE)
write_csv(source, panel_file, na = "")
stopifnot(any(is.finite(source$mean_percent)))
scale_max <- max(scale_step, ceiling(max(source$mean_percent, na.rm = TRUE) / scale_step) * scale_step)
figure_width <- panel_width * nrow(condition_order) + annotation_width

theme_set(theme_classic(base_family = font_family, base_size = font_size) +
            theme(text = element_text(size = font_size, family = font_family),
                  axis.text = element_text(size = font_size, color = "black"),
                  axis.title = element_text(size = font_size),
                  axis.line = element_blank(), axis.ticks = element_blank(),
                  strip.text = element_text(size = font_size), strip.background = element_blank(),
                  legend.position = "right", legend.text = element_text(size = font_size),
                  legend.title = element_text(size = font_size),
                  legend.key.height = grid::unit(0.15, "in")))
plot_data <- source %>%
  mutate(condition = factor(condition, levels = condition_order$condition),
         upstream_base = factor(upstream_base, levels = bases),
         downstream_base = factor(downstream_base, levels = bases))
p <- ggplot(plot_data, aes(downstream_base, upstream_base, fill = mean_percent)) +
  geom_tile() +
  geom_text(aes(label = if_else(is.na(mean_percent), "", sprintf("%.0f", mean_percent))),
            color = number_color, family = font_family, size = font_size / ggplot2::.pt) +
  facet_wrap(~condition, nrow = 1, drop = FALSE) +
  scale_fill_gradient(low = fill_low, high = fill_high, na.value = fill_missing,
                      limits = c(0, scale_max), name = legend_title) +
  scale_x_discrete(drop = FALSE, limits = bases, labels = rna_labels) +
  scale_y_discrete(drop = FALSE, limits = rev(bases), labels = rev(rna_labels)) +
  coord_fixed() +
  labs(x = "3′ flanking base", y = "5′ flanking base")
ggsave(png_file, p, width = figure_width, height = figure_height, dpi = png_dpi)
ggsave(pdf_file, p, width = figure_width, height = figure_height, device = cairo_pdf)
if (file.info(pdf_file)$size > 2 * file.info(png_file)$size) {
  message("PDF exceeds twice the PNG size; both files have been retained.")
}
