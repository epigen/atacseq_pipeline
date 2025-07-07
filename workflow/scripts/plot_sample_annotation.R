#### libraries ####
library(data.table)
library(tidyverse)
library(patchwork)
library(ggplot2)
library(ggnewscale)
library(stringr)
# for interactive plotting
library(plotly)
library(htmlwidgets)
# set plot base size
theme_set(theme_minimal(base_size = 6))

#### configs ####
# input
sample_annotation_path <- snakemake@input[["sample_annotation"]]
sample_annotation_w_QC_path <- snakemake@input[["sample_annotation_w_QC"]]

# output
sample_annotation_plot_path <- snakemake@output[["sample_annotation_plot"]]
sample_annotation_html_path <- snakemake@output[["sample_annotation_html"]]

#### load & prepare data ####
# load data
sample_annotation <- data.table::fread(file.path(sample_annotation_path), header = TRUE)
sample_annotation <- data.frame(sample_annotation[!duplicated(sample_annotation[[1]]), ], row.names = 1, check.names = FALSE)

anno <- data.frame(fread(file.path(sample_annotation_w_QC_path), header=TRUE), row.names=1)

# determine QC (pipeline provided) columns
names(sample_annotation) <- gsub(" +", "_", names(sample_annotation)) # replace empty space ` ` with underscore `_`
qc_cols <- setdiff(names(anno), names(sample_annotation))

# determine metadata (user provided) columns by removing non-numeric columns that are unique for each row (e.g., bam_file)
sample_annotation <- sample_annotation %>% select(where(\(.x){ n <- dplyr::n_distinct(na.omit(.x)); (is.numeric(.x) || n < nrow(sample_annotation)) }))
meta_cols <- names(sample_annotation)

# drop columns with no variation (e.g., read_type)
has_var <- function(x) length(unique(na.omit(x))) > 1
qc_cols  <- keep(qc_cols,  ~ has_var(anno[[.x]]))
meta_cols <- keep(meta_cols, ~ has_var(anno[[.x]]))

# collapse "duplicates" (ie redundant columns) only among non-numeric
num_cols  <- meta_cols[vapply(anno[meta_cols], is.numeric, logical(1))]
cat_cols  <- setdiff(meta_cols, num_cols)
sig <- vapply(cat_cols,
              \(col) paste(match(anno[[col]], unique(anno[[col]])), collapse = "|"),
              character(1))
meta_cols <- c(cat_cols[!duplicated(sig)], num_cols)

#### Z-score & cluster QC data ####
qc_mat <- anno |> select(all_of(qc_cols)) |> scale() |> as.matrix()

row_ord <- hclust(dist(qc_mat))$order
col_ord <- hclust(dist(t(qc_mat)))$order

qc_long <- as_tibble(qc_mat[row_ord, col_ord], rownames = "sample") |>
           pivot_longer(-sample, names_to = "metric", values_to = "z")

#### prepare metadata for plotting ####
meta_long <- anno[row_ord, meta_cols]                                     %>% 
  mutate(sample = rownames(anno)[row_ord])                                %>% 
  mutate(across(-sample, as.character))                                   %>% 
  pivot_longer(-sample, names_to = "meta", values_to = "value")           %>% 
  group_by(meta)                                                          %>% 
  mutate(num_val = suppressWarnings(as.numeric(value)),
         type    = if (all(!is.na(num_val))) "numeric" else "factor",
         col     = if (type[1] == "numeric") {
                      scales::col_numeric("plasma",
                                          domain = range(num_val, na.rm = TRUE))(num_val)
                    } else {
                      pal <- scales::hue_pal(l = 65)(n_distinct(value))
                      setNames(pal, sort(unique(value)))[value]
                    }
         )                                                   %>% 
  ungroup()                                                               %>% 
  select(-num_val)

#### embed ALL metadata into a tooltip string of interactive plot ####
meta_txt <- anno[row_ord, meta_cols] %>% mutate(across(everything(), as.character))
meta_txt <- pmap_chr(meta_txt, \(...) {
               vals <- c(...)
               paste(paste(names(vals), vals, sep = ": "), collapse = "<br>")
           })
names(meta_txt) <- rownames(anno)[row_ord]

#### plot heatmaps ####

#### QC heatmap ####

# add (un-scaled) metric values (raw) for tooltip of interactive plot
qc_raw_long <- anno[row_ord, qc_cols] %>%
               mutate(sample = rownames(anno)[row_ord]) %>% 
               pivot_longer(-sample, names_to = "metric", values_to = "raw")

qc_long <- qc_long %>% left_join(qc_raw_long, by = c("sample", "metric"))

# keep the clustered order in the plot and add hover-tooltip
qc_long <- qc_long %>% 
  mutate(sample = factor(sample,  levels = rownames(qc_mat)[row_ord]), # row-order
         metric = factor(metric,  levels = colnames(qc_mat)[col_ord]), # col-order
         hover  = paste0("Sample: ", sample,
                         "<br>Metric: ", metric,
                         "<br>Value: ",   signif(raw, 4),
                         "<br>", meta_txt[as.character(sample)])
        )

# plot
p_qc <- ggplot(qc_long, aes(x = metric,
                            y = sample,
                            fill = z,
                            text = hover)) +
        geom_tile() +
        scale_x_discrete(limits = colnames(qc_mat)[col_ord]) +           # enforce col order
        scale_y_discrete(limits = rownames(qc_mat)[row_ord]) +           # enforce row order
        scale_fill_gradient2(low = "blue", mid = "white", high = "red", name = "z-score",
                     guide = guide_colourbar(barheight = 2,  # thinner
                                             barwidth  = 0.15)) +
        labs(x = NULL, y = NULL, title = "QC metrics (scaled)") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              panel.grid = element_blank())

#### metadata heatmap as "annotation" ####
p_meta <- NULL
if(length(meta_cols) > 0){
    # order columns (x) exactly like the QC heatmap
    meta_long <- meta_long %>% mutate(sample = factor(sample,  levels = rownames(qc_mat)[row_ord]))   # row-order
    meta_levels <- unique(meta_long$meta)
    meta_long   <- meta_long %>% mutate(meta = factor(meta, levels = meta_levels))
    
    p_meta <- ggplot() +
              scale_y_discrete(limits = levels(qc_long$sample)) +
              scale_x_discrete(limits = meta_levels) +
                labs(x = NULL, y = NULL, title = "Metadata") +
              theme(axis.text.x = element_text(angle = 45, hjust = 1),
                    axis.text.y = element_blank(),
                    axis.title   = element_blank(),
                    panel.grid   = element_blank())
    
    for (v in meta_levels) {
        dat <- dplyr::filter(meta_long, meta == v)
    
        p_meta <- p_meta + ggnewscale::new_scale_fill()   # reset “fill” for this column
    
        if (dat$type[1] == "numeric") {                  # continuous legend
            p_meta <- p_meta +
                geom_tile(data = dat, aes(x = meta, y = sample, fill = as.numeric(value)), colour = "grey60", linewidth = 0.1) +
                scale_fill_viridis_c(name = v, option = "plasma",
                     guide = guide_colourbar(barheight = 2,  # thinner
                                             barwidth  = 0.15)) 
        } else {                                         # categorical legend
            pal <- setNames(dat$col, dat$value)

            # reduce legend in case of more than 10 levels
            max_items   <- min(10, length(unique(dat$value)))
            all_levels  <- unique(names(pal))
            show_levels <- all_levels[1:max_items]
            
            p_meta <- p_meta +
                geom_tile(data = dat, aes(x = meta, y = sample, fill = value), colour = "grey60", linewidth = 0.1) +
                scale_fill_manual(values = pal, 
                                  # name = v,
                                  breaks = show_levels,
                                  guide = guide_legend(keywidth  = 0.25,
                                                       keyheight = 0.4,
                                                       ncol=1,
                                                       byrow = TRUE,
                                                       title = ifelse(
                                                           length(all_levels) <= max_items,
                                                           v,
                                                           paste0(v, " (showing ", max_items, "/", length(all_levels), ")")
                                                           )
                                                      )
                                 )
        }
    }
}

#### combine and save plots ####
p_combined <- if (is.null(p_meta)) p_qc else (p_qc | p_meta) + plot_layout(widths = c(length(qc_cols), length(meta_cols)), guides = "collect") & theme(legend.position = "right")

# determine sizes
n_rows <- nrow(qc_mat)
n_cols <- length(qc_cols) + length(meta_cols)
max_row_label <- max(nchar(rownames(anno)))
max_col_label <- max(nchar(c(qc_cols, meta_cols)))

height_in <- n_rows * 0.08 + max_col_label * 0.05 + 1
width_in  <- n_cols * 0.10 + max_row_label * 0.05 + 2

# options(repr.plot.width = width_in, repr.plot.height = height_in)
# p_combined

ggsave(sample_annotation_plot_path, plot = p_combined, width = width_in, height = height_in, units = "in", dpi = 300)

#### interactive plot ####
# determine sizes in pixels
width_px  <- round((length(qc_cols) * 0.10 + max_row_label * 0.05 + 2) * 96)
height_px <- round(height_in * 96)

p_qc_interactive  <- plotly::ggplotly(p_qc,  tooltip = "text", width = width_px, height = height_px)

# p_qc_interactive

htmlwidgets::saveWidget(p_qc_interactive, sample_annotation_html_path, selfcontained = TRUE, title = "Sample annotation")
