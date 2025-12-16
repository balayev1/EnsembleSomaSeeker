################### Analysis of NGSCheckMate results

## load libraries
library(data.table)
library(ggplot2)
library(ggnewscale)

## load results
corr.table.ngscheckmate <- fread("/scratch.global/balay011/NGSCheckMate/CHORDBAI/ngscheckmate_out/ngscheckmate_out_output_corr_matrix.txt")

## load SRR to subject conversion table
conv.table <- fread("/projects/standard/aventeic/balay011/scripts/WGS_chordomaBai_sampleinfo.txt")

## convert
conv.table$add_subject_id <- paste0("T", conv.table$submitted_subject_id)
conv.table <- conv.table[, add_subject_id := paste0(add_subject_id,
    ifelse(histological_type == "recurrent tumor frozen tissue", "R_tumor",
        ifelse(histological_type == "frozen tissue", "_tumor",
            ifelse(histological_type == "peripheral blood", "_blood", ""))))]

id_map <- setNames(conv.table$add_subject_id, conv.table$Run)
corr.table.ngscheckmate <- corr.table.ngscheckmate[, sample_ID := id_map[sample_ID]]
setnames(corr.table.ngscheckmate, names(corr.table.ngscheckmate)[-1], id_map[names(corr.table.ngscheckmate)[-1]])

## save the converted output
write.table(corr.table.ngscheckmate, "/projects/standard/aventeic/balay011/scripts/ngscheckmate_out_converted_sampleid.txt",
    sep = "\t", row.names = FALSE, quote = FALSE)
## make correlation heatmap
sample_order <- names(corr.table.ngscheckmate)[-1]
corr.dt.long <- melt(corr.table.ngscheckmate,
                     id.vars = "sample_ID",
                     variable.name = "Sample2",
                     value.name = "Correlation")
setnames(corr.dt.long, "sample_ID", "Sample1")
corr.dt.long[, Sample1 := factor(Sample1, levels = sample_order)]
corr.dt.long[, Sample2 := factor(Sample2, levels = sample_order)]

axis_dt <- unique(corr.dt.long[, .(Sample = Sample1)])
setorder(axis_dt, Sample)

### extract Subject ID
axis_dt[, Subject_ID := gsub("R?_(blood|tumor)$", "", as.character(Sample))]

### extract Tissue Type
axis_dt[, Tissue_Type := fifelse(grepl("_blood", Sample), "Blood",
                                 fifelse(grepl("R_tumor", Sample), "Recurrence",
                                         "Tumor"))]

tissue_colors <- c("Blood" = "#1F78B4", "Tumor" = "#E31A1C", "Recurrence" = "#FF7F00")

### define a color palette for Subject IDsunique_subjects <- unique(axis_dt$Subject_ID)
unique_subjects <- unique(axis_dt$Subject_ID)
num_subjects <- length(unique_subjects)
subject_color_vector <- scales::hue_pal()(num_subjects)
set.seed(42)
shuffled_subjects <- sample(unique_subjects)
subject_palette <- setNames(subject_color_vector, shuffled_subjects)

ANNOTATION_SIZE <- 3.0
TOTAL_ANNOTATION_WIDTH <- 2 * ANNOTATION_SIZE # 6.0 units total extent
heatmap_plot <- ggplot(corr.dt.long, aes(x = Sample1, y = Sample2)) +
    geom_tile(aes(fill = Correlation), color = "black", linewidth = 0.1) +
        scale_fill_gradient(low = "white", high = "steelblue",
                        limit = c(0, 1), 
                        name = "Correlation\nValue",
                        guide = guide_colorbar(order = 1)) +
    ggnewscale::new_scale_fill() +
    geom_tile(data = axis_dt, aes(x = Sample, y = 0.5 - ANNOTATION_SIZE/2, fill = Tissue_Type), # Adjust y-center
              height = ANNOTATION_SIZE, 
              inherit.aes = FALSE, color = "black", linewidth = 0.1) +
    geom_tile(data = axis_dt, aes(x = 0.5 - ANNOTATION_SIZE/2, y = Sample, fill = Tissue_Type), # Adjust x-center
              width = ANNOTATION_SIZE, 
              inherit.aes = FALSE, color = "black", linewidth = 0.1) +
    scale_fill_manual(values = tissue_colors, 
                      name = "Tissue Type",
                      guide = guide_legend(order = 2)) +
    ggnewscale::new_scale_fill() +
    geom_tile(data = axis_dt, aes(x = Sample, y = 0.5 - (ANNOTATION_SIZE * 1.5), fill = Subject_ID),
              height = ANNOTATION_SIZE, 
              inherit.aes = FALSE, color = "black", linewidth = 0.1) +
    geom_tile(data = axis_dt, aes(x = 0.5 - (ANNOTATION_SIZE * 1.5), y = Sample, fill = Subject_ID),
              width = ANNOTATION_SIZE, 
              inherit.aes = FALSE, color = "black", linewidth = 0.1) +
    scale_fill_manual(values = subject_palette, 
                      name = "Subject ID",
                      guide = guide_legend(order = 3)) +
    scale_x_discrete(expand = expansion(add = 0.5)) +
    scale_y_discrete(expand = expansion(add = 0.5)) +
    coord_cartesian(xlim = c(0.5 - TOTAL_ANNOTATION_WIDTH, N_samples + 0.5), 
                    ylim = c(0.5 - TOTAL_ANNOTATION_WIDTH, N_samples + 0.5), 
                    expand = FALSE) +
    labs(title = "Sample Correlation CHORDBAI", x = NULL, y = NULL) +
    theme_minimal() +
    theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
        plot.margin = margin(5, 5, 5, 5))

### save the heatmap
ggsave(
    "sample_correlation_chordbai_ngscheckmate.pdf",
    plot = heatmap_plot,
    width = 15,
    height = 15,
    units = "in")

ggsave(
    "sample_correlation_chordbai_ngscheckmate.png",
    plot = heatmap_plot,
    width = 15,
    height = 15,
    units = "in")