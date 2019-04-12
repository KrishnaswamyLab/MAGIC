# To run this file:
# - Set the working directory to 'R/tests'.

library(Rmagic)
library(ggplot2)
library(readr)
library(viridis)

seurat_obj <- function() {
  # load data
  data <- read.csv('../../data/HMLE_TGFb_day_8_10.csv.gz')
  
  seurat_raw_data <- t(data)
  rownames(seurat_raw_data) <- colnames(data)
  colnames(seurat_raw_data) <- rownames(data)
  seurat_obj <- Seurat::CreateSeuratObject(raw.data=seurat_raw_data)

  # run MAGIC
  data_MAGIC <- magic(data)
  seurat_MAGIC <- magic(seurat_obj)
  stopifnot(all(data_MAGIC$result == t(seurat_MAGIC@data)))

  # plot
  p <- ggplot(data) +
         geom_point(aes(VIM, CDH1, colour=ZEB1)) +
         scale_colour_viridis(option="B")
  ggsave('EMT_data_R_before_magic.png', plot=p, width=5, height=5)

  p_m <- ggplot(data_MAGIC) +
           geom_point(aes(VIM, CDH1, colour=ZEB1)) +
           scale_colour_viridis(option="B")
  ggsave('EMT_data_R_after_magic.png', plot=p_m, width=5, height=5)
}
