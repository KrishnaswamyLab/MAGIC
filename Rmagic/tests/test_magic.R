# To run this file:
# - Set the working directory to 'R/tests'.

library(Rmagic)
library(ggplot2)
library(readr)
library(viridis)

test_magic <- function() {
  # load data
  data <- read.csv('../../data/HMLE_TGFb_day_8_10.csv.gz')

  # run MAGIC
  data_MAGIC <- magic(data)

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

test_seurat <- function() {
  data(magic_testdata)

  # seurat_obj <- Seurat::CreateSeuratObject(raw.data=t(magic_testdata))
  seurat_obj <- Seurat::CreateSeuratObject(counts = t(x = magic_testdata))

  # run MAGIC
  data_MAGIC <- magic(data = magic_testdata, seed = 42)
  seurat_obj <- magic(data = seurat_obj, seed = 42)
  stopifnot(all(data_MAGIC$result == t(x = Seurat::GetAssayData(object = seurat_obj, slot = 'data', assay = 'MAGIC'))))
}
