library(readr)
magic_testdata <- read_csv("../../data/HMLE_TGFb_day_8_10.csv.gz")
set.seed(42)
keep_cols <- colSums(magic_testdata > 0) > 10
keep_rows <- rowSums(magic_testdata) > 2000
magic_testdata <- magic_testdata[keep_rows,keep_cols]
magic_testdata <- Rmagic::library.size.normalize(magic_testdata)
magic_testdata <- sqrt(magic_testdata)
select_cols <- c(colnames(magic_testdata)[ceiling(runif(200) * nrow(magic_testdata))],
                 c("VIM", "CDH1", "ZEB1"))
magic_testdata <- magic_testdata[,colnames(magic_testdata) %in% select_cols]
select_rows <- ceiling(runif(500) * nrow(magic_testdata))
magic_testdata <- magic_testdata[select_rows,]
write_csv(magic_testdata, "../../data/test_data.csv")
usethis::use_data(magic_testdata)
