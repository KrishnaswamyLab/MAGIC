# To run this file:
# - Set the working directory to 'R/tests'.
# - Unzip the data files in 'data'.

source('../R/run_magic.R')

test_magic <- function() {
  require(ggplot2)

	# load data
	data <- read.csv('../data/HMLE_TGFb_day_8_10.csv', header=TRUE, sep=',')

	# run MAGIC
	data_MAGIC <- run_magic(data, rescale_percent=0.99)

	# plot
	p <- ggplot(data)
	p + geom_point(aes(VIM, CDH1, colour=ZEB1)) + scale_colour_gradient(low = 'purple', high='yellow')
  ggsave('EMT_data_R_before_magic.png', width=5, height=5)

  p_m <- ggplot(data_MAGIC)
  p_m + geom_point(aes(VIM, CDH1, colour=ZEB1)) + scale_colour_gradient(low = 'purple', high='yellow')
  ggsave('EMT_data_R_after_magic.png', width=5, height=5)
}
