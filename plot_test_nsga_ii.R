# to plot the estimated pareto fronts from the NSGA-II implementation, to compare with the plot in the paper.
library(tidyverse)

plot_one_res <- function(res_file, out_plot_name) {
  fitnesses <- read_csv(res_file)
  # all are biobjective problems
  p <- ggplot(fitnesses) +
    geom_point(aes(x = objective0, y = objective1))
  ggsave(filename = out_plot_name, plot = p)
}

plot_one_res("sch-test-fronts.csv", "sch-test-plot.png")
plot_one_res("kur-test-fronts.csv", "kur-test-plot.png")
plot_one_res("zdt2-test-fronts.csv", "zdt2-test-plot.png")
plot_one_res("zdt4-test-fronts.csv", "zdt4-test-plot.png")
