library(magrittr)
library(ggplot2)
library(gridExtra)
source("load_data.R")
source("implementations.R")
dataset <- import_data()
CIR_fn <- CIR_functions()
L_CIR_fn <- Lamperti_CIR_functions()
DW_fn <- DW_functions()
log_CIR_fn <- Log_CIR_functions(mean(-log(dataset$Ca2)))
exp_DW_fn <- exp_DW_functions(mean(-log(dataset$Ca2)))

sim_par <- optimizer(
  par = DW_fn$inverse_transform(c(1, -1 ,1, 1)),
  fn = local_linearization,
  fs = DW_fn,
  x = dataset$logCa2,
  dt = 0.02
) %>% DW_fn$transform()


double_well_stationary_distribution <- function(x, par = sim_par){
  shape_dens <- function(x){
    exp(2 / par[4]^2 * (-par[3] / 4 * x^4 + par[2] / 2 * x^2 + par[1] * x))
  }
  Z <- integrate(shape_dens, lower = -Inf, upper = Inf)$value
  1/Z * shape_dens(x)
}

double_well_stationary_distribution_original <- function(x, par = sim_par, C = mean(-log(dataset$Ca2))){
  shape_dens <- function(x){
    double_well_stationary_distribution(-log(x) - C) / x
  }
  Z <- integrate(shape_dens, lower = 1, upper = Inf)$value
  shape_dens(x) / Z
}


p1 <- ggplot() + 
  # geom_histogram(data = dataset, aes(x = logCa2, after_stat(density)), col = "black", fill = "grey85", bins = 25) + 
  geom_density(data = dataset, aes(x = logCa2, after_stat(density)), bw = "bcv", kernel = "epanechnikov") +
  stat_function(fun = double_well_stationary_distribution, linewidth = 1, linetype = "dashed") + 
  ggthemes::theme_base() + theme(plot.background = element_blank())


ggsave("Stationary distribution of Double-Well.jpeg", p1)


p2 <- ggplot() + 
  # geom_histogram(data = dataset, aes(x = Ca2, after_stat(density)), col = "black", fill = "grey85", bins = 25) +
  geom_density(data = dataset, aes(x = Ca2, after_stat(density)), bw = bw.bcv(dataset$Ca2, lower = 0.01), kernel = "epanechnikov") +
  stat_function(fun = double_well_stationary_distribution_original, linewidth = 1, linetype = "dashed") + 
  ggthemes::theme_base() + theme(plot.background = element_blank())

ggsave("Stationary distribution of transformed Double-Well.jpeg", p2)

grid.arrange(p1,p2, ncol = 2)

