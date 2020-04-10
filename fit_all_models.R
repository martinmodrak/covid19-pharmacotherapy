library(tidyverse)
library(tidybayes)
library(rstan)
library(cmdstanr)
library(readxl)
library(here)
library(cowplot)
library(bayesplot)
theme_set(theme_cowplot())

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)


source("R/hmm_tools.R")
source("R/models_tools.R")
source("R/model_defs.R")
source("R/load_data.R")


gautret <- load_gautret_data()
model <- rstan::stan_model(file = here::here("hmm.stan"))

fit_all_models(model, gautret)
