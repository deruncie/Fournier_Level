library(microbenchmark)
library(MCMCpack)
library(BSFG)


# load data
data = read.delim('Fournier-Leveletal_DataSEM.txt')[,1:16]

data = data[,c('UniqueID','Family','Loc','PTUtB','PTUtF','LfL','Branch','SIL_num')]
data$Survival = data$SIL_num > 0 & data$Branch > 0

data$log_Branch = log(data$Branch + 1)
data$log_SIL_num = log(data$SIL_num + 1)

Y = data[,c('PTUtB','PTUtF','LfL','log_Branch','log_SIL_num')] #,'Survival'

# create a folder for holding the posterior samples of the current chain (multiple folders could be used for different chains)
rep = "1"
folder = sprintf('Rep_%s',rep)
try(dir.create(folder))
setwd(folder)

print('Initializing')


run_parameters = BSFG_control(
  sampler = 'fast_BSFG',
  # sampler = 'general_BSFG',
  scale_Y = TRUE,
  simulation = FALSE,
  h2_divisions = 20,
  h2_step_size = 0.3,
  burn = 100,
  k_init = 5
)

priors = BSFG_priors(
  fixed_var = list(V = 1,     nu = 3),
  # tot_Y_var = list(V = 0.5,   nu = 3),
  tot_Y_var = list(V = 0.5,   nu = 10),
  tot_F_var = list(V = 18/20, nu = 20),
  h2_priors_resids_fun = function(h2s,n) 1,#pmax(pmin(ddirichlet(c(h2s,1-sum(h2s)),rep(2,length(h2s)+1)),10),1e-10),
  h2_priors_factors_fun = function(h2s,n) 1,#ifelse(h2s == 0,n,n/(n-1))
  Lambda_prior = list(
    sampler = sample_Lambda_prec_ARD,
    Lambda_df = 3,
    delta_1   = list(shape = 2.1,  rate = 1/20),
    delta_2   = list(shape = 3, rate = 1)
  ),
  # Lambda_prior = list(
  #   sampler = sample_Lambda_prec_ARD_v2,
  #   Lambda_df = 3,
  #   delta_1   = list(shape = 2,  rate = 40),
  #   delta_2   = list(shape = 2, rate = 3)
  # ),
  # Lambda_prior = list(
  #   sampler = sample_Lambda_prec_TPB,
  #   Lambda_A      = .5,
  #   Lambda_B      = .5,
  #   delta_1   = list(shape = 2.1,  rate = 1/20),
  #   delta_2   = list(shape = 3, rate = 1)
  # ),
  B_prior = list(
    sampler = sample_B_prec_ARD,
    B_df      = 3,
    B_F_df    = 3
  )
  # B_prior = list(
  #   sampler = sample_B_prec_TPB,
  #   B_A      = .5,
  #   B_B      = .5,
  #   B_omega  = 1/10,
  #   B_F_A      = .5,
  #   B_F_B      = .5,
  #   B_F_omega  = 1/10
  # )
)

BSFG_state = BSFG_init(Y, model=~0+Loc+(0+Loc|UniqueID), data,
                       run_parameters=run_parameters,
                       priors=priors)


BSFG_state = reorder_factors(BSFG_state)
BSFG_state = clear_Posterior(BSFG_state)
n_samples = 100;
for(i  in 1:20) {
  if(i == 10){
    BSFG_state = reorder_factors(BSFG_state)
    BSFG_state = clear_Posterior(BSFG_state)
  }
  print(sprintf('Run %d',i))
  BSFG_state = sample_BSFG(BSFG_state,n_samples,grainSize=1,ncores = 1)
  if(BSFG_state$Posterior$total_samples>0) trace_plot(BSFG_state$Posterior$tot_F_prec[,1,])
  if(BSFG_state$Posterior$total_samples>0) trace_plot(log(BSFG_state$Posterior$delta[,1,]))
  if(BSFG_state$current_state$nrun < BSFG_state$run_parameters$burn) {
    BSFG_state = reorder_factors(BSFG_state)
    # BSFG_state$current_state = update_k(BSFG_state)
    BSFG_state = clear_Posterior(BSFG_state)
    BSFG_state$run_parameters$burn = max(c(BSFG_state$run_parameters$burn,BSFG_state$current_state$nrun+100))
    print(BSFG_state$run_parameters$burn)
  }
  BSFG_state = save_posterior_chunk(BSFG_state)
  print(BSFG_state)
  plot(BSFG_state)
}


BSFG_state$Posterior = reload_Posterior(BSFG_state)
