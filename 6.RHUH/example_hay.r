library(virosolver)
library(lazymcmc)
library(tidyverse)
library(ggplot2)
library(arrow)
library(extraDistr)
library(deSolve)
library(doSNOW)
library(doParallel)
library(doMPI)


## Attach simulated data
leb_ct_data <- arrow::read_feather("Data_Ct.feather")
ct_detect <- 40

## Plot only detectable Ct values
p_ct_data <- ggplot(leb_ct_data %>% filter(ct < 40)) + 
  geom_violin(aes(x=t,group=t,y=ct),scale="width",fill="grey70",draw_quantiles=c(0.025,0.5,0.975)) + 
  geom_jitter(aes(x=t,y=ct),size=0.1,width=2,height=0) + 
  scale_y_continuous(trans="reverse") +
  theme_bw() +
  ylab("Ct value") +
  xlab("Observation time") +
  ggtitle("Observed Ct values < 40 (the limit of detection) over time")
p_ct_data

p_detectable_data <- leb_ct_data %>% 
  mutate(detect=ct < 30) %>% 
  group_by(t) %>% 
  summarize(prev=sum(detect)/n()) %>% 
  ggplot() + geom_point(aes(x=t,y=prev)) + 
  theme_bw() + 
  scale_y_continuous(limits=c(0,0.6)) + 
  ylab("Proportion detectable") +
  ggtitle("Proportion of samples with Ct values < 40") +
  xlab("Observation time") 
## `summarise()` ungrouping output (override with `.groups` argument)
p_detectable_data

data(example_gp_partab)
pars <- example_gp_partab$values
names(pars) <- example_gp_partab$names

## Solve the Ct model over a range of times since infection (referred to as "ages")
test_ages <- seq(1,50,by=1)

## This gives the modal Ct value
cts <- viral_load_func(pars, test_ages)

p_ct_model <- ggplot(data.frame(ct=c(40,cts),t=c(0,test_ages))) + 
  geom_line(aes(x=t,y=ct)) + 
  scale_y_continuous(trans="reverse",
                     limits=c(40,10)) +
  theme_bw() +
  ylab("Modal Ct value") +
  xlab("Days since infection")

## Note that this model does not solve for t=0, 
## as it is always assumed that no one is detectable 0 days post infection
prop_detect <- prop_detectable(test_ages,pars, cts)
p_ct_model_detectable <- ggplot(data.frame(p=c(0,prop_detect),t=c(0,test_ages))) + 
  geom_line(aes(x=t,y=p)) + 
  theme_bw() +
  ylab("Proportion of infections\n still detectable") +
  xlab("Days since infection")
p_ct_model
p_ct_model_detectable

## Draw samples from probabilistic Ct model
sim_cts <- simulate_viral_loads_example(test_ages, pars,N=200)
print(head(sim_cts))
## # A tibble: 6 x 3
##     age i        ct
##   <dbl> <chr> <dbl>
## 1     1 1      40  
## 2     1 2      40  
## 3     1 3      35.3
## 4     1 4      40  
## 5     1 5      32.5
## 6     1 6      35.2
p_sim_cts_age <- ggplot(sim_cts %>% filter(ct < 40)) +
  geom_point(aes(x=age,y=ct),alpha=0.25) +
  scale_y_continuous(trans="reverse",limits=c(40,10)) +
  theme_bw() +
  ylab("Ct value") +
  xlab("Days since infection") +
  ggtitle("Simulated detectable Ct values on each day post infection")
p_sim_cts_age

############################################################
# Define Markov Chain Monte Carlo parameters
data(example_gp_partab)
head(example_gp_partab)
##       values        names fixed lower_bound upper_bound steps lower_start
## 1  0.5000000 overall_prob     0           0           1   0.1         0.0
## 2  0.0000000       tshift     1           0           3   0.1         0.0
## 3  5.0000000 desired_mode     1           0           7   0.1         0.0
## 4 19.7359875   viral_peak     0           0          40   0.1        15.0
## 5  5.0000000       obs_sd     0           0          25   0.1         1.0
## 6  0.7888288       sd_mod     1           0           1   0.1         0.4
##   upper_start
## 1         1.0
## 2        10.0
## 3        10.0
## 4        25.0
## 5        10.0
## 6         0.6
## Illustration -- set the `viral_peak` parameter to be estimated during the procedure, and the `intercept` parameter to be fixed
example_gp_partab <- example_gp_partab %>% filter(names == "viral_peak") %>% mutate(fixed=0)
example_gp_partab <- example_gp_partab %>% filter(names == "intercept") %>% mutate(fixed=1)

## Priors on incidence

## Read in the SEIR model parameter control table
data(example_seir_partab)
## Pull out the current values for each parameter, and set these as the prior means
means <- example_seir_partab$values
names(means) <- example_seir_partab$names
## Set standard deviations of prior distribution
sds_seir <- c("obs_sd"=0.5,"viral_peak"=2,
         "wane_rate2"=1,"t_switch"=3,"level_switch"=1,
         "prob_detect"=0.03,
         "incubation"=0.25, "infectious"=0.5)

## Define a function that returns the log prior probability for a given vector of parameter
## values in `pars`, given the prior means and standard deviations set above.
prior_func_seir <- function(pars,...){
    ## Ct model priors
    obs_sd_prior <- dnorm(pars["obs_sd"], means[which(names(means) == "obs_sd")], sds_seir["obs_sd"],log=TRUE)
    viral_peak_prior <- dnorm(pars["viral_peak"], means[which(names(means) == "viral_peak")], sds_seir["viral_peak"],log=TRUE)
    wane_2_prior <- dnorm(pars["wane_rate2"],means[which(names(means) == "wane_rate2")],sds_seir["wane_rate2"],log=TRUE)
    tswitch_prior <- dnorm(pars["t_switch"],means[which(names(means) == "t_switch")],sds_seir["t_switch"],log=TRUE)
    level_prior <- dnorm(pars["level_switch"],means[which(names(means) == "level_switch")],sds_seir["level_switch"],log=TRUE)
    ## Beta prior on the prob_detect parameter to ensure between 0 and 1
    beta1_mean <- means[which(names(means) == "prob_detect")]
    beta1_sd <- sds_seir["prob_detect"]
    beta_alpha <- ((1-beta1_mean)/beta1_sd^2 - 1/beta1_mean)*beta1_mean^2
    beta_beta <- beta_alpha*(1/beta1_mean - 1)
    beta_prior <- dbeta(pars["prob_detect"],beta_alpha,beta_beta,log=TRUE)

    ## SEIR model priors
    incu_prior <- dlnorm(pars["incubation"],log(means[which(names(means) == "incubation")]), sds_seir["incubation"], TRUE)
    infectious_prior <- dlnorm(pars["infectious"],log(means[which(names(means) == "infectious")]),sds_seir["infectious"],TRUE)

    ## Sum up
    obs_sd_prior + viral_peak_prior + 
        wane_2_prior + tswitch_prior + level_prior + beta_prior +
        incu_prior + infectious_prior
}

## Likelihood and posterior function

## Point to a function that expects a vector of named parameters and returns a vector of daily infection probabilities/incidence
incidence_function <- solveSEIRModel_lsoda_wrapper

## Use the example parameter table
data(example_seir_partab)

## Create the posterior function used in the MCMC framework
posterior_func <- create_posterior_func(parTab=example_seir_partab,
                                        data=leb_ct_data,
                                        PRIOR_FUNC=prior_func_seir,
                                        INCIDENCE_FUNC=incidence_function,
                                        use_pos=FALSE) ## Important argument, see text

## Test with default parameters to find the log likelihood
posterior_func(example_seir_partab$values)
##    obs_sd 
## -6006.231

## Filter to Ct values on day 111 (fifth timepoint)
ct_data_use <- leb_ct_data %>% filter(t == 30)
p_ct_use <- ggplot(ct_data_use %>% filter(ct < 40)) + geom_histogram(aes(x=ct)) + theme_bw()
p_ct_use
## `stat_bin()` using `bins = 30`. Pick better value with `binwidth

## Function from the virosolver package to generate random starting parameter values that return a finite likelihood
start_tab <- generate_viable_start_pars(example_seir_partab,ct_data_use,
                                        create_posterior_func,
                                        incidence_function,
                                        prior_func_seir)

covMat <- diag(nrow(start_tab))
mvrPars <- list(covMat,2.38/sqrt(nrow(start_tab[start_tab$fixed==0,])),w=0.8)
mcmc_pars <- c("iterations"=200000,"popt"=0.234,"opt_freq"=2000,
              "thin"=1000,"adaptive_period"=100000,"save_block"=100)

dir.create("mcmc_chains/readme_single_cross_section",recursive=TRUE)

##################################
## RUN THE MCMC FRAMEWORK
## Run 3 MCMC chains. Note that it is possible to parallelize this loop with foreach and doPar
## Note the `use_pos` argument needs to be set here too
nchains <- 3
res <- foreach(chain_no=1:nchains,.packages = c("virosolver","lazymcmc","extraDistr","tidyverse","patchwork")) %dopar% {
  outputs <- run_MCMC(parTab=start_tab,
                      data=ct_data_use,
                      INCIDENCE_FUNC=incidence_function,
                      PRIOR_FUNC=prior_func_seir,
                      mcmcPars=mcmc_pars,
                      filename=paste0("mcmc_chains/readme_single_cross_section/readme_seir_",chain_no),
                      CREATE_POSTERIOR_FUNC=create_posterior_func,
                      mvrPars=mvrPars,
                      use_pos=FALSE) ## Important argument
}

## Read in the MCMC chains
chains <- load_mcmc_chains(location="mcmc_chains/readme_single_cross_section",
                           parTab=start_tab,
                           burnin=mcmc_pars["adaptive_period"],
                           chainNo=TRUE,
                           unfixed=TRUE,
                           multi=TRUE)
## [1] "mcmc_chains/readme_single_cross_section/readme_seir_1_multivariate_chain.csv"
## [2] "mcmc_chains/readme_single_cross_section/readme_seir_2_multivariate_chain.csv"
## [3] "mcmc_chains/readme_single_cross_section/readme_seir_3_multivariate_chain.csv"
## [[1]]
## [1] 201
## 
## [[2]]
## [1] 201
## 
## [[3]]
## [1] 201
## Reshape for plotting
chains_melted <- chains$chain %>% as_tibble %>% group_by(chain) %>% mutate(sampno=1:n()) %>% pivot_longer(-c(sampno,chain))
## Look at trace plots
p_trace <- ggplot(chains_melted) + 
  geom_line(aes(x=sampno,y=value,col=as.factor(chain))) + 
  facet_wrap(~name,scales="free_y") + 
  scale_color_viridis_d(name="Chain") + 
  theme_bw() +
  xlab("Iteration") +
  ylab("Value")
p_trace


## Load in MCMC chains again, but this time read in the fixed parameters too 
## to ensure that the posterior draws are compatible with the model functions
chains <- load_mcmc_chains(location="mcmc_chains/readme_single_cross_section",
                           parTab=start_tab,
                           burnin=mcmc_pars["adaptive_period"],
                           chainNo=FALSE,
                           unfixed=FALSE,
                           multi=TRUE)
## [1] "mcmc_chains/readme_single_cross_section/readme_seir_1_multivariate_chain.csv"
## [2] "mcmc_chains/readme_single_cross_section/readme_seir_2_multivariate_chain.csv"
## [3] "mcmc_chains/readme_single_cross_section/readme_seir_3_multivariate_chain.csv"
## [[1]]
## [1] 201
## 
## [[2]]
## [1] 201
## 
## [[3]]
## [1] 201
## Do some reshaping to allow correct subsampling (we need each sampno to correspond to one unique posterior draw)
chain_comb <- chains$chain %>% as_tibble() %>% mutate(sampno=1:n()) %>% as.data.frame()

# Get Lebanon incidence curve
leb_incidence_data <- arrow::read_feather("Data_incidence.feather")

## Load in true incidence curve to compare to our prediction
data(example_seir_incidence)
predictions <- plot_prob_infection(chain_comb, 
                                   nsamps=100, 
                                   INCIDENCE_FUNC=incidence_function,
                                   solve_times=0:max(ct_data_use$t),
                                   obs_dat=ct_data_use,
                                   true_prob_infection=leb_incidence_data)
p_incidence_prediction <- predictions$plot + scale_x_continuous(limits=c(0,200))
p_incidence_prediction
