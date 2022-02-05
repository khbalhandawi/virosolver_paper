########################################
## SCRIPT TO PERFORM SIM-RECOVERY OF MASSACHUSSETTS EPIDEMIC ON THE CLUSTER
## James Hay 6 November 2020
##  - Simulation-recovery of nursing homes with different population sizes and prior strengths
##  - First part of the script simulates an outbreak in a nursing home with specified population size and R0
##  - Then, simulates observed Ct values at the specified observation times. Here, use 3 obs times for growth, peak and decline phases
##  - Finally, re-fit the SEIR model to these simulations
##  - CLUSTER NOTE: this reads in pre-simulted Ct data and SEIR dynamics and fits the Ct 
##    model with specified settings to whatever timepoint was chosen

########################################
## 1. Headers
########################################
library(tidyverse)
library(ggthemes)
library(ggpubr)
library(data.table)
library(patchwork)
library(fitdistrplus)
library(deSolve)
library(lazymcmc) ## devtools::install_github("jameshay218/lazymcmc")
library(doParallel)
#HOME_WD <- "~/Documents/GitHub/"
HOME_WD <- Sys.getenv("GITDIR")

# devtools::load_all(paste0(HOME_WD,"/lazymcmc")) # parallel tempering branch
devtools::load_all(paste0(HOME_WD,"/virosolver"))

rerun_mcmc <- TRUE
use_pt <- FALSE
cap_patients <- 1000
p_patients <- 0.1
adaptive_period <- 100000
thin <- 50

start_date <- "2020-04-15"
start_plot <- "2020-04-15"
end_date <- "2021-04-01"
end_plot <- "2021-04-15"

# start_date <- "2020-06-01"
# start_plot <- "2020-04-15"
# end_date <- "2021-04-01"
# end_plot <- "2021-04-15"

## use parallel tempering
if (use_pt){
  devtools::load_all(paste0(HOME_WD,"/lazymcmc")) # parallel tempering branch
  
  ## MCMC parameters for Ct model fits if using parallel tempering branch
  n_temperatures <- 10
  # mcmcPars_ct <- list("iterations"=200000,"popt"=0.44,"opt_freq"=1000,
  #                     "thin"=100,"adaptive_period"=100000,"save_block"=100,
  #                     "temperature" = seq(1,101,length.out=n_temperatures),
  #                     "parallel_tempering_iter" = 5,"max_adaptive_period" = 100000, 
  #                     "adaptiveLeeway" = 0.2, "max_total_iterations" = 200000)
  # 
  mcmcPars_ct <- list("iterations"=50000,"popt"=0.44,"opt_freq"=1000,
                      "thin"=5,"adaptive_period"=50000,"save_block"=100,
                      "temperature" = seq(1,101,length.out=n_temperatures),
                      "parallel_tempering_iter" = 5,"max_adaptive_period" = 50000, 
                      "adaptiveLeeway" = 0.2, "max_total_iterations" = 50000)
} else {
  ## MCMC parameters for Ct model fits
  mcmcPars_ct <- c("iterations"=100000,"popt"=0.44,"opt_freq"=1000,
                   "thin"=thin,"adaptive_period"=adaptive_period,"save_block"=100)
}

## Arguments for this run, controlled by a csv file
control_table <- read_csv(paste0(HOME_WD,"/virosolver_paper/pars/lebanon/sim_leb_control_use.csv"))

## Priors for all models - EDIT THIS FILE TO CHANGE PRIORS!
source(paste0(HOME_WD,"/virosolver_paper/code/priors.R"))
source(paste0(HOME_WD,"/virosolver_paper/code/plot_funcs.R"))

# read command line arguments
args = commandArgs(trailingOnly=TRUE)

## Set random seed
simno <- strtoi(args[1])
set.seed(simno)

########################################
## 2. Read in run control parameters
########################################
## Name of this run
runname <- control_table$run_name[simno] 
## Which subfolder to save in? Usually different, but should remain unchanged if running multiple chains
run_index <- control_table$run_index[simno] 
## File for the Ct values
data_file_cts <- control_table$data_file_cts[simno] 
## Using all preceding times or just the most recent time?
cumu_times <- control_table$cumu_times[simno] 
## Using all Cts or just positives?
use_pos <- control_table$use_pos[simno] 
## Is the GP prior fixed or estimated? (if applicable)
gp_fixed <- control_table$gp_fixed[simno] 
## Which incidence model are we assuming?
model_version <- control_table$model_version[simno] 
## Where to store MCMC chains?
top_chainwd <- control_table$chainwd[simno] 
## Where to save output plots?
top_plotwd <- control_table$plotwd[simno] 
## In the vector of timepoints, which observation time are we at?
timepoint <- control_table$timepoint[simno] 
## If running repeats, what's the MCMC chain index?
chainno <- control_table$chainno[simno]
## Pointer to prior function to use
prior_func_use <- get(control_table$prior_func[simno])
## Pointer to incidence function to use
inc_func_use <- get(control_table$inc_func[simno])
## Maximum allowable age of infection?
age_max <- control_table$age_max[simno]

## Most importantly, where is the parameter control table stored?
parTab <- read.csv(control_table$parTab_file[simno], stringsAsFactors=FALSE)

########################################
## 3. Final admin
########################################
## Fixed control parameters
n_samp <- adaptive_period/thin

## CHANGE TO MAIN WD and manage save locations
chainwd <- paste0(top_chainwd,runname,"/",run_index,"/timepoint_",timepoint)
plot_wd <- paste0(top_plotwd,runname,"/",run_index)
setwd(HOME_WD)

if(!file.exists(chainwd)) dir.create(chainwd,recursive = TRUE)
if(!file.exists(plot_wd)) dir.create(plot_wd,recursive = TRUE)

########################################
## 4. Model parameters
########################################
## GP model parameters for fitting
pars <- parTab$values
names(pars) <- parTab$names

## Means for priors
means <- parTab$values
names(means) <- parTab$names

## Read in Ct data
obs_dat_all <- read_csv(data_file_cts) %>% rename(panther_Ct=Ct) %>%
  mutate(platform="Panther",first_pos=1) %>%
  mutate(id=1:n()) %>%
  filter(panther_Ct < 40)

obs_dat <-  obs_dat_all %>% 
  filter(platform=="Panther" &
           first_pos %in% c(1,0)) %>%
  filter(Date > start_date & Date < end_date) %>% ## After biased symptomatic sampling time
  rename(date=Date) %>%
  left_join(epi_calendar) %>%
  dplyr::select(first_day,  panther_Ct, id) %>%
  mutate(first_day = as.numeric(first_day)) %>%
  mutate(first_day = first_day - min(first_day) + 35) %>% ## Start 35 days before first sample
  arrange(first_day) %>%
  rename(t = first_day, ct=panther_Ct) %>%
  group_by(t) %>%
  filter(row_number() < cap_patients) %>% 
  do(top_n(.,as.integer(ceiling(p_patients*n())),id))

## Get possible observation times
obs_times <- unique(obs_dat$t)
n_days <- last(obs_dat$t)

## Cannot observe after the last observation time
timepoint <- min(timepoint, length(obs_times))

## Filter data for either all preceding time points or only this time point
if(cumu_times) {
  obs_dat_use <- obs_dat %>% filter(t <= obs_times[timepoint])
} else {
  obs_dat_use <- obs_dat %>% filter(t == obs_times[timepoint])
}
## If only using positive Cts
if(use_pos) {
  obs_dat_use <- obs_dat_use %>% filter(ct < pars["intercept"])
}

## Re-center observation times if we have a maximum age of infection
if(!is.na(age_max)){
  obs_dat_use <- obs_dat_use %>% mutate(t = t - min(t),t = t + max_age)
}

frac_days <- last(obs_dat_use$t) / last(obs_dat$t)

## Vectors of times/infection ages for simulation
ages <- 1:max(obs_dat_use$t)
times <- 0:max(obs_dat_use$t)

## This is for the GP version
mat <- matrix(rep(times, each=length(times)),ncol=length(times))
t_dist <- abs(apply(mat, 2, function(x) x-times))
if(model_version == "gp"){
  parTab <- bind_rows(parTab[parTab$names != "prob",], parTab[parTab$names == "prob",][1:length(times),])
  pars <- parTab$values
  names(pars) <- parTab$names
  
  ## Means for priors
  means <- parTab$values
  names(means) <- parTab$names
  
  if(gp_fixed == FALSE){
    parTab[parTab$names %in% c("nu","rho"),"fixed"] <- 0
  }
  
}

## Epidemic cannot start after first observation time
parTab[parTab$names == "t0",c("upper_bound","upper_start")] <- min(obs_dat_use$t)

## Test that model works with these settings
f <- create_posterior_func(parTab, obs_dat_use, prior_func_use, 
                           inc_func_use,solve_ver="likelihood",
                           use_pos=use_pos,
                           t_dist=t_dist)
f(pars)

########################################
## 5. Run MCMC
########################################
## Generate viable starting parameters
if(rerun_mcmc){
  startTab <- generate_viable_start_pars(parTab=parTab,
                                        obs_dat=obs_dat_use,
                                        CREATE_POSTERIOR_FUNC=create_posterior_func,
                                        INCIDENCE_FUNC=inc_func_use,
                                        PRIOR_FUNC=prior_func_use,
                                        use_pos=use_pos,
                                        t_dist=t_dist)
  
  covMat <- diag(nrow(startTab))
  mvrPars <- list(covMat,2.38/sqrt(nrow(startTab[startTab$fixed==0,])),w=0.8)
      
  output <- run_MCMC(parTab=startTab,
                      data=obs_dat_use,
                      INCIDENCE_FUNC=inc_func_use,
                      PRIOR_FUNC = prior_func_use,
                      solve_likelihood=TRUE,
                      mcmcPars=mcmcPars_ct,
                      filename=paste0(chainwd,"/",runname,"_",run_index,"_",timepoint,"_chainno_",chainno),
                      CREATE_POSTERIOR_FUNC=create_posterior_func,
                      mvrPars=NULL,
                      OPT_TUNING=0.2,
                      use_pos=use_pos,
                      t_dist=t_dist)
}
########################################
## 6. Check MCMC run
########################################    
## Read in chain and remove burn in period
chain <- read.csv(paste0(chainwd,"/",runname,"_",run_index,"_",timepoint,"_chainno_",chainno,"_univariate_chain.csv"))
chain <- chain[chain$sampno > mcmcPars_ct["adaptive_period"],]
chain_comb <- chain

## Get smoothed growth rates and raw incidence projections
samps <- sample(unique(chain_comb$sampno),n_samp)
trajs_smoothed <- matrix(0, nrow=n_samp,ncol=length(times))
trajs_normal <- matrix(0, nrow=n_samp,ncol=length(times))
for(ii in seq_along(samps)){
  ## Ensure these are all positive
  tmp_smoothed <- pmax(smooth.spline(inc_func_use(get_index_pars(chain_comb, samps[ii]),times))$y,0.0000001)
  tmp_normal <- pmax(inc_func_use(get_index_pars(chain_comb, samps[ii]),times),0.0000001)
  
  trajs_smoothed[ii,] <- (tmp_smoothed/sum(tmp_smoothed)) * frac_days
  trajs_normal[ii,] <- (tmp_normal/sum(tmp_normal)) * frac_days
}

## Get daily growth rate from smoothed version, find quantiles
trajs1 <- t(apply(trajs_smoothed, 1, function(x) log(x[2:length(x)]/x[1:(length(x)-1)])))
trajs1_quants <- t(apply(trajs1, 2, function(x) quantile(x,c(0.025,0.5,0.975))))
trajs1_quants <- as.data.frame(trajs1_quants)
trajs1_quants$t <- 1:nrow(trajs1_quants)
colnames(trajs1_quants) <- c("lower","median","upper","t")

## Growth rate plot
p_gr <- ggplot(trajs1_quants) + geom_ribbon(aes(x=t,ymin=lower,ymax=upper),alpha=0.25) + 
  geom_line(aes(x=t,y=median)) + 
  coord_cartesian(ylim=c(-0.5,0.5))

## Get daily incidence from unsmoothed version, find quantiles
trajs_quants <- t(apply(trajs_normal, 2, function(x) quantile(x,c(0.025,0.5,0.975))))
trajs_quants <- as.data.frame(trajs_quants)
trajs_quants$t <- 1:nrow(trajs_quants)
colnames(trajs_quants) <- c("lower","median","upper","t")

## Ct distribution plot
p_dat <- ggplot(obs_dat_use %>% filter(ct < pars["intercept"])) + 
  geom_violin(aes(x=t,group=t,y=ct),scale="width",fill="grey70",draw_quantiles=c(0.025,0.5,0.975)) + 
  scale_y_continuous(trans="reverse") +
  export_theme +
  scale_x_continuous(limits=c(0,n_days+5))

## Incidence rate plot
p_inc <- ggplot(trajs_quants) + geom_ribbon(aes(x=t,ymin=lower,ymax=upper),alpha=0.25) + 
  geom_line(aes(x=t,y=median)) + 
  geom_line(data=tibble(t=times,y=inc_func_use(get_best_pars(chain_comb),times)),aes(x=t,y=y),col="green") +
  #geom_line(data=tibble(t=1:,n_days,y=(seir_dynamics$incidence/population_n)[1:,n_days]),aes(x=t,y=y),col="red") +
  export_theme +
  ylab("Per capita incidence") +
  xlab("Days since start") +
  scale_x_continuous(limits=c(0,n_days+5)) +
  coord_cartesian(ylim=c(0,0.025))

## Get viral load trajectories
ages1 <- 1:50
vl_trajs <-  matrix(0, nrow=n_samp,ncol=length(ages1))
for(ii in 1:n_samp){
  tmp_pars <- get_index_pars(chain_comb, samps[ii])
  tmp <- viral_load_func(tmp_pars,ages1,FALSE)
  tmp1 <- extraDistr::rgumbel(length(tmp),tmp, tmp_pars["obs_sd"])
  vl_trajs[ii,] <- tmp1
}
vl_trajs1_quants <- t(apply(vl_trajs, 2, function(x) quantile(x,c(0.025,0.5,0.975))))
vl_trajs1_quants <- as.data.frame(vl_trajs1_quants)
vl_trajs1_quants$t <- 1:nrow(vl_trajs1_quants)
colnames(vl_trajs1_quants) <- c("lower","median","upper","t")

## Viral load trajectory plot
p_vl <- ggplot(vl_trajs1_quants) + geom_ribbon(aes(x=t,ymin=lower,ymax=upper),alpha=0.25) + 
  geom_line(aes(x=t,y=median))

## Trace plot
chain1 <- chain
if("prob" %in% colnames(chain)){
  use_cols <- c(which(colnames(chain) != "prob"), which(colnames(chain) == "prob")[1])
  chain1 <- chain[,use_cols]
}
chain_comb <- chain_comb[,colnames(chain_comb) != "chain"]

p_trace <- chain1[,c("sampno",unique(parTab[which(parTab$fixed == 0),"names"]))] %>%
  pivot_longer(-sampno) %>%
  ggplot() +
  geom_line(aes(x=sampno,y=value)) +
  facet_wrap(~name,scales="free_y")+
  scale_x_continuous(breaks=seq(min(chain$sampno),max(chain$sampno),by=max(chain$sampno)/5)) +
  export_theme

plotname_tmp <- paste0(plot_wd,"/",runname,"_",run_index,"_",timepoint,"_chainno_",chainno)
ggsave(paste0(plotname_tmp,"_trace.png"),p_trace,width=7,height=4,dpi=150)
ggsave(paste0(plotname_tmp,"_vl.png"),p_vl,width=7,height=4,dpi=150)
ggsave(paste0(plotname_tmp,"_growth_rate.png"),p_gr,width=7,height=4,dpi=150)
ggsave(paste0(plotname_tmp,"_fits.png"),p_dat/p_inc,width=7,height=6,dpi=150)

print("Job completed successfully")

