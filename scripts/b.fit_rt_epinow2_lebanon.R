## Pre-script 2: Rt estimations based on MA case data
library(tidyverse)
library(patchwork)
library(EpiEstim)
library(zoo)
library(ggsci)
library(EpiNow2)
#HOME_WD <- "~/"
HOME_WD <- "C:/Users/Khalil/Desktop/repos"

## Where the MCMC chains are stored
main_wd <- paste0(HOME_WD,"/virosolver_paper/")

setwd(main_wd)

rerun <- TRUE

## Create an epiweek calendar
dates <- seq(as.Date("2020-01-01"),as.Date("2020-12-31"),by="1 day")
epiweeks <- lubridate::epiweek(dates)
epi_calendar <- tibble(date=dates,week=epiweeks)
epi_calendar <- epi_calendar %>% group_by(week) %>% mutate(first_day=min(date))

## LEB data
leb_dat <- read_csv("data/RHUH_cases_data.csv")

## Plot case counts
p1 <- leb_dat %>% 
  ## Highlight anomalous cases
  #mutate(new_cases1=ifelse(new_cases > 3500, lag(new_cases,1), new_cases)) %>%
  ## Get rolling mean
  mutate(roll_mean=rollmean(new_cases, 7,fill=NA,align="right")) %>%
  ggplot()+ 
  geom_bar(aes(x=date,y=new_cases),stat="identity",alpha=0.25) +
  #geom_ribbon(aes(x=date,ymax=roll_mean,ymin=0),fill="grey70",alpha=0.5,col="grey20") +
  #scale_fill_npg() +
  scale_fill_manual(values=c("grey40","darkred")) +
  scale_y_continuous(expand=c(0,0),limits=c(0,2000)) +
  scale_x_date(limits=as.Date(c("2020-04-01", "2020-12-01"), "%Y-%m-%d"), breaks="7 days",
               expand=c(0,0)) +
  theme_classic()+
  theme(legend.position="none",
        panel.grid.minor=element_blank(),
        #axis.text.x=element_text(angle=45,hjust=1),
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.ticks.x=element_blank()) +
  xlab("Date") +
  ylab("New infections") +
  labs(tag="A")
p1

## Rt estimation on these case counts
rt_dat <- leb_dat %>% 
  ungroup() %>%
  dplyr::select(date, new_cases) %>%
  drop_na() %>%
  rename(I=new_cases,
         dates=date)

reporting_delay <- EpiNow2::bootstrapped_dist_fit(rlnorm(100, log(3), 0.5))
## Set max allowed delay to 30 days to truncate computation
reporting_delay$max <- 30

generation_time <- get_generation_time(disease = "SARS-CoV-2", source = "ganyani")
incubation_period <- get_incubation_period(disease = "SARS-CoV-2", source = "lauer")

reported_cases <- EpiNow2::example_confirmed[1:50]
reported_cases <- rt_dat %>% rename(confirm=I, date=dates)

if(rerun){
  
  estimates <- epinow(reported_cases = reported_cases,
                      generation_time = generation_time,
                      delays = delay_opts(incubation_period, reporting_delay),
                      stan = stan_opts(samples = 4000, warmup = 1000, cores = 4, chains = 4, verbose = TRUE,
                      control = list(adapt_delta = 0.95)),
                      horizon = 7)

  saveRDS(estimates,"results/RHUH_rt_fit.RData")
}
  