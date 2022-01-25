## IMPORTANT NOTE: the single cross section SEIR models are indexed from 1st February, so t[0] = 1st February there.
## The exponential and GP models are indexed from 35 days prior to the first sampling time. So all SEIR model times are 36
## days after the exp and GP models. HOWEVER, we index obs_dat based on the exponential and GP models. This is why
## you will see some times as 35, and others as 71 referencing the same day.

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

AAAS_palette <- c("blue1"="#3B4992FF","red1"="#EE0000FF","green1"="#008B45FF",
                  "purple1"="#631879FF","teal1"="#008280FF","red2"="#BB0021FF",
                  "purple2"="#5F559BFF","purple3"="#A20056FF",
                  "grey1"="#808180FF","black"="1B1919FF")

#HOME_WD <- "~/"
HOME_WD <- "C:/Users/Khalil/Desktop/repos"
devtools::load_all(paste0(HOME_WD,"/virosolver"))

## Where the MCMC chains are stored
main_wd <- paste0(HOME_WD,"/virosolver_paper/")

## Set flag to write plots to disk or not
save_plots <- FALSE

## CHANGE TO MAIN WD
setwd(main_wd)
source("code/plot_funcs.R")

export_theme <- export_theme + theme(plot.margin=unit(c(0,0,0,0),units="cm"),plot.tag=element_text(size=10,face="bold"))


## Arguments for this run
set.seed(1)
n_samp <- 1000

chainwd_gp <- paste0(HOME_WD,"/virosolver_paper/mcmc_chains/4.real_ma_ct/RHUH_gp/")

## MCMC parameters for Ct model fits
mcmcPars_ct_exp <- c("adaptive_period"=30000)

## Create an epiweek calendar
dates <- seq(as.Date("2020-01-01"),as.Date("2020-12-31"),by="1 day")
epiweeks <- lubridate::epiweek(dates)
epi_calendar <- tibble(date=dates,week=epiweeks)
epi_calendar <- epi_calendar %>% group_by(week) %>% mutate(first_day=min(date))

########################################
## 3. MA incidence plot and Rt
########################################
## LEB data
leb_dat <- read_csv("data/RHUH_cases_data.csv")

## Rt fit
#estimates <- readRDS("~/Documents/GitHub/ct_dynamics_preprint/results/ma_rt_fit.RData")
estimates <- readRDS("results/RHUH_rt_fit.RData")
rt_dat <- estimates$estimates$summarised %>% filter(variable %in% c("growth_rate") & type=="estimate")

## Split based on above or below Rt = 1
rt_dat_top <- rt_dat %>%
  mutate(bottom=pmax(0,lower_90),
         top=pmax(0,upper_90))
rt_dat_bot <- rt_dat %>%
  mutate(bottom=pmin(0,lower_90),
         top=pmin(0,upper_90))
## Same thing for median line
rt_median <- rt_dat %>% mutate(is_grow=ifelse(median>0,"Growing","Declining"))
index <- 1
rt_median$index <- index
for(i in 2:nrow(rt_median)){
  if(rt_median$is_grow[i] != rt_median$is_grow[i-1]) {
    index <- index + 1
  }
  rt_median$index[i] <- index
}
coeff <- 2500
yshift <- 500
p1 <-  ggplot(leb_dat)+ 
  geom_bar(aes(x=date,y=new_cases),stat="identity",alpha=0.5,fill=AAAS_palette["grey1"]) +
  #geom_ribbon(data=dat_infections,aes(x=date,ymin=lower_90,ymax=upper_90),
  #            fill=AAAS_palette["grey1"],alpha=0.5) +
  #geom_line(data=dat_infections,aes(x=date,y=mean),col=AAAS_palette["grey1"]) +
  geom_hline(yintercept=yshift,linetype="dashed",col=AAAS_palette["grey1"]) +
  geom_ribbon(data=rt_dat_top,aes(x=date,ymin=lower_90*coeff + yshift,ymax=upper_90*coeff + yshift),
              fill="darkorange",alpha=0.25) +
  geom_ribbon(data=rt_dat_bot,aes(x=date,ymin=lower_90*coeff + yshift,ymax=upper_90*coeff + yshift),
              fill=AAAS_palette["green1"],alpha=0.25) +
  geom_line(data=rt_median,aes(x=date,y=median*coeff + yshift,col=is_grow,group=index)) +
  scale_y_continuous(expand=c(0,0),limits=c(0,2000),
                     sec.axis=sec_axis(~.*1/coeff - yshift/coeff, name="Growth rate")) +
  scale_color_manual(values=c("Growing"="darkorange","Declining"=as.character(AAAS_palette["green1"])))+
  scale_x_date(limits=as.Date(c("2020-03-01", "2020-12-01"), "%Y-%m-%d"), breaks="1 month",
               expand=c(0,0)) +
  export_theme+
  theme(legend.position="none",
        panel.grid.minor=element_blank(),
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.line.x=element_blank(),
        axis.ticks.x=element_blank()) +
  xlab("Date") +
  ylab("New cases") +
  labs(tag="A")
p1



########################################
## 4. RHUH Ct data
########################################
obs_dat_all <- read_csv(paste0(main_wd,"/data/RHUH_Ct_data.csv")) %>% rename(panther_Ct=Ct) %>%
  mutate(platform="Panther",first_pos=1) %>%
  mutate(id=1:n()) %>%
  drop_na() %>%
  filter(panther_Ct < 40) %>%
  filter(Date > "2020-04-15")

obs_dat1 <-  obs_dat_all %>% 
  filter(platform=="Panther" &
           first_pos %in% c(1,0)) %>%
  filter(Date > "2020-04-15") %>% ## After biased symptomatic sampling time
  rename(date=Date) %>%
  left_join(epi_calendar) %>%
  dplyr::select(first_day,  panther_Ct, id) %>%
  mutate(first_day = as.numeric(first_day)) %>%
  mutate(first_day = first_day - min(first_day) + 35) %>% ## Start 35 days before first sample
  #mutate(earliest_date=as.Date("2020-02-01")) %>%
  #mutate(first_day = as.numeric(first_day - earliest_date)) %>% ## Assume 1st February is start date
  arrange(first_day) %>%
  rename(t = first_day, ct=panther_Ct)

obs_dat_all <- obs_dat_all %>% 
  filter(platform=="Panther" &
           first_pos %in% c(1,0)) %>%
  filter(Date > "2020-04-15") %>% ## After biased symptomatic sampling time
  rename(date=Date) %>%
  left_join(epi_calendar) %>%
  dplyr::select(first_day,  panther_Ct, id) %>%
  arrange(first_day) %>%
  rename(date = first_day, ct=panther_Ct)

comb_dat <- left_join(obs_dat1, obs_dat_all)
date_key <- distinct(comb_dat %>% dplyr::select(t, date))

date_min_date <- min(date_key$date)
date_min_t <- min(date_key$t)

date_max_date <- max(date_key$date)
date_max_t <- max(date_key$t)

integer_seq_times <- seq(0, date_max_t)
date_seq_times <- seq(date_min_date-date_min_t, date_max_date,by="1 day")
date_key <- tibble(t=integer_seq_times,date=date_seq_times)

p_dat <- ggplot(obs_dat_all) + 
  geom_violin(aes(x=date,group=date,y=ct),scale="width",fill=AAAS_palette["grey1"],
              alpha=0.5,color="black",size=0.1,
              draw_quantiles=c(0.025,0.5,0.975)) + 
  geom_dotplot(aes(x=date, y=ct,group=date),binaxis="y",
               binwidth=1,stackdir="center",binpositions="all",dotsize=0.1) +
  geom_smooth(data=obs_dat_all %>% group_by(date) %>% summarize(median_ct=median(ct)),
              aes(x=date,y=median_ct),col=AAAS_palette["blue1"],se=FALSE) +
  scale_y_continuous(trans="reverse",limits=c(42, 10),expand=c(0,0)) +
  geom_hline(yintercept=40,linetype="dashed") +
  export_theme +
  theme(axis.text.x=element_blank(), axis.title.x=element_blank(),
        axis.line.x = element_blank(), axis.ticks.x = element_blank()) +
  scale_x_date(limits=as.Date(c("2020-04-15","2020-12-15")),breaks="1 month",expand=c(0,0)) +
  xlab("Date of sample") +
  ylab("Ct value") +
  labs(tag="C")
p_dat

########################################
## 5. Rt scatterplot
########################################
## Plot Ct data
rhuh_data <- read_csv(paste0(main_wd,"/data/RHUH_Ct_data.csv")) %>% rename(panther_Ct=Ct) %>%
  mutate(platform="Panther",first_pos=1) %>%
  mutate(id=1:n())

rhuh_data_use <-  rhuh_data %>% 
  filter(platform=="Panther" &
           first_pos %in% c(1,0)) %>%
  rename(date=Date) %>%
  left_join(epi_calendar)

rhuh_data_week <- rhuh_data_use %>% 
  group_by(date) %>% 
  summarize(median_ct=median(panther_Ct),
            skew_ct=moments::skewness(panther_Ct),
            n=n())

combined_dat <- rhuh_data_week %>% 
  full_join(estimates$estimates$summarised %>% as_tibble() %>% 
              filter(variable == "R") %>%
              mutate(date = date)
  ) %>%
  arrange(date) %>% 
  rename(mean_rt=mean) %>%
  dplyr::select(date, median_ct,skew_ct, n, mean_rt) %>%
  filter(n >= 10)

combined_dat1 <- combined_dat %>% 
  pivot_longer(-c(date,n)) %>%
  left_join(epi_calendar) %>%
  group_by(first_day, name) %>% 
  summarize(mean=mean(value)) %>%
  pivot_wider(values_from=mean,names_from=name)
combined_dat1 <- combined_dat

median_cts <- smooth(combined_dat$median_ct)
skew_cts <- combined_dat$skew_ct
Rt <- combined_dat$mean_rt
ccf(median_cts, Rt,lag.max = 25)
ccf(skew_cts, Rt,lag.max = 25)

p_rt <- ggplot(combined_dat1 %>% rename(Rt=mean_rt)) +
  geom_point(aes(x=skew_ct,y=median_ct,col=Rt),alpha=0.9,size=2) +
  scale_color_gradient2(low="green",mid="blue",high="red",midpoint=1)+
  scale_y_continuous(trans="reverse") +
  xlab("Skewness of Ct distribution") +
  ylab("Median of Ct distribution") +
  export_theme + 
  theme(legend.position=c(0.2,0.7)) + 
  labs(tag="B")

########################################
## 7. GP fit
########################################
chains <- lazymcmc::load_mcmc_chains(chainwd_gp, parTab,FALSE,1,mcmcPars_ct_exp["adaptive_period"],
                                     multi=FALSE,chainNo=TRUE,PTchain = FALSE)
chain <- as.data.frame(chains$chain)
chain$sampno <- 1:nrow(chain)
chain_comb <- chain
chain_comb$sampno <- 1:nrow(chain_comb)

times <- 0:max(obs_dat1$t)
## Get smoothed growth rates
samps <- sample(unique(chain_comb$sampno),n_samp)
gp_trajs <- matrix(0, nrow=n_samp,ncol=length(times))
for(ii in seq_along(samps)){
  tmp <- pmax(smooth.spline(gaussian_process_model(get_index_pars(chain_comb, samps[ii]),times))$y,0)
  gp_trajs[ii,] <- tmp/sum(tmp)
  #trajs[ii,] <- pmax(inc_func_use(get_index_pars(chain_comb, samps[ii]),times),0.0000001)
}

gp_trajs1 <- t(apply(gp_trajs, 1, function(x) log(x[2:length(x)]/x[1:(length(x)-1)])))
gp_trajs1_quants <- t(apply(gp_trajs1, 2, function(x) quantile(x,c(0.025,0.5,0.975),na.rm=TRUE)))
gp_trajs1_quants <- as.data.frame(gp_trajs1_quants)
gp_trajs1_quants$t <- 1:nrow(gp_trajs1_quants)
colnames(gp_trajs1_quants) <- c("lower","median","upper","t")

## Growth rate plot
p_gr <- ggplot(trajs1_quants) + geom_ribbon(aes(x=t,ymin=lower,ymax=upper),alpha=0.25) + 
  geom_line(aes(x=t,y=median)) + 
  coord_cartesian(ylim=c(-0.5,0.5))


gp_trajs_quants <- t(apply(gp_trajs, 2, function(x) quantile(x,c(0.025,0.25,0.5,0.75,0.975),na.rm=TRUE)))
gp_trajs_quants <- as.data.frame(gp_trajs_quants)
gp_trajs_quants$t <- 1:nrow(gp_trajs_quants)
colnames(gp_trajs_quants) <- c("lower","mid_lower","median","mid_upper","upper","t")

## Growth rate plot
p_inc <- ggplot(gp_trajs_quants %>% left_join(date_key)) + 
  geom_ribbon(aes(x=date,ymin=lower,ymax=upper),alpha=0.25,fill=AAAS_palette["red1"]) + 
  geom_ribbon(aes(x=date,ymin=mid_lower,ymax=mid_upper),alpha=0.5,fill=AAAS_palette["red1"]) + 
  geom_line(aes(x=date,y=median),col=AAAS_palette["red1"]) + 
  #geom_line(data=tibble(t=times,y=inc_func_use(get_best_pars(chain_comb),times)),aes(x=t,y=y),col="green") +
  #geom_line(data=tibble(t=1:200,y=(seir_dynamics$incidence/population_n)[1:200]),aes(x=t,y=y),col="red") +
  export_theme +
  ylab("Scaled probability of infection") +
  xlab("Date") +
  scale_x_date(limits=as.Date(c("2020-04-01","2020-11-15")),breaks="1 month",expand=c(0,0)) +
  coord_cartesian(ylim=c(-0.0001,0.025)) +
  scale_y_continuous(expand=c(0,0)) +
  labs(tag="E")


########################################
## 8. Growth rate comparison
########################################
gp_trajs1_quants_joined <- left_join(gp_trajs1_quants, date_key)
gr_dat <- estimates$estimates$summarised %>% filter(variable == "growth_rate") %>% 
  filter(type=="estimate") %>%
  dplyr::select(date, mean,upper_90,lower_90) %>%
  rename(gr_mean=mean,
         gr_lower=lower_90,
         gr_upper=upper_90)

gr_dat_top <- gr_dat %>%
  mutate(gr_lower=pmax(0,gr_lower),
         gr_upper=pmax(0,gr_upper))
gr_dat_bot <- gr_dat %>%
  mutate(gr_lower=pmin(0,gr_lower),
         gr_upper=pmin(0,gr_upper))

gr_dat_median <- gr_dat %>% mutate(is_grow=ifelse(gr_mean>0,"R(t), growing","R(t), declining"))
index <- 1
gr_dat_median$index <- index
for(i in 2:nrow(gr_dat_median)){
  if(gr_dat_median$is_grow[i] != gr_dat_median$is_grow[i-1]) {
    index <- index + 1
  }
  gr_dat_median$index[i] <- index
}


trajs_dat_top <- gp_trajs1_quants_joined %>%
  mutate(lower=pmax(0,lower),
         upper=pmax(0,upper))
trajs_dat_bot <- gp_trajs1_quants_joined %>%
  mutate(lower=pmin(0,lower),
         upper=pmin(0,upper))

trajs_dat_median <- gp_trajs1_quants_joined %>% mutate(is_grow=ifelse(median>0,"Ct estimate","Ct estimate"))
index <- 1
trajs_dat_median$index <- index
for(i in 2:nrow(trajs_dat_median)){
  if(trajs_dat_median$is_grow[i] != trajs_dat_median$is_grow[i-1]) {
    index <- index + 1
  }
  trajs_dat_median$index[i] <- index
}

p_grs <- ggplot() + 
  geom_ribbon(data=gr_dat_top,aes(x=date,ymin=gr_lower,ymax=gr_upper),alpha=0.25,fill="darkorange") +
  geom_ribbon(data=gr_dat_bot,aes(x=date,ymin=gr_lower,ymax=gr_upper),alpha=0.25,fill=AAAS_palette["green1"]) +
  geom_line(data=gr_dat_median,aes(x=date,y=gr_mean,col=is_grow,group=index)) +
  geom_ribbon(data=trajs_dat_top,aes(x=date,ymin=lower,ymax=upper),alpha=0.25,fill=AAAS_palette["purple1"]) +
  geom_ribbon(data=trajs_dat_bot,aes(x=date,ymin=lower,ymax=upper),alpha=0.25,fill=AAAS_palette["purple1"]) +
  geom_line(data=trajs_dat_median,aes(x=date,y=median,col=is_grow,group=index)) +
  coord_cartesian(ylim=c(-0.125,0.25)) +
  scale_x_date(expand=c(0,0),limits=as.Date(c("2020-03-09", "2020-12-01"), "%Y-%m-%d"), 
               breaks=as.Date(c("2020-04-01","2020-07-01","2020-10-01"))) +
  geom_hline(yintercept=0,linetype="dashed",col=AAAS_palette["grey1"]) +
  scale_color_manual(values=c("R(t), growing"="darkorange",
                              "R(t), declining"=as.character(AAAS_palette["green1"]),
                              "Ct estimate"=as.character(AAAS_palette["purple1"])))+
  export_theme +
  theme(legend.position=c(0.5,0.9),
        legend.title = element_blank(),
        axis.text.x=element_text(size=6),
        legend.text=element_text(size=6),
        plot.margin=unit(c(0,0,0,0),units="cm")) +
  xlab("Date") +
  ylab("Growth rate") +
  labs(tag="F")

########################################
## 9. Pull together
########################################
if(TRUE){
  pdf("figures/Figure5.pdf",height=6,width=9)
  (((p1/p_dat/((ps_all_daily | plot_spacer()))/p_inc) + plot_layout(heights=c(4,3.5,0.5,4))) | (((plot_spacer()|p_rt)+plot_layout(widths=c(1,50)))/p_use_dial/p_grs)) + plot_layout(widths=c(2,1))
  dev.off()
  
  png("figures/Figure5.png",height=6,width=9,units="in",res=300)
  (((p1/p_dat/((ps_all_daily | plot_spacer()))/p_inc) + plot_layout(heights=c(4,3.5,0.5,4))) | (((plot_spacer()|p_rt)+plot_layout(widths=c(1,50)))/p_use_dial/p_grs)) + plot_layout(widths=c(2,1))
  dev.off()
}
if(save_plots){
  pdf("figures/Figure5.pdf",height=6,width=9)
  ((p1/p_dat/p_inc)| (((plot_spacer()|p_rt)+plot_layout(widths=c(1,50)))/p_use_dial/p_grs)) + plot_layout(widths=c(2,1))
  dev.off()
  
  png("figures/Figure5.png",height=6,width=9,units="in",res=300)
  ((p1/p_dat/p_inc)| (((plot_spacer()|p_rt)+plot_layout(widths=c(1,50)))/p_use_dial/p_grs)) + plot_layout(widths=c(2,1))
  dev.off()
}