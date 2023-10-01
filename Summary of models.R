rm(list=ls())
setwd("~/Dropbox (Personal)/RA/COVID Rt/Model Summary")
library(tidyverse)
library(EpiEstim)
library(tidyverse)
library(data.table)
library(doMC)
library(ggpubr)

#### function set up stochastic SIR model
SIR_model <- function(
  time_steps, # number of total timesteps
  pop.size, # population size
  seeds, # number of initial infections
  recovered, # number of recovered
  trans_prob, # probability of transmission given contact: beta
  inf_days, # average days of infectiousness*
  T0, # time to intervention
  trans_prob2 # trans_prob after intervention
  ){
  
  beta = trans_prob  

  # track states
  S = rep(0, time_steps)
  I = rep(0, time_steps)
  inc = rep(0, time_steps)
  R = rep(0, time_steps)
  C = rep(0, time_steps)
  mean = rep(0, time_steps)
  
  # initial conditions
  S[1] = pop.size - seeds - recovered
  I[1] = seeds
  R[1] = recovered
  
  for(i in 2:time_steps){
    # change beta at the change time point
    if(i == T0 + 1) beta = trans_prob2
    
    # set up random draw
    mean[i] = beta*I[i-1]/pop.size*S[i-1]
    trans_t = rpois(1, lambda = mean[i])
    
    # susceptible
    S[i] = S[i-1] - trans_t
    
    # infectious
    I[i] = (1-1/inf_days)*I[i-1] + trans_t
    inc[i] = trans_t
    C[i] = C[i-1] + trans_t
    
    # recovered
    R[i] = R[i-1] + 1/inf_days*I[i-1]
  }
  
  # return data table
  d = data.table(pop.size, S, I, R, inc, C, mean) %>% mutate(t = row_number())
  return(d)
}

run_model <- function(pop.size1, pop.size2, N = 2, N1 = 1, T.total, T1, seed1 = 50, seed2 = 50, R1 = 0, R2 = 0,
                      burnin, trans_prob.base1, trans_prob.base2, inf_days, eff.add1 = 0, eff.multi1 = 1,
                      eff.multi2 = 1, eff.add2 = 0, model_fxn = SIR_model) {
  T0 <- T.total + burnin - T1
  trt.IDs <- 1:N1
  
  out <- lapply(1:N, function(ind) { 
    if (ind %in% trt.IDs) {
      model_fxn(pop.size=pop.size1, time_steps=(T.total+burnin), trans_prob=trans_prob.base1, 
                T0 = T0, trans_prob2 = trans_prob.base1*eff.multi1 + eff.add1, 
                inf_days = inf_days, seeds = seed1, recovered = R1)
    } else {
      model_fxn(pop.size=pop.size2, time_steps=(T.total+burnin), trans_prob=trans_prob.base2, 
                T0 = T0, trans_prob2=trans_prob.base2*eff.multi2 + eff.add2, 
                inf_days = inf_days, seeds = seed2, recovered = R2)
    }})
  out.df <- rbindlist(out)
  
  df <- out.df %>% 
    mutate(unit = rep(1:N, each=(T.total+burnin)),
           trt.unit = unit %in% trt.IDs,
           trt.time = t > T0,
           trt_post = trt.unit & trt.time)
  
  return(df)
}
##############################################################################################################################
# set up variable inputs
T.total <- 80
burnin <- 25
beta <- 0.105 # baseline beta
beta2 <- 0.20 # to change for different transmission
beta3 <- 0.16 # to change for non-trivial susceptible depletion
inf_days <- 10
pop <- 1e+08 # make it large enough so there is NO susceptible depletion in other cases
pop2 <- 1e+04 # to force a non-trivial susceptible depletion where needed
agg <- 2 # days of aggregation / smoothing for log growth
##############################################################################################################################
# parallelization run test scenarios
doMC::registerDoMC(cores = detectCores()-5)
foreach::getDoParWorkers()
test1.tmp <- foreach(k=1:500, .errorhandling = "pass", .combine=function(x,y) rbindlist(list(x,y))) %dopar% {
  run_model(pop.size1=pop, pop.size2=pop, T.total=T.total, burnin=burnin, T1=0, inf_days=inf_days,
            trans_prob.base1=beta, trans_prob.base2=beta)
}
test2.tmp <- foreach(k=1:500, .errorhandling = "pass", .combine=function(x,y) rbindlist(list(x,y))) %dopar% {
  run_model(seed1=200, pop.size1=pop, pop.size2=pop, T.total=T.total, burnin=burnin, T1=0, inf_days=inf_days,
            trans_prob.base1=beta, trans_prob.base2=beta)
}
test3.tmp <- foreach(k=1:500, .errorhandling = "pass", .combine=function(x,y) rbindlist(list(x,y))) %dopar% {
  run_model(R1=pop*0.05, pop.size1=pop, pop.size2=pop, T.total=T.total, burnin=burnin, T1=0, inf_days=inf_days,
            trans_prob.base1=beta, trans_prob.base2=beta)
}
test4.tmp <- foreach(k=1:3000, .errorhandling = "pass", .combine=function(x,y) rbindlist(list(x,y))) %dopar% {
  run_model(trans_prob.base1=beta2, pop.size1=pop, pop.size2=pop, T.total=T.total, burnin=burnin, 
            T1=30,
            inf_days=inf_days, trans_prob.base2=beta,
            eff.multi1 = 2, eff.multi2 = 2)
}
test5.tmp <- foreach(k=1:3000, .errorhandling = "pass", .combine=function(x,y) rbindlist(list(x,y))) %dopar% {
  run_model(pop.size1=pop2, pop.size2=pop, T.total=T.total, burnin=burnin, T1=0, inf_days=inf_days,
            trans_prob.base1=beta3, trans_prob.base2=beta)
}
##############################################################################################################################
# Summarize / Average over simulations
test1 <- test1.tmp %>%
  group_by(t, unit, pop.size) %>%
  reframe(S = mean(S), I = mean(I), R = mean(R),
          inc = mean(inc), C = mean(C), mean = mean(mean),
          id = "'(a) All parameters match'", 
          lab = "")
test2 <- test2.tmp %>%
  group_by(t, unit, pop.size) %>%
  reframe(S = mean(S), I = mean(I), R = mean(R),
          inc = mean(inc), C = mean(C), mean = mean(mean),
          id = "'(b) Different initial conditions'", lab = "I[0]~`=`~200")
test3 <- test3.tmp %>%
  group_by(t, unit, pop.size) %>%
  reframe(S = mean(S), I = mean(I), R = mean(R),
          inc = mean(inc), C = mean(C), mean = mean(mean),
          id = "'(c) Different susceptible\npopulations'", lab = "S[0]/N~`=`~0.95")
test4 <- test4.tmp %>%
  group_by(t, unit, pop.size) %>%
  reframe(S = mean(S), I = mean(I), R = mean(R),
          inc = mean(inc), C = mean(C), mean = mean(mean),
          id = "'(d) Different transmission'", lab = 'beta[t]~"doubled in both groups at time 45"')
test5 <- test5.tmp %>%
  group_by(t, unit, pop.size) %>%
  reframe(S = mean(S), I = mean(I), R = mean(R), inc = mean(inc), 
          C = mean(C), mean = mean(mean),
          id = "'(e) Non-trivial susceptible\ndepletion'", 
          lab = paste0("N*`=`*10*`,`*0*0*0~`&`~beta[t]~`=`~", beta3))
##############################################################################################################################
# calculate the outcomes in each model
inc_df <- rbind(test1, test2, test3, test4, test5) %>%  
  mutate(y = inc, 
         y = ifelse(id=="'(d) Different transmission'", y/1500, y), # scale down; to make figure look better; otherwise all other lines are shown as a flat line at 0
         y = ifelse(id=="'(e) Non-trivial susceptible\ndepletion'", y/3, y),
         model = "incidence") %>%
  dplyr::select(model, unit, id, lab, y, t)

loginc_df <- rbind(test1, test2, test3, test4, test5) %>% 
  mutate(y = ifelse(inc > 0, log(inc), 0),
         y = ifelse(id=="'(d) Different transmission'", y/2, y),
         model = "'log incidence'") %>%
  dplyr::select(model, unit, id, lab, y, t)

growth_df <- rbind(test1, test2, test3, test4, test5) %>%
  mutate(week = ceiling(t/agg)) %>%
  group_by(unit, id, lab, week) %>%
  summarize(inc = sum(inc)) %>% 
  mutate(y = log(inc) - log(lag(inc, 1)),
         y = ifelse(id=="'(d) Different transmission'", y/2, y),
         t = week*agg,
         model = "'log growth'") %>%
  dplyr::select(model, unit, id, lab, y, t)

Rt_df <- rbind(test1, test2, test3, test4, test5) %>%
  group_by(id, lab, unit) %>%
  reframe(t = 3:(T.total+burnin),
          y = wallinga_teunis(round(inc), method = "parametric_si",
                              config = list(
                                t_start = t-1,
                                t_end = t,
                                method = "parametric_si", 
                                mean_si = inf_days, std_si = inf_days,
                                n_sim = 0))$R$`Mean(R)`) %>%
  mutate(y = ifelse(id=="'(d) Different transmission'", y/1.5, y),
         y = log(y),
         model = "log~R[t]") %>%
  dplyr::select(model, unit, id, lab, y, t)

beta_df <- rbind(test1, test2, test3, test4, test5) %>%
  mutate(sus_frac = S/pop.size) %>%
  merge(Rt_df, by = c("id", "lab", "unit", "t")) %>%
  mutate(y = exp(y)/sus_frac,
         y = log(y),
         model = "log~\u03B2[t]") %>%
  dplyr::select(model, unit, id, lab, y, t)
##############################################################################################################################
# merge all model data together for graphing
p.df <- inc_df %>% rbind(loginc_df) %>% rbind(beta_df) %>% rbind(growth_df) %>% rbind(Rt_df) %>% 
  mutate(model = factor(model, levels = c("incidence", "'log incidence'", "'log growth'", 
                                          "log~R[t]", "log~\u03B2[t]"))) %>%
  filter(t <= T.total,
         t > (burnin + inf_days*0.5)) %>%
  mutate(t = t - (burnin + inf_days*0.5)) %>%
  group_by(id, model, unit) %>%
  mutate(first.time = ifelse(t==min(t), 1, 0))
##############################################################################################################################
ggplot(p.df, aes(x = t, y = y, group = unit, col = factor(unit))) + geom_line() + 
  facet_grid(model~id, scales = "free_y", labeller = label_parsed) + 
  scale_y_continuous(expand = c(0.2, 0.2), n.breaks = 4) +
  geom_text(data = p.df %>% filter(unit==1 & first.time==1), size = 4,
            aes(label = lab,  x = 0, y = Inf), parse = T, hjust = 0, vjust = 1.2) +
  labs(x = "", y = "") + 
  scale_color_brewer(guide = "none", name = "", palette = "Set1") + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.text = element_text(size = 14))
ggsave("Figures/Summary of Models.png", width = 15, height = 10)
