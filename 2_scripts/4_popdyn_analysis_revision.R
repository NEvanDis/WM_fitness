# Analysis of Field data: population numbers ####
# How much variation in population growth be explained by timing asynchrony?

# Geert Jan (GJ) Sieperda explored models at Area level and at Site/Tree level
# Final model that we settled on was at Area level, so that's what I will use here
# Similarly GJ tried both GAM and lm models, and we settled on lm models
# Scripts GJ not included in deposition


# Open R project in main folder

# Load packages
#-----------------------------------
library(tidyverse)
library(cowplot)
library(tseries) #kpss test
library(forecast) # AR models
library(lme4)
library(lmerTest)
library(rms) # for backward elimination
library(effectsize) # for omega variance measure
theme_set(theme_cowplot()) #white background instead of grey -> don't load if want grey grid


# Load data ####
#-----------------------------------
# All data = per year per Area (Site and Tree info not included)

d_popnr <- read.csv("1_data/1994-2021_Field_PopNum_cor.csv") # produced by script <prep_FieldData.R>
d_popnr <- d_popnr %>% mutate(FemLogAbund=log(NoFemCor + 1),# correct for zero values
                              Area=ifelse(AreaShortName=="DO", "Doorwerth", 
                                          ifelse(AreaShortName=="HV", "Hoge Veluwe", 
                                                 ifelse(AreaShortName=="OH", "Oosterhout", "Warnsborn")))) 
head(d_popnr) # Population numbers corrected for trapping effort, grp=for plotting data gaps correctly

popnr_list <- d_popnr %>% group_split(AreaShortName) # store as list
names(popnr_list) <- c("DO", "HV", "OH", "WA")

# Data gaps ####
lapply(popnr_list, function(X) table(X$grp))
# Hoge Veluwe has a large data gap

# Nr of years per area ####
lapply(popnr_list, nrow)

# What Mismatch measure want to use? ####
# --First calculated mismatch per tree from individuals, then calculated Mismatch.mean per area from that + Mismatch.sd

d_mismatch <- read.csv("1_data/1994-2021_rel_timing_werror.csv") # produced by script <prep_FieldData.R>
head(d_mismatch) 

# Eggs
d_eggs <- read_csv("1_data/1994-2021_Field_Eggs_deposit.csv")
head(d_eggs)

# Want mean Eggs per female per year per Area
Eggs <-aggregate(Eggs~ YearCatch + AreaName, data=d_eggs, mean)
Eggs$Lag1 <- Eggs$YearCatch+1 # explanatory variable = mean number of eggs of YearX-1 that could explain population growth in YearX
head(Eggs)


#-----------------------------------
# Determine cyclicity ####
#-----------------------------------

# Spectral analysis
#-----------------------------------
# Example: https://web.stanford.edu/class/earthsys214/notes/series.html
?spectrum

l_spectr <- lapply(popnr_list, function(X) {
  
  if(X$AreaShortName[1]=="HV"){
    X <- filter(X, grp==2)
  }# HV has a data gap, disregard first two years
  
  spect <- spectrum(X$FemLogAbund, log='no', plot=FALSE) # log transform values beforehand to deal with zeros
  
  Area <- unique(X$Area)
  
  d_spec <- data.frame(spec_x=1/spect$freq, #Freq = calculated in terms of cycles per year, so 0.1 cycles per year = 1 cycle every 10 years
                       spec_y=spect$spec*2, #multiply spectral density by 2 so that area under periodogram equals variance of time series
                       Area=Area)
  
  d_spec$peak <- d_spec[d_spec$spec_y==max(d_spec$spec_y), "spec_x"]
  
  d_spec # return
})
lapply(l_spectr, head)

res_spect <- data.table::rbindlist(l_spectr, use.names=TRUE)
head(res_spect)

# Results ####
periods <- res_spect[!duplicated(res_spect[,c("Area", "peak")]),] %>% select(Area, peak)
# write.csv(periods, file="results/Periodicity.csv", row.names=F)
periods # Period for HV is equal to number of years used for spectral analysis!
lapply(popnr_list, nrow)

p_spec <- ggplot(data=res_spect, aes(x=spec_x, y=spec_y)) +
  facet_wrap(~Area) +
  geom_line(size=1, col="grey27") +
  labs(x="Periodicity in years", y="Spectral Density") +
  geom_vline(aes(xintercept=peak), linetype="dashed", size=1.5, col="black")+
  scale_x_continuous(breaks=seq(0,50,by=5)) + scale_y_continuous(lim=c(0,25), breaks=seq(0,30,by=5)) +
  theme(strip.text=element_text(size=18), axis.text = element_text(size=15), axis.title=element_text(size=18))+
  background_grid(major = "xy") #get grid lines
p_spec
# ggsave(filename="results/SpectrDens.png", plot=p_spec , device="png", width=220, height=150, units="mm", dpi="print")

rm(l_spectr, res_spect) # clean up


# Autocorrelation function (ACF)
#-----------------------------------
ci_lev <- 0.95 # confidence intervals to plot

l_ACF <- lapply(popnr_list, function(X) {
  
  if(X$AreaShortName[1]=="HV"){
    X <- filter(X, grp==2)
  }# HV has a data gap, disregard first two years
  
  ACF <- acf(X$FemLogAbund, lag.max = 20, # Log abundance
             type = c("correlation"),
             plot = FALSE, na.action = na.pass, demean = TRUE)
  
  #plot(ACF, main = 'Doorwerth') # where are the confidence intervals coming from?
  
  # Confidence values calculated like this
  ci <- qnorm((1 + ci_lev)/2)/sqrt(ACF$n.used)
  # Source: https://stackoverflow.com/questions/14266333/extract-confidence-interval-values-from-acf-correlogram
  
  d_acf <- data.frame(ACF=ACF$acf, Lag=ACF$lag, ci=ci, Area=unique(X$Area))
  d_acf$sign <- ifelse(d_acf$ACF>0, d_acf$ACF>d_acf$ci, d_acf$ACF<d_acf$ci*-1)
  d_acf
})
lapply(l_ACF, head)

res_ACF <- data.table::rbindlist(l_ACF, use.names=TRUE)
head(res_ACF)
# write.csv(res_ACF, file="results/ACF.csv", row.names=F)

filter(res_ACF, sign==TRUE) # lags that pass the confidence interval and whose ACF are thus significantly different from zero P<0.05?

p_ACF <- ggplot(data=res_ACF, aes(x=Lag, y=ACF)) +
  facet_wrap(~Area)+
  geom_bar(stat="identity", width=0.1, alpha=0.8, fill="grey27") +
  geom_hline(yintercept=0)+
  geom_hline(aes(yintercept=ci), col="black", linetype="dashed", size=1)+ 
  geom_hline(aes(yintercept=-ci), col="black", linetype="dashed", size=1)+ # Confidence intervals
  scale_x_continuous(breaks=seq(0,20,by=5)) + scale_y_continuous(lim=c(min(res_ACF$ACF), 1), breaks=seq(-1,1,by=0.5))+
  theme(strip.text=element_text(size=18), axis.text = element_text(size=15), axis.title=element_text(size=18))
p_ACF
# ggsave(filename="results/ACF.png", plot=p_ACF , device="png", width=220, height=150, units="mm", dpi="print")

rm(l_ACF, ci_lev) # clean up


#---------------------------------------------
# Identifying order of feedback structures
#---------------------------------------------

# AR models ####
# AR = log-linear autoregressive model
# ARIMA = AutoRegressive Integrated Moving Average, uses AR together with moving averages

#kpss test for stationarity
kpss.test <- lapply(popnr_list, function(X) {
  
  if(X$AreaShortName[1]=="HV"){
    X <- filter(X, grp==2)
  }# HV has a data gap, disregard first two years
  
  p <- tseries::kpss.test(X$FemLogAbund, null="Trend")$p.value
  return(data.frame(Area=X$Area[1], KPSS.p=p))
  })
kpss.test <- data.table::rbindlist(kpss.test) # all populations are stationary

forecast::tsdisplay(popnr_list[[1]]$FemLogAbund) #display LogAbund over time, ACF and PACF

lapply(popnr_list, function(X) { # Run AR models
  # ARIMA model(p,d,q): 
  # p=order of the AR model
  # d=order of differencing (a.k.a. detrending)
  # q=order of moving average (a.k.a. smoothing)
  # if d & q zero, then basically fitting a normal AR model
  # source: https://towardsdatascience.com/time-series-analysis-with-auto-arima-in-r-2b220b20e8ab
  
  if(X$AreaShortName[1]=="HV"){
    X <- filter(X, grp==2)
  }# HV has a data gap, disregard first two years
  
  print(X$Area[1])
  print(nrow(X))
  fit <- forecast::auto.arima(X$FemLogAbund, trace = TRUE, stepwise = FALSE, d=0, max.q=0, approximation=FALSE, ic='aicc')
  print(fit)
  cat("\n\n")
  
  # Visualize
  png(paste("results/AR_", X$Area[1], ".png", sep=""), width=250, height=150, units="mm", res=300)
  plot <- ggplot(data=data.frame(Fitted=as.numeric(fit$fitted), LogAbund=X$FemLogAbund, YearCatch=X$YearCatch),
                 aes(x=YearCatch, y=Fitted))+
    geom_line(size=1, linetype="dashed", col="grey50") +
    geom_line(aes(y=LogAbund), size=1) +
    geom_point(aes(y=LogAbund), size=3) +
    labs(y="Female density (log-scale)") +
    scale_x_continuous(lim=c(1993.5,2021.5), breaks=seq(1990,2020,by=5)) + 
    scale_y_continuous(lim=c(-0.4,4.6), breaks=seq(0,5,by=0.5))+
    theme(axis.text = element_text(size=15), axis.title.x=element_blank(), axis.title=element_text(size=18)) +
    background_grid(major = "xy") #get grid lines
  print(plot)
  dev.off()
}) %>% capture.output(file="results/AR_output.txt") # capture directly into .txt file because a lot of output printed on screen


# PRCF function ####
# Partial rate correlation function, based on example script by Luis Cayuela

PRCF <- function(df_Abund, log=TRUE, lagmax=15){
  
  if(df_Abund$AreaShortName[1]=="HV"){
    df_Abund <- filter(df_Abund, grp==2)
  }# HV has a large data gap, disregard first two years
  
  # Calculate R (log per capita growth rate)
  percapR <- vector(mode="numeric") # Initialize
  tmax <- nrow(df_Abund)-1
  
  for(year in 1:tmax){ # loop through years
    
    # For timeseries that are not Log transformed
    if(log==TRUE) { 
      d_Abund <- df_Abund$NoFemCor
      df_Abund$FemLogAbund <- log(df_Abund$NoFemCor + 1)
      
      d_Abund2 <- d_Abund +1 # +1 to correct for zero values in Log transformation
      percapRi<-log(d_Abund2[year+1]/d_Abund2[year])   # Calculate per capita growth rate
      percapR<-c(percapR,percapRi)   # Add per capita growth rate [i] to vector
      
      rm(d_Abund2)
      
    } else { # For timeseries that are log transformed                   
      LogAbund <- df_Abund$FemLogAbund
      
      percapRi <- LogAbund[year+1] - LogAbund[year]
      percapR <- c(percapR,percapRi)
    }
  }
  rm(year, percapRi)
  class(percapR)
  
  # Calculate correlations
  df <- data.frame(R=percapR, Abund=df_Abund$NoFemCor[c(1:tmax)]) # needs to be in same order
  df <- na.omit(df) # drop any NA values
  
  cors <-cor(df$R, df$Abund) # Correlate R with population numbers per year
  
  # Partial autocorrelation function, but has wrong null model for biological populations [Berryman & Turchin 2001, Oikos(92)]
  PACF <- pacf(df_Abund$FemLogAbund, lag.max = lagmax, plot = FALSE, na.action = na.pass)$acf 
  PRCF <- PACF
  PRCF[1]<-cors # Replace first value PACF with calculated correlation coefficient to correct null model
  
  # Results
  res_PRCF <- data.frame(Lag=seq(1, lagmax, by=1), PRCF=PRCF, Area=df_Abund$Area[1])
  res_PRCF$barlett <- 2/sqrt(length(na.omit(df_Abund$NoFemCor))) # Calculate confidence intervals
  
  return(res_PRCF)
}

list_PRCF <- lapply(popnr_list, function(X) PRCF(X, log=FALSE, lagmax=15))
res_PRCF <- data.table::rbindlist(list_PRCF, use.names = TRUE) %>% 
  mutate(sign=ifelse(PRCF<0, PRCF<barlett*-1, PRCF>barlett))
head(res_PRCF)
# write.csv(res_PRCF, file="results/PRCF.csv", row.names=F)
rm(list_PRCF)

# Visualize
p_PRCF <- ggplot(data=res_PRCF, aes(x=Lag, y=PRCF)) +
  facet_wrap(~Area)+
  geom_bar(stat="identity", col="black", fill="grey")+
  geom_hline(yintercept=0)+
  geom_hline(aes(yintercept=barlett), col="black", linetype="dashed", size=1)+ 
  geom_hline(aes(yintercept=-barlett), col="black", linetype="dashed", size=1)+ # Confidence intervals
  scale_x_continuous(lim=c(0.5,10.5), breaks=seq(1,10,by=1)) + # Plot only until lag10
  scale_y_continuous(lim=c(-1, 1), breaks=seq(-1,1,by=0.5))+
  theme(strip.text=element_text(size=18), axis.text = element_text(size=15), axis.title=element_text(size=18))
p_PRCF
# ggsave(filename="results/PRCF.png", plot=p_PRCF , device="png", width=220, height=150, units="mm", dpi="print")


# Conclusion: Best feedback order depends on Population, but always <=4 orders ####


#-----------------------------------
# Population growth models ####
#-----------------------------------

# Model type: linear model (LM)
# Dependent variable = R = realized per-capita growth rate = ln(Nt)-ln(Nt-1) = ln(Nt/Nt-1)
# Full model: R ~ ln(Nt-1) + ln(Nt-2) + ln(Nt-3) + ln(Nt-4) + Mismatch.mean + Mismatch.sd + Eggs + Area
# double check NAs

# Prep data ####
#-----------------------------------
lapply(popnr_list, head)

# Add lag variables
popnrlag <- list()
for(pop in 1:length(popnr_list)){
  Area <- popnr_list[[pop]]$AreaShortName[1]
  
  # Deal with HV data gap
  if(Area=="HV"){
    d_lags1 <- popnr_list[[pop]] %>% filter(grp==1) %>% # first two years
      mutate(LogAbund_1=NA, LogAbund_2=NA, LogAbund_3=NA, LogAbund_4=NA) # create lag variables
    d_lags1$LogAbund_1[2:nrow(d_lags1)] <- d_lags1$FemLogAbund[1:(nrow(d_lags1)-1)]
    
    d_lags2 <- popnr_list[[pop]] %>% filter(grp==2) %>% # years after data gap
      mutate(LogAbund_1=NA, LogAbund_2=NA, LogAbund_3=NA, LogAbund_4=NA) # create lag variables
    d_lags2$LogAbund_1[2:nrow(d_lags2)] <- d_lags2$FemLogAbund[1:(nrow(d_lags2)-1)]
    d_lags2$LogAbund_2[3:nrow(d_lags2)] <- d_lags2$FemLogAbund[1:(nrow(d_lags2)-2)]
    d_lags2$LogAbund_3[4:nrow(d_lags2)] <- d_lags2$FemLogAbund[1:(nrow(d_lags2)-3)]
    d_lags2$LogAbund_4[5:nrow(d_lags2)] <- d_lags2$FemLogAbund[1:(nrow(d_lags2)-4)]
    
    d_lags <- rbind(d_lags1, d_lags2)
    rm(d_lags1, d_lags2)
    
  } else{
    d_lags <- popnr_list[[pop]] %>% mutate(LogAbund_1=NA, LogAbund_2=NA, LogAbund_3=NA, LogAbund_4=NA) # create lag variables
    d_lags$LogAbund_1[2:nrow(d_lags)] <- d_lags$FemLogAbund[1:(nrow(d_lags)-1)]
    d_lags$LogAbund_2[3:nrow(d_lags)] <- d_lags$FemLogAbund[1:(nrow(d_lags)-2)]
    d_lags$LogAbund_3[4:nrow(d_lags)] <- d_lags$FemLogAbund[1:(nrow(d_lags)-3)]
    d_lags$LogAbund_4[5:nrow(d_lags)] <- d_lags$FemLogAbund[1:(nrow(d_lags)-4)]
  }

  d_lags <- d_lags %>% mutate(R=as.numeric(FemLogAbund-LogAbund_1))
  
  popnrlag[[pop]] <- d_lags
}
popnrlag <- data.table::rbindlist(popnrlag)
rm(d_lags, Area)

# Add explanatory variables
popnrlag <- merge(popnrlag, d_mismatch[,c("YearHatch", "AreaShortName", "Mismatch.Tree", "sd.Tree")], 
                   by.x=c("YearCatch", "AreaShortName"), by.y=c("YearHatch", "AreaShortName"), all.x=TRUE)
popnrlag <- merge(popnrlag, Eggs, by.x=c("YearCatch", "Area"), by.y=c("Lag1", "AreaName"), all.x=TRUE)

# Keep only necessary variables
popnrlag <- popnrlag %>% mutate(Mismatch.mean=Mismatch.Tree, Mismatch.sd=sd.Tree, Eggs_1=Eggs) %>%
  select(YearCatch, Area, grp, FemLogAbund, R, LogAbund_1, LogAbund_2, LogAbund_3, LogAbund_4, Mismatch.mean, Mismatch.sd, Eggs_1)
head(popnrlag)
str(popnrlag)
table(popnrlag$Area, popnrlag$YearCatch, is.na(popnrlag$Eggs_1)) # see years which are dropped when including Eggs_1 as covariate in the model


# Fit models ####
#-----------------------------------
lm_full<- lm(R ~ (LogAbund_1 + LogAbund_2 + LogAbund_3 + LogAbund_4)*Area + Mismatch.mean + Mismatch.sd + Eggs_1, data=popnrlag)
anovf <- car::Anova(lm_full) %>% as.data.frame() # none of interactions significant
anovf$model <- "full"
anovf$N <- nrow(popnrlag) - (nrow(popnrlag) - length(lm_full$residuals))
summary(lm_full) # 46 observations dropped (almost half!)
AICcmodavg::AICc(lm_full)

# Backward elimination ####
lm_fit <- ols(R ~ -1 + (LogAbund_1 + LogAbund_2 + LogAbund_3 + LogAbund_4)*Area + Mismatch.mean + Mismatch.sd + Eggs_1, 
              data=popnrlag, x=T, y=T) #x & y = store data for later use (e.g. bootstrapping)
lm_bw <- rms::fastbw(lm_fit, rule="aic") # fast backward elimination, doesn't have AICc
lm_bw # only LogAbund_4 kept

lm_bw_val <- rms::validate(lm_fit, B=500, bw=T) # validate backward elimination with bootstrapping
print(lm_bw_val, B=10) # but not very robust, number of factors retained very variable


# Model comparisons by hand ####
lm1 <- update(lm_full, .~. - LogAbund_1:Area -LogAbund_2:Area -LogAbund_3:Area -LogAbund_4:Area) # simplify model
anova1 <- car::Anova(lm1)
anova1$model <- "lm1"
anova1$N <- nrow(popnrlag) - (nrow(popnrlag) - length(lm1$residuals))
summary(lm1) # 46 observations dropped (almost half!)

lm2 <- update(lm1, .~. - Eggs_1) # simplify model: remove Eggs_1 (not variable of interest, not significant, and many missing values)
anova2 <- car::Anova(lm2)
anova2$model <- "lm2"
anova2$N <- nrow(popnrlag) - (nrow(popnrlag) - length(lm2$residuals))
summary(lm2) # only 25 observations dropped, probably why now Mismatch.mean significant, more data = more power

lm3 <- update(lm2, .~. - Mismatch.sd) # remove Mismatch.sd (not variable of interest and not significant)
anova3 <- car::Anova(lm3)
anova3$model <- "lm3"
anova3$N <-nrow(popnrlag) - (nrow(popnrlag) - length(lm3$residuals))
summary(lm3) # removing Mismatch.sd does not influence estimate for Mismatch.mean that much


# Final model ####
#-----------------------------------
lm_final <- lm(R ~ LogAbund_1 + LogAbund_2 + LogAbund_3 + LogAbund_4 + Mismatch.mean + Area, data=popnrlag)
summary(lm_final)
car::Anova(lm_final)
AICcmodavg::AICc(lm_final)

varexpl <- effectsize::omega_squared(lm_final) %>% as.data.frame # measure for how much variance each estimate explains

coef <- summary(lm_final)$coefficients %>% as.data.frame()
coef$Parameter <- row.names(coef)
coef <- coef %>% mutate(Area=ifelse(grepl("Area", Parameter)==T, gsub("Area(\\w+)","\\1",Parameter), NA),
                                  Parameter=ifelse(grepl("Area", Parameter)==T, "Area", Parameter))
coef <- merge(coef, varexpl, by="Parameter", all.x=T)
coef

plot(lm_final) # residuals check: look ok

res_anova <- rbind(anovf, anova1, anova2, anova3)
res_anova

# Mismatch effect size
# How to back transform: https://data.library.virginia.edu/interpreting-log-transformations-in-a-linear-model/
(exp(lm_final$coefficients["Mismatch.mean"])-1)*100 # effect size mismatch
(exp(filter(coef, Parameter=="Mismatch.mean")$`Std. Error`)-1)*100 # SE of that effect size

coef$EffectSize <- NA
coef$EffectSize[2:nrow(coef)] <- (exp(coef$Estimate[2:nrow(coef)])-1)*100

# write.csv(res_anova, file="results/anova.csv", row.names=T)
# write.csv(coef, file="results/output.csv", row.names=F)


# Plot effect
popnrlag$row <- row.names(popnrlag)
d <- na.omit(popnrlag[,c("row","LogAbund_1", "LogAbund_2", "LogAbund_3", "LogAbund_4", "Area", "Mismatch.mean","R")])
d <- merge(d, popnrlag[,c("row", "YearCatch")], by="row") %>% arrange(YearCatch, Area) # HV 1995 dropped anyway
# write.csv(select(d, -row), file="results/Observations_included.csv", row.names=F)

aggregate(row ~ Area,data=d, length)# number of observations included per Area

m_res <- lm(R ~ -1 + LogAbund_1 + LogAbund_2 + LogAbund_3 + LogAbund_4 + Area, data=d)
d$resid <- m_res$residuals # get residuals from model with all fixed effects except Mismatch.mean

m_pred <- lm(resid ~ Mismatch.mean, data=d) # reestimate Mismatch.mean effect on residuals
summary(m_pred) # bit smaller, but similar effect size as final model
d$pred <- predict(m_pred, new.data=d)
head(d)

p_eff <- ggplot(data=d, aes(x=Mismatch.mean, y=resid)) + # plot Mismatch effect, filter(d, Mismatch.mean>-21 & Mismatch.mean<0)
  geom_vline(xintercept=0, linetype="dashed", size=1, col="grey50")+
  geom_point(size=5, shape=21, fill="dodgerblue2", col="black", alpha=0.5) + #"grey27", alpha=0.7
  geom_line(aes(y=pred), size=1.5, col="black")+
  #geom_smooth(method=lm)+ # to check what happens to slope if leave out some extremes, slope becames a bit less steep
  labs(x="Mismatch with oak budburst date (days)", y="Residual R")+
  scale_x_continuous(breaks=seq(-30,10,by=5)) + scale_y_continuous(lim=c(-1.5, 1.5),breaks=seq(-1.5,1.5,by=0.5))+
  theme(strip.text=element_text(size=18), axis.text = element_text(size=15), axis.title=element_text(size=18))
p_eff
# ggsave(filename="_results/Mismatch_effect_rev.png", plot=p_eff , device="png", width=220, height=150, units="mm", dpi="print")


sessionInfo() %>% capture.output(file="src/env_PopDyn_analysis.txt")
