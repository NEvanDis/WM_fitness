# Analysis of Phenological mismatch experiment 2021 ####
# Manipulated timing of egg hatching of eggs from wild Mothers caught in 2020
# Either hatching on day of budburst (Day0), before (Day-4 to -1), or after (Day+1 to +5)
# Disentangle effects of photoperiod and food quality: photoperiod treatment (changing or constant)


# Open R project in main folder

# Load packages
#-----------------------------------
library(tidyverse)
library(readxl)
library(cowplot)
theme_set(theme_cowplot()) #white background instead of grey -> don't load if want grey grid
library(coxme)
library(lme4)
library(lmerTest)


# Load data ####
#-----------------------------------
d <- read.csv("1_data/CatFood2021_deposit.csv")
head(d)

length(unique(d$TubeID)) # should be 22 mothers
table(d$Treatment) # photoperiod and mismatch treatment coded in one variable


# Descriptives
#-----------------------------------

# N per Area
table(d[!duplicated(d$TubeID), "AreaShortName"])


#---------------------------------------------------------------------------------------------------------------------------
# Fitness curve ####
#----------------------------------
# RQ1: What are the fitness consequences of day to day timing (a)synchrony with budburst? ####

# Survival data ####
d_surv <- d %>% mutate(PhotoTreat=gsub("(\\w+)Day.+","\\1",Treatment), MismTreat=gsub("\\w+(Day.+)","\\1",Treatment)) %>% 
  select(TubeID, Treatment, PhotoTreat, MismTreat, CaterpillarID, DeadAprilDay, PupationAprilDay) %>%
  pivot_longer(cols=c(DeadAprilDay, PupationAprilDay), names_to="Info", values_to="TimeOfEvent") %>%
  filter(!is.na(TimeOfEvent)) %>%
  mutate(Event=ifelse(Info=="DeadAprilDay", 1, 0), Treatment=as.factor(Treatment), PhotoTreat=as.factor(ifelse(PhotoTreat=="Chang", "Changing", "Constant")), MismTreatf=as.factor(MismTreat), TubeID=as.factor(TubeID)) %>%
  mutate(Treatment=factor(Treatment, levels=c("ChangDay-4", "ChangDay-3", "ChangDay-2", "ChangDay-1", "ChangDay0", "ChangDay+1", "ChangDay+2", "ChangDay+3", "ChangDay+4",
                                              "ChangDay+5",  "ConstDay-4", "ConstDay-2", "ConstDay0", "ConstDay+2", "ConstDay+4")), 
    MismTreatf=factor(MismTreat, levels=c("Day-4", "Day-3", "Day-2", "Day-1", "Day0", "Day+1", "Day+2", "Day+3", "Day+4", "Day+5")),
    MismTreat=as.numeric(gsub("Day(.+)","\\1",MismTreat))) %>%
  mutate(MismTreat1=MismTreat+5, # no negatives to be able to fit squared term
         MismTreat2=(MismTreat+5)^2) # squared term to add in model
# Vidisha coded it as TimeOfEvent=DeadAprilDay or PupationAprilDay, with event=Died or Survived
head(d_surv)
str(d_surv)
table(d_surv$MismTreat2)

length(unique(d_surv$CaterpillarID)) # should be 976


#-----------------------------------
# Survival analysis ####
#-----------------------------------
levels(d_surv$Treatment)
levels(d_surv$PhotoTreat)
levels(d_surv$MismTreatf) # as factor or not? Marcel thinks not ####
levels(d_surv$TubeID)
table(d_surv$TimeOfEvent)
table(d_surv$Event)

# Visualize survival probabilities ####
head(d_surv)

surv_probs <- aggregate(Event~MismTreat + PhotoTreat + TubeID, d_surv, sum) # per mother
surv_probs$samplesize <- aggregate(Info~MismTreat + PhotoTreat + TubeID, d_surv, length)$Info
surv_probs$probs <- 100 - (surv_probs$Event/surv_probs$samplesize*100) # event = death
head(surv_probs)

surv_avg <- Rmisc::summarySE(surv_probs, measurevar="probs", groupvars=c("MismTreat")) # average of two photoperiod treatments
surv_avg$samplesize <- aggregate(Info~MismTreat, d_surv, length)$Info
surv_avg

raw_surv <- ggplot(data=surv_avg, aes(x=MismTreat, y=probs))+
  scale_colour_manual(values=c("grey27", "orangered2"))+ #"dodgerblue4"
  geom_jitter(data=surv_probs, aes(col=PhotoTreat), alpha=0.3, size=3, height=0.5, width=0.25)+
  geom_point(size=5, col="black") +
  geom_errorbar(aes(ymax = probs+se, ymin=probs-se), width=0.3, col="black") +
  geom_text(aes(label=samplesize, y=probs+8.3), col="black", size=4,fontface="bold")+ # N caterpillars in each treatment
  labs(y="Survival (%)", x="Mismatch with oak budburst date (days)")+
  scale_y_continuous(breaks=seq(0,100, by=10))+ scale_x_continuous(breaks=seq(-4, 5, by=1))+
  theme(axis.title.y=element_text(size=18, vjust=2), axis.title.x=element_text(size=18, vjust=-0.5),
        axis.text=element_text(size=16), legend.text = element_text(size=16), legend.title=element_text(size=17))
raw_surv 
# ggsave(filename="_results/Survival_raw.png", plot=raw_surv , device="png", width=200, height=150, units="mm", dpi="print")


# Fit binomial model ####
#-----------------------------------
head(d_surv) # test if probability of survival differs between treatments

glm1 <- glmer(Event ~ (MismTreat1 + MismTreat2)*PhotoTreat + (1|TubeID), family=binomial, data=d_surv,
              na.action="na.fail", control=glmerControl(calc.derivs=F)) # helps convergence
anova1 <- drop1(glm1,test="Chi") %>% as.data.frame # interaction not significant
anova1$mod <- "glm1"

glm2 <- update(glm1, ~ . -MismTreat1:PhotoTreat - MismTreat2:PhotoTreat) # simplify model
anova2 <- drop1(glm2,test="Chi") %>% as.data.frame #no effect of PhotoTreatment, but effect of MismTreat and MismTreat^2
anova2$mod <- "glm2"

# Final model ####
glm_final <- glm2
summary(glm_final)# Estimates are log odds
glm_res <- summary(glm_final)$coefficients %>% as.data.frame

# write.csv(glm_res, file="_results/output_Surv_glmer.csv", row.names=T)
# write.csv(rbind(anova1, anova2), file="_results/anova_Surv_glmer.csv", row.names=T)

# Get predictions ####
glm.pred <- d_surv[!duplicated(d_surv[,c("TubeID", "Treatment")]),] # each replicate assigned same prediction, so remove duplicates
glm.pred$pred <- predict(glm_final, newdata=glm.pred, type="response") # predictions are probability of dying now
glm.pred$survprob <- (1-glm.pred$pred)
aggregate(survprob~MismTreat, data=glm.pred, mean) # peak at Day2
glm.pred$rel <- glm.pred$survprob/mean(filter(glm.pred, MismTreat==1)$survprob) # expressive relative to peak
head(glm.pred)

# Visualize predictions ####
pred <- Rmisc::summarySE(glm.pred, measurevar="survprob", groupvars=c("MismTreat")) # average of two photoperiod treatments
pred$samplesize <- aggregate(CaterpillarID~MismTreat, data=d_surv, length)$CaterpillarID

# Add predictions to raw data figure
p_surv <- raw_surv + #geom_line(data=pred, aes(y=survprob*100)) +
  geom_smooth(data=pred, aes(y=survprob*100), se=F, col="red3")
p_surv
# ggsave(filename="_results/Survival_wpred_rev.png", plot=p_surv, device="png", width=200, height=150, units="mm", dpi="print")

rm(anova1, anova2, glm_res, glm1, glm2, pred, surv_probs, surv_avg, raw_surv) #cleanup


#-----------------------------------
# Pupation weight analysis ####
#-----------------------------------
head(d)

d_pupa <- d %>% mutate(PhotoTreat=gsub("(\\w+)Day.+","\\1",Treatment), MismTreat=gsub("\\w+(Day.+)","\\1",Treatment), PupaWeight=PupaWeight_ingrams*1000) %>%
  select(ExperimentName, TubeID, Treatment, PhotoTreat, MismTreat, CaterpillarID, PupationAprilDay, PupaWeight) %>%
  filter(!is.na(PupationAprilDay)) %>%
  mutate(Treatment=as.factor(Treatment), PhotoTreat=as.factor(ifelse(PhotoTreat=="Chang", "Changing", "Constant")), MismTreatf=as.factor(MismTreat), TubeID=as.factor(TubeID)) %>%
  mutate(Treatment=factor(Treatment, levels=c("ChangDay-4", "ChangDay-3", "ChangDay-2", "ChangDay-1", "ChangDay0", "ChangDay+1", "ChangDay+2", "ChangDay+3", "ChangDay+4",
                                              "ChangDay+5",  "ConstDay-4", "ConstDay-2", "ConstDay0", "ConstDay+2", "ConstDay+4")), 
         MismTreatf=factor(MismTreat, levels=c("Day-4", "Day-3", "Day-2", "Day-1", "Day0", "Day+1", "Day+2", "Day+3", "Day+4", "Day+5")),
         MismTreat=as.numeric(gsub("Day(.+)","\\1",MismTreat))) %>%
  mutate(MismTreat1=MismTreat+5, # no negatives
         MismTreat2=(MismTreat+5)^2) # squared term to add in model
head(d_pupa) # 346 individuals survived until pupation
nrow(d_pupa)/nrow(d)*100 # ~35%


# Visualize ####
weight <- Rmisc::summarySE(d_pupa, measurevar="PupaWeight", groupvars=c("MismTreat", "PhotoTreat"))
weight$pos <- ifelse(is.na(weight$se)==T, 0, weight$se) # position of sample size labels
weight

raw_weight <- ggplot(data=weight, aes(x=MismTreat, y=PupaWeight, col=PhotoTreat, fill=PhotoTreat))+
  scale_colour_manual(values=c("grey27", "orangered2"))+ #"dodgerblue4"
  scale_fill_manual(values=c("grey27", "orangered2"))+
  geom_jitter(data=d_pupa, aes(col=PhotoTreat), alpha=0.3, size=3, height=0, width=0.25)+ #alpha=0.3, size=2, height=0, width=0.25
  geom_errorbar(data=filter(weight, PhotoTreat=="Changing"), aes(ymax = PupaWeight+se, ymin=PupaWeight-se), width=0.3, col="black") +
  geom_errorbar(data=filter(weight, PhotoTreat=="Constant"), aes(ymax = PupaWeight+se, ymin=PupaWeight-se), width=0.3, col="orangered4") +
  geom_point(size=5, shape=21, col="black")+
  #geom_text(data=filter(weight, PhotoTreat=="Changing"),aes(label=N, y=PupaWeight-pos-2.3), col="black", size=4, fontface="bold")+
  #geom_text(data=filter(weight, PhotoTreat=="Constant"),aes(label=N, y=PupaWeight+pos+2.3), col="black", size=4, fontface="bold")+
  labs(y="Weight at pupation (mg)", x="Mismatch with oak budburst date (days)")+
  scale_y_continuous(breaks=seq(15,75, by=10))+ scale_x_continuous(breaks=seq(-4, 5, by=1))+
  theme(axis.title.y=element_text(size=18, vjust=2), axis.title.x=element_text(size=18, vjust=-0.5),
        axis.text=element_text(size=16), legend.text = element_text(size=16), legend.title=element_text(size=17))
raw_weight  
# ggsave(filename="_results/PupWeight_raw.png", plot=raw_weight , device="png", width=200, height=150, units="mm", dpi="print")


# Fit linear mixed model ####
#-----------------------------------
lm1 <- lmer(PupaWeight ~ (MismTreat1 + MismTreat2)*PhotoTreat + (1|TubeID), data=d_pupa)
anova1 <- anova(lm1) %>% as.data.frame() # interaction not significant
anova1$mod <- "lm1"

lm2 <- update(lm1, ~ . - MismTreat1:PhotoTreat - MismTreat2:PhotoTreat) # simplify model
anova2 <- anova(lm2) %>% as.data.frame() # Squared mismatch not significant
anova2$mod <- "lm2"

lm3 <- update(lm2, ~ . - MismTreat2) # simplify model
anova3 <- anova(lm3) %>% as.data.frame() # PhotoTreat and MismTreat significant
anova3$mod <- "lm3"

# Still there if exclude first time point with low sample size?
lm4 <- lmer(PupaWeight ~ -1 + MismTreat1 + PhotoTreat + (1|TubeID), data=filter(d_pupa, MismTreat!=-4))
anova(lm4) # yes


# Final model ####
lm_final <- lm3
summary(lm_final)
lm_res <- summary(lm_final)$coefficients %>% as.data.frame

plot(lm_final) #equal variance? ok
qqnorm(resid(lm_final)) #normally distributed? ok
qqline(resid(lm_final))

# write.csv(lm_res, file="_results/output_PupaWeight_lmer.csv", row.names=T)
# write.csv(rbind(anova1, anova2, anova3), file="_results/anova_PupaWeight_lmer.csv", row.names=T)

# Get predictions ####
lm.pred <- d_pupa[!duplicated(d_pupa[,c("TubeID", "Treatment")]),] # each replicate assigned same prediction, so remove duplicates
lm.pred$pred <- predict(lm_final, newdata=lm.pred, type="response")
head(lm.pred)

# Visualize predictions ####
pred1 <- Rmisc::summarySE(lm.pred, measurevar="pred", groupvars=c("MismTreat", "PhotoTreat")) # significant effect of photoperiod, so show separate means
pred1$samplesize <- weight$N
pred1$pos <- ifelse(is.na(pred1$se)==T, 0, pred1$se) # position of sample size labels

# add predictions to raw data figure
p_weight <- raw_weight + #geom_line(data=pred1, aes(y=pred)) +
  geom_smooth(data=pred1, aes(y=pred, col=PhotoTreat), se=F, method=lm)+
  geom_text(data=filter(weight, PhotoTreat=="Changing"),aes(label=N, y=PupaWeight-pos-1.5), col="black", size=4, fontface="bold")+
  geom_text(data=filter(weight, PhotoTreat=="Constant"),aes(label=N, y=PupaWeight+pos+2.3), col="black", size=4, fontface="bold")
p_weight
# ggsave(filename="_results/PupWeight_wpred_rev.png", plot=p_weight, device="png", width=200, height=150, units="mm", dpi="print")

rm(anova1, anova2, anova3, lm1, lm2, lm3, lm4, lm_res, raw_weight, weight, pred1) # clean up


#--------------------------------------------
# Get fitness curve ####
#--------------------------------------------

# Don't care about PhotoTreat effect, drop from models ####
glm_fit <- glmer(Event ~ MismTreat1 + MismTreat2 + (1 | TubeID), family="binomial", data=d_surv)
lm_fit <- lmer(PupaWeight ~ MismTreat1 + (1 | TubeID), data=d_pupa)

# Get predictions to use for curve ####
glm.fit <- d_surv[!duplicated(d_surv[,c("TubeID", "MismTreatf")]),] # each replicate assigned same prediction, so remove duplicates
glm.fit$pred <- predict(glm_fit, newdata=glm.fit, type="response") # predictions are probability of dying now
glm.fit$survpred <- (1-glm.fit$pred)

lm.fit <- d_pupa[!duplicated(d_pupa[,c("TubeID", "MismTreatf")]),] # each replicate assigned same prediction, so remove duplicates
lm.fit$pred <- predict(lm_fit, newdata=lm.fit, type="response")

head(glm.fit) # pred = probability of dying, survpred=1-pred, 220 observations = 22 mothers * 10 MismTreat groups
head(lm.fit) # pred=predicted weight from lmer, only 154 observations


# Fit curve to absolute fitness ####
#-----------------------------------
RelFit <- merge(glm.fit[,c("TubeID", "MismTreat", "survpred")], lm.fit[,c("TubeID", "MismTreat", "pred")], by=c("TubeID", "MismTreat"), all=T)
colnames(RelFit)[c(3,4)] <- c("survpred", "pupwpred")
RelFit$Fit <- RelFit$survpred*RelFit$pupwpred # multiply absolute values
head(RelFit) # can only do for 154 observations, clutches with >=1 caterpillar surviving until pupation
table(RelFit$MismTreat, is.na(RelFit$Fit)) # for the other clutches, fitness = 0
RelFit$Fit2 <- ifelse(is.na(RelFit$Fit)==T, 0, RelFit$Fit)

# loess model to describe the curve ####
loess_mod <- loess(Fit2~ -1 + MismTreat,  data=RelFit) 
summary(loess_mod)

curve <- RelFit[!duplicated(RelFit[,c("MismTreat")]),] %>% select(MismTreat)
curve$pred <- predict(loess_mod, newdata=curve)
curve <- arrange(curve, MismTreat)
curve # peak at day2
curve$rel <- curve$pred/filter(curve, MismTreat==2)$pred # expressive relative to peak

RelFit$rel <- RelFit$Fit2/mean(filter(RelFit, MismTreat==2)$Fit2)

RelFit_means <- Rmisc::summarySE(RelFit, measurevar="rel", groupvars=c("MismTreat"))
#RelFit_means$samplesize <- aggregate(TubeID~MismTreat, data=d_pupa, length)$TubeID # number of caterpillars curve is based on
RelFit_means$samplesize <- aggregate(TubeID~MismTreat, data=d_surv, length)$TubeID # number of caterpillars curve is based on = all
RelFit_means$curve <- curve$rel
head(RelFit_means)
# write.csv(RelFit_means, file="_results/RelFitness_rev.csv", row.names=F)

p_relfit <- ggplot(data=RelFit, aes(x=MismTreat, y=rel)) +
  geom_jitter(size=3, alpha=0.4, height=0, width=0.25, shape=21, fill="dodgerblue2", col="dodgerblue4")+
  geom_point(data=RelFit_means, size=5) +
  geom_errorbar(data=RelFit_means, aes(ymax = rel+se, ymin=rel-se), width=0.3) +
  #geom_line(data=curve, aes(y=pred), col="darkred", size=1)+
  #geom_smooth(data=curve, se=F, col="red3", size=1)+ # makes it look more like a smooth line, but otherwise exactly the same as using geom_line
  geom_text(data=RelFit_means,aes(label=samplesize, y=rel+0.12), col="black", size=5, fontface="bold")+ # number of caterpillars
  geom_hline(yintercept=1, linetype="dashed")+
  labs(y="Relative fitness", x="Mismatch with oak budburst date (days)")+
  scale_y_continuous(lim=c(0,1.21), breaks=seq(0,1.6, by=0.2))+ scale_x_continuous(breaks=seq(-4,5, by=1))+ #lim=c(-0.1,1.3)
  theme(legend.position="none")+
  theme(axis.title.y=element_text(size=18, vjust=2), axis.title.x=element_text(size=18, vjust=-0.5),
        axis.text=element_text(size=16), legend.text = element_text(size=16), legend.title=element_text(size=17))
p_relfit
# ggsave(filename="_results/FitnessCurve_rev.png", plot=p_relfit, device="png", width=200, height=150, units="mm", dpi="print")

rm(curve, RelFit_means) # clean up



#---------------------------------------------------------------------------------------------------------------------------
# Timing of life stages ####
#----------------------------------
# RQ2: Can food quality affect the timing of life stages? ####

# Timing data ####
d_tim <- d %>% 
  select(YearCatch, YearHatch, Treatment, TubeID, ClutchID, CaterpillarID, HatchAprilDay, PupationAprilDay, PupaWeight_ingrams) %>%
  filter(!is.na(PupationAprilDay)) %>%
  mutate(PhotoTreat=gsub("(\\w+)Day.+","\\1",Treatment), MismTreat=gsub("\\w+(Day.+)","\\1",Treatment), Event=1, PupaWeight=PupaWeight_ingrams) %>% # all included have pupated
  mutate(Treatment=as.factor(Treatment), PhotoTreat=as.factor(ifelse(PhotoTreat=="Chang", "Changing", "Constant")), MismTreatf=as.factor(MismTreat), TubeID=as.factor(TubeID)) %>%
  mutate(Treatment=factor(Treatment, levels=c("ChangDay-4", "ChangDay-3", "ChangDay-2", "ChangDay-1", "ChangDay0", "ChangDay+1", "ChangDay+2", "ChangDay+3", "ChangDay+4",
                                              "ChangDay+5",  "ConstDay-4", "ConstDay-2", "ConstDay0", "ConstDay+2", "ConstDay+4")), 
         MismTreatf=factor(MismTreat, levels=c("Day-4", "Day-3", "Day-2", "Day-1", "Day0", "Day+1", "Day+2", "Day+3", "Day+4", "Day+5")),
         MismTreat=as.numeric(gsub("Day(.+)","\\1",MismTreat))) %>%
  mutate(MismTreat1=MismTreat+5, # no negatives
         MismTreat2=(MismTreat+5)^2) %>% # squared term to add in model
  mutate(larv_devtime=PupationAprilDay-HatchAprilDay, CaterpillarID=as.factor(CaterpillarID)) #, pup_devtime=(AdultNovDate+214)-PupationAprilDay)
head(d_tim)
str(d_tim)


# Visualize ####
summary(d_tim$larv_devtime) # range 29-49

# get for each larv_devtime proportion of unpupated larvae per treatment
larv_prop <- aggregate(CaterpillarID~larv_devtime + Treatment, d_tim, length)
colnames(larv_prop)[ncol(larv_prop)] <- "Nr_cats"

tot_larv <- aggregate(CaterpillarID~ Treatment, d_tim, length) # sample sizes
colnames(tot_larv)[ncol(tot_larv)] <- "Total"

larv_prop <- merge(larv_prop, tot_larv, by="Treatment")

#add 0 timepoint
min(larv_prop$larv_devtime) # add 25 ok for plotting

add0 <- larv_prop[!duplicated(larv_prop$Treatment),] %>% 
  mutate(larv_devtime=25, Nr_cats=0)
add0

larv_prop <- rbind(larv_prop, add0)
larv_prop <- larv_prop %>% 
  mutate(PhotoTreat=as.factor(gsub("(\\w+)Day.+","\\1",Treatment)), MismTreat=as.factor(gsub("\\w+(Day.+)","\\1",Treatment))) %>%
  mutate(PhotoTreat=ifelse(PhotoTreat=="Chang", "Changing", "Constant"),
    Treatment=factor(Treatment, levels=c("ChangDay-4", "ChangDay-3", "ChangDay-2", "ChangDay-1", "ChangDay0", "ChangDay+1", "ChangDay+2", "ChangDay+3", "ChangDay+4",
                                              "ChangDay+5",  "ConstDay-4", "ConstDay-2", "ConstDay0", "ConstDay+2", "ConstDay+4")), 
         MismTreat=factor(MismTreat, levels=c("Day-4", "Day-3", "Day-2", "Day-1", "Day0", "Day+1", "Day+2", "Day+3", "Day+4", "Day+5")))
head(larv_prop)

treats <- as.list(unique(larv_prop$Treatment))

for(treat in 1:length(treats)){
  larv_sub <- filter(larv_prop, Treatment==treats[[treat]])
  
  larv_sub$accum <- NA
  larv_sub$accum[1] <- larv_sub$Nr_cats[1]
  
  for (row in 2:nrow(larv_sub)){
    larv_sub$accum[row] <- larv_sub$accum[row-1] + larv_sub$Nr_cats[row]
    larv_sub$accum[row] <- ifelse(is.na(larv_sub$accum[row]), larv_sub$accum[row-1], larv_sub$accum[row])
  }
  
  treats[[treat]] <- larv_sub
  
}
larv_prop <- data.table::rbindlist(treats) %>% arrange(larv_devtime, Treatment)

larv_prop$prop <- 1-(larv_prop$accum/larv_prop$Tot)
larv_prop$prop <- ifelse(larv_prop$larv_devtime==25, 1, larv_prop$prop)
head(larv_prop)

rm(treats, add0, larv_sub, tot_larv, row, treat) #clean up


# Plot larv_devtime ####
#----------------------------------

# trajectories
p_larvtim <- ggplot(data=larv_prop, aes(x=larv_devtime, y=prop, col=MismTreat, linetype=PhotoTreat))+
  scale_colour_manual(values=c("red4", "tomato2", "red", "orangered", "springgreen4", "steelblue2", "royalblue3", "blue4", "grey39", "black"))+
  geom_step(size=1)+
  scale_x_continuous(lim=c(25,55), breaks=seq(5,75,by=10))+ scale_y_continuous(breaks=seq(0,1,by=0.1))+
  labs(y="Proportion of larvae not yet pupated", x="Larval developmental time (days)")+
  theme(legend.title= element_blank(), legend.text=element_text(size=18), axis.title = element_text(size=18), 
        axis.text=element_text(size=16))
p_larvtim  
# ggsave(filename="_results/LarvDevTimeTraj.png", plot=p_larvtim, device="png", width=200, height=150, units="mm", dpi="print")

# Averages
devtime <- Rmisc::summarySE(d_tim, measurevar="larv_devtime", groupvars=c("MismTreat", "PhotoTreat"))
devtime$samplesize <- aggregate(CaterpillarID~PhotoTreat + MismTreat, d_tim, length)$CaterpillarID
devtime$pos <- ifelse(is.na(devtime$se)==T, 0, devtime$se) # position of sample size labels

p_devtime <- ggplot(data=devtime, aes(y=larv_devtime, x=MismTreat, col=PhotoTreat, fill=PhotoTreat)) +
  scale_colour_manual(values=c("grey27", "orangered2"))+ #"dodgerblue4"
  scale_fill_manual(values=c("grey27", "orangered2"))+
  geom_jitter(data=d_tim, alpha=0.3, size=3, height=0, width=0.25)+
  geom_point(size=5, shape=21, col="black")+
  geom_errorbar(aes(ymax = larv_devtime+se, ymin=larv_devtime-se), width=0.3) +
  geom_text(data=filter(devtime, PhotoTreat=="Changing"),aes(label=samplesize, y=larv_devtime+pos+0.5), col="black", size=4, fontface="bold")+
  geom_text(data=filter(devtime, PhotoTreat=="Constant"),aes(label=samplesize, y=larv_devtime-pos-0.5), col="black", size=4, fontface="bold")+
  labs(y="Larval development time (days)", x="Mismatch with oak budburst date (days)")+
  scale_x_continuous(breaks=seq(-4, 5, by=1)) + scale_y_continuous(lim=c(29,50),breaks=seq(20,50, by=5))+
  theme(axis.title = element_text(size=18), axis.text=element_text(size=16), 
        legend.text = element_text(size=16), legend.title=element_text(size=17))
p_devtime
# ggsave(filename="_results/LarvDevTime_raw.png", plot=p_devtime, device="png", width=200, height=150, units="mm", dpi="print")


#-----------------------------------
# Linear model ####
#----------------------------------
timlm1 <- lmer(larv_devtime ~ (MismTreat1 + MismTreat2)*PhotoTreat + (1|TubeID), data=d_tim)
anova1 <- anova(timlm1) # interactions not significant
anova1$mod <- "timlm1"

timlm2 <- update(timlm1, .~. -MismTreat1:PhotoTreat - MismTreat2:PhotoTreat)
anova2 <- anova(timlm2) # all fixed effects significant
anova2$mod <- "timlm2"

# Still there if remove -4 singularities?
timlm3 <- lmer(larv_devtime ~ MismTreat1 + MismTreat2 + PhotoTreat + (1|TubeID), data=filter(d_tim, MismTreat!=-4))
anova(timlm3) # yes

# Final model
timlm_final <- timlm2
summary(timlm_final)
timlm_res <- summary(timlm_final)$coefficients %>% as.data.frame

plot(timlm_final) #equal variance? ok
qqnorm(resid(timlm_final)) #normally distributed? okish
qqline(resid(timlm_final))

# write.csv(timlm_res, file="_results/output_LarvDevTime_lmer.csv", row.names=T)
# write.csv(rbind(anova1, anova2), file="_results/anova_LarvDevTime_lmer.csv", row.names=T)

# Get predictions ####
timlm.pred <- d_tim[!duplicated(d_tim[,c("TubeID", "Treatment")]),] # each replicate assigned same prediction, so remove duplicates
timlm.pred$pred <- predict(timlm_final, newdata=timlm.pred, type="response")
head(timlm.pred)

# Visualize predictions ####
pred <- Rmisc::summarySE(timlm.pred, measurevar="pred", groupvars=c("MismTreat", "PhotoTreat")) # significant effect of photoperiod, so show separate means
samplesize <- aggregate(CaterpillarID~MismTreat + PhotoTreat, data=d_tim, length) %>% arrange(MismTreat, PhotoTreat) %>% select(CaterpillarID)
pred$samplesize <- samplesize$CaterpillarID
rm(samplesize)

p_larvdev <- p_devtime + geom_line(data=filter(pred, PhotoTreat=="Constant"), aes(y=pred), size=1)+
  geom_smooth(data=filter(pred, PhotoTreat=="Changing"), aes(y=pred), se=F)
p_larvdev
# ggsave(filename="_results/LarvDevTime_wpred_rev.png", plot=p_larvdev, device="png", width=200, height=150, units="mm", dpi="print")

# Post-hoc test ####
summary(timlm2)
str(d_tim)
d_tim$TreatGroup <- as.factor(ifelse((d_tim$MismTreat1-5)<1, "BEFORE", # split treatment into hatching before, on, or after fitness peak groups
                           ifelse((d_tim$MismTreat1-5)>2, "AFTER", "ON")))
d_tim$TreatGroup <- factor(d_tim$TreatGroup, levels=c("ON", "BEFORE", "AFTER")) # day0 = ref level
table(d_tim$TreatGroup)
levels(d_tim$TreatGroup)

timlm_posth <- lmer(larv_devtime ~ TreatGroup + PhotoTreat + (1|TubeID), data=d_tim) # take final model and replace linear MismTreat terms for factor variable
anova(timlm_posth)
summary(timlm_posth)
# NB: if do not group mismatch+1 and mismatch+2 together do not find effect of TreatGroup ####

# emm_group <- emmeans::emmeans(timlm_posth, "TreatGroup", by=c("PhotoTreat")) #only compare treatments within PhotoTreat
emm_group <- emmeans::emmeans(timlm_posth, "TreatGroup") # don't have to do for PhotoTreat separately because final model does not include interaction effect
posthoc_res <- pairs(emm_group)
# write.csv(posthoc_res, file="_results/LarvDevTime_posthoc.csv", row.names=T)

rm(anova1, anova2, devtime, larv_prop, p_devtime, p_larvtim, pred, timlm1, timlm2, timlm3, timlm_res, timlm_posth, emm_group, posthoc_res)


#-----------------------------------
# Visualize adult data ####
#----------------------------------
d_adult <- d %>% 
  mutate(PhotoTreat=gsub("(\\w+)Day.+","\\1",Treatment), MismTreat=gsub("\\w+(Day.+)","\\1",Treatment), Event=1) %>% # all included have emerged
  mutate(Treatment=as.factor(Treatment), PhotoTreat=as.factor(ifelse(PhotoTreat=="Chang", "Changing", "Constant")), 
         MismTreatf=as.factor(MismTreat), TubeID=as.factor(TubeID)) %>%
  mutate(MismTreat=as.numeric(gsub("Day(.+)","\\1",MismTreat)),
    Treatment=factor(Treatment, levels=c("ChangDay-4", "ChangDay-3", "ChangDay-2", "ChangDay-1", "ChangDay0", "ChangDay+1", "ChangDay+2", "ChangDay+3", "ChangDay+4",
                                              "ChangDay+5",  "ConstDay-4", "ConstDay-2", "ConstDay0", "ConstDay+2", "ConstDay+4")), 
         MismTreatPlot=factor(MismTreatf, levels=c("Day-4", "Day-3", "Day-2", "Day-1", "Day0", "Day+1", "Day+2", "Day+3", "Day+4", "Day+5")),
         MismTreatf=factor(MismTreatf, levels=c("Day+1", "Day-4", "Day-3", "Day-2", "Day-1", "Day0", "Day+2", "Day+3", "Day+4", "Day+5"))) %>% # fitness peak = reference level
  mutate(larv_devtime=PupationAprilDay-HatchAprilDay, 
         pup_devtime=as.integer(as.Date(AdultNovDate, origin=as.Date("2021-10-31"))-as.Date(PupationAprilDay, origin=as.Date("2021-3-31"))), 
         CaterpillarID=as.factor(CaterpillarID), PupaWeight=PupaWeight*1000, AdultWeight=AdultWeight*1000) %>% 
  select(CaterpillarID, Treatment, PhotoTreat, MismTreat, MismTreatPlot, larv_devtime, PupaWeight, AdultNovDate, pup_devtime, AdultWeight, Sex) %>%
  arrange(PhotoTreat, MismTreatPlot) %>%
  filter(!is.na(AdultNovDate))
d_adult 


# Sex proportions ####
#-----------------------------------
sex <- table(d_adult$PhotoTreat, d_adult$MismTreatPlot, d_adult$Sex) %>% as.data.frame
colnames(sex) <- c("PhotoTreat", "MismTreat", "Sex", "Freq")
sex <- filter(sex, Freq>0) %>% droplevels() %>% arrange(PhotoTreat, Sex, MismTreat)
sex$N <- sex$Freq
sex$Pos <- c(0.85, 0.68, 0.6, 0.5, 0.35, 0.23, 0.13, 0.05, 0.85, 0.5, 0.18)

sex1 <- table(d_adult$PhotoTreat, d_adult$Sex) %>% as.data.frame
colnames(sex1) <- c("PhotoTreat", "Sex", "N")
sex1$Freq <- c(0.75, 0.85,0.3,0.35)
sex1

# Per Photoperiod
p_sex <- ggplot(data=sex, aes(y=Freq, x=PhotoTreat)) +
  #scale_colour_manual(values=c("black", "dodgerblue4"))+
  geom_bar(position="fill", stat="identity", size=1, alpha=0.6, aes(fill=Sex)) +
  geom_text(data = sex1, aes(label = N), size = 5, fontface="bold") +
  labs(y="Proportion")+
  theme(axis.title.y = element_text(size=18), axis.title.x=element_blank(), axis.text=element_text(size=16), 
        legend.text = element_text(size=16), legend.title=element_text(size=17))
p_sex
# ggsave(filename="_results/SexProp_raw.png", plot=p_sex, device="png", width=200, height=150, units="mm", dpi="print")

# Per Photoperiod + MismTreat
p_sex1 <- ggplot(data=sex, aes(y=Freq, x=PhotoTreat)) + 
  scale_colour_manual(values=c("grey45", "orangered4", "blue4", "white", "black"))+
  geom_bar(position="fill", stat="identity", size=1, alpha=0.6, aes(fill=Sex, col=MismTreat)) +
  geom_text(aes(y=Pos, label = N), size = 5, fontface="bold") +
  labs(y="Proportion")+
  theme(axis.title.y = element_text(size=18), axis.title.x=element_blank(), axis.text=element_text(size=16), 
        legend.text = element_text(size=16), legend.title=element_text(size=17),
        panel.background = element_rect(fill = 'lightsteelblue1'))
p_sex1
# ggsave(filename="_results/SexProp_raw_MismTr.png", plot=p_sex1, device="png", width=200, height=150, units="mm", dpi="print")

# Piecharts for MismTreat per Photoperiod ####
list <- list(filter(sex, Sex=="female" & PhotoTreat=="Changing"), 
             filter(sex, Sex=="female" & PhotoTreat=="Constant"), 
             filter(sex, Sex=="male" & PhotoTreat=="Constant"),
             filter(sex, Sex=="male" & PhotoTreat=="Changing"))

levels(list[[1]]$MismTreat)
cols <- data.frame(MismTreat=levels(list[[1]]$MismTreat), 
                   col=colorRampPalette(RColorBrewer::brewer.pal(5, "Paired"))(5)[1:5]) # make sure always same color used for each MismTreat

for(item in 1:length(list)){
  d <- list[[item]]
  
  file <- paste0("_results/PieChart_", unique(d$PhotoTreat), "_", unique(d$Sex) ,".png")
  print(file)
  
  plot <- ggplot(data=d, aes(x="", y=Freq, fill=MismTreat)) +
    scale_fill_manual(values=c(filter(cols, MismTreat %in% unique(d$MismTreat))$col)) +
    #scale_fill_brewer(palette="Paired") +
    geom_bar(stat="identity", width=1, col="black") +
    coord_polar("y", start=0) +
    theme_classic() +
    theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title=element_blank(),
          legend.text = element_text(size=20), legend.title = element_text(size=22),
          plot.title=element_text(size=25, face="bold"))
  png(file, width=300, height=250, units="mm", res=300)
  print(plot)
  dev.off()
}
rm(item, d, list, file)


# Development time ####
#-----------------------------------
adevtime <- Rmisc::summarySE(d_adult, measurevar="pup_devtime", groupvars=c("MismTreat", "PhotoTreat"))
adevtime <- adevtime %>% mutate(MismTreat=as.numeric(gsub("Day(.+)","\\1",MismTreatPlot)))
adevtime <- merge(devtime[,c("MismTreatPlot","PhotoTreat")], adevtime, all.x=TRUE, by=c("MismTreatPlot","PhotoTreat")) # want all treatments in there, even if no adults emerged
#devtime$pos <- ifelse(is.na(devtime$se)==T, 0, devtime$se) # position of sample size labels

p_adevtime <- ggplot(data=adevtime, aes(y=pup_devtime, x=MismTreat, col=PhotoTreat)) +
  scale_colour_manual(values=c("grey27", "dodgerblue4"))+
  geom_jitter(data=d_adult, size=4, alpha=0.7, height=0, width=0.3)+
  geom_text(data=filter(adevtime, PhotoTreat=="Changing"),aes(label=N, y=pup_devtime+2), col="black", size=4, fontface="bold")+
  geom_text(data=filter(adevtime, PhotoTreat=="Constant"),aes(label=N, y=pup_devtime+2), col="black", size=4, fontface="bold")+
  labs(y="Pupal development time (days)", x="Mismatch with oak budburst date (days)")+
  scale_x_continuous(lim=c(-4,5), breaks=seq(-4, 5, by=1)) + scale_y_continuous(lim=c(160,220),breaks=seq(150,220, by=10))+
  theme(axis.title = element_text(size=18), axis.text=element_text(size=16), 
        legend.text = element_text(size=16), legend.title=element_text(size=17))
p_adevtime
# ggsave(filename="_results/PupDevTime_raw.png", plot=p_adevtime, device="png", width=200, height=150, units="mm", dpi="print")



#---------------------------------------------------------------------------------------------------------------------------
# PupaWeight good proxy for fitness? ####
#-----------------------------------------

# Load data 
#-----------------------------------------
# NB: data not included in deposition ####

# Field data
dbField <- read_xlsx("_data/db_Field_WeightEggs.xlsx") # get data CatchYear2021 as well ####
colnames(dbField)[9] <- "AdultWeight"
table(dbField$AreaName) # Want at least 10 obs

dbField <- dbField %>% mutate(AdultWeight=AdultWeight*1000) %>% 
  filter(AreaName %in% c("Doorwerth", "Hoge Veluwe", "Oosterhout", "Warnsborn")) # only keep main areas
head(dbField)

# Experimental data
dbfood <- read_xlsx("_data/db_FeedExps_PupWeight_Eggs.xlsx") # 3 experiments: Margriet van Asch 2003, Lucia Salis 2013, Renske Jongen 2015

db1 <- read_xlsx("_data/db_MultGen1.xlsx", col_types="text")
db1 <- filter(db1, !is.na(PupaWeight) & (!is.na(AdultWeight) | !is.na(Eggs))) %>% filter(Sex!=2)
head(db1)
table(db1$ExperimentName)

db2 <- read_xlsx("_data/db_MultGen2.xlsx")
db2 <- filter(db2, !is.na(PupaWeight) & (!is.na(AdultWeight) | !is.na(Eggs)))
table(db2$ExperimentName) # what are these? also multigeneration experiments?
head(db2)

db <- rbind(db1[,colnames(db1) %in% colnames(db2)], db2[,colnames(db2) %in% colnames(db1)]) %>%
  mutate(PupaWeight=as.numeric(PupaWeight)*1000, AdultWeight=as.numeric(AdultWeight)*1000, Eggs=as.numeric(Eggs))
rm(db1, db2)

table(db$ExperimentName) # Why has FoodPhot2015 been dropped in Query?

# Add FoodPhot2015 (AdultWeights are missing in db)
dbfood_sub1 <- read_xlsx("_data/AllData_Access_PhotoperiodFoodExperiment2015.xlsx", sheet=2) # raw data sheet from Renske
dbfood_sub1 <- dbfood_sub1 %>% mutate(AdultWeight=as.numeric(`Adult weight(mg)`), PupaWeight=`Pupal weight(mg)`, YearCatch=2014, YearHatch=2016, ExperimentName="FoodPhot2015",
         Treatment=paste(`Treatment Photoperiod`, ifelse(`Treatment Food`=="High quality", "HQF", "LQF"), sep=" "), 
         CaterpillarID=NA, TubeID=NA, ClutchID=NA, Eggs=NA)
dbfood_sub1$CaterpillarID <- paste("Cat", seq(1,nrow(dbfood_sub1), 1), sep="") # make up a CaterpillarID for now
head(dbfood_sub1) # get AdultWeight

dbfood_sub2 <- filter(dbfood, ExperimentName=="FoodPhot2015") %>% mutate(PupaWeight=PupaWeight*1000, AdultWeight=NA)

db <- rbind(db, dbfood_sub1[,colnames(dbfood_sub1) %in% colnames(db)])
db <- rbind(db, dbfood_sub2[,colnames(dbfood_sub2) %in% colnames(db)]) # every individual added double now, because can't easily connect CaterpillarID to raw data!
rm(dbfood, dbfood_sub1, dbfood_sub2)

db$Treatment %>% unique # want to include all treatments? ####
head(db)

# What are treatments Exp2003 Margriet? L1, L2, L3, T1, T2, V1, V2
# Two food qualities: 0 days, +5 days = T and L, then -5=V
# Numbers = Replicate experiments?

# Treatment CT = control
# Treatment 0duplicate = control?
# Treatment CG = outside fieldshed?

#db_filt <- filter(db, Treatment %in% c("0 weeks", "control HQF", "control LQF")) # filter out treatments? Maybe test first if necessary
doubles <- db[duplicated(db[,c("CaterpillarID")]),] # multiple clutches of eggs for some caterpillars

adults <- db[!duplicated(db[,c("CaterpillarID")]),]
eggs <- aggregate(Eggs~CaterpillarID, data=db[!duplicated(db[c("ClutchID"),])], sum)

adults <- merge(adults[,!colnames(adults) %in% "Eggs"], eggs, by="CaterpillarID", all.x=T)
head(adults)


#-----------------------------------------
# Visualize and test relationships
#-----------------------------------------

# Field data
#-----------------------------------------
summary(dbField$Eggs) # never more than 505 Eggs
summary(dbField$AdultWeight) # some extreme outliers
cutoff <- mean(dbField$AdultWeight, na.rm=T) + sd(dbField$AdultWeight, na.rm=T)*3 # exclude extremes
table(dbField$AdultWeight<=cutoff)
filter(dbField, AdultWeight>cutoff) # probably entered as mg instead of g? doesn't look like it though

p_adulteggs <- ggplot(data=dbField, aes(x=AdultWeight, y=Eggs)) +
  #scale_color_manual(values=c("royalblue", "red3", "black"))+
  facet_wrap(~AreaName)+
  #facet_wrap(~YearCatch)+
  geom_point(size=2) + #aes(col=AreaName)
  geom_smooth(method=lm, se=F)+
  labs(x="Adult weight (mg)") +
  scale_y_continuous(lim=c(0,510))+ scale_x_continuous(lim=c(0,cutoff), breaks=seq(0,50, by=10))+
  theme(axis.title = element_text(size=18), axis.text=element_text(size=16), 
        strip.text = element_text(size=16))
p_adulteggs


# Test relationship
lm1 <- lm(Eggs ~ AdultWeight*AreaName + YearCatch, data=filter(dbField, AdultWeight<=cutoff)) # Year=numeric
car::Anova(lm1) # interaction with Area significant and YearCatch as well
summary(lm1) # for every 1mg increase in weight, adult lays 6 eggs extra
plot(lm1) # residuals are off, probably because zero inflated?


# Experimental data
#-----------------------------------------

# PupaWeight and AdultWeight
table(adults$ExperimentName) # Want at least 10 obs
adults_weight <- filter(adults, !(ExperimentName %in% c("CatFoodExp2021")) & AdultWeight>0)
# zeros in RN2009, probably missing values (because weight can't be zero if it's been measured)

p_pupadult <- ggplot(data=filter(adults, AdultWeight>0), aes(x=PupaWeight, y=AdultWeight)) +
  #scale_color_manual(values=c("royalblue", "red3", "black"))+
  facet_wrap(~ExperimentName)+
  geom_point(size=2)+ 
  geom_smooth(method=lm, se=F)+
  labs(x="Weight at pupation (mg)", y="Adult weight (mg)")+
  theme(axis.title = element_text(size=18), axis.text=element_text(size=16), 
        strip.text = element_text(size=16))
p_pupadult # clear relationship
# ggsave(filename="_results/Exp_Pup_Adult_Weight.png", plot=p_pupadult, device="png", width=200, height=150, units="mm", dpi="print")

# Test relationship
explm1a <- lm(AdultWeight ~ PupaWeight*ExperimentName, data=adults_weight) # Year=numeric
car::Anova(explm1a) # interaction significant
summary(explm1a) # relationship depends on Experiment,  why? Differences in treatments used maybe why pupa weight develops into adult weight in different ways?
plot(explm1a) #ok

# PupaWeight estimates per Experiment
est <- summary(explm1a)$coefficients %>% as.data.frame()
est$parameter <- rownames(est)
est <- est[grepl('PupaWeight:', rownames(est))==T,]
est <- est %>% mutate(PupaWeight=Estimate + explm1a$coefficients["PupaWeight"])
est[nrow(est)+1,] <- c(summary(explm1a)$coefficients[2,], "PupaWeight",explm1a$coefficients["PupaWeight"])
rownames(est) <- seq(1, nrow(est), by=1)
est
# write.csv(est, "_results/output_PupAdulWeight.csv", row.names=F)

# Eggs and PupaWeight
summary(adults$Eggs) # extreme outliers, keep same range as Field data?
table(adults$Eggs<510, useNA="always")
filter(adults, Eggs>=1000) # something off with how some data is entered in the database it seems for RN2009
filter(db, CaterpillarID=="2406") # For same CaterpillarID multiple PupaWeights, shouldn't be possible
filter(db, CaterpillarID=="3612") # This one from CatPhot2013 seems ok, just insane number of eggs
adults_eggs <- filter(adults, !is.na(Eggs))
table(adults_eggs$ExperimentName) # want at least 10 obs and exclude RN2009
adults_eggs <- filter(adults_eggs, !(ExperimentName %in% c("PupaeTemp2013", "RN2009")))

p_pupeggs <- ggplot(data=adults_eggs, aes(x=PupaWeight, y=Eggs)) + 
  #scale_color_manual(values=c("royalblue", "red3", "black"))+
  facet_wrap(~ExperimentName)+
  geom_point(size=2) + 
  #geom_smooth(method=lm, se=F)+
  labs(x="Weight at pupation (mg)") +
  theme(axis.title = element_text(size=18), axis.text=element_text(size=16), 
        strip.text = element_text(size=16))
p_pupeggs

# Test relationship
explm1 <- lm(Eggs ~ PupaWeight*ExperimentName, data=adults_eggs) # Year=numeric
car::Anova(explm1) # interaction not significant

explm2 <- update(explm1, .~. - PupaWeight:ExperimentName)
car::Anova(explm2)
summary(explm2) # for every 1mg increase in weight, adult lays ~7 eggs extra
plot(explm2) #ok

adults_eggs$pred <- predict(explm2, newdata=adults_eggs)

pred_pupegg <- p_pupeggs + geom_line(data=adults_eggs, aes(y=pred), size=1, col="blue")
pred_pupegg
# ggsave(filename="_results/Exp_PupWeight_Eggs.png", plot=pred_pupegg, device="png", width=200, height=150, units="mm", dpi="print")


sessionInfo() %>% capture.output(file="_src/env_CatFoodExp2021_analysis.txt")
