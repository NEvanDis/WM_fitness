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
library(lme4)
library(lmerTest)


# Load data ####
#-----------------------------------
d <- read.csv("1_data/CatFood2021_deposit.csv")
head(d)

length(unique(d$TubeID)) # should be 22 mothers
table(d$Treatment) # photoperiod and mismatch treatment coded in one variable


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

rm(anova1, anova2, larv_prop, p_devtime, p_larvtim, pred, timlm1, timlm2, timlm3, timlm_res, timlm_posth, emm_group, posthoc_res)


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
         CaterpillarID=as.factor(CaterpillarID), PupaWeight=PupaWeight_ingrams*1000, AdultWeight=AdultWeight_ingrams*1000) %>% 
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
adevtime <- Rmisc::summarySE(d_adult, measurevar="pup_devtime", groupvars=c("MismTreatPlot", "PhotoTreat"))
adevtime <- adevtime %>% mutate(MismTreat=as.numeric(gsub("Day(.+)","\\1",MismTreatPlot)))
adevtime <- merge(devtime[,c("MismTreat","PhotoTreat")], adevtime, all.x=TRUE, by=c("MismTreat","PhotoTreat")) # want all treatments in there, even if no adults emerged
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


