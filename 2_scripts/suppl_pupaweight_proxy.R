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

