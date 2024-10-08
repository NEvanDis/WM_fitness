# Visualize population numbers long-term winter moth monitoring

# Open R project file in top folder

# Load required packages
#------------------------------------------
library(tidyverse)


# Load data ####
#------------------------------------------
d_popnum <- read.csv("1_data/1994-2021_Field_PopNum_raw_deposit.csv") # Population numbers
d_trapsnr <- read.csv(file="1_data/1994-2021_Field_nr-traps_deposit.csv") # Trapping effort
timing <- read.csv("1_data/1994-2021_Field_timing_deposit.csv") # Seasonal timing


# Descriptives ####
#------------------------------------------
# Number of trees each year
aggregate(Traps~AreaShortName, data=d_trapsnr, min)
aggregate(Traps~AreaShortName, data=d_trapsnr, max)
aggregate(Traps~AreaShortName, data=d_trapsnr, mean)


# Prep data ####
#------------------------------------------

# Population numbers ####
pop_names <- unique(d_popnum$AreaShortName)
pops <- list(subset(d_popnum, AreaShortName==pop_names[1]),subset(d_popnum, AreaShortName==pop_names[2]), 
             subset(d_popnum, AreaShortName==pop_names[3]), subset(d_popnum, AreaShortName==pop_names[4]))
names(pops) <- pop_names
lapply(pops, head)

#-- Get number of adults per year, correcting for Trapping effort
popnr <- list()
for(pop in 1:length(pops)){ 
  d <- pops[[pop]]
  
  # Get number of adults per year
  nr <- aggregate(NoFemale ~YearCatch,data=d, sum)
  NoMale <- aggregate(NoMale~YearCatch,data=d, sum)
  nr <- merge(nr, NoMale, by="YearCatch", all=T)
  nr <- nr %>% mutate(NoFemale=ifelse(is.na(NoFemale)==T, 0, NoFemale), NoMale=ifelse(is.na(NoMale)==T, 0, NoMale)) %>%
    mutate(All=NoFemale+NoMale)
  nr
  
  # Correct for nr traps per year
  nr <- merge(nr, filter(d_trapsnr, AreaShortName==names(pops)[pop]), by="YearCatch", all.x=T)
  nr <- nr %>% mutate(NoFemCor=NoFemale/Traps, NoMalCor=NoMale/Traps, AllCor=All/Traps)
  
  nr$grp <- "1" # for plotting later to make sure data gaps are not connected  by lines
  grp <- 1 # initialize
  
  for(row in 2:nrow(nr)){
    grp <- ifelse((nr$YearCatch[row] - nr$YearCatch[row-1])==1, grp, grp+1)
    nr$grp[row] <- as.character(grp)
  }
  
  head(nr)
  
  popnr[[pop]] <- nr
  names(popnr)[pop] <- names(pops)[pop]
}

lapply(popnr, head)
popnr <- data.table::rbindlist(popnr, use.names=T)
head(popnr)

#-- Data gaps
table(popnr$AreaShortName, popnr$grp)
# Hoge Veluwe has a large data gap

# Warnsborn has a data gap of 1 year, but I don't think this is an actual data gap
# No mention in Fieldbook of WA, but probably only because no adults caught that year
# Otherwise would have mentioned no traps were placed that year!!!
popnr <- rbind(popnr,data.frame(YearCatch=2013, NoFemale=0, NoMale=0, All=0, AreaShortName="WA", Traps=NA,
                                NoFemCor=0, NoMalCor=0, AllCor=0, grp=1)) # add year
popnr <- popnr %>% mutate(grp=ifelse(AreaShortName=="WA", 1, grp)) %>% arrange(AreaShortName, YearCatch)
# write.csv(popnr, file="1_data/1994-2021_Field_PopNum_cor.csv", row.names=F)

rm(d, NoMale, nr, p, pop, grp, row) # clean-up


# Seasonal timing ####

#-- First take average per Tree, then calculate average per Area
timing <- timing %>% mutate(Tree_un=paste(AreaShortName, Site, Tree, sep="_"))
timing_tree <- Rmisc::summarySE(timing, measurevar="mismatch", groupvars=c("YearHatch", "Tree_un"), na.rm=TRUE) %>% 
  mutate(AreaShortName=gsub("(\\w+)_\\d+_\\d+","\\1",Tree_un))
timing_avg <- Rmisc::summarySE(timing_tree[,c("YearHatch", "AreaShortName", "mismatch")], measurevar="mismatch", groupvars=c("YearHatch", "AreaShortName"), na.rm=TRUE)
timing_avg <- timing_avg %>% mutate(Mismatch.Tree=mismatch, sd.Tree=sd) %>% select(YearHatch, AreaShortName, N, Mismatch.Tree, sd.Tree)
write.csv(timing_avg, file="1_data/1994-2021_rel_timing_werror.csv", row.names=FALSE)

head(timing_avg)
rm(timing_tree)
