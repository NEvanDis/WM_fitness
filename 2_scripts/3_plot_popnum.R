# Visualize population numbers long-term winter moth monitoring

# Open R project file in top folder

# Load required packages
#------------------------------------------
library(tidyverse)
library(cowplot) #remove for grey backgrounds
theme_set(theme_cowplot()) #change default ggplot2 theme to white background


# Load data ####
#-----------------------------------
# All data = per year per Area (Site and Tree info not included)
d_popnr <- read.csv("1_data/1994-2021_Field_PopNum_cor.csv") # produced by script <prep_FieldData.R>
timing <- read.csv("1_data/1994-2021_rel_timing_werror.csv") # produced by script <prep_FieldData.R>


# Visualize ####
#------------------------------------------

# Figure with 4 panels ####
pop_names <- data.frame(AreaShortName=c("HV", "OH", "DO", "WA"), Area=c("Hoge Veluwe", "Oosterhout", "Doorwerth", "Warnsborn"))
d_popnr <- merge(d_popnr, pop_names, by="AreaShortName", all.x=T)

p_popnr <- ggplot(d_popnr, aes(x=YearCatch, y=AllCor/10))+ # Number of females (NoFemale) = used for population models; pattern looks very similar to AllCor
  facet_wrap(~Area)+
  #geom_hline(yintercept=0)+ # add for timing graph below
  geom_line(aes(group=grp))+ # to not draw line for data gaps
  geom_point(size=2)+ #col="darkred"
  labs(y="Number of Adults")+
  #labs(y="Number of Females", x="Catch Year")+ #for Number of females plot
  scale_x_continuous(limits=c(1994, 2021), breaks=seq(1995,2020,by=5)) + 
  scale_y_continuous(lim=c(0,22), breaks=seq(0,35,by=5)) + # for plotting all
  #scale_y_continuous(lim=c(0,100), breaks=seq(0,100,by=25)) + # for plotting number of females
  theme(axis.title.y=element_text(size=20), axis.title.x=element_blank(), axis.text = element_text(size=15), 
        strip.text = element_text(size=20))+
  background_grid(major = "xy") #get grid lines
p_popnr
# ggsave(filename="results/PopNrs.png", plot=p_popnr , device="png", width=280, height=150, units="mm", dpi="print")
# ggsave(filename="results/FemNrs.png", plot=p_popnr , device="png", width=280, height=150, units="mm", dpi="print")


# Add amount of mismatch ####
timing1 <- merge(timing, d_popnr[,c("YearCatch","AreaShortName","grp")], # add grp info for plotting
                 by.x=c("YearHatch", "AreaShortName"), by.y=c("YearCatch", "AreaShortName"), all=T) %>% 
  arrange(AreaShortName, YearHatch)
timing1 <- merge(timing1, pop_names, by="AreaShortName")
summary(timing1$Mismatch.Tree)

timing1$Mismatch.Tree2 <- (timing1$Mismatch.Tree + 30) # turn into positive number for plotting
timing1[is.na(timing1$grp)==TRUE, "grp"] <- 1 # fix NA
head(timing1) # for HV 1995 only Mismatch.Area, no Mismatch.Tree available; same for WA 2013 (winter moth catch trees were no budburst scored)

ylabels <- c("","","", "", "0", "50", "100", "150", "200", "250", "300")
ylabels2 <- c("-30","-25", "-20", "-15", "-10", "-5", "0")

# Figure with 4 panels 
p_time <- p_popnr + geom_hline(yintercept=30, linetype="dashed", size=1) +
  geom_line(data=timing1, aes(x=YearHatch, y=Mismatch.Tree2, group=grp), col="grey50") +
  geom_point(data=timing1, aes(x=YearHatch, y=Mismatch.Tree2), size=2, col="grey50", alpha=0.9)+ #add mismatch, col="red"
  scale_y_continuous(sec.axis = sec_axis(~ ., breaks=seq(0,30,by=5), labels=ylabels2, name="Mismatch with oak budburst (days)\n"), breaks=seq(-20,30,by=5), labels=ylabels)
p_time
# ggsave(filename="results/PopNrs_wmismatch.png", plot=p_time , device="png", width=280, height=150, units="mm", dpi="print")
