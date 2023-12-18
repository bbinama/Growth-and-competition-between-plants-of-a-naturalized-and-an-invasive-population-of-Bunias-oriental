
#Packages
library(lme4)#fitting LMMs
library(DHARMa)#model checking
library(lsmeans)#least sqare means
library(car)#Anova cmd 
library(sjmisc)#plot
library(ggplot2)#plot
library(dplyr)  # For data manipulation
library(ggplot2)  # For data visualisation
library(car) # mutiple comparison and variance test
library(ggpubr)#plot scatter
library(MASS)#mean comparison
library(lattice)
library(tidyverse)
library(rstatix)
library(stats)
library(moonBook)
library(boot)
library(sjmisc)
library(moonBook)
library(ggiraph)
library(ggiraphExtra)
library(plyr)
library(ggeffects)
library(multcompView)
library(emmeans)
library(multcomp)
library(rcompanion)
library(lme4)
library(DHARMa)
library(ggplot2)
library(car)
library(ggpubr)
library(MASS)
library(glmmTMB)
library(dplyr)
library(ggpubr)
library(reshape2) 
library(tidyr)
library(survival)
library(permute)
library(lattice)
library(vegan)
library(DescTools)# DunneTest
library(pbkrtest)
library (patchwork)


setwd ("/Volumes/BBlaise/PhD/Data_CDRIVE/R Data/Competition") # for mac
#dat_1<-read.table(file="Single.txt", header=TRUE,sep="\t",dec = ".")
#Check if groups (single, two, three)was confounded with fertilization differences

#JENA
data_f<-read.table(file="Harvest1_JE_Fert.txt", header=TRUE,sep="\t",dec = ".")
str(data_f)
library("ggpubr")
ggboxplot(data_f, x = "Fert_Effect", y = "JE_Leaf_Count", 
          color = "Fert_Effect", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
        
          ylab = "JE_Leaf_Count", xlab = "Treatment")

#or

boxplot(JE_Leaf_Count ~Mixture_Code*Fert_Effect,data=data_f)

kruskal.test(JE_Leaf_Count ~ Fert_Effect, data = data_f)

#Kruskal-Wallis rank sum test

#data:  JE_Leaf_Count by Fert_Effect
#Kruskal-Wallis chi-squared = 13.449, df = 2, p-value = 0.001201

# Pairwise comparison
pairwise.wilcox.test(data_f$JE_Leaf_Count, data_f$Fert_Effect,
                     p.adjust.method = "BH")

#Pairwise comparisons using Wilcoxon rank sum test with continuity correction 

#data:  data_f$JE_Leaf_Count and data_f$Fert_Effect 

#    a20    b40   
#b40 0.0477 -     
#c60 0.0016 0.0477

#P value adjustment method: BH 
######################################

#HARVEST ONE COMPETITION IN TWO (JE VS GO)
#DRY BIOMASS
dat_1<-read.table(file="JE_GO_Competition_H1a.txt", header=TRUE,sep="\t",dec = ".")

str(dat_1)

#Use lmer with RMEL=F
M1 <- lmer(JE_Biomass ~ Population+(1 | Mixture_Code)+Initial_biomass, REML=F, data = dat_1)
M1a <- lmer(sqrt(JE_Biomass) ~ Population+(1 | Mixture_Code)+Initial_biomass, REML=F, data = dat_1)
M1b <- lmer(log(JE_Biomass) ~ Population+(1 | Mixture_Code)+Initial_biomass, REML=F, data = dat_1)

#Check themodel(DHARMarequired)
sim1a<-simulateResiduals(fittedModel= M1,n=1000)#fit better
plot(sim1a, quantreg= FALSE)

sim1a<-simulateResiduals(fittedModel= M1a,n=1000) #
plot(sim1a, quantreg= FALSE)

sim1a<-simulateResiduals(fittedModel= M1b,n=1000)
plot(sim1a, quantreg= FALSE)

summary(M1a)
#sqrt transformation is better
#Significant test

Anova(M1,type="III",test="Chisq")

##

#DRY BIOMASS THREE
dat_1a<-read.table(file="JE_GO_Competition_H1ab.txt", header=TRUE,sep="\t",dec = ".")

str(dat_1a)

#Use lmer with RMEL=F
M1k <- lmer(GO_Biomass ~ Population+(1 | Mixture_Code)+Initial_biomass, REML=F, data = dat_1a)
M1a <- lmer(sqrt(GO_Biomass) ~ Population+(1 | Mixture_Code)+Initial_biomass, REML=F, data = dat_1a)
M1b <- lmer(log(GO_Biomass) ~ Population+(1 | Mixture_Code)+Initial_biomass, REML=F, data = dat_1a)

#Check themodel(DHARMarequired)
sim1a<-simulateResiduals(fittedModel= M1,n=1000)#fit better
plot(sim1a, quantreg= FALSE)

sim1a<-simulateResiduals(fittedModel= M1a,n=1000) #
plot(sim1a, quantreg= FALSE)

sim1a<-simulateResiduals(fittedModel= M1b,n=1000)
plot(sim1a, quantreg= FALSE)

summary(M1a)
#sqrt transformation is better
#Significant test

Anova(M1k,type="III",test="Chisq")

#TOTAL NUMBER OF LEAVES
dat_1<-read.table(file="JE_GO_Competition_H1a.txt", header=TRUE,sep="\t",dec = ".")

str(dat_1)

#Use lmer with RMEL=F
M1 <- lmer(JE_Leaf_Count ~ Population+(1 | Mixture_Code)+Initial_biomass, REML=F, data = dat_1)
M1a <- lmer(sqrt(JE_Leaf_Count) ~ Population+(1 | Mixture_Code)+Initial_biomass, REML=F, data = dat_1)
M1b <- lmer(log(JE_Leaf_Count) ~ Population+(1 | Mixture_Code)+Initial_biomass, REML=F, data = dat_1)

#Check themodel(DHARMarequired)
sim1a<-simulateResiduals(fittedModel= M1,n=1000)#fit better
plot(sim1a, quantreg= FALSE)

sim1a<-simulateResiduals(fittedModel= M1a,n=1000) # FIT BETTER
plot(sim1a, quantreg= FALSE)

sim1a<-simulateResiduals(fittedModel= M1b,n=1000)
plot(sim1a, quantreg= FALSE)

summary(M1a)
#sqrt transformation is better
#Significant test

Anova(M1a,type="III",test="Chisq")

##

#LEAF COUNT THREE
dat_1a<-read.table(file="JE_GO_Competition_H1ab.txt", header=TRUE,sep="\t",dec = ".")

str(dat_1a)

#Use lmer with RMEL=F
M1k <- lmer(GO_Leaf_Count  ~ Population+(1 | Mixture_Code)+Initial_biomass, REML=F, data = dat_1a)
M1a <- lmer(sqrt(GO_Leaf_Count ) ~ Population+(1 | Mixture_Code)+Initial_biomass, REML=F, data = dat_1a)
M1b <- lmer(log(GO_Leaf_Count ) ~ Population+(1 | Mixture_Code)+Initial_biomass, REML=F, data = dat_1a)

#Check themodel(DHARMarequired)
sim1a<-simulateResiduals(fittedModel= M1k,n=1000)#fit better
plot(sim1a, quantreg= FALSE)

sim1a<-simulateResiduals(fittedModel= M1a,n=1000) #
plot(sim1a, quantreg= FALSE)

sim1a<-simulateResiduals(fittedModel= M1b,n=1000)
plot(sim1a, quantreg= FALSE)

summary(M1a)
#sqrt transformation is better
#Significant test

Anova(M1k,type="III",test="Chisq")



#LEAF LENGTH
dat_1<-read.table(file="JE_GO_Competition_H1a.txt", header=TRUE,sep="\t",dec = ".")

str(dat_1)

#Use lmer with RMEL=F
M1 <- lmer(JE_LLL ~ Population+(1 | Mixture_Code)+Initial_biomass, REML=F, data = dat_1)
M1a <- lmer(sqrt(JE_LLL) ~ Population+(1 | Mixture_Code)+Initial_biomass, REML=F, data = dat_1)
M1b <- lmer(log(JE_LLL) ~ Population+(1 | Mixture_Code)+Initial_biomass, REML=F, data = dat_1)

#Check themodel(DHARMarequired)
sim1a<-simulateResiduals(fittedModel= M1,n=1000)#
plot(sim1a, quantreg= FALSE)

sim1a<-simulateResiduals(fittedModel= M1a,n=1000) # FIT BETTER
plot(sim1a, quantreg= FALSE)

sim1a<-simulateResiduals(fittedModel= M1b,n=1000)
plot(sim1a, quantreg= FALSE)

summary(M1a)
#sqrt transformation is better
#Significant test

Anova(M1a,type="III",test="Chisq")

##

#LEAF LENGTH THREE
dat_1a<-read.table(file="JE_GO_Competition_H1ab.txt", header=TRUE,sep="\t",dec = ".")

str(dat_1a)

#Use lmer with RMEL=F
M1k <- lmer(GO_LLL  ~ Population+(1 | Mixture_Code)+Initial_biomass, REML=F, data = dat_1a)
M1a <- lmer(sqrt(GO_LLL ) ~ Population+(1 | Mixture_Code)+Initial_biomass, REML=F, data = dat_1a)
M1b <- lmer(log(GO_LLL ) ~ Population+(1 | Mixture_Code)+Initial_biomass, REML=F, data = dat_1a)

#Check themodel(DHARMarequired)
sim1a<-simulateResiduals(fittedModel= M1k,n=1000)#fit better
plot(sim1a, quantreg= FALSE)

sim1a<-simulateResiduals(fittedModel= M1a,n=1000) #
plot(sim1a, quantreg= FALSE)

sim1a<-simulateResiduals(fittedModel= M1b,n=1000)
plot(sim1a, quantreg= FALSE)

summary(M1a)
#sqrt transformation is better
#Significant test

Anova(M1k,type="III",test="Chisq")

## BOXPLOT FOR LEAF LENGTH
data_H<-read.table(file="JE_GO_Competition_H1Box.txt", header=TRUE,sep="\t",dec = ".")
str(data_H)

boxplot(JE_LLL  ~Population,data=data_H)

position<-c(1,2,4,5)
boxplot(JE_LLL ~Population,data=data_H,axes=F,at=position ,border="gray40",xlab= "Treatment", ylab="Length of the longest leaf (cm)",las=2, 
        col=c("#FFFF00","#0033FF","#FFFF00","#0033FF"))
axis(2,at=seq(15,35,5), seq(15,35,5), las=1)
axis(1, at=c(1,2,4,5), labels=c("GO2","JE2","GO3","JE3"), las=1)
box()



#add mean for Competition in two
Mixture_Code1<-as.factor(dat_1$Population)
str(dat_1)
M1a <- lmer(sqrt(JE_LLL) ~ Mixture_Code1+(1 | Mixture_Code)+Initial_biomass, REML=F, data = dat_1)														

lsm<-emmeans(M1a,~Mixture_Code1+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Mixture_Code1+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(1,2)
points(position,plot$response,pch=4,cex=2)

#add mean for GO incompetition with JE
Mixture_Code1<-as.factor(dat_1a$Population)
str(dat_1a)
M1k <- lmer(GO_LLL  ~ Mixture_Code1+(1 | Mixture_Code)+Initial_biomass, REML=F, data = dat_1a)														

lsm<-emmeans(M1k,~Mixture_Code1+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm
pan.emm <- emmeans(regrid(lsm), list(~Mixture_Code1+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(4,5)
points(position,plot$emmean,pch=4,cex=2)

###
#SPECIFIC LEAF AREA FOR TWO

dat_1<-read.table(file="JE_GO_Competition_H1a.txt", header=TRUE,sep="\t",dec = ".")

str(dat_1)

#Use lmer with RMEL=F
M1 <- lmer(JE_SLA  ~ Population+(1 | Mixture_Code)+Initial_biomass, REML=F, data = dat_1)
M1a <- lmer(sqrt(JE_SLA ) ~ Population+(1 | Mixture_Code)+Initial_biomass, REML=F, data = dat_1)
M1b <- lmer(log(JE_SLA ) ~ Population+(1 | Mixture_Code)+Initial_biomass, REML=F, data = dat_1)

#Check themodel(DHARMarequired)
sim1a<-simulateResiduals(fittedModel= M1,n=1000)# FIT BETTER
plot(sim1a, quantreg= FALSE)

sim1a<-simulateResiduals(fittedModel= M1a,n=1000) # 
plot(sim1a, quantreg= FALSE)

sim1a<-simulateResiduals(fittedModel= M1b,n=1000)
plot(sim1a, quantreg= FALSE)

summary(M1a)
#sqrt transformation is better
#Significant test

Anova(M1,type="III",test="Chisq")

##

# CARBON TO NITROGEN RATIO TWO
dat_1<-read.table(file="JE_GO_Competition_H1a.txt", header=TRUE,sep="\t",dec = ".")

str(dat_1)

#Use lmer with RMEL=F
M1 <- lmer(C_N_.  ~ Population+(1 | Mixture_Code)+Initial_biomass, REML=F, data = dat_1)
M1a <- lmer(sqrt(C_N_. ) ~ Population+(1 | Mixture_Code)+Initial_biomass, REML=F, data = dat_1)
M1b <- lmer(log(C_N_. ) ~ Population+(1 | Mixture_Code)+Initial_biomass, REML=F, data = dat_1)

#Check themodel(DHARMarequired)
sim1a<-simulateResiduals(fittedModel= M1,n=1000)# 
plot(sim1a, quantreg= FALSE)

sim1a<-simulateResiduals(fittedModel= M1a,n=1000) # 
plot(sim1a, quantreg= FALSE)

sim1a<-simulateResiduals(fittedModel= M1b,n=1000) #fit better
plot(sim1a, quantreg= FALSE)

summary(M1a)
#sqrt transformation is better
#Significant test

Anova(M1b,type="III",test="Chisq")

##
# NITROGEN TWO
dat_1<-read.table(file="JE_GO_Competition_H1a.txt", header=TRUE,sep="\t",dec = ".")

str(dat_1)

#Use lmer with RMEL=F
M1 <- lmer(N_.  ~ Population+(1 | Mixture_Code)+Initial_biomass, REML=F, data = dat_1)
M1a <- lmer(sqrt(N_. ) ~ Population+(1 | Mixture_Code)+Initial_biomass, REML=F, data = dat_1)
M1b <- lmer(log(N_. ) ~ Population+(1 | Mixture_Code)+Initial_biomass, REML=F, data = dat_1)

#Check themodel(DHARMarequired)
sim1a<-simulateResiduals(fittedModel= M1,n=1000)# 
plot(sim1a, quantreg= FALSE)

sim1a<-simulateResiduals(fittedModel= M1a,n=1000) # 
plot(sim1a, quantreg= FALSE)

sim1a<-simulateResiduals(fittedModel= M1b,n=1000) #fit better
plot(sim1a, quantreg= FALSE)

summary(M1a)
#sqrt transformation is better
#Significant test

Anova(M1b,type="III",test="Chisq")

#
# CARBON TWO
dat_1<-read.table(file="JE_GO_Competition_H1a.txt", header=TRUE,sep="\t",dec = ".")

str(dat_1)

#Use lmer with RMEL=F
M1 <- lmer(C_.  ~ Population+(1 | Mixture_Code)+Initial_biomass, REML=F, data = dat_1)
M1a <- lmer(sqrt(C_. ) ~ Population+(1 | Mixture_Code)+Initial_biomass, REML=F, data = dat_1)
M1b <- lmer(log(C_. ) ~ Population+(1 | Mixture_Code)+Initial_biomass, REML=F, data = dat_1)

#Check themodel(DHARMarequired)
sim1a<-simulateResiduals(fittedModel= M1,n=1000)# 
plot(sim1a, quantreg= FALSE)

sim1a<-simulateResiduals(fittedModel= M1a,n=1000) # 
plot(sim1a, quantreg= FALSE)

sim1a<-simulateResiduals(fittedModel= M1b,n=1000) #fit better
plot(sim1a, quantreg= FALSE)

summary(M1a)
#sqrt transformation is better
#Significant test

Anova(M1b,type="III",test="Chisq")
#LEAF WATER CONTENT TWO
dat_1<-read.table(file="JE_GO_Competition_H1a.txt", header=TRUE,sep="\t",dec = ".")

str(dat_1)

#Use lmer with RMEL=F
M1 <- lmer(JE_LWC ~ Population+(1 | Mixture_Code)+Initial_biomass, REML=F, data = dat_1)
M1a <- lmer(sqrt(JE_LWC) ~ Population+(1 | Mixture_Code)+Initial_biomass, REML=F, data = dat_1)
M1b <- lmer(log(JE_LWC) ~ Population+(1 | Mixture_Code)+Initial_biomass, REML=F, data = dat_1)

#Check themodel(DHARMarequired)
sim1a<-simulateResiduals(fittedModel= M1,n=1000)# fit better
plot(sim1a, quantreg= FALSE)

sim1a<-simulateResiduals(fittedModel= M1a,n=1000) # FIT BETTER
plot(sim1a, quantreg= FALSE)

sim1a<-simulateResiduals(fittedModel= M1b,n=1000)
plot(sim1a, quantreg= FALSE)

summary(M1a)
#sqrt transformation is better
#Significant test

Anova(M1,type="III",test="Chisq")

##

#WATER CONTENT THREE
dat_1a<-read.table(file="JE_GO_Competition_H1ab.txt", header=TRUE,sep="\t",dec = ".")

str(dat_1a)

#Use lmer with RMEL=F
M1k <- lmer(GO_LWC  ~ Population+(1 | Mixture_Code)+Initial_biomass, REML=F, data = dat_1a)
M1a <- lmer(sqrt(GO_LWC ) ~ Population+(1 | Mixture_Code)+Initial_biomass, REML=F, data = dat_1a)
M1b <- lmer(log(GO_LWC ) ~ Population+(1 | Mixture_Code)+Initial_biomass, REML=F, data = dat_1a)

#Check themodel(DHARMarequired)
sim1a<-simulateResiduals(fittedModel= M1k,n=1000)#fit better
plot(sim1a, quantreg= FALSE)

sim1a<-simulateResiduals(fittedModel= M1a,n=1000) #
plot(sim1a, quantreg= FALSE)

sim1a<-simulateResiduals(fittedModel= M1b,n=1000)
plot(sim1a, quantreg= FALSE)

summary(M1a)
#sqrt transformation is better
#Significant test

Anova(M1k,type="III",test="Chisq")

###
## BOXPLOT FOR LEAF WATER CONTENT
data_H<-read.table(file="JE_GO_Competition_H1Box.txt", header=TRUE,sep="\t",dec = ".")
str(data_H)

boxplot(JE_LWC  ~Population,data=data_H)

position<-c(1,2,4,5)
boxplot(JE_LWC ~Population,data=data_H,axes=F,at=position ,border="gray40",xlab= "Treatment", ylab="Water content (%)",las=2, 
        col=c("#FFFF00","#0033FF","#FFFF00","#0033FF"))
axis(2,at=seq(82,92,2), seq(82,92,2), las=1)
axis(1, at=c(1.5,4.5), labels=c(" Competition for two","Competition for three"), las=1)
box()



#add mean for Competition in two
Mixture_Code1<-as.factor(dat_1$Population)
str(dat_1)
M1 <- lmer(JE_LWC ~ Mixture_Code1+(1 | Mixture_Code)+Initial_biomass, REML=F, data = dat_1)													

lsm<-emmeans(M1,~Mixture_Code1+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Mixture_Code1+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(1,2)
points(position,plot$emmean,pch=4,cex=2)

#add mean for GO incompetition with JE
Mixture_Code1<-as.factor(dat_1a$Population)
str(dat_1a)
M1k <- lmer(GO_LWC  ~ Mixture_Code1+(1 | Mixture_Code)+Initial_biomass, REML=F, data = dat_1a)														

lsm<-emmeans(M1k,~Mixture_Code1+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm
pan.emm <- emmeans(regrid(lsm), list(~Mixture_Code1+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(4,5)
points(position,plot$emmean,pch=4,cex=2)


#######################################################
#SECOND HARVEST TWO

#DRY BIOMASS
dat_1<-read.table(file="JE_GO_Competition_H2a.txt", header=TRUE,sep="\t",dec = ".")

str(dat_1)

#Use lmer with RMEL=F
M1 <- lmer(JE_Biomass ~ Population+(1 | Mixture_Code)+Initial_biomass, REML=F, data = dat_1)
M1a <- lmer(sqrt(JE_Biomass) ~ Population+(1 | Mixture_Code)+Initial_biomass, REML=F, data = dat_1)
M1b <- lmer(log(JE_Biomass) ~ Population+(1 | Mixture_Code)+Initial_biomass, REML=F, data = dat_1)

#Check themodel(DHARMarequired)
sim1a<-simulateResiduals(fittedModel= M1,n=1000)#fit better
plot(sim1a, quantreg= FALSE)

sim1a<-simulateResiduals(fittedModel= M1a,n=1000) #
plot(sim1a, quantreg= FALSE)

sim1a<-simulateResiduals(fittedModel= M1b,n=1000)
plot(sim1a, quantreg= FALSE)

summary(M1a)
#sqrt transformation is better
#Significant test

Anova(M1,type="III",test="Chisq")

##

#DRY BIOMASS THREE
dat_1a<-read.table(file="JE_GO_Competition_H2b.txt", header=TRUE,sep="\t",dec = ".")

str(dat_1a)

#Use lmer with RMEL=F
M1k <- lmer(GO_Biomass ~ Population+(1 | Mixture_Code)+Initial_biomass, REML=F, data = dat_1a)
M1a <- lmer(sqrt(GO_Biomass) ~ Population+(1 | Mixture_Code)+Initial_biomass, REML=F, data = dat_1a)
M1b <- lmer(log(GO_Biomass) ~ Population+(1 | Mixture_Code)+Initial_biomass, REML=F, data = dat_1a)

#Check themodel(DHARMarequired)
sim1a<-simulateResiduals(fittedModel= M1,n=1000)#
plot(sim1a, quantreg= FALSE)

sim1a<-simulateResiduals(fittedModel= M1a,n=1000) # fit better
plot(sim1a, quantreg= FALSE)

sim1a<-simulateResiduals(fittedModel= M1b,n=1000)
plot(sim1a, quantreg= FALSE)

summary(M1a)
#sqrt transformation is better
#Significant test

Anova(M1a,type="III",test="Chisq")

## BOXPLOT FOR LEAF LENGTH
data_H<-read.table(file="JE_GO_Competition_H2Box.txt", header=TRUE,sep="\t",dec = ".")
str(data_H)

boxplot(JE_Biomass  ~Population,data=data_H)

position<-c(1,2,4,5)
boxplot(JE_Biomass ~Population,data=data_H,axes=F,at=position ,border="gray40",xlab= "Treatment", ylab="Aboveground biomass (g DW)",las=2, 
        col=c("#FFFF00","#0033FF","#FFFF00","#0033FF"))
axis(2,at=seq(0.0,1.2,0.2), seq(0.0,1.2,0.2), las=1)
axis(1, at=c(1.5,4.5), labels=c("Group of two","Group of three"), las=1)
box()



#add mean for Competition in two
Mixture_Code1<-as.factor(dat_1$Population)
str(dat_1)
M1a <- lmer(JE_Biomass ~ Mixture_Code1+(1 | Mixture_Code)+Initial_biomass, REML=F, data = dat_1)														

lsm<-emmeans(M1a,~Mixture_Code1+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Mixture_Code1+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(1,2)
points(position,plot$emmean,pch=4,cex=2)

#add mean for three incompetition 
Mixture_Code1<-as.factor(dat_1a$Population)
str(dat_1a)
M1k <- lmer(sqrt(GO_Biomass)  ~ Mixture_Code1+(1 | Mixture_Code)+Initial_biomass, REML=F, data = dat_1a)														

lsm<-emmeans(M1k,~Mixture_Code1+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm
pan.emm <- emmeans(regrid(lsm), list(~Mixture_Code1+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(4,5)
points(position,plot$response,pch=4,cex=2)


#TOTAL NUMBER OF LEAVES
dat_1<-read.table(file="JE_GO_Competition_H2a.txt", header=TRUE,sep="\t",dec = ".")

str(dat_1)

#Use lmer with RMEL=F
M1 <- lmer(JE_Leaf_Count ~ Population+(1 | Mixture_Code)+Initial_biomass, REML=F, data = dat_1)
M1a <- lmer(sqrt(JE_Leaf_Count) ~ Population+(1 | Mixture_Code)+Initial_biomass, REML=F, data = dat_1)
M1b <- lmer(log(JE_Leaf_Count) ~ Population+(1 | Mixture_Code)+Initial_biomass, REML=F, data = dat_1)

#Check themodel(DHARMarequired)
sim1a<-simulateResiduals(fittedModel= M1,n=1000)# fit better
plot(sim1a, quantreg= FALSE)

sim1a<-simulateResiduals(fittedModel= M1a,n=1000) # 
plot(sim1a, quantreg= FALSE)

sim1a<-simulateResiduals(fittedModel= M1b,n=1000)
plot(sim1a, quantreg= FALSE)

summary(M1a)
#sqrt transformation is better
#Significant test

Anova(M1a,type="III",test="Chisq")

##
## BOXPLOT FOR LEAF COUNT
data_H<-read.table(file="JE_GO_Competition_H2Box.txt", header=TRUE,sep="\t",dec = ".")
str(data_H)

boxplot(JE_Leaf_Count  ~Population,data=data_H)

position<-c(1,2,4,5)
boxplot(JE_Leaf_Count ~Population,data=data_H,axes=F,at=position ,border="gray40",xlab= "Competition treatment", ylab="Total number of leaves",las=2, 
        col=c("#FFFF00","#0033FF","#FFFF00","#0033FF"))
axis(2,at=seq(4,12,2), seq(4,12,2), las=1)
axis(1, at=c(1.5,4.5), labels=c("Group of two","Group of three"), las=1)
box()



#add mean for Competition in two
Mixture_Code1<-as.factor(dat_1$Population)
str(dat_1)
M1a <- lmer(JE_Leaf_Count ~ Mixture_Code1+(1 | Mixture_Code)+Initial_biomass, REML=F, data = dat_1)														

lsm<-emmeans(M1a,~Mixture_Code1+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Mixture_Code1+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(1,2)
points(position,plot$emmean,pch=4,cex=2)

#add mean for three incompetition 
Mixture_Code1<-as.factor(dat_1a$Population)
str(dat_1a)
M1k <- lmer(GO_Leaf_Count   ~ Mixture_Code1+(1 | Mixture_Code)+Initial_biomass, REML=F, data = dat_1a)														

lsm<-emmeans(M1k,~Mixture_Code1+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm
pan.emm <- emmeans(regrid(lsm), list(~Mixture_Code1+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(4,5)
points(position,plot$emmean,pch=4,cex=2)


#LEAF COUNT THREE
dat_1a<-read.table(file="JE_GO_Competition_H2b.txt", header=TRUE,sep="\t",dec = ".")

str(dat_1a)

#Use lmer with RMEL=F
M1k <- lmer(GO_Leaf_Count  ~ Population+(1 | Mixture_Code)+Initial_biomass, REML=F, data = dat_1a)
M1a <- lmer(sqrt(GO_Leaf_Count ) ~ Population+(1 | Mixture_Code)+Initial_biomass, REML=F, data = dat_1a)
M1b <- lmer(log(GO_Leaf_Count ) ~ Population+(1 | Mixture_Code)+Initial_biomass, REML=F, data = dat_1a)

#Check themodel(DHARMarequired)
sim1a<-simulateResiduals(fittedModel= M1k,n=1000)#fit better
plot(sim1a, quantreg= FALSE)

sim1a<-simulateResiduals(fittedModel= M1a,n=1000) 
plot(sim1a, quantreg= FALSE)

sim1a<-simulateResiduals(fittedModel= M1b,n=1000)
plot(sim1a, quantreg= FALSE)

summary(M1a)
# transformation is better
#Significant test

Anova(M1k,type="III",test="Chisq")



#LEAF LENGTH
dat_1<-read.table(file="JE_GO_Competition_H2a.txt", header=TRUE,sep="\t",dec = ".")

str(dat_1)

#Use lmer with RMEL=F
M1 <- lmer(JE_LLL ~ Population+(1 | Mixture_Code)+Initial_biomass, REML=F, data = dat_1)
M1a <- lmer(sqrt(JE_LLL) ~ Population+(1 | Mixture_Code)+Initial_biomass, REML=F, data = dat_1)
M1b <- lmer(log(JE_LLL) ~ Population+(1 | Mixture_Code)+Initial_biomass, REML=F, data = dat_1)

#Check themodel(DHARMarequired)
sim1a<-simulateResiduals(fittedModel= M1,n=1000)# fit better
plot(sim1a, quantreg= FALSE)

sim1a<-simulateResiduals(fittedModel= M1a,n=1000) #
plot(sim1a, quantreg= FALSE)

sim1a<-simulateResiduals(fittedModel= M1b,n=1000)
plot(sim1a, quantreg= FALSE)

summary(M1a)
# transformation is better
#Significant test

Anova(M1a,type="III",test="Chisq")

##

#LEAF LENGTH THREE
dat_1a<-read.table(file="JE_GO_Competition_H2b.txt", header=TRUE,sep="\t",dec = ".")

str(dat_1a)

#Use lmer with RMEL=F
M1k <- lmer(GO_LLL  ~ Population+(1 | Mixture_Code)+Initial_biomass, REML=F, data = dat_1a)
M1aa <- lmer(sqrt(GO_LLL ) ~ Population+(1 | Mixture_Code)+Initial_biomass, REML=F, data = dat_1a)
M1b <- lmer(log(GO_LLL ) ~ Population+(1 | Mixture_Code)+Initial_biomass, REML=F, data = dat_1a)

#Check themodel(DHARMarequired)
sim1a<-simulateResiduals(fittedModel= M1k,n=1000)#fit better
plot(sim1a, quantreg= FALSE)

sim1a<-simulateResiduals(fittedModel= M1a,n=1000) #
plot(sim1a, quantreg= FALSE)

sim1a<-simulateResiduals(fittedModel= M1b,n=1000)
plot(sim1a, quantreg= FALSE)

summary(M1a)
# transformation is better
#Significant test

Anova(M1aa,type="III",test="Chisq")

## BOXPLOT FOR LEAF LENGTH
data_H<-read.table(file="JE_GO_Competition_H2Box.txt", header=TRUE,sep="\t",dec = ".")
str(data_H)

boxplot(JE_LLL  ~Population,data=data_H)

position<-c(1,2,4,5)
boxplot(JE_LLL ~Population,data=data_H,axes=F,at=position ,border="gray40",xlab= "Treatment", ylab="Length of the longest leaf (cm)",las=2, 
        col=c("#FFFF00","#0033FF","#FFFF00","#0033FF"))
axis(2,at=seq(5,20,5), seq(5,20,5), las=1)
axis(1, at=c(1.5,4.5), labels=c("Competition for two","Competition for three"), las=1)
box()



#add mean for Competition in two
Mixture_Code1<-as.factor(dat_1$Population)
str(dat_1)
M1 <- lmer(JE_LLL ~ Mixture_Code1+(1 | Mixture_Code)+Initial_biomass, REML=F, data = dat_1)														

lsm<-emmeans(M1,~Mixture_Code1+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Mixture_Code1+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(1,2)
points(position,plot$emmean,pch=4,cex=2)

#add mean for GO incompetition with JE
Mixture_Code1<-as.factor(dat_1a$Population)
str(dat_1a)
M1aa <- lmer(sqrt(GO_LLL)  ~ Mixture_Code1+(1 | Mixture_Code)+Initial_biomass, REML=F, data = dat_1a)														

lsm<-emmeans(M1aa,~Mixture_Code1+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm
pan.emm <- emmeans(regrid(lsm), list(~Mixture_Code1+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(4,5)
points(position,plot$response,pch=4,cex=2)

###
#SPECIFIC LEAF AREA FOR TWO

dat_1<-read.table(file="JE_GO_Competition_H2a.txt", header=TRUE,sep="\t",dec = ".")

str(dat_1)

#Use lmer with RMEL=F
M1 <- lmer(SLA  ~ Population+(1 | Mixture_Code)+Initial_biomass, REML=F, data = dat_1)
M1a <- lmer(sqrt(SLA ) ~ Population+(1 | Mixture_Code)+Initial_biomass, REML=F, data = dat_1)
M1b <- lmer(log(SLA ) ~ Population+(1 | Mixture_Code)+Initial_biomass, REML=F, data = dat_1)

#Check themodel(DHARMarequired)
sim1a<-simulateResiduals(fittedModel= M1,n=1000)# 
plot(sim1a, quantreg= FALSE)

sim1a<-simulateResiduals(fittedModel= M1a,n=1000) # fit better
plot(sim1a, quantreg= FALSE)

sim1a<-simulateResiduals(fittedModel= M1b,n=1000)
plot(sim1a, quantreg= FALSE)

summary(M1a)
#sqrt transformation is better
#Significant test

Anova(M1a,type="III",test="Chisq")

##

# CARBON TO NITROGEN RATIO TWO
dat_1<-read.table(file="JE_GO_Competition_H2a.txt", header=TRUE,sep="\t",dec = ".")

str(dat_1)

#Use lmer with RMEL=F
M1 <- lmer(C_N.  ~ Population+(1 | Mixture_Code)+Initial_biomass, REML=F, data = dat_1)
M1a <- lmer(sqrt(C_N. ) ~ Population+(1 | Mixture_Code)+Initial_biomass, REML=F, data = dat_1)
M1b <- lmer(log(C_N. ) ~ Population+(1 | Mixture_Code)+Initial_biomass, REML=F, data = dat_1)

#Check themodel(DHARMarequired)
sim1a<-simulateResiduals(fittedModel= M1,n=1000)# 
plot(sim1a, quantreg= FALSE)

sim1a<-simulateResiduals(fittedModel= M1a,n=1000) # 
plot(sim1a, quantreg= FALSE)

sim1a<-simulateResiduals(fittedModel= M1b,n=1000) #fit better
plot(sim1a, quantreg= FALSE)

summary(M1a)
#log transformation is better
#Significant test

Anova(M1b,type="III",test="Chisq")

##
# NITROGEN TWO
dat_1<-read.table(file="JE_GO_Competition_H2a.txt", header=TRUE,sep="\t",dec = ".")

str(dat_1)

#Use lmer with RMEL=F
M1 <- lmer(N_.  ~ Population+(1 | Mixture_Code)+Initial_biomass, REML=F, data = dat_1)
M1a <- lmer(sqrt(N_. ) ~ Population+(1 | Mixture_Code)+Initial_biomass, REML=F, data = dat_1)
M1b <- lmer(log(N_. ) ~ Population+(1 | Mixture_Code)+Initial_biomass, REML=F, data = dat_1)

#Check themodel(DHARMarequired)
sim1a<-simulateResiduals(fittedModel= M1,n=1000)# 
plot(sim1a, quantreg= FALSE)

sim1a<-simulateResiduals(fittedModel= M1a,n=1000) # 
plot(sim1a, quantreg= FALSE)

sim1a<-simulateResiduals(fittedModel= M1b,n=1000) #fit better
plot(sim1a, quantreg= FALSE)

summary(M1a)
#log transformation is better
#Significant test

Anova(M1b,type="III",test="Chisq")

#
# CARBON TWO
dat_1<-read.table(file="JE_GO_Competition_H2a.txt", header=TRUE,sep="\t",dec = ".")

str(dat_1)

#Use lmer with RMEL=F
M1 <- lmer(C_.  ~ Population+(1 | Mixture_Code)+Initial_biomass, REML=F, data = dat_1)
M1a <- lmer(sqrt(C_. ) ~ Population+(1 | Mixture_Code)+Initial_biomass, REML=F, data = dat_1)
M1b <- lmer(log(C_. ) ~ Population+(1 | Mixture_Code)+Initial_biomass, REML=F, data = dat_1)

#Check themodel(DHARMarequired)
sim1a<-simulateResiduals(fittedModel= M1,n=1000)# fit better
plot(sim1a, quantreg= FALSE)

sim1a<-simulateResiduals(fittedModel= M1a,n=1000) # 
plot(sim1a, quantreg= FALSE)

sim1a<-simulateResiduals(fittedModel= M1b,n=1000) #
plot(sim1a, quantreg= FALSE)

summary(M1a)
#sqrt transformation is better
#Significant test

Anova(M1,type="III",test="Chisq")


## BOXPLOT FOR CARBON
dat_1<-read.table(file="JE_GO_Competition_H2a.txt", header=TRUE,sep="\t",dec = ".")
str(dat_1)

boxplot(C_.  ~Population,data=dat_1)

position<-c(1,3)
boxplot(C_. ~Population,data=dat_1,axes=F,at=position ,border="gray40",xlab= "Treatment", ylab="Carbon content (%)",las=2, 
        col=c("#FFFF00","#0033FF","#FFFF00","#0033FF"))
axis(2,at=seq(28,40,2), seq(28,40,2), las=1)
axis(1, at=c(2), labels=c("Competition for two"), las=1)
box()

#add mean for Competition in two
Mixture_Code1<-as.factor(dat_1$Population)
str(dat_1)
M1 <- lmer(C_. ~ Mixture_Code1+(1 | Mixture_Code)+Initial_biomass, REML=F, data = dat_1)													

lsm<-emmeans(M1,~Mixture_Code1+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Mixture_Code1+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(1,3)
points(position,plot$emmean,pch=4,cex=2)

#LEAF WATER CONTENT TWO
dat_1<-read.table(file="JE_GO_Competition_H2a.txt", header=TRUE,sep="\t",dec = ".")

str(dat_1)

#Use lmer with RMEL=F
M1 <- lmer(JE_LWC ~ Population+(1 | Mixture_Code)+Initial_biomass, REML=F, data = dat_1)
M1a <- lmer(sqrt(JE_LWC) ~ Population+(1 | Mixture_Code)+Initial_biomass, REML=F, data = dat_1)
M1b <- lmer(log(JE_LWC) ~ Population+(1 | Mixture_Code)+Initial_biomass, REML=F, data = dat_1)

#Check themodel(DHARMarequired)
sim1a<-simulateResiduals(fittedModel= M1,n=1000)# 
plot(sim1a, quantreg= FALSE)

sim1a<-simulateResiduals(fittedModel= M1a,n=1000) # 
plot(sim1a, quantreg= FALSE)

sim1a<-simulateResiduals(fittedModel= M1b,n=1000)
plot(sim1a, quantreg= FALSE)

summary(M1a)
#transformation is better
#Significant test

Anova(M1b,type="III",test="Chisq")

##

#WATER CONTENT THREE
dat_1a<-read.table(file="JE_GO_Competition_H2b.txt", header=TRUE,sep="\t",dec = ".")

str(dat_1a)

#Use lmer with RMEL=F
M1k <- lmer(LWC  ~ Population+(1 | Mixture_Code)+Initial_biomass, REML=F, data = dat_1a)
M1a <- lmer(sqrt(LWC ) ~ Population+(1 | Mixture_Code)+Initial_biomass, REML=F, data = dat_1a)
M1b <- lmer(log(LWC ) ~ Population+(1 | Mixture_Code)+Initial_biomass, REML=F, data = dat_1a)

#Check themodel(DHARMarequired)
sim1a<-simulateResiduals(fittedModel= M1k,n=1000)#fit better
plot(sim1a, quantreg= FALSE)

sim1a<-simulateResiduals(fittedModel= M1a,n=1000) #
plot(sim1a, quantreg= FALSE)

sim1a<-simulateResiduals(fittedModel= M1b,n=1000)
plot(sim1a, quantreg= FALSE)

summary(M1a)
#transformation is better
#Significant test

Anova(M1k,type="III",test="Chisq")

###
## BOXPLOT FOR LEAF WATER CONTENT
data_H<-read.table(file="JE_GO_Competition_H2Box.txt", header=TRUE,sep="\t",dec = ".")
str(data_H)

boxplot(JE_LWC  ~Population,data=data_H)

position<-c(1,2,4,5)
boxplot(JE_LWC ~Population,data=data_H,axes=F,at=position ,border="gray40",xlab= "Treatment", ylab="Water content (%)",las=2, 
        col=c("#FFFF00","#0033FF","#FFFF00","#0033FF"))
axis(2,at=seq(70,100,5), seq(70,100,5), las=1)
axis(1, at=c(1.5,4.5), labels=c(" Competition for two","Competition for three"), las=1)
box()



#add mean for Competition in two
Mixture_Code1<-as.factor(dat_1$Population)
str(dat_1)
M1 <- lmer(JE_LWC ~ Mixture_Code1+(1 | Mixture_Code)+Initial_biomass, REML=F, data = dat_1)													

lsm<-emmeans(M1,~Mixture_Code1+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Mixture_Code1+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(1,2)
points(position,plot$emmean,pch=4,cex=2)

#add mean for GO incompetition with JE
Mixture_Code1<-as.factor(dat_1a$Population)
str(dat_1a)
M1k <- lmer(LWC  ~ Mixture_Code1+(1 | Mixture_Code)+Initial_biomass, REML=F, data = dat_1a)														

lsm<-emmeans(M1k,~Mixture_Code1+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm
pan.emm <- emmeans(regrid(lsm), list(~Mixture_Code1+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(4,5)
points(position,plot$emmean,pch=4,cex=2)


#####################################


#PAIRWISE COMPARISON FIRST HARVEST







#GO
data_b<-read.table(file="Harvest1_GO_Fert.txt", header=TRUE,sep="\t",dec = ".")
str(data_b)

#or

boxplot(GO_Leaf_Count ~Fert_Effect,data=data_b)

kruskal.test(GO_Leaf_Count ~ Fert_Effect, data = data_b)

#Kruskal-Wallis rank sum test

#data:  GO_Leaf_Count by Fert_Effect
#Kruskal-Wallis chi-squared = 19.419, df = 2, p-value = 6.069e-05
# Pairwise comparison
pairwise.wilcox.test(data_b$GO_Leaf_Count, data_b$Fert_Effect,
                     p.adjust.method = "BH")


#Pairwise comparisons using Wilcoxon rank sum test with continuity correction 

#data:  data_b$GO_Leaf_Count and data_b$Fert_Effect 

#    a20    b40   
#b40 0.0759 -     
#c60 8e-05  0.0077

#P value adjustment method: BH 

#POPULATION
#single (JE VS GO)BIOMASS
data_a<-read.table(file="Single1.txt", header=TRUE,sep="\t",dec = ".")
str(data_a)

# use lm for JE incompetition GO  
str(data_a)

# Test homogeneity of variance
leveneTest(GO_Biomass~Population,data=data_a)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.25


#Check normality with shapiro test
shapiro.test(data_a$GO_Biomass) 

#P VALUE > 0.05 the assumption is valid, p =0.11

#transformation
Ms<-lm(GO_Biomass ~  Population +Initial_biomass,data=data_a)

Anova(Ms,type="III",test="F")

#Single (JE VS GO)GO_Leaf_Count

# Test homogeneity of variance
leveneTest(GO_Leaf_Count~Population,data=data_a)
#Assumption is valid if pvalue >.05
# assumption is not valid, p =0.02


#Check normality with shapiro test
shapiro.test(data_a$GO_Leaf_Count) 

#P VALUE > 0.05 the assumption is not valid, p =0.001

#transformation
Msa<-lm(GO_Leaf_Count ~  Population +Initial_biomass,data=data_a)
Msb<-lm(sqrt(GO_Leaf_Count) ~  Population +Initial_biomass,data=data_a)
Msc<-lm(log(GO_Leaf_Count) ~  Population +Initial_biomass,data=data_a)


plot(Msa) # 
hist(resid(Msa))

plot(Msb) # 
hist(resid(Msb)) #

plot(Msc) # 
hist(resid(Msc))
summary(m4)

#sqrt transformation fit better
Anova(Msc,type="III",test="F")

## Single (JE VS GO) GO_LLL

# Test homogeneity of variance
leveneTest(GO_LLL~Population,data=data_a)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.25


#Check normality with shapiro test
shapiro.test(data_a$GO_LLL) 

#P VALUE > 0.05 the assumption is valid, p =0.82

#transformation
Msa<-lm(GO_LLL ~  Population +Initial_biomass,data=data_a)

Anova(Msa,type="III",test="F")

#Single (JE VS GO) GO_LWC

# Test homogeneity of variance
leveneTest(GO_LWC~Population,data=data_a)
#Assumption is valid if pvalue >.05
# assumption is not valid, p =0.027


#Check normality with shapiro test
shapiro.test(data_a$GO_LWC) 

#P VALUE > 0.05 the assumption is valid, p =0.81

Msa<-lm(LWC ~  Population +Initial_biomass,data=data_a)

Msb<-lm(sqrt(LWC) ~  Population +Initial_biomass,data=data_a)
Msc<-lm(log(LWC) ~  Population +Initial_biomass,data=data_a)


plot(Msa) # 
hist(resid(Msa))

plot(Msb) # 
hist(resid(Msb)) #

plot(Msc) # 
hist(resid(Msc)) # fit better
summary(m4)

#transformation 
Anova(Msc,type="III",test="F")

#Single (JE vs GO) GO_SLA

# Test homogeneity of variance
leveneTest(GO_SLA~Population,data=data_a)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.201


#Check normality with shapiro test
shapiro.test(data_a$GO_SLA) 

#P VALUE > 0.05 the assumption is valid, p =0.71

Msa<-lm(GO_SLA ~  Population +Initial_biomass,data=data_a)

#transformation 
Anova(Msa,type="III",test="F")

## Single(JE VS GO) C/N
# Test homogeneity of variance
leveneTest(C_N.~Population,data=data_a)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.26


#Check normality with shapiro test
shapiro.test(data_a$C_N.) 

#P VALUE > 0.05 the assumption is NOT valid, p =0.022

Msa<-lm(C_N. ~  Population +Initial_biomass,data=data_a)

Msb<-lm(sqrt(C_N.) ~  Population +Initial_biomass,data=data_a)
Msc<-lm(log(C_N.) ~  Population +Initial_biomass,data=data_a)


plot(Msa) # 
hist(resid(Msa))

plot(Msb) # 
hist(resid(Msb)) #

plot(Msc) # 
hist(resid(Msc)) # fit better
summary(m4)

#transformation 
Anova(Msc,type="III",test="F")

## Single(JE VS GO) N
str(data_a)
# Test homogeneity of variance
leveneTest(N_.~Population,data=data_a)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.28


#Check normality with shapiro test
shapiro.test(data_a$N_.) 

#P VALUE > 0.05 the assumption is valid, p =0.56

Msa<-lm(N_. ~  Population +Initial_biomass,data=data_a)

Msb<-lm(sqrt(N_.) ~  Population +Initial_biomass,data=data_a)
Msc<-lm(log(N_.) ~  Population +Initial_biomass,data=data_a)


plot(Msa) # 
hist(resid(Msa)) # fit better

plot(Msb) # 
hist(resid(Msb)) #

plot(Msc) # 
hist(resid(Msc)) # 
summary(m4)

#transformation 
Anova(Msa,type="III",test="F")

## Single(JE VS GO) C
str(data_a)
# Test homogeneity of variance
leveneTest(C_.~Population,data=data_a)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.83


#Check normality with shapiro test
shapiro.test(data_a$C_.) 

#P VALUE > 0.05 the assumption is valid, p =0.17

Msa<-lm(C_. ~  Population +Initial_biomass,data=data_a)

Msb<-lm(sqrt(C_.) ~  Population +Initial_biomass,data=data_a)
Msc<-lm(log(C_.) ~  Population +Initial_biomass,data=data_a)


plot(Msa) # 
hist(resid(Msa))

plot(Msb) # 
hist(resid(Msb)) # fit better

plot(Msc) # 
hist(resid(Msc)) # 
summary(m4)

#transformation 
Anova(Msb,type="III",test="F")

##INTER1_TWO POPULATIONS
data_b<-read.table(file="Inter1_Two.txt", header=TRUE,sep="\t",dec = ".")
str(data_b)
#BIOMASS
# Test homogeneity of variance
leveneTest(Biomass~Population,data=data_b)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.89


#Check normality with shapiro test
shapiro.test(data_b$Biomass) 

#P VALUE > 0.05 the assumption is valid, p =0.74

Msa<-lm(Biomass ~  Population +Initial_biomass,data=data_b)

Msb<-lm(sqrt(Biomass) ~  Population +Initial_biomass,data=data_b)
Msc<-lm(log(Biomass) ~  Population +Initial_biomass,data=data_b)


plot(Msa) # 
hist(resid(Msa)) # fit better

plot(Msb) # 
hist(resid(Msb)) #

plot(Msc) # 
hist(resid(Msc)) # 
summary(m4)

#transformation 
Anova(Msa,type="III",test="F")

#TOTAL LEAVES
# Test homogeneity of variance
leveneTest(Leaf_Count~Population,data=data_b)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.44


#Check normality with shapiro test
shapiro.test(data_b$Leaf_Count) 

#P VALUE > 0.05 the assumption is not valid, p =0.001

Msa<-lm(Leaf_Count ~  Population +Initial_biomass,data=data_b)

Msb<-lm(sqrt(Leaf_Count) ~  Population +Initial_biomass,data=data_b)
Msc<-lm(log(Leaf_Count) ~  Population +Initial_biomass,data=data_b)


plot(Msa) # 
hist(resid(Msa))

plot(Msb) # 
hist(resid(Msb)) #

plot(Msc) # 
hist(resid(Msc)) # fit better
summary(m4)

#transformation 
Anova(Msc,type="III",test="F")

#LENGTH LONGEST LEAF

#Test homogeneity of variance
leveneTest(LLL ~Population,data=data_b)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.13


#Check normality with shapiro test
shapiro.test(data_b$LLL) 

#P VALUE > 0.05 the assumption is NOT valid, p =0.002

Msa<-lm(LLL ~  Population +Initial_biomass,data=data_b)

Msb<-lm(sqrt(LLL) ~  Population +Initial_biomass,data=data_b)
Msc<-lm(log(LLL) ~  Population +Initial_biomass,data=data_b)


plot(Msa) # 
hist(resid(Msa)) # fit better

plot(Msb) # 
hist(resid(Msb)) #

plot(Msc) # 
hist(resid(Msc)) # 
summary(m4)

#transformation 
Anova(Msa,type="III",test="F")

#SPECIFIC LEAF AREA
#Test homogeneity of variance
leveneTest(SLA ~Population,data=data_b)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.93


#Check normality with shapiro test
shapiro.test(data_b$SLA) 

#P VALUE > 0.05 the assumption is valid, p =0.19

Msa<-lm(SLA ~  Population +Initial_biomass,data=data_b)

Msb<-lm(sqrt(SLA) ~  Population +Initial_biomass,data=data_b)
Msc<-lm(log(SLA) ~  Population +Initial_biomass,data=data_b)


plot(Msa) # 
hist(resid(Msa))

plot(Msb) # 
hist(resid(Msb)) # fit better

plot(Msc) # 
hist(resid(Msc)) 
summary(m4)

#transformation 
Anova(Msb,type="III",test="F")

#LEAF WATER CONTENT
#Test homogeneity of variance
leveneTest(LWC ~Population,data=data_b)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.24


#Check normality with shapiro test
shapiro.test(data_b$LWC) 

#P VALUE > 0.05 the assumption is valid, p =0.65

Msa<-lm(LWC ~  Population +Initial_biomass,data=data_b)

Msb<-lm(sqrt(LWC) ~  Population +Initial_biomass,data=data_b)
Msc<-lm(log(LWC) ~  Population +Initial_biomass,data=data_b)


plot(Msa) # 
hist(resid(Msa))

plot(Msb) # 
hist(resid(Msb)) #

plot(Msc) # 
hist(resid(Msc)) # fit better
summary(m4)

#transformation 
Anova(Msc,type="III",test="F")

# CARBON TO NITROGEN CONTENT
#Test homogeneity of variance
leveneTest(C_N. ~Population,data=data_b)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.26


#Check normality with shapiro test
shapiro.test(data_b$C_N.) 

#P VALUE > 0.05 the assumption is NOT valid, p =0.014

Msa<-lm(C_N. ~  Population +Initial_biomass,data=data_b)

Msb<-lm(sqrt(C_N.) ~  Population +Initial_biomass,data=data_b)
Msc<-lm(log(C_N.) ~  Population +Initial_biomass,data=data_b)


plot(Msa) # 
hist(resid(Msa)) #

plot(Msb) # 
hist(resid(Msb)) #

plot(Msc) # 
hist(resid(Msc)) # 
summary(m4)

#transformation 
Anova(Msa,type="III",test="F")

## INTER 1(JE VS GO) N
str(data_b)
# Test homogeneity of variance
leveneTest(N_.~Population,data=data_b)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.58


#Check normality with shapiro test
shapiro.test(data_b$N_.) 

#P VALUE > 0.05 the assumption is NOT valid, p =0.0001

Msa<-lm(N_. ~  Population +Initial_biomass,data=data_b)

Msb<-lm(sqrt(N_.) ~  Population +Initial_biomass,data=data_b)
Msc<-lm(log(N_.) ~  Population +Initial_biomass,data=data_b)


plot(Msa) # 
hist(resid(Msa))

plot(Msb) # 
hist(resid(Msb)) # fit better

plot(Msc) # 
hist(resid(Msc)) # fit better
summary(m4)

#transformation 
Anova(Msc,type="III",test="F")

## INTER1 TWO (JE VS GO) C
str(data_b)
# Test homogeneity of variance
leveneTest(C_.~Population,data=data_b)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.59


#Check normality with shapiro test
shapiro.test(data_b$C_.) 

#P VALUE > 0.05 the assumption is NOT valid, p =0.0001

Msa<-lm(C_. ~  Population +Initial_biomass,data=data_b)

Msb<-lm(sqrt(C_.) ~  Population +Initial_biomass,data=data_b)
Msc<-lm(log(C_.) ~  Population +Initial_biomass,data=data_b)


plot(Msa) # 
hist(resid(Msa))

plot(Msb) # 
hist(resid(Msb)) #

plot(Msc) # 
hist(resid(Msc)) # fit better
summary(m4)

#transformation 
Anova(Msc,type="III",test="F")

##
#INTRA1_TWO
data_c<-read.table(file="Intra1_Two.txt", header=TRUE,sep="\t",dec = ".")
str(data_c)

#BIOMASS
# Test homogeneity of variance
leveneTest(Biomass~Population,data=data_c)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.88


#Check normality with shapiro test
shapiro.test(data_c$Biomass) 

#P VALUE > 0.05 the assumption is valid, p =0.98

Msa<-lm(Biomass ~  Population +Initial_biomass,data=data_c)

Msb<-lm(sqrt(Biomass) ~  Population +Initial_biomass,data=data_c)
Msc<-lm(log(Biomass) ~  Population +Initial_biomass,data=data_c)


plot(Msa) # 
hist(resid(Msa)) # fit better

plot(Msb) # 
hist(resid(Msb)) #

plot(Msc) # 
hist(resid(Msc)) # 
summary(m4)

#transformation 
Anova(Msa,type="III",test="F")
#Posthoc
Population1<-as.factor(data_c$Population)
str(data_c)
Msa<-lm(Biomass ~  Population1 +Initial_biomass,data=data_c)
gltt_Msa<-glht(Msa,linfct=mcp(Population1="Tukey"))
gltt_Msa
summary(gltt_Msa)
post<-cld(gltt_Msa, level = 0.05, decreasing = F)
post
boxplot(Biomass~Population, data=data_c)

#TOTAL LEAVES
# Test homogeneity of variance
leveneTest(Leaf_Count~Population,data=data_c)
#Assumption is valid if pvalue >.05
# assumption is valid, p =1


#Check normality with shapiro test
shapiro.test(data_c$Leaf_Count) 

#P VALUE > 0.05 the assumption is valid, p =0.058

Msa<-lm(Leaf_Count ~  Population +Initial_biomass,data=data_c)

Msb<-lm(sqrt(Leaf_Count) ~  Population +Initial_biomass,data=data_c)
Msc<-lm(log(Leaf_Count) ~  Population +Initial_biomass,data=data_c)


plot(Msa) # 
hist(resid(Msa))

plot(Msb) # 
hist(resid(Msb)) #

plot(Msc) # 
hist(resid(Msc)) # fit better
summary(m4)

#transformation 
Anova(Msc,type="III",test="F")

#mutliple comparison
Population1<-as.factor(data_c$Population)
str(data_c)
Msc<-lm(Leaf_Count ~  Population1 +Initial_biomass,data=data_c)
gltt_Msc<-glht(Msa,linfct=mcp(Population1="Tukey"))
gltt_Msc
summary(gltt_Msc)
post<-cld(gltt_Msc, level = 0.05, decreasing = F)
post
boxplot(Leaf_Count~Population, data=data_c)

#LENGTH LONGEST LEAF

#Test homogeneity of variance
leveneTest(LLL ~Population,data=data_c)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.64


#Check normality with shapiro test
shapiro.test(data_c$LLL) 

#P VALUE > 0.05 the assumption is NOT valid, p =0.09

Msa<-lm(LLL ~  Population +Initial_biomass,data=data_c)

Msb<-lm(sqrt(LLL) ~  Population +Initial_biomass,data=data_c)
Msc<-lm(log(LLL) ~  Population +Initial_biomass,data=data_c)


plot(Msa) # 
hist(resid(Msa)) # fit better

plot(Msb) # 
hist(resid(Msb)) #

plot(Msc) # 
hist(resid(Msc)) # 
summary(m4)

#transformation 
Anova(Msa,type="III",test="F")

#SPECIFIC LEAF AREA
#Test homogeneity of variance
leveneTest(SLA ~Population,data=data_c)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.301


#Check normality with shapiro test
shapiro.test(data_c$SLA) 

#P VALUE > 0.05 the assumption is valid, p =0.14

Msa<-lm(SLA ~  Population +Initial_biomass,data=data_c)

Msb<-lm(sqrt(SLA) ~  Population +Initial_biomass,data=data_c)
Msc<-lm(log(SLA) ~  Population +Initial_biomass,data=data_c)


plot(Msa) # 
hist(resid(Msa)) #fit better

plot(Msb) # 
hist(resid(Msb)) #

plot(Msc) # 
hist(resid(Msc)) 
summary(m4)

#transformation 
Anova(Msa,type="III",test="F")
#mutliple comparison
Population1<-as.factor(data_c$Population)
str(data_c)
Msc<-lm(SLA ~  Population1 +Initial_biomass,data=data_c)
gltt_Msc<-glht(Msa,linfct=mcp(Population1="Tukey"))
gltt_Msc
summary(gltt_Msc)
post<-cld(gltt_Msc, level = 0.05, decreasing = F)
post
boxplot(SLA~Population, data=data_c)

#LEAF WATER CONTENT
#Test homogeneity of variance
leveneTest(LWC ~Population,data=data_c)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.12


#Check normality with shapiro test
shapiro.test(data_c$LWC) 

#P VALUE > 0.05 the assumption is not valid, p =0.03

Msa<-lm(LWC ~  Population +Initial_biomass,data=data_c)

Msb<-lm(sqrt(LWC) ~  Population +Initial_biomass,data=data_c)
Msc<-lm(log(LWC) ~  Population +Initial_biomass,data=data_c)


plot(Msa) # 
hist(resid(Msa))

plot(Msb) # 
hist(resid(Msb)) #

plot(Msc) # 
hist(resid(Msc)) # fit better
summary(m4)

#transformation 
Anova(Msc,type="III",test="F")

# CARBON TO NITROGEN CONTENT
#Test homogeneity of variance
leveneTest(C_N. ~Population,data=data_c)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.37


#Check normality with shapiro test
shapiro.test(data_c$C_N.) 

#P VALUE > 0.05 the assumption is valid, p =0.25

Msa<-lm(C_N. ~  Population +Initial_biomass,data=data_c)

Msb<-lm(sqrt(C_N.) ~  Population +Initial_biomass,data=data_c)
Msc<-lm(log(C_N.) ~  Population +Initial_biomass,data=data_c)


plot(Msa) # 
hist(resid(Msa)) # fit better

plot(Msb) # 
hist(resid(Msb)) #

plot(Msc) # 
hist(resid(Msc)) # 
summary(m4)

#transformation 
Anova(Msa,type="III",test="F")

## INTRA1 TWO (JE VS GO) N
str(data_c)
# Test homogeneity of variance
leveneTest(N_.~Population,data=data_c)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.62


#Check normality with shapiro test
shapiro.test(data_c$N_.) 

#P VALUE > 0.05 the assumption is NOT valid, p =0.003

Msa<-lm(N_. ~  Population +Initial_biomass,data=data_c)

Msb<-lm(sqrt(N_.) ~  Population +Initial_biomass,data=data_c)
Msc<-lm(log(N_.) ~  Population +Initial_biomass,data=data_c)


plot(Msa) # 
hist(resid(Msa))

plot(Msb) # 
hist(resid(Msb)) #

plot(Msc) # 
hist(resid(Msc)) # fit better
summary(m4)

#transformation 
Anova(Msc,type="III",test="F")

## INTRA1 TWO (JE VS GO) C
str(data_c)
# Test homogeneity of variance
leveneTest(C_.~Population,data=data_c)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.65


#Check normality with shapiro test
shapiro.test(data_c$C_.) 

#P VALUE > 0.05 the assumption is NOT valid, p =0.0001

Msa<-lm(C_. ~  Population +Initial_biomass,data=data_c)

Msb<-lm(sqrt(C_.) ~  Population +Initial_biomass,data=data_c)
Msc<-lm(log(C_.) ~  Population +Initial_biomass,data=data_c)


plot(Msa) # 
hist(resid(Msa))

plot(Msb) # 
hist(resid(Msb)) # fit better

plot(Msc) # 
hist(resid(Msc)) # 
summary(m4)

#transformation 
Anova(Msb,type="III",test="F")

##
#INTER1_TWO-THREE

data_d<-read.table(file="Inter1_Two_One.txt", header=TRUE,sep="\t",dec = ".")
str(data_d)

#BIOMASS
# Test homogeneity of variance
leveneTest(Biomass~Population,data=data_d)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.14


#Check normality with shapiro test
shapiro.test(data_d$Biomass) 

#P VALUE > 0.05 the assumption is valid, p =0.61

Msa<-lm(Biomass ~  Population +Initial_biomass,data=data_d)

Msb<-lm(sqrt(Biomass) ~  Population +Initial_biomass,data=data_d)
Msc<-lm(log(Biomass) ~  Population +Initial_biomass,data=data_d)


plot(Msa) # 
hist(resid(Msa)) # fit better

plot(Msb) # 
hist(resid(Msb)) #

plot(Msc) # 
hist(resid(Msc)) # 
summary(m4)

#transformation 
Anova(Msa,type="III",test="F")
#Posthoc
Population1<-as.factor(data_d$Population)
str(data_d)
Msa<-lm(Biomass ~  Population1 +Initial_biomass,data=data_d)
gltt_Msa<-glht(Msa,linfct=mcp(Population1="Tukey"))
gltt_Msa
summary(gltt_Msa)
post<-cld(gltt_Msa, level = 0.05, decreasing = F)
post
boxplot(Biomass~Population, data=data_d)

#TOTAL LEAVES
# Test homogeneity of variance
leveneTest(Leaf_Count~Population,data=data_d)
#Assumption is valid if pvalue >.05
# assumption is valid, p =1


#Check normality with shapiro test
shapiro.test(data_d$Leaf_Count) 

#P VALUE > 0.05 the assumption is valid, p =0.07

Msa<-lm(Leaf_Count ~  Population +Initial_biomass,data=data_d)

Msb<-lm(sqrt(Leaf_Count) ~  Population +Initial_biomass,data=data_d)
Msc<-lm(log(Leaf_Count) ~  Population +Initial_biomass,data=data_d)


plot(Msa) # 
hist(resid(Msa)) #fit better

plot(Msb) # 
hist(resid(Msb)) #

plot(Msc) # 
hist(resid(Msc)) # 
summary(m4)

#transformation 
Anova(Msa,type="III",test="F")

#LENGTH LONGEST LEAF

#Test homogeneity of variance
leveneTest(LLL ~Population,data=data_d)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.30


#Check normality with shapiro test
shapiro.test(data_c$LLL) 

#P VALUE > 0.05 the assumption is valid, p =0.09

Msa<-lm(LLL ~  Population +Initial_biomass,data=data_d)

Msb<-lm(sqrt(LLL) ~  Population +Initial_biomass,data=data_d)
Msc<-lm(log(LLL) ~  Population +Initial_biomass,data=data_d)


plot(Msa) # 
hist(resid(Msa)) # fit better

plot(Msb) # 
hist(resid(Msb)) #

plot(Msc) # 
hist(resid(Msc)) # 
summary(m4)

#transformation 
Anova(Msa,type="III",test="F")


#LEAF WATER CONTENT
#Test homogeneity of variance
leveneTest(LWC ~Population,data=data_d)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.22


#Check normality with shapiro test
shapiro.test(data_d$LWC) 

#P VALUE > 0.05 the assumption is valid, p =0.33

Msa<-lm(LWC ~  Population +Initial_biomass,data=data_d)

Msb<-lm(sqrt(LWC) ~  Population +Initial_biomass,data=data_d)
Msc<-lm(log(LWC) ~  Population +Initial_biomass,data=data_d)


plot(Msa) # 
hist(resid(Msa)) #fit better

plot(Msb) # 
hist(resid(Msb)) #

plot(Msc) # 
hist(resid(Msc)) # 
summary(m4)

#transformation 
Anova(Msa,type="III",test="F")

##INTER1_One_TWO
data_e<-read.table(file="Inter1_One_Two.txt", header=TRUE,sep="\t",dec = ".")
str(data_e)

#BIOMASS
# Test homogeneity of variance
leveneTest(Biomass~Population,data=data_e)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.78


#Check normality with shapiro test
shapiro.test(data_e$Biomass) 

#P VALUE > 0.05 the assumption is valid, p =0.56

Msa<-lm(Biomass ~  Population +Initial_biomass,data=data_e)

Msb<-lm(sqrt(Biomass) ~  Population +Initial_biomass,data=data_e)
Msc<-lm(log(Biomass) ~  Population +Initial_biomass,data=data_e)


plot(Msa) # 
hist(resid(Msa)) # 

plot(Msb) # 
hist(resid(Msb)) # fit better

plot(Msc) # 
hist(resid(Msc)) # 
summary(m4)

#transformation 
Anova(Msb,type="III",test="F")

boxplot(Biomass~Population, data=data_e)

#TOTAL LEAVES
# Test homogeneity of variance
leveneTest(Leaf_Count~Population,data=data_e)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.09


#Check normality with shapiro test
shapiro.test(data_e$Leaf_Count) 

#P VALUE > 0.05 the assumption is valid, p =0.14

Msa<-lm(Leaf_Count ~  Population +Initial_biomass,data=data_e)

Msb<-lm(sqrt(Leaf_Count) ~  Population +Initial_biomass,data=data_e)
Msc<-lm(log(Leaf_Count) ~  Population +Initial_biomass,data=data_e)


plot(Msa) # 
hist(resid(Msa)) #

plot(Msb) # 
hist(resid(Msb)) #

plot(Msc) # 
hist(resid(Msc)) # fit better
summary(m4)

#transformation 
Anova(Msc,type="III",test="F")

#LENGTH LONGEST LEAF

#Test homogeneity of variance
leveneTest(LLL ~Population,data=data_e)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.33


#Check normality with shapiro test
shapiro.test(data_e$LLL) 

#P VALUE > 0.05 the assumption is valid, p =0.059

Msa<-lm(LLL ~  Population +Initial_biomass,data=data_e)

Msb<-lm(sqrt(LLL) ~  Population +Initial_biomass,data=data_e)
Msc<-lm(log(LLL) ~  Population +Initial_biomass,data=data_e)


plot(Msa) # 
hist(resid(Msa)) # fit better

plot(Msb) # 
hist(resid(Msb)) #

plot(Msc) # 
hist(resid(Msc)) # 
summary(m4)

#transformation 
Anova(Msa,type="III",test="F")


#LEAF WATER CONTENT
#Test homogeneity of variance
leveneTest(LWC ~Population,data=data_e)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.16


#Check normality with shapiro test
shapiro.test(data_e$LWC) 

#P VALUE > 0.05 the assumption is not valid, p =0.02

Msa<-lm(LWC ~  Population +Initial_biomass,data=data_e)

Msb<-lm(sqrt(LWC) ~  Population +Initial_biomass,data=data_e)
Msc<-lm(log(LWC) ~  Population +Initial_biomass,data=data_e)


plot(Msa) # 
hist(resid(Msa)) #

plot(Msb) # 
hist(resid(Msb)) # fit better

plot(Msc) # 
hist(resid(Msc)) # 
summary(m4)

#transformation 
Anova(Msb,type="III",test="F")

#Intra1_Three (JE VS GO)
data_i3<-read.table(file="Intra1_Three.txt", header=TRUE,sep="\t",dec = ".")
str(data_i3)
## BIOMASS
# Test homogeneity of variance
leveneTest(Biomass~Population,data=data_i3)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.99


#Check normality with shapiro test
shapiro.test(data_i3$Biomass) 

#P VALUE > 0.05 the assumption is valid, p =0.052

Mi3a<-lm(Biomass ~  Population +Initial_biomass,data=data_i3)

Mi3b<-lm(sqrt(Biomass) ~  Population +Initial_biomass,data=data_i3)
Mi3c<-lm(log(Biomass) ~  Population +Initial_biomass,data=data_i3)


plot(Mi3a) # 
hist(resid(Mi3a))

plot(Mi3b) # 
hist(resid(Mi3b)) #

plot(Mi3c) # 
hist(resid(Mi3c)) # fit better
summary(Mi3c)

#transformation 
Anova(Mi3c,type="III",test="F")

boxplot(Biomass ~Population,data=data_i3)
## TOTAL NUMBER OF LEAVES
# Test homogeneity of variance
leveneTest(Leaf_Count~Population,data=data_i3)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.16


#Check normality with shapiro test
shapiro.test(data_i3$Leaf_Count) 

#P VALUE > 0.05 the assumption is valid, p =0.34

Msa<-lm(Leaf_Count ~  Population +Initial_biomass,data=data_i3)

Msb<-lm(sqrt(Leaf_Count) ~  Population +Initial_biomass,data=data_i3)
Msc<-lm(log(Leaf_Count) ~  Population +Initial_biomass,data=data_i3)


plot(Msa) # 
hist(resid(Msa))

plot(Msb) # 
hist(resid(Msb)) #

plot(Msc) # 
hist(resid(Msc)) # fit better
summary(m4)

#transformation 
Anova(Msc,type="III",test="F")

#LEAF LONGEST LEAF
# Test homogeneity of variance
leveneTest(LLL~Population,data=data_i3)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.68


#Check normality with shapiro test
shapiro.test(data_i3$LLL ) 

#P VALUE > 0.05 the assumption is valid, p =0.28

Msa<-lm(LLL  ~  Population +Initial_biomass,data=data_i3)

Msb<-lm(sqrt(LLL) ~  Population +Initial_biomass,data=data_i3)
Msc<-lm(log(LLL) ~  Population +Initial_biomass,data=data_i3)


plot(Msa) # 
hist(resid(Msa))

plot(Msb) # 
hist(resid(Msb)) #

plot(Msc) # 
hist(resid(Msc)) # fit better
summary(m4)

#transformation 
Anova(Msc,type="III",test="F")

##LEAF WATER CONTENT
#Test homogeneity of variance
leveneTest(LWC~Population,data=data_i3)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.85


#Check normality with shapiro test
shapiro.test(data_i3$LWC ) 

#P VALUE > 0.05 the assumption is valid, p =0.11

Msa<-lm(LWC  ~  Population +Initial_biomass,data=data_i3)

Msb<-lm(sqrt(LWC) ~  Population +Initial_biomass,data=data_i3)
Msc<-lm(log(LWC) ~  Population +Initial_biomass,data=data_i3)


plot(Msa) # 
hist(resid(Msa))

plot(Msb) # 
hist(resid(Msb)) #

plot(Msc) # 
hist(resid(Msc)) # fit better
summary(m4)

#transformation 
Anova(Msc,type="III",test="F")

## SINGLE TO COMPETITION JENA
data_f<-read.table(file="Harvest1_JE_Fert.txt", header=TRUE,sep="\t",dec = ".")
str(data_f)

#ABOVEGROUND BIOMASS

boxplot(JE_Biomass ~Mixture_Code ,data=data_f)
#Use lmer with RMEL=F
Mixture_Code1<-as.factor(data_f$Mixture_Code)
str(data_f)

# Test homogeneity of variance
leveneTest(JE_Biomass~Mixture_Code1,data=data_f)
#Assumption is valid if pvalue >.05
# assumption is not valid, p =0.003


#Check normality with shapiro test
shapiro.test(data_f$JE_Biomass) 

#P VALUE > 0.05 the assumption is not valid, p =0.0004

#transformation
M2a<-lm(JE_Biomass ~  Mixture_Code1 +Initial_biomass,data=data_f)
M2b<-lm(sqrt(JE_Biomass) ~  Mixture_Code1 +Initial_biomass,data=data_f)
M2c<-lm(log(JE_Biomass) ~  Mixture_Code1 +Initial_biomass,data=data_f)


plot(M2a)  
hist(resid(M2a)) 

plot(M2b)  
hist(resid(M2b)) #fit better

plot(M2c) # 
hist(resid(M2c)) # 

#SQRT transformation fit better
Anova(M2b,type="III",test="F")
##
DunnettTest(sqrt(JE_Biomass) ~  Mixture_Code1,data=data_f)
#multiple comparaison

gltt_M2b<-glht(M2b,linfct=mcp(Mixture_Code1="Tukey"))
gltt_M2b
summary(gltt_M2b)
post<-cld(gltt_M2b, level = 0.05, decreasing = T)
post

#DIFFERENCE WITHIN REPLACEMENT SERIES GO incompetition JE
data_g<-read.table(file="Harvest1_GO_Fert.txt", header=TRUE,sep="\t",dec = ".")
Mixture_Code1<-as.factor(data_g$Mixture_Code)
str(data_g)

# Test homogeneity of variance
leveneTest(GO_Biomass~Mixture_Code1,data=data_g)
#Assumption is valid if pvalue >.05
# assumption is not valid, p =0.001


#Check normality with shapiro test
shapiro.test(data_g$GO_Biomass) 

#P VALUE > 0.05 the assumption is not valid, p =0.004

#transformation
M2a<-lm(GO_Biomass ~  Mixture_Code1 +Initial_biomass,data=data_g)
M2b<-lm(sqrt(GO_Biomass) ~  Mixture_Code1 +Initial_biomass,data=data_g)
M2c<-lm(log(GO_Biomass) ~  Mixture_Code1 +Initial_biomass,data=data_g)


plot(M2a)  
hist(resid(M2a)) #fit better

plot(M2b)  
hist(resid(M2b)) #

plot(M2c) # 
hist(resid(M2c)) # 

# transformation fit better
Anova(M2b,type="III",test="F")
##
DunnettTest(GO_Biomass ~  Mixture_Code1,data=data_g)
#multiple comparaison

gltt_M2a<-glht(M2b,linfct=mcp(Mixture_Code1="Tukey"))
gltt_M2a
summary(gltt_M2a)
post<-cld(gltt_M2a, level = 0.05, decreasing  = T)
post

#Data for Box_plot 
data_box<-read.table(file="Harvest1_JE_GO_Fert.txt", header=TRUE,sep="\t",dec = ".")
str(data_box)

boxplot(Biomass ~Mixture_Code,data=data_box)


position<-c(1,2,3,4,5,6,8,9,10,11,12,13)
boxplot(Biomass ~Mixture_Code,data=data_box,axes=F,at=position ,border="gray40",xlab= "Population", ylab="Aboveground biomass (g DW)",las=2, 
        col=c("#FFFF00","#FFFF00","#FFFF00","#FFFF00","#FFFF00","#FFFF00","#0033FF","#0033FF","#0033FF","#0033FF","#0033FF","#0033FF"))
axis(2,at=seq(0,5,1), seq(0,5,1), las=1)
axis(1, at=c(1,2,3,4,5,6,8,9,10,11,12,13), labels=c("0:1","1:1","0:2","2:1","1:2","0:3","0:1","1:1","0:2","2:1","1:2","0:3"), las=1)
box()

#add overall mean for JE and GO
#model predicted means

lsm<-emmeans(M1b,~Population+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Population+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(6.5,7.5)
points(position,plot$response,pch=4,cex=2)

#add mean for JE incompetition with GO
Mixture_Code1<-as.factor(data_f$Mixture_Code)
str(data_f)
M2b<-lm(sqrt(JE_Biomass) ~  Mixture_Code1 +Initial_biomass,data=data_f)

lsm<-emmeans(M2b,~Mixture_Code1+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Mixture_Code1+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(8,9,10,11,12,13)
points(position,plot$response,pch=4,cex=2)

#add mean for GO incompetition with JE
Mixture_Code1<-as.factor(data_g$Mixture_Code)
str(data_g)
M2a<-lm(GO_Biomass ~  Mixture_Code1 +Initial_biomass,data=data_g)

lsm<-emmeans(M2a,~Mixture_Code1+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm
pan.emm <- emmeans(regrid(lsm), list(~Mixture_Code1+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(1,2,3,4,5,6)
points(position,plot$emmean,pch=4,cex=2)

 #TOTAL NUMBER OF LEAVES
data_f<-read.table(file="Harvest1_JE_Fert.txt", header=TRUE,sep="\t",dec = ".")
boxplot(JE_Leaf_Count ~Mixture_Code ,data=data_f)
#Use lmer with RMEL=F
Mixture_Code1<-as.factor(data_f$Mixture_Code)
str(data_f)

# Test homogeneity of variance
leveneTest(JE_Leaf_Count~Mixture_Code1,data=data_f)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.11


#Check normality with shapiro test
shapiro.test(data_f$JE_Leaf_Count) 

#P VALUE > 0.05 the assumption is not valid, p =0.0002

#transformation
M2a<-lm(JE_Leaf_Count ~  Mixture_Code1 +Initial_biomass,data=data_f)
M2b<-lm(sqrt(JE_Leaf_Count) ~  Mixture_Code1 +Initial_biomass,data=data_f)
M2c<-lm(log(JE_Leaf_Count) ~  Mixture_Code1 +Initial_biomass,data=data_f)


plot(M2a)  
hist(resid(M2a)) 

plot(M2b)  
hist(resid(M2b)) #

plot(M2c) # 
hist(resid(M2c)) # fit better

#SQRT transformation fit better
Anova(M2c,type="III",test="F")
##
DunnettTest(log(JE_Leaf_Count) ~  Mixture_Code1,data=data_f)
#multiple comparaison

gltt_M2b<-glht(M2c,linfct=mcp(Mixture_Code1="Tukey"))
gltt_M2b
summary(gltt_M2b)
post<-cld(gltt_M2b, level = 0.05, decreasing = T)
post

#GO incompetition JE
data_g<-read.table(file="Harvest1_GO_Fert.txt", header=TRUE,sep="\t",dec = ".")
Mixture_Code1<-as.factor(data_g$Mixture_Code)
str(data_g)

# Test homogeneity of variance
leveneTest(GO_Leaf_Count~Mixture_Code1,data=data_g)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.46


#Check normality with shapiro test
shapiro.test(data_g$GO_Leaf_Count) 

#P VALUE > 0.05 the assumption is not valid, p =0.02

#transformation
M3a<-lm(GO_Leaf_Count ~  Mixture_Code1 +Initial_biomass,data=data_g)
M3b<-lm(sqrt(GO_Leaf_Count) ~  Mixture_Code1 +Initial_biomass,data=data_g)
M3c<-lm(log(GO_Leaf_Count) ~  Mixture_Code1 +Initial_biomass,data=data_g)


plot(M3a)  
hist(resid(M3a)) #

plot(M3b)  
hist(resid(M3b)) #

plot(M3c) # 
hist(resid(M3c)) # fit better

# transformation fit better
Anova(M3c,type="III",test="F")
##
DunnettTest(log(GO_Leaf_Count) ~  Mixture_Code1,data=data_g)
#multiple comparaison

gltt_M2a<-glht(M3c,linfct=mcp(Mixture_Code1="Tukey"))
gltt_M2a
summary(gltt_M2a)
post<-cld(gltt_M2a, level = 0.05, decreasing = T)
post

#Data for Box_plot 
data_box<-read.table(file="Harvest1_JE_GO_Fert.txt", header=TRUE,sep="\t",dec = ".")
str(data_box)

boxplot(Leaf_Count ~Mixture_Code,data=data_box)


position<-c(1,2,3,4,5,6,8,9,10,11,12,13)
boxplot(Leaf_Count ~Mixture_Code,data=data_box,axes=F,at=position ,border="gray40",xlab= "Population", ylab="Total number of leaves",las=2, 
        col=c("#FFFF00","#FFFF00","#FFFF00","#FFFF00","#FFFF00","#FFFF00","#0033FF","#0033FF","#0033FF","#0033FF","#0033FF","#0033FF"))
axis(2,at=seq(6,14,2), seq(6,14,2), las=1)
axis(1, at=c(1,2,3,4,5,6,8,9,10,11,12,13), labels=c("0:1","1:1","0:2","2:1","1:2","0:3","0:1","1:1","0:2","2:1","1:2","0:3"), las=1)
box()


#add mean for JE incompetition with GO
Mixture_Code1<-as.factor(data_f$Mixture_Code)
str(data_f)
M2c<-lm(log(JE_Leaf_Count) ~  Mixture_Code1 +Initial_biomass,data=data_f)

lsm<-emmeans(M2c,~Mixture_Code1+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Mixture_Code1+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(8,9,10,11,12,13)
points(position,plot$response,pch=4,cex=2)

#add mean for GO incompetition with JE
Mixture_Code1<-as.factor(data_g$Mixture_Code)
str(data_g)
M3c<-lm(log(GO_Leaf_Count) ~  Mixture_Code1 +Initial_biomass,data=data_g)

lsm<-emmeans(M3c,~Mixture_Code1+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm
pan.emm <- emmeans(regrid(lsm), list(~Mixture_Code1+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(1,2,3,4,5,6)
points(position,plot$response,pch=4,cex=2)


#LEAF LENGTH

boxplot(JE_LLL ~Mixture_Code ,data=data_f)
#Use lmer with RMEL=F
Mixture_Code1<-as.factor(data_f$Mixture_Code)
str(data_f)

# Test homogeneity of variance
leveneTest(JE_LLL~Mixture_Code1,data=data_f)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.52


#Check normality with shapiro test
shapiro.test(data_f$JE_LLL) 

#P VALUE > 0.05 the assumption is valid, p =0.83

#transformation
M2a<-lm(JE_LLL ~  Mixture_Code1 +Initial_biomass,data=data_f)
M2b<-lm(sqrt(JE_LLL) ~  Mixture_Code1 +Initial_biomass,data=data_f)
M2c<-lm(log(JE_LLL) ~  Mixture_Code1 +Initial_biomass,data=data_f)


plot(M2a)  
hist(resid(M2a)) # fit better

plot(M2b)  
hist(resid(M2b)) #

plot(M2c) # 
hist(resid(M2c)) # 

# untransformed fit better
Anova(M2a,type="III",test="F")
##
DunnettTest(JE_LLL ~  Mixture_Code1,data=data_f)
#multiple comparaison

gltt_M2b<-glht(M2a,linfct=mcp(Mixture_Code1="Tukey"))
gltt_M2b
summary(gltt_M2b)
post<-cld(gltt_M2b, level = 0.05, decreasing = T)
post

#GO incompetition JE
data_g<-read.table(file="Harvest1_GO_Fert.txt", header=TRUE,sep="\t",dec = ".")
Mixture_Code1<-as.factor(data_g$Mixture_Code)
str(data_g)

# Test homogeneity of variance
leveneTest(GO_LLL~Mixture_Code1,data=data_g)
#Assumption is valid if pvalue >.05
# assumption is not valid, p =0.02


#Check normality with shapiro test
shapiro.test(data_g$GO_LLL) 

#P VALUE > 0.05 the assumption is valid, p =0.07

#transformation
M3a<-lm(GO_LLL ~  Mixture_Code1 +Initial_biomass,data=data_g)
M3b<-lm(sqrt(GO_LLL) ~  Mixture_Code1 +Initial_biomass,data=data_g)
M3c<-lm(log(GO_LLL) ~  Mixture_Code1 +Initial_biomass,data=data_g)


plot(M3a)  
hist(resid(M3a)) #

plot(M3b)  
hist(resid(M3b)) # fit better

plot(M3c) # 
hist(resid(M3c)) # 

# sqrt transformation fit better
Anova(M3b,type="III",test="F")
##
DunnettTest(log(GO_LLL) ~  Mixture_Code1,data=data_g)
#multiple comparaison

gltt_M2a<-glht(M3b,linfct=mcp(Mixture_Code1="Tukey"))
gltt_M2a
summary(gltt_M2a)
post<-cld(gltt_M2a, level = 0.05, decreasing = T)
post

#Data for Box_plot 
data_box<-read.table(file="Harvest1_JE_GO_Fert.txt", header=TRUE,sep="\t",dec = ".")
str(data_box)

boxplot(LLL ~Mixture_Code,data=data_box)


position<-c(1,2,3,4,5,6,8,9,10,11,12,13)
boxplot(LLL ~Mixture_Code,data=data_box,axes=F,at=position ,border="gray40",xlab= "Population", ylab="Length of the longest leaf (cm)",las=2, 
        col=c("#FFFF00","#FFFF00","#FFFF00","#FFFF00","#FFFF00","#FFFF00","#0033FF","#0033FF","#0033FF","#0033FF","#0033FF","#0033FF"))
axis(2,at=seq(15,40,5), seq(15,40,5), las=1)
axis(1, at=c(1,2,3,4,5,6,8,9,10,11,12,13), labels=c("0:1","1:1","0:2","2:1","1:2","0:3","0:1","1:1","0:2","2:1","1:2","0:3"), las=1)
box()



#add mean for JE incompetition with GO
Mixture_Code1<-as.factor(data_f$Mixture_Code)
str(data_f)
M2a<-lm(JE_LLL ~  Mixture_Code1 +Initial_biomass,data=data_f)

lsm<-emmeans(M2a,~Mixture_Code1+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Mixture_Code1+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(8,9,10,11,12,13)
points(position,plot$emmean,pch=4,cex=2)

#add mean for GO incompetition with JE
Mixture_Code1<-as.factor(data_g$Mixture_Code)
str(data_g)
M3b<-lm(sqrt(GO_LLL)~  Mixture_Code1 +Initial_biomass,data=data_g)

lsm<-emmeans(M3b,~Mixture_Code1+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm
pan.emm <- emmeans(regrid(lsm), list(~Mixture_Code1+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(1,2,3,4,5,6)
points(position,plot$response,pch=4,cex=2)

#SPECIFIC LEAF AREA
data_SP<-read.table(file="JE_LMA.txt", header=TRUE,sep="\t",dec = ".")
boxplot(JE_SLA ~Mixture_Code ,data=data_SP)
#Use lmer with RMEL=F
Mixture_Code1<-as.factor(data_f$Mixture_Code)
str(data_f)

# Test homogeneity of variance
leveneTest(JE_SLA~Mixture_Code1,data=data_f)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.11


#Check normality with shapiro test
shapiro.test(data_f$JE_SLA) 

#P VALUE > 0.05 the assumption is not valid, p =0.0002

#transformation
M2a<-lm(JE_SLA ~  Mixture_Code1 +Initial_biomass,data=data_SP)
M2b<-lm(sqrt(JE_SLA) ~  Mixture_Code1 +Initial_biomass,data=data_SP)
M2c<-lm(log(JE_SLA) ~  Mixture_Code1 +Initial_biomass,data=data_SP)


plot(M2a)  
hist(resid(M2a)) #fit better

plot(M2b)  
hist(resid(M2b)) #

plot(M2c) # 
hist(resid(M2c)) # 

#UN transformed fit better
Anova(M2a,type="III",test="F")
##
DunnettTest(JE_SLA ~  Mixture_Code1,data=data_SP)
#multiple comparaison

gltt_M2b<-glht(M2a,linfct=mcp(Mixture_Code1="Tukey"))
gltt_M2b
summary(gltt_M2b)
post<-cld(gltt_M2b, level = 0.05, decreasing = T)
post

#GO incompetition JE
data_SPB<-read.table(file="GO_LMA.txt", header=TRUE,sep="\t",dec = ".")
Mixture_Code1<-as.factor(data_SPB$Mixture_Code)
str(data_g)

# Test homogeneity of variance
leveneTest(GO_SLA~Mixture_Code1,data=data_SPB)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.46


#Check normality with shapiro test
shapiro.test(data_SPB$GO_SLA) 

#P VALUE > 0.05 the assumption is not valid, p =0.02

#transformation
M3a<-lm(GO_SLA ~  Mixture_Code1 +Initial_biomass,data=data_SPB)
M3b<-lm(sqrt(GO_SLA) ~  Mixture_Code1 +Initial_biomass,data=data_SPB)
M3c<-lm(log(GO_SLA) ~  Mixture_Code1 +Initial_biomass,data=data_SPB)


plot(M3a)  
hist(resid(M3a)) #

plot(M3b)  
hist(resid(M3b)) # FIT BETTER

plot(M3c) # 
hist(resid(M3c)) # 

# transformation fit better
Anova(M3b,type="III",test="F")
##
DunnettTest(log(GO_SLA) ~  Mixture_Code1,data=data_g)
#multiple comparaison

gltt_M2a<-glht(M3b,linfct=mcp(Mixture_Code1="Tukey"))
gltt_M2a
summary(gltt_M2a)
post<-cld(gltt_M2a, level = 0.05, decreasing = T)
post

#Data for Box_plot 
data_box<-read.table(file="Box_LMA.txt", header=TRUE,sep="\t",dec = ".")
str(data_box)

boxplot(SLA  ~Mixture_Code,data=data_box)


position<-c(1,2,3,5,6,7)
boxplot(SLA  ~Mixture_Code,data=data_box,axes=F,at=position ,border="gray40",xlab= "Population", ylab="Specific leaf area (mm2 mg-1)",las=2, 
        col=c("#FFFF00","#FFFF00","#FFFF00","#0033FF","#0033FF","#0033FF"))
axis(2,at=seq(6,14,2), seq(6,14,2), las=1)
axis(1, at=c(1,2,3,5,6,7), labels=c("0:1","1:1","0:2","0:1","1:1","0:2"), las=1)
box()


#add mean for JE incompetition with GO
Mixture_Code1<-as.factor(data_f$Mixture_Code)
str(data_f)
M2c<-lm(JE_SLA ~  Mixture_Code1 +Initial_biomass,data=data_f)

lsm<-emmeans(M2a,~Mixture_Code1+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Mixture_Code1+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(5,6,7)
points(position,plot$response,pch=4,cex=2)

#add mean for GO incompetition with JE
Mixture_Code1<-as.factor(data_g$Mixture_Code)
str(data_g)
M3b<-lm(sqrt(GO_SLA) ~  Mixture_Code1 +Initial_biomass,data=data_g)

lsm<-emmeans(M3b,~Mixture_Code1+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm
pan.emm <- emmeans(regrid(lsm), list(~Mixture_Code1+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(1,2,3)
points(position,plot$response,pch=4,cex=2)

#CHECK CORRELATION OF 2:1 LLL VS 2:1 SLA
### SPEARMAN's CORRELATION
library(ggplot2)
library(lmerTest)
#SLA-LLL-2:1

data_Cor_2.1<-read.table(file="SLA_LLL_2.1.txt", header=TRUE,sep="\t",dec = ".")

cor.test(data_Cor_2.1$JE_SLA,data_Cor_2.1$JE_LLL,method = "spearman")
ggscatter(data_Cor_2.1, x = "JE_SLA", y = "JE_LLL", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Specific leaf area (mm2 mg-1)", ylab = "Length of the longest leaf (cm)")


data_Cor$predlm<-predict(m5a)
ggplot(data_Cor,aes(x=JE_SLA,y=JE_LLL, color=Mixture_Code))+
  geom_point()+ geom_smooth(method = "lm") +theme_bw() + 
  ylab("\nLength of the longest leaf (cm)")+
  xlab("\nSpecific leaf area (mm2 mg-1)")+
  theme(axis.text.x =element_text(vjust = 1, hjust = 1),
        panel.grid = element_blank())+stat_cor(aes(color=Mixture_Code),method = "pearson" )

#LEAF WATER CONTENT

data_f<-read.table(file="Harvest1_JE_Fert.txt", header=TRUE,sep="\t",dec = ".")
str(data_f)

boxplot(JE_LWC ~Mixture_Code ,data=data_f)
#Use lmer with RMEL=F
Mixture_Code1<-as.factor(data_f$Mixture_Code)
str(data_f)

# Test homogeneity of variance
leveneTest(JE_LWC ~Mixture_Code1,data=data_f)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.09


#Check normality with shapiro test
shapiro.test(data_f$JE_LWC) 

#P VALUE > 0.05 the assumption is not valid, p =0.001

#transformation
M2a<-lm(JE_LWC ~  Mixture_Code1 +Initial_biomass,data=data_f)
M2b<-lm(sqrt(JE_LWC) ~  Mixture_Code1 +Initial_biomass,data=data_f)
M2c<-lm(log(JE_LWC) ~  Mixture_Code1 +Initial_biomass,data=data_f)


plot(M2a)  
hist(resid(M2a)) 

plot(M2b)  
hist(resid(M2b)) #fit better

plot(M2c) # 
hist(resid(M2c)) # 

#SQRT transformation fit better
Anova(M2b,type="III",test="F")
##
DunnettTest(sqrt(JE_LWC) ~  Mixture_Code1,data=data_f)
#multiple comparaison

gltt_M2b<-glht(M2b,linfct=mcp(Mixture_Code1="Tukey"))
gltt_M2b
summary(gltt_M2b)
post<-cld(gltt_M2b, level = 0.05, decreasing = T)
post

#GO incompetition JE
data_g<-read.table(file="Harvest1_GO_Fert.txt", header=TRUE,sep="\t",dec = ".")
Mixture_Code1<-as.factor(data_g$Mixture_Code)
str(data_g)

# Test homogeneity of variance
leveneTest(GO_LWC~Mixture_Code1,data=data_g)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.23


#Check normality with shapiro test
shapiro.test(data_g$GO_LWC) 

#P VALUE > 0.05 the assumption is valid, p =0.67

#transformation
M2a<-lm(GO_LWC ~  Mixture_Code1 +Initial_biomass,data=data_g)
M2b<-lm(sqrt(GO_LWC) ~  Mixture_Code1 +Initial_biomass,data=data_g)
M2c<-lm(log(GO_LWC) ~  Mixture_Code1 +Initial_biomass,data=data_g)


plot(M2a)  
hist(resid(M2a)) # FIT BETTER

plot(M2b)  
hist(resid(M2b)) #

plot(M2c) # 
hist(resid(M2c)) # 

# transformation fit better
Anova(M2a,type="III",test="F")
##
DunnettTest(GO_LWC ~  Mixture_Code1,data=data_g)
#multiple comparaison

gltt_M2a<-glht(M2a,linfct=mcp(Mixture_Code1="Tukey"))
gltt_M2a
summary(gltt_M2a)
post<-cld(gltt_M2a, level = 0.05, decreasing = T)
post

#Data for Box_plot 
data_box<-read.table(file="Harvest1_JE_GO_Fert.txt", header=TRUE,sep="\t",dec = ".")
str(data_box)

boxplot(LWC ~Mixture_Code,data=data_box)


position<-c(1,2,3,4,5,6,8,9,10,11,12,13)
boxplot(LWC ~Mixture_Code,data=data_box,axes=F,at=position ,border="gray40",xlab= "Population", ylab="Water content (%)",las=2, 
        col=c("#FFFF00","#FFFF00","#FFFF00","#FFFF00","#FFFF00","#FFFF00","#0033FF","#0033FF","#0033FF","#0033FF","#0033FF","#0033FF"))
axis(2,at=seq(82,92,2), seq(82,92,2), las=1)
axis(1, at=c(1,2,3,4,5,6,8,9,10,11,12,13), labels=c("0:1","1:1","0:2","2:1","1:2","0:3","0:1","1:1","0:2","2:1","1:2","0:3"), las=1)
box()


#add mean for JE incompetition with GO
Mixture_Code1<-as.factor(data_f$Mixture_Code)
str(data_f)
M2b<-lm(sqrt(JE_LWC) ~  Mixture_Code1 +Initial_biomass,data=data_f)

lsm<-emmeans(M2b,~Mixture_Code1+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Mixture_Code1+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(8,9,10,11,12,13)
points(position,plot$response,pch=4,cex=2)

#add mean for GO incompetition with JE
Mixture_Code1<-as.factor(data_g$Mixture_Code)
str(data_g)
M2a<-lm(GO_LWC ~  Mixture_Code1 +Initial_biomass,data=data_g)

lsm<-emmeans(M2a,~Mixture_Code1+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm
pan.emm <- emmeans(regrid(lsm), list(~Mixture_Code1+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(1,2,3,4,5,6)
points(position,plot$emmean,pch=4,cex=2)




#RELATIVE COMPETITION INTENSITY BASED ON CONTROL
#Biomass was considered as measure of yield

#RCI_CONT_TWO
dat_1a<-read.table(file="RCI_Cont_Two.txt", header=TRUE,sep="\t",dec = ".")

str(dat_1a)
#Use lmer with RMEL=F
MC <- lmer(Competitive_ability_JE  ~ Population+(1 | Mixture_Code)+Initial_biomass, REML=F, data = dat_1a)
summary(MC)


#Check themodel(DHARMarequired)
sim1a<-simulateResiduals(fittedModel= MC,n=1000) #fit better
plot(sim1a, quantreg= FALSE)


summary(MC)
#untransformed is better
#Significant test
Anova(MC,type="III",test="Chisq")

#multiple comparaison
library(multcomp)
gltt_M1<-glht(MC,linfct=mcp(Population="Tukey"))
gltt_M1
summary(gltt_M1)
post<-cld(gltt_M1, level = 0.05, decreasing = F)
post
#

#RCI_CONT_THREE
dat_1b<-read.table(file="RCI_Cont_Three.txt", header=TRUE,sep="\t",dec = ".")

str(dat_1b)
#Use lmer with RMEL=F
MB <- lmer(Competitive_ability_JE  ~ Population+(1 | Mixture_Code)+Initial_biomass, REML=F, data = dat_1b)
summary(MC)


#Check themodel(DHARMarequired)
sim1a<-simulateResiduals(fittedModel= MB,n=1000) #fit better
plot(sim1a, quantreg= FALSE)


summary(MB)
#untransformed is better
#Significant test
Anova(MB,type="III",test="Chisq")

#multiple comparaison
library(multcomp)
gltt_M1<-glht(MB,linfct=mcp(Population="Tukey"))
gltt_M1
summary(gltt_M1)
post<-cld(gltt_M1, level = 0.05, decreasing = T)
post
#


#RCI BOX PLOT TWO AND THREE (JE VS GO) IN SUPPLEMENT
#Data for Box_plot 
dat_4b<-read.table(file="RCI_Box.txt", header=TRUE,sep="\t",dec = ".")
str(dat_4b)

boxplot(Competitive_ability_JE ~Population,data=dat_4b)


position<-c(1,2,4,5)
boxplot(Competitive_ability_JE ~Population,data=dat_4b,axes=F,at=position ,border="gray40",xlab= "Population", ylab="Relative competitive intensity",las=2, 
        col=c("#FFFF00","#0033FF","#FFFF00","#0033FF"))
axis(2,at=seq(0.0,0.8,0.2), seq(-0.0,0.8,0.2), las=1)
axis(1, at=c(1,2,4,5), labels=c("2_GO","2_JE","3_GO","3_JE"), las=1)
box()



#add mean for two competition

MC <- lmer(Competitive_ability_JE  ~ Population+(1 | Mixture_Code)+Initial_biomass, REML=F, data = dat_1a)

lsm<-emmeans(MC,~Population+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Population+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(1,2)
points(position,plot$emmean,pch=4,cex=2)

#add mean for three competition

str(dat_2b)
MB <- lmer(Competitive_ability_JE  ~ Population+(1 | Mixture_Code)+Initial_biomass, REML=F, data = dat_1b)

lsm<-emmeans(MB,~Population+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Population+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(4,5)
points(position,plot$emmean,pch=4,cex=2)
##


#RCI_CONT (Compare GO vs JE based on similar control)
dat_1c<-read.table(file="RCI_Cont.txt", header=TRUE,sep="\t",dec = ".")

str(dat_1c)
#Use lmer with RMEL=F
MD <- lmer(Competitive_ability_JE  ~ Population+(1 | Mixture_Code)+Initial_biomass, REML=F, data = dat_1c)
summary(MD)


#Check themodel(DHARMarequired)
sim1a<-simulateResiduals(fittedModel= MD,n=1000) #fit better
plot(sim1a, quantreg= FALSE)

#untransformed is better
#Significant test
Anova(MD,type="III",test="Chisq")


#COMPETITION CONTROL JE

# use lm   

dat_2a<-read.table(file="Competition_Cont_JE.txt", header=TRUE,sep="\t",dec = ".")
Mixture_Code1<-as.factor(dat_2a$Mixture_Code)
str(dat_2a)
# lm


# Test homogeneity of variance
leveneTest(Competitive_ability_JE~Mixture_Code,data=dat_2a)
#Assumption is valid if pvalue >.05
# assumption is not valid, p =0.0004


#Check normality with shapiro test
shapiro.test(dat_2a$Competitive_ability_JE) 

#P VALUE > 0.05 the assumption is valid, p =0.976
M2C<-lm(Competitive_ability_JE ~  Mixture_Code1 +Initial_biomass,data=dat_2a)
M2B<-lm(log(Competitive_ability_JE) ~  Mixture_Code1 +Initial_biomass,data=dat_2a)
M2D<-lm(sqrt(Competitive_ability_JE) ~  Mixture_Code1 +Initial_biomass,data=dat_2a)
#LM
plot(M2C)  
hist(resid(M2C)) 

plot(M2B)  
hist(resid(M2B))

plot(M2D)  
hist(resid(M2D)) # fit better


#untransformed fit better
Anova(M2D,type="III",test="F")
##
#Multiple comparison
gltt_M2a<-glht(M2D,linfct=mcp(Mixture_Code1="Tukey"))
gltt_M2a
summary(gltt_M2a)
post<-cld(gltt_M2a, level = 0.05, decreasing = F)
post

#COMPETITION CONTROL GO

# use lm   

dat_2b<-read.table(file="Competition_Cont_GO.txt", header=TRUE,sep="\t",dec = ".")
Mixture_Code1<-as.factor(dat_2b$Mixture_Code)
str(dat_2b)
# lm


# Test homogeneity of variance
leveneTest(Competitive_ability_JE~Mixture_Code,data=dat_2b)
#Assumption is valid if pvalue >.05
# assumption is not valid, p =0.0001


#Check normality with shapiro test
shapiro.test(dat_2b$Competitive_ability_JE) 

#P VALUE > 0.05 the assumption is valid, p =0.867
M2C<-lm(Competitive_ability_JE ~  Mixture_Code1 +Initial_biomass,data=dat_2b)
M2B<-lm(log(Competitive_ability_JE) ~  Mixture_Code1 +Initial_biomass,data=dat_2b)
M2D<-lm(sqrt(Competitive_ability_JE) ~  Mixture_Code1 +Initial_biomass,data=dat_2b)
#LM
plot(M2C)  
hist(resid(M2C)) #fit better

plot(M2B)  
hist(resid(M2B))

plot(M2D)  
hist(resid(M2D)) # 


#untransformed fit better
Anova(M2C,type="III",test="F")
##
#Multiple comparison
gltt_M2a<-glht(M2C,linfct=mcp(Mixture_Code1="Tukey"))
gltt_M2a
summary(gltt_M2a)
post<-cld(gltt_M2a, level = 0.05, decreasing = F)
post

#COMPETITION_CONTROL_1:1

# use lm   

dat_2c<-read.table(file="Competition_Cont_1.1.txt", header=TRUE,sep="\t",dec = ".")
Mixture_Code1<-as.factor(dat_2c$Mixture_Code)
str(dat_2c)
# lm


# Test homogeneity of variance
leveneTest(Competitive_ability_JE~Mixture_Code,data=dat_2c)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.8


#Check normality with shapiro test
shapiro.test(dat_2c$Competitive_ability_JE) 

#P VALUE > 0.05 the assumption is valid, p =0.80
M2C<-lm(Competitive_ability_JE ~  Mixture_Code1 +Initial_biomass,data=dat_2c)
M2B<-lm(log(Competitive_ability_JE) ~  Mixture_Code1 +Initial_biomass,data=dat_2c)
M2D<-lm(sqrt(Competitive_ability_JE) ~  Mixture_Code1 +Initial_biomass,data=dat_2c)
#LM
plot(M2C)  
hist(resid(M2C)) #fit better

plot(M2B)  
hist(resid(M2B))

plot(M2D)  
hist(resid(M2D)) # 


#untransformed fit better
Anova(M2C,type="III",test="F")
##
#Multiple comparison
gltt_M2a<-glht(M2C,linfct=mcp(Mixture_Code1="Tukey"))
gltt_M2a
summary(gltt_M2a)
post<-cld(gltt_M2a, level = 0.05, decreasing = T)
post

# COMPETITION CONTROL 0:2
dat_2d<-read.table(file="Competition_Cont_0.2.txt", header=TRUE,sep="\t",dec = ".")
Mixture_Code1<-as.factor(dat_2d$Mixture_Code)
str(dat_2d)
# lm


# Test homogeneity of variance
leveneTest(Competitive_ability_JE~Mixture_Code,data=dat_2d)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.95


#Check normality with shapiro test
shapiro.test(dat_2d$Competitive_ability_JE) 

#P VALUE > 0.05 the assumption is valid, p =0.93
M2C<-lm(Competitive_ability_JE ~  Mixture_Code1 +Initial_biomass,data=dat_2d)
M2B<-lm(log(Competitive_ability_JE) ~  Mixture_Code1 +Initial_biomass,data=dat_2d)
M2D<-lm(sqrt(Competitive_ability_JE) ~  Mixture_Code1 +Initial_biomass,data=dat_2d)
#LM
plot(M2C)  
hist(resid(M2C)) #fit better

plot(M2B)  
hist(resid(M2B))

plot(M2D)  
hist(resid(M2D)) # 


#untransformed fit better
Anova(M2C,type="III",test="F")
##
#Multiple comparison
gltt_M2a<-glht(M2C,linfct=mcp(Mixture_Code1="Tukey"))
gltt_M2a
summary(gltt_M2a)
post<-cld(gltt_M2a, level = 0.05, decreasing = T)
post

#COMPETITION CONTROL 2:1
dat_2e<-read.table(file="Competition_Cont_2.1.txt", header=TRUE,sep="\t",dec = ".")
Mixture_Code1<-as.factor(dat_2e$Mixture_Code)
str(dat_2e)
# lm


# Test homogeneity of variance
leveneTest(Competitive_ability_JE~Mixture_Code,data=dat_2e)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.11


#Check normality with shapiro test
shapiro.test(dat_2e$Competitive_ability_JE) 

#P VALUE > 0.05 the assumption is valid, p =0.58
M2C<-lm(Competitive_ability_JE ~  Mixture_Code1 +Initial_biomass,data=dat_2e)
M2B<-lm(log(Competitive_ability_JE) ~  Mixture_Code1 +Initial_biomass,data=dat_2e)
M2D<-lm(sqrt(Competitive_ability_JE) ~  Mixture_Code1 +Initial_biomass,data=dat_2e)
#LM
plot(M2C)  
hist(resid(M2C)) #

plot(M2B)  
hist(resid(M2B))

plot(M2D)  
hist(resid(M2D)) # fit better


#untransformed fit better
Anova(M2D,type="III",test="F")
##
#Multiple comparison
gltt_M2a<-glht(M2D,linfct=mcp(Mixture_Code1="Tukey"))
gltt_M2a
summary(gltt_M2a)
post<-cld(gltt_M2a, level = 0.05, decreasing = T)
post

#COMPETITION CONTROL_1:2

dat_2f<-read.table(file="Competition_Cont_1.2.txt", header=TRUE,sep="\t",dec = ".")
Mixture_Code1<-as.factor(dat_2f$Mixture_Code)
str(dat_2f)
# lm


# Test homogeneity of variance
leveneTest(Competitive_ability_JE~Mixture_Code,data=dat_2f)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.87


#Check normality with shapiro test
shapiro.test(dat_2f$Competitive_ability_JE) 

#P VALUE > 0.05 the assumption is valid, p =0.59
M2C<-lm(Competitive_ability_JE ~  Mixture_Code1 +Initial_biomass,data=dat_2f)
M2B<-lm(log(Competitive_ability_JE) ~  Mixture_Code1 +Initial_biomass,data=dat_2f)
M2D<-lm(sqrt(Competitive_ability_JE) ~  Mixture_Code1 +Initial_biomass,data=dat_2f)
#LM
plot(M2C)  
hist(resid(M2C)) # fit better

plot(M2B)  
hist(resid(M2B))

plot(M2D)  
hist(resid(M2D)) # 

#untransformed fit better
Anova(M2C,type="III",test="F")
##
#Multiple comparison
gltt_M2a<-glht(M2D,linfct=mcp(Mixture_Code1="Tukey"))
gltt_M2a
summary(gltt_M2a)
post<-cld(gltt_M2a, level = 0.05, decreasing = T)
post

#COMPETITION CONTROL_0:3

dat_2g<-read.table(file="Competition_Cont_0.3.txt", header=TRUE,sep="\t",dec = ".")
Mixture_Code1<-as.factor(dat_2g$Mixture_Code)
str(dat_2g)
# lm


# Test homogeneity of variance
leveneTest(Competitive_ability_JE~Mixture_Code,data=dat_2g)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.96


#Check normality with shapiro test
shapiro.test(dat_2g$Competitive_ability_JE) 

#P VALUE > 0.05 the assumption is valid, p =0.054
M2C<-lm(Competitive_ability_JE ~  Mixture_Code1 +Initial_biomass,data=dat_2g)
M2B<-lm(log(Competitive_ability_JE) ~  Mixture_Code1 +Initial_biomass,data=dat_2g)
M2D<-lm(sqrt(Competitive_ability_JE) ~  Mixture_Code1 +Initial_biomass,data=dat_2g)
#LM
plot(M2C)  
hist(resid(M2C)) # 

plot(M2B)  
hist(resid(M2B))

plot(M2D)  
hist(resid(M2D)) # fit better

#untransformed fit better
Anova(M2D,type="III",test="F")
##
#Multiple comparison
gltt_M2a<-glht(M2D,linfct=mcp(Mixture_Code1="Tukey"))
gltt_M2a
summary(gltt_M2a)
post<-cld(gltt_M2a, level = 0.05, decreasing = T)
post

#Data for Box_plot 
dat_4a<-read.table(file="Box_Cont.txt", header=TRUE,sep="\t",dec = ".")
str(dat_4a)

boxplot(Competitive_ability_JE ~Mixture_Code,data=dat_4a)


position<-c(1,2,3,4,5,7,8,9,10,11)
boxplot(Competitive_ability_JE ~Mixture_Code,data=dat_4a,axes=F,at=position ,border="gray40",xlab= "Population", ylab="Relative competitive intensity",las=2, 
        col=c("#FFFF00","#FFFF00","#FFFF00","#FFFF00","#FFFF00","#0033FF","#0033FF","#0033FF","#0033FF","#0033FF"))
axis(2,at=seq(0.0,0.8,0.2), seq(-0.0,0.8,0.2), las=1)
axis(1, at=c(1,2,3,4,5,7,8,9,10,11), labels=c("1:1","0:2","2:1","1:2","0:3","1:1","0:2","2:1","1:2","0:3"), las=1)
box()



#add mean for JE incompetition with GO
Mixture_Code1<-as.factor(dat_2a$Mixture_Code)
str(dat_2a)
M2D<-lm(sqrt(Competitive_ability_JE) ~  Mixture_Code1 +Initial_biomass,data=dat_2a)			

lsm<-emmeans(M2D,~Mixture_Code1+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Mixture_Code1+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(7,8,9,10,11)
points(position,plot$response,pch=4,cex=2)

#add mean for GO incompetition with JE
Mixture_Code1<-as.factor(dat_2b$Mixture_Code)
str(dat_2b)
M2C<-lm(Competitive_ability_JE ~  Mixture_Code1 +Initial_biomass,data=dat_2b)			


lsm<-emmeans(M2C,~Mixture_Code1+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Mixture_Code1+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(1,2,3,4,5)
points(position,plot$emmean,pch=4,cex=2)
##

#RELATIVE COMPETITIVE ABILITY BASED ON MONOCULTURE

#Biomass was considered as measure of yield

#BASED ON MONOCULTURE OF TWO
dat_2c2<-read.table(file="Competition_Mon_JE_GO_1.1.txt", header=TRUE,sep="\t",dec = ".")

str(dat_2c2)

# use lm for JE incompetition GO  

str(dat_2c2)
# lm



# Test homogeneity of variance
leveneTest(Competitive_ability_JE~Population,data=dat_2c2)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.64


#Check normality with shapiro test
shapiro.test(dat_2c2$Competitive_ability_JE) 

#P VALUE > 0.05 the assumption is valid, p =0.76

M2C<-lm(Competitive_ability_JE ~  Population +Initial_biomass,data=dat_2c2)
M2B<-lm(log(Competitive_ability_JE) ~  Population +Initial_biomass,data=dat_2c2)
M2D<-lm(sqrt(Competitive_ability_JE) ~  Population +Initial_biomass,data=dat_2c2)
#LM
plot(M2C)  
hist(resid(M2C)) #fit better
plot(M2B)  
hist(resid(M2B))
plot(M2D)  
hist(resid(M2D))


#untransformed fit better
Anova(M2C,type="III",test="F")
##
#multiple comparaison
Mixture_Code1<-as.factor(dat_2c2$Mixture_Code)
str(dat_2c2)
M2C<-lm(Competitive_ability_JE ~  Mixture_Code1 +Initial_biomass,data=dat_2c2)

gltt_M2a<-glht(M2C,linfct=mcp(Mixture_Code1="Tukey"))
gltt_M2a
summary(gltt_M2a)
post<-cld(gltt_M2a, level = 0.05, decreasing = F)
post
#BASED ON MONOCULTURE OF THREE _JE
dat_3c2<-read.table(file="JE_Mon_Competition.txt", header=TRUE,sep="\t",dec = ".")
Mixture_Code1<-as.factor(dat_3c2$Mixture_series)
str(dat_3c2)
M3C<-lm(Competitive_ability_JE  ~  Mixture_Code1 +Initial_biomass,data=dat_3c2)
# Test homogeneity of variance
leveneTest(Competitive_ability_JE ~Mixture_Code,data=dat_3c2)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.13

#Check normality with shapiro test
shapiro.test(dat_3c2$Competitive_ability_JE ) 

#P VALUE > 0.05 the assumption is valid, p =0.42


plot(M3C) # 
hist(resid(M3C)) #


summary(M3)

#transformed data fit better
Anova(M3C,type="III",test="F")
##

#multiple comparaison
Mixture_Code1<-as.factor(dat_3c2$Mixture_Code)
str(dat_3c2)
M3C<-lm(Competitive_ability_JE  ~  Mixture_Code1 +Initial_biomass,data=dat_3c2)

gltt_M3C<-glht(M3C,linfct=mcp(Mixture_Code1="Tukey"))
gltt_M3C
summary(gltt_M3C)
post<-cld(gltt_M3C, level = 0.05, decreasing = F)
post

##BASED ON MONOCULTURE OF THREE _GO
dat_3c3<-read.table(file="GO_Mon_Competition.txt", header=TRUE,sep="\t",dec = ".")
Mixture_Code1<-as.factor(dat_3c3$Mixture_series)
str(dat_3c3)

# Test homogeneity of variance
leveneTest(Competitive_ability_JE ~Mixture_Code,data=dat_3c3)
#Assumption is valid if pvalue >.05
# assumption is NOT valid, p =0.0006

#Check normality with shapiro test
shapiro.test(dat_3c3$Competitive_ability_JE ) 

#P VALUE > 0.05 the assumption is valid, p =0.11

M3C<-lm(Competitive_ability_JE  ~  Mixture_Code1 +Initial_biomass,data=dat_3c3)
plot(M3C) # 
hist(resid(M3C)) #


summary(M3)

#transformed data fit better
Anova(M3C,type="III",test="F")
##

#multiple comparaison
Mixture_Code1<-as.factor(dat_3c3$Mixture_Code)
str(dat_3c3)
M3C<-lm(Competitive_ability_JE  ~  Mixture_Code1 +Initial_biomass,data=dat_3c3)

gltt_M3C<-glht(M3C,linfct=mcp(Mixture_Code1="Tukey"))
gltt_M3C
summary(gltt_M3C)
post<-cld(gltt_M3C, level = 0.05, decreasing = F)
post

##BASED ON MONOCULTURE OF THREE _JE_GO_2:1
dat_3c4<-read.table(file="Competition_Mon_JE_GO_2.1.txt", header=TRUE,sep="\t",dec = ".")
Mixture_Code1<-as.factor(dat_3c4$Mixture_series)
str(dat_3c4)

# Test homogeneity of variance
leveneTest(Competitive_ability_JE ~Mixture_Code,data=dat_3c4)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.08

#Check normality with shapiro test
shapiro.test(dat_3c4$Competitive_ability_JE ) 

#P VALUE > 0.05 the assumption is valid, p =0.58

M3C<-lm(Competitive_ability_JE  ~  Population +Initial_biomass,data=dat_3c4)
plot(M3C) # 
hist(resid(M3C)) #


summary(M3)

#transformed data fit better
Anova(M3C,type="III",test="F")
##

#multiple comparaison
Mixture_Code1<-as.factor(dat_3c4$Population)
str(dat_3c4)
M3C<-lm(Competitive_ability_JE  ~  Mixture_Code1 +Initial_biomass,data=dat_3c4)

gltt_M3C<-glht(M3C,linfct=mcp(Mixture_Code1="Tukey"))
gltt_M3C
summary(gltt_M3C)
post<-cld(gltt_M3C, level = 0.05, decreasing = T)
post


##BASED ON MONOCULTURE OF THREE _JE_GO_1:2
dat_3c5<-read.table(file="Competition_Mon_JE_GO_1.2.txt", header=TRUE,sep="\t",dec = ".")
Mixture_Code1<-as.factor(dat_3c5$Mixture_Code)
str(dat_3c5)

# Test homogeneity of variance
leveneTest(Competitive_ability_JE ~Mixture_Code,data=dat_3c5)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.95

#Check normality with shapiro test
shapiro.test(dat_3c5$Competitive_ability_JE ) 

#P VALUE > 0.05 the assumption is valid, p =0.63

M3C<-lm(Competitive_ability_JE  ~  Population +Initial_biomass,data=dat_3c5)
plot(M3C) # 
hist(resid(M3C)) #


summary(M3)

#transformed data fit better
Anova(M3C,type="III",test="F")
##

#multiple comparaison
Mixture_Code1<-as.factor(dat_3c5$Population)
str(dat_3c4)
M3C<-lm(Competitive_ability_JE  ~  Mixture_Code1 +Initial_biomass,data=dat_3c5)

gltt_M3C<-glht(M3C,linfct=mcp(Mixture_Code1="Tukey"))
gltt_M3C
summary(gltt_M3C)
post<-cld(gltt_M3C, level = 0.05, decreasing = T)
post



#Data for Box_plot 
dat_4ac<-read.table(file="Box_Mon.txt", header=TRUE,sep="\t",dec = ".")
str(dat_4ac)

boxplot(Competitive_ability_JE ~Mixture_Code,data=dat_4ac)


position<-c(1,2,3.5,4.5,6.5,7.5)
boxplot(Competitive_ability_JE ~Mixture_Code,data=dat_4ac,axes=F,at=position ,border="gray40",xlab= "Population", ylab="Relative competitive intensity",las=2, 
        col=c("#FFFF00","#FFFF00","#FFFF00","#0033FF","#0033FF","#0033FF"))
axis(2,at=seq(-1.5,0.5,0.5), seq(-1.5,0.5,0.5), las=1)
axis(1, at=c(1,2,3.5,4.5,6.5,7.5), labels=c("1:1","2:1","1:2","1:1","2:1","1:2"), las=1)
box()


#add mean for Monoculture in two competition
Mixture_Code1<-as.factor(dat_2c2$Mixture_Code)
str(dat_2c2)
M2C<-lm(Competitive_ability_JE ~  Mixture_Code1 +Initial_biomass,data=dat_2c2)						


lsm<-emmeans(M2C,~Mixture_Code1+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Mixture_Code1+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(1,4.5)
points(position,plot$emmean,pch=4,cex=2)

#add mean for Monoculture for GO_Mon_Competition
Mixture_Code1<-as.factor(dat_3c2$Mixture_Code)
str(dat_3c2)
M3C<-lm(Competitive_ability_JE  ~  Mixture_Code1 +Initial_biomass,data=dat_3c2)						


lsm<-emmeans(M3C,~Mixture_Code1+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Mixture_Code1+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(2,3.5)
points(position,plot$emmean,pch=4,cex=2)
##

#add mean for Monoculture for JE_Mon_Competition
Mixture_Code1<-as.factor(dat_3c3$Mixture_Code)
str(dat_3c3)
M3C<-lm(Competitive_ability_JE  ~  Mixture_Code1 +Initial_biomass,data=dat_3c3)						


lsm<-emmeans(M3C,~Mixture_Code1+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Mixture_Code1+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(6.5,7.5)
points(position,plot$emmean,pch=4,cex=2)
##

###COMPETITION RCI MONOCULTURE THREE (JE VS GO)
dat_4ad<-read.table(file="RCI_Mon_Three.txt", header=TRUE,sep="\t",dec = ".")
str(dat_4ad)
#Use lmer with RMEL=F
MC <- lmer(Competitive_ability_JE  ~ Population+(1 | Mixture_Code)+Initial_biomass, REML=F, data = dat_4ad)

#Check themodel(DHARMarequired)
sim1a<-simulateResiduals(fittedModel= MC,n=1000) #fit better
plot(sim1a, quantreg= FALSE)


summary(MC)
#untransformed is better
#Significant test
Anova(MC,type="III",test="Chisq")

#multiple comparaison
library(multcomp)
gltt_M1<-glht(MC,linfct=mcp(Population="Tukey"))
gltt_M1
summary(gltt_M1)
post<-cld(gltt_M1, level = 0.05, decreasing = T)
post
#

#Data for Box_plot 
dat_4ae<-read.table(file="Box_Mon_Two_Three.txt", header=TRUE,sep="\t",dec = ".")
str(dat_4ad)

boxplot(Competitive_ability_JE ~Population,data=dat_4ae)


position<-c(1,2,4,5)
boxplot(Competitive_ability_JE ~Population,data=dat_4ae,axes=F,at=position ,border="gray40",xlab= "Population", ylab="Relative competitive intensity",las=2, 
        col=c("#FFFF00","#0033FF","#FFFF00","#0033FF"))
axis(2,at=seq(-1.5,0.5,0.5), seq(-1.5,0.5,0.5), las=1)
axis(1, at=c(1,2,4,5), labels=c("GO_2","JE_2","GO_3","JE_3"), las=1)
box()

#add mean for Monoculture in two competition JE VS GO
Mixture_Code1<-as.factor(dat_2c2$Population)
str(dat_2c2)
M2C<-lm(Competitive_ability_JE ~  Mixture_Code1 +Initial_biomass,data=dat_2c2)			


lsm<-emmeans(M2C,~Mixture_Code1+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Mixture_Code1+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(1,2)
points(position,plot$emmean,pch=4,cex=2)

#add mean for Monoculture in three competition JE VS GO
Mixture_Code1<-as.factor(dat_4ad$Population)
str(dat_4ad)
MC <- lmer(Competitive_ability_JE  ~ Mixture_Code1+(1 | Mixture_Code)+Initial_biomass, REML=F, data = dat_4ad)


lsm<-emmeans(MC,~Mixture_Code1+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Mixture_Code1+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(4,5)
points(position,plot$emmean,pch=4,cex=2)



#############################################################
#SECOND HARVEST


#JENA
data_f<-read.table(file="Harvest2_JE_Fert.txt", header=TRUE,sep="\t",dec = ".")
str(data_f)

#GO
data_b<-read.table(file="Harvest2_GO_Fert.txt", header=TRUE,sep="\t",dec = ".")
str(data_b)

#POPULATION
#single (JE VS GO)BIOMASS
data_a<-read.table(file="Single2.txt", header=TRUE,sep="\t",dec = ".")
str(data_a)

# use lm for JE incompetition GO  
str(data_a)

# Test homogeneity of variance
leveneTest(Biomass~Population,data=data_a)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.075


#Check normality with shapiro test
shapiro.test(data_a$Biomass) 

#P VALUE > 0.05 the assumption is valid, p =0.39

#transformation
Ms<-lm(Biomass ~  Population +Initial_biomass,data=data_a)

Anova(Ms,type="III",test="F")

#Single (JE VS GO)GO_Leaf_Count

# Test homogeneity of variance
leveneTest(Leaf_Count~Population,data=data_a)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.22


#Check normality with shapiro test
shapiro.test(data_a$Leaf_Count) 

#P VALUE > 0.05 the assumption is not valid, p =0.0001

#transformation
Msa<-lm(Leaf_Count ~  Population +Initial_biomass,data=data_a)
Msb<-lm(sqrt(Leaf_Count) ~  Population +Initial_biomass,data=data_a)
Msc<-lm(log(Leaf_Count) ~  Population +Initial_biomass,data=data_a)


plot(Msa) # 
hist(resid(Msa))

plot(Msb) # 
hist(resid(Msb)) #

plot(Msc) # 
hist(resid(Msc)) #fit better
summary(m4)

#sqrt transformation fit better
Anova(Msc,type="III",test="F")

## Single (JE VS GO) GO_LLL

# Test homogeneity of variance
leveneTest(LLL~Population,data=data_a)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.16


#Check normality with shapiro test
shapiro.test(data_a$LLL) 

#P VALUE > 0.05 the assumption is not valid, p =0.0003

#transformation
Msa<-lm(LLL ~  Population +Initial_biomass,data=data_a)
Msb<-lm(sqrt(LLL) ~  Population +Initial_biomass,data=data_a)
Msc<-lm(log(LLL) ~  Population +Initial_biomass,data=data_a)


plot(Msa) # 
hist(resid(Msa))

plot(Msb) # 
hist(resid(Msb)) #

plot(Msc) # 
hist(resid(Msc)) #fit better
summary(m4)
Anova(Msc,type="III",test="F")

#Single (JE VS GO) LWC

# Test homogeneity of variance
leveneTest(LWC~Population,data=data_a)
#Assumption is valid if pvalue >.05
# assumption is not valid, p =0.77


#Check normality with shapiro test
shapiro.test(data_a$LWC) 

#P VALUE > 0.05 the assumption is valid, p =0.001

Msa<-lm(LWC ~  Population +Initial_biomass,data=data_a)

Msb<-lm(sqrt(LWC) ~  Population +Initial_biomass,data=data_a)
Msc<-lm(log(LWC) ~  Population +Initial_biomass,data=data_a)


plot(Msa) # 
hist(resid(Msa))

plot(Msb) # 
hist(resid(Msb)) #

plot(Msc) # 
hist(resid(Msc)) # fit better
summary(m4)

#transformation 
Anova(Msc,type="III",test="F")

#Single (JE vs GO) SLA

# Test homogeneity of variance
leveneTest(SLA~Population,data=data_a)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.85


#Check normality with shapiro test
shapiro.test(data_a$SLA) 

#P VALUE > 0.05 the assumption is valid, p =0.27

Msa<-lm(SLA ~  Population +Initial_biomass,data=data_a)
Msb<-lm(sqrt(SLA) ~  Population +Initial_biomass,data=data_a)
Msc<-lm(log(SLA) ~  Population +Initial_biomass,data=data_a)


plot(Msa) # 
hist(resid(Msa))#fit better

plot(Msb) # 
hist(resid(Msb)) #

plot(Msc) # 
hist(resid(Msa)) # 
#transformation 
Anova(Msa,type="III",test="F")

## Single(JE VS GO) C/N
# Test homogeneity of variance
leveneTest(C_N.~Population,data=data_a)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.35


#Check normality with shapiro test
shapiro.test(data_a$C_N.) 

#P VALUE > 0.05 the assumption is NOT valid, p =0.002

Msa<-lm(C_N. ~  Population +Initial_biomass,data=data_a)

Msb<-lm(sqrt(C_N.) ~  Population +Initial_biomass,data=data_a)
Msc<-lm(log(C_N.) ~  Population +Initial_biomass,data=data_a)


plot(Msa) # 
hist(resid(Msa))

plot(Msb) # 
hist(resid(Msb)) #

plot(Msc) # 
hist(resid(Msc)) # fit better
summary(m4)

#transformation 
Anova(Msc,type="III",test="F")

## Single(JE VS GO) N
str(data_a)
# Test homogeneity of variance
leveneTest(N_.~Population,data=data_a)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.62


#Check normality with shapiro test
shapiro.test(data_a$N_.) 

#P VALUE > 0.05 the assumption is not valid, p =0.01

Msa<-lm(N_. ~  Population +Initial_biomass,data=data_a)

Msb<-lm(sqrt(N_.) ~  Population +Initial_biomass,data=data_a)
Msc<-lm(log(N_.) ~  Population +Initial_biomass,data=data_a)


plot(Msa) # 
hist(resid(Msa)) # 

plot(Msb) # 
hist(resid(Msb)) #

plot(Msc) # 
hist(resid(Msc)) # fit better
summary(m4)

#transformation 
Anova(Msc,type="III",test="F")

## Single(JE VS GO) C
str(data_a)
# Test homogeneity of variance
leveneTest(C_.~Population,data=data_a)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.68


#Check normality with shapiro test
shapiro.test(data_a$C_.) 

#P VALUE > 0.05 the assumption is valid, p =0.29

Msa<-lm(C_. ~  Population +Initial_biomass,data=data_a)

Msb<-lm(sqrt(C_.) ~  Population +Initial_biomass,data=data_a)
Msc<-lm(log(C_.) ~  Population +Initial_biomass,data=data_a)


plot(Msa) # 
hist(resid(Msa)) #fit better

plot(Msb) # 
hist(resid(Msb)) # 

plot(Msc) # 
hist(resid(Msc)) # 
summary(m4)

#transformation 
Anova(Msa,type="III",test="F")

##INTER1_TWO POPULATIONS
data_b<-read.table(file="Inter2_Two.txt", header=TRUE,sep="\t",dec = ".")
str(data_b)
#BIOMASS
# Test homogeneity of variance
leveneTest(Biomass~Population,data=data_b)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.72


#Check normality with shapiro test
shapiro.test(data_b$Biomass) 

#P VALUE > 0.05 the assumption is valid, p =0.47

Msa<-lm(Biomass ~  Population +Initial_biomass,data=data_b)

Msb<-lm(sqrt(Biomass) ~  Population +Initial_biomass,data=data_b)
Msc<-lm(log(Biomass) ~  Population +Initial_biomass,data=data_b)


plot(Msa) # 
hist(resid(Msa)) # FIT BETTER

plot(Msb) # 
hist(resid(Msb)) #

plot(Msc) # 
hist(resid(Msc)) # 
summary(m4)

#transformation 
Anova(Msa,type="III",test="F")

#TOTAL LEAVES
# Test homogeneity of variance
leveneTest(Leaf_Count~Population,data=data_b)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.16


#Check normality with shapiro test
shapiro.test(data_b$Leaf_Count) 

#P VALUE > 0.05 the assumption is not valid, p =0.001

Msa<-lm(Leaf_Count ~  Population +Initial_biomass,data=data_b)

Msb<-lm(sqrt(Leaf_Count) ~  Population +Initial_biomass,data=data_b)
Msc<-lm(log(Leaf_Count) ~  Population +Initial_biomass,data=data_b)


plot(Msa) # 
hist(resid(Msa))

plot(Msb) # 
hist(resid(Msb)) #

plot(Msc) # 
hist(resid(Msc)) # fit better
summary(m4)

#transformation 
Anova(Msc,type="III",test="F")

#LENGTH LONGEST LEAF

#Test homogeneity of variance
leveneTest(LLL ~Population,data=data_b)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.43


#Check normality with shapiro test
shapiro.test(data_b$LLL) 

#P VALUE > 0.05 the assumption is NOT valid, p =0.03

Msa<-lm(LLL ~  Population +Initial_biomass,data=data_b)

Msb<-lm(sqrt(LLL) ~  Population +Initial_biomass,data=data_b)
Msc<-lm(log(LLL) ~  Population +Initial_biomass,data=data_b)


plot(Msa) # 
hist(resid(Msa)) # 

plot(Msb) # 
hist(resid(Msb)) # fit better

plot(Msc) # 
hist(resid(Msc)) # 
summary(m4)

#transformation 
Anova(Msb,type="III",test="F")

#SPECIFIC LEAF AREA
#Test homogeneity of variance
leveneTest(SLA ~Population,data=data_b)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.07


#Check normality with shapiro test
shapiro.test(data_b$SLA) 

#P VALUE > 0.05 the assumption is valid, p =0.84

Msa<-lm(SLA ~  Population +Initial_biomass,data=data_b)

Msb<-lm(sqrt(SLA) ~  Population +Initial_biomass,data=data_b)
Msc<-lm(log(SLA) ~  Population +Initial_biomass,data=data_b)


plot(Msa) # 
hist(resid(Msa))

plot(Msb) # 
hist(resid(Msb)) # 

plot(Msc) # 
hist(resid(Msc)) #fit better
summary(m4)

#transformation 
Anova(Msc,type="III",test="F")

#LEAF WATER CONTENT
#Test homogeneity of variance
leveneTest(LWC ~Population,data=data_b)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.34


#Check normality with shapiro test
shapiro.test(data_b$LWC) 

#P VALUE > 0.05 the assumption is not valid, p =0.002

Msa<-lm(LWC ~  Population +Initial_biomass,data=data_b)

Msb<-lm(sqrt(LWC) ~  Population +Initial_biomass,data=data_b)
Msc<-lm(log(LWC) ~  Population +Initial_biomass,data=data_b)


plot(Msa) # 
hist(resid(Msa))

plot(Msb) # 
hist(resid(Msb)) #

plot(Msc) # 
hist(resid(Msc)) # fit better
summary(m4)

#transformation 
Anova(Msc,type="III",test="F")

# CARBON TO NITROGEN CONTENT
#Test homogeneity of variance
leveneTest(C_N. ~Population,data=data_b)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.14


#Check normality with shapiro test
shapiro.test(data_b$C_N.) 

#P VALUE > 0.05 the assumption is valid, p =0.35

Msa<-lm(C_N. ~  Population +Initial_biomass,data=data_b)

Msb<-lm(sqrt(C_N.) ~  Population +Initial_biomass,data=data_b)
Msc<-lm(log(C_N.) ~  Population +Initial_biomass,data=data_b)


plot(Msa) # 
hist(resid(Msa)) # fit better

plot(Msb) # 
hist(resid(Msb)) #

plot(Msc) # 
hist(resid(Msc)) # 
summary(m4)

#transformation 
Anova(Msa,type="III",test="F")

## INTER 1(JE VS GO) N
str(data_b)
# Test homogeneity of variance
leveneTest(N_.~Population,data=data_b)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.21


#Check normality with shapiro test
shapiro.test(data_b$N_.) 

#P VALUE > 0.05 the assumption is NOT valid, p =0.01

Msa<-lm(N_. ~  Population +Initial_biomass,data=data_b)

Msb<-lm(sqrt(N_.) ~  Population +Initial_biomass,data=data_b)
Msc<-lm(log(N_.) ~  Population +Initial_biomass,data=data_b)


plot(Msa) # 
hist(resid(Msa))

plot(Msb) # 
hist(resid(Msb)) # 

plot(Msc) # 
hist(resid(Msc)) # fit better
summary(m4)

#transformation 
Anova(Msc,type="III",test="F")

## INTER1 TWO (JE VS GO) C
str(data_b)
# Test homogeneity of variance
leveneTest(C_.~Population,data=data_b)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.16


#Check normality with shapiro test
shapiro.test(data_b$C_.) 

#P VALUE > 0.05 the assumption is not valid, p =0.0001

Msa<-lm(C_. ~  Population +Initial_biomass,data=data_b)

Msb<-lm(sqrt(C_.) ~  Population +Initial_biomass,data=data_b)
Msc<-lm(log(C_.) ~  Population +Initial_biomass,data=data_b)


plot(Msa) # 
hist(resid(Msa))

plot(Msb) # 
hist(resid(Msb)) # fit better

plot(Msc) # 
hist(resid(Msc)) # 
summary(m4)

#transformation 
Anova(Msb,type="III",test="F")

##
#INTRA1_TWO
data_c<-read.table(file="Intra2_Two.txt", header=TRUE,sep="\t",dec = ".")
str(data_c)

#BIOMASS
# Test homogeneity of variance
leveneTest(Biomass~Population,data=data_c)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.62


#Check normality with shapiro test
shapiro.test(data_c$Biomass) 

#P VALUE > 0.05 the assumption is valid, p =0.13

Msa<-lm(Biomass ~  Population +Initial_biomass,data=data_c)

Msb<-lm(sqrt(Biomass) ~  Population +Initial_biomass,data=data_c)
Msc<-lm(log(Biomass) ~  Population +Initial_biomass,data=data_c)


plot(Msa) # 
hist(resid(Msa)) # fit better

plot(Msb) # 
hist(resid(Msb)) #

plot(Msc) # 
hist(resid(Msc)) # 
summary(m4)

#transformation 
Anova(Msa,type="III",test="F")


#TOTAL LEAVES
# Test homogeneity of variance
leveneTest(Leaf_Count~Population,data=data_c)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.5


#Check normality with shapiro test
shapiro.test(data_c$Leaf_Count) 

#P VALUE > 0.05 the assumption is valid, p =0.059

Msa<-lm(Leaf_Count ~  Population +Initial_biomass,data=data_c)

Msb<-lm(sqrt(Leaf_Count) ~  Population +Initial_biomass,data=data_c)
Msc<-lm(log(Leaf_Count) ~  Population +Initial_biomass,data=data_c)


plot(Msa) # 
hist(resid(Msa)) #fit better

plot(Msb) # 
hist(resid(Msb)) #

plot(Msc) # 
hist(resid(Msc)) # 
summary(m4)

#transformation 
Anova(Msa,type="III",test="F")

#mutliple comparison
Population1<-as.factor(data_c$Population)
str(data_c)
Msa<-lm(Leaf_Count ~  Population1 +Initial_biomass,data=data_c)
gltt_Msc<-glht(Msa,linfct=mcp(Population1="Tukey"))
gltt_Msc
summary(gltt_Msc)
post<-cld(gltt_Msc, level = 0.05, decreasing = T)
post
boxplot(Leaf_Count~Population, data=data_c)

#LENGTH LONGEST LEAF

#Test homogeneity of variance
leveneTest(LLL ~Population,data=data_c)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.09


#Check normality with shapiro test
shapiro.test(data_c$LLL) 

#P VALUE > 0.05 the assumption is NOT valid, p =0.004

Msa<-lm(LLL ~  Population +Initial_biomass,data=data_c)

Msb<-lm(sqrt(LLL) ~  Population +Initial_biomass,data=data_c)
Msc<-lm(log(LLL) ~  Population +Initial_biomass,data=data_c)


plot(Msa) # 
hist(resid(Msa)) # 

plot(Msb) # 
hist(resid(Msb)) # FIT BETER

plot(Msc) # 
hist(resid(Msc)) # 
summary(m4)

#transformation 
Anova(Msb,type="III",test="F")

#SPECIFIC LEAF AREA
#Test homogeneity of variance
leveneTest(SLA ~Population,data=data_c)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.36


#Check normality with shapiro test
shapiro.test(data_c$SLA) 

#P VALUE > 0.05 the assumption is not valid, p =0.0002

Msa<-lm(SLA ~  Population +Initial_biomass,data=data_c)

Msb<-lm(sqrt(SLA) ~  Population +Initial_biomass,data=data_c)
Msc<-lm(log(SLA) ~  Population +Initial_biomass,data=data_c)


plot(Msa) # 
hist(resid(Msa)) #

plot(Msb) # 
hist(resid(Msb)) #

plot(Msc) # 
hist(resid(Msc)) #fit better
summary(m4)

#transformation 
Anova(Msc,type="III",test="F")


#LEAF WATER CONTENT
#Test homogeneity of variance
leveneTest(LWC ~Population,data=data_c)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.38


#Check normality with shapiro test
shapiro.test(data_c$LWC) 

#P VALUE > 0.05 the assumption is not valid, p =0.001

Msa<-lm(LWC ~  Population +Initial_biomass,data=data_c)

Msb<-lm(sqrt(LWC) ~  Population +Initial_biomass,data=data_c)
Msc<-lm(log(LWC) ~  Population +Initial_biomass,data=data_c)


plot(Msa) # 
hist(resid(Msa))

plot(Msb) # 
hist(resid(Msb)) #

plot(Msc) # 
hist(resid(Msc)) # fit better
summary(m4)

#transformation 
Anova(Msc,type="III",test="F")

# CARBON TO NITROGEN CONTENT
#Test homogeneity of variance
leveneTest(C_N. ~Population,data=data_c)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.43


#Check normality with shapiro test
shapiro.test(data_c$C_N.) 

#P VALUE > 0.05 the assumption is valid, p =0.55

Msa<-lm(C_N. ~  Population +Initial_biomass,data=data_c)

Msb<-lm(sqrt(C_N.) ~  Population +Initial_biomass,data=data_c)
Msc<-lm(log(C_N.) ~  Population +Initial_biomass,data=data_c)


plot(Msa) # 
hist(resid(Msa)) # 

plot(Msb) # 
hist(resid(Msb)) #

plot(Msc) # 
hist(resid(Msc)) # fit better
summary(m4)

#transformation 
Anova(Msc,type="III",test="F")

## INTRA1 TWO (JE VS GO) N
str(data_c)
# Test homogeneity of variance
leveneTest(N_.~Population,data=data_c)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.09


#Check normality with shapiro test
shapiro.test(data_c$N_.) 

#P VALUE > 0.05 the assumption is valid, p =0.28

Msa<-lm(N_. ~  Population +Initial_biomass,data=data_c)

Msb<-lm(sqrt(N_.) ~  Population +Initial_biomass,data=data_c)
Msc<-lm(log(N_.) ~  Population +Initial_biomass,data=data_c)


plot(Msa) # 
hist(resid(Msa)) #fit better

plot(Msb) # 
hist(resid(Msb)) #

plot(Msc) # 
hist(resid(Msc)) # 
summary(m4)

#transformation 
Anova(Msa,type="III",test="F")

## INTRA1 TWO (JE VS GO) C
str(data_c)
# Test homogeneity of variance
leveneTest(C_.~Population,data=data_c)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.64


#Check normality with shapiro test
shapiro.test(data_c$C_.) 

#P VALUE > 0.05 the assumption is valid, p =0.86

Msa<-lm(C_. ~  Population +Initial_biomass,data=data_c)

Msb<-lm(sqrt(C_.) ~  Population +Initial_biomass,data=data_c)
Msc<-lm(log(C_.) ~  Population +Initial_biomass,data=data_c)


plot(Msa) # 
hist(resid(Msa))

plot(Msb) # 
hist(resid(Msb)) # 

plot(Msc) # 
hist(resid(Msc)) # fit better
summary(m4)

#transformation 
Anova(Msc,type="III",test="F")

##
#INTER1_TWO-THREE

data_d<-read.table(file="Inter2_Two_One.txt", header=TRUE,sep="\t",dec = ".")
str(data_d)

#BIOMASS
# Test homogeneity of variance
leveneTest(Biomass~Population,data=data_d)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.18


#Check normality with shapiro test
shapiro.test(data_d$Biomass) 

#P VALUE > 0.05 the assumption is valid, p =0.31

Msa<-lm(Biomass ~  Population +Initial_biomass,data=data_d)

Msb<-lm(sqrt(Biomass) ~  Population +Initial_biomass,data=data_d)
Msc<-lm(log(Biomass) ~  Population +Initial_biomass,data=data_d)


plot(Msa) # 
hist(resid(Msa)) # 

plot(Msb) # 
hist(resid(Msb)) # fit better

plot(Msc) # 
hist(resid(Msc)) # 
summary(m4)

#transformation 
Anova(Msb,type="III",test="F")
#Posthoc
Population1<-as.factor(data_d$Population)
str(data_d)
Msb<-lm(Biomass ~  Population1 +Initial_biomass,data=data_d)
gltt_Msa<-glht(Msb,linfct=mcp(Population1="Tukey"))
gltt_Msa
summary(gltt_Msa)
post<-cld(gltt_Msa, level = 0.05, decreasing = F)
post
boxplot(Biomass~Population, data=data_d)

#TOTAL LEAVES
# Test homogeneity of variance
leveneTest(Leaf_Count~Population,data=data_d)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.87


#Check normality with shapiro test
shapiro.test(data_d$Leaf_Count) 

#P VALUE > 0.05 the assumption is not valid, p =0.007

Msa<-lm(Leaf_Count ~  Population +Initial_biomass,data=data_d)

Msb<-lm(sqrt(Leaf_Count) ~  Population +Initial_biomass,data=data_d)
Msc<-lm(log(Leaf_Count) ~  Population +Initial_biomass,data=data_d)


plot(Msa) # 
hist(resid(Msa)) #

plot(Msb) # 
hist(resid(Msb)) #

plot(Msc) # 
hist(resid(Msc)) # fit better
summary(m4)

#transformation 
Anova(Msc,type="III",test="F")

#LENGTH LONGEST LEAF

#Test homogeneity of variance
leveneTest(LLL ~Population,data=data_d)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.17


#Check normality with shapiro test
shapiro.test(data_c$LLL) 

#P VALUE > 0.05 the assumption is valid, p =0.004

Msa<-lm(LLL ~  Population +Initial_biomass,data=data_d)

Msb<-lm(sqrt(LLL) ~  Population +Initial_biomass,data=data_d)
Msc<-lm(log(LLL) ~  Population +Initial_biomass,data=data_d)


plot(Msa) # 
hist(resid(Msa)) # 

plot(Msb) # 
hist(resid(Msb)) # fit better

plot(Msc) # 
hist(resid(Msc)) # 
summary(m4)

#transformation 
Anova(Msc,type="III",test="F")


#LEAF WATER CONTENT
#Test homogeneity of variance
leveneTest(LWC ~Population,data=data_d)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.722


#Check normality with shapiro test
shapiro.test(data_d$LWC) 

#P VALUE > 0.05 the assumption is not valid, p =0.02

Msa<-lm(LWC ~  Population +Initial_biomass,data=data_d)

Msb<-lm(sqrt(LWC) ~  Population +Initial_biomass,data=data_d)
Msc<-lm(log(LWC) ~  Population +Initial_biomass,data=data_d)


plot(Msa) # 
hist(resid(Msa)) #

plot(Msb) # 
hist(resid(Msb)) #

plot(Msc) # 
hist(resid(Msc)) # fit better
summary(m4)

#transformation 
Anova(Msc,type="III",test="F")

##INTER1_One_TWO
data_e<-read.table(file="Inter2_One_Two.txt", header=TRUE,sep="\t",dec = ".")
str(data_e)

#BIOMASS
# Test homogeneity of variance
leveneTest(Biomass~Population,data=data_e)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.41


#Check normality with shapiro test
shapiro.test(data_e$Biomass) 

#P VALUE > 0.05 the assumption is valid, p =0.49

Msa<-lm(Biomass ~  Population +Initial_biomass,data=data_e)

Msb<-lm(sqrt(Biomass) ~  Population +Initial_biomass,data=data_e)
Msc<-lm(log(Biomass) ~  Population +Initial_biomass,data=data_e)


plot(Msa) # 
hist(resid(Msa)) # 

plot(Msb) # 
hist(resid(Msb)) # 

plot(Msc) # 
hist(resid(Msc)) # FIT BETTER
summary(m4)

#transformation 
Anova(Msc,type="III",test="F")

boxplot(Biomass~Population, data=data_e)

#TOTAL LEAVES
# Test homogeneity of variance
leveneTest(Leaf_Count~Population,data=data_e)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.07


#Check normality with shapiro test
shapiro.test(data_e$Leaf_Count) 

#P VALUE > 0.05 the assumption is valid, p =0.15

Msa<-lm(Leaf_Count ~  Population +Initial_biomass,data=data_e)

Msb<-lm(sqrt(Leaf_Count) ~  Population +Initial_biomass,data=data_e)
Msc<-lm(log(Leaf_Count) ~  Population +Initial_biomass,data=data_e)


plot(Msa) # 
hist(resid(Msa)) #

plot(Msb) # 
hist(resid(Msb)) #

plot(Msc) # 
hist(resid(Msc)) # fit better
summary(m4)

#transformation 
Anova(Msc,type="III",test="F")

#LENGTH LONGEST LEAF

#Test homogeneity of variance
leveneTest(LLL ~Population,data=data_e)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.53


#Check normality with shapiro test
shapiro.test(data_e$LLL) 

#P VALUE > 0.05 the assumption is valid, p =0.46

Msa<-lm(LLL ~  Population +Initial_biomass,data=data_e)

Msb<-lm(sqrt(LLL) ~  Population +Initial_biomass,data=data_e)
Msc<-lm(log(LLL) ~  Population +Initial_biomass,data=data_e)


plot(Msa) # 
hist(resid(Msa)) # 

plot(Msb) # 
hist(resid(Msb)) #

plot(Msc) # 
hist(resid(Msc)) # fit better
summary(m4)

#transformation 
Anova(Msc,type="III",test="F")


#LEAF WATER CONTENT
#Test homogeneity of variance
leveneTest(LWC ~Population,data=data_e)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.09


#Check normality with shapiro test
shapiro.test(data_e$LWC) 

#P VALUE > 0.05 the assumption is valid, p =0.11

Msa<-lm(LWC ~  Population +Initial_biomass,data=data_e)

Msb<-lm(sqrt(LWC) ~  Population +Initial_biomass,data=data_e)
Msc<-lm(log(LWC) ~  Population +Initial_biomass,data=data_e)


plot(Msa) # 
hist(resid(Msa)) #

plot(Msb) # 
hist(resid(Msb)) #

plot(Msc) # 
hist(resid(Msc)) # fit better
summary(m4)

#transformation 
Anova(Msc,type="III",test="F")
boxplot(LWC~Population, data=data_e)

#Intra1_Three (JE VS GO)
data_i3<-read.table(file="Intra2_Three.txt", header=TRUE,sep="\t",dec = ".")
str(data_i3)
## BIOMASS
# Test homogeneity of variance
leveneTest(Biomass~Population,data=data_i3)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.81


#Check normality with shapiro test
shapiro.test(data_i3$Biomass) 

#P VALUE > 0.05 the assumption is valid, p =0.39

Mi3a<-lm(Biomass ~  Population +Initial_biomass,data=data_i3)

Mi3b<-lm(sqrt(Biomass) ~  Population +Initial_biomass,data=data_i3)
Mi3c<-lm(log(Biomass) ~  Population +Initial_biomass,data=data_i3)


plot(Mi3a) # 
hist(resid(Mi3a))

plot(Mi3b) # 
hist(resid(Mi3b)) #

plot(Mi3c) # 
hist(resid(Mi3c)) # fit better
summary(Mi3c)

#transformation 
Anova(Mi3c,type="III",test="F")

boxplot(Biomass ~Population,data=data_i3)
## TOTAL NUMBER OF LEAVES
# Test homogeneity of variance
leveneTest(Leaf_Count~Population,data=data_i3)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.20


#Check normality with shapiro test
shapiro.test(data_i3$Leaf_Count) 

#P VALUE > 0.05 the assumption is not valid, p =0.01

Msa<-lm(Leaf_Count ~  Population +Initial_biomass,data=data_i3)

Msb<-lm(sqrt(Leaf_Count) ~  Population +Initial_biomass,data=data_i3)
Msc<-lm(log(Leaf_Count) ~  Population +Initial_biomass,data=data_i3)


plot(Msa) # 
hist(resid(Msa))

plot(Msb) # 
hist(resid(Msb)) #

plot(Msc) # 
hist(resid(Msc)) # fit better
summary(m4)

#transformation 
Anova(Msc,type="III",test="F")

#LEAF LONGEST LEAF
# Test homogeneity of variance
leveneTest(LLL~Population,data=data_i3)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.46


#Check normality with shapiro test
shapiro.test(data_i3$LLL ) 

#P VALUE > 0.05 the assumption is valid, p =0.74

Msa<-lm(LLL  ~  Population +Initial_biomass,data=data_i3)

Msb<-lm(sqrt(LLL) ~  Population +Initial_biomass,data=data_i3)
Msc<-lm(log(LLL) ~  Population +Initial_biomass,data=data_i3)


plot(Msa) # 
hist(resid(Msa))

plot(Msb) # 
hist(resid(Msb)) #

plot(Msc) # 
hist(resid(Msc)) # fit better
summary(m4)

#transformation 
Anova(Msb,type="III",test="F")

##LEAF WATER CONTENT
#Test homogeneity of variance
leveneTest(LWC~Population,data=data_i3)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.41


#Check normality with shapiro test
shapiro.test(data_i3$LWC ) 

#P VALUE > 0.05 the assumption is valid, p =0.29

Msa<-lm(LWC  ~  Population +Initial_biomass,data=data_i3)

Msb<-lm(sqrt(LWC) ~  Population +Initial_biomass,data=data_i3)
Msc<-lm(log(LWC) ~  Population +Initial_biomass,data=data_i3)


plot(Msa) # 
hist(resid(Msa))

plot(Msb) # 
hist(resid(Msb)) #

plot(Msc) # 
hist(resid(Msc)) # fit better
summary(m4)

#transformation 
Anova(Msc,type="III",test="F")

## SINGLE TO COMPETITION JENA
data_f<-read.table(file="Harvest2_JE_Fert.txt", header=TRUE,sep="\t",dec = ".")
str(data_f)


boxplot(JE_Biomass~Mixture_Code ,data=data_f)
#ABOVEGROUND BIOMASS

#Use lmer with RMEL=F
Mixture_Code1<-as.factor(data_f$Mixture_Code)
str(data_f)

# Test homogeneity of variance
leveneTest(JE_Biomass~Mixture_Code1,data=data_f)
#Assumption is valid if pvalue >.05
# assumption is not valid, p =0.009


#Check normality with shapiro test
shapiro.test(data_f$JE_Biomass) 

#P VALUE > 0.05 the assumption is not valid, p =0.001

#transformation
M2a<-lm(JE_Biomass ~  Mixture_Code1 +Initial_biomass,data=data_f)
M2b<-lm(sqrt(JE_Biomass) ~  Mixture_Code1 +Initial_biomass,data=data_f)
M2c<-lm(log(JE_Biomass) ~  Mixture_Code1 +Initial_biomass,data=data_f)


plot(M2a)  
hist(resid(M2a)) 

plot(M2b)  
hist(resid(M2b)) #fit better

plot(M2c) # 
hist(resid(M2c)) # 

#SQRT transformation fit better
Anova(M2b,type="III",test="F")
##
DunnettTest(sqrt(JE_Biomass) ~  Mixture_Code1,data=data_f)
#multiple comparaison

gltt_M2b<-glht(M2b,linfct=mcp(Mixture_Code1="Tukey"))
gltt_M2b
summary(gltt_M2b)
post<-cld(gltt_M2b, level = 0.05, decreasing = T)
post

#GO incompetition JE
data_g<-read.table(file="Harvest2_GO_Fert.txt", header=TRUE,sep="\t",dec = ".")
Mixture_Code1<-as.factor(data_g$Mixture_Code)
str(data_g)

# Test homogeneity of variance
leveneTest(GO_Biomass~Mixture_Code1,data=data_g)
#Assumption is valid if pvalue >.05
# assumption is not valid, p =0.01


#Check normality with shapiro test
shapiro.test(data_g$GO_Biomass) 

#P VALUE > 0.05 the assumption is valid, p =0.15

#transformation
M2a<-lm(GO_Biomass ~  Mixture_Code1 +Initial_biomass,data=data_g)
M2b<-lm(sqrt(GO_Biomass) ~  Mixture_Code1 +Initial_biomass,data=data_g)
M2c<-lm(log(GO_Biomass) ~  Mixture_Code1 +Initial_biomass,data=data_g)


plot(M2a)  
hist(resid(M2a)) #

plot(M2b)  
hist(resid(M2b)) # fit better

plot(M2c) # 
hist(resid(M2c)) # 

# transformation fit better
Anova(M2b,type="III",test="F")
##
DunnettTest(sqrt(GO_Biomass) ~  Mixture_Code1,data=data_g)
#multiple comparaison

gltt_M2a<-glht(M2b,linfct=mcp(Mixture_Code1="Tukey"))
gltt_M2a
summary(gltt_M2a)
post<-cld(gltt_M2a, level = 0.05, decreasing = T)
post

#Data for Box_plot 
data_box<-read.table(file="Harvest2_JE_GO_Fert.txt", header=TRUE,sep="\t",dec = ".")
str(data_box)

boxplot(Biomass ~Mixture_Code,data=data_box)


position<-c(1,2,3,4,5,6,8,9,10,11,12,13)
boxplot(Biomass ~Mixture_Code,data=data_box,axes=F,at=position ,border="gray40",xlab= "Population", ylab="Aboveground biomass (g DW)",las=2, 
        col=c("#FFFF00","#FFFF00","#FFFF00","#FFFF00","#FFFF00","#FFFF00","#0033FF","#0033FF","#0033FF","#0033FF","#0033FF","#0033FF"))
axis(2,at=seq(0,1.5,0.5), seq(0,1.5,0.5), las=1)
axis(1, at=c(1,2,3,4,5,6,8,9,10,11,12,13), labels=c("0:1","1:1","0:2","2:1","1:2","0:3","0:1","1:1","0:2","2:1","1:2","0:3"), las=1)
box()


#add mean for JE incompetition with GO
Mixture_Code1<-as.factor(data_f$Mixture_Code)
str(data_f)
M2b<-lm(sqrt(JE_Biomass) ~  Mixture_Code1 +Initial_biomass,data=data_f)								


lsm<-emmeans(M2b,~Mixture_Code1+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Mixture_Code1+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(8,9,10,11,12,13)
points(position,plot$response,pch=4,cex=2)

#add mean for GO incompetition with JE
Mixture_Code1<-as.factor(data_g$Mixture_Code)
str(data_g)
M2b<-lm(sqrt(GO_Biomass) ~  Mixture_Code1 +Initial_biomass,data=data_g)								


lsm<-emmeans(M2b,~Mixture_Code1+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm
pan.emm <- emmeans(regrid(lsm), list(~Mixture_Code1+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(1,2,3,4,5,6)
points(position,plot$response,pch=4,cex=2)

#TOTAL NUMBER OF LEAVES
boxplot(JE_Leaf_Count ~Mixture_Code ,data=data_f)
#Use lmer with RMEL=F
Mixture_Code1<-as.factor(data_f$Mixture_Code)
str(data_f)

# Test homogeneity of variance
leveneTest(JE_Leaf_Count~Mixture_Code1,data=data_f)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.89


#Check normality with shapiro test
shapiro.test(data_f$JE_Leaf_Count) 

#P VALUE > 0.05 the assumption is valid, p =0.056

#transformation
M2a<-lm(JE_Leaf_Count ~  Mixture_Code1 +Initial_biomass,data=data_f)
M2b<-lm(sqrt(JE_Leaf_Count) ~  Mixture_Code1 +Initial_biomass,data=data_f)
M2c<-lm(log(JE_Leaf_Count) ~  Mixture_Code1 +Initial_biomass,data=data_f)


plot(M2a)  
hist(resid(M2a)) 

plot(M2b)  
hist(resid(M2b)) #fit better

plot(M2c) # 
hist(resid(M2c)) # 

#SQRT transformation fit better
Anova(M2b,type="III",test="F")
##
DunnettTest(log(JE_Leaf_Count) ~  Mixture_Code1,data=data_f)
#multiple comparaison

gltt_M2b<-glht(M2c,linfct=mcp(Mixture_Code1="Tukey"))
gltt_M2b
summary(gltt_M2b)
post<-cld(gltt_M2b, level = 0.05, decreasing = T)
post

#GO incompetition JE
data_g<-read.table(file="Harvest2_GO_Fert.txt", header=TRUE,sep="\t",dec = ".")
Mixture_Code1<-as.factor(data_g$Mixture_Code)
str(data_g)

# Test homogeneity of variance
leveneTest(GO_Leaf_Count~Mixture_Code1,data=data_g)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.17


#Check normality with shapiro test
shapiro.test(data_g$GO_Leaf_Count) 

#P VALUE > 0.05 the assumption is not valid, p =0.02

#transformation
M3a<-lm(GO_Leaf_Count ~  Mixture_Code1 +Initial_biomass,data=data_g)
M3b<-lm(sqrt(GO_Leaf_Count) ~  Mixture_Code1 +Initial_biomass,data=data_g)
M3c<-lm(log(GO_Leaf_Count) ~  Mixture_Code1 +Initial_biomass,data=data_g)


plot(M3a)  
hist(resid(M3a)) #

plot(M3b)  
hist(resid(M3b)) #

plot(M3c) # 
hist(resid(M3c)) # fit better

# transformation fit better
Anova(M3c,type="III",test="F")
##
DunnettTest(log(GO_Leaf_Count) ~  Mixture_Code1,data=data_g)
#multiple comparaison

gltt_M2a<-glht(M3c,linfct=mcp(Mixture_Code1="Tukey"))
gltt_M2a
summary(gltt_M2a)
post<-cld(gltt_M2a, level = 0.05, decreasing = T)
post

#Data for Box_plot 
data_box<-read.table(file="Harvest2_JE_GO_Fert.txt", header=TRUE,sep="\t",dec = ".")
str(data_box)

boxplot(Leaf_Count ~Mixture_Code,data=data_box)


position<-c(1,2,3,4,5,6,8,9,10,11,12,13)
boxplot(Leaf_Count ~Mixture_Code,data=data_box,axes=F,at=position ,border="gray40",xlab= "Population", ylab="Total number of leaves",las=2, 
        col=c("#FFFF00","#FFFF00","#FFFF00","#FFFF00","#FFFF00","#FFFF00","#0033FF","#0033FF","#0033FF","#0033FF","#0033FF","#0033FF"))
axis(2,at=seq(6,14,2), seq(6,14,2), las=1)
axis(1, at=c(1,2,3,4,5,6,8,9,10,11,12,13), labels=c("0:1","1:1","0:2","2:1","1:2","0:3","0:1","1:1","0:2","2:1","1:2","0:3"), las=1)
box()


#add mean for JE incompetition with GO
Mixture_Code1<-as.factor(data_f$Mixture_Code)
str(data_f)
M2b<-lm(sqrt(JE_Leaf_Count) ~  Mixture_Code1 +Initial_biomass,data=data_f)									


lsm<-emmeans(M2b,~Mixture_Code1+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Mixture_Code1+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(8,9,10,11,12,13)
points(position,plot$response,pch=4,cex=2)

#add mean for GO incompetition with JE
Mixture_Code1<-as.factor(data_g$Mixture_Code)
str(data_g)
M3c<-lm(log(GO_Leaf_Count) ~  Mixture_Code1 +Initial_biomass,data=data_g)									

lsm<-emmeans(M3c,~Mixture_Code1+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm
pan.emm <- emmeans(regrid(lsm), list(~Mixture_Code1+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(1,2,3,4,5,6)
points(position,plot$response,pch=4,cex=2)


#LEAF LENGTH

boxplot(JE_LLL ~Mixture_Code ,data=data_f)
#Use lmer with RMEL=F
Mixture_Code1<-as.factor(data_f$Mixture_Code)
str(data_f)

# Test homogeneity of variance
leveneTest(JE_LLL~Mixture_Code1,data=data_f)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.22


#Check normality with shapiro test
shapiro.test(data_f$JE_LLL) 

#P VALUE > 0.05 the assumption is NOT valid, p =0.002

#transformation
M2a<-lm(JE_LLL ~  Mixture_Code1 +Initial_biomass,data=data_f)
M2b<-lm(sqrt(JE_LLL) ~  Mixture_Code1 +Initial_biomass,data=data_f)
M2c<-lm(log(JE_LLL) ~  Mixture_Code1 +Initial_biomass,data=data_f)


plot(M2a)  
hist(resid(M2a)) # 

plot(M2b)  
hist(resid(M2b)) # fit better

plot(M2c) # 
hist(resid(M2c)) # 

# untransformed fit better
Anova(M2b,type="III",test="F")
##
DunnettTest(JE_LLL ~  Mixture_Code1,data=data_f)
#multiple comparaison

gltt_M2b<-glht(M2b,linfct=mcp(Mixture_Code1="Tukey"))
gltt_M2b
summary(gltt_M2b)
post<-cld(gltt_M2b, level = 0.05, decreasing = T)
post

#GO incompetition JE
data_g<-read.table(file="Harvest1_GO_Fert.txt", header=TRUE,sep="\t",dec = ".")
Mixture_Code1<-as.factor(data_g$Mixture_Code)
str(data_g)

# Test homogeneity of variance
leveneTest(GO_LLL~Mixture_Code1,data=data_g)
#Assumption is valid if pvalue >.05
# assumption is not valid, p =0.004


#Check normality with shapiro test
shapiro.test(data_g$GO_LLL) 

#P VALUE > 0.05 the assumption is not valid, p =0.002

#transformation
M3a<-lm(GO_LLL ~  Mixture_Code1 +Initial_biomass,data=data_g)
M3b<-lm(sqrt(GO_LLL) ~  Mixture_Code1 +Initial_biomass,data=data_g)
M3c<-lm(log(GO_LLL) ~  Mixture_Code1 +Initial_biomass,data=data_g)


plot(M3a)  
hist(resid(M3a)) #

plot(M3b)  
hist(resid(M3b)) # fit better

plot(M3c) # 
hist(resid(M3c)) # 

# sqrt transformation fit better
Anova(M3c,type="III",test="F")
##
DunnettTest(log(GO_LLL) ~  Mixture_Code1,data=data_g)
#multiple comparaison

gltt_M2a<-glht(M3c,linfct=mcp(Mixture_Code1="Tukey"))
gltt_M2a
summary(gltt_M2a)
post<-cld(gltt_M2a, level = 0.05, decreasing = T)
post

#Data for Box_plot 
data_box<-read.table(file="Harvest1_JE_GO_Fert.txt", header=TRUE,sep="\t",dec = ".")
str(data_box)

boxplot(LLL ~Mixture_Code,data=data_box)


position<-c(1,2,3,4,5,6,8,9,10,11,12,13)
boxplot(LLL ~Mixture_Code,data=data_box,axes=F,at=position ,border="gray40",xlab= "Population", ylab="Length of the longest leaf (cm)",las=2, 
        col=c("#FFFF00","#FFFF00","#FFFF00","#FFFF00","#FFFF00","#FFFF00","#0033FF","#0033FF","#0033FF","#0033FF","#0033FF","#0033FF"))
axis(2,at=seq(5,20,5), seq(5,20,5), las=1)
axis(1, at=c(1,2,3,4,5,6,8,9,10,11,12,13), labels=c("0:1","1:1","0:2","2:1","1:2","0:3","0:1","1:1","0:2","2:1","1:2","0:3"), las=1)
box()



#add mean for JE incompetition with GO
Mixture_Code1<-as.factor(data_f$Mixture_Code)
str(data_f)
M2b<-lm(sqrt(JE_LLL) ~  Mixture_Code1 +Initial_biomass,data=data_f)								

lsm<-emmeans(M2b,~Mixture_Code1+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Mixture_Code1+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(8,9,10,11,12,13)
points(position,plot$response,pch=4,cex=2)

#add mean for GO incompetition with JE
Mixture_Code1<-as.factor(data_g$Mixture_Code)
str(data_g)
M3c<-lm(log(GO_LLL) ~  Mixture_Code1 +Initial_biomass,data=data_g)								

lsm<-emmeans(M3c,~Mixture_Code1+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm
pan.emm <- emmeans(regrid(lsm), list(~Mixture_Code1+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(1,2,3,4,5,6)
points(position,plot$response,pch=4,cex=2)

#SPECIFIC LEAF AREA
data_SP<-read.table(file="JE_LMA.txt", header=TRUE,sep="\t",dec = ".")
boxplot(JE_SLA ~Mixture_Code ,data=data_SP)
#Use lmer with RMEL=F
Mixture_Code1<-as.factor(data_f$Mixture_Code)
str(data_f)

# Test homogeneity of variance
leveneTest(JE_SLA~Mixture_Code1,data=data_f)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.11


#Check normality with shapiro test
shapiro.test(data_f$JE_SLA) 

#P VALUE > 0.05 the assumption is not valid, p =0.0002

#transformation
M2a<-lm(JE_SLA ~  Mixture_Code1 +Initial_biomass,data=data_SP)
M2b<-lm(sqrt(JE_SLA) ~  Mixture_Code1 +Initial_biomass,data=data_SP)
M2c<-lm(log(JE_SLA) ~  Mixture_Code1 +Initial_biomass,data=data_SP)


plot(M2a)  
hist(resid(M2a)) #fit better

plot(M2b)  
hist(resid(M2b)) #

plot(M2c) # 
hist(resid(M2c)) # 

#UN transformed fit better
Anova(M2a,type="III",test="F")
##
DunnettTest(JE_SLA ~  Mixture_Code1,data=data_SP)
#multiple comparaison

gltt_M2b<-glht(M2a,linfct=mcp(Mixture_Code1="Tukey"))
gltt_M2b
summary(gltt_M2b)
post<-cld(gltt_M2b, level = 0.05, decreasing = T)
post

#GO incompetition JE
data_SPB<-read.table(file="GO_LMA.txt", header=TRUE,sep="\t",dec = ".")
Mixture_Code1<-as.factor(data_SPB$Mixture_Code)
str(data_g)

# Test homogeneity of variance
leveneTest(GO_SLA~Mixture_Code1,data=data_SPB)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.46


#Check normality with shapiro test
shapiro.test(data_SPB$GO_SLA) 

#P VALUE > 0.05 the assumption is not valid, p =0.02

#transformation
M3a<-lm(GO_SLA ~  Mixture_Code1 +Initial_biomass,data=data_SPB)
M3b<-lm(sqrt(GO_SLA) ~  Mixture_Code1 +Initial_biomass,data=data_SPB)
M3c<-lm(log(GO_SLA) ~  Mixture_Code1 +Initial_biomass,data=data_SPB)


plot(M3a)  
hist(resid(M3a)) #

plot(M3b)  
hist(resid(M3b)) # FIT BETTER

plot(M3c) # 
hist(resid(M3c)) # 

# transformation fit better
Anova(M3b,type="III",test="F")
##
DunnettTest(log(GO_SLA) ~  Mixture_Code1,data=data_g)
#multiple comparaison

gltt_M2a<-glht(M3b,linfct=mcp(Mixture_Code1="Tukey"))
gltt_M2a
summary(gltt_M2a)
post<-cld(gltt_M2a, level = 0.05, decreasing = T)
post

#Data for Box_plot 
data_box<-read.table(file="Box_LMA.txt", header=TRUE,sep="\t",dec = ".")
str(data_box)

boxplot(SLA  ~Mixture_Code,data=data_box)


position<-c(1,2,3,5,6,7)
boxplot(SLA  ~Mixture_Code,data=data_box,axes=F,at=position ,border="gray40",xlab= "Population", ylab="Specific leaf area (mm2 mg-1)",las=2, 
        col=c("#FFFF00","#FFFF00","#FFFF00","#0033FF","#0033FF","#0033FF"))
axis(2,at=seq(6,14,2), seq(6,14,2), las=1)
axis(1, at=c(1,2,3,5,6,7), labels=c("0:1","1:1","0:2","0:1","1:1","0:2"), las=1)
box()


#add mean for JE incompetition with GO
Mixture_Code1<-as.factor(data_f$Mixture_Code)
str(data_f)
M2c<-lm(JE_SLA ~  Mixture_Code1 +Initial_biomass,data=data_f)

lsm<-emmeans(M2a,~Mixture_Code1+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Mixture_Code1+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(5,6,7)
points(position,plot$response,pch=4,cex=2)

#add mean for GO incompetition with JE
Mixture_Code1<-as.factor(data_g$Mixture_Code)
str(data_g)
M3b<-lm(sqrt(GO_SLA) ~  Mixture_Code1 +Initial_biomass,data=data_g)

lsm<-emmeans(M3b,~Mixture_Code1+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm
pan.emm <- emmeans(regrid(lsm), list(~Mixture_Code1+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(1,2,3)
points(position,plot$response,pch=4,cex=2)

#CHECK CORRELATION OF 2:1 LLL VS 2:1 SLA
### SPEARMAN's CORRELATION
library(ggplot2)
library(lmerTest)
#SLA-LLL-2:1

data_Cor_2.1<-read.table(file="SLA_LLL_2.1.txt", header=TRUE,sep="\t",dec = ".")

cor.test(data_Cor_2.1$JE_SLA,data_Cor_2.1$JE_LLL,method = "spearman")
ggscatter(data_Cor_2.1, x = "JE_SLA", y = "JE_LLL", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Specific leaf area (mm2 mg-1)", ylab = "Length of the longest leaf (cm)")


data_Cor$predlm<-predict(m5a)
ggplot(data_Cor,aes(x=JE_SLA,y=JE_LLL, color=Mixture_Code))+
  geom_point()+ geom_smooth(method = "lm") +theme_bw() + 
  ylab("\nLength of the longest leaf (cm)")+
  xlab("\nSpecific leaf area (mm2 mg-1)")+
  theme(axis.text.x =element_text(vjust = 1, hjust = 1),
        panel.grid = element_blank())+stat_cor(aes(color=Mixture_Code),method = "pearson" )

#LEAF WATER CONTENT

data_f<-read.table(file="Harvest1_JE_Fert.txt", header=TRUE,sep="\t",dec = ".")
str(data_f)

boxplot(JE_LWC ~Mixture_Code ,data=data_f)
#Use lmer with RMEL=F
Mixture_Code1<-as.factor(data_f$Mixture_Code)
str(data_f)

# Test homogeneity of variance
leveneTest(JE_LWC ~Mixture_Code1,data=data_f)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.09


#Check normality with shapiro test
shapiro.test(data_f$JE_LWC) 

#P VALUE > 0.05 the assumption is not valid, p =0.001

#transformation
M2a<-lm(JE_LWC ~  Mixture_Code1 +Initial_biomass,data=data_f)
M2b<-lm(sqrt(JE_LWC) ~  Mixture_Code1 +Initial_biomass,data=data_f)
M2c<-lm(log(JE_LWC) ~  Mixture_Code1 +Initial_biomass,data=data_f)


plot(M2a)  
hist(resid(M2a)) 

plot(M2b)  
hist(resid(M2b)) #fit better

plot(M2c) # 
hist(resid(M2c)) # 

#SQRT transformation fit better
Anova(M2b,type="III",test="F")
##
DunnettTest(sqrt(JE_LWC) ~  Mixture_Code1,data=data_f)
#multiple comparaison

gltt_M2b<-glht(M2b,linfct=mcp(Mixture_Code1="Tukey"))
gltt_M2b
summary(gltt_M2b)
post<-cld(gltt_M2b, level = 0.05, decreasing = T)
post

#GO incompetition JE
data_g<-read.table(file="Harvest1_GO_Fert.txt", header=TRUE,sep="\t",dec = ".")
Mixture_Code1<-as.factor(data_g$Mixture_Code)
str(data_g)

# Test homogeneity of variance
leveneTest(GO_LWC~Mixture_Code1,data=data_g)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.23


#Check normality with shapiro test
shapiro.test(data_g$GO_LWC) 

#P VALUE > 0.05 the assumption is valid, p =0.67

#transformation
M2a<-lm(GO_LWC ~  Mixture_Code1 +Initial_biomass,data=data_g)
M2b<-lm(sqrt(GO_LWC) ~  Mixture_Code1 +Initial_biomass,data=data_g)
M2c<-lm(log(GO_LWC) ~  Mixture_Code1 +Initial_biomass,data=data_g)


plot(M2a)  
hist(resid(M2a)) # FIT BETTER

plot(M2b)  
hist(resid(M2b)) #

plot(M2c) # 
hist(resid(M2c)) # 

# transformation fit better
Anova(M2a,type="III",test="F")
##
DunnettTest(GO_LWC ~  Mixture_Code1,data=data_g)
#multiple comparaison

gltt_M2a<-glht(M2a,linfct=mcp(Mixture_Code1="Tukey"))
gltt_M2a
summary(gltt_M2a)
post<-cld(gltt_M2a, level = 0.05, decreasing = T)
post

#Data for Box_plot 
data_box<-read.table(file="Harvest1_JE_GO_Fert.txt", header=TRUE,sep="\t",dec = ".")
str(data_box)

boxplot(LWC ~Mixture_Code,data=data_box)


position<-c(1,2,3,4,5,6,8,9,10,11,12,13)
boxplot(LWC ~Mixture_Code,data=data_box,axes=F,at=position ,border="gray40",xlab= "Population", ylab="Water content (%)",las=2, 
        col=c("#FFFF00","#FFFF00","#FFFF00","#FFFF00","#FFFF00","#FFFF00","#0033FF","#0033FF","#0033FF","#0033FF","#0033FF","#0033FF"))
axis(2,at=seq(82,92,2), seq(82,92,2), las=1)
axis(1, at=c(1,2,3,4,5,6,8,9,10,11,12,13), labels=c("0:1","1:1","0:2","2:1","1:2","0:3","0:1","1:1","0:2","2:1","1:2","0:3"), las=1)
box()


#add mean for JE incompetition with GO
Mixture_Code1<-as.factor(data_f$Mixture_Code)
str(data_f)
M2b<-lm(sqrt(JE_LWC) ~  Mixture_Code1 +Initial_biomass,data=data_f)

lsm<-emmeans(M2b,~Mixture_Code1+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Mixture_Code1+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(8,9,10,11,12,13)
points(position,plot$response,pch=4,cex=2)

#add mean for GO incompetition with JE
Mixture_Code1<-as.factor(data_g$Mixture_Code)
str(data_g)
M2a<-lm(GO_LWC ~  Mixture_Code1 +Initial_biomass,data=data_g)

lsm<-emmeans(M2a,~Mixture_Code1+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm
pan.emm <- emmeans(regrid(lsm), list(~Mixture_Code1+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(1,2,3,4,5,6)
points(position,plot$emmean,pch=4,cex=2)




#RELATIVE COMPETITION INTENSITY BASED ON CONTROL
#Biomass was considered as measure of yield

#RCI_CONT_TWO
dat_1a<-read.table(file="RCI_Cont_Two.txt", header=TRUE,sep="\t",dec = ".")

str(dat_1a)
#Use lmer with RMEL=F
MC <- lmer(Competitive_ability_JE  ~ Population+(1 | Mixture_Code)+Initial_biomass, REML=F, data = dat_1a)
summary(MC)


#Check themodel(DHARMarequired)
sim1a<-simulateResiduals(fittedModel= MC,n=1000) #fit better
plot(sim1a, quantreg= FALSE)


summary(MC)
#untransformed is better
#Significant test
Anova(MC,type="III",test="Chisq")

#multiple comparaison
library(multcomp)
gltt_M1<-glht(MC,linfct=mcp(Population="Tukey"))
gltt_M1
summary(gltt_M1)
post<-cld(gltt_M1, level = 0.05, decreasing = F)
post
#

#RCI_CONT_THREE
dat_1b<-read.table(file="RCI_Cont_Three.txt", header=TRUE,sep="\t",dec = ".")

str(dat_1b)
#Use lmer with RMEL=F
MB <- lmer(Competitive_ability_JE  ~ Population+(1 | Mixture_Code)+Initial_biomass, REML=F, data = dat_1b)
summary(MC)


#Check themodel(DHARMarequired)
sim1a<-simulateResiduals(fittedModel= MB,n=1000) #fit better
plot(sim1a, quantreg= FALSE)


summary(MB)
#untransformed is better
#Significant test
Anova(MB,type="III",test="Chisq")

#multiple comparaison
library(multcomp)
gltt_M1<-glht(MB,linfct=mcp(Population="Tukey"))
gltt_M1
summary(gltt_M1)
post<-cld(gltt_M1, level = 0.05, decreasing = T)
post
#


#RCI BOX PLOT TWO AND THREE (JE VS GO) IN SUPPLEMENT
#Data for Box_plot 
dat_4b<-read.table(file="RCI_Box.txt", header=TRUE,sep="\t",dec = ".")
str(dat_4b)

boxplot(Competitive_ability_JE ~Population,data=dat_4b)


position<-c(1,2,4,5)
boxplot(Competitive_ability_JE ~Population,data=dat_4b,axes=F,at=position ,border="gray40",xlab= "Population", ylab="Relative competitive intensity",las=2, 
        col=c("#FFFF00","#0033FF","#FFFF00","#0033FF"))
axis(2,at=seq(0.0,0.8,0.2), seq(-0.0,0.8,0.2), las=1)
axis(1, at=c(1,2,4,5), labels=c("2_GO","2_JE","3_GO","3_JE"), las=1)
box()



#add mean for two competition

MC <- lmer(Competitive_ability_JE  ~ Population+(1 | Mixture_Code)+Initial_biomass, REML=F, data = dat_1a)

lsm<-emmeans(MC,~Population+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Population+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(1,2)
points(position,plot$emmean,pch=4,cex=2)

#add mean for three competition

str(dat_2b)
MB <- lmer(Competitive_ability_JE  ~ Population+(1 | Mixture_Code)+Initial_biomass, REML=F, data = dat_1b)

lsm<-emmeans(MB,~Population+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Population+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(4,5)
points(position,plot$emmean,pch=4,cex=2)
##


#RCI_CONT (Compare GO vs JE based on similar control)
dat_1c<-read.table(file="RCI_Cont.txt", header=TRUE,sep="\t",dec = ".")

str(dat_1c)
#Use lmer with RMEL=F
MD <- lmer(Competitive_ability_JE  ~ Population+(1 | Mixture_Code)+Initial_biomass, REML=F, data = dat_1c)
summary(MD)


#Check themodel(DHARMarequired)
sim1a<-simulateResiduals(fittedModel= MD,n=1000) #fit better
plot(sim1a, quantreg= FALSE)

#untransformed is better
#Significant test
Anova(MD,type="III",test="Chisq")


#COMPETITION CONTROL JE

# use lm   

dat_2a<-read.table(file="Competition_Cont_JE.txt", header=TRUE,sep="\t",dec = ".")
Mixture_Code1<-as.factor(dat_2a$Mixture_Code)
str(dat_2a)
# lm


# Test homogeneity of variance
leveneTest(Competitive_ability_JE~Mixture_Code,data=dat_2a)
#Assumption is valid if pvalue >.05
# assumption is not valid, p =0.0004


#Check normality with shapiro test
shapiro.test(dat_2a$Competitive_ability_JE) 

#P VALUE > 0.05 the assumption is valid, p =0.976
M2C<-lm(Competitive_ability_JE ~  Mixture_Code1 +Initial_biomass,data=dat_2a)
M2B<-lm(log(Competitive_ability_JE) ~  Mixture_Code1 +Initial_biomass,data=dat_2a)
M2D<-lm(sqrt(Competitive_ability_JE) ~  Mixture_Code1 +Initial_biomass,data=dat_2a)
#LM
plot(M2C)  
hist(resid(M2C)) 

plot(M2B)  
hist(resid(M2B))

plot(M2D)  
hist(resid(M2D)) # fit better


#untransformed fit better
Anova(M2D,type="III",test="F")
##
#Multiple comparison
gltt_M2a<-glht(M2D,linfct=mcp(Mixture_Code1="Tukey"))
gltt_M2a
summary(gltt_M2a)
post<-cld(gltt_M2a, level = 0.05, decreasing = F)
post

#COMPETITION CONTROL GO

# use lm   

dat_2b<-read.table(file="Competition_Cont_GO.txt", header=TRUE,sep="\t",dec = ".")
Mixture_Code1<-as.factor(dat_2b$Mixture_Code)
str(dat_2b)
# lm


# Test homogeneity of variance
leveneTest(Competitive_ability_JE~Mixture_Code,data=dat_2b)
#Assumption is valid if pvalue >.05
# assumption is not valid, p =0.0001


#Check normality with shapiro test
shapiro.test(dat_2b$Competitive_ability_JE) 

#P VALUE > 0.05 the assumption is valid, p =0.867
M2C<-lm(Competitive_ability_JE ~  Mixture_Code1 +Initial_biomass,data=dat_2b)
M2B<-lm(log(Competitive_ability_JE) ~  Mixture_Code1 +Initial_biomass,data=dat_2b)
M2D<-lm(sqrt(Competitive_ability_JE) ~  Mixture_Code1 +Initial_biomass,data=dat_2b)
#LM
plot(M2C)  
hist(resid(M2C)) #fit better

plot(M2B)  
hist(resid(M2B))

plot(M2D)  
hist(resid(M2D)) # 


#untransformed fit better
Anova(M2C,type="III",test="F")
##
#Multiple comparison
gltt_M2a<-glht(M2C,linfct=mcp(Mixture_Code1="Tukey"))
gltt_M2a
summary(gltt_M2a)
post<-cld(gltt_M2a, level = 0.05, decreasing = F)
post

#COMPETITION_CONTROL_1:1

# use lm   

dat_2c<-read.table(file="Competition_Cont_1.1.txt", header=TRUE,sep="\t",dec = ".")
Mixture_Code1<-as.factor(dat_2c$Mixture_Code)
str(dat_2c)
# lm


# Test homogeneity of variance
leveneTest(Competitive_ability_JE~Mixture_Code,data=dat_2c)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.8


#Check normality with shapiro test
shapiro.test(dat_2c$Competitive_ability_JE) 

#P VALUE > 0.05 the assumption is valid, p =0.80
M2C<-lm(Competitive_ability_JE ~  Mixture_Code1 +Initial_biomass,data=dat_2c)
M2B<-lm(log(Competitive_ability_JE) ~  Mixture_Code1 +Initial_biomass,data=dat_2c)
M2D<-lm(sqrt(Competitive_ability_JE) ~  Mixture_Code1 +Initial_biomass,data=dat_2c)
#LM
plot(M2C)  
hist(resid(M2C)) #fit better

plot(M2B)  
hist(resid(M2B))

plot(M2D)  
hist(resid(M2D)) # 


#untransformed fit better
Anova(M2C,type="III",test="F")
##
#Multiple comparison
gltt_M2a<-glht(M2C,linfct=mcp(Mixture_Code1="Tukey"))
gltt_M2a
summary(gltt_M2a)
post<-cld(gltt_M2a, level = 0.05, decreasing = T)
post

# COMPETITION CONTROL 0:2
dat_2d<-read.table(file="Competition_Cont_0.2.txt", header=TRUE,sep="\t",dec = ".")
Mixture_Code1<-as.factor(dat_2d$Mixture_Code)
str(dat_2d)
# lm


# Test homogeneity of variance
leveneTest(Competitive_ability_JE~Mixture_Code,data=dat_2d)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.95


#Check normality with shapiro test
shapiro.test(dat_2d$Competitive_ability_JE) 

#P VALUE > 0.05 the assumption is valid, p =0.93
M2C<-lm(Competitive_ability_JE ~  Mixture_Code1 +Initial_biomass,data=dat_2d)
M2B<-lm(log(Competitive_ability_JE) ~  Mixture_Code1 +Initial_biomass,data=dat_2d)
M2D<-lm(sqrt(Competitive_ability_JE) ~  Mixture_Code1 +Initial_biomass,data=dat_2d)
#LM
plot(M2C)  
hist(resid(M2C)) #fit better

plot(M2B)  
hist(resid(M2B))

plot(M2D)  
hist(resid(M2D)) # 


#untransformed fit better
Anova(M2C,type="III",test="F")
##
#Multiple comparison
gltt_M2a<-glht(M2C,linfct=mcp(Mixture_Code1="Tukey"))
gltt_M2a
summary(gltt_M2a)
post<-cld(gltt_M2a, level = 0.05, decreasing = T)
post

#COMPETITION CONTROL 2:1
dat_2e<-read.table(file="Competition_Cont_2.1.txt", header=TRUE,sep="\t",dec = ".")
Mixture_Code1<-as.factor(dat_2e$Mixture_Code)
str(dat_2e)
# lm


# Test homogeneity of variance
leveneTest(Competitive_ability_JE~Mixture_Code,data=dat_2e)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.11


#Check normality with shapiro test
shapiro.test(dat_2e$Competitive_ability_JE) 

#P VALUE > 0.05 the assumption is valid, p =0.58
M2C<-lm(Competitive_ability_JE ~  Mixture_Code1 +Initial_biomass,data=dat_2e)
M2B<-lm(log(Competitive_ability_JE) ~  Mixture_Code1 +Initial_biomass,data=dat_2e)
M2D<-lm(sqrt(Competitive_ability_JE) ~  Mixture_Code1 +Initial_biomass,data=dat_2e)
#LM
plot(M2C)  
hist(resid(M2C)) #

plot(M2B)  
hist(resid(M2B))

plot(M2D)  
hist(resid(M2D)) # fit better


#untransformed fit better
Anova(M2D,type="III",test="F")
##
#Multiple comparison
gltt_M2a<-glht(M2D,linfct=mcp(Mixture_Code1="Tukey"))
gltt_M2a
summary(gltt_M2a)
post<-cld(gltt_M2a, level = 0.05, decreasing = T)
post

#COMPETITION CONTROL_1:2

dat_2f<-read.table(file="Competition_Cont_1.2.txt", header=TRUE,sep="\t",dec = ".")
Mixture_Code1<-as.factor(dat_2f$Mixture_Code)
str(dat_2f)
# lm


# Test homogeneity of variance
leveneTest(Competitive_ability_JE~Mixture_Code,data=dat_2f)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.87


#Check normality with shapiro test
shapiro.test(dat_2f$Competitive_ability_JE) 

#P VALUE > 0.05 the assumption is valid, p =0.59
M2C<-lm(Competitive_ability_JE ~  Mixture_Code1 +Initial_biomass,data=dat_2f)
M2B<-lm(log(Competitive_ability_JE) ~  Mixture_Code1 +Initial_biomass,data=dat_2f)
M2D<-lm(sqrt(Competitive_ability_JE) ~  Mixture_Code1 +Initial_biomass,data=dat_2f)
#LM
plot(M2C)  
hist(resid(M2C)) # fit better

plot(M2B)  
hist(resid(M2B))

plot(M2D)  
hist(resid(M2D)) # 

#untransformed fit better
Anova(M2C,type="III",test="F")
##
#Multiple comparison
gltt_M2a<-glht(M2D,linfct=mcp(Mixture_Code1="Tukey"))
gltt_M2a
summary(gltt_M2a)
post<-cld(gltt_M2a, level = 0.05, decreasing = T)
post

#COMPETITION CONTROL_0:3

dat_2g<-read.table(file="Competition_Cont_0.3.txt", header=TRUE,sep="\t",dec = ".")
Mixture_Code1<-as.factor(dat_2g$Mixture_Code)
str(dat_2g)
# lm


# Test homogeneity of variance
leveneTest(Competitive_ability_JE~Mixture_Code,data=dat_2g)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.96


#Check normality with shapiro test
shapiro.test(dat_2g$Competitive_ability_JE) 

#P VALUE > 0.05 the assumption is valid, p =0.054
M2C<-lm(Competitive_ability_JE ~  Mixture_Code1 +Initial_biomass,data=dat_2g)
M2B<-lm(log(Competitive_ability_JE) ~  Mixture_Code1 +Initial_biomass,data=dat_2g)
M2D<-lm(sqrt(Competitive_ability_JE) ~  Mixture_Code1 +Initial_biomass,data=dat_2g)
#LM
plot(M2C)  
hist(resid(M2C)) # 

plot(M2B)  
hist(resid(M2B))

plot(M2D)  
hist(resid(M2D)) # fit better

#untransformed fit better
Anova(M2D,type="III",test="F")
##
#Multiple comparison
gltt_M2a<-glht(M2D,linfct=mcp(Mixture_Code1="Tukey"))
gltt_M2a
summary(gltt_M2a)
post<-cld(gltt_M2a, level = 0.05, decreasing = T)
post

#Data for Box_plot 
dat_4a<-read.table(file="Box_Cont.txt", header=TRUE,sep="\t",dec = ".")
str(dat_4a)

boxplot(Competitive_ability_JE ~Mixture_Code,data=dat_4a)


position<-c(1,2,3,4,5,7,8,9,10,11)
boxplot(Competitive_ability_JE ~Mixture_Code,data=dat_4a,axes=F,at=position ,border="gray40",xlab= "Population", ylab="Relative competitive intensity",las=2, 
        col=c("#FFFF00","#FFFF00","#FFFF00","#FFFF00","#FFFF00","#0033FF","#0033FF","#0033FF","#0033FF","#0033FF"))
axis(2,at=seq(0.0,0.8,0.2), seq(-0.0,0.8,0.2), las=1)
axis(1, at=c(1,2,3,4,5,7,8,9,10,11), labels=c("1:1","0:2","2:1","1:2","0:3","1:1","0:2","2:1","1:2","0:3"), las=1)
box()



#add mean for JE incompetition with GO
Mixture_Code1<-as.factor(dat_2a$Mixture_Code)
str(dat_2a)
M2D<-lm(sqrt(Competitive_ability_JE) ~  Mixture_Code1 +Initial_biomass,data=dat_2a)			

lsm<-emmeans(M2D,~Mixture_Code1+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Mixture_Code1+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(7,8,9,10,11)
points(position,plot$response,pch=4,cex=2)

#add mean for GO incompetition with JE
Mixture_Code1<-as.factor(dat_2b$Mixture_Code)
str(dat_2b)
M2C<-lm(Competitive_ability_JE ~  Mixture_Code1 +Initial_biomass,data=dat_2b)			


lsm<-emmeans(M2C,~Mixture_Code1+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Mixture_Code1+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(1,2,3,4,5)
points(position,plot$emmean,pch=4,cex=2)
##

#RELATIVE COMPETITIVE ABILITY BASED ON MONOCULTURE

#Biomass was considered as measure of yield

#BASED ON MONOCULTURE OF TWO
dat_2c2<-read.table(file="Competition_Mon_JE_GO_1.1.txt", header=TRUE,sep="\t",dec = ".")

str(dat_2c2)

# use lm for JE incompetition GO  

str(dat_2c2)
# lm



# Test homogeneity of variance
leveneTest(Competitive_ability_JE~Population,data=dat_2c2)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.64


#Check normality with shapiro test
shapiro.test(dat_2c2$Competitive_ability_JE) 

#P VALUE > 0.05 the assumption is valid, p =0.76

M2C<-lm(Competitive_ability_JE ~  Population +Initial_biomass,data=dat_2c2)
M2B<-lm(log(Competitive_ability_JE) ~  Population +Initial_biomass,data=dat_2c2)
M2D<-lm(sqrt(Competitive_ability_JE) ~  Population +Initial_biomass,data=dat_2c2)
#LM
plot(M2C)  
hist(resid(M2C)) #fit better
plot(M2B)  
hist(resid(M2B))
plot(M2D)  
hist(resid(M2D))


#untransformed fit better
Anova(M2C,type="III",test="F")
##
#multiple comparaison
Mixture_Code1<-as.factor(dat_2c2$Mixture_Code)
str(dat_2c2)
M2C<-lm(Competitive_ability_JE ~  Mixture_Code1 +Initial_biomass,data=dat_2c2)

gltt_M2a<-glht(M2C,linfct=mcp(Mixture_Code1="Tukey"))
gltt_M2a
summary(gltt_M2a)
post<-cld(gltt_M2a, level = 0.05, decreasing = F)
post
#BASED ON MONOCULTURE OF THREE _JE
dat_3c2<-read.table(file="JE_Mon_Competition.txt", header=TRUE,sep="\t",dec = ".")
Mixture_Code1<-as.factor(dat_3c2$Mixture_series)
str(dat_3c2)
M3C<-lm(Competitive_ability_JE  ~  Mixture_Code1 +Initial_biomass,data=dat_3c2)
# Test homogeneity of variance
leveneTest(Competitive_ability_JE ~Mixture_Code,data=dat_3c2)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.13

#Check normality with shapiro test
shapiro.test(dat_3c2$Competitive_ability_JE ) 

#P VALUE > 0.05 the assumption is valid, p =0.42


plot(M3C) # 
hist(resid(M3C)) #


summary(M3)

#transformed data fit better
Anova(M3C,type="III",test="F")
##

#multiple comparaison
Mixture_Code1<-as.factor(dat_3c2$Mixture_Code)
str(dat_3c2)
M3C<-lm(Competitive_ability_JE  ~  Mixture_Code1 +Initial_biomass,data=dat_3c2)

gltt_M3C<-glht(M3C,linfct=mcp(Mixture_Code1="Tukey"))
gltt_M3C
summary(gltt_M3C)
post<-cld(gltt_M3C, level = 0.05, decreasing = F)
post

##BASED ON MONOCULTURE OF THREE _GO
dat_3c3<-read.table(file="GO_Mon_Competition.txt", header=TRUE,sep="\t",dec = ".")
Mixture_Code1<-as.factor(dat_3c3$Mixture_series)
str(dat_3c3)

# Test homogeneity of variance
leveneTest(Competitive_ability_JE ~Mixture_Code,data=dat_3c3)
#Assumption is valid if pvalue >.05
# assumption is NOT valid, p =0.0006

#Check normality with shapiro test
shapiro.test(dat_3c3$Competitive_ability_JE ) 

#P VALUE > 0.05 the assumption is valid, p =0.11

M3C<-lm(Competitive_ability_JE  ~  Mixture_Code1 +Initial_biomass,data=dat_3c3)
plot(M3C) # 
hist(resid(M3C)) #


summary(M3)

#transformed data fit better
Anova(M3C,type="III",test="F")
##

#multiple comparaison
Mixture_Code1<-as.factor(dat_3c3$Mixture_Code)
str(dat_3c3)
M3C<-lm(Competitive_ability_JE  ~  Mixture_Code1 +Initial_biomass,data=dat_3c3)

gltt_M3C<-glht(M3C,linfct=mcp(Mixture_Code1="Tukey"))
gltt_M3C
summary(gltt_M3C)
post<-cld(gltt_M3C, level = 0.05, decreasing = F)
post

##BASED ON MONOCULTURE OF THREE _JE_GO_2:1
dat_3c4<-read.table(file="Competition_Mon_JE_GO_2.1.txt", header=TRUE,sep="\t",dec = ".")
Mixture_Code1<-as.factor(dat_3c4$Mixture_series)
str(dat_3c4)

# Test homogeneity of variance
leveneTest(Competitive_ability_JE ~Mixture_Code,data=dat_3c4)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.08

#Check normality with shapiro test
shapiro.test(dat_3c4$Competitive_ability_JE ) 

#P VALUE > 0.05 the assumption is valid, p =0.58

M3C<-lm(Competitive_ability_JE  ~  Population +Initial_biomass,data=dat_3c4)
plot(M3C) # 
hist(resid(M3C)) #


summary(M3)

#transformed data fit better
Anova(M3C,type="III",test="F")
##

#multiple comparaison
Mixture_Code1<-as.factor(dat_3c4$Population)
str(dat_3c4)
M3C<-lm(Competitive_ability_JE  ~  Mixture_Code1 +Initial_biomass,data=dat_3c4)

gltt_M3C<-glht(M3C,linfct=mcp(Mixture_Code1="Tukey"))
gltt_M3C
summary(gltt_M3C)
post<-cld(gltt_M3C, level = 0.05, decreasing = T)
post


##BASED ON MONOCULTURE OF THREE _JE_GO_1:2
dat_3c5<-read.table(file="Competition_Mon_JE_GO_1.2.txt", header=TRUE,sep="\t",dec = ".")
Mixture_Code1<-as.factor(dat_3c5$Mixture_Code)
str(dat_3c5)

# Test homogeneity of variance
leveneTest(Competitive_ability_JE ~Mixture_Code,data=dat_3c5)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.95

#Check normality with shapiro test
shapiro.test(dat_3c5$Competitive_ability_JE ) 

#P VALUE > 0.05 the assumption is valid, p =0.63

M3C<-lm(Competitive_ability_JE  ~  Population +Initial_biomass,data=dat_3c5)
plot(M3C) # 
hist(resid(M3C)) #


summary(M3)

#transformed data fit better
Anova(M3C,type="III",test="F")
##

#multiple comparaison
Mixture_Code1<-as.factor(dat_3c5$Population)
str(dat_3c4)
M3C<-lm(Competitive_ability_JE  ~  Mixture_Code1 +Initial_biomass,data=dat_3c5)

gltt_M3C<-glht(M3C,linfct=mcp(Mixture_Code1="Tukey"))
gltt_M3C
summary(gltt_M3C)
post<-cld(gltt_M3C, level = 0.05, decreasing = T)
post



#Data for Box_plot 
dat_4ac<-read.table(file="Box_Mon.txt", header=TRUE,sep="\t",dec = ".")
str(dat_4ac)

boxplot(Competitive_ability_JE ~Mixture_Code,data=dat_4ac)


position<-c(1,2,3.5,4.5,6.5,7.5)
boxplot(Competitive_ability_JE ~Mixture_Code,data=dat_4ac,axes=F,at=position ,border="gray40",xlab= "Population", ylab="Relative competitive intensity",las=2, 
        col=c("#FFFF00","#FFFF00","#FFFF00","#0033FF","#0033FF","#0033FF"))
axis(2,at=seq(-1.5,0.5,0.5), seq(-1.5,0.5,0.5), las=1)
axis(1, at=c(1,2,3.5,4.5,6.5,7.5), labels=c("1:1","2:1","1:2","1:1","2:1","1:2"), las=1)
box()


#add mean for Monoculture in two competition
Mixture_Code1<-as.factor(dat_2c2$Mixture_Code)
str(dat_2c2)
M2C<-lm(Competitive_ability_JE ~  Mixture_Code1 +Initial_biomass,data=dat_2c2)						


lsm<-emmeans(M2C,~Mixture_Code1+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Mixture_Code1+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(1,4.5)
points(position,plot$emmean,pch=4,cex=2)

#add mean for Monoculture for GO_Mon_Competition
Mixture_Code1<-as.factor(dat_3c2$Mixture_Code)
str(dat_3c2)
M3C<-lm(Competitive_ability_JE  ~  Mixture_Code1 +Initial_biomass,data=dat_3c2)						


lsm<-emmeans(M3C,~Mixture_Code1+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Mixture_Code1+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(2,3.5)
points(position,plot$emmean,pch=4,cex=2)
##

#add mean for Monoculture for JE_Mon_Competition
Mixture_Code1<-as.factor(dat_3c3$Mixture_Code)
str(dat_3c3)
M3C<-lm(Competitive_ability_JE  ~  Mixture_Code1 +Initial_biomass,data=dat_3c3)						


lsm<-emmeans(M3C,~Mixture_Code1+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Mixture_Code1+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(6.5,7.5)
points(position,plot$emmean,pch=4,cex=2)
##

###COMPETITION RCI MONOCULTURE THREE (JE VS GO)
dat_4ad<-read.table(file="RCI_Mon_Three.txt", header=TRUE,sep="\t",dec = ".")
str(dat_4ad)
#Use lmer with RMEL=F
MC <- lmer(Competitive_ability_JE  ~ Population+(1 | Mixture_Code)+Initial_biomass, REML=F, data = dat_4ad)

#Check themodel(DHARMarequired)
sim1a<-simulateResiduals(fittedModel= MC,n=1000) #fit better
plot(sim1a, quantreg= FALSE)


summary(MC)
#untransformed is better
#Significant test
Anova(MC,type="III",test="Chisq")

#multiple comparaison
library(multcomp)
gltt_M1<-glht(MC,linfct=mcp(Population="Tukey"))
gltt_M1
summary(gltt_M1)
post<-cld(gltt_M1, level = 0.05, decreasing = T)
post
#

#Data for Box_plot 
dat_4ae<-read.table(file="Box_Mon_Two_Three.txt", header=TRUE,sep="\t",dec = ".")
str(dat_4ad)

boxplot(Competitive_ability_JE ~Population,data=dat_4ae)


position<-c(1,2,4,5)
boxplot(Competitive_ability_JE ~Population,data=dat_4ae,axes=F,at=position ,border="gray40",xlab= "Population", ylab="Relative competitive intensity",las=2, 
        col=c("#FFFF00","#0033FF","#FFFF00","#0033FF"))
axis(2,at=seq(-1.5,0.5,0.5), seq(-1.5,0.5,0.5), las=1)
axis(1, at=c(1,2,4,5), labels=c("GO_2","JE_2","GO_3","JE_3"), las=1)
box()

#add mean for Monoculture in two competition JE VS GO
Mixture_Code1<-as.factor(dat_2c2$Population)
str(dat_2c2)
M2C<-lm(Competitive_ability_JE ~  Mixture_Code1 +Initial_biomass,data=dat_2c2)			


lsm<-emmeans(M2C,~Mixture_Code1+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Mixture_Code1+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(1,2)
points(position,plot$emmean,pch=4,cex=2)

#add mean for Monoculture in three competition JE VS GO
Mixture_Code1<-as.factor(dat_4ad$Population)
str(dat_4ad)
MC <- lmer(Competitive_ability_JE  ~ Mixture_Code1+(1 | Mixture_Code)+Initial_biomass, REML=F, data = dat_4ad)


lsm<-emmeans(MC,~Mixture_Code1+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Mixture_Code1+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(4,5)
points(position,plot$emmean,pch=4,cex=2)

#ROOT BIOMASS RATIO (Below biomass / Aboveground biomass 2 harvest)


# use lm for JE ONLY  
str(dat_2r)
dat_2r<-read.table(file="Root_JE_H2.txt", header=TRUE,sep="\t",dec = ".")
Mixture_Code1<-as.factor(dat_2r$Mixture_series)
str(dat_2r)
# lm
M2C<-lm(Root_biomass_fraction_. ~  Mixture_Code1 +Initial_biomass,data=dat_2r)
M2Cb<-lm(sqrt(Root_biomass_fraction_.) ~  Mixture_Code1 +Initial_biomass,data=dat_2r)
M2Cc<-lm(log(Root_biomass_fraction_.) ~  Mixture_Code1 +Initial_biomass,data=dat_2r)


# Test homogeneity of variance
leveneTest(Root_biomass_fraction_.~Mixture_Code,data=dat_2r)
#Assumption is valid if pvalue >.05
# assumption is not valid, p =0.04


#Check normality with shapiro test
shapiro.test(dat_2r$Root_biomass_fraction_.) 

#P VALUE > 0.05 the assumption is valid, p =0.003

#LM
plot(M2C)  
hist(resid(M2C)) 

plot(M2Cb)  
hist(resid(M2Cb))

plot(M2Cc)  
hist(resid(M2Cc))

#log transformed fit better
Anova(M2Cc,type="III",test="F")
##
#multiple comparaison
Mixture_Code1<-as.factor(dat_2r$Mixture_Code)
str(dat_2r)
M2Cc<-lm(log(Root_biomass_fraction_.) ~  Mixture_Code1 +Initial_biomass,data=dat_2r)

gltt_M2a<-glht(M2Cc,linfct=mcp(Mixture_Code1="Tukey"))
gltt_M2a
summary(gltt_M2a)
post<-cld(gltt_M2a, level = 0.05, decreasing = F)
post
#USE LM GO ONLY
dat_3r<-read.table(file="Root_GO_H2.txt", header=TRUE,sep="\t",dec = ".")
Mixture_Code1<-as.factor(dat_3r$Mixture_series)
str(dat_3r)
M3C<-lm(Root_biomass_fraction_.  ~  Mixture_Code1 +Initial_biomass,data=dat_3r)
M3Cb<-lm(sqrt(Root_biomass_fraction_.)  ~  Mixture_Code1 +Initial_biomass,data=dat_3r)
M3Cc<-lm(log(Root_biomass_fraction_.)  ~  Mixture_Code1 +Initial_biomass,data=dat_3r)



# Test homogeneity of variance
leveneTest(Root_biomass_fraction_. ~Mixture_Code,data=dat_3r)
#Assumption is valid if pvalue >.05
# assumption is not valid, p =0.57

#Check normality with shapiro test
shapiro.test(dat_3r$Root_biomass_fraction_. ) 

#P VALUE > 0.05 the assumption is valid, p =0.005


plot(M3C) # 
hist(resid(M3C)) #

plot(M3Cb) # 
hist(resid(M3Cb)) 
plot(M3Cc) # 
hist(resid(M3Cc)) #fit

summary(M3)

#log transformed data fit better
Anova(M3Cc,type="III",test="F")
##

#multiple comparaison
Mixture_Code1<-as.factor(dat_3r$Mixture_Code)
str(dat_3r)
M3Cc<-lm(log(Root_biomass_fraction_.)  ~  Mixture_Code1 +Initial_biomass,data=dat_3r)

gltt_M3C<-glht(M3Cc,linfct=mcp(Mixture_Code1="Tukey"))
gltt_M3C
summary(gltt_M3C)
post<-cld(gltt_M3C, level = 0.05, decreasing = F)
post

#Single (JE VS GO)

dat_2r<-read.table(file="Root_JE_H2.txt", header=TRUE,sep="\t",dec = ".")
Mixture_Code1<-as.factor(dat_2r$Mixture_series)
str(dat_2r)
# lm
M2C<-lm(Root_biomass_fraction_. ~  Mixture_Code1 +Initial_biomass,data=dat_2r)
M2Cb<-lm(sqrt(Root_biomass_fraction_.) ~  Mixture_Code1 +Initial_biomass,data=dat_2r)
M2Cc<-lm(log(Root_biomass_fraction_.) ~  Mixture_Code1 +Initial_biomass,data=dat_2r)


# Test homogeneity of variance
leveneTest(Root_biomass_fraction_.~Mixture_Code,data=dat_2r)
#Assumption is valid if pvalue >.05
# assumption is not valid, p =0.04


#Check normality with shapiro test
shapiro.test(dat_2r$Root_biomass_fraction_.) 

#P VALUE > 0.05 the assumption is valid, p =0.003

#LM
plot(M2C)  
hist(resid(M2C)) 

plot(M2Cb)  
hist(resid(M2Cb))

plot(M2Cc)  
hist(resid(M2Cc))

#log transformed fit better
Anova(M2Cc,type="III",test="F")
##
#multiple comparaison
Mixture_Code1<-as.factor(dat_2r$Mixture_Code)
str(dat_2r)
M2Cc<-lm(log(Root_biomass_fraction_.) ~  Mixture_Code1 +Initial_biomass,data=dat_2r)

gltt_M2a<-glht(M2Cc,linfct=mcp(Mixture_Code1="Tukey"))
gltt_M2a
summary(gltt_M2a)
post<-cld(gltt_M2a, level = 0.05, decreasing = F)
post

#INTRA 0:2 (JE VS GO)
dat_3r<-read.table(file="Root_GO_H2.txt", header=TRUE,sep="\t",dec = ".")
Mixture_Code1<-as.factor(dat_3r$Mixture_series)
str(dat_3r)
M3C<-lm(Root_biomass_fraction_.  ~  Mixture_Code1 +Initial_biomass,data=dat_3r)
M3Cb<-lm(sqrt(Root_biomass_fraction_.)  ~  Mixture_Code1 +Initial_biomass,data=dat_3r)
M3Cc<-lm(log(Root_biomass_fraction_.)  ~  Mixture_Code1 +Initial_biomass,data=dat_3r)



# Test homogeneity of variance
leveneTest(Root_biomass_fraction_. ~Mixture_Code,data=dat_3r)
#Assumption is valid if pvalue >.05
# assumption is not valid, p =0.57

#Check normality with shapiro test
shapiro.test(dat_3r$Root_biomass_fraction_. ) 

#P VALUE > 0.05 the assumption is valid, p =0.005


plot(M3C) # 
hist(resid(M3C)) #

plot(M3Cb) # 
hist(resid(M3Cb)) 
plot(M3Cc) # 
hist(resid(M3Cc)) #fit

summary(M3)

#log transformed data fit better
Anova(M3Cc,type="III",test="F")
##

#Data for Box_plot 
dat_4r<-read.table(file="Root_Box_H2.txt", header=TRUE,sep="\t",dec = ".")
str(dat_4r)

boxplot(Root_biomass_fraction_. ~Mixture_Code,data=dat_4r)


position<-c(1,2,3.5,4.5,6,7)
boxplot(Root_biomass_fraction_. ~Mixture_Code,data=dat_4r,axes=F,at=position ,border="gray40",xlab= "Population", ylab="Root biomass fraction [%]",las=2, 
        col=c("#0033FF","#0033FF","#0033FF","#FFFF00","#FFFF00","#FFFF00"))
axis(2,at=seq(30,80,10), seq(30,80,10), las=1)
axis(1, at=c(1,2,3.5,4.5,6,7), labels=c("0:1","0:2","JE","GO","0:1","0:2"), las=1)
box()

#add overall mean for JE and GO
#model predicted means

lsm<-emmeans(MCc,~Population+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Population+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(3.5,4.5)
points(position,plot$response,pch=4,cex=2)

#add mean for JE incompetition with GO
M2Cc<-lm(log(Root_biomass_fraction_.) ~  Mixture_Code1 +Initial_biomass,data=dat_2r)












#MIXTURE SERIES COMBINED HARVEST 1

dat_1<-read.table(file="Series_Combined.txt", header=TRUE,sep="\t",dec = ",")

dat_1<-read.table(file="Combined1.txt", header=TRUE,sep="\t",dec = ".")# add fert effect

str(dat_1)
#DRY BIOMASS
#Use lmer with RMEL=F
M1 <- lmer(Biomass_DW ~ Population+(1 | Mixture_series)+(1 | Incompetition)+Initial_biomass, REML=F, data = dat_1)
M1a <- lmer(sqrt(Biomass_DW) ~ Population+(1 | Mixture_series)+(1 | Incompetition)+Initial_biomass, REML=F, data = dat_1)
M1b <- lmer(log(Biomass_DW) ~ Population+(1 | Mixture_series)+(1 | Incompetition)+Initial_biomass, REML=F, data = dat_1)

#Check themodel(DHARMarequired)
sim1a<-simulateResiduals(fittedModel= M1,n=1000)
plot(sim1a, quantreg= FALSE)

sim1a<-simulateResiduals(fittedModel= M1a,n=1000) #fit better
plot(sim1a, quantreg= FALSE)

sim1a<-simulateResiduals(fittedModel= M1b,n=1000)
plot(sim1a, quantreg= FALSE)

summary(M1a)
#sqrt transformation is better
#Significant test
Anova(M1a,type="III",test="Chisq")


# use lm for JE incompetition GO  
str(dat_2)
dat_2<-read.table(file="JE.txt", header=TRUE,sep="\t",dec = ".")
Mixture_Code1<-as.factor(dat_2$Mixture_Code)

## LEAF COUNT
str(dat_1)
#Use lmer with RMEL=F
M1 <- lmer(Leaf_Count ~ Population+(1 | Mixture_series)+(1 | Incompetition)+Initial_biomass, REML=F, data = dat_1)
M1a <- lmer(sqrt(Leaf_Count) ~ Population+(1 | Mixture_series)+(1 | Incompetition)+Initial_biomass, REML=F, data = dat_1)
M1b <- lmer(log(Leaf_Count) ~ Population+(1 | Fertilization)+(1 | Mixture_series)+(1 | Incompetition)+Initial_biomass, REML=F, data = dat_1)

#check model with Rf
M1b <- lmer(log(Leaf_Count) ~ Population+(1 | Fertilization)+(1 | Mixture_series)+(1 | Incompetition)+Initial_biomass, REML=F, data = dat_1)
M1ba <- lmer(log(Leaf_Count) ~ Population+(1 | Fertilization)+(1 | Mixture_series)+Initial_biomass, REML=F, data = dat_1)
M1bc <- lmer(log(Leaf_Count) ~ Population+(1 | Fertilization)+Initial_biomass, REML=F, data = dat_1)
M1bd <- lmer(log(Leaf_Count) ~ Population+(1 | Mixture_series)+Initial_biomass, REML=F, data = dat_1)

anova(M1b,M1ba, M1bc, M1bd)
#M1ba is better

#Check themodel(DHARMarequired)
sim1a<-simulateResiduals(fittedModel= M1,n=1000)
plot(sim1a, quantreg= FALSE)

sim1a<-simulateResiduals(fittedModel= M1a,n=1000)
plot(sim1a, quantreg= FALSE)

sim1a<-simulateResiduals(fittedModel= M1b,n=1000)#fit better
plot(sim1a, quantreg= FALSE)

summary(M1b)
#log transformation is better
#Significant test
Anova(M1b,type="III",test="Chisq")

#multiple comparaison
library(multcomp)
gltt_M1a<-glht(M1a,linfct=mcp(Population="Tukey"))
gltt_M1a
summary(gltt_M1a)
post<-cld(gltt_M1a, level = 0.05, decreasing = F)
post
#

# use lm for JE incompetition GO  
str(dat_2)
dat_2<-read.table(file="JE.txt", header=TRUE,sep="\t",dec = ".")


M2a<-lm(JE_Leaf_Count ~  Mixture_series +Initial_biomass,data=dat_2)

# Test homogeneity of variance
leveneTest(JE_Leaf_Count~Mixture_Code,data=dat_2)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.17


#Check normality with shapiro test
shapiro.test(dat_2$JE_Leaf_Count) 

#P VALUE > 0.05 the assumption is valid, p =0.0008

#transformation
M2b<-lm(sqrt(JE_Leaf_Count) ~  Mixture_series +Initial_biomass,data=dat_2)
M2c<-lm(log(JE_Leaf_Count) ~  Mixture_series +Initial_biomass,data=dat_2)

plot(M2a)  
hist(resid(M2a))

plot(M2b)  
hist(resid(M2b)) 

plot(M2c) # 
hist(resid(M2c)) # log fit better
summary(m4)

#log transformation fit better
Anova(M2c,type="III",test="F")
##
DunnettTest(log(JE_Leaf_Count) ~  Mixture_series,data=dat_2)
#multiple comparaison
Mixture_Code1<-as.factor(dat_2$Mixture_Code)
str(dat_2)
M2c<-lm(log(JE_Leaf_Count) ~  Mixture_Code1 +Initial_biomass,data=dat_2)

gltt_M2c<-glht(M2c,linfct=mcp(Mixture_Code1="Tukey"))
gltt_M2c
summary(gltt_M2c)
post<-cld(gltt_M2c, level = 0.05, decreasing = F)
post
#GO incompetition JE
dat_3<-read.table(file="GO.txt", header=TRUE,sep="\t",dec = ".")
Mixture_Code1<-as.factor(dat_3$Mixture_Code)
str(dat_3)
M3<-lm(GO_Leaf_Count ~  Mixture_Code1 +Initial_biomass,data=dat_3)
M3a<-lm(sqrt(GO_Leaf_Count) ~  Mixture_Code1 +Initial_biomass,data=dat_3)
M3b<-lm(log(GO_Leaf_Count) ~  Mixture_Code1 +Initial_biomass,data=dat_3)
# Test homogeneity of variance
leveneTest(GO_Leaf_Count~Mixture_Code1,data=dat_3)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.4

#Check normality with shapiro test
shapiro.test(dat_3$GO_Leaf_Count) 

#P VALUE > 0.05 the assumption is valid, p =0.008


plot(M3) # 
hist(resid(M3))
plot(M3a) # 
hist(resid(M3a)) 
plot(M3b) # 
hist(resid(M3b)) #log fit better
summary(M3)

#log transformed data fit better
Anova(M3b,type="III",test="F")
##
DunnettTest(log(GO_Leaf_Count)~  Mixture_series,data=dat_3)
#multiple comparaison
Mixture_Code1<-as.factor(dat_3$Mixture_Code)
str(dat_3)
M3b<-lm(log(GO_Leaf_Count) ~  Mixture_Code1 +Initial_biomass,data=dat_3)

gltt_M3b<-glht(M3b,linfct=mcp(Mixture_Code1="Tukey"))
gltt_M3b
summary(gltt_M3b)
post<-cld(gltt_M3b, level = 0.05, decreasing = F)
post
#Data for Box_plot 
dat_4<-read.table(file="Data_Box.txt", header=TRUE,sep="\t",dec = ",")
str(dat_4)

boxplot(Leaf_Count ~Mixture_Code,data=dat_4)


position<-c(1,2,3,4,5,6.5,7.5,9,10,11,12,13)
boxplot(Leaf_Count ~Mixture_Code,data=dat_4,axes=F,at=position ,border="gray40",xlab= "Population", ylab="Total number of leaves",las=2, 
        col=c("#0033FF","#0033FF","#0033FF","#0033FF","#0033FF","#0033FF","#FFFF00","#FFFF00","#FFFF00","#FFFF00","#FFFF00","#FFFF00"))
axis(2,at=seq(5,14,2), seq(5,14,2), las=1)
axis(1, at=c(1,2,3,4,5,6.5,7.5,9,10,11,12,13), labels=c("0:1","1:1","2:1","1:2","0:2","JE","GO","0:1","1:1","2:1","1:2","0:2"), las=1)
box()

#add overall mean for JE and GO
#model predicted means

lsm<-emmeans(M1b,~Population+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Population+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(6.5,7.5)
points(position,plot$response,pch=4,cex=2)

#add mean for JE incompetition with GO
M2c<-lm(log(JE_Leaf_Count) ~  Mixture_Code1 +Initial_biomass,data=dat_2)

lsm<-emmeans(M2c,~Mixture_Code1+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Mixture_Code1+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(1,2,3,4,5)
points(position,plot$response,pch=4,cex=2)
#means<-tapply(dat_2$JE_Biomass,dat_2$Mixture_Code,mean)
#means
#points(means,col="black",pch=4, at=position)

#add mean for GO incompetition with JE
#add mean for JE incompetition with GO
M3b<-lm(log(GO_Leaf_Count) ~  Mixture_Code1 +Initial_biomass,data=dat_3)

lsm<-emmeans(M3b,~Mixture_Code1+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Mixture_Code1+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(9,10,11,12,13)
points(position,plot$response,pch=4,cex=2)

## LENGTH OF LONGEST LEAF
str(dat_1)
#Use lmer with RMEL=F
M1 <- lmer(LLL ~ Population+(1 | Mixture_series)+(1 | Incompetition)+Initial_biomass, REML=F, data = dat_1)
M1a <- lmer(sqrt(LLL) ~ Population+(1 | Mixture_series)+(1 | Incompetition)+Initial_biomass, REML=F, data = dat_1)
M1b <- lmer(log(LLL) ~ Population+(1 | Mixture_series)+(1 | Incompetition)+Initial_biomass, REML=F, data = dat_1)

#Check themodel(DHARMarequired)
sim1a<-simulateResiduals(fittedModel= M1,n=1000) #fit better
plot(sim1a, quantreg= FALSE)

sim1a<-simulateResiduals(fittedModel= M1a,n=1000)
plot(sim1a, quantreg= FALSE)

sim1a<-simulateResiduals(fittedModel= M1b,n=1000)
plot(sim1a, quantreg= FALSE)

summary(M1a)
#untransformed is better
#Significant test
Anova(M1,type="III",test="Chisq")

#multiple comparaison
library(multcomp)
gltt_M1<-glht(M1,linfct=mcp(Population="Tukey"))
gltt_M1
summary(gltt_M1)
post<-cld(gltt_M1, level = 0.05, decreasing = F)
post
#

# use lm for JE incompetition GO  
str(dat_2)
dat_2<-read.table(file="JE.txt", header=TRUE,sep="\t",dec = ".")

Mixture_Code1<-as.factor(dat_2$Mixture_series)
# lm
M2a<-lm(JE_LLL ~  Mixture_Code1 +Initial_biomass,data=dat_2)

# Test homogeneity of variance
leveneTest(JE_LLL~Mixture_Code1,data=dat_2)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.8


#Check normality with shapiro test
shapiro.test(dat_2$JE_LLL) 

#P VALUE > 0.05 the assumption is valid, p =0.7

#transformation
M2b<-lm(sqrt(JE_LLL) ~  Mixture_Code1 +Initial_biomass,data=dat_2)
M2c<-lm(log(JE_LLL) ~  Mixture_Code1 +Initial_biomass,data=dat_2)

plot(M2a)  
hist(resid(M2a)) #fit better

plot(M2b)  
hist(resid(M2b)) 

plot(M2c) # 
hist(resid(M2c)) 
summary(m4)

#untransformation fit better
Anova(M2a,type="III",test="F")
##
DunnettTest(log(JE_LLL) ~  Mixture_series,data=dat_2)
#multiple comparaison
Mixture_Code1<-as.factor(dat_2$Mixture_Code)
str(dat_2)
M2a<-lm(JE_LLL ~  Mixture_Code1 +Initial_biomass,data=dat_2)

gltt_M2a<-glht(M2a,linfct=mcp(Mixture_Code1="Tukey"))
gltt_M2a
summary(gltt_M2a)
post<-cld(gltt_M2a, level = 0.05, decreasing = F)
post
#GO incompetition JE
dat_3<-read.table(file="GO.txt", header=TRUE,sep="\t",dec = ".")
Mixture_Code1<-as.factor(dat_3$Mixture_series)
str(dat_3)
M3<-lm(GO_LLL ~  Mixture_Code1 +Initial_biomass,data=dat_3)
M3a<-lm(sqrt(GO_LLL) ~  Mixture_Code1 +Initial_biomass,data=dat_3)
M3b<-lm(log(GO_LLL) ~  Mixture_Code1 +Initial_biomass,data=dat_3)
# Test homogeneity of variance
leveneTest(GO_Leaf_Count~Mixture_Code1,data=dat_3)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.4

#Check normality with shapiro test
shapiro.test(dat_3$GO_LLL) 

#P VALUE > 0.05 the assumption is valid, p =0.02


plot(M3) # 
hist(resid(M3)) 
plot(M3a) # 
hist(resid(M3a)) #fit better
plot(M3b) # 
hist(resid(M3b)) 
summary(M3)

#sqrt transformed data fit better
Anova(M3a,type="III",test="F")
##
DunnettTest(sqrt(GO_LLL)~  Mixture_series,data=dat_3)
#multiple comparaison
Mixture_Code1<-as.factor(dat_3$Mixture_Code)
str(dat_3)
M3a<-lm(sqrt(GO_LLL) ~  Mixture_Code1 +Initial_biomass,data=dat_3)

gltt_M3a<-glht(M3a,linfct=mcp(Mixture_Code1="Tukey"))
gltt_M3a
summary(gltt_M3a)
post<-cld(gltt_M3a, level = 0.05, decreasing = F)
post
#Data for Box_plot 
dat_4<-read.table(file="Data_Box.txt", header=TRUE,sep="\t",dec = ",")
str(dat_4)

boxplot(LLL ~Mixture_Code,data=dat_4)


position<-c(1,2,3,4,5,6.5,7.5,9,10,11,12,13)
boxplot(LLL ~Mixture_Code,data=dat_4,axes=F,at=position ,border="gray40",xlab= "Population", ylab="Length of the longest leaf [cm]",las=2, 
        col=c("#0033FF","#0033FF","#0033FF","#0033FF","#0033FF","#0033FF","#FFFF00","#FFFF00","#FFFF00","#FFFF00","#FFFF00","#FFFF00"))
axis(2,at=seq(15,40,5), seq(15,40,5), las=1)
axis(1, at=c(1,2,3,4,5,6.5,7.5,9,10,11,12,13), labels=c("0:1","1:1","2:1","1:2","0:2","JE","GO","0:1","1:1","2:1","1:2","0:2"), las=1)
box()

#add overall mean for JE and GO
#model predicted means

lsm<-emmeans(M1,~Population+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Population+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(6.5,7.5)
points(position,plot$emmean,pch=4,cex=2)

#add mean for JE incompetition with GO
M2a<-lm(JE_LLL ~  Mixture_Code1 +Initial_biomass,data=dat_2)

lsm<-emmeans(M2a,~Mixture_Code1+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Mixture_Code1+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(1,2,3,4,5)
points(position,plot$emmean,pch=4,cex=2)

#add mean for GO incompetition with JE
M3a<-lm(sqrt(GO_LLL) ~  Mixture_Code1 +Initial_biomass,data=dat_3)

lsm<-emmeans(M3a,~Mixture_Code1+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Mixture_Code1+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(9,10,11,12,13)
points(position,plot$response,pch=4,cex=2)

### LEAF WATER CONTENT
str(dat_1)
#Use lmer with RMEL=F
M1 <- lmer(LWC ~ Population+(1 | Mixture_series)+(1 | Incompetition)+Initial_biomass, REML=F, data = dat_1)
M1a <- lmer(sqrt(LWC) ~ Population+(1 | Mixture_series)+(1 | Incompetition)+Initial_biomass, REML=F, data = dat_1)
M1b <- lmer(log(LWC) ~ Population+(1 | Mixture_series)+(1 | Incompetition)+Initial_biomass, REML=F, data = dat_1)

#Check themodel(DHARMarequired)
sim1a<-simulateResiduals(fittedModel= M1,n=1000) #fit better
plot(sim1a, quantreg= FALSE)

sim1a<-simulateResiduals(fittedModel= M1a,n=1000)
plot(sim1a, quantreg= FALSE)

sim1a<-simulateResiduals(fittedModel= M1b,n=1000)
plot(sim1a, quantreg= FALSE)

summary(M1a)
#untransformed is better
#Significant test
Anova(M1,type="III",test="Chisq")

#multiple comparaison
library(multcomp)
gltt_M1<-glht(M1,linfct=mcp(Population="Tukey"))
gltt_M1
summary(gltt_M1)
post<-cld(gltt_M1, level = 0.05, decreasing = F)
post
#

# use lm for JE incompetition GO  
str(dat_2)
dat_2<-read.table(file="JE.txt", header=TRUE,sep="\t",dec = ".")

Mixture_Code1<-as.factor(dat_2$Mixture_series)
str(dat_2)
# lm
M2a<-lm(JE_LWC ~  Mixture_Code1 +Initial_biomass,data=dat_2)

# Test homogeneity of variance
leveneTest(JE_LWC~Mixture_Code1,data=dat_2)
#Assumption is valid if pvalue >.05
# assumption is not valid, p =0.03 


#Check normality with shapiro test
shapiro.test(dat_2$JE_LWC) 

#P VALUE > 0.05 the assumption is valid, p =0.0003

#transformation
M2b<-lm(sqrt(JE_LWC) ~  Mixture_Code1 +Initial_biomass,data=dat_2)
M2c<-lm(log(JE_LWC) ~  Mixture_Code1 +Initial_biomass,data=dat_2)
#LM
plot(M2a)  
hist(resid(M2a)) 

plot(M2b)  
hist(resid(M2b)) 

plot(M2c) # 
hist(resid(M2c)) 

summary(m4)

#SQRTtransformation fit better
Anova(M2b,type="III",test="F")
##
DunnettTest(sqrt(JE_LWC) ~  Mixture_Code1,data=dat_2)
#multiple comparaison
Mixture_Code1<-as.factor(dat_2$Mixture_Code)
str(dat_2)
M2b<-lm(sqrt(JE_LWC )~  Mixture_Code1 +Initial_biomass,data=dat_2)

gltt_M2b<-glht(M2b,linfct=mcp(Mixture_Code1="Tukey"))
gltt_M2b
summary(gltt_M2b)
post<-cld(gltt_M2b, level = 0.05, decreasing = F)
post
#GO incompetition JE
dat_3<-read.table(file="GO.txt", header=TRUE,sep="\t",dec = ".")
Mixture_Code1<-as.factor(dat_3$Mixture_series)
str(dat_3)
M3<-lm(GO_LWC ~  Mixture_Code1 +Initial_biomass,data=dat_3)
M3a<-lm(sqrt(GO_LWC) ~  Mixture_Code1 +Initial_biomass,data=dat_3)
M3b<-lm(log(GO_LWC) ~  Mixture_Code1 +Initial_biomass,data=dat_3)
# Test homogeneity of variance
leveneTest(GO_LWC~Mixture_Code1,data=dat_3)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.13

#Check normality with shapiro test
shapiro.test(dat_3$GO_LWC) 

#P VALUE > 0.05 the assumption is valid, p =0.78


plot(M3) # 
hist(resid(M3)) # fit better
plot(M3a) # 
hist(resid(M3a)) 
plot(M3b) # 
hist(resid(M3b)) 
summary(M3)

#untransformed data fit better
Anova(M3,type="III",test="F")
##
DunnettTest(GO_LWC~  Mixture_series,data=dat_3)
#multiple comparaison
Mixture_Code1<-as.factor(dat_3$Mixture_Code)
str(dat_3)
M3<-lm(GO_LWC ~  Mixture_Code1 +Initial_biomass,data=dat_3)

gltt_M3<-glht(M3,linfct=mcp(Mixture_Code1="Tukey"))
gltt_M3
summary(gltt_M3)
post<-cld(gltt_M3, level = 0.05, decreasing = F)
post
#Data for Box_plot 
dat_4<-read.table(file="Data_Box.txt", header=TRUE,sep="\t",dec = ",")
str(dat_4)

boxplot(LWC ~Mixture_Code,data=dat_4)


position<-c(1,2,3,4,5,6.5,7.5,9,10,11,12,13)
boxplot(LWC ~Mixture_Code,data=dat_4,axes=F,at=position ,border="gray40",xlab= "Population", ylab="Leaf water content [%]",las=2, 
        col=c("#0033FF","#0033FF","#0033FF","#0033FF","#0033FF","#0033FF","#FFFF00","#FFFF00","#FFFF00","#FFFF00","#FFFF00","#FFFF00"))
axis(2,at=seq(82,92,2), seq(82,92,2), las=1)
axis(1, at=c(1,2,3,4,5,6.5,7.5,9,10,11,12,13), labels=c("0:1","1:1","2:1","1:2","0:2","JE","GO","0:1","1:1","2:1","1:2","0:2"), las=1)
box()

#add overall mean for JE and GO
#model predicted means

lsm<-emmeans(M1,~Population+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Population+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(6.5,7.5)
points(position,plot$emmean,pch=4,cex=2)

#add mean for JE incompetition with GO
M2b<-lm(sqrt(JE_LWC) ~  Mixture_Code1 +Initial_biomass,data=dat_2)

lsm<-emmeans(M2b,~Mixture_Code1+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Mixture_Code1+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(1,2,3,4,5)
points(position,plot$response,pch=4,cex=2)

#add mean for GO incompetition with JE
M3<-lm(GO_LWC ~  Mixture_Code1 +Initial_biomass,data=dat_3)

lsm<-emmeans(M3,~Mixture_Code1+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Mixture_Code1+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(9,10,11,12,13)
points(position,plot$emmean,pch=4,cex=2)

## LEAF MASS PER AREA (LMA)
dat_1h1<-read.table(file="Combine_LMA.txt", header=TRUE,sep="\t",dec = ".")

str(dat_1h1)
#Use lmer with RMEL=F
M1 <- lmer(LMA  ~ Population+(1 | Mixture_series)+(1 | Incompetition)+Initial_biomass, REML=F, data = dat_1a)
M1a <- lmer(sqrt(LMA ) ~ Population+(1 | Mixture_series)+(1 | Incompetition)+Initial_biomass, REML=F, data = dat_1a)
M1b <- lmer(log(LMA ) ~ Population+(1 | Mixture_series)+(1 | Incompetition)+Initial_biomass, REML=F, data = dat_1a)

#Check themodel(DHARMarequired)
sim1a<-simulateResiduals(fittedModel= M1,n=1000) #fit better
plot(sim1a, quantreg= FALSE)

sim1a<-simulateResiduals(fittedModel= M1a,n=1000)
plot(sim1a, quantreg= FALSE)

sim1a<-simulateResiduals(fittedModel= M1b,n=1000)
plot(sim1a, quantreg= FALSE)

summary(M1a)
#no much difference btwn untransformed and transformed
#Significant test
Anova(M1,type="III",test="Chisq")

#multiple comparaison
library(multcomp)
gltt_M1<-glht(M1,linfct=mcp(Population="Tukey"))
gltt_M1
summary(gltt_M1)
post<-cld(gltt_M1, level = 0.05, decreasing = F)
post
#

# use lm for JE incompetition GO  
str(dat_2h1)
dat_2h1<-read.table(file="JE_LMA.txt", header=TRUE,sep="\t",dec = ".")
Mixture_Code1<-as.factor(dat_2h1$Mixture_Code)
str(dat_2h1)
# lm
M2a<-lm(LMA ~  Mixture_Code1 +Initial_biomass,data=dat_2h1)
M2b<-lm(sqrt(LMA) ~  Mixture_Code1 +Initial_biomass,data=dat_2h1)
M2c<-lm(log(LMA) ~  Mixture_Code1 +Initial_biomass,data=dat_2h1)

# Test homogeneity of variance
leveneTest(LMA~Mixture_Code1,data=dat_2h1)
#Assumption is valid if pvalue >.05
# assumption is not valid, p =0.03 


#Check normality with shapiro test
shapiro.test(dat_2h1$LMA) 

#P VALUE > 0.05 the assumption is valid, p =0.47

#LM
plot(M2a)  
hist(resid(M2a)) 
plot(M2b)  
hist(resid(M2b))
plot(M2c)  
hist(resid(M2c))
summary(m4)

# sqrt ransformed fit better
Anova(M2b,type="III",test="F")
##
DunnettTest(LMA ~  Mixture_series,data=dat_2h1)
#multiple comparaison
Mixture_Code1<-as.factor(dat_2h1$Mixture_Code)
str(dat_2h1)
M2a<-lm(LMA ~  Mixture_Code1 +Initial_biomass,data=dat_2h1)

gltt_M2a<-glht(M2a,linfct=mcp(Mixture_Code1="Tukey"))
gltt_M2a
summary(gltt_M2a)
post<-cld(gltt_M2a, level = 0.05, decreasing = F)
post
#GO incompetition JE
dat_3h1<-read.table(file="GO_LMA.txt", header=TRUE,sep="\t",dec = ".")
str(dat_3h1)
Mixture_Code1<-as.factor(dat_3h1$Mixture_Code) # Mixture_Code1= Mixture_series
str(dat_3h1)
M3<-lm(LMA  ~  Mixture_Code1 +Initial_biomass,data=dat_3h1)
M3a<-lm(sqrt(LMA ) ~  Mixture_Code1 +Initial_biomass,data=dat_3h1)
M3b<-lm(log(LMA) ~  Mixture_Code1 +Initial_biomass,data=dat_3h1)
# Test homogeneity of variance
leveneTest(LMA ~Mixture_Code,data=dat_3h1)
#Assumption is valid if pvalue >.05
# assumption is not valid, p =0.03

#Check normality with shapiro test
shapiro.test(dat_3h1$LMA ) 

#P VALUE > 0.05 the assumption is valid, p =0.04


plot(M3) # 
hist(resid(M3)) #
plot(M3a) # 
hist(resid(M3a)) #
plot(M3b) # 
hist(resid(M3b)) #
summary(M3)

#no difference, untransformed data fit better
Anova(M3,type="III",test="F")
##
DunnettTest(LMA ~  Mixture_series,data=dat_3h1)
#multiple comparaison
Mixture_Code1<-as.factor(dat_3h1$Mixture_Code)
str(dat_3)
M3<-lm(LMA  ~  Mixture_Code1 +Initial_biomass,data=dat_3h1)

gltt_M3b<-glht(M3,linfct=mcp(Mixture_Code1="Tukey"))
gltt_M3b
summary(gltt_M3b)
post<-cld(gltt_M3b, level = 0.05, decreasing = F)
post
#Data for Box_plot 
dat_4<-read.table(file="Box_LMA.txt", header=TRUE,sep="\t",dec = ".")
str(dat_4)

boxplot(LMA~Mixture_Code,data=dat_4)


position<-c(1,2,3,4.5,5.5,7,8,9)
boxplot(LMA ~Mixture_Code,data=dat_4,axes=F,at=position ,border="gray40",xlab= "Population", ylab="Leaf mass per area [mg mm-2]",las=2, 
        col=c("#0033FF","#0033FF","#0033FF","#0033FF","#FFFF00","#FFFF00","#FFFF00","#FFFF00"))
axis(2,at=seq(0.020,0.040,0.005), seq(0.020,0.040,0.005), las=1)
axis(1, at=c(1,2,3,4.5,5.5,7,8,9), labels=c("0:1","1:1","0:2","JE","GO","0:1","1:1","0:2"), las=1)
box()

#add overall mean for JE and GO
#model predicted means

lsm<-emmeans(M1,~Population+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Population+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(4.5,5.5)
points(position,plot$emmean,pch=4,cex=2)

#add mean for JE incompetition with GO
M2a<-lm(LMA ~  Mixture_Code1 +Initial_biomass,data=dat_2a)

lsm<-emmeans(M2a,~Mixture_Code1+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Mixture_Code1+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(1,2,3)
points(position,plot$emmean,pch=4,cex=2)

#add mean for GO incompetition with JE
M3a<-lm(LMA ~  Mixture_Code1 +Initial_biomass,data=dat_3)

lsm<-emmeans(M3a,~Mixture_Code1+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Mixture_Code1+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(7,8,9)
points(position,plot$emmean,pch=4,cex=2)
### Specific Leaf Area (SLA)
str(dat_1h1)
#Use lmer with RMEL=F
M1 <- lmer(JE_SLA  ~ Population+(1 | Mixture_series)+(1 | Incompetition)+Initial_biomass, REML=F, data = dat_1h1)
M1a <- lmer(sqrt(JE_SLA ) ~ Population+(1 | Mixture_series)+(1 | Incompetition)+Initial_biomass, REML=F, data = dat_1h1)
M1b <- lmer(log(JE_SLA ) ~ Population+(1 | Mixture_series)+(1 | Incompetition)+Initial_biomass, REML=F, data = dat_1h1)

#Check themodel(DHARMarequired)
sim1a<-simulateResiduals(fittedModel= M1,n=1000) #fit better
plot(sim1a, quantreg= FALSE)

sim1a<-simulateResiduals(fittedModel= M1a,n=1000)
plot(sim1a, quantreg= FALSE)

sim1a<-simulateResiduals(fittedModel= M1b,n=1000)
plot(sim1a, quantreg= FALSE)

summary(M1)
#untransformed is better
#Significant test
Anova(M1,type="III",test="Chisq")

#multiple comparaison
library(multcomp)
gltt_M1<-glht(M1,linfct=mcp(Population="Tukey"))
gltt_M1
summary(gltt_M1)
post<-cld(gltt_M1, level = 0.05, decreasing = F)
post
#

# use lm for JE incompetition GO  
str(dat_2a)
dat_2a<-read.table(file="Competitive_JE.txt", header=TRUE,sep="\t",dec = ".")


# lm
M2a<-lm(Competitive_ability_JE ~  Mixture_series +Initial_biomass,data=dat_2a)

# Test homogeneity of variance
leveneTest(Competitive_ability_JE~Mixture_Code,data=dat_2a)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.17 


#Check normality with shapiro test
shapiro.test(dat_2a$Competitive_ability_JE) 

#P VALUE > 0.05 the assumption is valid, p =0.6

#LM
plot(M2a)  
hist(resid(M2a)) 

summary(m4)

#untransformed fit better
Anova(M2a,type="III",test="F")
##
DunnettTest(Competitive_ability_JE ~  Mixture_series,data=dat_2a)
#multiple comparaison
Mixture_Code1<-as.factor(dat_2h1$Mixture_Code)
str(dat_2h1)
M2a<-lm(JE_SLA ~  Mixture_Code1 +Initial_biomass,data=dat_2h1)
M2b<-lm(JE_SLA ~  Mixture_Code1 +Initial_biomass,data=dat_2h1)
M2c<-lm(JE_SLA ~  Mixture_Code1 +Initial_biomass,data=dat_2h1)
# Test homogeneity of variance
leveneTest(JE_SLA ~Mixture_Code1,data=dat_2h1)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.26

#Check normality with shapiro test
shapiro.test(dat_2h1$JE_SLA ) 

#P VALUE > 0.05 the assumption is valid, p =0.45


plot(M2a) # 
hist(resid(M2a)) # fit better
plot(M2b) # 
hist(resid(M2b)) 
plot(M2b) # 
hist(resid(M2c)) 
summary(M2c)

#untransformed data fit better
Anova(M2a,type="III",test="F")
##
DunnettTest(JE_SLA ~  Mixture_series,data=dat_2h1)
gltt_M2a<-glht(M2a,linfct=mcp(Mixture_Code1="Tukey"))
gltt_M2a
summary(gltt_M2a)
post<-cld(gltt_M2a, level = 0.05, decreasing = F)
post
#GO incompetition JE

dat_3<-read.table(file="GO.txt", header=TRUE,sep="\t",dec = ".")
Mixture_Code1<-as.factor(dat_3h1$Mixture_Code)
str(dat_3)
M3<-lm(GO_SLA  ~  Mixture_Code1 +Initial_biomass,data=dat_3h1)
M3a<-lm(sqrt(GO_SLA ) ~  Mixture_Code1 +Initial_biomass,data=dat_3h1)
M3b<-lm(log(GO_SLA ) ~  Mixture_Code1 +Initial_biomass,data=dat_3h1)
# Test homogeneity of variance
leveneTest(GO_SLA ~Mixture_Code1,data=dat_3h1)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.1

#Check normality with shapiro test
shapiro.test(dat_3h1$GO_SLA ) 

#P VALUE > 0.05 the assumption is valid, p =0.2


plot(M3) # 
hist(resid(M3)) #
plot(M3a) # 
hist(resid(M3a)) 
plot(M3b) # 
hist(resid(M3b)) #fit better
summary(M3)

#transformed data fit better
Anova(M3b,type="III",test="F")
##
DunnettTest(GO_SLA ~  Mixture_Code1,data=dat_3h1)# not needed
#multiple comparaison
Mixture_Code1<-as.factor(dat_3$Mixture_Code)
str(dat_3)
M3b<-lm(log(GO_SLA)  ~  Mixture_Code1 +Initial_biomass,data=dat_3)

gltt_M3b<-glht(M3b,linfct=mcp(Mixture_Code1="Tukey"))
gltt_M3b
summary(gltt_M3b)
post<-cld(gltt_M3b, level = 0.05, decreasing = F)
post
#Data for Box_plot 
dat_4h1<-read.table(file="Box_LMA.txt", header=TRUE,sep="\t",dec = ".")
str(dat_4h1)
Mixture_Code1<-as.factor(dat_3h1$Mixture_Code) # Mixture_Code = Mixture_series
boxplot(JE_SLA ~Mixture_Code,data=dat_4h1)


position<-c(1,2,3,4.5,5.5,7,8,9)
boxplot(JE_SLA ~Mixture_Code,data=dat_4h1,axes=F,at=position ,border="gray40",xlab= "Population", ylab="Specific leaf area [mm2mg-1]",las=2, 
        col=c("#0033FF","#0033FF","#0033FF","#0033FF","#FFFF00","#FFFF00","#FFFF00","#FFFF00"))
axis(2,at=seq(25,55,5), seq(25,55,5), las=1)
axis(1, at=c(1,2,3,4.5,5.5,7,8,9), labels=c("0:1","1:1","0:2","JE","GO","0:1","1:1","0:2"), las=1)
box()

#add overall mean for JE and GO
#model predicted means

lsm<-emmeans(M1,~Population+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Population+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(4.5,5.5)
points(position,plot$emmean,pch=4,cex=2)

#add mean for JE incompetition with GO
M2a<-lm(JE_SLA ~  Mixture_Code1 +Initial_biomass,data=dat_2)

lsm<-emmeans(M2a,~Mixture_Code1+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Mixture_Code1+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(1,2,3)
points(position,plot$emmean,pch=4,cex=2)

#add mean for GO incompetition with JE
M3b<-lm(log(GO_SLA) ~  Mixture_Code1 +Initial_biomass,data=dat_3)

lsm<-emmeans(M3b,~Mixture_Code1+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Mixture_Code1+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(7,8,9)
points(position,plot$response,pch=4,cex=2)
##
## NITROGEN TO CARBON CONTENT
dat_1<-read.table(file="Combined_CN.txt", header=TRUE,sep="\t",dec = ".")

str(dat_1)
#NITROGEN
M1 <- lmer(N_. ~ Population+(1 | Mixture_series)+(1 | Incompetition)+Initial_biomass, REML=F, data = dat_1)
M1a <- lmer(sqrt(N_.) ~ Population+(1 | Mixture_series)+(1 | Incompetition)+Initial_biomass, REML=F, data = dat_1)
M1b <- lmer(log(N_.) ~ Population+(1 | Mixture_series)+(1 | Incompetition)+Initial_biomass, REML=F, data = dat_1)

#Check themodel(DHARMarequired)
sim1a<-simulateResiduals(fittedModel= M1,n=1000) #fit better
plot(sim1a, quantreg= FALSE)

sim1a<-simulateResiduals(fittedModel= M1a,n=1000)
plot(sim1a, quantreg= FALSE)

sim1a<-simulateResiduals(fittedModel= M1b,n=1000)
plot(sim1a, quantreg= FALSE)

summary(M1a)
#LOG transformed is better
#Significant test
Anova(M1b,type="III",test="Chisq")

#multiple comparaison
library(multcomp)
gltt_M1<-glht(M1b,linfct=mcp(Population="Tukey"))
gltt_M1
summary(gltt_M1)
post<-cld(gltt_M1, level = 0.05, decreasing = F)
post
#

# use lm for JE incompetition GO  
str(dat_2)
dat_2<-read.table(file="CN_JE.txt", header=TRUE,sep="\t",dec = ".")

Mixture_Code1<-as.factor(dat_2$Mixture_series)
str(dat_2)
# lm
M2a<-lm(N_. ~  Mixture_Code1 +Initial_biomass,data=dat_2)

# Test homogeneity of variance
leveneTest(N_.~Mixture_Code1,data=dat_2)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.67 


#Check normality with shapiro test
shapiro.test(dat_2$N_.) 

#P VALUE > 0.05 the assumption is not valid, p =0.0001

#transformation
M2b<-lm(sqrt(N_.) ~  Mixture_Code1 +Initial_biomass,data=dat_2)
M2c<-lm(log(N_.) ~  Mixture_Code1 +Initial_biomass,data=dat_2)
#LM
plot(M2a)  
hist(resid(M2a)) 

plot(M2b)  
hist(resid(M2b)) 

plot(M2c) # 
hist(resid(M2c)) 

summary(m4)

#log transformation fit better
Anova(M2c,type="III",test="F")
##
DunnettTest(log(N_.) ~  Mixture_Code1,data=dat_2)
#multiple comparaison
Mixture_Code1<-as.factor(dat_2$Mixture_Code)
str(dat_2)
M2c<-lm(log(N_. )~  Mixture_Code1 +Initial_biomass,data=dat_2)

gltt_M2b<-glht(M2c,linfct=mcp(Mixture_Code1="Tukey"))
gltt_M2b
summary(gltt_M2b)
post<-cld(gltt_M2b, level = 0.05, decreasing = F)
post
#GO incompetition JE
dat_3<-read.table(file="CN_GO.txt", header=TRUE,sep="\t",dec = ".")
Mixture_Code1<-as.factor(dat_3$Mixture_Code)
str(dat_3)
M3<-lm(N_.  ~  Mixture_Code1 +Initial_biomass,data=dat_3)
M3a<-lm(sqrt(N_. ) ~  Mixture_Code1 +Initial_biomass,data=dat_3)
M3b<-lm(log(N_. ) ~  Mixture_Code1 +Initial_biomass,data=dat_3)
# Test homogeneity of variance
leveneTest(N_.~Mixture_Code1,data=dat_3)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.84

#Check normality with shapiro test
shapiro.test(dat_3$N_.) 

#P VALUE > 0.05 the assumption is valid, p =0.11


plot(M3) # 
hist(resid(M3)) # 
plot(M3a) # 
hist(resid(M3a)) 
plot(M3b) # 
hist(resid(M3b)) 
summary(M3)

#log transformed data fit better
Anova(M3b,type="III",test="F")
##
DunnettTest(N_.~  Mixture_series,data=dat_3)
#multiple comparaison
Mixture_Code1<-as.factor(dat_3$Mixture_Code)
str(dat_3)
M3b<-lm(log(N_.) ~  Mixture_Code1 +Initial_biomass,data=dat_3)

gltt_M3<-glht(M3b,linfct=mcp(Mixture_Code1="Tukey"))
gltt_M3
summary(gltt_M3)
post<-cld(gltt_M3, level = 0.05, decreasing = F)
post
#Data for Box_plot 
dat_4<-read.table(file="CN_Box.txt", header=TRUE,sep="\t",dec = ".")
str(dat_4)

boxplot(N_. ~Mixture_Code,data=dat_4)

position<-c(1,2,3,4.5,5.5,7,8,9)
boxplot(N_. ~Mixture_Code,data=dat_4,axes=F,at=position ,border="gray40",xlab= "Population", ylab="Nitrogen content [%]",las=2, 
        col=c("#0033FF","#0033FF","#0033FF","#0033FF","#FFFF00","#FFFF00","#FFFF00","#FFFF00"))
axis(2,at=seq(2.0,5.5,0.5), seq(2.0,5.5,0.5), las=1)
axis(1, at=c(1,2,3,4.5,5.5,7,8,9), labels=c("0:1","1:1","0:2","JE","GO","0:1","1:1","0:2"), las=1)
box()

#add overall mean for JE and GO
#model predicted means

lsm<-emmeans(M1b,~Population+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Population+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(4.5,5.5)
points(position,plot$response,pch=4,cex=2)

#add mean for JE incompetition with GO
M2b<-lm(log(N_.) ~  Mixture_Code1 +Initial_biomass,data=dat_2)

lsm<-emmeans(M2b,~Mixture_Code1+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Mixture_Code1+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(1,2,3)
points(position,plot$response,pch=4,cex=2)

#add mean for GO incompetition with JE
M3a<-lm(log(N_.) ~  Mixture_Code1 +Initial_biomass,data=dat_3)

lsm<-emmeans(M3a,~Mixture_Code1+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Mixture_Code1+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(7,8,9)
points(position,plot$response,pch=4,cex=2)


#CARBON CONTENT

dat_1<-read.table(file="Combined_CN.txt", header=TRUE,sep="\t",dec = ".")

str(dat_1)

M1 <- lmer(C_. ~ Population+(1 | Mixture_series)+(1 | Incompetition)+Initial_biomass, REML=F, data = dat_1)
M1a <- lmer(sqrt(C_.) ~ Population+(1 | Mixture_series)+(1 | Incompetition)+Initial_biomass, REML=F, data = dat_1)
M1b <- lmer(log(C_.) ~ Population+(1 | Mixture_series)+(1 | Incompetition)+Initial_biomass, REML=F, data = dat_1)
#install.packages('TMB', type = 'source')

#Check themodel(DHARMarequired)
sim1a<-simulateResiduals(fittedModel= M1,n=1000) 
plot(sim1a, quantreg= FALSE)

sim1a<-simulateResiduals(fittedModel= M1a,n=1000)
plot(sim1a, quantreg= FALSE)

sim1a<-simulateResiduals(fittedModel= M1b,n=1000)
plot(sim1a, quantreg= FALSE)

summary(M1b)
#untransformed is better
#Significant test
Anova(M1,type="III",test="Chisq")

#multiple comparaison
library(multcomp)
gltt_M1<-glht(M1b,linfct=mcp(Population="Tukey"))
gltt_M1
summary(gltt_M1)
post<-cld(gltt_M1, level = 0.05, decreasing = F)
post
#

# use lm for JE incompetition GO  

dat_2<-read.table(file="CN_JE.txt", header=TRUE,sep="\t",dec = ".")

Mixture_Code1<-as.factor(dat_2$Mixture_series)
str(dat_2)
# lm
M2a<-lm(C_. ~  Mixture_COde1 +Initial_biomass,data=dat_2)

# Test homogeneity of variance
leveneTest(C_.~Mixture_Code1,data=dat_2)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.46 


#Check normality with shapiro test
shapiro.test(dat_2$C_.) 

#P VALUE > 0.05 the assumption is not valid, p =0.0001

#transformation
M2b<-lm(sqrt(C_.) ~  Mixture_Code1 +Initial_biomass,data=dat_2)
M2c<-lm(log(C_.) ~  Mixture_Code1 +Initial_biomass,data=dat_2)
#LM
plot(M2a)  
hist(resid(M2a)) 

plot(M2b)  
hist(resid(M2b)) 

plot(M2c) # 
hist(resid(M2c)) 

summary(m4)

#log transformation is better
Anova(M2c,type="III",test="F")
##
DunnettTest(log(C_.) ~  Mixture_series,data=dat_2)
#multiple comparaison
Mixture_Code1<-as.factor(dat_2$Mixture_Code)
str(dat_2)
M2c<-lm(log(C_. )~  Mixture_Code1 +Initial_biomass,data=dat_2)

gltt_M2b<-glht(M2c,linfct=mcp(Mixture_Code1="Tukey"))
gltt_M2b
summary(gltt_M2b)
post<-cld(gltt_M2b, level = 0.05, decreasing = F)
post
#GO incompetition JE
dat_3<-read.table(file="CN_GO.txt", header=TRUE,sep="\t",dec = ".")
Mixture_Code1<-as.factor(dat_3$Mixture_series)
str(dat_3)
M3<-lm(C_.  ~  Mixture_series +Initial_biomass,data=dat_3)
M3a<-lm(sqrt(C_. ) ~  Mixture_series +Initial_biomass,data=dat_3)
M3b<-lm(log(C_. ) ~  Mixture_series +Initial_biomass,data=dat_3)
# Test homogeneity of variance
leveneTest(C_.~Mixture_Code1,data=dat_3)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.48

#Check normality with shapiro test
shapiro.test(dat_3$C_.) 

#P VALUE > 0.05 the assumption is valid, p =0.0001


plot(M3) # 
hist(resid(M3)) # 
plot(M3a) # 
hist(resid(M3a)) 
plot(M3b) # 
hist(resid(M3b)) 
summary(M3)

#log transformed data fit better
Anova(M3b,type="III",test="F")
##
DunnettTest(C_.~  Mixture_series,data=dat_3)
#multiple comparaison
Mixture_Code1<-as.factor(dat_3$Mixture_Code)
str(dat_3)
M3b<-lm(log(C_.) ~  Mixture_Code1 +Initial_biomass,data=dat_3)

gltt_M3<-glht(M3b,linfct=mcp(Mixture_Code1="Tukey"))
gltt_M3
summary(gltt_M3)
post<-cld(gltt_M3, level = 0.05, decreasing = F)
post
#Data for Box_plot 
dat_4<-read.table(file="CN_Box.txt", header=TRUE,sep="\t",dec = ".")
str(dat_4)

boxplot(C_. ~Mixture_Code,data=dat_4)

position<-c(1,2,3,4.5,5.5,7,8,9)
boxplot(C_. ~Mixture_Code,data=dat_4,axes=F,at=position ,border="gray40",xlab= "Population", ylab="Carbon content [%]",las=2, 
        col=c("#0033FF","#0033FF","#0033FF","#0033FF","#FFFF00","#FFFF00","#FFFF00","#FFFF00"))
axis(2,at=seq(25,55,5), seq(25,55,5), las=1)
axis(1, at=c(1,2,3,4.5,5.5,7,8,9), labels=c("0:1","1:1","0:2","JE","GO","0:1","1:1","0:2"), las=1)
box()

#add overall mean for JE and GO
#model predicted means

lsm<-emmeans(M1b,~Population+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Population+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(4.5,5.5)
points(position,plot$response,pch=4,cex=2)

#add mean for JE incompetition with GO
M2c<-lm(log(C_.) ~  Mixture_Code1 +Initial_biomass,data=dat_2)

lsm<-emmeans(M2c,~Mixture_Code1+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Mixture_Code1+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(1,2,3)
points(position,plot$response,pch=4,cex=2)

#add mean for GO incompetition with JE
M3b<-lm(log(C_.) ~  Mixture_Code1 +Initial_biomass,data=dat_3)

lsm<-emmeans(M3b,~Mixture_Code1+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Mixture_Code1+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(7,8,9)
points(position,plot$response,pch=4,cex=2)

#CARBON TO NITROGEN RATIO

dat_1<-read.table(file="Combined_CN.txt", header=TRUE,sep="\t",dec = ".")

str(dat_1)

M1 <- lmer(C_N. ~ Population+(1 | Mixture_series)+(1 | Incompetition)+Initial_biomass, REML=F, data = dat_1)
M1a <- lmer(sqrt(C_N.) ~ Population+(1 | Mixture_series)+(1 | Incompetition)+Initial_biomass, REML=F, data = dat_1)
M1b <- lmer(log(C_N.) ~ Population+(1 | Mixture_series)+(1 | Incompetition)+Initial_biomass, REML=F, data = dat_1)

#Check themodel(DHARMarequired)
sim1a<-simulateResiduals(fittedModel= M1,n=1000) 
plot(sim1a, quantreg= FALSE)

sim1a<-simulateResiduals(fittedModel= M1a,n=1000)
plot(sim1a, quantreg= FALSE)

sim1a<-simulateResiduals(fittedModel= M1b,n=1000)
plot(sim1a, quantreg= FALSE)

summary(M1b)
#untransformed is better
#Significant test
Anova(M1,type="III",test="Chisq")

#multiple comparaison
library(multcomp)
gltt_M1<-glht(M1b,linfct=mcp(Population="Tukey"))
gltt_M1
summary(gltt_M1)
post<-cld(gltt_M1, level = 0.05, decreasing = F)
post
#

# use lm for JE incompetition GO  

dat_2<-read.table(file="CN_JE.txt", header=TRUE,sep="\t",dec = ".")
Mixture_Code1<-as.factor(dat_2$Mixture_series)
str(dat_2)
# lm
M2a<-lm(C_N_. ~  Mixture_Code1 +Initial_biomass,data=dat_2)

# Test homogeneity of variance
leveneTest(C_N_.~Mixture_Code1,data=dat_2)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.77 


#Check normality with shapiro test
shapiro.test(dat_2$C_N_.) 

#P VALUE > 0.05 the assumption is not valid, p =0.02

#transformation
M2b<-lm(sqrt(C_N_.) ~  Mixture_Code1 +Initial_biomass,data=dat_2)
M2c<-lm(log(C_N_.) ~  Mixture_Code1 +Initial_biomass,data=dat_2)
#LM
plot(M2a)  
hist(resid(M2a)) 

plot(M2b)  
hist(resid(M2b)) 

plot(M2c) # 
hist(resid(M2c)) 

summary(m4)

#log transformation is better
Anova(M2c,type="III",test="F")
##
DunnettTest(log(C_N_.) ~  Mixture_series,data=dat_2)
#multiple comparaison
Mixture_Code1<-as.factor(dat_2$Mixture_Code)
str(dat_2)
M2c<-lm(log(C_N_. )~  Mixture_Code1 +Initial_biomass,data=dat_2)

gltt_M2b<-glht(M2c,linfct=mcp(Mixture_Code1="Tukey"))
gltt_M2b
summary(gltt_M2b)
post<-cld(gltt_M2b, level = 0.05, decreasing = F)
post
#GO incompetition JE
dat_3<-read.table(file="CN_GO.txt", header=TRUE,sep="\t",dec = ".")
Mixture_Code1<-as.factor(dat_3$Mixture_series)
str(dat_3)
M3<-lm(C_N.  ~  Mixture_Code1 +Initial_biomass,data=dat_3)
M3a<-lm(sqrt(C_N. ) ~  Mixture_Code1 +Initial_biomass,data=dat_3)
M3b<-lm(log(C_N. ) ~  Mixture_Code1 +Initial_biomass,data=dat_3)
# Test homogeneity of variance
leveneTest(C_N.~Mixture_Code,data=dat_3)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.98

#Check normality with shapiro test
shapiro.test(dat_3$C_N.) 

#P VALUE > 0.05 the assumption is valid, p =0.01


plot(M3) # 
hist(resid(M3)) # 
plot(M3a) # 
hist(resid(M3a)) 
plot(M3b) # 
hist(resid(M3b)) 
summary(M3)

#log transformed data fit better
Anova(M3b,type="III",test="F")
##
DunnettTest(C_N.~  Mixture_series,data=dat_3)
#multiple comparaison
Mixture_Code1<-as.factor(dat_3$Mixture_Code)
str(dat_3)
M3b<-lm(log(N_.) ~  Mixture_Code1 +Initial_biomass,data=dat_3)

gltt_M3<-glht(M3b,linfct=mcp(Mixture_Code1="Tukey"))
gltt_M3
summary(gltt_M3)
post<-cld(gltt_M3, level = 0.05, decreasing = F)
post
#Data for Box_plot 
dat_4<-read.table(file="CN_Box.txt", header=TRUE,sep="\t",dec = ".")
str(dat_4)

boxplot(C_N.  ~Mixture_Code,data=dat_4)

position<-c(1,2,3,4.5,5.5,7,8,9)
boxplot(C_N.  ~Mixture_Code,data=dat_4,axes=F,at=position ,border="gray40",xlab= "Population", ylab="Carbon to nitrogen ratio [%]",las=2, 
        col=c("#0033FF","#0033FF","#0033FF","#0033FF","#FFFF00","#FFFF00","#FFFF00","#FFFF00"))
axis(2,at=seq(10,18,2), seq(10,18,2), las=1)
axis(1, at=c(1,2,3,4.5,5.5,7,8,9), labels=c("0:1","1:1","0:2","JE","GO","0:1","1:1","0:2"), las=1)
box()

#add overall mean for JE and GO
#model predicted means

lsm<-emmeans(M1,~Population+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Population+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(4.5,5.5)
points(position,plot$emmean,pch=4,cex=2)

#add mean for JE incompetition with GO
M2c<-lm(log(C_N_.) ~  Mixture_Code1 +Initial_biomass,data=dat_2)

lsm<-emmeans(M2c,~Mixture_Code1+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Mixture_Code1+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(1,2,3)
points(position,plot$response,pch=4,cex=2)

#add mean for GO incompetition with JE
M3b<-lm(log(C_N.) ~  Mixture_Code1 +Initial_biomass,data=dat_3)

lsm<-emmeans(M3b,~Mixture_Code1+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Mixture_Code1+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(7,8,9)
points(position,plot$response,pch=4,cex=2)



#MIXTURE SERIES COMBINED HARVEST 2

dat_1<-read.table(file="Combine_Series_H2.txt", header=TRUE,sep="\t",dec = ".")
Mixture_Code1<-as.factor(dat_1$Mixture_series)
str(dat_1)
#Mixed-effect analysis of variance (ANOVA)
#M1 <- lmer(Biomass ~ Population+(1 |Mixture_series )+(1 | Incompetition), data = dat_1,control = lmerControl(optimizer ="Nelder_Mead"))
#DRY BIOMASS
#Use lmer with RMEL=F
M1 <- lmer(Biomass ~ Population+(1 | Mixture_Code1)+(1 | Incompetition)+Initial_biomass, REML=F, data = dat_1)
M1a <- lmer(sqrt(Biomass) ~ Population+(1 | Mixture_Code1)+(1 | Incompetition)+Initial_biomass, REML=F, data = dat_1)
M1b <- lmer(log(Biomass) ~ Population+(1 | Mixture_Code1)+(1 | Incompetition)+Initial_biomass, REML=F, data = dat_1)

#Check themodel(DHARMarequired)
sim1a<-simulateResiduals(fittedModel= M1,n=1000)
plot(sim1a, quantreg= FALSE)

sim1a<-simulateResiduals(fittedModel= M1a,n=1000) #fit better
plot(sim1a, quantreg= FALSE)

sim1a<-simulateResiduals(fittedModel= M1b,n=1000)
plot(sim1a, quantreg= FALSE)

summary(M1a)
#sqrt transformation is better
#Significant test
Anova(M1a,type="III",test="Chisq")


# use lm for JE incompetition GO  

dat_2<-read.table(file="JE_H2.txt", header=TRUE,sep="\t",dec = ".")
Mixture_Code1<-as.factor(dat_2$Mixture_series)
str(dat_2)
# or lm
M2a<-lm(JE_Biomass ~  Mixture_Code1 +Initial_biomass,data=dat_2)

# Test homogeneity of variance
leveneTest(JE_Biomass~Mixture_Code,data=dat_2)
#Assumption is valid if pvalue >.05
# assumption is not valid, p =0.007


#Check normality with shapiro test
shapiro.test(dat_2$JE_Biomass) 

#P VALUE > 0.05 the assumption is not valid, p =0.005

#transformation
M2b<-lm(sqrt(JE_Biomass) ~  Mixture_Code1 +Initial_biomass,data=dat_2)
M2c<-lm(log(JE_Biomass) ~  Mixture_Code1 +Initial_biomass,data=dat_2)

plot(M2a) # 
hist(resid(M2a))

plot(M2b) # 
hist(resid(M2b)) 

plot(M2c) # 
hist(resid(M2c))
summary(m4)

#sqrt transformation fit better
Anova(M2b,type="III",test="F")
##
DunnettTest(sqrt(JE_Biomass) ~  Mixture_series,data=dat_2)

#multiple comparaison
Mixture_Code1<-as.factor(dat_2$Mixture_Code)
str(dat_2)
M2b<-lm(sqrt(JE_Biomass) ~  Mixture_Code1 +Initial_biomass,data=dat_2)

gltt_M2b<-glht(M2b,linfct=mcp(Mixture_Code1="Tukey"))
gltt_M2b
summary(gltt_M2b)
post<-cld(gltt_M2b, level = 0.05, decreasing = F)
post

#GO incompetition JE
dat_3<-read.table(file="GO_H2.txt", header=TRUE,sep="\t",dec = ".")
Mixture_Code1<-as.factor(dat_3$Mixture_series)
str(dat_3)
M3<-lm(GO_Biomass ~  Mixture_Code1 +Initial_biomass,data=dat_3)
M3a<-lm(sqrt(GO_Biomass) ~  Mixture_Code1 +Initial_biomass,data=dat_3)
M3b<-lm(log(GO_Biomass) ~  Mixture_series +Initial_biomass,data=dat_3)
# Test homogeneity of variance
leveneTest(GO_Biomass~Mixture_Code,data=dat_3)
#Assumption is valid if pvalue >.05
# assumption is not valid, p =0.01

#Check normality with shapiro test
shapiro.test(dat_3$GO_Biomass) 

#P VALUE > 0.05 the assumption is valid, p =0.18


plot(M3) # 
hist(resid(M3))
plot(M3a) # 
hist(resid(M3a)) #sqrt fit better
plot(M3b) # 
hist(resid(M3b))
summary(M3)

#sqrt transformed data fit better
Anova(M3a,type="III",test="F")
##
DunnettTest(sqrt(GO_Biomass)~  Mixture_series,data=dat_3)
#multiple comparaison
Mixture_Code1<-as.factor(dat_3$Mixture_Code)
str(dat_3)
M3a<-lm(sqrt(GO_Biomass) ~  Mixture_Code1 +Initial_biomass,data=dat_3)

gltt_M3a<-glht(M3a,linfct=mcp(Mixture_Code1="Tukey"))
gltt_M3a
summary(gltt_M3a)
post<-cld(gltt_M3a, level = 0.05, decreasing = F)
post
#Data for Box_plot 
dat_4b<-read.table(file="Data_Box_H2.txt", header=TRUE,sep="\t",dec = ".")
str(dat_4)

boxplot(Biomass ~Mixture_Code,data=dat_4)


position<-c(1,2,3,4,5,6.5,7.5,9,10,11,12,13)
Harvest2<-boxplot(Biomass ~Mixture_Code,data=dat_4b,axes=F,at=position ,border="gray40",xlab= "Population", ylab="Aboveground biomass [g DW]",las=2, 
        col=c("#0033FF","#0033FF","#0033FF","#0033FF","#0033FF","#0033FF","#FFFF00","#FFFF00","#FFFF00","#FFFF00","#FFFF00","#FFFF00"))
axis(2,at=seq(0,5,0.5), seq(0,5,0.5), las=1)
axis(1, at=c(1,2,3,4,5,6.5,7.5,9,10,11,12,13), labels=c("0:1","1:1","2:1","1:2","0:2","JE","GO","0:1","1:1","2:1","1:2","0:2"), las=1) +
box()

Harvest1a + Harvest2
#add overall mean for JE and GO
#model predicted means

lsm<-emmeans(M1a,~Population+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Population+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(6.5,7.5)
points(position,plot$response,pch=4,cex=2)

#add mean for JE incompetition with GO
M2a<-lm(JE_Biomass ~  Mixture_Code1 +Initial_biomass,data=dat_2)

lsm<-emmeans(M2a,~Mixture_Code+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Mixture_Code+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(1,2,3,4,5)
points(position,plot$emmean,pch=4,cex=2)
#means<-tapply(dat_2$JE_Biomass,dat_2$Mixture_Code,mean)
#means
#points(means,col="black",pch=4, at=position)

#add mean for GO incompetition with JE
#add mean for JE incompetition with GO
M3a<-lm(sqrt(GO_Biomass) ~  Mixture_Code +Initial_biomass,data=dat_3)

lsm<-emmeans(M3a,~Mixture_Code+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Mixture_Code+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(9,10,11,12,13)
points(position,plot$response,pch=4,cex=2)


## LEAF COUNT
str(dat_1)
#Use lmer with RMEL=F
M1 <- lmer(Leaf_Count ~ Population+(1 | Mixture_Code1)+(1 | Incompetition)+Initial_biomass, REML=F, data = dat_1)
M1a <- lmer(sqrt(Leaf_Count) ~ Population+(1 | Mixture_Code1)+(1 | Incompetition)+Initial_biomass, REML=F, data = dat_1)
M1b <- lmer(log(Leaf_Count) ~ Population+(1 | Mixture_Code1)+(1 | Incompetition)+Initial_biomass, REML=F, data = dat_1)

#Check themodel(DHARMarequired)
sim1a<-simulateResiduals(fittedModel= M1,n=1000) #fit better
plot(sim1a, quantreg= FALSE)

sim1a<-simulateResiduals(fittedModel= M1a,n=1000)
plot(sim1a, quantreg= FALSE)

sim1a<-simulateResiduals(fittedModel= M1b,n=1000)
plot(sim1a, quantreg= FALSE)

summary(M1a)
# untransformation is better
#Significant test
Anova(M1,type="III",test="Chisq")# no difference

#multiple comparaison
library(multcomp)
gltt_M1a<-glht(M1a,linfct=mcp(Population="Tukey"))
gltt_M1a
summary(gltt_M1a)
post<-cld(gltt_M1a, level = 0.05, decreasing = F)
post
#

# use lm for JE incompetition GO  
str(dat_2)
dat_2<-read.table(file="JE_H2.txt", header=TRUE,sep="\t",dec = ".")
Mixture_Code1<-as.factor(dat_2$Mixture_series)
str(dat_2)
M2a<-lm(JE_Leaf_Count ~  Mixture_Code1 +Initial_biomass,data=dat_2)

# Test homogeneity of variance
leveneTest(JE_Leaf_Count~Mixture_Code,data=dat_2)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.85


#Check normality with shapiro test
shapiro.test(dat_2$JE_Leaf_Count) 

#P VALUE > 0.05 the assumption is valid, p =0.07

#transformation
M2b<-lm(sqrt(JE_Leaf_Count) ~  Mixture_Code1 +Initial_biomass,data=dat_2)
M2c<-lm(log(JE_Leaf_Count) ~  Mixture_Code1 +Initial_biomass,data=dat_2)

plot(M2a)  
hist(resid(M2a))

plot(M2b)  
hist(resid(M2b)) 

plot(M2c) # 
hist(resid(M2c)) 
summary(m4)

#untransformed fit better
Anova(M2b,type="III",test="F")
##
DunnettTest(sqrt(JE_Leaf_Count) ~  Mixture_series,data=dat_2)
#multiple comparaison
Mixture_Code1<-as.factor(dat_2$Mixture_Code)
str(dat_2)
M2a<-lm(JE_Leaf_Count ~  Mixture_Code1 +Initial_biomass,data=dat_2)

gltt_M2b<-glht(M2a,linfct=mcp(Mixture_Code1="Tukey"))
gltt_M2b
summary(gltt_M2b)
post<-cld(gltt_M2b, level = 0.05, decreasing = F)
post
#GO incompetition JE
dat_3<-read.table(file="GO_H2.txt", header=TRUE,sep="\t",dec = ".")
Mixture_Code1<-as.factor(dat_3$Mixture_series)
str(dat_3)
M3<-lm(GO_Leaf_Count ~  Mixture_Code1 +Initial_biomass,data=dat_3)
M3a<-lm(sqrt(GO_Leaf_Count) ~  Mixture_Code1 +Initial_biomass,data=dat_3)
M3b<-lm(log(GO_Leaf_Count) ~  Mixture_Code1 +Initial_biomass,data=dat_3)
# Test homogeneity of variance
leveneTest(GO_Leaf_Count~Mixture_Code,data=dat_3)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.21

#Check normality with shapiro test
shapiro.test(dat_3$GO_Leaf_Count) 

#P VALUE > 0.05 the assumption is not valid, p =0.026


plot(M3) # 
hist(resid(M3))
plot(M3a) # 
hist(resid(M3a)) 
plot(M3b) # 
hist(resid(M3b)) #log fit better
summary(M3)

#log transformed data fit better
Anova(M3b,type="III",test="F")
##
DunnettTest(log(GO_Leaf_Count)~  Mixture_series,data=dat_3)
#multiple comparaison
Mixture_Code1<-as.factor(dat_3$Mixture_Code)
str(dat_3)
M3b<-lm(log(GO_Leaf_Count) ~  Mixture_Code1 +Initial_biomass,data=dat_3)

gltt_M3b<-glht(M3b,linfct=mcp(Mixture_Code1="Tukey"))
gltt_M3b
summary(gltt_M3b)
post<-cld(gltt_M3b, level = 0.05, decreasing = F)
post
#Data for Box_plot 
dat_4<-read.table(file="Data_Box_H2.txt", header=TRUE,sep="\t",dec = ".")
str(dat_4)

boxplot(Leaf_Count ~Mixture_Code,data=dat_4)


position<-c(1,2,3,4,5,6.5,7.5,9,10,11,12,13)
boxplot(Leaf_Count ~Mixture_Code,data=dat_4,axes=F,at=position ,border="gray40",xlab= "Population", ylab="Total number of leaves",las=2, 
        col=c("#0033FF","#0033FF","#0033FF","#0033FF","#0033FF","#0033FF","#FFFF00","#FFFF00","#FFFF00","#FFFF00","#FFFF00","#FFFF00"))
axis(2,at=seq(4,12,2), seq(4,12,2), las=1)
axis(1, at=c(1,2,3,4,5,6.5,7.5,9,10,11,12,13), labels=c("0:1","1:1","2:1","1:2","0:2","JE","GO","0:1","1:1","2:1","1:2","0:2"), las=1)
box()

#add overall mean for JE and GO
#model predicted means
M1 <- lmer(Leaf_Count ~ Population+(1 | Mixture_Code1)+(1 | Incompetition)+Initial_biomass, REML=F, data = dat_1)

lsm<-emmeans(M1,~Population+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Population+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(6.5,7.5)
points(position,plot$emmean,pch=4,cex=2)

#add mean for JE incompetition with GO
M2b<-lm(sqrt(JE_Leaf_Count) ~  Mixture_Code1 +Initial_biomass,data=dat_2)

lsm<-emmeans(M2b,~Mixture_Code1+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Mixture_Code1+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(1,2,3,4,5)
points(position,plot$response,pch=4,cex=2)


#add mean for GO incompetition with JE

M3b<-lm(log(GO_Leaf_Count) ~  Mixture_Code1 +Initial_biomass,data=dat_3)

lsm<-emmeans(M3b,~Mixture_Code1+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Mixture_Code1+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(9,10,11,12,13)
points(position,plot$response,pch=4,cex=2)

## LENGTH OF LONGEST LEAF
str(dat_1)
Mixture_Code1<-as.factor(dat_1$Mixture_series)
#Use lmer with RMEL=F
M1 <- lmer(LLL ~ Population+(1 | Mixture_Code1)+(1 | Incompetition)+Initial_biomass, REML=F, data = dat_1)
M1a <- lmer(sqrt(LLL) ~ Population+(1 | Mixture_Code1)+(1 | Incompetition)+Initial_biomass, REML=F, data = dat_1)
M1b <- lmer(log(LLL) ~ Population+(1 | Mixture_Code1)+(1 | Incompetition)+Initial_biomass, REML=F, data = dat_1)

#Check themodel(DHARMarequired)
sim1a<-simulateResiduals(fittedModel= M1,n=1000) #fit better
plot(sim1a, quantreg= FALSE)

sim1a<-simulateResiduals(fittedModel= M1a,n=1000)
plot(sim1a, quantreg= FALSE)

sim1a<-simulateResiduals(fittedModel= M1b,n=1000)
plot(sim1a, quantreg= FALSE)

summary(M1a)
#untransformed is better
#Significant test
Anova(M1,type="III",test="Chisq")

#multiple comparaison
library(multcomp)
gltt_M1<-glht(M1,linfct=mcp(Population="Tukey"))
gltt_M1
summary(gltt_M1)
post<-cld(gltt_M1, level = 0.05, decreasing = F)
post
#

# use lm for JE incompetition GO  
str(dat_2)
dat_2<-read.table(file="JE_H2.txt", header=TRUE,sep="\t",dec = ".")
Mixture_Code1<-as.factor(dat_2$Mixture_series)

# lm
M2a<-lm(JE_LLL ~  Mixture_Code1 +Initial_biomass,data=dat_2)

# Test homogeneity of variance
leveneTest(JE_LLL~Mixture_Code1,data=dat_2)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.28


#Check normality with shapiro test
shapiro.test(dat_2$JE_LLL) 

#P VALUE > 0.05 the assumption is not valid, p =0.005

#transformation
M2b<-lm(sqrt(JE_LLL) ~  Mixture_Code1 +Initial_biomass,data=dat_2)
M2c<-lm(log(JE_LLL) ~  Mixture_Code1 +Initial_biomass,data=dat_2)

plot(M2a)  
hist(resid(M2a)) 

plot(M2b)  
hist(resid(M2b)) 

plot(M2c) # 
hist(resid(M2c)) 
summary(m4)

#sqrt transformation fit better
Anova(M2b,type="III",test="F")
##
DunnettTest(sqrt(JE_LLL) ~  Mixture_series,data=dat_2)
#multiple comparaison
Mixture_Code1<-as.factor(dat_2$Mixture_Code)
str(dat_2)
M2b<-lm(sqrt(JE_LLL) ~  Mixture_Code1 +Initial_biomass,data=dat_2)

gltt_M2b<-glht(M2b,linfct=mcp(Mixture_Code1="Tukey"))
gltt_M2b
summary(gltt_M2b)
post<-cld(gltt_M2b, level = 0.05, decreasing = F)
post
#GO incompetition JE
dat_3<-read.table(file="GO_H2.txt", header=TRUE,sep="\t",dec = ".")
Mixture_Code1<-as.factor(dat_3$Mixture_series)
str(dat_3)
M3<-lm(GO_LLL ~  Mixture_Code1 +Initial_biomass,data=dat_3)
M3a<-lm(sqrt(GO_LLL) ~  Mixture_Code1 +Initial_biomass,data=dat_3)
M3b<-lm(log(GO_LLL) ~  Mixture_Code1 +Initial_biomass,data=dat_3)
# Test homogeneity of variance
leveneTest(GO_LLL~Mixture_Code,data=dat_3)
#Assumption is valid if pvalue >.05
# assumption is not valid, p =0.008

#Check normality with shapiro test
shapiro.test(dat_3$GO_LLL) 

#P VALUE > 0.05 the assumption is not valid, p =0.002


plot(M3) # 
hist(resid(M3)) 
plot(M3a) # 
hist(resid(M3a)) #fit better
plot(M3b) # 
hist(resid(M3b)) 
summary(M3)

#sqrt transformed data fit better
Anova(M3a,type="III",test="F")
##
DunnettTest(sqrt(GO_LLL)~  Mixture_series,data=dat_3)
#multiple comparaison
Mixture_Code1<-as.factor(dat_3$Mixture_Code)
str(dat_3)
M3a<-lm(sqrt(GO_LLL) ~  Mixture_Code1 +Initial_biomass,data=dat_3)

gltt_M3a<-glht(M3a,linfct=mcp(Mixture_Code1="Tukey"))
gltt_M3a
summary(gltt_M3a)
post<-cld(gltt_M3a, level = 0.05, decreasing = F)
post
#Data for Box_plot 
dat_4<-read.table(file="Data_Box_H2.txt", header=TRUE,sep="\t",dec = ".")
str(dat_4)

boxplot(LLL ~Mixture_Code,data=dat_4)


position<-c(1,2,3,4,5,6.5,7.5,9,10,11,12,13)
boxplot(LLL ~Mixture_Code,data=dat_4,axes=F,at=position ,border="gray40",xlab= "Population", ylab="Length of the longest leaf [cm]",las=2, 
        col=c("#0033FF","#0033FF","#0033FF","#0033FF","#0033FF","#0033FF","#FFFF00","#FFFF00","#FFFF00","#FFFF00","#FFFF00","#FFFF00"))
axis(2,at=seq(5,20,5), seq(5,20,5), las=1)
axis(1, at=c(1,2,3,4,5,6.5,7.5,9,10,11,12,13), labels=c("0:1","1:1","2:1","1:2","0:2","JE","GO","0:1","1:1","2:1","1:2","0:2"), las=1)
box()

#add overall mean for JE and GO
#model predicted means

lsm<-emmeans(M1,~Population+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Population+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(6.5,7.5)
points(position,plot$emmean,pch=4,cex=2)

#add mean for JE incompetition with GO
M2b<-lm(sqrt(JE_LLL) ~  Mixture_Code1 +Initial_biomass,data=dat_2)

lsm<-emmeans(M2b,~Mixture_Code1+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Mixture_Code1+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(1,2,3,4,5)
points(position,plot$response,pch=4,cex=2)

#add mean for GO incompetition with JE
M3a<-lm(sqrt(GO_LLL) ~  Mixture_Code1 +Initial_biomass,data=dat_3)

lsm<-emmeans(M3a,~Mixture_Code1+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Mixture_Code1+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(9,10,11,12,13)
points(position,plot$response,pch=4,cex=2)

### LEAF WATER CONTENT
Mixture_Code1<-as.factor(dat_1$Mixture_series)
str(dat_1)
#Use lmer with RMEL=F
M1 <- lmer(LWC ~ Population+(1 | Mixture_Code1)+(1 | Incompetition)+Initial_biomass, REML=F, data = dat_1)
M1a <- lmer(sqrt(LWC) ~ Population+(1 | Mixture_Code1)+(1 | Incompetition)+Initial_biomass, REML=F, data = dat_1)
M1b <- lmer(log(LWC) ~ Population+(1 | Mixture_Code1)+(1 | Incompetition)+Initial_biomass, REML=F, data = dat_1)

#Check themodel(DHARMarequired)
sim1a<-simulateResiduals(fittedModel= M1,n=1000) 
plot(sim1a, quantreg= FALSE)

sim1a<-simulateResiduals(fittedModel= M1a,n=1000)
plot(sim1a, quantreg= FALSE)

sim1a<-simulateResiduals(fittedModel= M1b,n=1000)
plot(sim1a, quantreg= FALSE)

summary(M1a)
#untransformed is better
#Significant test
Anova(M1,type="III",test="Chisq")

#multiple comparaison
library(multcomp)
gltt_M1<-glht(M1,linfct=mcp(Population="Tukey"))
gltt_M1
summary(gltt_M1)
post<-cld(gltt_M1, level = 0.05, decreasing = F)
post
#

# use lm for JE incompetition GO  

dat_2<-read.table(file="JE_H2.txt", header=TRUE,sep="\t",dec = ".")
Mixture_Code1<-as.factor(dat_2$Mixture_series)
str(dat_2)
# lm
M2a<-lm(JE_LWC ~  Mixture_Code1 +Initial_biomass,data=dat_2)

# Test homogeneity of variance
leveneTest(JE_LWC~Mixture_Code,data=dat_2)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.6


#Check normality with shapiro test
shapiro.test(dat_2$JE_LWC) 

#P VALUE > 0.05 the assumption is not valid, p =0.0003

#transformation
M2b<-lm(sqrt(JE_LWC) ~  Mixture_Code1 +Initial_biomass,data=dat_2)
M2c<-lm(log(JE_LWC) ~  Mixture_Code1 +Initial_biomass,data=dat_2)
#LM
plot(M2a)  
hist(resid(M2a)) 

plot(M2b)  
hist(resid(M2b)) #fit better

plot(M2c) # 
hist(resid(M2c)) 

summary(m4)

#SQRTtransformation fit better
Anova(M2b,type="III",test="F")
##
DunnettTest(sqrt(JE_LWC) ~  Mixture_series,data=dat_2)
#multiple comparaison
Mixture_Code1<-as.factor(dat_2$Mixture_Code)
str(dat_2)
M2b<-lm(sqrt(JE_LWC )~  Mixture_Code1 +Initial_biomass,data=dat_2)

gltt_M2b<-glht(M2b,linfct=mcp(Mixture_Code1="Tukey"))
gltt_M2b
summary(gltt_M2b)
post<-cld(gltt_M2b, level = 0.05, decreasing = F)
post
#GO incompetition JE
dat_3<-read.table(file="GO_H2.txt", header=TRUE,sep="\t",dec = ".")
Mixture_Code1<-as.factor(dat_3$Mixture_series)
str(dat_3)
M3<-lm(GO_LWC ~  Mixture_Code1 +Initial_biomass,data=dat_3)
M3a<-lm(sqrt(GO_LWC) ~  Mixture_Code1 +Initial_biomass,data=dat_3)
M3b<-lm(log(GO_LWC) ~  Mixture_Code1 +Initial_biomass,data=dat_3)
# Test homogeneity of variance
leveneTest(GO_LWC~Mixture_Code,data=dat_3)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.52

#Check normality with shapiro test
shapiro.test(dat_3$GO_LWC) 

#P VALUE > 0.05 the assumption is not valid, p =0.0001


plot(M3) # 
hist(resid(M3)) # fit better
plot(M3a) # 
hist(resid(M3a)) 
plot(M3b) # 
hist(resid(M3b)) 
summary(M3)

#untransformed data fit better
Anova(M3,type="III",test="F")
##
DunnettTest(GO_LWC~  Mixture_series,data=dat_3)
#multiple comparaison
Mixture_Code1<-as.factor(dat_3$Mixture_Code)
str(dat_3)
M3<-lm(GO_LWC ~  Mixture_Code1 +Initial_biomass,data=dat_3)

gltt_M3<-glht(M3,linfct=mcp(Mixture_Code1="Tukey"))
gltt_M3
summary(gltt_M3)
post<-cld(gltt_M3, level = 0.05, decreasing = F)
post
#Data for Box_plot 
dat_4<-read.table(file="Data_Box_H2.txt", header=TRUE,sep="\t",dec = ".")
str(dat_4)

boxplot(LWC ~Mixture_Code,data=dat_4)


position<-c(1,2,3,4,5,6.5,7.5,9,10,11,12,13)
boxplot(LWC ~Mixture_Code,data=dat_4,axes=F,at=position ,border="gray40",xlab= "Population", ylab="Leaf water content [%]",las=2, 
        col=c("#0033FF","#0033FF","#0033FF","#0033FF","#0033FF","#0033FF","#FFFF00","#FFFF00","#FFFF00","#FFFF00","#FFFF00","#FFFF00"))
axis(2,at=seq(60,100,10), seq(60,100,10), las=1)
axis(1, at=c(1,2,3,4,5,6.5,7.5,9,10,11,12,13), labels=c("0:1","1:1","2:1","1:2","0:2","JE","GO","0:1","1:1","2:1","1:2","0:2"), las=1)
box()

#add overall mean for JE and GO
#model predicted means

lsm<-emmeans(M1,~Population+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Population+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(6.5,7.5)
points(position,plot$emmean,pch=4,cex=2)

#add mean for JE incompetition with GO
M2b<-lm(sqrt(JE_LWC) ~  Mixture_Code1 +Initial_biomass,data=dat_2)

lsm<-emmeans(M2b,~Mixture_Code1+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Mixture_Code1+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(1,2,3,4,5)
points(position,plot$response,pch=4,cex=2)

#add mean for GO incompetition with JE
M3<-lm(GO_LWC ~  Mixture_Code1 +Initial_biomass,data=dat_3)

lsm<-emmeans(M3,~Mixture_Code1+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Mixture_Code1+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(9,10,11,12,13)
points(position,plot$emmean,pch=4,cex=2)

## LEAF MASS PER AREA (LMA)
dat_1a<-read.table(file="Combine_LMA_H2.txt", header=TRUE,sep="\t",dec = ".")

str(dat_1a)
#Use lmer with RMEL=F
M1 <- lmer(LMA  ~ Population+(1 | Mixture_series)+(1 | Incompetition)+Initial_biomass, REML=F, data = dat_1a)
M1a <- lmer(sqrt(LMA ) ~ Population+(1 | Mixture_series)+(1 | Incompetition)+Initial_biomass, REML=F, data = dat_1a)
M1b <- lmer(log(LMA ) ~ Population+(1 | Mixture_series)+(1 | Incompetition)+Initial_biomass, REML=F, data = dat_1a)

#Check themodel(DHARMarequired)
sim1a<-simulateResiduals(fittedModel= M1,n=1000) 
plot(sim1a, quantreg= FALSE)

sim1a<-simulateResiduals(fittedModel= M1a,n=1000)
plot(sim1a, quantreg= FALSE)

sim1a<-simulateResiduals(fittedModel= M1b,n=1000) # fit better
plot(sim1a, quantreg= FALSE)

summary(M1a)
#log transformed fit better
#Significant test
Anova(M1b,type="III",test="Chisq")

#multiple comparaison
library(multcomp)
gltt_M1<-glht(M1,linfct=mcp(Population="Tukey"))
gltt_M1
summary(gltt_M1)
post<-cld(gltt_M1, level = 0.05, decreasing = F)
post
#

# use lm for JE incompetition GO  
str(dat_2a)
dat_2a<-read.table(file="JE_LMA_H2.txt", header=TRUE,sep="\t",dec = ".")


# lm
M2a<-lm(LMA ~  Mixture_series +Initial_biomass,data=dat_2a)
M2b<-lm(sqrt(LMA) ~  Mixture_series +Initial_biomass,data=dat_2a)
M2c<-lm(log(LMA) ~  Mixture_series +Initial_biomass,data=dat_2a)

# Test homogeneity of variance
leveneTest(LMA~Mixture_Code,data=dat_2a)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.89


#Check normality with shapiro test
shapiro.test(dat_2a$LMA) 

#P VALUE > 0.05 the assumption is valid, p =0.12

#LM
plot(M2a)  
hist(resid(M2a)) 
plot(M2b)  
hist(resid(M2b)) #fit better
plot(M2c)  
hist(resid(M2c))
summary(m4)

# sqrt transformed fit better
Anova(M2b,type="III",test="F")
##
DunnettTest(sqrt(LMA) ~  Mixture_series,data=dat_2a)
#multiple comparaison
Mixture_Code1<-as.factor(dat_2a$Mixture_Code)
str(dat_2a)
M2b<-lm(sqrt(LMA) ~  Mixture_Code1 +Initial_biomass,data=dat_2a)

gltt_M2b<-glht(M2b,linfct=mcp(Mixture_Code1="Tukey"))
gltt_M2b
summary(gltt_M2b)
post<-cld(gltt_M2b, level = 0.05, decreasing = F)
post
#GO incompetition JE
dat_3<-read.table(file="GO_LMA_H2.txt", header=TRUE,sep="\t",dec = ".")
str(dat_3)
M3<-lm(LMA  ~  Mixture_series +Initial_biomass,data=dat_3)
M3a<-lm(sqrt(LMA ) ~  Mixture_series +Initial_biomass,data=dat_3)
M3b<-lm(log(LMA) ~  Mixture_series +Initial_biomass,data=dat_3)
# Test homogeneity of variance
leveneTest(LMA ~Mixture_Code,data=dat_3)
#Assumption is valid if pvalue >.05
# assumption is not valid, p =0.48

#Check normality with shapiro test
shapiro.test(dat_3$LMA ) 

#P VALUE > 0.05 the assumption is valid, p =0.0001


plot(M3) # 
hist(resid(M3)) #
plot(M3a) # 
hist(resid(M3a)) #
plot(M3b) # 
hist(resid(M3b)) # log fit better
summary(M3)

# log transformed data fit better
Anova(M3b,type="III",test="F")
##
DunnettTest(LMA ~  Mixture_series,data=dat_3)
#multiple comparaison
Mixture_Code1<-as.factor(dat_3$Mixture_Code)
str(dat_3)
M3b<-lm(log(LMA)  ~  Mixture_Code1 +Initial_biomass,data=dat_3)

gltt_M3b<-glht(M3b,linfct=mcp(Mixture_Code1="Tukey"))
gltt_M3b
summary(gltt_M3b)
post<-cld(gltt_M3b, level = 0.05, decreasing = F)
post
#Data for Box_plot 
dat_4<-read.table(file="Box_LMA_H2.txt", header=TRUE,sep="\t",dec = ".")
str(dat_4)

boxplot(LMA~Mixture_Code,data=dat_4)


position<-c(1,2,3,4.5,5.5,7,8,9)
boxplot(LMA ~Mixture_Code,data=dat_4,axes=F,at=position ,border="gray40",xlab= "Population", ylab="Leaf mass per area [mg mm-2]",las=2, 
        col=c("#0033FF","#0033FF","#0033FF","#0033FF","#FFFF00","#FFFF00","#FFFF00","#FFFF00"))
axis(2,at=seq(0.02,0.16,0.02), seq(0.02,0.16,0.02), las=1)
axis(1, at=c(1,2,3,4.5,5.5,7,8,9), labels=c("0:1","1:1","0:2","JE","GO","0:1","1:1","0:2"), las=1)
box()

#add overall mean for JE and GO
#model predicted means

lsm<-emmeans(M1b,~Population+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Population+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(4.5,5.5)
points(position,plot$response,pch=4,cex=2)

#add mean for JE incompetition with GO
M2b<-lm(sqrt(LMA) ~  Mixture_Code1 +Initial_biomass,data=dat_2a)

lsm<-emmeans(M2b,~Mixture_Code1+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Mixture_Code1+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(1,2,3)
points(position,plot$response,pch=4,cex=2)

#add mean for GO incompetition with JE
M3b<-lm(log(LMA) ~  Mixture_Code1 +Initial_biomass,data=dat_3)

lsm<-emmeans(M3b,~Mixture_Code1+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Mixture_Code1+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(7,8,9)
points(position,plot$response,pch=4,cex=2)
### Specific Leaf Area (SLA)
str(dat_1a)
#Use lmer with RMEL=F
M1 <- lmer(SLA  ~ Population+(1 | Mixture_series)+(1 | Incompetition)+Initial_biomass, REML=F, data = dat_1a)
M1a <- lmer(sqrt(SLA ) ~ Population+(1 | Mixture_series)+(1 | Incompetition)+Initial_biomass, REML=F, data = dat_1a)
M1b <- lmer(log(SLA ) ~ Population+(1 | Mixture_series)+(1 | Incompetition)+Initial_biomass, REML=F, data = dat_1a)

#Check themodel(DHARMarequired)
sim1a<-simulateResiduals(fittedModel= M1,n=1000) #fit better
plot(sim1a, quantreg= FALSE)

sim1a<-simulateResiduals(fittedModel= M1a,n=1000)
plot(sim1a, quantreg= FALSE)

sim1a<-simulateResiduals(fittedModel= M1b,n=1000)
plot(sim1a, quantreg= FALSE)

summary(M1a)
# sqrt transformed is better
#Significant test
Anova(M1a,type="III",test="Chisq")

#multiple comparaison
library(multcomp)
gltt_M1<-glht(M1a,linfct=mcp(Population="Tukey"))
gltt_M1
summary(gltt_M1)
post<-cld(gltt_M1, level = 0.05, decreasing = F)
post
#
# use lm for JE incompetition GO  
str(dat_2a)
dat_2a<-read.table(file="JE_LMA_H2.txt", header=TRUE,sep="\t",dec = ".")
Mixture_Code1<-as.factor(dat_2a$Mixture_series)
str(dat_2a)
# lm
M2a<-lm(SLA ~  Mixture_Code1 +Initial_biomass,data=dat_2a)
M2b<-lm(sqrt(SLA) ~  Mixture_Code1 +Initial_biomass,data=dat_2a)
M2c<-lm(log(SLA) ~  Mixture_Code1 +Initial_biomass,data=dat_2a)

# Test homogeneity of variance
leveneTest(SLA~Mixture_Code,data=dat_2a)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.66


#Check normality with shapiro test
shapiro.test(dat_2a$SLA) 

#P VALUE > 0.05 the assumption is valid, p =0.51

#LM
plot(M2a)  
hist(resid(M2a)) # untransformed fit better
plot(M2b)  
hist(resid(M2b)) 
plot(M2c)  
hist(resid(M2c))
summary(m4)

# untransformed fit better
Anova(M2a,type="III",test="F")
##
DunnettTest(SLA ~  Mixture_series,data=dat_2a)
#multiple comparaison
Mixture_Code1<-as.factor(dat_2a$Mixture_Code)
str(dat_2a)
M2a<-lm(SLA~  Mixture_Code1 +Initial_biomass,data=dat_2a)

gltt_M2b<-glht(M2a,linfct=mcp(Mixture_Code1="Tukey"))
gltt_M2b
summary(gltt_M2b)
post<-cld(gltt_M2b, level = 0.05, decreasing = F)
post
#GO incompetition JE
dat_3<-read.table(file="GO_LMA_H2.txt", header=TRUE,sep="\t",dec = ".")
Mixture_Code1<-as.factor(dat_3$Mixture_series)
str(dat_3)
M3<-lm(GO_SLA  ~  Mixture_Code1 +Initial_biomass,data=dat_3)
M3a<-lm(sqrt(GO_SLA ) ~  Mixture_Code1 +Initial_biomass,data=dat_3)
M3b<-lm(log(GO_SLA) ~  Mixture_Code1 +Initial_biomass,data=dat_3)
# Test homogeneity of variance
leveneTest(GO_SLA ~Mixture_Code,data=dat_3)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.37

#Check normality with shapiro test
shapiro.test(dat_3$GO_SLA ) 

#P VALUE > 0.05 the assumption is valid, p =0.0001


plot(M3) # 
hist(resid(M3)) #
plot(M3a) # 
hist(resid(M3a)) #
plot(M3b) # 
hist(resid(M3b)) # log fit better
summary(M3)

# log transformed data fit better
Anova(M3b,type="III",test="F")
##
DunnettTest(GO_SLA ~  Mixture_series,data=dat_3)
#multiple comparaison
Mixture_Code1<-as.factor(dat_3$Mixture_Code)
str(dat_3)
M3b<-lm(log(SLA)  ~  Mixture_Code1 +Initial_biomass,data=dat_3)

gltt_M3b<-glht(M3b,linfct=mcp(Mixture_Code1="Tukey"))
gltt_M3b
summary(gltt_M3b)
post<-cld(gltt_M3b, level = 0.05, decreasing = F)
post
#Data for Box_plot 
dat_4<-read.table(file="Box_LMA_H2.txt", header=TRUE,sep="\t",dec = ".")
str(dat_4)

boxplot(SLA~Mixture_Code,data=dat_4)


position<-c(1,2,3,4.5,5.5,7,8,9)
boxplot(SLA ~Mixture_Code,data=dat_4,axes=F,at=position ,border="gray40",xlab= "Population", ylab="Leaf mass per area [mg mm-2]",las=2, 
        col=c("#0033FF","#0033FF","#0033FF","#0033FF","#FFFF00","#FFFF00","#FFFF00","#FFFF00"))
axis(2,at=seq(10,60,10), seq(10,60,10), las=1)
axis(1, at=c(1,2,3,4.5,5.5,7,8,9), labels=c("0:1","1:1","0:2","JE","GO","0:1","1:1","0:2"), las=1)
box()

#add overall mean for JE and GO
#model predicted means

lsm<-emmeans(M1a,~Population+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Population+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(4.5,5.5)
points(position,plot$response,pch=4,cex=2)

#add mean for JE incompetition with GO
M2a<-lm(SLA ~  Mixture_Code1 +Initial_biomass,data=dat_2a)

lsm<-emmeans(M2a,~Mixture_Code1+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Mixture_Code1+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(1,2,3)
points(position,plot$emmean,pch=4,cex=2)

#add mean for GO incompetition with JE
M3b<-lm(log(GO_SLA) ~  Mixture_Code1 +Initial_biomass,data=dat_3)

lsm<-emmeans(M3b,~Mixture_Code1+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Mixture_Code1+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(7,8,9)
points(position,plot$response,pch=4,cex=2)


# use lm for JE incompetition GO  
str(dat_2a)
dat_2a<-read.table(file="Competitive_JE.txt", header=TRUE,sep="\t",dec = ".")


# lm
M2a<-lm(Competitive_ability_JE ~  Mixture_series +Initial_biomass,data=dat_2a)

# Test homogeneity of variance
leveneTest(Competitive_ability_JE~Mixture_Code,data=dat_2a)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.17 


#Check normality with shapiro test
shapiro.test(dat_2a$Competitive_ability_JE) 

#P VALUE > 0.05 the assumption is valid, p =0.6

#LM
plot(M2a)  
hist(resid(M2a)) 

summary(m4)

#untransformed fit better
Anova(M2a,type="III",test="F")
##
DunnettTest(Competitive_ability_JE ~  Mixture_series,data=dat_2a)
#multiple comparaison
Mixture_Code1<-as.factor(dat_2$Mixture_Code)
str(dat_2)
M2a<-lm(JE_SLA ~  Mixture_Code1 +Initial_biomass,data=dat_2)

gltt_M2a<-glht(M2a,linfct=mcp(Mixture_Code1="Tukey"))
gltt_M2a
summary(gltt_M2a)
post<-cld(gltt_M2a, level = 0.05, decreasing = F)
post
#GO incompetition JE
dat_3<-read.table(file="GO.txt", header=TRUE,sep="\t",dec = ".")
str(dat_3)
M3<-lm(GO_SLA  ~  Mixture_series +Initial_biomass,data=dat_3)
M3a<-lm(sqrt(GO_SLA ) ~  Mixture_series +Initial_biomass,data=dat_3)
M3b<-lm(log(GO_SLA ) ~  Mixture_series +Initial_biomass,data=dat_3)
# Test homogeneity of variance
leveneTest(GO_SLA ~Mixture_Code,data=dat_3)
#Assumption is valid if pvalue >.05
# assumption is not valid, p =0.04

#Check normality with shapiro test
shapiro.test(dat_3$GO_SLA ) 

#P VALUE > 0.05 the assumption is valid, p =0.04


plot(M3) # 
hist(resid(M3)) #
plot(M3a) # 
hist(resid(M3a)) 
plot(M3b) # 
hist(resid(M3b)) #fit better
summary(M3)

#transformed data fit better
Anova(M3b,type="III",test="F")
##
DunnettTest(GO_SLA ~  Mixture_series,data=dat_3)# not needed
#multiple comparaison
Mixture_Code1<-as.factor(dat_3$Mixture_Code)
str(dat_3)
M3b<-lm(log(GO_SLA)  ~  Mixture_Code1 +Initial_biomass,data=dat_3)

gltt_M3b<-glht(M3b,linfct=mcp(Mixture_Code1="Tukey"))
gltt_M3b
summary(gltt_M3b)
post<-cld(gltt_M3b, level = 0.05, decreasing = F)
post
#Data for Box_plot 
dat_4<-read.table(file="Data_Box.txt", header=TRUE,sep="\t",dec = ",")
str(dat_4)

boxplot(SLA ~Mixture_Code,data=dat_4)


position<-c(1,2,3,4,5,6.5,7.5,9,10,11,12,13)
boxplot(SLA ~Mixture_Code,data=dat_4,axes=F,at=position ,border="gray40",xlab= "Population", ylab="Specific leaf area [mm2mg-1]",las=2, 
        col=c("#0033FF","#0033FF","#0033FF","#0033FF","#0033FF","#0033FF","#FFFF00","#FFFF00","#FFFF00","#FFFF00","#FFFF00","#FFFF00"))
axis(2,at=seq(20,60,10), seq(20,60,10), las=1)
axis(1, at=c(1,2,3,4,5,6.5,7.5,9,10,11,12,13), labels=c("0:1","1:1","2:1","1:2","0:2","JE","GO","0:1","1:1","2:1","1:2","0:2"), las=1)
box()

#add overall mean for JE and GO
#model predicted means

lsm<-emmeans(M1,~Population+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Population+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(6.5,7.5)
points(position,plot$emmean,pch=4,cex=2)

#add mean for JE incompetition with GO
M2a<-lm(JE_SLA ~  Mixture_Code1 +Initial_biomass,data=dat_2)

lsm<-emmeans(M2a,~Mixture_Code1+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Mixture_Code1+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(1,2,3,4,5)
points(position,plot$emmean,pch=4,cex=2)

#add mean for GO incompetition with JE
M3b<-lm(log(GO_SLA) ~  Mixture_Code1 +Initial_biomass,data=dat_3)

lsm<-emmeans(M3b,~Mixture_Code1+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Mixture_Code1+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(9,10,11,12,13)
points(position,plot$response,pch=4,cex=2)
##
## NITROGEN TO CARBON CONTENT
dat_1<-read.table(file="Combine_CN_H2.txt", header=TRUE,sep="\t",dec = ".")


#NITROGEN
Mixture_Code1<-as.factor(dat_1$Mixture_series)
str(dat_1)
M1 <- lmer(N_. ~ Population+(1 | Mixture_series)+(1 | Incompetition)+Initial_biomass, REML=F, data = dat_1)
M1a <- lmer(sqrt(N_.) ~ Population+(1 | Mixture_series)+(1 | Incompetition)+Initial_biomass, REML=F, data = dat_1)
M1b <- lmer(log(N_.) ~ Population+(1 | Mixture_series)+(1 | Incompetition)+Initial_biomass, REML=F, data = dat_1)

#Check themodel(DHARMarequired)
sim1a<-simulateResiduals(fittedModel= M1,n=1000) 
plot(sim1a, quantreg= FALSE)

sim1a<-simulateResiduals(fittedModel= M1a,n=1000)
plot(sim1a, quantreg= FALSE)

sim1a<-simulateResiduals(fittedModel= M1b,n=1000) # log fit better
plot(sim1a, quantreg= FALSE)

summary(M1a)
#LOG transformed is better
#Significant test
Anova(M1b,type="III",test="Chisq")

#multiple comparaison
library(multcomp)
gltt_M1<-glht(M1b,linfct=mcp(Population="Tukey"))
gltt_M1
summary(gltt_M1)
post<-cld(gltt_M1, level = 0.05, decreasing = F)
post
#

# use lm for JE incompetition GO  
str(dat_2)
dat_2<-read.table(file="CN_JE_H2.txt", header=TRUE,sep="\t",dec = ".")
Mixture_Code1<-as.factor(dat_2$Mixture_series)

# lm
M2a<-lm(N_. ~  Mixture_Code1 +Initial_biomass,data=dat_2)

# Test homogeneity of variance
leveneTest(N_.~Mixture_Code,data=dat_2)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.19 


#Check normality with shapiro test
shapiro.test(dat_2$N_.) 

#P VALUE > 0.05 the assumption is not valid, p =0.98

#transformation
M2b<-lm(sqrt(N_.) ~  Mixture_Code1 +Initial_biomass,data=dat_2)
M2c<-lm(log(N_.) ~  Mixture_Code1 +Initial_biomass,data=dat_2)
#LM
plot(M2a)  
hist(resid(M2a)) 

plot(M2b)  
hist(resid(M2b)) #fit better

plot(M2c) # 
hist(resid(M2c)) 

summary(m4)

#sqrt transformation fit better
Anova(M2b,type="III",test="F")
##
DunnettTest(sqrt(N_.) ~  Mixture_series,data=dat_2)
#multiple comparaison
Mixture_Code1<-as.factor(dat_2$Mixture_Code)
str(dat_2)
M2b<-lm(sqrt(N_. )~  Mixture_Code1 +Initial_biomass,data=dat_2)

gltt_M2b<-glht(M2b,linfct=mcp(Mixture_Code1="Tukey"))
gltt_M2b
summary(gltt_M2b)
post<-cld(gltt_M2b, level = 0.05, decreasing = F)
post
#GO incompetition JE
dat_3<-read.table(file="CN_GO_H2.txt", header=TRUE,sep="\t",dec = ".")
Mixture_Code1<-as.factor(dat_3$Mixture_series)
str(dat_3)
M3<-lm(N_.  ~  Mixture_Code1 +Initial_biomass,data=dat_3)
M3a<-lm(sqrt(N_. ) ~  Mixture_Code1 +Initial_biomass,data=dat_3)
M3b<-lm(log(N_. ) ~  Mixture_Code1 +Initial_biomass,data=dat_3)
# Test homogeneity of variance
leveneTest(N_.~Mixture_Code,data=dat_3)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.34

#Check normality with shapiro test
shapiro.test(dat_3$N_.) 

#P VALUE > 0.05 the assumption is not valid, p =0.01


plot(M3) # 
hist(resid(M3)) # 
plot(M3a) # 
hist(resid(M3a)) 
plot(M3b) # 
hist(resid(M3b)) 
summary(M3)

#log transformed data fit better
Anova(M3b,type="III",test="F")
##
DunnettTest(N_.~  Mixture_series,data=dat_3)
#multiple comparaison
Mixture_Code1<-as.factor(dat_3$Mixture_Code)
str(dat_3)
M3b<-lm(log(N_.) ~  Mixture_Code1 +Initial_biomass,data=dat_3)

gltt_M3<-glht(M3b,linfct=mcp(Mixture_Code1="Tukey"))
gltt_M3
summary(gltt_M3)
post<-cld(gltt_M3, level = 0.05, decreasing = F)
post
#Data for Box_plot 
dat_4<-read.table(file="CN_Box_H2.txt", header=TRUE,sep="\t",dec = ".")
str(dat_4)

boxplot(N_. ~Mixture_Code,data=dat_4)

position<-c(1,2,3,4.5,5.5,7,8,9)
boxplot(N_. ~Mixture_Code,data=dat_4,axes=F,at=position ,border="gray40",xlab= "Population", ylab="Nitrogen content [%]",las=2, 
        col=c("#0033FF","#0033FF","#0033FF","#0033FF","#FFFF00","#FFFF00","#FFFF00","#FFFF00"))
axis(2,at=seq(1.0,3.5,0.5), seq(1.0,3.5,0.5), las=1)
axis(1, at=c(1,2,3,4.5,5.5,7,8,9), labels=c("0:1","1:1","0:2","JE","GO","0:1","1:1","0:2"), las=1)
box()

#add overall mean for JE and GO
#model predicted means

lsm<-emmeans(M1b,~Population+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Population+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(4.5,5.5)
points(position,plot$response,pch=4,cex=2)

#add mean for JE incompetition with GO
M2b<-lm(sqrt(N_.) ~  Mixture_Code1 +Initial_biomass,data=dat_2)

lsm<-emmeans(M2b,~Mixture_Code1+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Mixture_Code1+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(1,2,3)
points(position,plot$response,pch=4,cex=2)

#add mean for GO incompetition with JE
M3b<-lm(log(N_.) ~  Mixture_Code1 +Initial_biomass,data=dat_3)

lsm<-emmeans(M3b,~Mixture_Code1+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Mixture_Code1+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(7,8,9)
points(position,plot$response,pch=4,cex=2)


#CARBON CONTENT

dat_1<-read.table(file="Combine_CN_H2.txt", header=TRUE,sep="\t",dec = ".")
Mixture_Code1<-as.factor(dat_3$Mixture_series)
str(dat_1)

M1 <- lmer(C_. ~ Population+(1 | Mixture_Code1)+(1 | Incompetition)+Initial_biomass, REML=F, data = dat_1)
M1a <- lmer(sqrt(C_.) ~ Population+(1 | Mixture_Code1)+(1 | Incompetition)+Initial_biomass, REML=F, data = dat_1)
M1b <- lmer(log(C_.) ~ Population+(1 | Mixture_Code1)+(1 | Incompetition)+Initial_biomass, REML=F, data = dat_1)

#Check themodel(DHARMarequired)
sim1a<-simulateResiduals(fittedModel= M1,n=1000) #fit better
plot(sim1a, quantreg= FALSE)

sim1a<-simulateResiduals(fittedModel= M1a,n=1000)
plot(sim1a, quantreg= FALSE)

sim1a<-simulateResiduals(fittedModel= M1b,n=1000)
plot(sim1a, quantreg= FALSE)

summary(M1b)
#untransformed is better
#Significant test
Anova(M1,type="III",test="Chisq")

#multiple comparaison
library(multcomp)
gltt_M1<-glht(M1,linfct=mcp(Population="Tukey"))
gltt_M1
summary(gltt_M1)
post<-cld(gltt_M1, level = 0.05, decreasing = F)
post
#

# use lm for JE incompetition GO  

dat_2x<-read.table(file="CN_JE_H2.txt", header=TRUE,sep="\t",dec = ".")
Mixture_Code1<-as.factor(dat_2x$Mixture_series)
str(dat_2x)
# lm
M2a<-lm(C_. ~  Mixture_Code1 +Initial_biomass,data=dat_2x)

# Test homogeneity of variance
leveneTest(C_.~Mixture_Code,data=dat_2x)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.18


#Check normality with shapiro test
shapiro.test(dat_2x$C_.) 

#P VALUE > 0.05 the assumption is not valid, p =0.0001

#transformation
M2b<-lm(sqrt(C_.) ~  Mixture_Code1 +Initial_biomass,data=dat_2x)
M2c<-lm(log(C_.) ~  Mixture_Code1 +Initial_biomass,data=dat_2x)
#LM
plot(M2a)  
hist(resid(M2a)) 

plot(M2b)  
hist(resid(M2b)) 

plot(M2c) # 
hist(resid(M2c)) 

summary(m4)

#log transformation is better
Anova(M2c,type="III",test="F")
##
DunnettTest(log(N_.) ~  Mixture_series,data=dat_2x)
#multiple comparaison
Mixture_Code1<-as.factor(dat_2x$Mixture_Code)
str(dat_2x)
M2c<-lm(log(N_. )~  Mixture_Code1 +Initial_biomass,data=dat_2x)

gltt_M2b<-glht(M2c,linfct=mcp(Mixture_Code1="Tukey"))
gltt_M2b
summary(gltt_M2b)
post<-cld(gltt_M2b, level = 0.05, decreasing = F)
post
#GO incompetition JE
dat_3<-read.table(file="CN_GO_H2.txt", header=TRUE,sep="\t",dec = ".")
Mixture_Code1<-as.factor(dat_3$Mixture_series)
str(dat_3)
M3<-lm(C_.  ~  Mixture_series +Initial_biomass,data=dat_3)
M3a<-lm(sqrt(C_. ) ~  Mixture_Code1 +Initial_biomass,data=dat_3)
M3b<-lm(log(C_. ) ~  Mixture_Code1 +Initial_biomass,data=dat_3)
# Test homogeneity of variance
leveneTest(C_.~Mixture_Code,data=dat_3)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.67

#Check normality with shapiro test
shapiro.test(dat_3$C_.) 

#P VALUE > 0.05 the assumption is valid, p =0.08


plot(M3) # 
hist(resid(M3)) # fit better
plot(M3a) # 
hist(resid(M3a)) 
plot(M3b) # 
hist(resid(M3b)) 
summary(M3)

#untransformed data fit better
Anova(M3,type="III",test="F")
##
DunnettTest(C_.~  Mixture_series,data=dat_3)
#multiple comparaison
Mixture_Code1<-as.factor(dat_3$Mixture_Code)
str(dat_3)
M3<-lm(N_. ~  Mixture_Code1 +Initial_biomass,data=dat_3)

gltt_M<-glht(M3,linfct=mcp(Mixture_Code1="Tukey"))
gltt_M
summary(gltt_M3)
post<-cld(gltt_M3, level = 0.05, decreasing = F)
post
#Data for Box_plot 
dat_4<-read.table(file="CN_Box_H2.txt", header=TRUE,sep="\t",dec = ".")
str(dat_4)

boxplot(C_. ~Mixture_Code,data=dat_4)

position<-c(1,2,3,4.5,5.5,7,8,9)
boxplot(C_. ~Mixture_Code,data=dat_4,axes=F,at=position ,border="gray40",xlab= "Population", ylab="Carbon content [%]",las=2, 
        col=c("#0033FF","#0033FF","#0033FF","#0033FF","#FFFF00","#FFFF00","#FFFF00","#FFFF00"))
axis(2,at=seq(28,40,2), seq(28,40,2), las=1)
axis(1, at=c(1,2,3,4.5,5.5,7,8,9), labels=c("0:1","1:1","0:2","JE","GO","0:1","1:1","0:2"), las=1)
box()

#add overall mean for JE and GO
#model predicted means

lsm<-emmeans(M1,~Population+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Population+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(4.5,5.5)
points(position,plot$emmean,pch=4,cex=2)

#add mean for JE incompetition with GO
M2c<-lm(log(C_.) ~  Mixture_Code1 +Initial_biomass,data=dat_2x)

lsm<-emmeans(M2c,~Mixture_Code1+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Mixture_Code1+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(1,2,3)
points(position,plot$response,pch=4,cex=2)

#add mean for GO incompetition with JE
M3<-lm(C_. ~  Mixture_Code1 +Initial_biomass,data=dat_3)

lsm<-emmeans(M3,~Mixture_Code1+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Mixture_Code1+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(7,8,9)
points(position,plot$emmean,pch=4,cex=2)

#CARBON TO NITROGEN RATIO

dat_1<-read.table(file="Combine_CN_H2.txt", header=TRUE,sep="\t",dec = ".")

str(dat_1)

M1 <- lmer(C_N. ~ Population+(1 | Mixture_series)+(1 | Incompetition)+Initial_biomass, REML=F, data = dat_1)
M1a <- lmer(sqrt(C_N.) ~ Population+(1 | Mixture_series)+(1 | Incompetition)+Initial_biomass, REML=F, data = dat_1)
M1b <- lmer(log(C_N.) ~ Population+(1 | Mixture_series)+(1 | Incompetition)+Initial_biomass, REML=F, data = dat_1)

#Check themodel(DHARMarequired)
sim1a<-simulateResiduals(fittedModel= M1,n=1000) #fit better
plot(sim1a, quantreg= FALSE)

sim1a<-simulateResiduals(fittedModel= M1a,n=1000)
plot(sim1a, quantreg= FALSE)

sim1a<-simulateResiduals(fittedModel= M1b,n=1000)
plot(sim1a, quantreg= FALSE)

summary(M1b)
#untransformed is better
#Significant test
Anova(M1,type="III",test="Chisq")

#multiple comparaison
library(multcomp)
gltt_M1<-glht(M1b,linfct=mcp(Population="Tukey"))
gltt_M1
summary(gltt_M1)
post<-cld(gltt_M1, level = 0.05, decreasing = F)
post
#

# use lm for JE incompetition GO  

dat_2x<-read.table(file="CN_JE_H2.txt", header=TRUE,sep="\t",dec = ".")
Mixture_Code1<-as.factor(dat_2x$Mixture_series)
str(dat_2x)
# lm
M2a<-lm(C_N. ~  Mixture_Code1 +Initial_biomass,data=dat_2x)

# Test homogeneity of variance
leveneTest(C_N.~Mixture_Code,data=dat_2x)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.08


#Check normality with shapiro test
shapiro.test(dat_2x$C_N.) 

#P VALUE > 0.05 the assumption is not valid, p =0.003

#transformation
M2b<-lm(sqrt(C_N.) ~  Mixture_Code1 +Initial_biomass,data=dat_2x)
M2c<-lm(log(C_N.) ~  Mixture_Code1 +Initial_biomass,data=dat_2x)
#LM
plot(M2a)  
hist(resid(M2a)) 

plot(M2b)  
hist(resid(M2b)) #fit better

plot(M2c) # 
hist(resid(M2c)) 

summary(m4)

#sqrt transformation is better
Anova(M2b,type="III",test="F")
##
DunnettTest(sqrt(C_N.) ~  Mixture_series,data=dat_2x)
#multiple comparaison
Mixture_Code1<-as.factor(dat_2x$Mixture_Code)
str(dat_2x)
M2b<-lm(sqrt(C_N. )~  Mixture_Code1 +Initial_biomass,data=dat_2x)

gltt_M2b<-glht(M2b,linfct=mcp(Mixture_Code1="Tukey"))
gltt_M2b
summary(gltt_M2b)
post<-cld(gltt_M2b, level = 0.05, decreasing = F)
post
#GO incompetition JE
dat_3<-read.table(file="CN_GO_H2.txt", header=TRUE,sep="\t",dec = ".")
Mixture_Code1<-as.factor(dat_3$Mixture_series)
str(dat_3)
M3<-lm(C_N.  ~  Mixture_Code1 +Initial_biomass,data=dat_3)
M3a<-lm(sqrt(C_N. ) ~  Mixture_Code1 +Initial_biomass,data=dat_3)
M3b<-lm(log(C_N. ) ~  Mixture_Code1 +Initial_biomass,data=dat_3)
# Test homogeneity of variance
leveneTest(C_N.~Mixture_Code,data=dat_3)
#Assumption is valid if pvalue >.05
# assumption is not valid, p =0.01

#Check normality with shapiro test
shapiro.test(dat_3$C_N.) 

#P VALUE > 0.05 the assumption is valid, p =0.16


plot(M3) # 
hist(resid(M3)) # 
plot(M3a) # 
hist(resid(M3a)) 
plot(M3b) # 
hist(resid(M3b)) 
summary(M3)

#un transformed data fit better
Anova(M3,type="III",test="F")
##
DunnettTest(C_N.~  Mixture_series,data=dat_3)
#multiple comparaison
Mixture_Code1<-as.factor(dat_3$Mixture_Code)
str(dat_3)
M3b<-lm(C_N. ~  Mixture_Code1 +Initial_biomass,data=dat_3)

gltt_M3<-glht(M3b,linfct=mcp(Mixture_Code1="Tukey"))
gltt_M3
summary(gltt_M3)
post<-cld(gltt_M3, level = 0.05, decreasing = F)
post
#Data for Box_plot 
dat_4<-read.table(file="CN_Box_H2.txt", header=TRUE,sep="\t",dec = ".")
str(dat_4)

boxplot(C_N.  ~Mixture_Code,data=dat_4)

position<-c(1,2,3,4.5,5.5,7,8,9)
boxplot(C_N.  ~Mixture_Code,data=dat_4,axes=F,at=position ,border="gray40",xlab= "Population", ylab="Carbon to nitrogen ratio [%]",las=2, 
        col=c("#0033FF","#0033FF","#0033FF","#0033FF","#FFFF00","#FFFF00","#FFFF00","#FFFF00"))
axis(2,at=seq(10,40,5), seq(10,40,5), las=1)
axis(1, at=c(1,2,3,4.5,5.5,7,8,9), labels=c("0:1","1:1","0:2","JE","GO","0:1","1:1","0:2"), las=1)
box()

#add overall mean for JE and GO
#model predicted means

lsm<-emmeans(M1,~Population+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Population+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(4.5,5.5)
points(position,plot$emmean,pch=4,cex=2)

#add mean for JE incompetition with GO
M2b<-lm(sqrt(C_N.) ~  Mixture_Code1 +Initial_biomass,data=dat_2x)

lsm<-emmeans(M2b,~Mixture_Code1+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Mixture_Code1+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(1,2,3)
points(position,plot$response,pch=4,cex=2)

#add mean for GO incompetition with JE
M3b<-lm(C_N. ~  Mixture_Code1 +Initial_biomass,data=dat_3)

lsm<-emmeans(M3b,~Mixture_Code1+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Mixture_Code1+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(7,8,9)
points(position,plot$emmean,pch=4,cex=2)


#RELATIVE COMPETITION INTENSITY BASED ON CONTROL
#Biomass was considered as measure of yield
dat_1a<-read.table(file="Combined_Competitive_Contr_H2.txt", header=TRUE,sep="\t",dec = ".")

str(dat_1a)
#Use lmer with RMEL=F
MC <- lmer(Competitive_ability_JE  ~ Population+(1 | Mixture_series)+(1 | Incompetition)+Initial_biomass, REML=F, data = dat_1a)
MCb <- lmer(sqrt(Competitive_ability_JE)  ~ Population+(1 | Mixture_series)+(1 | Incompetition)+Initial_biomass, REML=F, data = dat_1a)

MCc <- lmer(log(Competitive_ability_JE)  ~ Population+(1 | Mixture_series)+(1 | Incompetition)+Initial_biomass, REML=F, data = dat_1a)


#Check themodel(DHARMarequired)
sim1a<-simulateResiduals(fittedModel= MC,n=1000) 
plot(sim1a, quantreg= FALSE)

sim1b<-simulateResiduals(fittedModel= MCb,n=1000) #fit better
plot(sim1b, quantreg= FALSE)

sim1ac<-simulateResiduals(fittedModel= MCc,n=1000) 
plot(sim1a, quantreg= FALSE)
summary(MC)
#SQRT transformed is better
#Significant test
Anova(MCb,type="III",test="Chisq")

#multiple comparaison
library(multcomp)
gltt_M1<-glht(MC,linfct=mcp(Population="Tukey"))
gltt_M1
summary(gltt_M1)
post<-cld(gltt_M1, level = 0.05, decreasing = F)
post
#

# use lm for JE incompetition GO  

dat_2<-read.table(file="Competitive_Contr_JE_H2.txt", header=TRUE,sep="\t",dec = ".")
dat_2a<-dat_2[!is.na(dat_2$Competitive_ability_JE),]# remove na
Mixture_Code1<-as.factor(dat_2a$Mixture_series)
str(dat_2)
# lm

M2C<-lm(Competitive_ability_JE ~  Mixture_Code1 +Initial_biomass,data=dat_2a)
M2Cb<-lm(sqrt(Competitive_ability_JE) ~  Mixture_Code1 +Initial_biomass,data=dat_2a)
M2Cc<-lm(log(Competitive_ability_JE) ~  Mixture_Code1 +Initial_biomass,data=dat_2a)

# Test homogeneity of variance
leveneTest(Competitive_ability_JE~Mixture_Code,data=dat_2a)
#Assumption is valid if pvalue >.05
# assumption is valid, p =0.88


#Check normality with shapiro test
shapiro.test(dat_2a$Competitive_ability_JE) 

#P VALUE > 0.05 the assumption is NOT valid, p =0.0001

#LM
plot(M2C)  
hist(resid(M2C)) 
plot(M2Cb)  
hist(resid(M2Cb)) # fit better
plot(M2Cc)  
hist(resid(M2Cc))
summary(M2C)

#sqrt transformed fit better
Anova(M2Cb,type="III",test="F")
##
DunnettTest(Competitive_ability_JE ~  Mixture_series,data=dat_2a)
#multiple comparaison
Mixture_Code1<-as.factor(dat_2a$Mixture_Code)
str(dat_2a)
M2Cb<-lm(sqrt(Competitive_ability_JE) ~  Mixture_Code1 +Initial_biomass,data=dat_2a)

gltt_M2a<-glht(M2Cb,linfct=mcp(Mixture_Code1="Tukey"))
gltt_M2a
summary(gltt_M2a)
post<-cld(gltt_M2a, level = 0.05, decreasing = F)
post
#COMPETITIVE INTENSITY OF GO incompetition JE
dat_3a<-read.table(file="Competitive_Contr_GO_H2.txt", header=TRUE,sep="\t",dec = ".")
Mixture_Code1<-as.factor(dat_3a$Mixture_series)
str(dat_3a)
M3C<-lm(Competitive_ability_GO  ~  Mixture_Code1 +Initial_biomass,data=dat_3a)
# Test homogeneity of variance
leveneTest(Competitive_ability_GO ~Mixture_Code,data=dat_3a)
#Assumption is valid if pvalue >.05
# assumption is NOT valid, p =0.01

#Check normality with shapiro test
shapiro.test(dat_3a$Competitive_ability_GO ) 

#P VALUE > 0.05 the assumption is not valid, p =0.0001

M3Cb<-lm(sqrt(Competitive_ability_GO)  ~  Mixture_Code1 +Initial_biomass,data=dat_3a)
M3Cc<-lm(log(Competitive_ability_GO)  ~  Mixture_Code1 +Initial_biomass,data=dat_3a)

plot(M3C) # 
hist(resid(M3C)) #

plot(M3Cb) # 
hist(resid(M3Cb)) #fit better

plot(M3Cc) # 
hist(resid(M3Cc))
summary(M3)

#sqrt transformed data fit better
Anova(M3C,type="III",test="F")
##
DunnettTest(Competitive_ability_GO ~  Mixture_series,data=dat_3a)
#multiple comparaison
Mixture_Code1<-as.factor(dat_3a$Mixture_Code)
str(dat_3a)
M3Cb<-lm(sqrt(Competitive_ability_GO)  ~  Mixture_Code1 +Initial_biomass,data=dat_3a)

gltt_M3C<-glht(M3Cb,linfct=mcp(Mixture_Code1="Tukey"))
gltt_M3C
summary(gltt_M3C)
post<-cld(gltt_M3C, level = 0.05, decreasing = F)
post
#Data for Box_plot 
dat_4a<-read.table(file="Competitive_Box_Contr_H2.txt", header=TRUE,sep="\t",dec = ".")
str(dat_4a)

boxplot(Competitive_ability ~Mixture_Code,data=dat_4a)


position<-c(1,2,3,4,6,7,9,10,11,12)
boxplot(Competitive_ability ~Mixture_Code,data=dat_4a,axes=F,at=position ,border="gray40",xlab= "Population", ylab="Relative competitive intensity",las=2, 
        col=c("#0033FF","#0033FF","#0033FF","#0033FF","#0033FF","#FFFF00","#FFFF00","#FFFF00","#FFFF00","#FFFF00"))
axis(2,at=seq(-6,2,2), seq(-6,2,2), las=1)
axis(1, at=c(1,2,3,4,6,7,9,10,11,12), labels=c("1:1","2:1","1:2","0:2","JE","GO","1:1","2:1","1:2","0:2"), las=1)
box()

#add overall mean for JE and GO
#model predicted means

lsm<-emmeans(MCb,~Population+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Population+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(6,7)
points(position,plot$response,pch=4,cex=2)

#add mean for JE incompetition with GO
M2Cb<-lm(sqrt(Competitive_ability_JE) ~  Mixture_Code1 +Initial_biomass,data=dat_2a)

lsm<-emmeans(M2Cb,~Mixture_Code1+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Mixture_Code1+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(1,2,3,4)
points(position,plot$response,pch=4,cex=2)

#add mean for GO incompetition with JE
M3Cb<-lm(sqrt(Competitive_ability_GO) ~  Mixture_Code1 +Initial_biomass,data=dat_3a)

lsm<-emmeans(M3Cb,~Mixture_Code1+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Mixture_Code1+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(9,10,11,12)
points(position,plot$response,pch=4,cex=2)
##

#RELATIVE COMPETITIVE ABILITY BASED ON MONOCULTURE

#Biomass was considered as measure of yield
dat_1c2<-read.table(file="Competitive_Combine_Mon_H2.txt", header=TRUE,sep="\t",dec = ".")

str(dat_1c2)
#Use lmer with RMEL=F
MC <- lmer(Competitive_ability  ~ Population+(1 | Mixture_series)+(1 | Incompetition)+Initial_biomass, REML=F, data = dat_1c2)
MCb <- lmer(sqrt(Competitive_ability)  ~ Population+(1 | Mixture_series)+(1 | Incompetition)+Initial_biomass, REML=F, data = dat_1c2)
MCc <- lmer(log(Competitive_ability)  ~ Population+(1 | Mixture_series)+(1 | Incompetition)+Initial_biomass, REML=F, data = dat_1c2)

#Check themodel(DHARMarequired)
sim1a<-simulateResiduals(fittedModel= MC,n=1000) #fit better
plot(sim1a, quantreg= FALSE)

sim1a<-simulateResiduals(fittedModel= MC,n=1000) 
plot(sim1a, quantreg= FALSE)

sim1a<-simulateResiduals(fittedModel= MC,n=1000) 
plot(sim1a, quantreg= FALSE)

summary(MC)
#untransformed is better
#Significant test
Anova(MC,type="III",test="Chisq")

#multiple comparaison
library(multcomp)
gltt_M1<-glht(MC,linfct=mcp(Population="Tukey"))
gltt_M1
summary(gltt_M1)
post<-cld(gltt_M1, level = 0.05, decreasing = F)
post
#

# use lm for JE incompetition GO  
str(dat_2c2)
dat_2c2<-read.table(file="Competitive_JE_Mon_H2.txt", header=TRUE,sep="\t",dec = ".")
Mixture_Code1<-as.factor(dat_2c2$Mixture_series)

# lm
M2C<-lm(Competitive_ability_JE ~  Mixture_Code1 +Initial_biomass,data=dat_2c2)
M2Cb<-lm(sqrt(Competitive_ability_JE) ~  Mixture_Code1 +Initial_biomass,data=dat_2c2)
M2Cc<-lm(log(Competitive_ability_JE) ~  Mixture_Code1 +Initial_biomass,data=dat_2c2)


# Test homogeneity of variance
leveneTest(Competitive_ability_JE~Mixture_Code,data=dat_2c2)
#Assumption is valid if pvalue >.05
# assumption is not valid, p =0.8


#Check normality with shapiro test
shapiro.test(dat_2c2$Competitive_ability_JE) 

#P VALUE > 0.05 the assumption is valid, p =0.0001

#LM
plot(M2C)  
hist(resid(M2C)) 

plot(M2Cb)  
hist(resid(M2Cb))

plot(M2Cc)  
hist(resid(M2Cc))

#untransformed fit better
Anova(M2C,type="III",test="F")
##
#multiple comparaison
Mixture_Code1<-as.factor(dat_2c2$Mixture_Code)
str(dat_2c2)
M2C<-lm(Competitive_ability_JE ~  Mixture_Code1 +Initial_biomass,data=dat_2c2)

gltt_M2a<-glht(M2C,linfct=mcp(Mixture_Code1="Tukey"))
gltt_M2a
summary(gltt_M2a)
post<-cld(gltt_M2a, level = 0.05, decreasing = F)
post
#COMPETITIVE INTENSITY OF GO incompetition JE
dat_3c2<-read.table(file="Competitive_GO_Mon_H2.txt", header=TRUE,sep="\t",dec = ".")
Mixture_Code1<-as.factor(dat_3c2$Mixture_series)
str(dat_3a)
M3C<-lm(Competitive_ability_GO  ~  Mixture_Code1 +Initial_biomass,data=dat_3c2)
M3Cb<-lm(sqrt(Competitive_ability_GO)  ~  Mixture_Code1 +Initial_biomass,data=dat_3c2)
M3Cc<-lm(log(Competitive_ability_GO)  ~  Mixture_Code1 +Initial_biomass,data=dat_3c2)



# Test homogeneity of variance
leveneTest(Competitive_ability_GO ~Mixture_Code,data=dat_3c2)
#Assumption is valid if pvalue >.05
# assumption is not valid, p =0.53

#Check normality with shapiro test
shapiro.test(dat_3c2$Competitive_ability_GO ) 

#P VALUE > 0.05 the assumption is not valid, p =0.0001


plot(M3C) # 
hist(resid(M3C)) #

plot(M3Cb) # 
hist(resid(M3Cb)) 
plot(M3Cc) # 
hist(resid(M3Cc)) 

summary(M3)

#sqrt transformed data fit better
Anova(M3Cb,type="III",test="F")
##

#multiple comparaison
Mixture_Code1<-as.factor(dat_3c2$Mixture_Code)
str(dat_3c2)
M3C<-lm(Competitive_ability_GO  ~  Mixture_Code1 +Initial_biomass,data=dat_3c2)

gltt_M3C<-glht(M3C,linfct=mcp(Mixture_Code1="Tukey"))
gltt_M3C
summary(gltt_M3C)
post<-cld(gltt_M3C, level = 0.05, decreasing = F)
post
#Data for Box_plot 
dat_4ac<-read.table(file="Competitive_Box_Mon_H2.txt", header=TRUE,sep="\t",dec = ".")
str(dat_4ac)

boxplot(Competitive_ability ~Mixture_Code,data=dat_4ac)


position<-c(1,2,3,5,6,8,9,10)
boxplot(Competitive_ability ~Mixture_Code,data=dat_4ac,axes=F,at=position ,border="gray40",xlab= "Population", ylab="Relative competitive intensity",las=2, 
        col=c("#0033FF","#0033FF","#0033FF","#0033FF","#FFFF00","#FFFF00","#FFFF00","#FFFF00"))
axis(2,at=seq(-3,1,1), seq(-3,1,1), las=1)
axis(1, at=c(1,2,3,5,6,8,9,10), labels=c("1:1","2:1","1:2","JE","GO","1:1","2:1","1:2"), las=1)
box()

#add overall mean for JE and GO
#model predicted means

lsm<-emmeans(MC,~Population+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Population+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(5,6)
points(position,plot$emmean,pch=4,cex=2)

#add mean for JE incompetition with GO
M2C<-lm(Competitive_ability_JE ~  Mixture_Code1 +Initial_biomass,data=dat_2c2)

lsm<-emmeans(M2C,~Mixture_Code1+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Mixture_Code1+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(1,2,3)
points(position,plot$emmean,pch=4,cex=2)

#add mean for GO incompetition with JE
M3Cb<-lm(sqrt(Competitive_ability_GO ) ~  Mixture_Code1 +Initial_biomass,data=dat_3c2)

lsm<-emmeans(M3Cb,~Mixture_Code1+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Mixture_Code1+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(8,9,10)
points(position,plot$response,pch=4,cex=2)
##

#ROOT BIOMASS FRACTION (Below biomass / Total biomass 1 and 2 harvest)


dat_1r<-read.table(file="Root_biomass_combine_H2.txt", header=TRUE,sep="\t",dec = ".")

str(dat_1r)
#Use lmer with RMEL=F
MC <- lmer(Root_biomass_fraction_.  ~ Population+(1 | Mixture_series)+Initial_biomass, REML=F, data = dat_1r)
MCb <- lmer(sqrt(Root_biomass_fraction_.)  ~ Population+(1 | Mixture_series)+Initial_biomass, REML=F, data = dat_1r)
MCc <- lmer(log(Root_biomass_fraction_.)  ~ Population+(1 | Mixture_series)+Initial_biomass, REML=F, data = dat_1r)

#Check themodel(DHARMarequired)
sim1a<-simulateResiduals(fittedModel= MC,n=1000) # 
plot(sim1a, quantreg= FALSE)

sim1a<-simulateResiduals(fittedModel= MCb,n=1000) 
plot(sim1a, quantreg= FALSE)

sim1a<-simulateResiduals(fittedModel= MCc,n=1000) #log fit better
plot(sim1a, quantreg= FALSE)

summary(MC)
#log transformed is better
#Significant test
Anova(MCc,type="III",test="Chisq")

gltt_M2a<-glht(MCc,linfct=mcp(Mixture_Code="Tukey"))
gltt_M2a
summary(gltt_M2a)
post<-cld(gltt_M2a, level = 0.05, decreasing = F)
post
# use lm for JE ONLY  
str(dat_2r)
dat_2r<-read.table(file="Root_JE_H2.txt", header=TRUE,sep="\t",dec = ".")
Mixture_Code1<-as.factor(dat_2r$Mixture_series)
str(dat_2r)
# lm
M2C<-lm(Root_biomass_fraction_. ~  Mixture_Code1 +Initial_biomass,data=dat_2r)
M2Cb<-lm(sqrt(Root_biomass_fraction_.) ~  Mixture_Code1 +Initial_biomass,data=dat_2r)
M2Cc<-lm(log(Root_biomass_fraction_.) ~  Mixture_Code1 +Initial_biomass,data=dat_2r)


# Test homogeneity of variance
leveneTest(Root_biomass_fraction_.~Mixture_Code,data=dat_2r)
#Assumption is valid if pvalue >.05
# assumption is not valid, p =0.04


#Check normality with shapiro test
shapiro.test(dat_2r$Root_biomass_fraction_.) 

#P VALUE > 0.05 the assumption is valid, p =0.003

#LM
plot(M2C)  
hist(resid(M2C)) 

plot(M2Cb)  
hist(resid(M2Cb))

plot(M2Cc)  
hist(resid(M2Cc))

#log transformed fit better
Anova(M2Cc,type="III",test="F")
##
#multiple comparaison
Mixture_Code1<-as.factor(dat_2r$Mixture_Code)
str(dat_2r)
M2Cc<-lm(log(Root_biomass_fraction_.) ~  Mixture_Code1 +Initial_biomass,data=dat_2r)

gltt_M2a<-glht(M2Cc,linfct=mcp(Mixture_Code1="Tukey"))
gltt_M2a
summary(gltt_M2a)
post<-cld(gltt_M2a, level = 0.05, decreasing = F)
post
#USE LM GO ONLY
dat_3r<-read.table(file="Root_GO_H2.txt", header=TRUE,sep="\t",dec = ".")
Mixture_Code1<-as.factor(dat_3r$Mixture_series)
str(dat_3r)
M3C<-lm(Root_biomass_fraction_.  ~  Mixture_Code1 +Initial_biomass,data=dat_3r)
M3Cb<-lm(sqrt(Root_biomass_fraction_.)  ~  Mixture_Code1 +Initial_biomass,data=dat_3r)
M3Cc<-lm(log(Root_biomass_fraction_.)  ~  Mixture_Code1 +Initial_biomass,data=dat_3r)



# Test homogeneity of variance
leveneTest(Root_biomass_fraction_. ~Mixture_Code,data=dat_3r)
#Assumption is valid if pvalue >.05
# assumption is not valid, p =0.57

#Check normality with shapiro test
shapiro.test(dat_3r$Root_biomass_fraction_. ) 

#P VALUE > 0.05 the assumption is valid, p =0.005


plot(M3C) # 
hist(resid(M3C)) #

plot(M3Cb) # 
hist(resid(M3Cb)) 
plot(M3Cc) # 
hist(resid(M3Cc)) #fit

summary(M3)

#log transformed data fit better
Anova(M3Cc,type="III",test="F")
##

#multiple comparaison
Mixture_Code1<-as.factor(dat_3r$Mixture_Code)
str(dat_3r)
M3Cc<-lm(log(Root_biomass_fraction_.)  ~  Mixture_Code1 +Initial_biomass,data=dat_3r)

gltt_M3C<-glht(M3Cc,linfct=mcp(Mixture_Code1="Tukey"))
gltt_M3C
summary(gltt_M3C)
post<-cld(gltt_M3C, level = 0.05, decreasing = F)
post
#Data for Box_plot 
dat_4r<-read.table(file="Root_Box_H2.txt", header=TRUE,sep="\t",dec = ".")
str(dat_4r)

boxplot(Root_biomass_fraction_. ~Mixture_Code,data=dat_4r)


position<-c(1,2,3.5,4.5,6,7)
boxplot(Root_biomass_fraction_. ~Mixture_Code,data=dat_4r,axes=F,at=position ,border="gray40",xlab= "Population", ylab="Root biomass fraction [%]",las=2, 
        col=c("#0033FF","#0033FF","#0033FF","#FFFF00","#FFFF00","#FFFF00"))
axis(2,at=seq(30,80,10), seq(30,80,10), las=1)
axis(1, at=c(1,2,3.5,4.5,6,7), labels=c("0:1","0:2","JE","GO","0:1","0:2"), las=1)
box()

#add overall mean for JE and GO
#model predicted means

lsm<-emmeans(MCc,~Population+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Population+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(3.5,4.5)
points(position,plot$response,pch=4,cex=2)

#add mean for JE incompetition with GO
M2Cc<-lm(log(Root_biomass_fraction_.) ~  Mixture_Code1 +Initial_biomass,data=dat_2r)

lsm<-emmeans(M2Cc,~Mixture_Code1+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Mixture_Code1+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(1,2)
points(position,plot$response,pch=4,cex=2)

#add mean for GO incompetition with JE
M3Cc<-lm(log(Root_biomass_fraction_. ) ~  Mixture_Code1 +Initial_biomass,data=dat_3r)

lsm<-emmeans(M3Cc,~Mixture_Code1+Initial_biomass)
pairs(lsm,adjust="tukey")
lsm

pan.emm <- emmeans(regrid(lsm), list(~Mixture_Code1+Initial_biomass))
pan.emm
plot<-as.data.frame(pan.emm)
plot
position<-c(6,7)
points(position,plot$response,pch=4,cex=2)
##
