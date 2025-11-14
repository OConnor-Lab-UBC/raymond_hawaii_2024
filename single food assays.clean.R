### Calla Hawaii Fish experiment No.1 : no choice feeding assays 

# 1. Load data and packages --------------------------------------------------

### load libraries
library(MuMIn)
library(nlme)
library(plyr)
library(tidyverse) 
library(broom)
library(reshape2)
library(gridExtra)
library(RColorBrewer)
library(viridis)
library(googlesheets4)
library(colorspace)
library(car)
library(emmeans)
library(sjPlot)
library(lme4)

## load data

data <- read_sheet('https://docs.google.com/spreadsheets/d/1kANQcC19OzdMkrYXXnrXRRlUR2ecax9yrbSCe9qfnaE/edit?gid=668485896#gid=668485896', sheet=2) 

View(data)

data$fish.weight.g <- as.numeric(data$fish.weight.g)

# make factors 
data$fish.present <- as.factor(data$fish.present)
data$fish.species <- as.factor(data$fish.species)
data$algal.latin <- as.factor(data$algal.latin)
data$algal.species <- as.factor(data$algal.species)

# Set reference variable to the one with the most observations ##
data$fish.present <- relevel(data$fish.present, ref = "0")
data$fish.species <- relevel(data$fish.species, ref = "Control")

# 2a. Create datafiles for each fish --------------------------------------

#Unicorn
BUdata <- data[(data$fish.species=="Bluespine Unicornfish" | data$fish.species=="Control"),]
#View(BUdata)

#Sailfin
SFdata <- data[(data$fish.species=="Sailfin tang" | data$fish.species=="Control"),]
#View(SFdata)

#Convict
CTdata <- data[(data$fish.species=="Convict Tang" | data$fish.species=="Control"),]
#View(CTdata)

#Yellowfin 
YFdata <- data[(data$fish.species=="Yellowfin Surgeonfish" | data$fish.species=="Control"),]
#View(YFdata)

#Bullethead
BHdata <- data[(data$fish.species=="Bullethead Parrotfish" | data$fish.species=="Control"),]
#View(BHdata)

#Palenose
PNdata <- data[(data$fish.species=="Palenose Parrotfish" | data$fish.species=="Control"),]
#View(PNdata)

# 2b. Create datafiles for each algae --------------------------------------
#Gorilla ogo
gorilla.ogo <- subset(data, algal.species == "Gorilla ogo")

#Halimeda 
halimeda <- subset(data, algal.species == "Halimeda")

#Prickly seaweed
prickly <- subset(data, algal.species == "Prickly seaweed")

#Smothering seaweed
smothering <- subset(data, algal.species == "Smothering seaweed")


# 3. Model consumption for each fish ------------------------------------

############################################# Unicorn ####
summary(BUdata$delta.weight.g)
shapiro.test(BUdata$delta.weight.g)
hist(BUdata$delta.weight.g)

shapiro.test(exp(BUdata$delta.weight.g))
hist(exp(BUdata$delta.weight.g))
shapiro.test(inverse(BUdata$delta.weight.g))

#now looks normal - though still not normal shapiro, but p value less than exponent - go with exponent 

delta.weight.g.exp <- (exp(BUdata$delta.weight.g))

mod1a <- lme(delta.weight.g ~ algal.species * fish.present, random = list(fish.id = ~ 1, date = ~ 1), method ="ML", data = BUdata, na.action = na.exclude)

mod1b <- lme(delta.weight.g ~ algal.species + fish.present, random = list(fish.id = ~ 1, date = ~ 1), method ="ML", data = BUdata, na.action = na.exclude)

mod2a <- lme(delta.weight.g ~ algal.species * fish.present, random = ~ 1 | fish.id, method ="ML", data = BUdata, na.action = na.exclude)

mod2b <- lme(delta.weight.g ~ algal.species + fish.present, random = ~ 1 | fish.id, method ="ML", data = BUdata, na.action = na.exclude)

mod3a <- lme(delta.weight.g ~ algal.species * fish.present, random = ~ 1 | date, method ="ML", data = BUdata, na.action = na.exclude)

mod3b <- lme(delta.weight.g ~ algal.species + fish.present, random = ~ 1 | date, method ="ML", data = BUdata, na.action = na.exclude)

mod4 <- lme(delta.weight.g ~ algal.species, random = ~ 1 | fish.id, method ="ML", data = BUdata, na.action = na.exclude)

mod5 <- lme(delta.weight.g ~ fish.present, random = ~ 1 | fish.id, method ="ML", data = BUdata, na.action = na.exclude)

mod6 <- lme(delta.weight.g ~ 1, random = ~ 1 | fish.id, method ="ML", data = BUdata, na.action = na.exclude)

model.sel(mod1a, mod1b, mod2a, mod2b, mod3a, mod3b, mod4, mod5, mod6)

BU.models <- as.data.frame(model.sel(mod1a, mod1b, mod2a, mod2b, mod3a, mod3b, mod4, mod5, mod6))
# write.csv(BU.models, "model selection - BU E1.csv")

##Best is mod 3a - re-run as ML 
mod3a <- lmer(delta.weight.g ~ algal.species * fish.present + (1 | date),  data = BUdata, na.action = na.exclude)

Anova(mod3a)
summary(mod3a)
tab_model(mod3a)

## extraBU estimated means 
EMM.BU <- emmeans(mod3a, ~ algal.species * fish.present)

## export results 
BU.model <- capture.output(summary(mod3a))
# write.csv(BU.model, "BU best model E1..csv")

EMM.BU.df <- as.data.frame(EMM.BU)
# write.csv(EMM.BU.df, "EMM.BU.df.raw.csv")

### Simple pairwise comparisons...
pairs(EMM.BU, simple = "algal.species")    # compare treatments for fish presence
EMM.BU.PW <- pairs(EMM.BU, simple = "fish.present")     # compare doses for each algal species
EMM.BU.PW

plot(EMM.BU.PW)

# write.csv(as.data.frame(EMM.BU.PW), "Pairwise_means_BU.raw.csv")

### Plot only estimated means

ggplot(EMM.BU.df, aes(x = emmean, y = algal.species, color = as.factor(fish.present))) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = lower.CL, xmax = upper.CL), height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  labs(
    y = "Algal Species",
    x = "Estimated Weight Change (g)",
    color = "Fish Presence") +
  theme_minimal()

# ggsave("BU_emmeans.raw.png", device = "png", path = './figures/', width = 7, height = 4)

############################################# Sailfin ####

# test for normality 
shapiro.test(SFdata$delta.weight.g)
shapiro.test(exp(SFdata$delta.weight.g))
# 1.253 -05
shapiro.test(10^(SFdata$delta.weight.g))
# worse 

hist(SFdata$delta.weight.g)

# not normal
SFmax <- max(SFdata$delta.weight.g)

#reflect data 
SF.reflect <- SFdata$delta.weight.g - SFmax

SF.reflect.p <- abs(SF.reflect)

hist(1/(SF.reflect.p + 1))
shapiro.test(1/(SF.reflect.p + 1))
## still worse that exp

#now looks normal - though still not normal shapiro, but p value less than exponent - go with exponent 

delta.weight.g.exp <- (exp(SFdata$delta.weight.g))

mod1a <- lme(delta.weight.g.exp ~ algal.species * fish.present, random = list(fish.id = ~ 1, date = ~ 1), method ="ML", data = SFdata, na.action = na.exclude)

mod1b <- lme(delta.weight.g.exp ~ algal.species + fish.present, random = list(fish.id = ~ 1, date = ~ 1), method ="ML", data = SFdata, na.action = na.exclude)

mod2a <- lme(delta.weight.g.exp ~ algal.species * fish.present, random = ~ 1 | fish.id, method ="ML", data = SFdata, na.action = na.exclude)

mod2b <- lme(delta.weight.g.exp ~ algal.species + fish.present, random = ~ 1 | fish.id, method ="ML", data = SFdata, na.action = na.exclude)

mod3a <- lme(delta.weight.g.exp ~ algal.species * fish.present, random = ~ 1 | date, method ="ML", data = SFdata, na.action = na.exclude)

mod3b <- lme(delta.weight.g.exp ~ algal.species + fish.present, random = ~ 1 | date, method ="ML", data = SFdata, na.action = na.exclude)

mod4 <- lme(delta.weight.g.exp ~ algal.species, random = ~ 1 | fish.id, method ="ML", data = SFdata, na.action = na.exclude)

mod5 <- lme(delta.weight.g.exp ~ fish.present, random = ~ 1 | fish.id, method ="ML", data = SFdata, na.action = na.exclude)

mod6 <- lme(delta.weight.g.exp ~ 1, random = ~ 1 | fish.id, method ="ML", data = SFdata, na.action = na.exclude)

model.sel(mod1a, mod1b, mod2a, mod2b, mod3a, mod3b, mod4, mod5, mod6)
SF.models <- as.data.frame(model.sel(mod1a, mod1b, mod2a, mod2b, mod3a, mod3b, mod4, mod5, mod6))
## write.csv(SF.models, "model selection - SF E1.csv")

##Best is mod 3a - re-run as REML 
mod3a <- lme(exp(delta.weight.g) ~ algal.species * fish.present, random = ~ 1 | date, method = "REML", data = SFdata, na.action = na.exclude)

Anova(mod3a)
summary(mod3a)
tab_model(mod3a)

## extract estimated means 
EMM.SF <- emmeans(mod3a, ~ algal.species * fish.present, type = "response")
EMM.SF

## export results 
SF.model <- capture.output(summary(mod3a))
# write.csv(SF.model, "SF best model E1.csv")

EMM.SF.df <- as.data.frame(EMM.SF)
# write.csv(EMM.SF.df, "EMM.SF.df.csv")

### Plot only estimated means
emm_sf_df <- as.data.frame(EMM.SF) 

ggplot(emm_sf_df, aes(x = response, y = algal.species, color = as.factor(fish.present))) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = lower.CL, xmax = upper.CL), height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  labs(
    y = "Algal Species",
    x = "Estimated Weight Change (g)",
    color = "Fish Presence") +
  theme_minimal()

# ggsave("SF_emmeans.png", device = "png", path = './figures/', width = 7, height = 4)

### Simple pairwise comparisons...
pairs(EMM.SF, simple = "algal.species")    # compare treatments for fish presence
SF.pairs <- as.data.frame(pairs(EMM.SF, simple = "fish.present"))     # compare doses for each algal species

plot(SF.pairs)
# write.csv(SF.pairs,"SF.pairwise.emmeans.csv")

#############################################  Manini  ####
shapiro.test(CTdata$delta.weight.g)
hist(CTdata$delta.weight.g)

shapiro.test(exp(CTdata$delta.weight.g))
hist(exp(CTdata$delta.weight.g))

#now looks normal - though still not normal shapiro, but p value less than exponent - go with exponent 

delta.weight.g.exp <- (exp(CTdata$delta.weight.g))

mod1a <- lme(delta.weight.g.exp ~ algal.species * fish.present, random = list(fish.id = ~ 1, date = ~ 1), method ="ML", data = CTdata, na.action = na.exclude)

mod1b <- lme(delta.weight.g.exp ~ algal.species + fish.present, random = list(fish.id = ~ 1, date = ~ 1), method ="ML", data = CTdata, na.action = na.exclude)

mod2a <- lme(delta.weight.g.exp ~ algal.species * fish.present, random = ~ 1 | fish.id, method ="ML", data = CTdata, na.action = na.exclude)

mod2b <- lme(delta.weight.g.exp ~ algal.species + fish.present, random = ~ 1 | fish.id, method ="ML", data = CTdata, na.action = na.exclude)

mod3a <- lme(delta.weight.g.exp ~ algal.species * fish.present, random = ~ 1 | date, method ="ML", data = CTdata, na.action = na.exclude)

mod3b <- lme(delta.weight.g.exp ~ algal.species + fish.present, random = ~ 1 | date, method ="ML", data = CTdata, na.action = na.exclude)

mod4 <- lme(delta.weight.g.exp ~ algal.species, random = ~ 1 | fish.id, method ="ML", data = CTdata, na.action = na.exclude)

mod5 <- lme(delta.weight.g.exp ~ fish.present, random = ~ 1 | fish.id, method ="ML", data = CTdata, na.action = na.exclude)

mod6 <- lme(delta.weight.g.exp ~ 1, random = ~ 1 | fish.id, method ="ML", data = CTdata, na.action = na.exclude)

model.sel(mod1a, mod1b, mod2a, mod2b, mod3a, mod3b, mod4, mod5, mod6)

CT.models <- as.data.frame(model.sel(mod1a, mod1b, mod2a, mod2b, mod3a, mod3b, mod4, mod5, mod6))
# write.csv(CT.models, "model selection - CT E1.csv")

##Best is mod 3a - re-run as ML 
mod3a <- lme(exp(delta.weight.g) ~ algal.species * fish.present, random = ~ 1 | date, method = "REML", data = CTdata, na.action = na.exclude)

Anova(mod3a)
summary(mod3a)
tab_model(mod3a)

## extract estimated means 
EMM.CT <- emmeans(mod3a, ~ algal.species * fish.present, type = "response")
EMM.CT

## export results 
CT.model <- capture.output(summary(mod3a))
# write.csv(CT.model, "CT best model E1.csv")

EMM.CT.df <- as.data.frame(EMM.CT)
# write.csv(EMM.CT.df, "EMM.CT.df.csv")

### Plot only estimated means

ggplot(EMM.CT.df, aes(x = response, y = algal.species, color = as.factor(fish.present))) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = lower.CL, xmax = upper.CL), height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  labs(
    y = "Algal Species",
    x = "Estimated Weight Change (g)",
    color = "Fish Presence") +
  theme_minimal()

# ggsave("CT_emmeans.png", device = "png", path = './figures/', width = 7, height = 4)

### Simple pairwise comparisons...
pairs(EMM.CT, simple = "algal.species")    # compare treats for fish presence
CT.pairs <- as.data.frame(pairs(EMM.CT, simple = "fish.present"))     # compare doses for each algal species

# write.csv(CT.pairs,"CT.pairwise.emmeans.csv")

############################################# Yellowfin ####
shapiro.test(YFdata$delta.weight.g)
hist(YFdata$delta.weight.g)

shapiro.test(exp(YFdata$delta.weight.g))
hist(log(YFdata$delta.weight.g))

#now looks normal - though still not normal shapiro, but p value less than exponent - go with exponent 

delta.weight.g.exp <- (exp(YFdata$delta.weight.g))

mod1a <- lme(delta.weight.g.exp ~ algal.species * fish.present, random = list(fish.id = ~ 1, date = ~ 1), method ="ML", data = YFdata, na.action = na.exclude)

mod1b <- lme(delta.weight.g.exp ~ algal.species + fish.present, random = list(fish.id = ~ 1, date = ~ 1), method ="ML", data = YFdata, na.action = na.exclude)

mod2a <- lme(delta.weight.g.exp ~ algal.species * fish.present, random = ~ 1 | fish.id, method ="ML", data = YFdata, na.action = na.exclude)

mod2b <- lme(delta.weight.g.exp ~ algal.species + fish.present, random = ~ 1 | fish.id, method ="ML", data = YFdata, na.action = na.exclude)

mod3a <- lme(delta.weight.g.exp ~ algal.species * fish.present, random = ~ 1 | date, method ="ML", data = YFdata, na.action = na.exclude)

mod3b <- lme(delta.weight.g.exp ~ algal.species + fish.present, random = ~ 1 | date, method ="ML", data = YFdata, na.action = na.exclude)

mod4 <- lme(delta.weight.g.exp ~ algal.species, random = ~ 1 | fish.id, method ="ML", data = YFdata, na.action = na.exclude)

mod5 <- lme(delta.weight.g.exp ~ fish.present, random = ~ 1 | fish.id, method ="ML", data = YFdata, na.action = na.exclude)

mod6 <- lme(delta.weight.g.exp ~ 1, random = ~ 1 | fish.id, method ="ML", data = YFdata, na.action = na.exclude)

model.sel(mod1a, mod1b, mod2a, mod2b, mod3a, mod3b, mod4, mod5, mod6)

YF.models <- as.data.frame(model.sel(mod1a, mod1b, mod2a, mod2b, mod3a, mod3b, mod4, mod5, mod6))
# write.csv(YF.models, "model seleYFion - YF E1.csv")

##Best is mod 3a - re-run as ML 
mod3a <- lme(exp(delta.weight.g) ~ algal.species * fish.present, random = ~ 1 | date, method = "REML", data = YFdata, na.action = na.exclude)

Anova(mod3a)
summary(mod3a)
tab_model(mod3a)

## extraYF estimated means 
EMM.YF <- emmeans(mod3a, ~ algal.species * fish.present, type = "response")
EMM.YF

## export results 
YF.model <- capture.output(summary(mod3a))
# # write.csv(YF.model, "YF best model E1.csv")

EMM.YF.df <- as.data.frame(EMM.YF)
# write.csv(EMM.YF.df, "EMM.YF.df.csv")

### Plot only estimated means

ggplot(EMM.YF.df, aes(x = response, y = algal.species, color = as.factor(fish.present))) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = lower.CL, xmax = upper.CL), height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  labs(
    y = "Algal Species",
    x = "Estimated Weight Change (g)",
    color = "Fish Presence") +
  theme_minimal()

ggsave("YF_emmeans.png", device = "png", path = './figures/', width = 7, height = 4)

### Simple pairwise comparisons...
pairs(EMM.YF, simple = "algal.species")    # compare treats for fish presence
YF.pairs <- as.data.frame(pairs(EMM.YF, simple = "fish.present"))     # compare doses for each algal species

# write.csv(YF.pairs,"YF.pairwise.emmeans.csv")

############################################# Bullethead ####
shapiro.test(BHdata$delta.weight.g)
hist(BHdata$delta.weight.g)

shapiro.test(exp(BHdata$delta.weight.g))
hist(exp(BHdata$delta.weight.g))

#now looks normal - though still not normal shapiro, but p value less than exponent - go with exponent 

delta.weight.g.exp <- (exp(BHdata$delta.weight.g))

mod1a <- lme(delta.weight.g.exp ~ algal.species * fish.present, random = list(fish.id = ~ 1, date = ~ 1), method ="ML", data = BHdata, na.action = na.exclude)

mod1b <- lme(delta.weight.g.exp ~ algal.species + fish.present, random = list(fish.id = ~ 1, date = ~ 1), method ="ML", data = BHdata, na.action = na.exclude)

mod2a <- lme(delta.weight.g.exp ~ algal.species * fish.present, random = ~ 1 | fish.id, method ="ML", data = BHdata, na.action = na.exclude)

mod2b <- lme(delta.weight.g.exp ~ algal.species + fish.present, random = ~ 1 | fish.id, method ="ML", data = BHdata, na.action = na.exclude)

mod3a <- lme(delta.weight.g.exp ~ algal.species * fish.present, random = ~ 1 | date, method ="ML", data = BHdata, na.action = na.exclude)

mod3b <- lme(delta.weight.g.exp ~ algal.species + fish.present, random = ~ 1 | date, method ="ML", data = BHdata, na.action = na.exclude)

mod4 <- lme(delta.weight.g.exp ~ algal.species, random = ~ 1 | fish.id, method ="ML", data = BHdata, na.action = na.exclude)

mod5 <- lme(delta.weight.g.exp ~ fish.present, random = ~ 1 | fish.id, method ="ML", data = BHdata, na.action = na.exclude)

mod6 <- lme(delta.weight.g.exp ~ 1, random = ~ 1 | fish.id, method ="ML", data = BHdata, na.action = na.exclude)

model.sel(mod1a, mod1b, mod2a, mod2b, mod3a, mod3b, mod4, mod5, mod6)

BH.models <- as.data.frame(model.sel(mod1a, mod1b, mod2a, mod2b, mod3a, mod3b, mod4, mod5, mod6))
# write.csv(BH.models, "model selection - BH E1.csv")

##Best is mod 3a - re-run as REML 
mod3a <- lme(exp(delta.weight.g) ~ algal.species * fish.present, random = ~ 1 | date, method = "REML", data = BHdata, na.action = na.exclude)

Anova(mod3a)
summary(mod3a)
tab_model(mod3a)

## extraBH estimated means 
EMM.BH <- emmeans(mod3a, ~ algal.species * fish.present, type = "response")
EMM.BH

## export results 
BH.model <- capture.output(summary(mod3a))
# write.csv(BH.model, "BH best model E1.csv")

EMM.BH.df <- as.data.frame(EMM.BH)
# write.csv(EMM.BH.df, "EMM.BH.df.csv")

### Plot only estimated means

ggplot(EMM.BH.df, aes(x = response, y = algal.species, color = as.factor(fish.present))) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = lower.CL, xmax = upper.CL), height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  labs(
    y = "Algal Species",
    x = "Estimated Weight Change (g)",
    color = "Fish Presence") +
  theme_minimal()

ggsave("BH_emmeans.png", device = "png", path = './figures/', width = 7, height = 4)

### Simple pairwise comparisons...
pairs(EMM.BH, simple = "algal.species")    # compare treats for fish presence
BH.pairs <- as.data.frame(pairs(EMM.BH, simple = "fish.present"))     # compare doses for each algal species

# write.csv(BH.pairs,"BH.pairwise.emmeans.csv")

### The default is to apply a separate Tukey adjustment to the P values in each by group (so if each group has just 2 means, no adjustment at all is applied). If you want to adjust the whole family combined, you need to undo the by variable and specify the desired adjustment (which can't be Tukey because that method is invalid when you have more than one set of pairwise comparisons.) For example

test(pairs(EMM.BH, by = "fish.present"), by = NULL, adjust = "mvt")

############################################# Palenose ####
shapiro.test(PNdata$delta.weight.g)
hist(PNdata$delta.weight.g)

shapiro.test(exp(PNdata$delta.weight.g))
hist(exp(PNdata$delta.weight.g))

#now looks normal - though still not normal shapiro, but p value less than exponent - go with exponent 

delta.weight.g.exp <- (exp(PNdata$delta.weight.g))

mod1a <- lme(delta.weight.g.exp ~ algal.species * fish.present, random = list(fish.id = ~ 1, date = ~ 1), method ="ML", data = PNdata, na.action = na.exclude)

mod1b <- lme(delta.weight.g.exp ~ algal.species + fish.present, random = list(fish.id = ~ 1, date = ~ 1), method ="ML", data = PNdata, na.action = na.exclude)

mod2a <- lme(delta.weight.g.exp ~ algal.species * fish.present, random = ~ 1 | fish.id, method ="ML", data = PNdata, na.action = na.exclude)

mod2b <- lme(delta.weight.g.exp ~ algal.species + fish.present, random = ~ 1 | fish.id, method ="ML", data = PNdata, na.action = na.exclude)

mod3a <- lme(delta.weight.g.exp ~ algal.species * fish.present, random = ~ 1 | date, method ="ML", data = PNdata, na.action = na.exclude)

mod3b <- lme(delta.weight.g.exp ~ algal.species + fish.present, random = ~ 1 | date, method ="ML", data = PNdata, na.action = na.exclude)

mod4 <- lme(delta.weight.g.exp ~ algal.species, random = ~ 1 | fish.id, method ="ML", data = PNdata, na.action = na.exclude)

mod5 <- lme(delta.weight.g.exp ~ fish.present, random = ~ 1 | fish.id, method ="ML", data = PNdata, na.action = na.exclude)

mod6 <- lme(delta.weight.g.exp ~ 1, random = ~ 1 | fish.id, method ="ML", data = PNdata, na.action = na.exclude)

model.sel(mod1a, mod1b, mod2a, mod2b, mod3a, mod3b, mod4, mod5, mod6)

PN.models <- as.data.frame(model.sel(mod1a, mod1b, mod2a, mod2b, mod3a, mod3b, mod4, mod5, mod6))
# write.csv(PN.models, "model selection - PN E1.csv")

##Best is mod 3b - re-run as ML 
mod3a <- lme(exp(delta.weight.g) ~ algal.species * fish.present, random = ~ 1 | date, method = "REML", data = PNdata, na.action = na.exclude)

Anova(mod3a)
summary(mod3a)
tab_model(mod3a)

## extraPN estimated means 
EMM.PN <- emmeans(mod3a, ~ algal.species * fish.present, type = "response")
EMM.PN

## export results 
PN.model <- capture.output(summary(mod3a))
# # write.csv(PN.model, "PN best model E1.csv")

EMM.PN.df <- as.data.frame(EMM.PN)
# write.csv(EMM.PN.df, "EMM.PN.df.csv")

### Plot only estimated means

ggplot(EMM.PN.df, aes(x = response, y = algal.species, color = as.factor(fish.present))) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = lower.CL, xmax = upper.CL), height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  labs(
    y = "Algal Species",
    x = "Estimated Weight Change (g)",
    color = "Fish Presence") +
  theme_minimal()

# ggsave("PN_emmeans.png", device = "png", path = './figures/', width = 7, height = 4)

### Simple pairwise comparisons...
pairs(EMM.PN, simple = "algal.species")    # compare treats for fish presence
PN.pairs <- as.data.frame(pairs(EMM.PN, simple = "fish.present"))     # compare doses for each algal species

# write.csv(PN.pairs,"PN.pairwise.emmeans.csv")

# 3. Plot results --------------------------------------------------
############################################# 3.1 Figure one ------------------------------------------------------

### Create labels 
#xlabs <- c("Gorilla ogo", "", "Halimeda", "", "Prickly seaweed", "", "Smothering seaweed", "")

### Plot consumption by algal species 

### Merge all estimates means 
  
  ## re-name column in BU data
  EMM.BU.df <- as.data.frame(EMM.BU) |>
    dplyr::rename(response = emmean)
  
# Combine data frames and filter for relevant species and fish presence
EMM_combined <- bind_rows(
  EMM.PN.df %>% mutate(fish.species = "Palenose Parrotfish"),
  EMM.BH.df %>% mutate(fish.species = "Bullethead Parrotfish"),
  EMM.BU.df %>% mutate(fish.species = "Bluespine Unicornfish"),
  EMM.CT.df %>% mutate(fish.species = "Convict Tang"),
  EMM.SF.df %>% mutate(fish.species = "Sailfin tang"),
  EMM.YF.df %>% mutate(fish.species = "Yellowfin Surgeonfish")) 

# write.csv(EMM_combined, "E1.estimated.means.summary.csv")

#EMM_combined <- read.csv("E1.estimated.means.summary.edit.csv")

EMM_combined_pres <- EMM_combined %>%
  filter(fish.present == 1)
EMM_combined_abs<- EMM_combined %>%
  filter(fish.present == 0) %>%
  mutate(fish_species=c("Control"))

# Remove leading "0." or "1." from algal.species in EMM_combined_pres and EMM_combined_abs

#EMM_combined_pres$algal.species <- sub("^[01]\\.", "", EMM_combined_pres$algal.species)
#EMM_combined_abs$algal.species <- sub("^[01]\\.", "", EMM_combined_abs$algal.species)


#EMM_combined_go <- EMM_combined %>% filter(algal.species == "1.Gorilla ogo")
#EMM_combined_ps <- EMM_combined %>% filter(algal.species == "1.Prickly seaweed")
#EMM_combined_ss <- EMM_combined %>% filter(algal.species == "1.Smothering seaweed")
#EMM_combined_h <- EMM_combined %>% filter(algal.species == "1.Halimeda")


#### FIGURE ONE 

library(rphylopic)

algae.labs <- c( "Halimeda"  = "italic(H.~discoidea)", 
                 "Prickly seaweed" = "italic(A.~spicifera)", 
                 "Smothering seaweed" = "italic(E.~denticulatum)", 
                 "Gorilla ogo" = "italic(G.~salicornia)")

fish.labs <- c("Convict Tang" = expression(italic("A.triostegus")), 
               "Palenose Parrotfish" = expression(italic("S.psittacus")),
               "Bullethead Parrotfish" = expression(italic("C.spilurus")), 
               "Bluespine Unicornfish" = expression(italic("N.unicornis")), 
               "Sailfin tang" = expression(italic("Z.velifer")), 
               "Yellowfin Surgeonfish" = expression(italic("A.xanthopterus")),
               "Control" = "Control")
ggplot() +
  geom_jitter(data = data, aes(y = delta.weight.g, x = fish.species, color = delta.weight.g > 0), alpha = 0.6, width = 0.2, size = 2) +
  facet_wrap(~ algal.species,
             labeller = labeller(
               algal.species = as_labeller(algae.labs, label_parsed)), ncol = 2,) +
  #scale_colour_manual(values=c("red", "blue")) +
  geom_pointrange(data = EMM_combined_pres, aes(y = response, x = fish.species, ymin = lower.CL, ymax = upper.CL, fill = response>0), 
                shape = 21, color = "black") +
  geom_pointrange(data=EMM_combined_abs, aes(y=response, x=fish_species, ymin=lower.CL, ymax=upper.CL, fill = response>0),
                  shape = 21, color = "black") +
  xlab("Fish species") +
  scale_x_discrete(labels = fish.labs) + 
  ylab("Change in weight (g)") +
  theme_minimal() +
  geom_hline(yintercept=0, linetype="dashed") +
  coord_flip(ylim = c(-3, 1.0)) +
  theme(legend.position = "none") +
  add_phylopic(uuid = "56957948-cd3e-4d28-8c60-b6616a782b54", #bluespine unicornfish
               x=2, y=0.85, height = 0.75, alpha=1, fill = "black") +
  add_phylopic(name = "Chlorurus spilurus", #Bullethead Parrotfish
               x=3, y=0.85, height = 0.6, alpha=1, fill = "black", horizontal =  TRUE) +
  add_phylopic(name = "Acanthurus triostegus", #Convict Tang
               x=4, y=0.85, height = 1.2, alpha=1, horizontal = TRUE, fill = "black") +
  add_phylopic(name = "Scarus psittacus", #Palenose Parrotfish
               x=5, y=0.85, height = 1.1, alpha=1, horizontal = TRUE, fill = "black") +
  add_phylopic(name = "Zebrasoma velifer", #Sailfin tang
               x=6, y=0.85, height = 1.1, alpha=1, horizontal = TRUE,  fill = "black") +
  add_phylopic(uuid = "6ddba087-d2fb-429e-b87e-4935f8f222b1", #Yellowfin Surgeonfish
               x=7, y=0.85, height = 0.75, alpha=1, fill = "black")

ggsave("allspecies_e1.jitter.Nov14.png", device = "png", path = './figures/', width = 7, height = 4)

############################################# 3.2 Extra figures - consumption by fish species  #####
#Unicorn 
BUdata <- arrange(BUdata, algal.species, fish.present)
EMM.BU.df$algal.species <- interaction(EMM.BU.df$fish.present, EMM.BU.df$algal.species)
EMM.BU.df <- arrange(EMM.BU.df, algal.species, fish.present)

BU.plot <- ggplot(BUdata, aes(x = interaction(fish.present, algal.species), algal.species, y = delta.weight.g, color = factor(fish.present))) +
  geom_jitter(width = 0.3) +
  theme_bw() +
  xlab("Algae Species") +
  ylab("Algae Consumed (g)") +
  scale_x_discrete(labels = xlabs) +
  scale_color_manual(values = c("azure3", "lightseagreen"))+
  labs(title = "Bluespine Unicorn fish", size = 3) +
  
  # add estimated means 
  geom_point(data = EMM.BU.df, aes(x = algal.species, y = emmean, fill = factor(fish.present)),
             size = 4, shape = 21, color = "black") + # Triangles for estimated means
  scale_fill_manual(values = c("azure3", "lightseagreen")) + 
  
  geom_linerange(data = EMM.BU.df, aes(x = algal.species, ymin = lower.CL, ymax = upper.CL),
                 size = 0.5, inherit.aes = FALSE) # Error bars for confidence intervals

BU.plot
ggsave("BU_e1.jitter.png", device = "png", path = './figures/', width = 7, height = 4)

#Sailfin
SFdata <- arrange(SFdata, algal.species, fish.present)
EMM.SF.df$algal.species <- interaction(EMM.SF.df$fish.present, EMM.SF.df$algal.species)
EMM.SF.df <- arrange(EMM.SF.df, algal.species, fish.present)

SF.plot <- ggplot(SFdata, aes(x = interaction(fish.present, algal.species), algal.species, y = delta.weight.g, color = factor(fish.present))) +
  geom_jitter(width = 0.3) +
  theme_bw() +
  xlab("Algae Species") +
  ylab("Algae Consumed (g)") +
  scale_x_discrete(labels = xlabs) +
  scale_color_manual(values = c("azure3", "lightseagreen"))+
  labs(title = "Sailfin Tang", size = 3) +
  
  # add estimated means 
  geom_point(data = EMM.SF.df, aes(x = algal.species, y = response, fill = factor(fish.present)),
             size = 4, shape = 21, color = "black") + # Triangles for estimated means
  scale_fill_manual(values = c("azure3", "lightseagreen")) + 
  
  geom_linerange(data = EMM.SF.df, aes(x = algal.species, ymin = lower.CL, ymax = upper.CL),
                 size = 0.5, inherit.aes = FALSE) # Error bars for confidence intervals

SF.plot
ggsave("SF_e1.jitter.png", device = "png", path = './figures/', width = 7, height = 4)

#Convict
CTdata <- arrange(CTdata, algal.species, fish.present)
EMM.CT.df$algal.species <- interaction(EMM.CT.df$fish.present, EMM.CT.df$algal.species)
EMM.CT.df <- arrange(EMM.CT.df, algal.species, fish.present)

CT.plot <- ggplot(CTdata, aes(x = interaction(fish.present, algal.species), algal.species, y = delta.weight.g, color = factor(fish.present))) +
  geom_jitter(width = 0.3) +
  theme_bw() +
  xlab("Algae Species") +
  ylab("Algae Consumed (g)") +
  scale_x_discrete(labels = xlabs) +
  scale_color_manual(values = c("azure3", "lightseagreen"))+
  labs(title = "Convict Tang", size = 3) +
  
  # add estimated means 
  geom_point(data = EMM.CT.df, aes(x = algal.species, y = response, fill = factor(fish.present)),
             size = 4, shape = 21, color = "black") + # Triangles for estimated means
  scale_fill_manual(values = c("azure3", "lightseagreen")) + 
  
  geom_linerange(data = EMM.CT.df, aes(x = algal.species, ymin = lower.CL, ymax = upper.CL),
                 size = 0.5, inherit.aes = FALSE) # Error bars for confidence intervals

CT.plot 
ggsave("CT_e1.jitter.png", device = "png", path = './figures/', width = 7, height = 4)

#Yellowfin 
YFdata <- arrange(YFdata, algal.species, fish.present)
EMM.YF.df$algal.species <- interaction(EMM.YF.df$fish.present, EMM.YF.df$algal.species)
EMM.YF.df <- arrange(EMM.YF.df, algal.species, fish.present)

YF.plot <- ggplot(YFdata, aes(x = interaction(fish.present, algal.species), algal.species, y = delta.weight.g, color = factor(fish.present))) +
  geom_jitter(width = 0.3) +
  theme_bw() +
  xlab("Algae Species") +
  ylab("Algae Consumed (g)") +
  scale_x_discrete(labels = xlabs) +
  scale_color_manual(values = c("azure3", "lightseagreen"))+
  labs(title = "Yellowfin surgeonfish", size = 3) +
  
  # add estimated means 
  geom_point(data = EMM.YF.df, aes(x = algal.species, y = response, fill = factor(fish.present)),
             size = 4, shape = 21, color = "black") + # Triangles for estimated means
  scale_fill_manual(values = c("azure3", "lightseagreen")) + 
  
  geom_linerange(data = EMM.YF.df, aes(x = algal.species, ymin = lower.CL, ymax = upper.CL),
                 size = 0.5, inherit.aes = FALSE) # Error bars for confidence intervals

YF.plot
ggsave("YF_e1.jitter.png", device = "png", path = './figures/', width = 7, height = 4)

#Bullethead
BHdata <- arrange(BHdata, algal.species, fish.present)
EMM.BH.df$algal.species <- interaction(EMM.BH.df$fish.present, EMM.BH.df$algal.species)
EMM.BH.df <- arrange(EMM.BH.df, algal.species, fish.present)

BH.plot <- ggplot(BHdata, aes(x = interaction(fish.present, algal.species), algal.species, y = delta.weight.g, color = factor(fish.present))) +
  geom_jitter(width = 0.3) +
  theme_bw() +
  xlab("Algae Species") +
  ylab("Algae Consumed (g)") +
  scale_x_discrete(labels = xlabs) +
  scale_color_manual(values = c("azure3", "lightseagreen"))+
  labs(title = "Bullethead parrot", size = 3) +
  
  # add estimated means 
  geom_point(data = EMM.BH.df, aes(x = algal.species, y = response, fill = factor(fish.present)),
             size = 4, shape = 21, color = "black") + # Triangles for estimated means
  scale_fill_manual(values = c("azure3", "lightseagreen")) + 
  
  geom_linerange(data = EMM.BH.df, aes(x = algal.species, ymin = lower.CL, ymax = upper.CL),
                 size = 0.5, inherit.aes = FALSE) # Error bars for confidence intervals

BH.plot
ggsave("BH_e1.jitter.png", device = "png", path = './figures/', width = 7, height = 4)

#Palenose
PNdata <- arrange(PNdata, algal.species, fish.present)
EMM.PN.df$algal.species <- interaction(EMM.PN.df$fish.present, EMM.PN.df$algal.species)
EMM.PN.df <- arrange(EMM.PN.df, algal.species, fish.present)

PN.plot <- ggplot(PNdata, aes(x = interaction(fish.present, algal.species), algal.species, y = delta.weight.g, color = factor(fish.present))) +
  geom_jitter(width = 0.3) +
  theme_bw() +
  xlab("Algae Species") +
  ylab("Algae Consumed (g)") +
  scale_x_discrete(labels = xlabs) +
  scale_color_manual(values = c("azure3", "lightseagreen"))+
  labs(title = "Palenose parrot fish", size = 3) +
  
  # add estimated means 
  geom_point(data = EMM.PN.df, aes(x = algal.species, y = response, fill = factor(fish.present)),
             size = 4, shape = 21, color = "black") + # Triangles for estimated means
  scale_fill_manual(values = c("azure3", "lightseagreen")) + 
  
  geom_linerange(data = EMM.PN.df, aes(x = algal.species, ymin = lower.CL, ymax = upper.CL),
                 size = 0.5, inherit.aes = FALSE) # Error bars for confidence intervals

PN.plot
ggsave("PN_e1.jitter.png", device = "png", path = './figures/', width = 7, height = 4)

### bubble plot of change in algae 
## Create labels for algae species 
algae_labels <- c("Gorilla ogo" = expression(italic("G.salicornia")),
                  "Halimeda" = expression(italic("H.discoidea")),
                  "Prickly seaweed" = expression(italic("A.spicifera")),
                  "Smothering seaweed" = expression(italic("E.denticulatum")))

# 4. Correlations --------------------------------------------------
fish_data <- data[data$fish.present != "0", ]
control_data <- data[data$fish.present != "1", ]

## correlation between fish weight and algae consumed 
plot(data$delta.weight.g ~ data$fish.weight.g, col = group_col,
     xlab = "Fish weight (g)",
     ylab = "Change in algae weight (g)",
     las = 1,
     pch = 16,
     xlim = c(0,100))

lm_model <- lm(delta.weight.g ~ fish.weight.g, data = fish_data)
abline(lm_model, col = "black", lwd = 2)  


shapiro.test(data$delta.weight.g)
shapiro.test(data$fish.weight.g)
## Both non-normal

cor.test(fish_data$delta.weight.g, fish_data$fish.weight.g, method = "kendall")

## signifigant but small tau (-0.12)

ggplot(fish_data, aes(x = fish.weight.g, y = delta.weight.g, color = fish.latin)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~ fish.latin)

## Correlation between algae eaten and length 
plot(data$delta.weight.g ~ data$fish.length.cm,
     xlab = "Fish length (cm)",
     ylab = "Change in algae weight (g)",
     las = 1,
     pch = 16,
     xlim = c(0,20))

lm_model <- lm(delta.weight.g ~ fish.length.cm, data = data)
abline(lm_model, col = "black", lwd = 2)  

shapiro.test(data$fish.length.cm)
#normal 
cor.test(data$delta.weight.g, data$fish.length.cm)

### not signifignat 

## weight differences between groups 
shapiro.test(log(data$fish.weight.g))

aov <- aov((log(data$fish.weight.g) ~ data$fish.species))

summary(aov)

TukeyHSD(aov)