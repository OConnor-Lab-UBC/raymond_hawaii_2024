### 1. Load data and packages --------------------------------------------------

library(googlesheets4)

data <- read_sheet('https://docs.google.com/spreadsheets/d/1kANQcC19OzdMkrYXXnrXRRlUR2ecax9yrbSCe9qfnaE/edit?gid=300783329#gid=300783329', sheet=3) 

# Set reference variable to the one with the most observations ##
data$treatment <- as.factor(data$treatment)
data$treatment <- relevel(data$treatment, ref = "control")
data$treatment <- factor(data$treatment, levels = c("control", "3 species", "6 species", "manini", "psittacus", "spilurus", "unicorn", "velifer", "xanthopterus"))

# make factors 
data$diversity <- as.factor(data$diversity)
data$fish.presence<- as.factor(data$fish.presence)

# create a new variable for monocultures 
data$grouped_treatment <- ifelse(data$treatment %in% c("manini", "psittacus", "spilurus", "unicorn", "velifer", "xanthopterus"), 
                                 "monocultures", 
                                 as.character(data$treatment))

data$grouped_treatment <- factor(data$grouped_treatment, levels = c("control", "3 species", "6 species", "monocultures"))

data$grouped_treatment <- relevel(data$grouped_treatment, ref = "control")

### load libraries ###
library(MuMIn)
library(nlme)
library(plyr)
library(tidyverse) 
library(broom)
library(reshape2)
library(gridExtra)
library(RColorBrewer)
library(viridis)
library(colorspace)
library(shinyjs)
library(vioplot)
library(emmeans)
library(sjPlot)
library(blme)
library(visreg)
library(car)


### 2.  Model response for all treatments ####

# test for normality 
shapiro.test(data$delta.weight.total.s)
# p = 0.101
leveneTest(delta.weight.total.s ~ treatment, data = data)
# p = 0.3128

# mod1 <- lmer(delta.weight.total.s ~ treatment + (1 | date),  data = data)

## using blmer because estimates of random effect went to 0 - blmer incorporates a penalty on the variance in the likelihood that reduces the chances that estimates will go to zero.
mod1 <- blmer(delta.weight.total.s ~ treatment + (1 | date),  data = data)

visreg(mod1, xvar = "treatment", whitespace = 0.4, 
       points.par = list(cex = 1.1, col = "red"))

shapiro.test(residuals(mod1))
# p = 0.7382

tab_model(mod1)

## Estimate means using emmeans 
EMM.BEF <- emmeans(mod1, ~ treatment)
EMM.BEF

EMM.BEF.df <- as.data.frame(EMM.BEF)

# write.csv(EMM.BEF.df, "BEF.estimated.means.met.scale.Jan3.csv")

## plot means with comparison arrows
plot(EMM.BEF, comparisons = TRUE)

## plot means with no arrows 
plot(EMM.BEF,
     ylab = "Species treatment",
     xlab = "Estimated Weight Change (g)")

pwpp(EMM.BEF)

## complete pairwise comparisions 
pairs.full <- pairs(EMM.BEF)

plot(pairs.full)

# write.csv(pairs.full, "pairwise comparison. bef met scale.Jan3.csv")

    ### Visualize consumption patters across treatments (Figure 4) -----
    ## create custom axis labels
    xlabs <- c("Control", "3 species", "6 species", 
               expression(italic("A.triostegus")), 
               expression(italic("S.psittacus")),
               expression(italic("C.spilurus")), 
               expression(italic("N.unicornis")), 
               expression(italic("Z.velifer")), 
               expression(italic("A.xanthopterus")))
    
    ### Algal consumption by all treatments -  scatter plot with overlapping boxes ---
    
    ggplot(data, aes(x = treatment, y = delta.weight.total.s, color = treatment)) +
      # Scatter plot
      geom_jitter(width = 0.25, size = 3, alpha = 0.6) + 
      scale_x_discrete(labels = xlabs) + # Adjust transparency to make the points less overwhelming
      # Boxplot
      geom_boxplot(aes(group = treatment),
                   color = "black", 
                   fill = NA,
                   alpha = 0,  # Make box plot transparent
                   outlier.size = 2) +
      # Labels and themes
      labs(y = "Scaled change in weight (g)", x = "Fish Treatment") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotates x-axis labels for readability
      
      # Customize point shapes and colors
      scale_shape_manual(values = c(15, 17)) +
      scale_color_manual(values = c("#177F97", "#00AFAB", "#2EC6AF", 
                                    "#72008D", "#AB1488", "#D24E71", "#E8853A" ,"#DD6157","#ECC000"),                
                         labels = xlabs) +
      # Legends
      guides(color = guide_legend(title = "Treatment"))
    
    #ggsave("algalconsumption.alltreatments.bef.met.scaled.jan.3.0.6.png", device = "png", path = './figures/', width = 7, height = 4)
    
    
    ### 3. Expected vs observed consumption - adjusted for fish biomass 
    
### 3a. Calculate mean consumption per treatment and calculate Dmax -----
means <- data %>%
  dplyr::group_by(treatment) %>%
  dplyr::summarise(mean_delta = mean(delta.weight.total.s, na.rm = TRUE))

print(means)

#write.csv(as.data.frame(means), "bef.means.Met.scaled.csv")

# Identify polyculture value
poly_mean_6 <- means %>% filter(treatment == "6 species") %>% pull(mean_delta)
poly_mean_3 <- means %>% filter(treatment == "3 species") %>% pull(mean_delta)

# Identify the best (most negative = most consumption) monoculture
monoculture_means <- means %>%
  filter(treatment %in% c("manini", "velifer", "xanthopterus", "spilurus", "psittacus", "naso"))

best_mono <- min(monoculture_means$mean_delta, na.rm = TRUE)

Dmax_6 <- (poly_mean_6 - best_mono) / (best_mono)
print(Dmax_6)

Dmax_3 <- (poly_mean_3 - best_mono) / (best_mono)
print(Dmax_3)

### 3b. Predict total consumption of polycultures based on monocultures. ------

# Step 1: Subset to only the monoculture treatments
mono_data <- data %>%
  filter(diversity %in% c(0,1))

# Step 2: Fit the per-capita regression model using the number of individuals of each species
model <- lm(delta.weight.total.s ~ manini + velifer + xanthopterus + spilurus + psittacus + naso, 
            data = mono_data)

# test for normality in residuals 
shapiro.test(residuals(model))

car::vif(model)
cooks <- cooks.distance(model)
plot(cooks)

# model output 
tab_model(model)

# Step 3: Use this model to predict expected consumption for all treatments (including polycultures)
data$expected_polyculture <- predict(model, newdata = data)

# Step 4: Calculate the difference between observed and expected (residual) and Dt
data$residual_consumption <- data$delta.weight.total.s - data$expected_polyculture

lm_check <- lm(residual_consumption ~ total.weight, data = data)
summary(lm_check)

# Dt 
data$dt_relative <- (data$delta.weight.total.s - data$expected_polyculture) / data$expected_polyculture

## write.csv(data, "dt.data.bef.met.scaled.jan3.csv")

## Step 5: test for differences between observed and expected values 
poly_data <- data %>%
  dplyr::filter(treatment %in% c("3 species", "6 species"))

# Subset the data for each polyculture
#poly_3 <- subset(poly_data, treatment == "3 species")
#poly_6 <- subset(poly_data, treatment == "6 species")

# Paired test (observed vs expected)

# test normality assumptions 
shapiro.test(poly_data$residual_consumption)
qqnorm(poly_data$residual_consumption)
qqline(poly_data$residual_consumption)

## Abnormal based on q-q plot - use wilcox test 
wilcox_result <- wilcox.test(poly_data$residual_consumption)
print(wilcox_result)

    ###  Visualize observed vs expected consumption for polyculture treatments only (Figure 5a) ----
    
    # Summarize observed and expected means by treatment
    summary_df <- data %>%
      dplyr::group_by(treatment) %>%
      dplyr::summarise(
        mean_observed = mean(delta.weight.total.s, na.rm = TRUE),
        se_observed = sd(delta.weight.total.s, na.rm = TRUE) / sqrt(n()),
        mean_expected = mean(expected_polyculture, na.rm = TRUE),
        se_expected = sd(expected_polyculture, na.rm = TRUE) / sqrt(n())
      ) %>%
      pivot_longer(cols = c(mean_observed, mean_expected), 
                   names_to = "type", values_to = "value") %>%
      mutate(
        se = ifelse(type == "mean_observed", se_observed, se_expected),
        type = dplyr::case_when(
          type == "mean_observed" ~ "Observed",
          type == "mean_expected" ~ "Expected",
          TRUE ~ type    )  )
    
    # scatter plot 
    ggplot(poly_data, aes(x = expected_polyculture, y = delta.weight.total.s)) +
      geom_point(aes(color = treatment), size = 3) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +  # 1:1 line
      scale_x_continuous(limits = c(-26, 0), expand = c(0, 0)) +
      scale_y_continuous(limits = c(-26, 0), expand = c(0, 0)) +
     # scale_x_reverse() +
     # scale_y_reverse() +
     # facet_wrap(~treatment) +
      labs(x = "Mass-scaled expected consumption (g)", y = "Mass-scaled observed consumption (g)",
           color = "Treatment") +
      theme_minimal()
    
 # ggsave("exp.observed.scatter.met.scaled.merged.feb9.zero.png", device = "png", path = './figures/', width = 7, height = 4)
    

### 3c. Predict consumption by algal spices ------

# Step 0: Calculate mass scaled consumption for each algae type - normalizing for weight like the overall consumption done in excel
data$gs.delta.weight.total.s <- ifelse(is.finite(data$gs.delta.weight.total.g/ data$total.weight.s * 94),
                                  data$gs.delta.weight.total.g / data$total.weight.s* 94,
                                  data$gs.delta.weight.total.g)
data$as.delta.weight.total.s <- ifelse(is.finite(data$as.delta.weight.total.g / data$total.weight.s * 94),
                                  data$as.delta.weight.total.g / data$total.weight.s * 94,
                                  data$as.delta.weight.total.g)
data$ed.delta.weight.total.s <- ifelse(is.finite(data$ed.delta.weight.total.g / data$total.weight.s * 94),
                                  data$ed.delta.weight.total.g / data$total.weight.s * 94,
                                  data$ed.delta.weight.total.g)

# Step 1: Re-create mono data with the new column 
mono_data <- data %>%
  filter(diversity %in% c(0,1))

# Step 2: Fit the per-capita regression model for each algal species 

# Gracilaria model
mod_gs <- lm(gs.delta.weight.total.s ~ manini + velifer + xanthopterus + spilurus + psittacus + naso, data = mono_data)
tab_model(mod_gs)

# Acanthophora model
mod_acanth <- lm(as.delta.weight.total.s ~  manini + velifer + xanthopterus + spilurus + psittacus + naso, data = mono_data)
tab_model(mod_acanth)

# Eucheuma model
mod_eucheuma <- lm(ed.delta.weight.total.s ~ manini + velifer + xanthopterus + spilurus + psittacus + naso, data = mono_data)
tab_model(mod_eucheuma)

# Step 3: Use models to predict expected consumption for all treatments for each algal species 

## add results as columns in poly_data 
data$expected_gs.s <- predict(mod_gs, newdata = data)
data$expected_acanth.s <- predict(mod_acanth, newdata = data)
data$expected_eucheuma.s <- predict(mod_eucheuma, newdata = data)

# Step 4: Calculate the difference between observed and expected (residual)
poly_data <- data %>%
  dplyr::filter(treatment %in% c("3 species", "6 species"))

poly_data$residual_consumption_gs.s <- poly_data$gs.delta.weight.total.s - poly_data$expected_gs.s
poly_data$residual_consumption_acanth.s <- poly_data$as.delta.weight.total.s - poly_data$expected_acanth.s
poly_data$residual_consumption_eucheuma.s <- poly_data$ed.delta.weight.total.s - poly_data$expected_eucheuma.s

## Step 5: test for differences between observed and expected values 

# paired t-tests

## Gracilaria 
shapiro.test(poly_data$residual_consumption_gs.s)
qqnorm(poly_data$residual_consumption_gs.s)
qqline(poly_data$residual_consumption_gs.s)

t.test(poly_data$residual_consumption_gs.s)

## acanth
shapiro.test(poly_data$residual_consumption_acanth.s)
qqnorm(poly_data$residual_consumption_acanth.s)
qqline(poly_data$residual_consumption_acanth.s)

t.test(poly_data$residual_consumption_acanth.s)

## eucheuma 
shapiro.test(poly_data$residual_consumption_eucheuma.s)
qqnorm(poly_data$residual_consumption_eucheuma.s)
qqline(poly_data$residual_consumption_eucheuma.s)

t.test(poly_data$residual_consumption_eucheuma.s)

    ### Visualize observed vs expected consumption for polyculture treatments only (Figure 5b)----
    
    # Start with poly_data and reshape it manually into long format
    plot_data <- poly_data %>%
      select(treatment, replicate,
             gs_obs = gs.delta.weight.total.s,
             gs_exp = expected_gs.s,
             acanth_obs = as.delta.weight.total.s,
             acanth_exp = expected_acanth.s,
             eucheuma_obs = ed.delta.weight.total.s,
             eucheuma_exp = expected_eucheuma.s) %>%
      pivot_longer(
        cols = c(gs_obs, gs_exp, acanth_obs, acanth_exp, eucheuma_obs, eucheuma_exp),
        names_to = c("algae", "type"),
        names_sep = "_",
        values_to = "value"
      ) %>%
      pivot_wider(names_from = type, values_from = value) %>%
      rename(observed = obs, expected = exp)
    
    ### scatter plot with lines 
    plot_data$algae <- factor(plot_data$algae, levels = c("acanth", "eucheuma", "gs"))
    
    latin_labels <- c(
      acanth   = "italic('A. spicifera')",
      eucheuma = "italic('E. denticulatum')",
      gs       = "italic('G. salicornia')")
    
    ggplot(plot_data, aes(x = expected, y = observed)) +
      geom_point(aes(color = treatment), size = 3) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +  # 1:1 line
      scale_x_continuous(limits = c(-13, 1), expand = c(0, 0)) +
      scale_y_continuous(limits = c(-13, 1), expand = c(0, 0)) +
      #scale_x_reverse() +
      #scale_y_reverse() +
      facet_wrap(~ algae, labeller = as_labeller(latin_labels, default = label_parsed)) +
      labs(x = "Mass-scaled expected consumption (g)", y = "Mass-scaled observed consumption (g)",
           color = "Treatment") +
      theme_minimal()
    
    #ggsave("observed.expected.scaled.by.algal.mass.labs.feb13.png", device = "png", path = './figures/', width = 7, height = 4)
    
    
### Additional plots ------

## scaled consumption in BEF - species x species (figure S.6)
    algae_long <- data %>%
      select(treatment,
             gs.delta.weight.total.s,
             as.delta.weight.total.s,
             ed.delta.weight.total.s) %>%
      pivot_longer(
        cols = ends_with(".s"),
        names_to = "alga",
        values_to = "consumption_scaled"
      ) %>%
      mutate(
        alga = dplyr::recode(
          alga,
          "as.delta.weight.total.s" = "A. spicifera",
          "ed.delta.weight.total.s" = "E. denticulatum",
          "gs.delta.weight.total.s" = "G. salicornia"
        )
      )
    
    xlabs <- c("Control", "3 species", "6 species", 
               expression(italic("A.triostegus")), 
               expression(italic("S.psittacus")),
               expression(italic("C.spilurus")), 
               expression(italic("N.unicornis")), 
               expression(italic("Z.velifer")), 
               expression(italic("A.xanthopterus")))
    
    cols <- c("#177F97", "#00AFAB", "#2EC6AF",
              "#72008D", "#AB1488", "#D24E71",
              "#E8853A", "#DD6157", "#ECC000")
    
    ggplot(algae_long, aes(x = treatment, y = consumption_scaled, color = treatment)) +
      geom_jitter(width = 0.25, size = 3, alpha = 0.8) +
      geom_boxplot(aes(group = treatment),
                   color = "black",
                   fill = NA,
                   alpha = 0,
                   outlier.shape = NA) +
      facet_wrap(~ alga, nrow = 1) +
      scale_x_discrete(labels = xlabs) +
      scale_color_manual(values = cols, labels = xlabs) +
      labs(y = "Mass-scaled change in algal mass (g)", x = "Fish Treatment") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      guides(color = "none")  
    
    ggsave("algalconsumption.alltreatments.bef.met.scaled.by.aglae.jan.3.png", device = "png", path = './figures/', width = 7, height = 4)
    

# Observed vs expected for each replicate in BEF experiment 
      poly_3 <- subset(poly_data, treatment == "3 species")
      
      poly_3_long <- poly_3 %>%
        select(ID, observed = delta.weight.total.s, expected = expected_polyculture, replicate) %>%
        pivot_longer(cols = c(observed, expected),
                     names_to = "type",
                     values_to = "consumption")
      
      
      ggplot(poly_3_long, aes(x = factor(replicate), y = consumption, color = type, group = ID)) +
        geom_point(position = position_dodge(width = 0.3), size = 3) +
        geom_line(aes(group = ID), color = "gray", linetype = "dashed", alpha = 0.5) +
        labs(
          title = "Observed vs. Expected Consumption.s (3-Species Treatment)",
          x = "Replicate ID",
          y = "Consumption (g)",
          color = "Value Type"
        ) +
        scale_color_manual(values = c("observed" = "#1b9e77", "expected" = "#d95f02")) +
        theme_minimal()
      
      
      poly_6_long <- poly_6 %>%
        select(ID, observed = delta.weight.total.s, expected = expected_polyculture, replicate) %>%
        pivot_longer(cols = c(observed, expected),
                     names_to = "type",
                     values_to = "consumption")
      
      ggplot(poly_6_long, aes(x = factor(replicate), y = consumption, color = type, group = ID)) +
        geom_point(position = position_dodge(width = 0.3), size = 3) +
        geom_line(aes(group = ID), color = "gray", linetype = "dashed", alpha = 0.5) +
        labs(
          title = "Observed vs. Expected Consumption (6-Species Treatment)",
          x = "Replicate ID",
          y = "Consumption (g)",
          color = "Value Type"
        ) +
        scale_color_manual(values = c("observed" = "#1b9e77", "expected" = "#d95f02")) +
        theme_minimal()

 ## Observed vs expected by algal species 
        
        library(tidyverse)
        
        # Filter to 3-species treatment
        poly_3 <- subset(poly_data, treatment == "3 species")
        
        # Properly pivot to long format while preserving replicate ID
        # Prepare data
        poly_3_long_sp <- poly_3 %>%
          select(
            replicate,
            observed_gs = gs.delta.weight.total.s,
            expected_gs = expected_gs.s,
            observed_acanth = as.delta.weight.total.s,
            expected_acanth = expected_acanth.s,
            observed_eucheuma = ed.delta.weight.total.s,
            expected_eucheuma = expected_eucheuma.s
          ) %>%
          pivot_longer(
            cols = -replicate,
            names_to = c("type", "algae"),
            names_sep = "_",
            values_to = "consumption"
          ) %>%
          mutate(
            algae = dplyr::recode(algae,
                                  "gs" = "G.salicornia",
                                  "acanth" = "A.spicifera",
                                  "eucheuma" = "E.denticulatum"),
            type = dplyr::recode(type, "observed" = "Observed", "expected" = "Expected")
          )
        
        # Plot observed vs. expected for each algae, across replicates
        ggplot(poly_3_long_sp, aes(x = algae, y = consumption, color = type)) +
          geom_point(size = 3) +
          geom_line(aes(group = interaction(replicate, algae)), 
                    color = "gray40", 
                    linetype = "dashed", 
                    alpha = 0.6) +
          facet_wrap(~ replicate, ncol = 3) +
          labs(title = "Observed vs. Expected Algal Consumption by replicate (3-Species Treatment)",
               x = "Replicate ID",
               y = "Consumption (g)",
               color = "Value Type") +
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          scale_color_manual(values = c("Observed" = "#1b9e77", "Expected" = "#d95f02")) 
        
       # ggsave("observed.expected.3sp.by.algae.png", device = "png", path = './figures/', width = 7, height = 4)
        
        ### 6 sp 
        poly_6 <- subset(poly_data, treatment == "6 species")
        
        poly_6_long_sp <- poly_6 %>%
          select(
            replicate,
            observed_gs = gs.delta.weight.total.s,
            expected_gs = expected_gs.s,
            observed_acanth = as.delta.weight.total.s,
            expected_acanth = expected_acanth.s,
            observed_eucheuma = ed.delta.weight.total.s,
            expected_eucheuma = expected_eucheuma.s
          ) %>%
          pivot_longer(
            cols = -replicate,
            names_to = c("type", "algae"),
            names_sep = "_",
            values_to = "consumption"
          ) %>%
          mutate(
            algae = dplyr::recode(algae,
                                  "gs" = "G.salicornia",
                                  "acanth" = "A.spicifera",
                                  "eucheuma" = "E.denticulatum"),
            type = dplyr::recode(type, "observed" = "Observed", "expected" = "Expected")
          )
        
        # Plot observed vs. expected for each algae, across replicates
        ggplot(poly_6_long_sp, aes(x = algae, y = consumption, color = type)) +
          geom_point(size = 3) +
          geom_line(aes(group = interaction(replicate, algae)), 
                    color = "gray40", 
                    linetype = "dashed", 
                    alpha = 0.6) +
          facet_wrap(~ replicate, ncol = 3) +
          labs(title = "Observed vs. Expected Algal Consumption by Replicate (6-Species Treatment)",
               x = "Replicate ID",
               y = "Consumption (g)",
               color = "Value Type") +
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          scale_color_manual(values = c("Observed" = "#1b9e77", "Expected" = "#d95f02")) 
        
       # ggsave("observed.expected.6sp.by.algae.png", device = "png", path = './figures/', width = 7, height = 4)

        ## 3 species polyculture - consumption by replicate (Figure S7)
        
        p3data <- data[(data$treatment == "3 species"), ]
        
        plot(p3data$delta.weight.total.s ~ p3data$replicate,
             xlab = "Polyculture replicate",
             ylab = "Mass-scaled change in algal mass (g)",
             las = 1)       
        
        
### Correlations ####
algae.col <- sequential_hcl(3, palette = "TealGrn")
fish_data <- data[data$treatment != "control", ]

## correlation between fish weight and algae consumed 
plot(fish_data$delta.weight.total.g ~ fish_data$total.weight, col = algae.col,
     xlab = "Total fish weight (g)",
     ylab = "Change in algae weight (g)",
     las = 1,
     pch = 16,
     xlim = c(0,500))

lm_model <- lm(delta.weight.total.g ~ total.weight, data = fish_data)
abline(lm_model, col = "black", lwd = 2)  

shapiro.test(fish_data$delta.weight.total.g)
# p = 0.2591

shapiro.test(fish_data$total.weight)
#p= 0.0604

##both non- normal
cor.test(fish_data$delta.weight.total.g, fish_data$total.weight)

## not signifigant 

## MET scale assumptions 
plot(log(abs(fish_data$delta.weight.total.g)) ~ log(fish_data$total.weight), col = algae.col,
     xlab = "log(fish weight (g))",
     ylab = "log(abs(Change in algae weight (g))",
     las = 1,
     pch = 16)
     #xlim = c(2,5),
    # ylim = c(-3,2))

lm_model <- lm(log(abs(delta.weight.total.g)) ~ log(total.weight), data = fish_data)
summary(lm_model)
tab_model(lm_model)
slope <- coef(lm_model)[2]
slope 
linearHypothesis(lm_model, "log(total.weight) = 0.75")

abline(lm_model, col = "black", lwd = 2)  
