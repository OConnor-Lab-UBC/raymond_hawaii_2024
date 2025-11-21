### Fish experiment 2: full choice assays

# 1. Load data and packages --------------------------------------------------

### load libraries
library(MuMIn)
library(tidyverse) 
library(reshape2)
library(gridExtra)
library(viridis)
library(googlesheets4)
library(ggplot2)
library(tidyr)
library(emmeans)
library(dplyr)
library(sjPlot)
library(emmeans)
library(lme4)
library(nlme)
library(dplyr)
library(purrr)
library(colorspace)
library(car)


## load data
data <- read_sheet('https://docs.google.com/spreadsheets/d/1kANQcC19OzdMkrYXXnrXRRlUR2ecax9yrbSCe9qfnaE/edit?gid=1713746961#gid=1713746961', sheet=1) 

# make factors 
data$fish.species <- as.factor(data$fish.species)
data$fish.species <- relevel(data$fish.species, ref = "Control")

# convert column to numeric 
columns_to_convert <- 9:28

for (i in columns_to_convert) {
  data[[i]] <- as.numeric(as.character(data[[i]]))
}

columns_numeric <- sapply(data[, columns_to_convert], is.numeric)

print(columns_numeric)

View(data)

# Define the columns representing algae consumption
algae_columns <- c("G.salicornia.delta.wt.g", 
                   "A.spicifera.delta.wt.g", 
                   "E.denticulatum.delta.wt.g", 
                   "H.discoidea.delta.wt.g")

# Gather algae consumption data into long format
data_long <- data %>%
  dplyr::select(fish.species, fish.id, fish.weight.g, date, all_of(algae_columns)) %>%
  tidyr::gather(key = "algal.species", value = "delta.weight.g",  -fish.species, - fish.id, - fish.weight.g, -date)

## set control to reference level 
data_long$fish.species <- relevel(data_long$fish.species, ref = "Control")

## make factors 
data_long$algal.species <- as.factor(data_long$algal.species)
data_long$algal.species <- factor(data_long$algal.species, levels=c("H.discoidea.delta.wt.g", "A.spicifera.delta.wt.g", 
                                                                    "E.denticulatum.delta.wt.g", "G.salicornia.delta.wt.g"))
data_long$date <- as.factor(data_long$date)

# Filter just the trial data (excluding controls)
trial_data <- data_long %>%
  filter(fish.species != "Control")

# Keep only control observations
control_data <- data_long %>%
  filter(fish.species == "Control")

# 2. Set up control matching pipeline -------

# For each trial row, store the vector of all matching control deltas - this lists all control values for that algal species on the given day
trial_with_ctrl <- trial_data %>%
  mutate(
    control_vals = map2(
      date, algal.species,
      ~ control_data %>%
        filter(date == .x, algal.species == .y) %>%
        pull(delta.weight.g)
    )
  )

# 3. Bootstrap adjustment function ----------

# Negative delta = consumption; positive = growth.
# Subtract a random control value and set positive adjusted values to 0 to account for any cases where algae grew more in trial than control 
adjust_once <- function(trial_with_ctrl, iter_id) {
  trial_with_ctrl %>%
    mutate(
      # sample one control value per row
      control_value = map_dbl(control_vals, ~ sample(.x, 1, replace = TRUE)),
      adjusted_weight = delta.weight.g - control_value,
      # if adjusted is POSITIVE (fish grew algae more than control), set to 0
      adjusted_weight = if_else(adjusted_weight > 0, 0, adjusted_weight)
    ) %>%
    group_by(fish.id) %>%
    mutate(
      total_consumption = sum(adjusted_weight, na.rm = TRUE),
      proportional_consumption = adjusted_weight / total_consumption,
      bootstrap_id = iter_id
    ) %>%
    ungroup()
}

# 4. Run bootstrap -------------------------------------------------------------------

n_boot <- 1000
set.seed(123)

boot_results <- map_dfr(1:n_boot, ~ adjust_once(trial_with_ctrl, .x))

# if any fish have total_consumption = 0, proportional_consumption will be NaN/Inf
boot_results <- boot_results %>%
  filter(!is.na(proportional_consumption))

boot_results_flat <- boot_results %>%
  dplyr::select(where(~ !is.list(.)))

# write.csv(boot_results_flat, "E2.bootresults.Nov.15.csv", row.names = FALSE)


# 5. Fit models and extract emmeans / pairwise comparisons for each bootstrap iteration  -----------

emm_list      <- list()
pairwise_list <- list()

boot_ids <- sort(unique(boot_results$bootstrap_id))

for (i in boot_ids) {
  message("Running model for bootstrap replicate ", i)
  
  boot_iter <- boot_results %>%
    filter(bootstrap_id == i) %>%
    filter(!is.na(proportional_consumption))
  
  # Safety: skip if somehow empty
  if (nrow(boot_iter) == 0) {
    message("  -> no data for replicate ", i, ", skipping.")
    next
  }
  
  # Fit the model; if it errors, print the message and skip
  mod <- tryCatch(
    {
      lmer(
        sqrt(proportional_consumption) ~ algal.species * fish.species + (1 | fish.id),
        data = boot_iter
      )
    },
    error = function(e) {
      message("  -> lmer error in replicate ", i, ": ", e$message)
      return(NULL)
    }
  )
  
  if (is.null(mod)) next
  
  # emmeans for this replicate
  emm   <- emmeans(mod, ~ algal.species * fish.species, type = "response")
  emm_df <- as.data.frame(emm)
  emm_df$bootstrap_id <- i
  emm_list[[length(emm_list) + 1]] <- emm_df
  
  # pairwise contrasts within fish.species
  pw <- pairs(emm, by = "fish.species") %>% as.data.frame()
  pw$bootstrap_id <- i
  pairwise_list[[length(pairwise_list) + 1]] <- pw
}

# Combine all non-NULL entries
all_emm_results <- bind_rows(emm_list)
all_pairs       <- bind_rows(pairwise_list)

dim(all_emm_results)
dim(all_pairs)

# 6. Summarise emmeans and pairwise across bootstraps -------------------------------------------------------------------

# Emmeans summaries
summary_emm <- all_emm_results %>%
  group_by(fish.species, algal.species) %>%
  summarise(
    mean_response = mean(response, na.rm = TRUE),
    lower_95      = quantile(response, 0.025, na.rm = TRUE),
    upper_95      = quantile(response, 0.975, na.rm = TRUE),
    .groups       = "drop"
  )

# write.csv(summary_emm, "E2.bootstrapped.emm.Nov.15.csv", row.names = FALSE)

# Pairwise summaries
summary_pairs <- all_pairs %>%
  group_by(fish.species, contrast) %>%
  summarise(
    mean_estimate    = mean(estimate, na.rm = TRUE),
    lower_95         = quantile(estimate, 0.025, na.rm = TRUE),
    upper_95         = quantile(estimate, 0.975, na.rm = TRUE),
    p_below_0        = mean(estimate < 0, na.rm = TRUE),
    p_above_0        = mean(estimate > 0, na.rm = TRUE),
    prop_significant = mean(p.value < 0.05, na.rm = TRUE),
    median_p_value   = median(p.value, na.rm = TRUE),
    .groups          = "drop"
  )

# write.csv(summary_pairs, "E2.bootstrapped.pairwise.Nov.15.csv", row.names = FALSE)

# 7. Model total proportional consumption by fish species --------------------
fish_data <- data[data$fish.present != "0", ]

shapiro.test((fish_data$delta.weight.total.s))
#normal

leveneTest(delta.weight.total.s  ~ fish.species, data = fish_data)
# vars equal 

best.mod <- lmer(delta.weight.total.s ~ fish.species + (1 | date), data = data, na.action = na.exclude)

tab_model(best.mod)
shapiro.test(residuals(best.mod))

## Estimate means 
EMM.e2 <- emmeans(best.mod, ~ fish.species)
EMM.e2

# write.csv(EMM.e2, "emmeans.E2.consumpton.met.scale.nov.15.csv")

## plot means with comparison arrows
plot(EMM.e2, comparisons = TRUE)

## complete pairwise comparisions 
pairs.e2 <- pairs(EMM.e2)

plot(pairs.e2)

# write.csv(pairs.e2, "pairwise.e2.rate.met.scale.nov.15.csv")

# 8. plot results ------------------------------------------------------

## create labels for fish and algae 

algae.labs <- c(expression(italic("H.discoidea")), expression(italic("A.spicifera")), expression(italic("E.denticulatum")), expression(italic("G.salicornia")))

fish.labs <- c("Convict Tang" = expression(italic("A.triostegus")), 
               "Palenose Parrotfish" = expression(italic("S.psittacus")),
               "Bullethead Parrotfish" = expression(italic("C.spilurus")), 
               "Bluespine Unicornfish" = expression(italic("N.unicornis")), 
               "Sailfin tang" = expression(italic("Z.velifer")), 
               "Yellowfin Surgeonfish" = expression(italic("A.xanthopterus")))

# 8.1 - Figure 2 (Connected box and whisker plot using estimated means from proportianal analysis) 
summary_emm$fish.species <- as.factor(summary_emm$fish.species)  
summary_emm$fish.species <- factor(summary_emm$fish.species, levels = c("Convict Tang","Palenose Parrotfish","Bullethead Parrotfish",
                                                                        "Bluespine Unicornfish","Sailfin tang", "Yellowfin Surgeonfish" ))

## plot bootstrapped data 
ggplot(summary_emm, aes(x = algal.species, y = mean_response, 
                        color = fish.species, group = fish.species, linetype = fish.species)) +
  geom_pointrange(aes(ymin = lower_95, ymax = upper_95),
                  position = position_dodge(width = 0.5), size = 0.4, fatten = 2.3) +
  geom_line(position = position_dodge(width = 0.5), size = 0.6) +
  theme_bw() +
  scale_x_discrete(labels = algae.labs) +
  xlab("Algae Species") +
  ylab("Proportional Consumption (mean ± 95% CI)") +
  ylim(-0.05, 1) +
  labs(title = "", color = "Fish Species") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_linetype_manual(values = c(
    "Bluespine Unicornfish" = "dashed",
    "Bullethead Parrotfish" = "solid",
    "Convict Tang" = "solid",
    "Palenose Parrotfish" = "solid",
    "Sailfin tang" = "solid",
    "Yellowfin Surgeonfish" = "solid")) +
  scale_color_manual(values = c("#72008D", "#AB1488", "#D24E71", "#E8853A", "#DD6157", "#ECC000"), 
                     labels = fish.labs) +
  guides(linetype = "none")

ggsave("dot.whisker.e2.bootstrapped.model.Nov.15.png", device = "png", path = './figures/', width = 7, height = 4)

## 8.2 Figure 3 - total adjusted algae consumed by species - where point size is un-adjusted fish weight (Fig) 
ggplot(data, aes(y = delta.weight.total.s, 
                 x = fish.species, 
                 col = fish.species)) + 
  geom_jitter(
    aes(size = ifelse(is.na(fish.weight.g) | fish.weight.g == 0, 2, fish.weight.g)), 
    width = 0.25, alpha = 0.8
  ) +  
  scale_x_discrete(labels = fish.labs) +  
  geom_boxplot(aes(group = fish.species),
               color = "black",   
               fill = NA,   
               alpha = 0,   
               outlier.size = 2) +  
  xlab("Fish Species") +
  ylab(expression("Scaled change in weight ("*g~d^{-1}*")")) +
  guides(
    x = guide_axis(angle = 45),
    color = "none",
    size = guide_legend(title = "Fish Weight (g)", override.aes = list(shape = 16))
  ) +
  scale_size_continuous(
    breaks = c(2, 30, 50, 70, 90),
    labels = c("Control","30g", "50g", "70g", "90g"),
    range = c(2, 6)
  ) +
  scale_color_manual(values = c("#177F97", "#72008D", "#AB1488", "#D24E71", "#E8853A", "#DD6157", "#ECC000")) +
  theme_minimal()

ggsave("algal.by.fish.met.scale.aug.1.png", device = "png", path = './figures/', width = 7, height = 4)

## 8.3 individual fish consumption from bootstrapped data (Figure S3) 
# Summarize across bootstraps per fish x algae
boot_summary <- boot_results %>%
  dplyr::group_by(fish.id, fish.species, algal.species) %>%
  dplyr::summarise(
    mean_pc = mean(proportional_consumption, na.rm = TRUE),
    lower_ci = quantile(proportional_consumption, 0.025, na.rm = TRUE),
    upper_ci = quantile(proportional_consumption, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

# Plot
ggplot(boot_summary, aes(x = algal.species, y = mean_pc, group = fish.id)) +
  geom_point(aes(color = fish.id)) +
  geom_line(aes(color = fish.id)) +
  # geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci, color = fish.id), width = 0.2) +
  theme_bw() +
  facet_wrap( ~ fish.species, ncol = 2) + 
  scale_x_discrete(labels = algae.labs) +
  xlab("Algae Species") +
  ylab("Mean Proportional Consumption (bootstrapped)") +
  labs(title = "Individual Fish Consumption (Bootstrapped Means)")

ggsave("mean.individual.consumption.png", device = "png", path = './figures/', width = 9, height = 6)


# 9. Correlation tests / check of MET assumptions ------

## 9.a correlation between log(fish weight) and log(abs(algae consumed))

fish_data <- data[data$fish.present != "0", ]
control_data <- data[data$fish.present != "1", ]

fish.col <- sequential_hcl(8, palette = "TealGrn")
algae.col <- sequential_hcl(4, palette = "TealGrn")

plot(log(abs(fish_data$delta.weight.total.g)) ~ log(fish_data$fish.weight.g), col = algae.col,
     xlab = "log(fish weight (g))",
     ylab = "log(abs(Change in algae weight (g))",
     las = 1,
     pch = 16,
     xlim = c(2,5),
     ylim = c(-3,2))

lm_model <- lm(log(abs(delta.weight.total.g)) ~ log(fish.weight.g), data = fish_data)
abline(lm_model, col = "black", lwd = 2)  

# 9.b correlation between fish weight and algae consumed 
plot(fish_data$delta.weight.total.g ~ fish_data$fish.weight.g, col = algae.col,
     xlab = "fish weight (g)",
     ylab = "Change in algae weight (g)",
     las = 1,
     pch = 16,
     ylim = c(-5,0.25),
     xlim = c(0,100))

lm_model <- lm(delta.weight.total.g ~ fish.weight.g, data = fish_data)
abline(lm_model, col = "black", lwd = 2)  

shapiro.test(fish_data$delta.weight.total.g)
# p = 0.4433

shapiro.test(fish_data$fish.weight.g)
#p= 0.1578

##both non- normal
cor.test(fish_data$delta.weight.total.g, fish_data$fish.weight.g, method = "kendall")

boxplot(data$fish.weight.g ~ data$fish.species,
        col = fish.col)

weight.aov <- aov(data$fish.weight.g ~ data$fish.species)
