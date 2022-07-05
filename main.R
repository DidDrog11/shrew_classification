library(caret)
library(here)
library(nnet)
library(readxl)
library(rsample)
library(tidyverse)
library(truncnorm)

seed(123)
# Prepare data from literature

shrew_id <- read_csv(here("data", "rodent_ids.csv")) %>%
  filter(group == "shrew") %>%
      mutate(genus = str_split(name, "_", simplify = TRUE)[, 1]) %>%
  rename(species = name)

# Calculate the SD of the weight of each shrew species from it's min and max weight

weight <- tibble(species = shrew_id$species,
                 weight_mean = shrew_id$weight_mean,
                 weight_sd = (shrew_id$weight_max - shrew_id$weight_min)/4,
                 weight_min = shrew_id$weight_min,
                 weight_max = shrew_id$weight_max)

# Other measurements are likely to be associated with the weight of the shrew
# Produce three values for each species and measurement associated with the mean, min and max weights
# Assumption is that lightest specimens were also the smallest

measurements <- tibble(species = rep_len(shrew_id$species, 3*length(shrew_id$species)),
                       weight = c(shrew_id$weight_mean, shrew_id$weight_min, shrew_id$weight_max),
                       hb =c(shrew_id$head_body_mean, shrew_id$head_body_min, shrew_id$head_body_max),
                       tail = c(shrew_id$tail_mean, shrew_id$tail_min, shrew_id$tail_max),
                       hind_foot = c(shrew_id$hind_foot_mean, shrew_id$hind_foot_min, shrew_id$hind_foot_max))

# Assuming a linear relationship between weight and head body length we can model the relationship and subsequently use this
# to produce HB for simulated weights of the shrews

hb_model <- glm(hb ~ weight * species, data =  measurements, family = "gaussian")
tail_model <- glm(tail ~ weight * species, data =  measurements, family = "gaussian")
hf_model <- glm(hind_foot ~ weight * species, data =  measurements, family = "gaussian")

# Check if these produce reasonable predictions for HB, tail and hind foot

measurements$pred_hb <- predict(hb_model, newdata = measurements)
measurements$pred_tail <- predict(tail_model, newdata = measurements)
measurements$pred_hf <- predict(hf_model, newdata = measurements)

ggplot(measurements) +
  geom_point(aes(x = hb, y = pred_hb, colour = species)) +
  geom_line(aes(x = hb, y = pred_hb, colour = species)) +
  labs(title = "Predicted head-body length from weight",
       x = "Head-body length",
       y = "Predicted head-body length") +
  theme_bw()

ggplot(measurements) +
  geom_point(aes(x = tail, y = pred_tail, colour = species)) +
  geom_line(aes(x = tail, y = pred_tail, colour = species)) +
  labs(title = "Predicted tail length from weight",
       x = "Tail length",
       y = "Predicted tail length") +
  theme_bw()

ggplot(measurements) +
  geom_point(aes(x = hind_foot, y = pred_hf, colour = species)) +
  geom_line(aes(x = hind_foot, y = pred_hf, colour = species)) +
  labs(title = "Predicted hind-foot length from weight",
       x = "Hind-foot length",
       y = "Predicted hind-foot length") +
  theme_bw()
# I think these are acceptable


# Simulating data from weights --------------------------------------------
# We use a truncated normal distribution to sample 5,000 values from the distribution of a shrews weight.
# We set the SD as 2 times the SD to produce a range of 95% of the expected values
# The min and max are set to 50% and 150% of the min and max respectively to account for the likely unsampled variability and limit biologically implausible values.
# The models produced above use the 

sim_measures <- list()

for(i in 1:length(unique(weight$species))) {
  
  sim_measures[[i]] <- tibble(species = rep(weight$species[i], each = 5000),
                             weight = rtruncnorm(n = 5000, mean = weight$weight_mean[i], sd = weight$weight_sd[i] * 2, a = weight$weight_min[i]*0.5, b = weight$weight_max[i]*1.5))
  
  sim_measures[[i]]$head_body <- predict(hb_model, newdata = sim_measures[[i]])
  sim_measures[[i]]$tail <- predict(tail_model, newdata = sim_measures[[i]])
  sim_measures[[i]]$hind_foot <- predict(hf_model, newdata = sim_measures[[i]])
  sim_measures[[i]]$hb_tail_ratio <- sim_measures[[i]]$tail/sim_measures[[i]]$head_body
  
}

simulated_morphology <- bind_rows(sim_measures) %>%
  mutate(species = factor(species))


# Training model on this simulated dataset --------------------------------

simulated_morphology$species <- droplevels(relevel(simulated_morphology$species, ref = "crocidura_olivieri"))

split_data <- initial_split(simulated_morphology, prop = 0.7, strata = "species")
train_data <- training(split_data)
test_data <- testing(split_data)

crocidurae_multinom_model_m1 <- multinom(species ~ weight + head_body + tail + hind_foot, data = train_data,
                                         maxit = 200)

crocidurae_multinom_model_m2 <- multinom(species ~ weight + head_body + tail + hb_tail_ratio + hind_foot, data = train_data,
                                         maxit = 200)

# As sensitivity analysis we will limit to the shrew species that have been detected in the region
sens_1 <- simulated_morphology %>%
  filter(species %in% c("crocidura_buettikoferi", "crocidura_theresae", "crocidura_grandiceps", "crocidura_olivieri"))

split_data_s1 <- initial_split(sens_1, prop = 0.7, strata = "species")
train_data_s1 <- training(split_data_s1)
test_data_s1 <- testing(split_data_s1)

crocidurae_multinom_model_s1 <- multinom(species ~ weight + head_body + tail + hind_foot, data = train_data_s1,
                                         maxit = 200)

# Check the accuracy of these models on withheld test data

test_data$m1_pred <-  predict(crocidurae_multinom_model_m1, newdata = test_data, response = "class")
test_data$m2_pred <-  predict(crocidurae_multinom_model_m2, newdata = test_data, response = "class")
test_data_s1$s1_pred <-  predict(crocidurae_multinom_model_s1, newdata = test_data_s1, response = "class")

tab_m1 <- table(test_data$species, test_data$m1_pred)
m1_accuracy <- round((sum(diag(tab_m1))/sum(tab_m1)) * 100, 2)
tab_m2 <- table(test_data$species, test_data$m2_pred)
m2_accuracy <- round((sum(diag(tab_m2))/sum(tab_m2)) * 100, 2)

tab_s1 <- test_data_s1 %>%
  mutate(species = as.character(species),
         s1_pred = as.character(s1_pred)) %>%
  janitor::tabyl(species, s1_pred)

# There is no important accuracy increase from including hb_tail_ratio, so this will be removed as it will be colinear with both hb and tail

crocidura_confusion <- bind_rows(test_data %>%
                                   rename(species_predicted = m1_pred) %>%
                                   group_by(species, species_predicted) %>%
                                   summarise(n = n()) %>%
                                   mutate(prop = n/sum(n)) %>%
                                   mutate(data = "test",
                                          species = factor(str_to_sentence(str_replace_all(species, "_", " "))),
                                          species_predicted = factor(str_to_sentence(str_replace_all(species_predicted, "_", " ")))))

crocidura_confusion$species <- factor(crocidura_confusion$species, levels = unique(sort(as.character(str_to_sentence(str_replace_all(simulated_morphology$species, "_", " "))))))
crocidura_confusion$species_predicted <- factor(crocidura_confusion$species_predicted, levels = unique(sort(as.character(str_to_sentence(str_replace_all(simulated_morphology$species, "_", " "))))))

crocidura_confusion_matrix <- ggplot(crocidura_confusion) +
  geom_tile(aes(x = species, y = species_predicted, fill = prop)) +
  scale_fill_viridis_c() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(title = paste("Crocidura, accuracy =", m1_accuracy),
       x = "Species",
       y = "Predicted species",
       fill = "Proportion")


# Use to make out of sample predictions -----------------------------------
# Prepare BNITM data
bnitm_data <- read_xlsx(here("data", "bnitm_data.xlsx")) %>%
  rename(village = Village,
         ID = 2,
         genus = Genus,
         species = sp,
         weight = Wt,
         head_body = HB,
         tail = "T",
         hind_foot = HF) %>%
  mutate(tail = as.numeric(tail),
         hb_tail_ratio = tail/head_body,
         species = paste0(str_to_lower(genus), "_", species)) %>%
  select(village, ID, species, weight, head_body, tail, hind_foot) %>%
  drop_na(weight, head_body, tail, hind_foot) # A single shrew is missing tail measurements

# Apply to BNITM data
bnitm_data$species_predicted <- predict(crocidurae_multinom_model_m1, newdata = bnitm_data, type = "class")
probabilities <- data.frame(round(predict(crocidurae_multinom_model_m1, newdata = bnitm_data, type = "probs"), 3))
sensitivity <- bind_cols(bnitm_data %>%
                           rename(m1_pred = species_predicted),
                         s1_pred = predict(crocidurae_multinom_model_s1, newdata = bnitm_data, type = "class"))

bnitm_data <- bind_cols(bnitm_data, probabilities)

correct_classification <- bnitm_data %>%
  filter(species == species_predicted)

bnitm_crocidura_confusion <- bnitm_data %>%
  group_by(species, species_predicted) %>%
  summarise(n = n()) %>%
  mutate(prop = n/sum(n)) %>%
  mutate(data = "bnitm") %>%
  mutate(species = fct_relevel(factor(str_to_sentence(str_replace_all(species, "_", " ")), levels = levels(crocidura_confusion$species))),
         species_predicted = fct_relevel(factor(str_to_sentence(str_replace_all(species_predicted, "_", " ")), levels = levels(crocidura_confusion$species))))

bnitm_confusion_matrix <- ggplot(bnitm_crocidura_confusion) +
  geom_tile(aes(x = species, y = species_predicted, fill = prop)) +
  scale_fill_viridis_c() +
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(title = "Crocidura")

# Only 5 individuals have been correctly classified

weight_plot <- simulated_morphology %>%
  ggplot() +
  geom_freqpoly(aes(x = weight)) +
  geom_point(data = bnitm_data, aes(x = weight, y = 0.1), colour = "red") +
  facet_wrap(~ species, ncol = 2) +
  labs(title = "Shrew weight")

hb_plot <- simulated_morphology %>%
  ggplot() +
  geom_freqpoly(aes(x = head_body)) +
  geom_point(data = bnitm_data, aes(x = head_body, y = 0.1), colour = "red") +
  facet_wrap(~ species, ncol = 2) +
  labs(title = "Shrew head-body")

tail_plot <- simulated_morphology %>%
  ggplot() +
  geom_freqpoly(aes(x = tail)) +
  geom_point(data = bnitm_data, aes(x = tail, y = 1), colour = "red") +
  facet_wrap(~ species, ncol = 2) +
  labs(title = "Shrew tail")

hf_plot <- simulated_morphology %>%
  ggplot() +
  geom_freqpoly(aes(x = hind_foot)) +
  geom_point(data = bnitm_data, aes(x = hind_foot, y = 1), colour = "red") +
  facet_wrap(~ species, ncol = 2) +
  labs(title = "Shrew hind foot")

report_input <- list(model = crocidurae_multinom_model_m1,
                     training_data = simulated_morphology,
                     training_confusion = crocidura_confusion_matrix,
                     bnitm_data = bnitm_data,
                     bnitm_confusion_matrix = bnitm_confusion_matrix,
                     plots = list(weight_plot,
                                  hb_plot,
                                  tail_plot,
                                  hf_plot))

write_rds(report_input, here("data", "data_for_report.rds"))
