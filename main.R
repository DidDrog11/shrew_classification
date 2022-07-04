library(caret)
library(here)
library(nnet)
library(readxl)
library(rsample)
library(tidyverse)
library(truncnorm)

# Prepare data from literature

shrew_id <- read_csv(here("data", "rodent_ids.csv")) %>%
  filter(group == "shrew") %>%
      mutate(genus = str_split(name, "_", simplify = TRUE)[, 1])
    
mean_sd <- shrew_id %>%
  mutate(weight_sd = (weight_max - weight_mean)/4,
         head_body_sd = (head_body_max - head_body_min)/4,
         tail_sd = (tail_max -tail_mean)/4,
         hind_foot_sd = (hind_foot_max - hind_foot_mean)/4,
         ear_sd = (ear_max - ear_min)/4,
         length_skull_sd = (length_skull_max - length_skull_mean)/4)

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
  select(village, ID, species, weight, head_body, hb_tail_ratio, tail, hind_foot)

# We use a truncated normal distribution to sample 10,000 values from these distributions. The min and max are set to 75% and 125% of the min and max respectively
simulated_morphology <- list()
    
for(i in 1:length(unique(mean_sd$name))) {
  
  species = rep(mean_sd$name[i], each = 2000)
  subfamily = rep(mean_sd$sub_family[i], each = 2000)
  weight = rtruncnorm(n = 2000, mean = mean_sd$weight_mean[i], sd = mean_sd$weight_sd[i], a = mean_sd$weight_min[i]*0.75, b = mean_sd$weight_max[i]*1.25)
  head_body = rtruncnorm(n = 2000, mean = mean_sd$head_body_mean[i], sd = mean_sd$head_body_sd[i], a = mean_sd$head_body_min[i]*0.75, b = mean_sd$head_body_max[i]*1.25)
  tail = rtruncnorm(n = 2000, mean = mean_sd$tail_mean[i], sd = mean_sd$tail_sd[i], a = mean_sd$tail_min[i]*0.75, b = mean_sd$tail_max[i]*1.25)
  hind_foot = rtruncnorm(n = 2000, mean = mean_sd$hind_foot_mean[i], sd = mean_sd$hind_foot_sd[i], a = mean_sd$hind_foot_min[i]*0.75, b = mean_sd$hind_foot_max[i]*1.25)
  ear = rtruncnorm(n = 2000, mean = mean_sd$ear_mean[i], sd = mean_sd$ear_sd[i], a = mean_sd$ear_min[i]*0.75, b = mean_sd$ear_max[i]*1.25)
  hb_tail_ratio = tail/head_body
  
  simulated_morphology[[i]] <- tibble(species, subfamily, weight, head_body, tail, hb_tail_ratio, hind_foot, ear)
}

# BNITM data does not include measurements for ear length so this will not be included in this model

## Multinomial model for Crocidurinae
crocidurae_multinom_train <- bind_rows(simulated_morphology) %>%
  mutate(species = factor(species)) %>%
  select(-ear)

crocidurae_multinom_train$species <- droplevels(relevel(crocidurae_multinom_train$species, ref = "crocidura_olivieri"))

split_data <- initial_split(crocidurae_multinom_train, prop = 0.7, strata = "species")
train_data <- training(split_data)
test_data <- testing(split_data)

crocidurae_multinom_model_m1 <- multinom(species ~ weight + head_body + tail + hb_tail_ratio + hind_foot, data = train_data,
                                         maxit = 200)

train_data$species_predicted <- predict(crocidurae_multinom_model_m1, newdata = train_data, response = "class")
crocidura_train <- train_data
tab <- table(crocidura_train$species, crocidura_train$species_predicted)
# Accuracy
crocidura_accuracy <- round((sum(diag(tab))/sum(tab)) * 100, 2)

test_data$species_predicted <-  predict(crocidurae_multinom_model_m1, newdata = test_data, response = "class")
crocidura_test <- test_data

crocidura_confusion <- bind_rows(crocidura_train %>%
                                   group_by(species, species_predicted) %>%
                                   summarise(n = n()) %>%
                                   mutate(prop = n/sum(n)) %>%
                                   mutate(data = "train"),
                                 crocidura_test %>%
                                   group_by(species, species_predicted) %>%
                                   summarise(n = n()) %>%
                                   mutate(prop = n/sum(n)) %>%
                                   mutate(data = "test"))

levels(crocidura_confusion$species) <- sort(as.character(levels(crocidura_confusion$species)))
levels(crocidura_confusion$species_predicted) <- sort(as.character(levels(crocidura_confusion$species_predicted)))

crocidura_confusion_matrix <- ggplot(crocidura_confusion) +
  geom_tile(aes(x = species, y = species_predicted, fill = prop)) +
  facet_grid(~ data) +
  scale_fill_viridis_c() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(title = paste("Crocidura, accuracy =", crocidura_accuracy))

# Apply to BNITM data
bnitm_data$species_predicted <- predict(crocidurae_multinom_model_m1, newdata = bnitm_data, type = "class")
probabilities <- data.frame(round(predict(crocidurae_multinom_model_m1, newdata = bnitm_data, type = "probs"), 3))

bnitm_data <- bind_cols(bnitm_data, probabilities)

correct_classification <- bnitm_data %>%
  filter(species == species_predicted)
tab <- bnitm_data %>%
  group_by(species, species_predicted) %>%
  summarise(n = n()) %>%
  pivot_wider(names_from = species_predicted, values_from = n) %>%
  ungroup()

bnitm_crocidura_confusion <- bnitm_data %>%
  group_by(species, species_predicted) %>%
  summarise(n = n()) %>%
  mutate(prop = n/sum(n)) %>%
  mutate(data = "bnitm") %>%
  mutate(species = factor(species, levels = levels(crocidurae_multinom_train$species))) %>%
  drop_na(species_predicted)

bnitm_confusion_matrix <- ggplot(bnitm_crocidura_confusion) +
  geom_tile(aes(x = species, y = species_predicted, fill = prop)) +
  scale_fill_viridis_c() +
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(title = "Crocidura")

# Only 14 individuals have been correctly classified

weight_plot <- crocidurae_multinom_train %>%
  ggplot() +
  geom_freqpoly(aes(x = weight)) +
  geom_point(data = bnitm_data, aes(x = weight, y = 1), colour = "red") +
  facet_wrap(~ species, ncol = 1) +
  labs(title = "Shrew weight")

hb_plot <- crocidurae_multinom_train %>%
  ggplot() +
  geom_freqpoly(aes(x = head_body)) +
  geom_point(data = bnitm_data, aes(x = head_body, y = 1), colour = "red") +
  facet_wrap(~ species, ncol = 1) +
  labs(title = "Shrew Head-Body")

tail_plot <- crocidurae_multinom_train %>%
  ggplot() +
  geom_freqpoly(aes(x = tail)) +
  geom_point(data = bnitm_data, aes(x = tail, y = 1), colour = "red") +
  facet_wrap(~ species, ncol = 1) +
  labs(title = "Shrew tail")

hbt_plot <- crocidurae_multinom_train %>%
  ggplot() +
  geom_freqpoly(aes(x = hb_tail_ratio)) +
  geom_point(data = bnitm_data, aes(x = hb_tail_ratio, y = 1), colour = "red") +
  facet_wrap(~ species, ncol = 1) +
  labs(title = "Shrew HB-Tail ratio")

hf_plot <- crocidurae_multinom_train %>%
  ggplot() +
  geom_freqpoly(aes(x = hind_foot)) +
  geom_point(data = bnitm_data, aes(x = hind_foot, y = 1), colour = "red") +
  facet_wrap(~ species, ncol = 1) +
  labs(title = "Shrew hind foot")

report_input <- list(model = crocidurae_multinom_model_m1,
                     training_data = crocidurae_multinom_train,
                     training_confusion = crocidura_confusion_matrix,
                     bnitm_data = bnitm_data,
                     plots = list(weight_plot,
                                  hb_plot,
                                  tail_plot,
                                  hbt_plot,
                                  hf_plot))

write_rds(report_input, here("data", "data_for_report.rds"))
