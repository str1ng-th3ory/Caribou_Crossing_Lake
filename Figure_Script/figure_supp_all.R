library(dplyr)
library(randomForest)
library(caret)
library(ROCR)
library(pROC)
library(MuMIn)
library(car)
library(glmnet)
library(MASS)
library(ggplot2)

#Figure S1 - Panels A
#SEE 09_fall_models.R
importance_data = as.data.frame(randomForest::importance(rf_fall))
importance_data$Variable = c(
  "Proportion of Time Spent in Lake Area",
  "APR at the Nearest Pixel",
  "APR Along the Potential Crossing Path",
  "APR of the Entire Lake Area",
  "Lake Width",
  "The Ratio of Circumnavigate Speeds Over Reference Speeds (log)",
  "The Ratio of Direct Speeds Over Reference Speeds (log)",
  "The Ratio of Circumnavigate Speeds Over Direct Speeds (log)",
  "Sex",
  "Migration Season"
)

p_rf = ggplot(importance_data, aes(x = reorder(Variable, MeanDecreaseGini), y = MeanDecreaseGini)) +
  geom_point(size = 3, color = "blue") +
  labs(x = "Variables", y = "Importance (Gini Decrease)") +
  coord_flip() +
  theme_classic()

print(p_rf)

#Figure S1 - Panels B
#SEE 09_fall_models.R
threshold_data = data.frame(
  log_ratio_cir_ref = RF_fall$log_ratio_cir_ref,
  predicted_probs = predicted_probs,
  location = case_when(
    RF_fall$response == 1 ~ "Known Crossing Event",
    RF_fall$response == 0 ~ "Known Circumnavigating Event"
  )
)

log_ratio_cir_ref_grid = seq(min(RF_fall$log_ratio_cir_ref), max(RF_fall$log_ratio_cir_ref), length = 100)
new_data_ratio = data.frame(log_ratio_cir_ref = log_ratio_cir_ref_grid)

predicted_probs_with_se = predict(
  logistic_model_fall,
  newdata = new_data_ratio,
  type = "response",
  se.fit = TRUE
)

new_data_ratio$lower_bound = pmax(0, predicted_probs_with_se$fit - 1.96 * predicted_probs_with_se$se.fit)
new_data_ratio$upper_bound = pmin(1, predicted_probs_with_se$fit + 1.96 * predicted_probs_with_se$se.fit)
new_data_ratio$predicted_probs = predicted_probs_with_se$fit

p_logit = ggplot(new_data_ratio, aes(x = log_ratio_cir_ref, y = predicted_probs)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower_bound, ymax = upper_bound), alpha = 0.2) +
  geom_point(data = threshold_data, aes(x = log_ratio_cir_ref, y = predicted_probs, color = location), size = 2) +
  labs(
    x = "The Ratio of Circumnavigate Speeds Over Reference Speeds (log)",
    y = "Predicted Probability of Crossing",
    color = "Event"
  ) +
  theme_classic()

print(p_logit)

#Figure S1 - Panels C
#SEE 07_combined_models.R

importance_data = as.data.frame(randomForest::importance(rf))
importance_data$Variable = c(
  "Proportion of Time Spent in Lake Area",
  "APR at the Nearest Pixel",
  "APR Along the Potential Crossing Path",
  "APR of the Entire Lake Area",
  "Lake Width",
  "The Ratio of Circumnavigate Speeds Over Reference Speeds (log)",
  "The Ratio of Direct Speeds Over Reference Speeds (log)",
  "The Ratio of Circumnavigate Speeds Over Direct Speeds (log)",
  "Sex",
  "Migration Season"
)

p_rf = ggplot(importance_data, aes(x = reorder(Variable, MeanDecreaseGini), y = MeanDecreaseGini)) +
  geom_point(size = 3, color = "blue") +
  labs(
    x = "Variables",
    y = "Importance (Gini Decrease)"
  ) +
  coord_flip() +
  theme_classic()

print(p_rf)




#Figure S1 - Panels D
#SEE 07_combined_models.R
threshold_data = data.frame(
  albedo_linear = RF_data$albedo_linear,
  lake_spent_proportion = RF_data$lake_spent_proportion,
  log_ratio_cir_ref = RF_data$log_ratio_cir_ref,
  predicted_probs = predicted_probs,
  location = case_when(
    RF_data$response == 1 ~ "Known Crossing Event",
    RF_data$response == 0 ~ "Known Circumnavigating Event"
  )
)

albedo_grid = seq(min(RF_data$albedo_linear), max(RF_data$albedo_linear), length = 100)
lake_spent_grid = seq(min(RF_data$lake_spent_proportion), max(RF_data$lake_spent_proportion), length = 100)
log_ratio_cir_ref_grid = seq(min(RF_data$log_ratio_cir_ref), max(RF_data$log_ratio_cir_ref), length = 100)

new_data_albedo = data.frame(
  albedo_linear = albedo_grid,
  lake_spent_proportion = median(RF_data$lake_spent_proportion),
  log_ratio_cir_ref = median(RF_data$log_ratio_cir_ref)
)

pred_albedo = predict(model_all_logistic, newdata = new_data_albedo, type = "response", se.fit = TRUE)

new_data_albedo$lower_bound = pmax(0, pred_albedo$fit - 1.96 * pred_albedo$se.fit)
new_data_albedo$upper_bound = pmin(1, pred_albedo$fit + 1.96 * pred_albedo$se.fit)
new_data_albedo$predicted_probs = pred_albedo$fit

p = ggplot(new_data_albedo, aes(x = albedo_linear, y = predicted_probs)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower_bound, ymax = upper_bound), alpha = 0.2) +
  geom_point(data = threshold_data, aes(x = albedo_linear, y = predicted_probs, color = location), size = 1.5) +
  labs(
    x = "Albedo Along the Potential Crossing Path",
    y = "Predicted Probability of Crossing",
    color = "Event"
  ) +
  theme_classic()

print(p)


#Figure S2 
#SEE 09_fall_models.R

unknown_data_fall = Unknown_event_fall %>%
  mutate(
    sex = as.factor(sex),
    season = as.factor(season)
  ) %>%
  dplyr::select(
    lake_spent_proportion,
    albedo_nearest_pixel,
    albedo_linear,
    albedo_whole_lake,
    lake_width,
    log_ratio_cir_ref,
    log_ratio_dir_ref,
    log_ratio_cir_dir,
    sex,
    season
  ) %>%
  filter(rowSums(is.na(.)) == 0)

unknown_predicted_probs_fall_rf = predict(rf_fall, unknown_data_fall, type = "prob")
unknown_data_fall$rf_predicted_probs = unknown_predicted_probs_fall_rf[, 2]

unknown_predicted_probs_fall_logit = predict(
  logistic_model_fall,
  newdata = unknown_data_fall,
  type = "response",
  se.fit = TRUE
)

unknown_data_fall$predicted_probs = unknown_predicted_probs_fall_logit$fit
unknown_data_fall$lower_bound = pmax(0, unknown_predicted_probs_fall_logit$fit - 1.96 * unknown_predicted_probs_fall_logit$se.fit)
unknown_data_fall$upper_bound = pmin(1, unknown_predicted_probs_fall_logit$fit + 1.96 * unknown_predicted_probs_fall_logit$se.fit)

p_unknown = ggplot(new_data_ratio, aes(x = log_ratio_cir_ref, y = predicted_probs)) +
  geom_line(color = "blue", linewidth = 1) +
  geom_point(data = threshold_data, aes(x = log_ratio_cir_ref, y = predicted_probs, color = location), size = 2) +
  geom_point(data = unknown_data_fall, aes(x = log_ratio_cir_ref, y = predicted_probs, color = "Unknown Event (Binomial Logistic Model)"), size = 2) +
  geom_point(data = unknown_data_fall, aes(x = log_ratio_cir_ref, y = rf_predicted_probs, color = "Unknown Event (Random Forest Model)"), size = 2) +
  labs(
    x = "The Ratio of Circumnavigate Speeds over Reference Speeds (log)",
    y = "Predicted Probability of Crossing",
    color = "Event"
  ) +
  scale_color_manual(
    values = c(
      "Known Crossing Event" = "#00BFC4",
      "Known Circumnavigating Event" = "#F8766D",
      "Unknown Event (Binomial Logistic Model)" = "orange",
      "Unknown Event (Random Forest Model)" = "#C77CFF"
    )
  ) +
  theme_classic()

print(p_unknown)

