# ----------------------------------------
# Figure6_spring_model_results.R
# Purpose:
#   Create Figure 6 for the spring-only dataset.
#
# Panels:
#   A. Random Forest variable importance
#   B. Binomial Logistic probability curve with
#      known-event observations and unknown-event predictions
#
# Prerequisite:
#   Run 08_spring_models.R first so that the following objects exist:
#   - rf_spring
#   - logistic_model_spring
#   - RF_spring
#   - Unknown_event_spring
#   - predicted_probs
#   - threshold_data
#   - new_data_ratio
# ----------------------------------------

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

# =========================
# Panel A. RF variable importance
# =========================
importance_data = as.data.frame(randomForest::importance(rf_spring))
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

p_fig6A = ggplot(
  importance_data,
  aes(x = reorder(Variable, MeanDecreaseGini), y = MeanDecreaseGini)
) +
  geom_point(size = 3, color = "blue") +
  coord_flip() +
  theme_classic() +
  labs(
    x = "Variables",
    y = "Importance (Gini Decrease)"
  )

print(p_fig6A)

# =========================
# Prepare unknown spring data for panel B
# =========================
unknown_data_spring = Unknown_event_spring %>%
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


unknown_data_spring$season = droplevels(unknown_data_spring$season)

# RF predictions
unknown_predicted_probs_spring_rf = predict(
  rf_spring,
  unknown_data_spring,
  type = "prob"
)
unknown_data_spring$rf_predicted_probs = unknown_predicted_probs_spring_rf[, 2]

# BL predictions
unknown_predicted_probs_spring_logit = predict(
  logistic_model_spring,
  newdata = unknown_data_spring,
  type = "response",
  se.fit = TRUE
)

unknown_data_spring$predicted_probs = unknown_predicted_probs_spring_logit$fit
unknown_data_spring$lower_bound = pmax(0, unknown_predicted_probs_spring_logit$fit - 1.96 * unknown_predicted_probs_spring_logit$se.fit)
unknown_data_spring$upper_bound = pmin(1, unknown_predicted_probs_spring_logit$fit + 1.96 * unknown_predicted_probs_spring_logit$se.fit)

# =========================
# Panel B. BL curve + known and unknown events
# =========================

p_fig6B = ggplot(new_data_ratio, aes(x = albedo_linear, y = predicted_probs)) +
  geom_line(color = "blue", linewidth = 0.6) +
  geom_ribbon(aes(ymin = lower_bound, ymax = upper_bound), alpha = 0.2) +
  geom_point(
    data = threshold_data,
    aes(x = albedo_linear, y = predicted_probs, color = location),
    size = 1.5
  ) +
  geom_point(
    data = unknown_data_spring,
    aes(
      x = albedo_linear,
      y = predicted_probs,
      color = "Unknown Event (Binomial Logistic Model)"
    ),
    size = 1.5
  ) +
  geom_point(
    data = unknown_data_spring,
    aes(
      x = albedo_linear,
      y = rf_predicted_probs,
      color = "Unknown Event (Random Forest Model)"
    ),
    size = 1.5
  ) +
  scale_color_manual(
    values = c(
      "Known Crossing Event" = "#00BFC4",
      "Known Circumnavigating Event" = "#F8766D",
      "Unknown Event (Binomial Logistic Model)" = "orange",
      "Unknown Event (Random Forest Model)" = "#C77CFF"
    )
  ) +
  theme_classic() +
  labs(
    x = "APR Along the Potential Crossing Path",
    y = "Predicted Probability of Crossing",
    color = "Event"
  )

print(p_fig6B)

# =========================
# Combine panels
# =========================
fig6_combined = ggarrange(
  p_fig6A,
  p_fig6B,
  labels = c("A", "B"),
  ncol = 2,
  nrow = 1,
  common.legend = FALSE
)

print(fig6_combined)

# Optional save
# ggsave("Figure6_spring_model_results.png", fig6_combined, width = 12, height = 6, dpi = 300)


