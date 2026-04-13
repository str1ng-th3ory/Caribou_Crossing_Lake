# ----------------------------------------
# 07_combined_models.R
# Purpose:
#   Fit combined spring + fall models for known events,
#   and predict unknown events using:
#   (1) Random Forest
#   (2) Binomial Logistic Regression
#
# Input:
#   - known_event_all.csv
#   - unknown_event_all.csv
#
# Output:
#   - combined_rf_test_predictions.csv
#   - combined_unknown_predictions.csv
#   - combined_model_input.csv
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
library(detectseparation)

# =========================
# 1. Load data
# =========================
Known_event_all = read.csv("known_event_all.csv", stringsAsFactors = FALSE)
Unknown_event_all = read.csv("unknown_event_all.csv", stringsAsFactors = FALSE)

# =========================
# 2. Prepare combined known-event modeling table
# =========================
RF_data = Known_event_all %>%
  mutate(
    response = factor(response, levels = c(0, 1)),
    sex = as.factor(sex),
    season = as.factor(season)
  ) %>%
  dplyr::select(
    response,
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
  )

# optional complete-case filter for modeling
RF_data = RF_data %>%
  filter(rowSums(is.na(.)) == 0)

write.csv(RF_data, "combined_model_input.csv", row.names = FALSE)

# =========================
# 3. Random Forest
# =========================

train_indices = sample(1:nrow(RF_data), 0.8 * nrow(RF_data), replace = FALSE)
train_data = RF_data[train_indices, ]
test_data  = RF_data[-train_indices, ]

class_counts = table(train_data$response)
total_samples = sum(class_counts)
num_classes = length(class_counts)

class_weights = total_samples / (num_classes * class_counts)
weight_for_0=class_weights[[1]]
weight_for_1=class_weights[[2]]

class_wt <- c(`0` = weight_for_0, `1` = weight_for_1)

rf = randomForest(
  response ~ .,
  data = train_data,
  mtry = 2,
  importance = TRUE,
  ntree = 1000,
  classwt = class_wt
)

print(rf)
varImpPlot(rf)

test_pred_class = predict(rf, test_data)
test_pred_prob  = predict(rf, test_data, type = "prob")

test_data$response = factor(test_data$response, levels = levels(test_pred_class))

conf_matrix = confusionMatrix(
  test_pred_class,
  test_data$response,
  mode = "everything",
  positive = "1"
)
print(conf_matrix)

f1_score = conf_matrix$byClass["F1"]
cat("RF F1 Score:", f1_score, "\n")

perf2 = prediction(test_pred_prob[, 2], test_data$response)
auc2 = performance(perf2, "auc")
cat("RF AUC:", auc2@y.values[[1]], "\n")

roc_curve_rf = roc(test_data$response, test_pred_prob[, 2])
plot(roc_curve_rf)
print(auc(roc_curve_rf))

rf_test_output = test_data
rf_test_output$pred_class = test_pred_class
rf_test_output$pred_prob = test_pred_prob[, 2]

write.csv(rf_test_output, "combined_rf_test_predictions.csv", row.names = FALSE)

# =========================
# plot
# =========================

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


# =========================
# 4. Logistic regression
# =========================
# check multi-collinearity before do this!

all_predictors = c(
  "lake_spent_proportion",
  "albedo_nearest_pixel",
  "albedo_linear",
  "albedo_whole_lake",
  "lake_width",
  "log_ratio_cir_ref",
  "log_ratio_dir_ref",
  "sex",
  "season"
)

formula_full = as.formula(
  paste("response ~", paste(all_predictors, collapse = " + "))
)

full_model = glm(formula_full, family = "binomial", data = RF_data, na.action = na.fail)

print(summary(full_model))

separation_check = detectseparation::detect_separation(full_model)
print(separation_check)

vif_values = vif(full_model)
print(vif_values)

# LASSO
x = model.matrix(formula_full, RF_data)[, -1]
y = RF_data$response

cv_fit = cv.glmnet(x, y, family = "binomial", alpha = 1)
best_lambda = cv_fit$lambda.min
lasso_model = glmnet(x, y, family = "binomial", alpha = 1, lambda = best_lambda)
print(coef(lasso_model))

# =========================
# Combined logistic model selection and diagnostics
# =========================

# reduced candidate model
adjust_model <- glm(
  response ~ lake_spent_proportion + albedo_linear + log_ratio_cir_ref,
  family = "binomial",
  data = RF_data
)

summary(adjust_model)

# -------------------------
# Stepwise regression
# -------------------------
step_model <- MASS::stepAIC(adjust_model, direction = "both", trace = FALSE)
summary(step_model)

# -------------------------
# Dredge on reduced model
# -------------------------
dredge_adjust_model <- glm(
  response ~ lake_spent_proportion + albedo_linear + log_ratio_cir_ref,
  family = "binomial",
  data = RF_data,
  na.action = na.fail
)

dredge_results_whole <- MuMIn::dredge(dredge_adjust_model)
print(dredge_results_whole)
summary(get.models(dredge_results_whole, 1)[[1]])

# -------------------------
# Check whether quadratic terms may be needed
# using binned observed crossing proportions
# -------------------------
k <- 12

tmp <- RF_data %>%
  mutate(bin = ntile(albedo_linear, k)) %>%
  group_by(bin) %>%
  summarise(
    albedo_mid = median(albedo_linear, na.rm = TRUE),
    p_cross = mean(response == 1, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(tmp, aes(albedo_mid, p_cross)) +
  geom_point() +
  geom_line() +
  labs(x = "Albedo", y = "Observed crossing proportion")

tmp <- RF_data %>%
  mutate(bin = ntile(lake_spent_proportion, k)) %>%
  group_by(bin) %>%
  summarise(
    lake_spent_proportion_mid = median(lake_spent_proportion, na.rm = TRUE),
    p_cross = mean(response == 1, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(tmp, aes(lake_spent_proportion_mid, p_cross)) +
  geom_point() +
  geom_line() +
  labs(x = "Lake spent proportion", y = "Observed crossing proportion")

tmp <- RF_data %>%
  mutate(bin = ntile(log_ratio_cir_ref, k)) %>%
  group_by(bin) %>%
  summarise(
    log_ratio_cir_ref_mid = median(log_ratio_cir_ref, na.rm = TRUE),
    p_cross = mean(response == 1, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(tmp, aes(log_ratio_cir_ref_mid, p_cross)) +
  geom_point() +
  geom_line() +
  labs(x = "log_ratio_cir_ref", y = "Observed crossing proportion")

# Visual inspection suggested monotonic relationships.
# No clear U-shaped or inverted-U pattern was observed,
# so quadratic terms were not included.

# -------------------------
# Interaction with season
# -------------------------
RF_data$season <- as.factor(RF_data$season)

cont_vars  <- c("lake_spent_proportion", "albedo_linear", "log_ratio_cir_ref")
factor_var <- "season"

rhs_main <- paste(c(cont_vars, factor_var), collapse = " + ")
rhs_int  <- paste(
  c("albedo_linear:season",
    "lake_spent_proportion:season"),
  collapse = " + "
)

formula_global <- as.formula(
  paste0("response ~ ", rhs_main, " + ", rhs_int)
)

global_fit <- glm(
  formula_global,
  family = binomial,
  data = RF_data,
  na.action = na.fail
)

subset_txt <- paste(
  "(!`albedo_linear:season` || (albedo_linear & season))",
  "(!`lake_spent_proportion:season` || (lake_spent_proportion & season))",
  sep = " & "
)
subset_call <- parse(text = subset_txt)[[1]]

dredge_tbl <- dredge(global_fit, subset = subset_call)
print(dredge_tbl)

best_mod <- get.models(dredge_tbl, 1)[[1]]
summary(best_mod)

# after excluding models with convergence problems and quasi-complete separation,
# the top supported model remained:
m_lin <- glm(
  response ~ lake_spent_proportion + albedo_linear + log_ratio_cir_ref,
  family = "binomial",
  data = RF_data
)

summary(m_lin)

# -------------------------
# Final combined logistic model performance
# -------------------------
model_all_logistic <- glm(
  response ~ albedo_linear + lake_spent_proportion + log_ratio_cir_ref,
  family = "binomial",
  data = RF_data
)

predicted_probs <- predict(model_all_logistic, type = "response")
predicted_classes <- ifelse(predicted_probs > 0.5, 1, 0)

table(Predicted = predicted_classes, Actual = RF_data$response)

accuracy <- mean(predicted_classes == RF_data$response)
print(paste("Accuracy:", accuracy))

conf_mat_logit <- confusionMatrix(
  as.factor(predicted_classes),
  RF_data$response,
  positive = "1",
  mode = "everything"
)
print(conf_mat_logit)
print(conf_mat_logit$byClass["F1"])

roc_curve <- roc(RF_data$response, predicted_probs)
plot(roc_curve)
auc(roc_curve)

# =========================
# 5. Predict unknown events
# =========================
Unknown_data = Unknown_event_all %>%
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

# RF predictions
rf_predicted_probs = predict(rf, Unknown_data, type = "prob")
Unknown_data$rf_predicted_probs = rf_predicted_probs[, 2]

# Logistic predictions
unknown_predicted_probs_logit = predict(
  model_all_logistic,
  newdata = Unknown_data,
  type = "response",
  se.fit = TRUE
)

Unknown_data$predicted_probs = unknown_predicted_probs_logit$fit

# Compare predictions
print(t.test(Unknown_data$predicted_probs, Unknown_data$rf_predicted_probs, paired = TRUE))

logistic_class = ifelse(Unknown_data$predicted_probs > 0.5, 1, 0)
rf_class = ifelse(Unknown_data$rf_predicted_probs > 0.5, 1, 0)

print(table(Logistic = logistic_class, RF = rf_class))

consistency = mean(logistic_class == rf_class)
print(paste("Classification Consistency: ", round(consistency * 100, 2), "%", sep = ""))

Unknown_data$logistic_class = logistic_class
Unknown_data$rf_class = rf_class

write.csv(Unknown_data, "combined_unknown_predictions.csv", row.names = FALSE)

# =========================
# 7. Logistic effect plots
# =========================
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
