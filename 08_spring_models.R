# ----------------------------------------
# 08_spring_models.R
# Purpose:
#   Fit spring-only models for known events,
#   and predict unknown spring events using:
#   (1) Random Forest
#   (2) Binomial Logistic Regression
#
# Input:
#   - known_event_spring.csv
#   - unknown_event_spring.csv
#
# Output:
#   - spring_rf_test_predictions.csv
#   - spring_unknown_predictions.csv
#   - spring_model_input.csv
# ----------------------------------------

library(dplyr)
library(randomForest)
library(caret)
library(ROCR)
library(pROC)
library(MuMIn)
library(car)
library(glmnet)
library(ggplot2)
library(svglite)

# =========================
# 1. Load data
# =========================
Known_event_spring = read.csv("known_event_spring.csv", stringsAsFactors = FALSE)
Unknown_event_spring = read.csv("unknown_event_spring.csv", stringsAsFactors = FALSE)

# =========================
# 2. Prepare spring known-event modeling table
# =========================
RF_spring = Known_event_spring %>%
  mutate(
    response = factor(response, levels = c(0, 1)),
    sex = as.factor(sex),
    season = as.factor(season)
  ) %>%
  select(
    response,
    crossing_duration,
    lake_spent_proportion,
    albedo_nearest_pixel,
    albedo_linear,
    albedo_whole_lake,
    lake_width,
    log_ratio_cir_ref,
    log_ratio_dir_ref,
    log_ratio_cir_dir,
    lake_spent,
    sex,
    season,
    crossing_date
  ) %>%
  filter(rowSums(is.na(.)) == 0)

write.csv(RF_spring, "spring_model_input.csv", row.names = FALSE)

# =========================
# 3. Spring Random Forest
# =========================
RF_data_spring_test = RF_spring %>%
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

set.seed(177)

train_indices = sample(1:nrow(RF_data_spring_test), 0.8 * nrow(RF_data_spring_test), replace = FALSE)
train_data = RF_data_spring_test[train_indices, ]
test_data = RF_data_spring_test[-train_indices, ]

class_counts = table(train_data$response)
total_samples = sum(class_counts)
num_classes = length(class_counts)
class_weights = total_samples / (num_classes * class_counts)

class_wt = c(
  `0` = class_weights[[1]],
  `1` = class_weights[[2]]
)

rf_spring = randomForest(
  response ~ .,
  data = train_data,
  mtry = 2,
  importance = TRUE,
  ntree = 1000,
  classwt = class_wt
)

print(rf_spring)
print(randomForest::importance(rf_spring))
varImpPlot(rf_spring)

# =========================
# 4. RF evaluation and variable importance
# =========================

test_pred_class = predict(rf_spring, test_data)
test_pred_prob = predict(rf_spring, test_data, type = "prob")

test_data$response = factor(test_data$response, levels = levels(test_pred_class))

conf_matrix = confusionMatrix(
  test_pred_class,
  test_data$response,
  mode = "everything",
  positive = "1"
)
print(conf_matrix)

perf_test = prediction(test_pred_prob[, 2], test_data$response)
auc_test = performance(perf_test, "auc")
cat("Spring RF test AUC:", auc_test@y.values[[1]], "\n")

roc_curve_rf = roc(test_data$response, test_pred_prob[, 2])
plot(roc_curve_rf)
print(auc(roc_curve_rf))

spring_rf_test_output = test_data
spring_rf_test_output$pred_class = test_pred_class
spring_rf_test_output$pred_prob = test_pred_prob[, 2]
write.csv(spring_rf_test_output, "spring_rf_test_predictions.csv", row.names = FALSE)

# =========================
# Additional RF evaluation
# =========================

# train AUC
pred1 = predict(rf_spring, train_data, type = "prob")
perf_train = prediction(pred1[, 2], train_data$response)
auc_train = performance(perf_train, "auc")
cat("Spring RF train AUC:", auc_train@y.values[[1]], "\n")

# train ROC
pred_train_roc = performance(perf_train, "tpr", "fpr")
plot(pred_train_roc, main = "ROC Curve for Spring Random Forest (Train)", col = 2, lwd = 2)
abline(a = 0, b = 1, lwd = 2, lty = 2, col = "gray")

# test class prediction
predictions = predict(rf_spring, test_data)

# test confusion matrix
confusion_matrix = table(test_data$response, predictions)
print(confusion_matrix)

# test accuracy
accuracy = sum(diag(confusion_matrix)) / sum(confusion_matrix)
cat("Spring RF test accuracy:", accuracy, "\n")

# test precision and recall
TN = confusion_matrix[1, 1]
FP = confusion_matrix[1, 2]
FN = confusion_matrix[2, 1]
TP = confusion_matrix[2, 2]

precision = TP / (TP + FP)
recall = TP / (TP + FN)

cat("Spring RF test precision:", precision, "\n")
cat("Spring RF test recall:", recall, "\n")

# test AUC
pred2 = predict(rf_spring, test_data, type = "prob")
perf_test = prediction(pred2[, 2], test_data$response)
auc_test = performance(perf_test, "auc")
cat("Spring RF test AUC:", auc_test@y.values[[1]], "\n")

# test ROC
pred_test_roc = performance(perf_test, "tpr", "fpr")
plot(pred_test_roc, main = "ROC Curve for Spring Random Forest (Test)", col = 2, lwd = 2)
abline(a = 0, b = 1, lwd = 2, lty = 2, col = "gray")


# =========================
# plot
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

p_rf = ggplot(importance_data, aes(x = reorder(Variable, MeanDecreaseGini), y = MeanDecreaseGini)) +
  geom_point(size = 3, color = "blue") +
  labs(x = "Variables", y = "Importance (Gini Decrease)") +
  coord_flip() +
  theme_classic()

print(p_rf)

# =========================
# 4. Spring logistic regression
# =========================
all_predictors = c(
  "lake_spent_proportion",
  "albedo_nearest_pixel",
  "albedo_linear",
  "albedo_whole_lake",
  "lake_width",
  "log_ratio_cir_ref",
  "log_ratio_dir_ref",
  "log_ratio_cir_dir",
  "sex"
)

RF_spring$response = droplevels(RF_spring$response)

formula_full = as.formula(paste("response ~", paste(all_predictors, collapse = " + ")))

full_model = glm(
  formula_full,
  family = "binomial",
  data = RF_spring,
  na.action = na.fail,
  control = list(maxit = 100)
)

print(summary(full_model))
print(vif(full_model))
print(alias(full_model))

# LASSO
x = model.matrix(formula_full, RF_spring)[, -1]
y = RF_spring$response

cv_fit = cv.glmnet(x, y, family = "binomial", alpha = 1)
best_lambda = cv_fit$lambda.min
lasso_model = glmnet(x, y, family = "binomial", alpha = 1, lambda = best_lambda)
print(coef(lasso_model))


# =========================
#  Check whether quadratic terms are needed
# using binned observed crossing proportions
# =========================
k = 12

tmp = RF_spring %>%
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

tmp = RF_spring %>%
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
  labs(x = "lake_spent_proportion", y = "Observed crossing proportion")

tmp = RF_spring %>%
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

# If the binned relationships appear approximately monotonic,
# quadratic terms are not retained in the candidate models.
# Visual inspection suggested no clear U-shaped or inverted-U pattern,
# so quadratic terms were not included in the spring candidate models.


# reduced model for the multi-collinearity #check VIF
adjust_predictors = c(
  "lake_spent_proportion",
  "albedo_linear",
  "log_ratio_cir_ref"
)

formula_adjust = as.formula(paste("response ~", paste(adjust_predictors, collapse = " + ")))
adjust_model = glm(
  formula_adjust,
  family = "binomial",
  data = RF_spring,
  na.action = na.fail
)

print(vif(adjust_model))

dredge_results_spring = dredge(adjust_model)
print(dredge_results_spring)


# after excluding models with convergence problems and quasi-complete separation,
# the top supported model remained:
logistic_model_spring = glm(
  response ~ albedo_linear,
  family = "binomial",
  data = RF_spring
)

print(summary(logistic_model_spring))

predicted_probs = predict(logistic_model_spring, type = "response")
predicted_classes = ifelse(predicted_probs > 0.5, 1, 0)

print(table(Predicted = predicted_classes, Actual = RF_spring$response))
print(mean(predicted_classes == RF_spring$response))

print(
  confusionMatrix(
    as.factor(predicted_classes),
    RF_spring$response,
    mode = "everything",
    positive = "1"
  )
)

roc_curve_logit = roc(RF_spring$response, predicted_probs)
plot(roc_curve_logit)
print(auc(roc_curve_logit))

# =========================
# 5. Spring logistic effect plot
# =========================
threshold_data = data.frame(
  albedo_linear = RF_spring$albedo_linear,
  predicted_probs = predicted_probs,
  lower_bound = NA_real_,
  upper_bound = NA_real_,
  location = case_when(
    RF_spring$response == 1 ~ "Known Crossing Event",
    RF_spring$response == 0 ~ "Known Circumnavigating Event"
  )
)

albedo_grid = seq(min(RF_spring$albedo_linear), max(RF_spring$albedo_linear), length = 100)
new_data_ratio = data.frame(albedo_linear = albedo_grid)

predicted_probs_with_se = predict(
  logistic_model_spring,
  newdata = new_data_ratio,
  type = "response",
  se.fit = TRUE
)

new_data_ratio$lower_bound = pmax(0, predicted_probs_with_se$fit - 1.96 * predicted_probs_with_se$se.fit)
new_data_ratio$upper_bound = pmin(1, predicted_probs_with_se$fit + 1.96 * predicted_probs_with_se$se.fit)
new_data_ratio$predicted_probs = predicted_probs_with_se$fit

p_logit = ggplot(new_data_ratio, aes(x = albedo_linear, y = predicted_probs)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower_bound, ymax = upper_bound), alpha = 0.2) +
  geom_point(data = threshold_data, aes(x = albedo_linear, y = predicted_probs, color = location), size = 2) +
  labs(
    x = "Albedo (the set of pixels along potential crossing paths)",
    y = "Predicted Probability of Crossing",
    color = "Event"
  ) +
  theme_classic()

print(p_logit)

# =========================
# 6. Predict unknown spring events
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
levels(unknown_data_spring$season)

unknown_predicted_probs_spring_rf = predict(rf_spring, unknown_data_spring, type = "prob")
unknown_data_spring$rf_predicted_probs = unknown_predicted_probs_spring_rf[, 2]

unknown_predicted_probs_spring_logit = predict(
  logistic_model_spring,
  newdata = unknown_data_spring,
  type = "response",
  se.fit = TRUE
)

unknown_data_spring$predicted_probs = unknown_predicted_probs_spring_logit$fit
unknown_data_spring$lower_bound = pmax(0, unknown_predicted_probs_spring_logit$fit - 1.96 * unknown_predicted_probs_spring_logit$se.fit)
unknown_data_spring$upper_bound = pmin(1, unknown_predicted_probs_spring_logit$fit + 1.96 * unknown_predicted_probs_spring_logit$se.fit)

p_unknown = ggplot(new_data_ratio, aes(x = albedo_linear, y = predicted_probs)) +
  geom_line(color = "blue", linewidth = 0.6) +
  geom_ribbon(aes(ymin = lower_bound, ymax = upper_bound), alpha = 0.2) +
  geom_point(data = threshold_data, aes(x = albedo_linear, y = predicted_probs, color = location), size = 1.5) +
  geom_point(data = unknown_data_spring, aes(x = albedo_linear, y = predicted_probs, color = "Unknown Event (Binomial Logistic Model)"), size = 1.5) +
  geom_point(data = unknown_data_spring, aes(x = albedo_linear, y = rf_predicted_probs, color = "Unknown Event (Random Forest Model)"), size = 1.5) +
  labs(
    x = "APR Along the Potential Crossing Path",
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

print(t.test(unknown_data_spring$predicted_probs, unknown_data_spring$rf_predicted_probs, paired = TRUE))

unknown_data_spring$logistic_class = ifelse(unknown_data_spring$predicted_probs > 0.5, 1, 0)
unknown_data_spring$rf_class = ifelse(unknown_data_spring$rf_predicted_probs > 0.5, 1, 0)

print(table(Logistic = unknown_data_spring$logistic_class, RF = unknown_data_spring$rf_class))
print(mean(unknown_data_spring$logistic_class == unknown_data_spring$rf_class))

write.csv(unknown_data_spring, "spring_unknown_predictions.csv", row.names = FALSE)



# t-test
t.test(unknown_data_spring$predicted_probs, unknown_data_spring$rf_predicted_probs, paired = TRUE)

# Compare predictions
unknown_data_spring$logistic_class <- ifelse(unknown_data_spring$predicted_probs > 0.5, 1, 0)
unknown_data_spring$rf_class <- ifelse(unknown_data_spring$rf_predicted_probs > 0.5, 1, 0)


table(Logistic = unknown_data_spring$logistic_class, RF = unknown_data_spring$rf_class)


consistency <- mean(unknown_data_spring$logistic_class == unknown_data_spring$rf_class)
print(paste("Classification Consistency: ", round(consistency * 100, 2), "%", sep = ""))

