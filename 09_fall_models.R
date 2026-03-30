# ----------------------------------------
# 09_fall_models.R
# Purpose:
#   Fit fall-only models for known events,
#   and predict unknown fall events using:
#   (1) Random Forest
#   (2) Binomial Logistic Regression
#
# Input:
#   - known_event_fall.csv
#   - unknown_event_fall.csv
#
# Output:
#   - fall_rf_test_predictions.csv
#   - fall_unknown_predictions.csv
#   - fall_model_input.csv
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
Known_event_fall = read.csv("known_event_fall.csv", stringsAsFactors = FALSE)
Unknown_event_fall = read.csv("unknown_event_fall.csv", stringsAsFactors = FALSE)

# =========================
# 2. Prepare fall known-event modeling table
# =========================
RF_fall = Known_event_fall %>%
  mutate(
    response = factor(response, levels = c(0, 1)),
    sex = as.factor(sex),
    season = as.factor(season)
  ) %>%
  dplyr::select(
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

write.csv(RF_fall, "fall_model_input.csv", row.names = FALSE)

# =========================
# 3. Fall Random Forest
# =========================
RF_data_fall_test = RF_fall %>%
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

set.seed(155)

train_indices = sample(1:nrow(RF_data_fall_test), 0.8 * nrow(RF_data_fall_test), replace = FALSE)
train_data = RF_data_fall_test[train_indices, ]
test_data = RF_data_fall_test[-train_indices, ]

class_counts = table(train_data$response)
total_samples = sum(class_counts)
num_classes = length(class_counts)
class_weights = total_samples / (num_classes * class_counts)

class_wt = c(
  `0` = class_weights[[1]],
  `1` = class_weights[[2]]
)

rf_fall = randomForest(
  response ~ .,
  data = train_data,
  mtry = 2,
  importance = TRUE,
  ntree = 1000,
  classwt = class_wt
)

print(rf_fall)
print(randomForest::importance(rf_fall))
varImpPlot(rf_fall)

# =========================
# 4. RF evaluation and variable importance
# =========================

test_pred_class = predict(rf_fall, test_data)
test_pred_prob = predict(rf_fall, test_data, type = "prob")

test_data$response = factor(test_data$response, levels = levels(test_pred_class))

conf_matrix = confusionMatrix(
  test_pred_class,
  test_data$response,
  mode = "everything",
  positive = "1"
)
print(conf_matrix)

perf2 = prediction(test_pred_prob[, 2], test_data$response)
auc2 = performance(perf2, "auc")
cat("Fall RF AUC:", auc2@y.values[[1]], "\n")

roc_curve_rf = roc(test_data$response, test_pred_prob[, 2])
plot(roc_curve_rf)
print(auc(roc_curve_rf))

fall_rf_test_output = test_data
fall_rf_test_output$pred_class = test_pred_class
fall_rf_test_output$pred_prob = test_pred_prob[, 2]

write.csv(fall_rf_test_output, "fall_rf_test_predictions.csv", row.names = FALSE)

# =========================
# Additional RF evaluation
# =========================

# train AUC
pred1 = predict(rf_fall, train_data, type = "prob")
perf_train = prediction(pred1[, 2], train_data$response)
auc_train = performance(perf_train, "auc")
cat("Fall RF train AUC:", auc_train@y.values[[1]], "\n")

# train ROC
pred_train_roc = performance(perf_train, "tpr", "fpr")
plot(pred_train_roc, main = "ROC Curve for Fall Random Forest (Train)", col = 2, lwd = 2)
abline(a = 0, b = 1, lwd = 2, lty = 2, col = "gray")

# test class prediction
predictions = predict(rf_fall, test_data)

# test confusion matrix
confusion_matrix = table(test_data$response, predictions)
print(confusion_matrix)

# test accuracy
accuracy = sum(diag(confusion_matrix)) / sum(confusion_matrix)
cat("Fall RF test accuracy:", accuracy, "\n")

# test precision and recall
TN = confusion_matrix[1, 1]
FP = confusion_matrix[1, 2]
FN = confusion_matrix[2, 1]
TP = confusion_matrix[2, 2]

precision = TP / (TP + FP)
recall = TP / (TP + FN)

cat('Fall RF test precision:', precision, "\n")
cat("Fall RF test recall:", recall, "\n")

# test AUC
pred2 = predict(rf_fall, test_data, type = "prob")
perf_test = prediction(pred2[, 2], test_data$response)
auc_test = performance(perf_test, "auc")
cat("Fall RF test AUC:", auc_test@y.values[[1]], "\n")

# test ROC
pred_test_roc = performance(perf_test, "tpr", "fpr")
plot(pred_test_roc, main = "ROC Curve for Fall Random Forest (Test)", col = 2, lwd = 2)
abline(a = 0, b = 1, lwd = 2, lty = 2, col = "gray")

# =========================
# plot
# =========================

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

# =========================
# 4. Fall logistic regression
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

RF_fall$response = droplevels(RF_fall$response)

formula_full = as.formula(paste("response ~", paste(all_predictors, collapse = " + ")))

full_model = glm(
  formula_full,
  family = "binomial",
  data = RF_fall,
  na.action = na.fail
)

print(summary(full_model))

#lasso regression
x = model.matrix(formula_full, RF_fall)[, -1]
y = RF_fall$response

cv_fit = cv.glmnet(x, y, family = "binomial", alpha = 1)
best_lambda = cv_fit$lambda.min
lasso_model = glmnet(x, y, family = "binomial", alpha = 1, lambda = best_lambda)
print(coef(lasso_model))

#ridge regression
cv_fit_ridge = cv.glmnet(x, y, family = "binomial", alpha = 0)
ridge_model = glmnet(x, y, family = "binomial", alpha = 0, lambda = cv_fit_ridge$lambda.min)
print(coef(ridge_model))

#remove variable
adjust_predictors = c(
  "lake_spent_proportion",
  "albedo_whole_lake",
  "log_ratio_cir_ref"
)

formula_adjust = as.formula(paste("response ~", paste(adjust_predictors, collapse = " + ")))
adjust_model = glm(formula_adjust, family = "binomial", data = RF_fall)
print(vif(adjust_model))


# =========================
#  Check whether quadratic terms are needed
# using binned observed crossing proportions
# =========================
k = 12

tmp = RF_fall %>%
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
  labs(x = "Time Spent", y = "Observed crossing proportion")

tmp = RF_fall %>%
  mutate(bin = ntile(albedo_whole_lake, k)) %>%
  group_by(bin) %>%
  summarise(
    albedo_whole_lake_mid = median(albedo_whole_lake, na.rm = TRUE),
    p_cross = mean(response == 1, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(tmp, aes(albedo_whole_lake_mid, p_cross)) +
  geom_point() +
  geom_line() +
  labs(x = "Albedo", y = "Observed crossing proportion")

tmp = RF_fall %>%
  mutate(bin = ntile(log_ratio_dir_ref, k)) %>%
  group_by(bin) %>%
  summarise(
    log_ratio_dir_ref_mid = median(log_ratio_dir_ref, na.rm = TRUE),
    p_cross = mean(response == 1, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(tmp, aes(log_ratio_dir_ref_mid, p_cross)) +
  geom_point() +
  geom_line() +
  labs(x = "log_ratio_dir_ref", y = "Observed crossing proportion")
# If the binned relationships appear approximately monotonic,
# quadratic terms are not retained in the candidate models.
# Visual inspection suggested no clear U-shaped or inverted-U pattern,
# so quadratic terms were not included in the spring candidate models.


#adjust_model <- glm(adjust_model, family = "binomial", data = RF_fall, na.action = na.fail,
#control = list(maxit = 100))
#dredge_results=dredge(adjust_model)
#only converage from 5 to 8 model 
#summary(get.models(dredge_results, 7)[[1]])

# after excluding models with convergence problems and quasi-complete separation,
# the top supported model remained:
logistic_model_fall = glm(
  response ~ log_ratio_cir_ref,
  family = "binomial",
  data = RF_fall
)

print(summary(logistic_model_fall))

predicted_probs = predict(logistic_model_fall, type = "response")
predicted_classes = ifelse(predicted_probs > 0.5, 1, 0)

print(table(Predicted = predicted_classes, Actual = RF_fall$response))
print(mean(predicted_classes == RF_fall$response))

print(
  confusionMatrix(
    as.factor(predicted_classes),
    RF_fall$response,
    mode = "everything",
    positive = "1"
  )
)

roc_curve_logit = roc(RF_fall$response, predicted_probs)
plot(roc_curve_logit)
print(auc(roc_curve_logit))


# =========================
# 6. Fall logistic effect plot
# =========================
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

# =========================
# 7. Predict unknown fall events
# =========================
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

print(t.test(unknown_data_fall$predicted_probs, unknown_data_fall$rf_predicted_probs, paired = TRUE))

unknown_data_fall$logistic_class = ifelse(unknown_data_fall$predicted_probs > 0.5, 1, 0)
unknown_data_fall$rf_class = ifelse(unknown_data_fall$rf_predicted_probs > 0.5, 1, 0)

print(table(Logistic = unknown_data_fall$logistic_class, RF = unknown_data_fall$rf_class))
print(mean(unknown_data_fall$logistic_class == unknown_data_fall$rf_class))

consistency <- mean(unknown_data_fall$logistic_class == unknown_data_fall$rf_class)
print(paste("Classification Consistency: ", round(consistency * 100, 2), "%", sep = ""))


write.csv(unknown_data_fall, "fall_unknown_predictions.csv", row.names = FALSE)


# optional confusion-style plot for agreement
confusion_matrix = table(
  Logistic = unknown_data_fall$logistic_class,
  RF = unknown_data_fall$rf_class
)

confusion_df = as.data.frame(as.table(confusion_matrix))
colnames(confusion_df) = c("Logistic", "RF", "Count")

ggplot(confusion_df, aes(x = RF, y = Logistic)) +
  geom_tile(aes(fill = Count), color = "black") +
  geom_text(aes(label = Count), size = 6) +
  scale_fill_gradient(low = "white", high = "blue") +
  labs(
    x = "Random Forest Model",
    y = "Logistic Model",
    title = "Agreement between Logistic and RF Predictions for Unknown Fall Events"
  ) +
  theme_minimal()
