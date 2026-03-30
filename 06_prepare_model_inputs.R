# ----------------------------------------
# 06_prepare_model_inputs.R
# Purpose:
#   Prepare cleaned focal-event datasets for downstream modeling.
#
# Input:
#   - all_focal_events_with_lake_spent.csv
#
# Output:
#   - feature_df.csv
#   - known_event_all.csv
#   - unknown_event_all.csv
#   - known_event_spring.csv
#   - known_event_fall.csv
#   - unknown_event_spring.csv
#   - unknown_event_fall.csv
# ----------------------------------------

library(dplyr)
library(readr)

# 1. Load focal-event table
# =========================
all_cases = read.csv("all_focal_events_with_lake_spent.csv", stringsAsFactors = FALSE)


# 2. Standardize data types
# =========================
all_cases = all_cases %>%
  rename(
    location = location,
    #albedo_linear = alebdo_linear
  ) %>%
  mutate(
    Year = as.numeric(Year),
    ID = as.character(ID),
    location = as.character(location),
    season = as.character(season),
    sex = as.character(sex),
    crossing_time = as.POSIXct(crossing_time, tz = "UTC"),
    crossing_duration = as.numeric(crossing_duration),
    crossing_date = as.numeric(crossing_date),
    lake_width = as.numeric(lake_width),
    circumvent_distance = as.numeric(circumvent_distance),
    straight_distance = as.numeric(straight_distance),
    circumvent_speed = as.numeric(circumvent_speed),
    straight_speed = as.numeric(straight_speed),
    albedo_linear = as.numeric(albedo_linear),
    albedo_nearest_pixel = as.numeric(albedo_nearest_pixel),
    albedo_whole_lake = as.numeric(albedo_whole_lake),
    average_speed_on_land = as.numeric(average_speed_on_land),
    max_speed_on_land = as.numeric(max_speed_on_land),
    lake_spent = as.numeric(lake_spent)
  )

# 3. Remove known problematic rows if needed
# =========================
# Example carried over from the original workflow:
# remove 2016 bat143 2016-07-16 12:00:00 because of NA in albedo_nearest_pixel
# all_cases = all_cases %>%
#   filter(
#     !(Year == 2016 &
#         ID == "bat143" &
#         crossing_time == as.POSIXct("2016-07-16 12:00:00", tz = "UTC")) )

# 4. Build feature-engineered table
# =========================
feature_df = all_cases %>%
  mutate(
    response = case_when(
      location == "Known Crossing Event" ~ 1,
      location == "Known Circumnavigating Event" ~ 0,
      TRUE ~ NA_real_
    ),
    response = as.factor(response),
    sex = as.factor(sex),
    season = as.factor(season),
    lake_spent_proportion = lake_spent / (crossing_duration * 1440),
    log_ratio_cir_ref = log(circumvent_speed / average_speed_on_land),
    log_ratio_dir_ref = log(straight_speed / average_speed_on_land),
    log_ratio_cir_dir = log(circumvent_speed / straight_speed)
  )

# 5. Split into known / unknown
# =========================
Known_event_all = feature_df %>%
  filter(location %in% c("Known Crossing Event", "Known Circumnavigating Event")) %>%
  filter(!is.na(average_speed_on_land))

Unknown_event_all = feature_df %>%
  filter(location == "Unknown Crossing Event") %>%
  filter(!is.na(average_speed_on_land))

# Optional: keep only rows with complete data for unknown-event prediction
# Unknown_event_all_complete = Unknown_event_all %>%
#   filter(rowSums(is.na(.)) == 0)

# 6. Split by season
# =========================
Known_event_spring = Known_event_all %>% filter(season == "spring")
Known_event_fall   = Known_event_all %>% filter(season == "fall")

Unknown_event_spring = Unknown_event_all %>% filter(season == "spring")
Unknown_event_fall   = Unknown_event_all %>% filter(season == "fall")

# 7. Quick checks
# =========================
cat("Total focal events:", nrow(all_cases), "\n")
cat("Known events:", nrow(Known_event_all), "\n")
cat("Unknown events:", nrow(Unknown_event_all), "\n")
cat("Known spring:", nrow(Known_event_spring), "\n")
cat("Known fall:", nrow(Known_event_fall), "\n")
cat("Unknown spring:", nrow(Unknown_event_spring), "\n")
cat("Unknown fall:", nrow(Unknown_event_fall), "\n")

# 9. Save outputs
# =========================
write.csv(feature_df, "feature_df.csv", row.names = FALSE)

write.csv(Known_event_all, "known_event_all.csv", row.names = FALSE)
write.csv(Unknown_event_all, "unknown_event_all.csv", row.names = FALSE)

write.csv(Known_event_spring, "known_event_spring.csv", row.names = FALSE)
write.csv(Known_event_fall, "known_event_fall.csv", row.names = FALSE)

write.csv(Unknown_event_spring, "unknown_event_spring.csv", row.names = FALSE)
write.csv(Unknown_event_fall, "unknown_event_fall.csv", row.names = FALSE)
