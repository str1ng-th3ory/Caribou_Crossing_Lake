# ----------------------------------------
# 05_prepare_modeling_dataset.R
# Purpose:
#   (1) Combine crossing_event, unknown_crossing_event, and circumnavigate_event
#   (2) Create event_id and derive reference-based variables
#   (3) Estimate lake_spent for focal events using ctmm occurrence
#   (4) Export combined event table and focal-row modeling table
# ----------------------------------------

library(dplyr)
library(lubridate)
library(sf)
library(terra)
library(raster)
library(ctmm)

# Load event tables
crossing = read.csv("crossing_event.csv", stringsAsFactors = FALSE)
unknown = read.csv("unknown_crossing_event.csv", stringsAsFactors = FALSE)
circum = read.csv("circumnavigate_event.csv", stringsAsFactors = FALSE)

# convert crossing_time
crossing$crossing_time = as.POSIXct(crossing$crossing_time, tz = "UTC")
unknown$crossing_time = as.POSIXct(unknown$crossing_time, tz = "UTC")
circum$crossing_time = as.POSIXct(circum$crossing_time, tz = "UTC")


# label event class
crossing$location = "Known Crossing Event"
unknown$location = "Unknown Crossing Event"
circum$location = "Known Circumnavigating Event"

# add source if helpful
crossing$source_table = "crossing_event"
unknown$source_table = "unknown_crossing_event"
circum$source_table = "circumnavigate_event"

# combine
all_event = bind_rows(crossing, unknown, circum)

# make sure type is numeric
all_event$type = as.numeric(all_event$type)

# Create event_id
# Each focal row (type == 0) starts a new event
# Do this within each source table to avoid cross-table mixing
# =========================
all_event = all_event %>%
  arrange(location, Year, ID, crossing_time) %>%
  group_by(location) %>%
  mutate(event_id_num = ceiling(row_number() / 7)) %>%
  # each focal event is expected to contribute 7 rows
  #	reference rows may be incomplete near seasonal boundaries, in which case users should verify grouping manually
  ungroup() %>%
  mutate(
    event_id = case_when(
      location == "Known Crossing Event" ~ paste0("known_crossing_", event_id_num),
      location == "Unknown Crossing Event" ~ paste0("unknown_crossing_", event_id_num),
      location == "Known Circumnavigating Event" ~ paste0("known_circumnavigate_", event_id_num)
    )) %>% dplyr::select(-event_id_num)


# =========================
all_event = all_event %>%
  dplyr::mutate(
    Year = as.numeric(Year),
    type = as.numeric(type),
    crossing_duration = as.numeric(crossing_duration),
    crossing_date = as.numeric(crossing_date),
    nearest_before_Lon = as.numeric(nearest_before_Lon),
    nearest_before_Lat = as.numeric(nearest_before_Lat),
    nearest_after_Lon = as.numeric(nearest_after_Lon),
    nearest_after_Lat = as.numeric(nearest_after_Lat),
    lake_width = as.numeric(lake_width),
    circumvent_distance = as.numeric(circumvent_distance),
    straight_distance = as.numeric(straight_distance),
    circumvent_speed = as.numeric(circumvent_speed),
    straight_speed = as.numeric(straight_speed),
    albedo_linear = as.numeric(albedo_linear),
    albedo_nearest_pixel = as.numeric(albedo_nearest_pixel),
    albedo_whole_lake = as.numeric(albedo_whole_lake)
  )

# =========================
# Derive reference-based variables
# These are assigned to focal rows only
# =========================
all_event = all_event %>%
  dplyr::group_by(event_id) %>%
  dplyr::mutate(
    average_speed_on_land = ifelse(
      type == 0,
      mean(straight_speed[type != 0], na.rm = TRUE),
      NA_real_
    ),
    max_speed_on_land = ifelse(
      type == 0,
      max(straight_speed[type != 0], na.rm = TRUE),
      NA_real_
    ),
    event_start = min(crossing_time, na.rm = TRUE),
    event_end = max(crossing_time, na.rm = TRUE)
  ) %>%
  dplyr::ungroup()

# focal rows only
all_cases = all_event %>%
  filter(type == 0)

# Load raw tracking data and lake polygon
# =========================
load("data/nwt_lakes.RData")
nwt_lakes = st_transform(nwt_lakes, 4326) %>% arrange(Time)

sfPolygons = st_read("data/Contwoyto_pure_water.shp", quiet = TRUE)
sfPolygons = st_transform(sfPolygons, 4326)

# attach sex from raw data if available
# =========================
if ("sex" %in% names(nwt_lakes)) {
  sex_lookup = nwt_lakes %>%
    st_drop_geometry() %>%
    distinct(ID, Year, sex)
  
  all_cases = all_cases %>%
    left_join(sex_lookup, by = c("ID", "Year"))
}

# =========================
# Function to estimate lake_spent for one focal event
# =========================
estimate_lake_spent_one_event = function(year, id, event_start, event_end, nwt_lakes, sfPolygons) {
  track_df = nwt_lakes %>%
    filter(Year == year, ID == id) %>%
    arrange(Time) %>%
    st_drop_geometry() %>%
    filter(Time >= event_start & Time <= event_end)
  
  if (nrow(track_df) < 2) return(NA_real_)
  
  test = track_df %>%
    mutate(
      x = Lon,
      y = Lat,
      time = Time,
      id = ID
    ) %>%
    as.data.frame()
  
  test$Time = as.POSIXct(test$Time, tz = "UTC")
  
  # ctmm telemetry object
  test_ctmm = as.telemetry(test)
  
  # fit occurrence model
  GUESS = ctmm.guess(test_ctmm, interactive = FALSE)
  FIT = ctmm.fit(test_ctmm, GUESS)
  test_occurrence = occurrence(test_ctmm, FIT)
  
  # convert occurrence to raster and reproject
  test_raster_terra = rast(raster(test_occurrence, DF = "PMF"))
  test_raster_wgs84 = terra::project(test_raster_terra, "EPSG:4326")
  
  # crop and mask by lake polygon
  cropped_raster = terra::crop(test_raster_wgs84, sfPolygons)
  masked_raster = terra::mask(cropped_raster, sfPolygons)
  
  # probability mass inside lake
  raster_sum = terra::global(masked_raster, fun = "sum", na.rm = TRUE)
  
  if (is.null(raster_sum) || nrow(raster_sum) == 0) return(NA_real_)
  
  # event duration in minutes
  time_difference_in_mins = as.numeric(difftime(test$Time[length(test$Time)], test$Time[1], units = "mins"))
  
  # expected time spent in lake
  time_in_lake = time_difference_in_mins * raster_sum[[1]]
  
  return(as.numeric(time_in_lake))
}

# =========================
# Estimate lake_spent for focal rows
# =========================
all_cases$lake_spent = NA_real_

for (i in seq_len(nrow(all_cases))) {
  all_cases$lake_spent[i] = estimate_lake_spent_one_event(
    year = all_cases$Year[i],
    id = all_cases$ID[i],
    event_start = all_cases$event_start[i],
    event_end = all_cases$event_end[i],
    nwt_lakes = nwt_lakes,
    sfPolygons = sfPolygons
  )
}

# =========================
# Join focal-level variables back to full event table
# =========================
focal_vars = all_cases %>%
  select(event_id, average_speed_on_land, max_speed_on_land, lake_spent, sex)

all_event_updated = all_event %>%
  select(-any_of(c("average_speed_on_land", "max_speed_on_land", "lake_spent", "sex"))) %>%
  left_join(focal_vars, by = "event_id")


# Build focal-row modeling table
# =========================
modeling_df = all_event_updated %>%
  filter(type == 0) %>%
  mutate(
    crossing_date = as.numeric(crossing_date),
    crossing_duration = as.numeric(crossing_duration),
    lake_width = as.numeric(lake_width),
    albedo_linear = as.numeric(albedo_linear),
    albedo_nearest_pixel = as.numeric(albedo_nearest_pixel),
    albedo_whole_lake = as.numeric(albedo_whole_lake),
    straight_distance = as.numeric(straight_distance),
    circumvent_distance = as.numeric(circumvent_distance),
    straight_speed = as.numeric(straight_speed),
    circumvent_speed = as.numeric(circumvent_speed),
    lake_spent_proportion = lake_spent / (crossing_duration * 1440),
    log_ratio_cir_ref = log(circumvent_speed / average_speed_on_land),
    log_ratio_dir_ref = log(straight_speed / average_speed_on_land),
    log_ratio_cir_dir = log(circumvent_speed / straight_speed)
  )


# Save outputs
# =========================
write.csv(all_event_updated, "all_event_combined.csv", row.names = FALSE)
write.csv(all_cases, "all_focal_events_with_lake_spent.csv", row.names = FALSE)
save(all_event, all_cases, file = "chapter1_cleaned_events.rda")

