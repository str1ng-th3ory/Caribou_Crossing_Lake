# ----------------------------------------
# Figure_4_known_event_maps.R
# Purpose:
#   Plot spatial distributions of known crossing and
#   known circumnavigating events in spring and fall.
#
# Panels:
#   A = Spring crossing, event steps (-3 to +3)
#   B = Spring crossing, full trajectories
#   C = Spring circumnavigating, event steps (-3 to +3)
#   D = Spring circumnavigating, full trajectories
#   E = Fall crossing, event steps (-3 to +3)
#   F = Fall crossing, full trajectories
#   G = Fall circumnavigating, event steps (-3 to +3)
#   H = Fall circumnavigating, full trajectories
# ----------------------------------------

library(dplyr)
library(lubridate)
library(sf)
library(terra)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)

# =========================
# 1. Load data
# =========================

# event table generated from your new workflow
# use whichever one you saved with event_id / event_start / event_end
all_event = read.csv("all_event_combined.csv", stringsAsFactors = FALSE)

# raw movement data
load("data/nwt_lakes.RData")

# lake polygon
lake = st_read("data/Contwoyto_pure_water.shp", quiet = TRUE)

# =========================
# 2. Standardize data types
# =========================
all_event = all_event %>%
  mutate(
    Year = as.numeric(Year),
    crossing_time = as.POSIXct(crossing_time, tz = "UTC"),
    event_start = as.POSIXct(event_start, tz = "UTC"),
    event_end = as.POSIXct(event_end, tz = "UTC"),
    nearest_before_Lon = as.numeric(nearest_before_Lon),
    nearest_before_Lat = as.numeric(nearest_before_Lat),
    nearest_after_Lon = as.numeric(nearest_after_Lon),
    nearest_after_Lat = as.numeric(nearest_after_Lat),
    type = as.numeric(type),
    season = as.character(season),
    location = as.character(location),
    event_id = as.character(event_id)
  )

nwt_lakes = st_transform(nwt_lakes, 4326) %>%
  arrange(Time)

lake = st_transform(lake, 4326)

# =========================
# 3. Keep only known events
# =========================
known_event_steps = all_event %>%
  filter(location %in% c("Known Crossing Event", "Known Circumnavigating Event"))

# focal rows only, used to recover full trajectories
known_event_focal = known_event_steps %>%
  filter(type == 0)

# =========================
# 4. Build full-trajectory table for each event
# =========================
build_full_trajectory_df = function(focal_df, track_df) {
  out_list = vector("list", nrow(focal_df))
  
  for (i in seq_len(nrow(focal_df))) {
    year_i = focal_df$Year[i]
    id_i = focal_df$ID[i]
    start_i = focal_df$event_start[i]
    end_i = focal_df$event_end[i]
    event_id_i = focal_df$event_id[i]
    
    tmp = track_df %>%
      filter(Year == year_i, ID == id_i) %>%
      arrange(Time) %>%
      filter(Time >= start_i & Time <= end_i) %>%
      st_drop_geometry() %>%
      mutate(event_id = event_id_i)
    
    out_list[[i]] = tmp
  }
  
  bind_rows(out_list)
}

known_full_traj = build_full_trajectory_df(
  focal_df = known_event_focal,
  track_df = nwt_lakes
)


# =========================
# 5. Split data for the 8 panels
# =========================

# event-step panels
steps_spring_cross = known_event_steps %>%
  filter(season == "spring", location == "Known Crossing Event")

steps_spring_circ = known_event_steps %>%
  filter(season == "spring", location == "Known Circumnavigating Event")

steps_fall_cross = known_event_steps %>%
  filter(season == "fall", location == "Known Crossing Event")

steps_fall_circ = known_event_steps %>%
  filter(season == "fall", location == "Known Circumnavigating Event")

# full-trajectory panels
traj_spring_cross = known_full_traj %>%
  semi_join(
    known_event_focal %>%
      filter(season == "spring", location == "Known Crossing Event") %>%
      select(event_id),
    by = "event_id"
  )

traj_spring_circ = known_full_traj %>%
  semi_join(
    known_event_focal %>%
      filter(season == "spring", location == "Known Circumnavigating Event") %>%
      select(event_id),
    by = "event_id"
  )

traj_fall_cross = known_full_traj %>%
  semi_join(
    known_event_focal %>%
      filter(season == "fall", location == "Known Crossing Event") %>%
      select(event_id),
    by = "event_id"
  )

traj_fall_circ = known_full_traj %>%
  semi_join(
    known_event_focal %>%
      filter(season == "fall", location == "Known Circumnavigating Event") %>%
      select(event_id),
    by = "event_id"
  )

# =========================
# 6. Plot helpers
# =========================
# event-step panels:
# show all step points from -3 to +3,
# and connect only the focal local segment (type 0 to 1)
plot_event_steps = function(step_df, panel_title = NULL) {
  ggplot(lake) +
    geom_sf(fill = "black", alpha = 0.5, lwd = 0) +
    geom_path(
      data = step_df %>% filter(type %in% c(0, 1)),
      aes(
        x = nearest_before_Lon,
        y = nearest_before_Lat,
        group = event_id,
        color = type
      )
    ) +
    geom_point(
      data = step_df,
      aes(
        x = nearest_before_Lon,
        y = nearest_before_Lat,
        color = type
      )
    ) +
    scale_colour_gradientn(
      colours = rev(brewer.pal(n = 11, name = "RdYlBu"))
    ) +
    theme_dark() +
    labs(
      title = panel_title,
      x = "Lon",
      y = "Lat",
      color = "Water Crossing Steps"
    )
}

# full-trajectory panels:
# show the complete trajectory segment during the event window
plot_full_trajectory = function(traj_df, panel_title = NULL, legend_title = legend_title) {
  ggplot(lake) +
    geom_sf(fill = "black", alpha = 0.5, lwd = 0) +
    geom_path(
      data = traj_df,
      aes(x = Lon, y = Lat, group = event_id, color = "Path")
    ) +
    geom_point(
      data = traj_df,
      aes(x = Lon, y = Lat, group = event_id, color = "GPS Location")
    ) +
    scale_color_manual(
      values = c(
        "Path" = "grey",
        "GPS Location" = "white"
      )
    ) +
    theme_dark() +
    labs(
      title = panel_title,
      x = "Lon",
      y = "Lat",
      color = legend_title
    )
}


# =========================
# 7. Build the 8 panels
# =========================

pA = plot_event_steps(
  steps_spring_cross,
  panel_title = "Spring crossing, event steps"
)

pB = plot_full_trajectory(
  traj_spring_cross,
  panel_title = "Spring crossing, full trajectories",
  legend_title = "Crossing Events"
)

pC = plot_event_steps(
  steps_spring_circ,
  panel_title = "Spring circumnavigating, event steps"
)

pD = plot_full_trajectory(
  traj_spring_circ,
  panel_title = "Spring circumnavigating, full trajectories",
  legend_title = "Circumnavigating Events"
)

pE = plot_event_steps(
  steps_fall_cross,
  panel_title = "Fall crossing, event steps"
)

pF = plot_full_trajectory(
  traj_fall_cross,
  panel_title = "Fall crossing, full trajectories",
  legend_title = "Crossing Events"
)

pG = plot_event_steps(
  steps_fall_circ,
  panel_title = "Fall circumnavigating, event steps"
)

pH = plot_full_trajectory(
  traj_fall_circ,
  panel_title = "Fall circumnavigating, full trajectories",
  legend_title = "Circumnavigating Events"
)

# =========================
# 8. Arrange final figure
# =========================
figure4 = ggarrange(
  pA, pB,
  pC, pD,
  pE, pF,
  pG, pH,
  labels = c("A", "B", "C", "D", "E", "F", "G", "H"),
  ncol = 2,
  nrow = 4,
  common.legend = FALSE
)

print(figure4)

# =========================
# 9. Save
# =========================
ggsave(
  filename = "Figure4_known_events_spatial_distribution.png",
  plot = figure4,
  width = 12,
  height = 20,
  dpi = 300
)
