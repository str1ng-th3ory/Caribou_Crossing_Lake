# ----------------------------------------
# Figure5_unknown_event_trajectories.R
# Purpose:
#   Plot spatial distribution of unknown transit events
#   during spring and fall migration.
#
# Panels:
#   A. Spring migration unknown events
#   B. Fall migration unknown events
#
# Required input files:
#   - unknown_event_spring.csv
#   - unknown_event_fall.csv
#   - Contwoyto_pure_water.shp
# ----------------------------------------

library(dplyr)
library(sf)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)

# =========================
# 1. Load data
# =========================

Unknown_event_spring = read.csv("unknown_event_spring.csv", stringsAsFactors = FALSE)
Unknown_event_fall = read.csv("unknown_event_fall.csv", stringsAsFactors = FALSE)

polys1_sf = st_read("Contwoyto_pure_water.shp", quiet = TRUE)
polys1_sf = st_transform(polys1_sf, 4326)
lake = polys1_sf

# =========================
# 2. Prepare unknown-event tables
# =========================
unknown_spring_plot = Unknown_event_spring %>%
  mutate(
    across(
      c(
        Year, crossing_duration, crossing_date,
        nearest_before_Lon, nearest_before_Lat,
        nearest_after_Lon, nearest_after_Lat,
        lake_width, type
      ),
      as.numeric
    )
  ) %>%
  filter(location == "Unknown Crossing Event")

unknown_fall_plot = Unknown_event_fall %>%
  mutate(
    across(
      c(
        Year, crossing_duration, crossing_date,
        nearest_before_Lon, nearest_before_Lat,
        nearest_after_Lon, nearest_after_Lat,
        lake_width, type
      ),
      as.numeric
    )
  ) %>%
  filter(location == "Unknown Crossing Event")

# =========================
# 3. Plot spring unknown events
# =========================
p_fig5A = ggplot(lake) +
  geom_sf(fill = "black", alpha = 0.5, lwd = 0) +
  geom_path(
    data = unknown_spring_plot %>% filter(type == 0 | type == 1),
    aes(
      x = nearest_before_Lon,
      y = nearest_before_Lat,
      group = event_id,
      color = type
    )
  ) +
  geom_point(
    data = unknown_spring_plot,
    aes(
      x = nearest_before_Lon,
      y = nearest_before_Lat,
      color = type
    )
  ) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu"))) +
  theme_dark() +
  labs(
    color = "Water Crossing Steps",
    x = "Longitude",
    y = "Latitude"
  )

print(p_fig5A)

# =========================
# 4. Plot fall unknown events
# =========================
p_fig5B = ggplot(lake) +
  geom_sf(fill = "black", alpha = 0.5, lwd = 0) +
  geom_path(
    data = unknown_fall_plot %>% filter(type == 0 | type == 1),
    aes(
      x = nearest_before_Lon,
      y = nearest_before_Lat,
      group = event_id,
      color = type
    )
  ) +
  geom_point(
    data = unknown_fall_plot,
    aes(
      x = nearest_before_Lon,
      y = nearest_before_Lat,
      color = type
    )
  ) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu"))) +
  theme_dark() +
  labs(
    color = "Transit Steps",
    x = "Longitude",
    y = "Latitude"
  )

print(p_fig5B)

# =========================
# 5. Combine panels
# =========================
fig5_combined = ggarrange(
  p_fig5A,
  p_fig5B,
  labels = c("A", "B"),
  ncol = 2,
  nrow = 1,
  common.legend = TRUE,
  legend = "top"
)

print(fig5_combined)

# Optional save
# ggsave("Figure5_unknown_event_trajectories.png", fig5_combined, width = 12, height = 6, dpi = 300)