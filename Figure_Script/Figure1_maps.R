# ----------------------------------------
# Figure1_maps.R
# Purpose:
#   Generate the three component maps used in Figure 1:
#   (1) inset map of North America with Canada highlighted
#   (2) regional map showing historical caribou GPS locations by herd
#   (3) study area map of Contwoyto Lake on satellite imagery
#
# Notes:
#   - Final panel assembly was completed in PowerPoint.
#   - The study area box and connecting lines were added in PowerPoint.
#   - The left main panel shows historical GPS locations, not polygon ranges.
#
# Required input objects/files:
#   - nwt_lakes.RData
#   - Contwoyto_pure_water.shp
# ----------------------------------------

library(dplyr)
library(sf)
library(ggplot2)
library(ggspatial)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggmap)

# =========================
# 1. Load movement data
# =========================
load("data/nwt_lakes.RData")
nwt_lakes = st_transform(nwt_lakes, 4326)

# lake polygon for study area map
polys1_sf = st_read("data/Contwoyto_pure_water.shp", quiet = TRUE)
polys1_sf = st_transform(polys1_sf, 4326)

# =========================
# 2. Base world data
# =========================
world = ne_countries(scale = "medium", returnclass = "sf")
north_america = world %>%
  filter(continent == "North America")

CAN = world %>%
  filter(admin == "Canada")

# =========================
# 3. Figure 1A: inset map
# =========================
fig1_inset = ggplot(data = north_america) +
  geom_sf(fill = "white", color = "black", linewidth = 0.2) +
  geom_sf(data = CAN, fill = "darkgrey", color = "black", linewidth = 0.2) +
  coord_sf(crs = st_crs(32612)) +
  theme_bw() +
  labs(
    x = "Longitude",
    y = "Latitude"
  )

print(fig1_inset)

# optional save
# ggsave("Figure1_inset_map.png", fig1_inset, width = 4, height = 4, dpi = 300)

# =========================
# 4. Figure 1B: regional map
#    historical GPS locations by herd
# =========================
fig1_regional = ggplot() +
  geom_sf(data = CAN, fill = "darkgrey") +
  geom_sf(data = nwt_lakes, aes(color = Herd), size = 1) +
  coord_sf(crs = st_crs(32612)) +
  theme_bw() +
  labs(
    x = "Longitude",
    y = "Latitude",
    color = "Caribou Herd"
  ) +
  theme(
    legend.position = c(0.85, 0.25),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.key = element_rect(fill = "white"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )

print(fig1_regional)

# optional save
# ggsave("Figure1_regional_map.png", fig1_regional, width = 8, height = 6, dpi = 300)

# =========================
# 5. Figure 1C: study area satellite map
# =========================

# bounding center used in your original script
center_lon = mean(c(-111.61050, -108.79082))
center_lat = mean(c(66.00000, 64.96828))

# IMPORTANT:
# replace with your own key or set it securely elsewhere
# register_google(key = "YOUR_GOOGLE_MAPS_KEY")

sat_map = get_googlemap(
  center = c(lon = center_lon, lat = center_lat),
  zoom = 8,
  size = c(640, 640),
  scale = 2,
  maptype = "satellite"
)

fig1_study_area = ggmap(sat_map) +
  geom_sf(
    data = polys1_sf,
    inherit.aes = FALSE,
    fill = "white",
    color = "white",
    linewidth = 2
  ) +
  annotation_scale(
    location = "bl",
    width_hint = 1,
    bar_cols = c("white", "black"),
    text_col = "white",
    line_col = "white",
    pad_y = unit(0.5, "cm"),
    unit = "km",
    style = "bar"
  ) +
  annotation_north_arrow(
    location = "tr",
    which_north = "true",
    style = north_arrow_fancy_orienteering(
      fill = c("white", "black"),
      text_col = "white"
    ),
    height = unit(2, "cm"),
    width = unit(2, "cm")
  ) +
  coord_sf(crs = 4326) +
  theme_minimal() +
  labs(
    x = "Longitude",
    y = "Latitude"
  )

print(fig1_study_area)

# optional save
# ggsave("Figure1_study_area_map.png", fig1_study_area, width = 8, height = 8, dpi = 300)

# =========================
# 6. Notes for final assembly
# =========================
# Final Figure 1 was assembled in PowerPoint using:
#   - Figure1_inset_map.png
#   - Figure1_regional_map.png
#   - Figure1_study_area_map.png
#
# In PowerPoint, add:
#   - the study area box on the regional map
#   - the connecting lines between the regional map and study area map
#   - final panel alignment and label placement











