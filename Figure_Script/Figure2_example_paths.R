# ----------------------------------------
# Figure2_example_paths.R
# Purpose:
#   Plot one illustrative unknown transit event showing:
#   (1) direct path
#   (2) circumnavigate path
#   (3) reference path (3 steps before and after)
#
# Example event: bev042, 2010
# ----------------------------------------

library(dplyr)
library(lubridate)
library(sf)
library(terra)
library(raster)
library(gdistance)
library(fasterize)
library(ggplot2)

# =========================
# 1. Load data
# =========================
load("data/nwt_lakes.RData")
nwt_lakes = st_transform(nwt_lakes, 4326) %>% arrange(Time)

polys1_sf = st_read("data/Contwoyto_pure_water.shp", quiet = TRUE)
polys1_sf = st_transform(polys1_sf, 4326)

# for plotting
lake = polys1_sf

# =========================
# 2. Select one example unknown event
# =========================
target_year = 2010
target_id = "bev042"

caribou_target = nwt_lakes %>%
  filter(year(Time) == target_year, ID == target_id) %>%
  arrange(Time) %>%
  st_drop_geometry()

# manually define focal event window
event_start = as.POSIXct("2010-05-22 15:00:00", tz = "UTC")
event_end   = as.POSIXct("2010-05-24 03:00:00", tz = "UTC")

cross_line = caribou_target %>%
  filter(Time >= event_start, Time <= event_end) %>%
  arrange(Time)

# make sure at least two points exist
stopifnot(nrow(cross_line) >= 2)

# focal start and end points
pt_start = cross_line[1, ]
pt_end   = cross_line[nrow(cross_line), ]

# =========================
# 3. Direct path
# =========================
cross_line_vct = vect(
  cross_line,
  geom = c("Lon", "Lat"),
  crs = "EPSG:4326"
)

cross_line_sf = st_as_sf(as.lines(cross_line_vct))

# =========================
# 4. Reference path
#    three steps before and after the event
# =========================
start_idx = which(caribou_target$Time == event_start)
end_idx   = which(caribou_target$Time == event_end)

# if exact matching fails, use nearest
if (length(start_idx) == 0) {
  start_idx = which.min(abs(as.numeric(caribou_target$Time) - as.numeric(event_start)))
}
if (length(end_idx) == 0) {
  end_idx = which.min(abs(as.numeric(caribou_target$Time) - as.numeric(event_end)))
}

ref_start = max(1, start_idx - 3)
ref_end   = min(nrow(caribou_target), end_idx + 3)

reference_df = caribou_target %>%
  slice(ref_start:ref_end) %>%
  arrange(Time)

reference_line_vct = vect(
  reference_df,
  geom = c("Lon", "Lat"),
  crs = "EPSG:4326"
)

reference_line_sf = st_as_sf(as.lines(reference_line_vct))

# =========================
# 5. Circumnavigate path
#    least-cost path over land
# =========================
lake_proj = st_transform(polys1_sf, 32612)

pt1 = st_as_sf(
  data.frame(Lon = pt_start$Lon, Lat = pt_start$Lat),
  coords = c("Lon", "Lat"),
  crs = 4326
) %>% st_transform(32612)

pt2 = st_as_sf(
  data.frame(Lon = pt_end$Lon, Lat = pt_end$Lat),
  coords = c("Lon", "Lat"),
  crs = 4326
) %>% st_transform(32612)

# build buffered domain around lake + points
bbox_proj = st_as_sfc(
  st_bbox(st_union(lake_proj, st_union(pt1, pt2)))
) %>%
  st_buffer(dist = 60000)

coast = st_difference(bbox_proj, st_union(lake_proj))
coast_bbox = st_bbox(coast)

r = raster(
  xmn = coast_bbox["xmin"],
  xmx = coast_bbox["xmax"],
  ymn = coast_bbox["ymin"],
  ymx = coast_bbox["ymax"],
  res = 100,
  crs = st_crs(lake_proj)$wkt
)
r[] = 1

lake_union = st_as_sf(st_union(lake_proj))
lake_union$val = 1

lake_raster = fasterize(lake_union, r, field = "val")
land_mask = raster::mask(r, lake_raster, inverse = TRUE)
land_mask[land_mask < 0] = NA

tr = transition(
  land_mask,
  transitionFunction = function(x) 1,
  directions = 16,
  symm = TRUE
)
tr = geoCorrection(tr, type = "c", scl = FALSE)

circumvent_line = shortestPath(
  tr,
  origin = st_coordinates(pt1),
  goal = st_coordinates(pt2),
  output = "SpatialLines"
)

circumvent_line_sf = st_as_sf(circumvent_line) %>%
  st_set_crs(32612) %>%
  st_transform(4326)

# =========================
# 6. GPS points for plotting
# =========================
caribou_points_sf = st_as_sf(
  reference_df,
  coords = c("Lon", "Lat"),
  crs = 4326,
  remove = FALSE
)

# =========================
# 7. Plot extent
# =========================
bounds = st_bbox(circumvent_line_sf)

# =========================
# 8. Plot
# =========================
p_fig2 = ggplot() +
  geom_sf(data = lake, aes(fill = "Lake Area"), color = NA, alpha = 0.5) +
  geom_sf(data = reference_line_sf, aes(color = "Reference Path"), linewidth = 1.2) +
  geom_sf(data = cross_line_sf, aes(color = "Direct Path"), linewidth = 1.2) +
  geom_sf(data = circumvent_line_sf, aes(color = "Circumnavigate Path"), linewidth = 1.2) +
  geom_sf(data = caribou_points_sf, aes(color = "Caribou GPS Locations"), size = 2) +
  scale_fill_manual(
    values = c("Lake Area" = "black"),
    name = "Feature Type"
  ) +
  scale_color_manual(
    values = c(
      "Direct Path" = "yellow",
      "Circumnavigate Path" = "red",
      "Reference Path" = "blue",
      "Caribou GPS Locations" = "white"
    ),
    name = "Path Type"
  ) +
  coord_sf(
    xlim = c(bounds["xmin"] - 0.8, bounds["xmax"] + 1),
    ylim = c(bounds["ymin"] - 0.5, bounds["ymax"] + 0.5),
    expand = FALSE
  ) +
  theme_dark() +
  labs(
    x = "Longitude",
    y = "Latitude"
  )

print(p_fig2)

# Optional save
# ggsave("Figure2_example_paths.png", p_fig2, width = 8, height = 6, dpi = 300)
