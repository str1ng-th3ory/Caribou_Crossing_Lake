# ----------------------------------------
# 01_visual_inspection.R
# Purpose:
#   Visualize annual caribou tracks relative to Contwoyto Lake
#   and manually record candidate on-lake crossing events.
#
# Note:
#   The caribou GPS dataset is not distributed with this repository because it is
#   subject to third-party data-sharing restrictions from GNWT-ECC.
#   Qualified researchers may request access directly from GNWT-ECC.
#
# Expected input object:
#   nwt_lakes: an sf object containing at least the following columns:
#   ID, Time, Lon, Lat
# ----------------------------------------
library(lubridate)
library(dplyr)
library(terra)
library(sf)
library(leaflet)

# Restricted input data are not included in this repository.
# Replace the path below with the local path to the GNWT-ECC-approved data file.
data_dir = "data"
load(file.path(data_dir, "nwt_lakes.RData"))

# Lake polygon
polys1 <- vect(file.path(data_dir, "Contwoyto_pure_water.shp"))
crs(polys1) = "EPSG:4326"

nwt_lakes = st_transform(nwt_lakes, crs = 4326)

# Pick a year, for example 2005
target_year = 2005
caribou_year = subset(nwt_lakes, year(nwt_lakes$Time) == target_year)
unique_ids = unique(caribou_year$ID)

process_id = function(id, data) {
  id_data = data[data$ID %in% id, ]
  whole_line = data.frame(id_data) %>% arrange(Time)
  if (nrow(whole_line) < 2) return(NULL)
  
  v = terra::vect(
    whole_line,
    geom = c("Lon", "Lat"),
    crs = "EPSG:4326"
  )
  terra::as.lines(v)
}

results = list()
for (id in unique_ids) {
  ln = process_id(id, caribou_year)   
  if (!is.null(ln)) results[[id]] = ln
}

results_sf = lapply(results, sf::st_as_sf)

# Visual inspection for the crossing
leaflet_map = leaflet() %>%
  addProviderTiles(providers$CartoDB.DarkMatterNoLabels) %>%
  addPolygons(data = polys1, fillColor = "#FFFAFA", color = "#FFFFFF",
              smoothFactor = 0.5, group = "Lake")

for (id in names(results_sf)) {
  leaflet_map = leaflet_map %>%
    addPolylines(data = results_sf[[id]], color = "red",
                 weight = 2, opacity = 0.8, group = id)
}

leaflet_map = leaflet_map %>%
  addLayersControl(
    overlayGroups = c("Lake", names(results_sf)),
    options = layersControlOptions(collapsed = FALSE)
  )

leaflet_map

# ===============================
# MANUAL RECORD (after visual check)
# ===============================

# for example
candidate_ids = c(
  "bat068",
  "bat089"
)

on_lake_candidates = data.frame(
  Year = target_year,
  ID   = candidate_ids,
  source = "visual_check",     
  stringsAsFactors = FALSE
)

#save
output_dir <- "outputs"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

write.csv(
  on_lake_candidates,
  file.path(output_dir, paste0("on_lake_candidates_", target_year, ".csv")),
  row.names = FALSE
)

# Optional: export to Google Sheets
# Requires the googlesheets4 package and authenticated access.
# library(googlesheets4)
# on_lake_candidates %>%
#   write_sheet(
#     ss = gs4_get("link"),
#     sheet = paste0("on_lake_candidates", target_year)
#   )

## optional
# Check the track and lake 
library(dplyr)
library(lubridate)
library(terra)
library(sf)

load("data/nwt_lakes.RData")
nwt_lakes = st_transform(nwt_lakes, 4326) %>% arrange(Time)

polys1_sf = st_read("data/Contwoyto_pure_water.shp", quiet = TRUE)
polys1_sf = st_transform(polys1_sf, 4326)
polygon_lake = terra::vect(polys1_sf)

data_i = nwt_lakes %>%
  filter(Year == 2010, ID == "bat080") %>% #need input the year and id number
  arrange(Time) %>%
  mutate(yday = yday(Time)) %>%
  filter(yday >= 182, yday <= 280) #need input day 

pts_vct = terra::vect(
  data.frame(data_i),
  geom = c("Lon", "Lat"),
  crs = "EPSG:4326",
  keepgeom = TRUE
)

inside_pts = terra::relate(pts_vct, polygon_lake, "within")

data_i$inside_lake = inside_pts
data_i %>% filter(inside_lake)
