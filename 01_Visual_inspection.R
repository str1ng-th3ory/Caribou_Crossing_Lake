
#check "on lake" points
library(lubridate)
library(dplyr)
library(terra)
library(sf)
library(leaflet)

# load caribou movement data
## Data files are not included due to privacy/size constraints. Paths are for demonstration.
load("data/nwt_lakes.RData")

# lake polygon
polys1 <- vect("Contwoyto_pure_water.shp")

crs(polys1) <- "EPSG:4326"

nwt_lakes_new <- st_transform(
  nwt_lakes,
  crs = '+proj=longlat +datum=WGS84 +no_defs'
)
nwt_lakes <- nwt_lakes_new

# pick a year, for example 2005
target_year <- 2005
caribou_year <- subset(nwt_lakes, year(nwt_lakes$Time) == target_year)
unique_ids <- unique(caribou_year$ID)

process_id <- function(id, data) {
  id_data <- data[data$ID %in% id, ]
  whole_line <- data.frame(id_data) %>% arrange(Time)
  if (nrow(whole_line) < 2) return(NULL)
  
  v <- terra::vect(
    whole_line,
    geom = c("Lon", "Lat"),
    crs = "+proj=longlat +datum=WGS84 +no_defs"
  )
  terra::as.lines(v)
}

results <- list()
for (id in unique_ids) {
  ln <- process_id(id, caribou_year)   
  if (!is.null(ln)) results[[id]] <- ln
}

results_sf <- lapply(results, sf::st_as_sf)

# visual inspection for the crossing
leaflet_map <- leaflet() %>%
  addProviderTiles(providers$CartoDB.DarkMatterNoLabels) %>%
  addPolygons(data = polys1, fillColor = "#FFFAFA", color = "#FFFFFF",
              smoothFactor = 0.5, group = "Lake")

for (id in names(results_sf)) {
  leaflet_map <- leaflet_map %>%
    addPolylines(data = results_sf[[id]], color = "red",
                 weight = 2, opacity = 0.8, group = id)
}

leaflet_map <- leaflet_map %>%
  addLayersControl(
    overlayGroups = c("Lake", names(results_sf)),
    options = layersControlOptions(collapsed = FALSE)
  )

leaflet_map

# ===============================
# MANUAL RECORD (after visual check)
# ===============================

target_year <- 2005

# for example
candidate_ids <- c(
  "bat068",
  "bat089"
)

on_lake_candidates <- data.frame(
  Year = target_year,
  ID   = candidate_ids,
  source = "visual_check",     
  stringsAsFactors = FALSE
)

#save
write.csv(
  on_lake_candidates,
  paste0("on_lake_candidates", target_year, ".csv"),
  row.names = FALSE
)

#or google sheet
on_lake_candidates %>%
  write_sheet(
    ss = gs4_get("link"),
    sheet = paste0("on_lake_candidates", target_year)
  )










