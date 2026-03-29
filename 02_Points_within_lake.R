# ----------------------------------------
# 02_extract_on_lake_points_and_albedo.R
# Purpose:
#   (1) Identify caribou GPS points that fall within the Contwoyto Lake polygon
#   during the analysis period (DOY 92–280), and
#   (2) extract corresponding pixel-level albedo percentile information.
#
# Notes:
#   - Raw caribou GPS data are not distributed with this repository because they
#     are subject to third-party data-sharing restrictions from GNWT-ECC.
#   - This script also requires pre-processed annual .mat files containing daily
#     albedo climatology outputs.
#
# Expected input objects/files:
#   - nwt_lakes.RData
#   - Contwoyto_pure_water.shp
#   - output_<year>.mat for each year analyzed
# ----------------------------------------

library(lubridate)
library(dplyr)
library(terra)
library(sf)
library(R.matlab)
library(raster)

##### Step1: Find out Points Within Lake

# figure out "on lake" points (fall within the lake polygon)

# load caribou movement data
## Data files are not included due to privacy/size constraints. Paths are for demonstration.
load("data/nwt_lakes.RData")

#older version
#require(rgdal)
#polys1 = readOGR(dsn = ".", layer = "Contwoyto_pure_water")

# lake polygon
polys1_sf = sf::st_read(  "data/Contwoyto_pure_water.shp", quiet = TRUE)
polys1_sf = sf::st_transform(polys1_sf, 4326)
polys1 = as(polys1_sf, "Spatial")
# option: lake = terra::vect(polys1_sf)  

#movement data
nwt_lakes = st_transform(
  nwt_lakes,
  crs = '+proj=longlat +datum=WGS84 +no_defs'
)
#nwt_lakes  = sf::st_transform(nwt_lakes, 4326)

years = 2000:2021

df = nwt_lakes %>%
  mutate(yday = yday(Time)) %>%
  filter(year(Time) %in% years, yday >= 92, yday <= 280) %>%
  arrange(Time)

df_no_geom = sf::st_drop_geometry(df)
pts = terra::vect(df_no_geom, geom = c("Lon", "Lat"), crs = "EPSG:4326", keepgeom = TRUE)

# find out points fall within the lake
lake = terra::vect(polys1_sf)
inside = relate(pts, lake, "within")
onlake_pts = pts[inside, ]

# transfer to dataframe
onlake_df = as.data.frame(onlake_pts)

# year-id list
onlake_year_id = onlake_df %>%
  mutate(Year = year(Time)) %>%
  distinct(Year, ID) %>%
  arrange(Year, ID)

# save
write.csv(
  onlake_year_id,
  paste0("on_lake_candidates", ".csv"),
  row.names = FALSE
)

# Optional: export candidate list to Google Sheets
# library(googlesheets4)
# onlake_year_id %>%
#   write_sheet(
#     ss = gs4_get("link"),
#     sheet = "on_lake_candidate"
#   )
#head(onlake_year_id)


##### Step 2: Obtain Albedo Data for Each Points on Lake  
# url = 'link'
# on_lake = gsheet2tbl(url)
on_lake=onlake_year_id

initial_list=as.list(on_lake)

# Creating a new list to store the results
result_list = list()

current_year = NA
percentile_rank = NULL
arr_climatology_brick = NULL
albedo_data = NULL

for (i in 1:length(initial_list$Year)) {
  year = initial_list$Year[i]
  id = initial_list$ID[i]
  
  if (is.na(current_year) || year != current_year) {
    current_year = year
    data_filename = paste0("output_", year, ".mat") #change year here
    
    climatology_data = readMat(data_filename)
    arr_climatology = array(unlist(climatology_data), dim = c(426, 724, 365))

    arr_climatology_brick = brick(
      arr_climatology, 
      xmn = -111.7125, xmx = -108.6958, 
      ymn = 64.375, ymx = 66.15, 
      crs = '+proj=longlat +datum=WGS84', 
      transpose = FALSE)
    

    arr_climatology_brick_crop = crop(arr_climatology_brick, extent(polys1))
    albedo_data = mask(arr_climatology_brick_crop, polys1)

    percentile_rank = calc(albedo_data[[92:280]], fun = function(x) {
      percentile = rank(x, na.last = "keep") / sum(!is.na(x))
      return (percentile)
      })
  }
  
  # Loading data for the specific year
  # Filtering data
  data = subset(nwt_lakes, year(nwt_lakes$Time) == year)
  id_data = data[data$ID %in% id,]
  whole_line = data.frame(id_data) %>% arrange(Time)
  
  if ("sfc" %in% class(whole_line) || any(grepl("geometry", names(whole_line), ignore.case = TRUE))) {
    whole_line = sf::st_drop_geometry(whole_line)
  }
  
  # Calculating intersection points
  cross_line_vct = terra::vect(whole_line, geom = c("Lon", "Lat"), 
                                crs = '+proj=longlat +datum=WGS84 +no_defs +type=crs', 
                                keepgeom = TRUE)
  
  polygon_lake = terra::vect(polys1)
  x = terra::intersect(polygon_lake, cross_line_vct)
  
  #check if there's intersection 
  if (is.empty(x)) {
    print(paste("No intersection for Year:", year, "ID:", id))
    next
  }
  
  if (nrow(x) == 0) {
    print(paste("No intersection for Year:", year, "ID:", id))
    next
  }
  
    # if have at least 1 intersection
    for (j in 1:nrow(x)) {
      
      # the jth intersected point
      x_j = x[j, ]
      day_on_lake =lubridate::yday(x_j$Time)
      
      # check if day_on_lake is within DOY 92 to 280 observation period
      if (day_on_lake < 92 || day_on_lake > 280) {
        print(paste("No spring migration on ice for Year:", year, "ID:", id, "Day:", day_on_lake))
        next
      }
      
      rank_number = day_on_lake - 91
      albedo_drop_percentile = percentile_rank[[rank_number]]
      
      test2 = raster(albedo_drop_percentile)
      values(test2) = values(albedo_drop_percentile)
      albedo_drop_percentile_rast = terra::rast(test2)
      
      intersected_pixels = terra::extract(albedo_drop_percentile_rast, x_j)
      albedo = intersected_pixels$layer
      
      # Adding the results to the new list
      result_list = append(result_list, list(list(
        Year = year, 
        ID = id, 
        Albedo = albedo,
        DayOfYear = day_on_lake, 
        Lon = x_j$Lon, 
        Lat = x_j$Lat, 
        Time = x_j$Time
      )))
    }
  if (i %% 50 == 0) message("Progress: ", i, "/", length(initial_list$Year))
  #save(result_list, file = "result_list.RData")
  #print(result_list)
  }

# organize the results

result_list_df = lapply(result_list, function(x) as.data.frame(t(unlist(x))))
final_data_frame = do.call(rbind, result_list_df)
final_data_frame$Time = as.POSIXct(as.numeric(final_data_frame$Time), origin = "1970-01-01", tz = "UTC")

# transfer the list to a dataframe
data_frames = lapply(result_list, function(x) {
  data.frame(
    Year = x$Year,
    ID = x$ID,
    Albedo=x$Albedo,
    Crossing_date = x$DayOfYear,
    Lon = x$Lon,
    Lat = x$Lat,
    Time = x$Time
  )
})

final_data_frame = bind_rows(data_frames)

#save
write.csv(
  final_data_frame,
  "on_lake.csv",
  row.names = FALSE
)

# Optional: export results to Google Sheets
# library(googlesheets4)
# final_data_frame %>%
#   write_sheet(
#     ss = gs4_get("link"),
#     sheet = "on_lake_result"
#   )


