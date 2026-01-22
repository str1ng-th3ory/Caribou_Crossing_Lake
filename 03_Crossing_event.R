library(lubridate)
library(dplyr)
library(terra)
library(sf)
library(R.matlab) 
library(sf)
library(gdistance)
library(fasterize)
library(usethis)
library(devtools)
library(recolorize)
library(mapview)
mapviewOptions()
mapviewOptions(method="ngb")
library(gdistance)
library(raster)
library(geosphere)
require(devtools)
library(googlesheets4)

# load caribou movement data
load("/Users/Qianr/Documents/GitHub/dessertation_phd/nwt_lakes.RData")
nwt_lakes  <- sf::st_transform(nwt_lakes, 4326)

#older version
#require(rgdal)
#polys1 <- readOGR(dsn = ".", layer = "Contwoyto_pure_water")

# lake polygon
polys1_sf = sf::st_read(  "/Users/Qianr/Documents/GitHub/dessertation_phd/Contwoyto_pure_water.shp", quiet = TRUE)
polys1_sf <- sf::st_transform(polys1_sf, 4326)

#new version
polys1 = as(polys1_sf, "Spatial")

# Find out Known Crossing Event
nwt_lakes_new <- st_transform(nwt_lakes, crs = '+proj=longlat +datum=WGS84 +no_defs')
nwt_lakes=nwt_lakes_new%>%arrange(Time)

# load information about crossing points
# option1: link google sheet
url <- 'link' # Replace the access link to the spreadsheets which contains crossing information: Year, ID, Albedo,Crossing_date,Lon,Lat,Time)
whole_list=gsheet2tbl(url) 

# option2: use document saved before
on_lake=final_data_frame #(see Point_within_lake. R) or load "on_lake.csv"

initial_list=as.list(on_lake)

pair_list <- list()

current_year <- NA
percentile_rank <- NULL
albedo_data <- NULL
arr_climatology_brick <- NULL

for (i in 1:length(initial_list$Year)) {
 
  year <- initial_list$Year[i]
  id <- initial_list$ID[i]
  
  #load albedo data based on Year
  if (is.na(current_year) || year != current_year) {
    
    current_year <- year
    data_filename <- paste0("output_", year, ".mat")
    
    climatology_data <- readMat(data_filename)
    arr_climatology <- array(unlist(climatology_data), dim = c(426, 724, 365))
    
    arr_climatology_brick <- brick(
      arr_climatology,
      xmn = -111.7125, xmx = -108.6958,
      ymn = 64.375, ymx = 66.15,
      crs = "+proj=longlat +datum=WGS84",
      transpose = FALSE
    )
  
    arr_climatology_brick_crop <- crop(arr_climatology_brick, extent(polys1))
    albedo_data <- mask(arr_climatology_brick_crop, polys1)
    
    percentile_rank <- calc(albedo_data[[92:280]], fun = function(x) {
      rank(x, na.last = "keep") / sum(!is.na(x))
    })
  }
  
  # obtain the current year's trajectory
  data <- subset(nwt_lakes, lubridate::year(nwt_lakes$Time) == year)
  id_data <- data[data$ID %in% id,]
  whole_line_year <- data.frame(id_data) %>% arrange(Time)

  whole_line_year$yday=yday(whole_line_year$Time)
  
  whole_line<- subset(whole_line_year, yday >= 92 & yday <= 280) 
  
  if (nrow(whole_line) < 2) next
  
  # find out index based on target time
  target_time <- as.POSIXct(initial_list$Time[i], tz = "UTC", format = "%Y-%m-%d %H:%M:%S")
  
  # match
  index <- which(whole_line$Time == target_time)
  # option:  index <- which.min(abs(as.numeric(whole_line$Time) - as.numeric(target_time)))
  
  before_index=index
  after_index=index+1
  if (after_index > nrow(whole_line)) next
  
  nearest_before <- whole_line[before_index, ]
  nearest_after <- whole_line[after_index, ]
  
  crossing_date=nearest_before$yday
  # migration season definition:
  # within yday 92–280:
  #   month <= 6  -> spring migration
  #   month >= 7  -> fall migration
  season_value <- ifelse(format(nearest_before$Time, "%m") > "06", "fall", "spring")
  
  #lake_width: The length of the intersection between the lake polygon and the line connecting consecutive sampling points.  
  point1 <- whole_line[before_index, c("Lon", "Lat")]
  point2 <- whole_line[after_index, c("Lon", "Lat")]
  line <- Line(rbind(point1, point2))
  lines <- Lines(list(line), ID="1")
  whole_line_vct_lines <- SpatialLines(list(lines))
  proj4string(whole_line_vct_lines) <- CRS('+proj=longlat +datum=WGS84 +no_defs')
  proj4string(polys1) <- CRS('+proj=longlat +datum=WGS84 +no_defs')
  
  sfLines <- st_as_sf(whole_line_vct_lines)
  sfPolygons <- st_as_sf(polys1)
  intersectionPoints <- st_intersection(sfLines,sfPolygons)
  lake_width=st_length(intersectionPoints)
  
  # straight distance/speed
  nearest_before_location <- c(nearest_before$Lon[[1]], nearest_before$Lat[[1]])
  nearest_after_location  <- c(nearest_after$Lon[[1]],  nearest_after$Lat[[1]])
  
  #straight distance/speed
  straight_distance=geosphere::distm(nearest_before_location,nearest_after_location, fun = distHaversine)[1,1]
  time_difference <- nearest_after$Time - nearest_before$Time
  time_difference_seconds <- as.numeric(time_difference, units = "secs")
  crossing_duration=as.numeric(time_difference, units = "days")
  
  straight_speed= straight_distance/time_difference_seconds
  
  # circumnavigate distance/speed ((projected to 32612)
  lake_proj <- sf::st_transform(polys1_sf, 32612)
  # make pt1 / pt2 (projected)
  pt1 <- sf::st_as_sf(
    data.frame(Lon = nearest_before$Lon[[1]], Lat = nearest_before$Lat[[1]]),
    coords = c("Lon", "Lat"),
    crs = 4326
  ) %>% sf::st_transform(32612)
  
  pt2 <- sf::st_as_sf(
    data.frame(Lon = nearest_after$Lon[[1]], Lat = nearest_after$Lat[[1]]),
    coords = c("Lon", "Lat"),
    crs = 4326
  ) %>% sf::st_transform(32612)
  
  # build coast (extent) based on lake bbox 
  bbox_proj <- sf::st_as_sfc(sf::st_bbox(lake_proj)) %>% sf::st_buffer(dist = 55000)
  coast <- sf::st_difference(bbox_proj, sf::st_union(lake_proj))
  
  # raster mask: land=1, lake=NA
  # optional: res can be different from 100m, 250m to 500m (same as albedo pixel)
  r <- raster::raster(raster::extent(coast), res = 100, crs = sf::st_crs(lake_proj)$wkt)
  r[] <- 1
  
  # rasterize lake onto r
  lake_raster <- fasterize::fasterize(sf::st_sf(geometry = sf::st_union(lake_proj)), r, field = NULL)
  # land_mask: lake cells become NA (blocked), land stays 1
  land_mask <- raster::mask(r, lake_raster, inverse = TRUE)
  land_mask[land_mask < 0] <- NA
  
  # transition + geoCorrection (projected CRS)
  tr <- gdistance::transition(land_mask,
                              transitionFunction = function(x) 1,
                              directions = 16,
                              symm = TRUE)
  tr <- gdistance::geoCorrection(tr, type = "c", scl = FALSE)
  
  #option: 
  # circumvent_path <- shortestPath(tr,
  #                                 origin = st_coordinates(pt1),
  #                                 goal = st_coordinates(pt2),
  #                                 output = "SpatialLines")
  
  # compute circumvent distance (meters)
  cost <- gdistance::costDistance(tr,
                                  sf::st_coordinates(pt1),
                                  sf::st_coordinates(pt2))
  circumvent_distance <- as.numeric(cost[1,1]) 
  circumvent_speed=circumvent_distance/time_difference_seconds

  # extract the albedo along the path
  rank_number <- crossing_date - 91 
  albedo_drop_percentile=percentile_rank[[(rank_number)]]
  test2=raster(albedo_drop_percentile)
  values(test2)=values(albedo_drop_percentile)
  albedo_drop_percentile_rast= terra:: rast(test2)
  
  line_before_after = whole_line %>% filter(between(whole_line$Time, 
                                                    nearest_before$Time, 
                                                    nearest_after$Time))
  
  line_before_after_vct=terra::vect(
    line_before_after, geom=c("Lon", "Lat"), 
    crs='+proj=longlat +datum=WGS84 +no_defs +type=crs',
    keepgeom=TRUE) %>% as.lines() 
  
  intersected_pixels=terra::extract(albedo_drop_percentile_rast, line_before_after_vct)
  
  #split dataframe results into multiple parts based on continious NA rows
  na_inds <- which(is.na(intersected_pixels$layer))
  start_inds <- c(1, na_inds[which(diff(na_inds) != 1)] + 1)
  end_inds <- c(na_inds[which(diff(na_inds) != 1)], nrow(intersected_pixels))
  split_dfs <- lapply(seq_along(start_inds), function(i) intersected_pixels[start_inds[i]:end_inds[i], ])
  
  # albedo percentile rank (APR) along the crossing path
  alebdo_linear=intersected_pixels$layer %>% mean(,na.rm=TRUE)

  # albedo percentile rank (APR) at the nearest pixels
  near_albedo=split_dfs[length(split_dfs)]%>%data.frame()%>%dplyr::select('layer') #option:split_dfs[[1]] %>% data.frame()
  albedo_nearest_pixel=near_albedo[(which(is.na(near_albedo))[1]-1),]
 
  # albedo percentile rank (APR) of the entire lake
  x.stats <- data.frame(x.mean=cellStats(albedo_data[[92:280]], "mean",na.rm=TRUE))
  albedo_whole_lake=(rank(x.stats$x.mean)/sum(!is.na(x.stats)))[(rank_number)]
  
  #add type for the crossing point 
  type=0
  
  #收集结果
  pair_list[[length(pair_list) + 1]] <- list(
    Year=year,
    ID= id,
    crossing_time=nearest_before$Time,
    crossing_duration=crossing_duration,
    crossing_date=nearest_before$yday,
    nearest_before_Lon=nearest_before$Lon[[1]],
    nearest_before_Lat=nearest_before$Lat[[1]],
    nearest_after_Lon=nearest_after$Lon[[1]],
    nearest_after_Lat=nearest_after$Lat[[1]],
    lake_width=lake_width,
    circumvent_distance=circumvent_distance,
    straight_distance=straight_distance,
    circumvent_speed=circumvent_speed,
    straight_speed= straight_speed,
    alebdo_linear=alebdo_linear,
    albedo_nearest_pixel=albedo_nearest_pixel,
    albedo_whole_lake=albedo_whole_lake,
    type=0,
    season=season_value)
  
  # add reference point for 
  reference_list=list()
  
  for (k in (c(
    index-1,
    index-2,
    index-3,
    index+1,
    index+2,
    index+3
  ))) 
    {
    if (k > 0 && k <= nrow(whole_line)) {
      nearest_location <- whole_line[k, ]
      nearest_location_next <- whole_line[min(k+1, nrow(whole_line)), ]  
      
      nearest_location_location <- c(nearest_location$Lon, nearest_location$Lat)
      displacement <- geosphere::distm(nearest_location_location, c(nearest_location_next$Lon, nearest_location_next$Lat), fun = distHaversine)[1,1]
      time_difference <- nearest_location_next$Time - nearest_location$Time
      time_difference_seconds <- as.numeric(time_difference, units = "secs")
      crossing_duration=as.numeric(time_difference, units = "days")
      average_speed <- ifelse(time_difference > 0, displacement / time_difference_seconds, NA)  
      season_value <- ifelse(format(nearest_location$Time, "%m") > "06", "fall", "spring")
      
      pair_list[[length(pair_list) + 1]] <- list(
        Year = Year,
        ID = id,
        crossing_time = nearest_location$Time,
        crossing_duration = crossing_duration,
        crossing_date = nearest_location$yday,
        nearest_before_Lon = nearest_location$Lon[[1]],
        nearest_before_Lat = nearest_location$Lat[[1]],
        nearest_after_Lon = nearest_location_next$Lon[[1]],
        nearest_after_Lat = nearest_location_next$Lat[[1]],
        lake_width = NA,
        circumvent_distance = NA,
        straight_distance = displacement,
        circumvent_speed = NA,
        straight_speed = average_speed,
        alebdo_linear = NA,
        albedo_nearest_pixel = NA,
        albedo_whole_lake = NA,
        season=season_value
        )
    }
  }
  pair_list <- append(pair_list, reference_list)
}

pair_list_df <- lapply(pair_list, function(x) as.data.frame(t(unlist(x))))
result_df <- do.call(rbind, pair_list_df)
result_df$crossing_time <- as.POSIXct(as.numeric(result_df$crossing_time), origin = "1970-01-01", tz = "UTC")

# assign "type" to distinguish the type of points in the trasit event(0 means crossing point)
result_df$type=c(0,-3,-2,-1, 1,2,3)

#save
write.csv(
  result_df,
  paste0("crossing_event", ".csv"),
  row.names = FALSE
)

# option
result_df %>%
  write_sheet(
    ss = gs4_get(
      "link" # Replace the access link to the spreadsheets
    ),
    sheet = "crossing_event")
