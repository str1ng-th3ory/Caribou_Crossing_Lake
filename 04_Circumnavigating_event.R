# ----------------------------------------
# 04_circumnavigating_event.R
# Purpose:
#   (1) Identify candidate circumnavigation events using spatial criteria,
#   (2) manually inspect candidate trajectories, and
#   (3) compute event-level movement and albedo/APR metrics for confirmed
#       circumnavigation events.
#
# Notes:
#   - Raw caribou GPS data are not distributed with this repository because they
#     are subject to third-party data-sharing restrictions from GNWT-ECC.
#   - This script also requires preprocessed annual .mat files containing daily
#     albedo climatology outputs.
#
# Expected input objects/files:
#   - nwt_lakes.RData
#   - Contwoyto_pure_water.shp
#   - output_<year>.mat for each year analyzed
#   - manually checked circumnavigation table with before_index / after_index
# ----------------------------------------

library(lubridate)
library(dplyr)
library(terra)
library(sf)
library(R.matlab)
library(gdistance)
library(fasterize)
library(mapview)
mapviewOptions(method = "ngb")
library(raster)
library(geosphere)

# load caribou movement data
load("data/nwt_lakes.RData")
nwt_lakes = st_transform(nwt_lakes, 4326) %>% arrange(Time)
# older version
# require(rgdal)
# polys1 = readOGR(dsn = ".", layer = "Contwoyto_pure_water")

# lake polygon
polys1_sf = sf::st_read("data/Contwoyto_pure_water.shp", quiet = TRUE)
polys1_sf = sf::st_transform(polys1_sf, 4326)
# new version
polys1 = as(polys1_sf, "Spatial")
polygon_lake = terra::vect(polys1_sf) 

# =========================
######Step 1: Find candidates
# an automated spatial pre-screening based on geometric criteria (lake buffer + axis intersection)
# =========================

# buffer (20 km guiding framework)
lake_utm = st_transform(polys1_sf, 32612)
lake_buffer_utm = st_buffer(lake_utm, 20000)
lake_buffer = st_transform(lake_buffer_utm, 4326)
polygon_lake_buffer = terra::vect(lake_buffer)

#calculate the median line along the lake’s long axis
lake_hull = st_convex_hull(polys1_sf)

extendLongestAxisPCA = function(polygon, extension_factor = 2.5) {
  coords = st_coordinates(polygon)
  pca = prcomp(coords[, 1:2])
  center = colMeans(coords[, 1:2])
  pc1 = pca$rotation[, 1] * sqrt(pca$sdev[1]) * extension_factor
  axis_start = center - pc1
  axis_end = center + pc1
  st_sfc(st_linestring(rbind(axis_start, axis_end)), crs = st_crs(polygon))
}

extended_axis = extendLongestAxisPCA(lake_hull, extension_factor = 2.5)
extended_axis_terra = terra::vect(extended_axis)

mapview(extended_axis)+mapview(polys1_sf)

#find out circumnavigating case candidates
#find out all possible year-id combos
year_id_combinations = expand.grid(
  year = unique(nwt_lakes$Year),
  id   = unique(nwt_lakes$ID),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

actual_combinations = year_id_combinations[
  with(year_id_combinations, paste(year, id) %in% paste(nwt_lakes$Year, nwt_lakes$ID)),
]

#2 season
season_windows = list(
  spring = c(92, 181), #spring migration season (DOY 92 - DOY 181)
  fall   = c(182, 280) #fall migration season (DOY 182 - DOY 280)
)

results = list()

for (i in seq_len(nrow(actual_combinations))) {
  
  year_i = actual_combinations$year[i]
  id_i   = as.character(actual_combinations$id[i])
  
  data_i = nwt_lakes %>%
    filter(Year == year_i, ID == id_i) %>%
    arrange(Time)
  
  if (nrow(data_i) == 0) next
  
  whole_line_year = data.frame(data_i) %>%
    arrange(Time) %>%
    mutate(yday = lubridate::yday(Time))
  
  for (ss in names(season_windows)) {
    
    y1 = season_windows[[ss]][1]
    y2 = season_windows[[ss]][2]
    
    whole_line = subset(whole_line_year, yday >= y1 & yday <= y2)
    
    if (nrow(whole_line) < 2) {next}
    
    # individual trajectory（terra）
    cross_line_vct = terra::vect(
      whole_line,
      geom = c("Lon", "Lat"),
      crs = "+proj=longlat +datum=WGS84 +no_defs",
      keepgeom = FALSE
    )
    whole_line_vct_line = terra::as.lines(cross_line_vct)
    
    # must meet three requirements
    # 1) intersect with lake_buffer（individual enter into the buffer）
    intersection_buffer = terra::intersect(polygon_lake_buffer, whole_line_vct_line)
    if (terra::is.empty(intersection_buffer) || nrow(intersection_buffer) == 0) next
    
    # 2) intersect with extended median axis (lake long axis)（"Candidate features" that cross from one side of the axis to the other.）
    intersection_axis = terra::intersect(extended_axis_terra, whole_line_vct_line)
    if (terra::is.empty(intersection_axis) || nrow(intersection_axis) == 0) next
    
    # 3) not intersect with lake polygon（not “crossing”，more likely "circumnavigating"）
    disjoint_lake = terra::intersect(polygon_lake, whole_line_vct_line)
    if (!terra::is.empty(disjoint_lake) && nrow(disjoint_lake) > 0) next
    
    # candidate
    results[[length(results) + 1]] = data.frame(
      Year   = year_i,
      ID     = id_i,
      season = ss,
      stringsAsFactors = FALSE
    )
  }
}

#option
# mapview(lake_buffer) +
#   mapview(extended_axis, color="pink") +
#   mapview(st_as_sf(whole_line_vct_line))


# save candidate list
results_df = if (length(results) == 0) {
  data.frame(Year = integer(), ID = character(), season = character())
} else {
  do.call(rbind, results)
}

# Optional: export candidate list to Google Sheets
# library(googlesheets4)
# results_df %>%
#   write_sheet(
#     ss = gs4_get("link"),
#     sheet = "Circumnavigate_initial"
#   )

results_df_circumnavigate = results_df

# =========================
######Step 2: Manual Inspection - MUST DO!!
# using trajectory visualization to confirm true circumnavigation and to define start/end points
# =========================

results_df_list = list()
for (i in 1:nrow(results_df)) {
  year_i   = results_df$Year[i]
  id_i     = as.character(results_df$ID[i])
  season_i = as.character(results_df$season[i]) 

  y1 = season_windows[[season_i]][1]
  y2 = season_windows[[season_i]][2]
  
  data = nwt_lakes %>%
    filter(lubridate::year(Time) == year_i, ID == id_i) %>%
    arrange(Time) %>%
    mutate(yday = yday(Time)) %>%
    filter(yday >= y1, yday <= y2)
  
  if (nrow(data) == 0) {
    results_df_list[[i]] = list(points = NULL, line = NULL)
    next
  }
  
  sf_data = st_as_sf(data, coords = c("Lon","Lat"), crs = 4326, remove = FALSE) %>%
    mutate(idx = dplyr::row_number())
  
  sf_data$popup = paste0(
    "idx: ", sf_data$idx,
    "<br>Time: ", format(sf_data$Time, "%Y-%m-%d %H:%M"),
    "<br>yday: ", sf_data$yday
  )
  
  line = NULL
  if (nrow(sf_data) > 1) {
    line = sf_data %>%
      summarize(geometry = st_combine(geometry)) %>%
      st_cast("LINESTRING")
  }
  
  results_df_list[[i]] = list(points = sf_data, line = line)
}

# Here, you need to change "i" and manually check each row in the candidate form (result_df)
i = 1 #need to change it each round 

geometry_list = results_df_list[[i]]

m = mapview(polys1_sf, layer.name="Lake") +
  mapview(lake_buffer, layer.name="20km buffer")

if (!is.null(geometry_list$line)) {
  m = m + mapview(geometry_list$line, layer.name=paste0("Line_", i))
}
if (!is.null(geometry_list$points)) {
  m = m + mapview(geometry_list$points, layer.name=paste0("Points_", i), popup = geometry_list$points$popup)
}

m 

# NEED TO DO:
# observe the spatial relationship between the lake polygon, lake buffer, trajectory and points
# check if it is a circumnavigation event
# need to record "idx" for the starting point and the ending point for each circumnavigation event
# update your "result_df" dataframe or on your google sheet, based on your manual inspection and filter out any case which is not a circumnavigation event
# new "result_df" should have 5 columns: "Year", "ID", "season", "before_index", "after_index"

# =========================
# Step 3: Compute metrics (APR, speed, etc) for confirmed circumnavigation events
# =========================
#load albedo data for one year
#check working directory at first
load_albedo_one_year = function(year, polys1_sp, mat_dir = "data") {
  data_filename = file.path(mat_dir, paste0("output_", year, ".mat"))
  if (!file.exists(data_filename)) {
    stop("Cannot find: ", data_filename)
  }
  
  climatology_data = readMat(data_filename)
  arr_climatology = array(unlist(climatology_data), dim = c(426, 724, 365))
  
  arr_climatology_brick = brick(
    arr_climatology,
    xmn = -111.7125, xmx = -108.6958,
    ymn = 64.375,   ymx = 66.15,
    crs = "+proj=longlat +datum=WGS84",
    transpose = FALSE
  )
  
  arr_climatology_brick_crop = crop(arr_climatology_brick, extent(polys1_sp))
  albedo_data = mask(arr_climatology_brick_crop, polys1_sp)
  
  return(albedo_data)  # RasterBrick with 365 layers
}

# loop
# results_df_confirmed must contain:
# Year, ID, season, before_index, after_index
  
pair_list = list()

current_year = NA_integer_
albedo_data = NULL
percentile_rank = NULL

# results_df_checked should contain manually confirmed circumnavigation events
# save it as circumnavigation_candidates_checked.csv
# with columns: Year, ID, season, before_index, after_index
results_df_checked = read.csv("circumnavigation_candidates_checked.csv", stringsAsFactors = FALSE)

for (i in seq_len(nrow(results_df_checked))) {
  
  year  = results_df_checked$Year[i]
  id    = as.character(results_df_checked$ID[i])
  season_i = as.character(results_df_checked$season[i])
  before_index = results_df_checked$before_index[i]
  after_index  = results_df_checked$after_index[i]
  
  # load albedo + calc percentile_rank 
  if (is.na(current_year) || year != current_year) {
    
    current_year = year

    # load albedo for this year 
    albedo_data = load_albedo_one_year(
      year = year,
      polys1_sp = polys1,     # your Spatial lake polygon
      mat_dir = "data"        # folder where output_YYYY.mat lives
    )
    
    percentile_rank = calc(albedo_data[[92:280]], fun = function(x) {
      n = sum(!is.na(x))
      if (n == 0) return(rep(NA_real_, length(x)))
      rank(x, na.last = "keep") / n
    })
  }
  
  data = nwt_lakes %>%
    filter(Year == year, ID == id) %>%
    arrange(Time) %>%
    mutate(yday = lubridate::yday(Time))
  
  y1 = season_windows[[season_i]][1]
  y2 = season_windows[[season_i]][2]
  whole_line = subset(data, yday >= y1 & yday <= y2)
  
  if (nrow(whole_line) < 2) next
  
  if (is.na(before_index) || is.na(after_index)) next
  if (before_index < 1 || after_index < 1) next
  if (before_index > nrow(whole_line) || after_index > nrow(whole_line)) next
  if (before_index >= after_index) next
  
  nearest_before = whole_line[before_index, ]
  nearest_after  = whole_line[after_index, ]
  crossing_date = nearest_before$yday
  
  #lake_width: The length of the intersection between the lake polygon and the line connecting consecutive sampling points.  
  point1 = whole_line[before_index, c("Lon", "Lat")]
  point2 = whole_line[after_index, c("Lon", "Lat")]
  
  line = Line(rbind(point1, point2))
  lines = Lines(list(line), ID="1")
  whole_line_vct_lines = SpatialLines(list(lines))
  
  proj4string(whole_line_vct_lines) = CRS('+proj=longlat +datum=WGS84 +no_defs')
  proj4string(polys1) = CRS('+proj=longlat +datum=WGS84 +no_defs')
  
  sfLines = st_as_sf(whole_line_vct_lines)
  sfPolygons = st_as_sf(polys1)
  
  intersectionPoints = st_intersection(sfLines,sfPolygons)
  lake_width = if (nrow(intersectionPoints) == 0) NA_real_ else as.numeric(st_length(intersectionPoints))
  
  # straight distance/speed
  nearest_before_location = c(nearest_before$Lon[[1]], nearest_before$Lat[[1]])
  nearest_after_location  = c(nearest_after$Lon[[1]],  nearest_after$Lat[[1]])
  straight_distance=geosphere::distm(nearest_before_location,nearest_after_location, fun = distHaversine)[1,1]
  time_difference = nearest_after$Time - nearest_before$Time
  time_difference_seconds = as.numeric(time_difference, units = "secs")
  crossing_duration=as.numeric(time_difference, units = "days")
  
  straight_speed= straight_distance/time_difference_seconds

  # circumnavigate distance/speed ((projected to 32612)
  lake_proj = sf::st_transform(polys1_sf, 32612)
  
  pt1 = sf::st_as_sf(
    data.frame(Lon = nearest_before$Lon[[1]], Lat = nearest_before$Lat[[1]]),
    coords = c("Lon", "Lat"),
    crs = 4326
  ) %>% sf::st_transform(32612)
  
  pt2 = sf::st_as_sf(
    data.frame(Lon = nearest_after$Lon[[1]], Lat = nearest_after$Lat[[1]]),
    coords = c("Lon", "Lat"),
    crs = 4326
  ) %>% sf::st_transform(32612)
  
  # compute domain = bbox(lake + pts) buffered
  bbox_proj = sf::st_as_sfc(sf::st_bbox(sf::st_union(lake_proj, sf::st_union(pt1, pt2)))) %>%
    sf::st_buffer(dist = 60000)
  coast = sf::st_difference(bbox_proj, sf::st_union(lake_proj))
  
  # raster mask: land=1, lake=NA
  # optional: res can be different from 100m, 250m to 500m (same as albedo pixel)
  r = raster::raster(raster::extent(coast), res = 100, crs = sf::st_crs(lake_proj)$wkt)
  r[] = 1
  
  # rasterize lake onto r
  lake_sf = sf::st_as_sf(sf::st_union(lake_proj))
  lake_sf$val = 1
  lake_raster = fasterize::fasterize(lake_sf, r, field = "val")
  
  # land_mask: lake cells become NA (blocked), land stays 1
  land_mask = raster::mask(r, lake_raster, inverse = TRUE)
  land_mask[land_mask < 0] = NA
  
  # transition + geoCorrection (projected CRS)
  tr = gdistance::transition(land_mask,
                             transitionFunction = function(x) 1,
                             directions = 16,
                             symm = TRUE)
  tr = gdistance::geoCorrection(tr, type = "c", scl = FALSE)
  
  # cost distance (meters)
  cost = gdistance::costDistance(tr, 
                                 sf::st_coordinates(pt1), 
                                 sf::st_coordinates(pt2))
  circumvent_distance = as.numeric(cost[1,1])
  
  circumvent_speed=circumvent_distance/time_difference_seconds
  
  # extract the albedo along the path
  rank_number = crossing_date - 91
  if (rank_number < 1 || rank_number > nlayers(percentile_rank)) next
  
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
  
  # albedo percentile rank (APR) at the nearest pixels
  # option1: find out the last lake pixel before the individual is back to the land
  layer = intersected_pixels$layer
  
  if (all(is.na(layer))) {
    albedo_nearest_pixel = NA_real_
  } else {
    # lake to land: from non NA to NA
    water_to_land = which(!is.na(layer[-length(layer)]) & is.na(layer[-1]))
    if (length(water_to_land) == 0) {
      albedo_nearest_pixel = tail(na.omit(layer), 1)
    } else {
      albedo_nearest_pixel = layer[max(water_to_land)]
    }
  }
  
  # option2: the first lake pixel after the individual enter into lake
  if (all(is.na(layer))) {
    albedo_first_water = NA_real_
  } else {
    # land to lake: from NA to non-NA
    land_to_water = which(is.na(layer[-length(layer)]) & !is.na(layer[-1]))
    
    if (length(land_to_water) == 0) {
      albedo_first_water = layer[which(!is.na(layer))[1]]
    } else {
      albedo_first_water = layer[min(land_to_water) + 1]
    }
  }
  
  #split data frame results into multiple parts based on continuous NA rows
  #na_inds = which(is.na(intersected_pixels$layer))
  #start_inds = c(1, na_inds[which(diff(na_inds) != 1)] + 1)
  #end_inds = c(na_inds[which(diff(na_inds) != 1)], nrow(intersected_pixels))
  #split_dfs = lapply(seq_along(start_inds), function(i) intersected_pixels[start_inds[i]:end_inds[i], ])
  # albedo percentile rank (APR) at the nearest pixels
  #near_albedo=split_dfs[length(split_dfs)]%>%data.frame()%>%dplyr::select('layer') #option:split_dfs[[1]] %>% data.frame()
  #albedo_nearest_pixel=near_albedo[(which(is.na(near_albedo))[1]-1),]
  
  # albedo percentile rank (APR) along the crossing path
  albedo_linear=intersected_pixels$layer %>% mean(,na.rm=TRUE)
  
  # albedo percentile rank (APR) of the entire lake
  x.stats = data.frame(x.mean=cellStats(albedo_data[[92:280]], "mean",na.rm=TRUE))
  albedo_whole_lake=(rank(x.stats$x.mean)/sum(!is.na(x.stats)))[(rank_number)]
  
  pair_list[[length(pair_list) + 1]] = list(
    Year = year,
    ID = id,
    crossing_time = nearest_before$Time,
    crossing_duration = crossing_duration,
    crossing_date = nearest_before$yday,
    nearest_before_Lon = nearest_before$Lon[[1]],
    nearest_before_Lat = nearest_before$Lat[[1]],
    nearest_after_Lon  = nearest_after$Lon[[1]],
    nearest_after_Lat  = nearest_after$Lat[[1]],
    lake_width = lake_width,
    circumvent_distance = circumvent_distance,
    straight_distance = straight_distance,
    circumvent_speed = circumvent_speed,
    straight_speed = straight_speed,
    albedo_linear = albedo_linear,
    albedo_nearest_pixel = albedo_nearest_pixel,
    albedo_whole_lake = albedo_whole_lake,
    type = 0,
    season = season_i
  )
  
  # reference rows: ±3 segments around the focal event
  reference_list = list()
  
  reference_offsets = c(-3, -2, -1, 1, 2, 3)
  
  for (offset in reference_offsets) {
    
    if (offset < 0) {
      k = before_index + offset
    } else {
      k = after_index + offset - 1
    }
    
    if (k > 0 && k < nrow(whole_line)) {
      nearest_location = whole_line[k, ]
      nearest_location_next = whole_line[k + 1, ]
      
      pA = c(nearest_location$Lon[[1]], nearest_location$Lat[[1]])
      pB = c(nearest_location_next$Lon[[1]], nearest_location_next$Lat[[1]])
      
      displacement = geosphere::distm(pA, pB, fun = geosphere::distHaversine)[1,1]
      td1 = nearest_location_next$Time - nearest_location$Time
      td1_seconds = as.numeric(td1, units = "secs")
      crossing_duration1 = as.numeric(td1, units = "days")
      
      average_speed = ifelse(td1_seconds > 0, displacement / td1_seconds, NA)
      
      reference_list[[length(reference_list) + 1]] = list(
        Year = year,
        ID = id,
        crossing_time = nearest_location$Time,
        crossing_duration = crossing_duration1,
        crossing_date = nearest_location$yday,
        nearest_before_Lon = nearest_location$Lon[[1]],
        nearest_before_Lat = nearest_location$Lat[[1]],
        nearest_after_Lon  = nearest_location_next$Lon[[1]],
        nearest_after_Lat  = nearest_location_next$Lat[[1]],
        lake_width = NA,
        circumvent_distance = NA,
        straight_distance = displacement,
        circumvent_speed = NA,
        straight_speed = average_speed,
        albedo_linear = NA,
        albedo_nearest_pixel = NA,
        albedo_whole_lake = NA,
        type = offset,
        season = season_i
      )
    }
  }
  
  pair_list = append(pair_list, reference_list)
}

pair_list_df = lapply(pair_list, function(x) as.data.frame(t(unlist(x))))
result_df_circumnavigate = do.call(rbind, pair_list_df)

result_df_circumnavigate$crossing_time = as.POSIXct(
  as.numeric(result_df_circumnavigate$crossing_time), 
  origin = "1970-01-01", tz = "UTC")

#save
write.csv(
  result_df_circumnavigate,
  "circumnavigate_event.csv",
  row.names = FALSE
)

# Optional: export results to Google Sheets
# library(googlesheets4)
# result_df_circumnavigate %>%
#   write_sheet(
#     ss = gs4_get("link"),
#     sheet = "circumnavigate_event"
#   )


# Optional: rebuild output after standardizing columns across list elements
# in case some rows are missing fields (e.g., albedo-related values for references)

# if the albedo_linear is missing or have NA
standard_columns = colnames(result_df_circumnavigate)  

pair_list_df_standardized = lapply(pair_list_df, function(df) {

  for (col in standard_columns) {
    if (!col %in% names(df)) {
      df[[col]] = NA  
    }
  }

  df = df[, standard_columns]
  return(df)
})

result_df1 = do.call(rbind, pair_list_df_standardized)

#result_df = do.call(rbind, pair_list_df)
result_df1$crossing_time = as.POSIXct(as.numeric(result_df1$crossing_time), origin = "1970-01-01", tz = "UTC")



