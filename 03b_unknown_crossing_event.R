# ----------------------------------------
# 03b_unknown_crossing_event.R
# Purpose:
#   (1) Identify candidate unknown crossing events based on lake-line intersection,
#   (2) manually review candidates,
#   (3) compute event-level movement and albedo/APR metrics for checked unknown events.
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
#   - manually checked candidate table with Year, ID, season
# ----------------------------------------

library(lubridate)
library(dplyr)
library(terra)
library(sf)
library(R.matlab)
library(gdistance)
library(fasterize)
library(raster)
library(geosphere)

# load caribou movement data
load("data/nwt_lakes.RData")
nwt_lakes = st_transform(nwt_lakes, 4326) %>% arrange(Time)

# lake polygon
polys1_sf = sf::st_read("data/Contwoyto_pure_water.shp", quiet = TRUE)
polys1_sf = sf::st_transform(polys1_sf, 4326)
polys1 = as(polys1_sf, "Spatial")
polygon_lake = terra::vect(polys1_sf)

season_windows = list(
  spring = c(92, 181), #spring migration season (DOY 92 - DOY 181)
  fall   = c(182, 280) #fall migration season (DOY 182 - DOY 280)
)


#### 这里记得加上mat_dir = "data"
# load one year of albedo
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
    ymn = 64.375, ymx = 66.15,
    crs = "+proj=longlat +datum=WGS84",
    transpose = FALSE
  )
  
  arr_climatology_brick_crop = crop(arr_climatology_brick, extent(polys1_sp))
  albedo_data = mask(arr_climatology_brick_crop, polys1_sp)
  
  return(albedo_data)
}

# =========================
# Step 1: Generate broad candidate table
# =========================
actual_combinations = nwt_lakes %>%
  st_drop_geometry() %>%
  distinct(Year, ID) %>%
  arrange(Year, ID)

results = list()

for (i in seq_len(nrow(actual_combinations))) {
  
  year_i = actual_combinations$Year[i]
  id_i   = as.character(actual_combinations$ID[i])
  
  data_i = nwt_lakes %>%
    filter(Year == year_i, ID == id_i) %>%
    arrange(Time)
  
  if (nrow(data_i) < 2) next
  
  whole_line_year = data.frame(data_i) %>%
    arrange(Time) %>%
    mutate(yday = lubridate::yday(Time))
  
  for (ss in names(season_windows)) {
    
    y1 = season_windows[[ss]][1]
    y2 = season_windows[[ss]][2]
    
    whole_line = subset(whole_line_year, yday >= y1 & yday <= y2)
    if (nrow(whole_line) < 2) next
    
    pts_vct = terra::vect(
      whole_line,
      geom = c("Lon", "Lat"),
      crs = "EPSG:4326",
      keepgeom = TRUE
    )
    
    whole_line_vct_line = terra::as.lines(pts_vct)
    intersection_lake = terra::intersect(polygon_lake, whole_line_vct_line)
    
    if (terra::is.empty(intersection_lake) || nrow(intersection_lake) == 0) next
    
    results[[length(results) + 1]] = data.frame(
      Year   = year_i,
      ID     = id_i,
      season = ss,
      stringsAsFactors = FALSE
    )
  }
}

lake_intersection_candidates = if (length(results) == 0) {
  data.frame(Year = integer(), ID = character(), season = character())
} else {
  do.call(rbind, results) %>%
    distinct(Year, ID, season) %>%
    arrange(Year, ID, season)
}

write.csv(
  lake_intersection_candidates,
  "lake_intersection_candidates.csv",
  row.names = FALSE
)


# =========================
# Manual inspection - MUST DO!
# =========================
# Review lake_intersection_candidates.csv and retain only true unknown-crossing cases.
# Record:
#   - Year
#   - ID
#   - season
#   - before_index
#   - after_index
#
# Save as:
#   unknown_crossing_candidates_checked.csv
#
# Note:
# before_index and after_index refer to row indices within the seasonal subset
# (whole_line), not within the full annual trajectory




# =========================
# Step 1.5: Manual review support
# Purpose:
#   Visualize a single Year-ID-season case, label seasonal points with idx,
#   and suggest candidate before_index / after_index values based on
#   trajectory-lake intersection geometry.
# =========================

library(mapview)

findNearestPointBefore = function(start_lon, start_lat, path) {
  path_xy = cbind(
    as.numeric(unlist(path$Lon)),
    as.numeric(unlist(path$Lat))
  )
  
  distances = geosphere::distm(
    x = matrix(c(start_lon, start_lat), ncol = 2),
    y = path_xy,
    fun = geosphere::distHaversine
  )
  
  min_index_start = which.min(distances)
  
  if (min_index_start > 1) {
    return(min_index_start)
  } else {
    return(NULL)
  }
}

findNearestPointAfter = function(end_lon, end_lat, path, nearest_before_index) {
  if (is.null(nearest_before_index)) return(NULL)
  
  path_xy = cbind(
    as.numeric(unlist(path$Lon)),
    as.numeric(unlist(path$Lat))
  )
  
  distances = geosphere::distm(
    x = matrix(c(end_lon, end_lat), ncol = 2),
    y = path_xy,
    fun = geosphere::distHaversine
  )
  
  min_index_end = which.min(distances)
  
  if (min_index_end < nrow(path) && min_index_end > nearest_before_index) {
    return(min_index_end)
  } else {
    return(NULL)
  }
}

# help check the location of the crossing event
suggest_crossing_indices = function(year_i, id_i, season_i,
                                    nwt_lakes, polygon_lake, polys1_sf,
                                    season_windows) {
  
  data_i = nwt_lakes %>%
    filter(Year == year_i, ID == id_i) %>%
    arrange(Time) %>%
    mutate(yday = yday(Time))
  
  if (nrow(data_i) < 2) {
    message("Not enough points.")
    return(NULL)
  }
  
  y1 = season_windows[[season_i]][1]
  y2 = season_windows[[season_i]][2]
  
  whole_line = subset(data_i, yday >= y1 & yday <= y2)
  if (nrow(whole_line) < 2) {
    message("Not enough points in this season window.")
    return(NULL)
  }
  
  whole_line_df = sf::st_drop_geometry(whole_line)
  
  pts_vct = terra::vect(
    whole_line_df,
    geom = c("Lon", "Lat"),
    crs = "EPSG:4326",
    keepgeom = TRUE
  )
  
  whole_line_vct_line = terra::as.lines(pts_vct)
  intersection_points = terra::intersect(polygon_lake, whole_line_vct_line)
  
  if (terra::is.empty(intersection_points) || nrow(intersection_points) == 0) {
    message("No lake intersection found.")
    return(NULL)
  }
  
  gmat = geom(intersection_points[1, ])
  total_segments = max(gmat[, "part"])
  
  suggestions = list()
  
  for (g in seq_len(total_segments)) {
    seg_rows = which(gmat[, "part"] == g)
    seg_geom = gmat[seg_rows, , drop = FALSE]
    
    start_lon = seg_geom[1, "x"]
    start_lat = seg_geom[1, "y"]
    end_lon   = seg_geom[nrow(seg_geom), "x"]
    end_lat   = seg_geom[nrow(seg_geom), "y"]
    
    before_index = findNearestPointBefore(start_lon, start_lat, whole_line_df)
    after_index  = findNearestPointAfter(end_lon, end_lat, whole_line_df, before_index)
    
    if (!is.null(before_index) && !is.null(after_index) && after_index > before_index) {
      suggestions[[length(suggestions) + 1]] = data.frame(
        Year = year_i,
        ID = id_i,
        season = season_i,
        segment = g,
        before_index = before_index,
        after_index = after_index,
        before_time = whole_line$Time[before_index],
        after_time = whole_line$Time[after_index],
        stringsAsFactors = FALSE
      )
    }
  }
  
  if (length(suggestions) == 0) {
    message("No valid before/after index pair found.")
    return(NULL)
  }
  
  suggestions_df = do.call(rbind, suggestions)
  
  # build point map with idx labels in popup
  whole_line_sf = st_as_sf(whole_line_df, coords = c("Lon", "Lat"), crs = 4326, remove = FALSE)
  whole_line_sf$idx = seq_len(nrow(whole_line_sf))
  whole_line_sf$popup = paste0(
    "idx: ", whole_line_sf$idx,
    "<br>Time: ", whole_line_sf$Time,
    "<br>yday: ", whole_line_sf$yday
  )
  
  line_sf = whole_line_sf %>%
    summarise(geometry = st_combine(geometry)) %>%
    st_cast("LINESTRING")
  
  map_obj = mapview(polys1_sf, col.regions = "white", alpha.regions = 0.3, layer.name = "Lake") +
    mapview(line_sf, color = "red", lwd = 2, layer.name = "Track") +
    mapview(whole_line_sf, popup = whole_line_sf$popup, layer.name = "Points")
  
  print(map_obj)
  print(suggestions_df)
  
  return(list(
    suggestions = suggestions_df,
    whole_line = whole_line,
    map = map_obj
  ))
}

# =========================
# Example manual review calls to inspect individual cases
# =========================
res_bat077 = suggest_crossing_indices(
  year_i = 2010,
  id_i = "bat077",
  season_i = "spring",
  nwt_lakes = nwt_lakes,
  polygon_lake = polygon_lake,
  polys1_sf = polys1_sf,
  season_windows = season_windows
)

# =========================
# Manual inspection required
# =========================
# Use suggest_crossing_indices() to inspect candidate cases one by one.
# Record the final confirmed indices in:
#   unknown_crossing_candidates_checked.csv
#
# Required columns:
#   Year, ID, season, before_index, after_index
#
# Note:
# before_index and after_index refer to row indices within the seasonal subset
# (whole_line), not within the full annual trajectory.


# =========================
# Part 2. Build unknown crossing event table
# =========================
initial_list = read.csv("unknown_crossing_candidates_checked.csv", stringsAsFactors = FALSE)

pair_list = list()

current_year = NA_integer_
albedo_data = NULL
percentile_rank = NULL

for (i in seq_len(nrow(initial_list))) {
  
  year = initial_list$Year[i]
  id = as.character(initial_list$ID[i])
  season_i = as.character(initial_list$season[i])
  before_index = initial_list$before_index[i]
  after_index = initial_list$after_index[i]
  
  # load annual albedo data once per year
  if (is.na(current_year) || year != current_year) {
    current_year = year
    
    albedo_data = load_albedo_one_year(
      year = year,
      polys1_sp = polys1,
      mat_dir = "data"
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
  
  whole_line_df = sf::st_drop_geometry(whole_line)
  
  if (is.na(before_index) || is.na(after_index)) next
  if (before_index < 1 || after_index < 1) next
  if (before_index > nrow(whole_line_df) || after_index > nrow(whole_line_df)) next
  if (before_index >= after_index) next
  
  nearest_before = whole_line_df[before_index, ]
  nearest_after  = whole_line_df[after_index, ]
  crossing_date = nearest_before$yday
  
  # focal season label
  season_value = season_i
  
  # lake width: intersection length between lake polygon and straight segment
  point1 = c(nearest_before$Lon[[1]], nearest_before$Lat[[1]])
  point2 = c(nearest_after$Lon[[1]], nearest_after$Lat[[1]])
  
  line = Line(rbind(point1, point2))
  lines = Lines(list(line), ID = "1")
  whole_line_vct_lines = SpatialLines(list(lines))
  
  proj4string(whole_line_vct_lines) = CRS("+proj=longlat +datum=WGS84 +no_defs")
  proj4string(polys1) = CRS("+proj=longlat +datum=WGS84 +no_defs")
  
  sfLines = st_as_sf(whole_line_vct_lines)
  sfPolygons = st_as_sf(polys1)
  
  intersectionPoints = st_intersection(sfLines, sfPolygons)
  lake_width = if (nrow(intersectionPoints) == 0) NA_real_ else as.numeric(st_length(intersectionPoints))
  
  # straight distance / speed
  nearest_before_location = c(nearest_before$Lon[[1]], nearest_before$Lat[[1]])
  nearest_after_location  = c(nearest_after$Lon[[1]], nearest_after$Lat[[1]])
  
  straight_distance = geosphere::distm(
    nearest_before_location,
    nearest_after_location,
    fun = geosphere::distHaversine
  )[1,1]
  
  time_difference = nearest_after$Time - nearest_before$Time
  time_difference_seconds = as.numeric(time_difference, units = "secs")
  crossing_duration = as.numeric(time_difference, units = "days")
  straight_speed = straight_distance / time_difference_seconds
  
  # circumnavigate distance / speed
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
  
  bbox_proj = sf::st_as_sfc(sf::st_bbox(sf::st_union(lake_proj, sf::st_union(pt1, pt2)))) %>%
    sf::st_buffer(dist = 60000)
  
  coast = sf::st_difference(bbox_proj, sf::st_union(lake_proj))
  
  coast_bbox = sf::st_bbox(coast)
  
  r = raster::raster(
    xmn = coast_bbox["xmin"],
    xmx = coast_bbox["xmax"],
    ymn = coast_bbox["ymin"],
    ymx = coast_bbox["ymax"],
    res = 100,
    crs = sf::st_crs(lake_proj)$wkt
  )
  r[] = 1
  
  lake_sf = sf::st_as_sf(sf::st_union(lake_proj))
  lake_sf$val = 1
  lake_raster = fasterize::fasterize(lake_sf, r, field = "val")
  
  land_mask = raster::mask(r, lake_raster, inverse = TRUE)
  land_mask[land_mask < 0] = NA
  
  tr = gdistance::transition(
    land_mask,
    transitionFunction = function(x) 1,
    directions = 16,
    symm = TRUE
  )
  tr = gdistance::geoCorrection(tr, type = "c", scl = FALSE)
  
  cost = gdistance::costDistance(
    tr,
    sf::st_coordinates(pt1),
    sf::st_coordinates(pt2)
  )
  circumvent_distance = as.numeric(cost[1,1])
  circumvent_speed = circumvent_distance / time_difference_seconds
  
  # APR extraction
  rank_number = crossing_date - 91
  if (rank_number < 1 || rank_number > nlayers(percentile_rank)) next
  
  albedo_drop_percentile = percentile_rank[[rank_number]]
  test2 = raster(albedo_drop_percentile)
  values(test2) = values(albedo_drop_percentile)
  albedo_drop_percentile_rast = terra::rast(test2)
  
  line_before_after = whole_line_df %>%
    filter(between(Time, nearest_before$Time, nearest_after$Time))
  
  line_before_after_vct = terra::vect(
    line_before_after,
    geom = c("Lon", "Lat"),
    crs = "+proj=longlat +datum=WGS84 +no_defs +type=crs",
    keepgeom = TRUE
  ) %>% as.lines()
  
  intersected_pixels = terra::extract(albedo_drop_percentile_rast, line_before_after_vct)
  layer_vals = intersected_pixels$layer
  
  # APR along the path
  albedo_linear = mean(layer_vals, na.rm = TRUE)
  
  # APR at the nearest lake pixel
  albedo_nearest_pixel = NA_real_
  if (!all(is.na(layer_vals))) {
    water_to_land = which(!is.na(layer_vals[-length(layer_vals)]) & is.na(layer_vals[-1]))
    
    if (length(water_to_land) > 0) {
      albedo_nearest_pixel = layer_vals[max(water_to_land)]
    } else {
      albedo_nearest_pixel = tail(na.omit(layer_vals), 1)
    }
  }
  
  # APR of the entire lake
  x.stats = data.frame(x.mean = cellStats(albedo_data[[92:280]], "mean", na.rm = TRUE))
  albedo_whole_lake = (rank(x.stats$x.mean) / sum(!is.na(x.stats)))[rank_number]
  
  # focal row
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
    season = season_value
  )
  
  # reference rows
  reference_offsets = c(-3, -2, -1, 1, 2, 3)
  
  for (offset in reference_offsets) {
    
    if (offset < 0) {
      k = before_index + offset
    } else {
      k = after_index + offset -1
    }
    
    if (k > 0 && k < nrow(whole_line_df)) {
      nearest_location = whole_line_df[k, ]
      nearest_location_next = whole_line_df[k + 1, ]
      
      pA = c(nearest_location$Lon[[1]], nearest_location$Lat[[1]])
      pB = c(nearest_location_next$Lon[[1]], nearest_location_next$Lat[[1]])
      
      displacement = geosphere::distm(pA, pB, fun = geosphere::distHaversine)[1,1]
      td1 = nearest_location_next$Time - nearest_location$Time
      td1_seconds = as.numeric(td1, units = "secs")
      crossing_duration_ref = as.numeric(td1, units = "days")
      
      average_speed = ifelse(td1_seconds > 0, displacement / td1_seconds, NA)
      
      pair_list[[length(pair_list) + 1]] = list(
        Year = year,
        ID = id,
        crossing_time = nearest_location$Time,
        crossing_duration = crossing_duration_ref,
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
        season = season_value
      )
    }
  }
}

pair_list_df = lapply(pair_list, function(x) as.data.frame(t(unlist(x))))
result_df_unknown = do.call(rbind, pair_list_df)

result_df_unknown$crossing_time = as.POSIXct(
  as.numeric(result_df_unknown$crossing_time),
  origin = "1970-01-01",
  tz = "UTC"
)

write.csv(
  result_df_unknown,
  "unknown_crossing_event.csv",
  row.names = FALSE
)
