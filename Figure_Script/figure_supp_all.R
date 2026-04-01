library(dplyr)
library(randomForest)
library(caret)
library(ROCR)
library(pROC)
library(MuMIn)
library(car)
library(glmnet)
library(MASS)
library(ggplot2)
library(sf)
library(terra)
library(raster)
library(R.matlab)
library(dplyr)
library(ggplot2)

# ========================================
# Appendix_S1_Figure_S1_albedo_scatter.R
# Purpose:
#   Plot pixel-wise albedo values from Contwoyto Lake in 2019
#   from day 92 (April 2) to day 280 (October 7).
# ========================================

# =========================
# 1. Load lake polygon
# =========================
polys1_sf = st_read("data/Contwoyto_pure_water.shp", quiet = TRUE)
polys1_sf = st_transform(polys1_sf, 4326)
polys1 = as(polys1_sf, "Spatial")

# =========================
# 2. Load 2019 albedo data
# =========================
year = 2019
data_filename = file.path("data", paste0("output_", year, ".mat"))

climatology_data = readMat(data_filename)
arr_climatology = array(unlist(climatology_data), dim = c(426, 724, 365))

arr_climatology_brick = brick(
  arr_climatology,
  xmn = -111.7125,
  xmx = -108.6958,
  ymn = 64.375,
  ymx = 66.15,
  crs = "+proj=longlat +datum=WGS84",
  transpose = FALSE
)

# crop and mask to lake
arr_climatology_brick_crop = crop(arr_climatology_brick, extent(polys1))
albedo_data = mask(arr_climatology_brick_crop, polys1)

# =========================
# 3. Keep only DOY 92-280
# =========================
albedo_subset = albedo_data[[92:280]]

# =========================
# 4. Convert raster layers to long dataframe
# =========================
albedo_list = vector("list", nlayers(albedo_subset))

for (i in seq_len(nlayers(albedo_subset))) {
  layer_data = raster::raster(albedo_subset, layer = i)
  layer_values = raster::values(layer_data)
  
  albedo_list[[i]] = data.frame(
    Day = i + 91,
    Albedo = layer_values
  )
}

albedo_df = bind_rows(albedo_list) %>%
  filter(!is.na(Albedo))

# =========================
# 5. Plot
# =========================
p_s1 = ggplot(albedo_df, aes(x = Day, y = Albedo)) +
  geom_point(alpha = 0.1, size = 0.5) +
  theme_classic() +
  labs(
    x = "Day of the Year",
    y = "Albedo"
  )

print(p_s1)

# ========================================
# Appendix S2 related exploratory code
# Convex hull approaches for shortest
# circumnavigation paths
#
# Main manuscript analyses use gdistance.
# Code below preserves alternative methods
# provided by Elie Gurarie and referenced
# in Appendix S2.
# ========================================


# ----------------------------------------
# Part A. Single-lake version
# Used to adapt the method to the Contwoyto
# Lake system after removing interior rings
# / islands.
# ----------------------------------------

getOuterPolygon = function(poly) {
  xy.poly = st_coordinates(poly) |> data.frame()
  st_sfc(
    st_polygon(
      xy.poly |> subset(L1 == 1) |> as.matrix() |> list()
    ),
    crs = st_crs(poly)
  )
}

circumnavigateLake = function(start, end, lake,
                               plotme = TRUE,
                               returntrack = FALSE,
                               verbose = FALSE) {
  if (inherits(start, "sf")) {
    start = st_coordinates(start)
    start.z = start[1] + 1i * start[2]
  } else {
    start.z = start
  }
  
  if (inherits(end, "sf")) {
    end = st_coordinates(end)
    end.z = end[1] + 1i * end[2]
  } else {
    end.z = end
  }
  
  if (inherits(lake, "sf") | inherits(lake, "sfc")) {
    lake = st_coordinates(lake)
  }
  
  if (start.z == end.z) {
    distance = NA
    return(distance)
  }
  
  z.track = c(start.z, end.z)
  z.lake = lake[, 1] + 1i * lake[, 2]
  z.lake.step = cbind(z.lake[-length(z.lake)], z.lake[-1])
  
  which.cross = which(apply(z.lake.step, 1, isIntersect, z.track))
  
  if (length(which.cross) == 0) {
    if (verbose) {
      cat("Straightest line does not touch the lake.\n")
    }
    if (returntrack) {
      return(list(
        distance = Mod(diff(z.track)),
        landtrack = cbind(x = Re(z.track), y = Im(z.track))
      ))
    } else {
      return(Mod(diff(z.track)))
    }
  } else {
    which.cross = which.cross[c(1, length(which.cross))]
    if (which.cross[1] > which.cross[2]) {
      which.cross = sort(which.cross)
      z.track = z.track[2:1]
    }
    
    z.tight.poly = c(z.track[1], z.lake[(which.cross[1] + 1):which.cross[2]], z.track[2])
    
    mcp.z.single = function(z.poly) {
      z.mcp = z.poly[cbind(Re(z.poly), Im(z.poly)) |> chull()]
      z.mcp
    }
    
    z.mcp = z.tight.poly |> mcp.z.single()
    
    which.start = which(z.mcp == start.z)
    which.end = which(z.mcp == end.z)
    which.1 = min(which.start, which.end)
    which.2 = max(which.start, which.end)
    
    if (is.infinite(which.1)) {
      which.1 = 0
    }
    if (is.infinite(which.2)) {
      which.2 = 0
    }
    
    z.track.tight = z.mcp[which.1:length(z.mcp)]
    
    if (plotme) {
      plot(
        lake,
        asp = 1,
        type = "l",
        ylim = range(Im(c(z.lake, z.track))),
        xlim = range(Re(c(z.lake, z.track)))
      )
      polygon(lake, col = "grey")
      lines(z.track.tight, col = 2, pch = 19, type = "o")
      lines(z.mcp, col = 2, lwd = 2)
      lines(z.tight.poly, col = "blue", lwd = 1)
      points(rbind(start, end), pch = 19)
    }
    
    landtrack = cbind(x = Re(z.track.tight), y = Im(z.track.tight))
    distance = sum(Mod(diff(z.track.tight)))
    
    if (returntrack) {
      return(list(
        distance = distance,
        landtrack = landtrack
      ))
    } else {
      return(distance)
    }
  }
}

isIntersect = function(z1, z2) {
  theta = Arg(z1[2] - z1[1])
  z1.t = (z1 - z1[1]) * exp(-1i * theta)
  z2.t = (z2 - z1[1]) * exp(-1i * theta)
  slope.r = (Re(z2.t[2]) - Re(z2.t[1])) / (Im(z2.t[2]) - Im(z2.t[1]))
  x.intercept = Re(z2.t[1]) - slope.r * Im(z2.t[1])
  intersect.x = x.intercept > 0 & x.intercept < Re(z1.t[2])
  intersect.x & (Im(z2.t[2]) * Im(z2.t[1])) < 0
}


# ----------------------------------------
# Part B. Multi-polygon version
# Exploratory version provided by EG for
# evaluating shortest paths around multiple
# irregular polygons and all combinations
# of clockwise / counterclockwise routing.
# ----------------------------------------

circumnavigatePoly = function(start, end, poly, clockwise = TRUE, ...) {
  if (inherits(start, "sf")) {
    start = st_coordinates(start)
    start.z = start[1] + 1i * start[2]
  } else {
    start.z = start
  }
  
  if (inherits(end, "sf")) {
    end = st_coordinates(end)
    end.z = end[1] + 1i * end[2]
  } else {
    end.z = end
  }
  
  if (inherits(poly, "sf") | inherits(poly, "sfc")) {
    poly = st_coordinates(poly)
    poly.z = poly[, 1] + 1i * poly[, 2]
  } else if (all(c("X", "Y") %in% colnames(poly))) {
    poly.z = poly[, "X"] + 1i * poly[, "Y"]
  } else {
    poly.z = poly
  }
  
  if (clockwise) {
    circumtrack = circumnavigatePoly_cc(start.z, end.z, poly.z, ...)
  } else {
    circumtrack.ccc = circumnavigatePoly_cc(end.z, start.z, poly.z, ...)
    circumtrack = circumtrack.ccc[nrow(circumtrack.ccc):1, ]
  }
  
  list(
    track = circumtrack,
    distance = sum(Mod(diff(circumtrack[, 1] + 1i * circumtrack[, 2])))
  )
}

circumnavigatePoly_cc = function(start.z, end.z, poly.z,
                                  plotme = FALSE,
                                  verbose = FALSE) {
  z.track = c(start.z, end.z)
  z.poly.step = cbind(poly.z[-length(poly.z)], poly.z[-1])
  
  which.cross = which(apply(z.poly.step, 1, is_intersection, z.track))
  
  if (length(which.cross) == 0) {
    if (verbose) {
      cat("Straightest line does not touch the poly.\n")
    }
    landtrack = cbind(x = Re(z.track), y = Im(z.track))
  } else {
    crossings = sapply(
      which.cross,
      function(i) find_intersection(z.poly.step[i, ], z.track)
    )
    
    crossing_order = order(Mod(crossings - start.z))
    
    enter = which.cross[crossing_order[1]]
    exit = which.cross[crossing_order[length(crossing_order)]]
    maxlen = length(poly.z)
    
    if (enter < exit) {
      z.tight.poly = c(z.track[1], poly.z[(enter + 1):exit], z.track[2])
    }
    if (exit < enter) {
      z.tight.poly = c(z.track[1], poly.z[c((enter + 1):maxlen, 1:exit)], z.track[2])
    }
    
    z.mcp = z.tight.poly |> mcp.z()
    z.mcp = c(z.mcp, z.mcp)
    
    which.start = min(which(z.mcp == start.z))
    which.end = which(z.mcp == end.z)
    which.end = which.end[which.end > which.start] |> min()
    
    z.track.tight = z.mcp[which.start:which.end]
    
    if (plotme) {
      plot(
        poly.z,
        asp = 1,
        type = "n",
        ylim = range(Im(c(poly.z, z.track))),
        xlim = range(Re(c(poly.z, z.track)))
      )
      polygon(poly.z, col = "grey", border = NA)
      lines(z.tight.poly, col = "blue", lwd = 1)
      lines(z.track.tight, col = "purple", lwd = 2, pch = 19, type = "o")
      points(rbind(start, end), pch = 19)
    }
    
    landtrack = cbind(x = Re(z.track.tight), y = Im(z.track.tight))
  }
  
  return(landtrack)
}

circumnavigateMultiPoly = function(start, end, polys,
                                    cc.sequence = rep(TRUE, length(polys)),
                                    plotme = TRUE) {
  if (inherits(start, "sf")) {
    start = st_coordinates(start)
    start.z = start[1] + 1i * start[2]
  } else {
    start.z = start
  }
  
  if (inherits(end, "sf")) {
    end = st_coordinates(end)
    end.z = end[1] + 1i * end[2]
  } else {
    end.z = end
  }
  
  if (inherits(polys, "sf") | inherits(polys, "sfc")) {
    polys.df = st_coordinates(polys) |> data.frame()
  }
  
  require(plyr)
  
  poly.list = dlply(polys.df, "L2")
  
  newtrack = circumnavigatePoly(start.z, end.z, polys[1, ], clockwise = cc.sequence[1])
  newtrack.z = newtrack$track[, 1] + 1i * newtrack$track[, 2]
  
  subtracks = list()
  subtracks[[1]] = newtrack.z[-length(newtrack.z)]
  
  for (i in 2:length(poly.list)) {
    lastpoly = st_coordinates(polys[length(subtracks), ])
    lastpoly.z = lastpoly[, "X"] + 1i * lastpoly[, "Y"]
    
    newstart = lastpoly.z[max(which(lastpoly.z %in% subtracks[[length(subtracks)]]))]
    newtrack = circumnavigatePoly(newstart, end, polys[i, ], clockwise = cc.sequence[i])
    
    newtrack.z = newtrack$track[, 1] + 1i * newtrack$track[, 2]
    
    if (nrow(newtrack$track) > 2) {
      subtracks[[length(subtracks) + 1]] = newtrack.z[-length(newtrack.z)]
    }
  }
  
  tighttrack = c(do.call("c", subtracks), end)
  
  z.mcp = tighttrack |> mcp.z()
  z.mcp = c(z.mcp, z.mcp)
  
  if (!any(!cc.sequence)) {
    which.start = min(which(z.mcp == start.z))
    which.end = which(z.mcp == end.z)
    which.end = which.end[which.end > which.start] |> min()
    circumtrack = z.mcp[which.start:which.end]
  } else if (!any(cc.sequence)) {
    which.start = which(z.mcp == start.z) |> max()
    which.end = which(z.mcp == end.z) |> min()
    circumtrack = z.mcp[which.start:which.end]
  } else {
    circumtrack = tighttrack
  }
  
  if (plotme) {
    plot(polys, border = NA, col = "grey", ylim = range(Im(c(start, end))))
    lines(rbind(start, end), type = "o")
    lines(tighttrack, col = "red", lwd = 2)
    lines(circumtrack, col = "purple", lwd = 2)
  }
  
  list(
    track = circumtrack,
    distance = sum(Mod(diff(circumtrack)))
  )
}

mcp.z = function(poly.z) {
  z.mcp = poly.z[cbind(Re(poly.z), Im(poly.z)) |> chull()]
  z.mcp
}

is_intersection = function(z1, z2) {
  theta = Arg(z1[2] - z1[1])
  z1.t = (z1 - z1[1]) * exp(-1i * theta)
  z2.t = (z2 - z1[1]) * exp(-1i * theta)
  slope.r = (Re(z2.t[2]) - Re(z2.t[1])) / (Im(z2.t[2]) - Im(z2.t[1]))
  x.intercept = Re(z2.t[1]) - slope.r * Im(z2.t[1])
  intersect.x = x.intercept > 0 & x.intercept < Re(z1.t[2])
  intersect.x & (Im(z2.t[2]) * Im(z2.t[1])) < 0
}

find_intersection = function(z1, z2) {
  theta = Arg(z1[2] - z1[1])
  z1.t = (z1 - z1[1]) * exp(-1i * theta)
  z2.t = (z2 - z1[1]) * exp(-1i * theta)
  
  crossing.t = Re(z2.t[1]) - (diff(Re(z2.t)) / diff(Im(z2.t))) * Im(z2.t[1])
  
  if (Re(crossing.t) > max(Re(z1.t)) |
      Re(crossing.t) < min(Re(z1.t)) |
      min(Im(z2.t)) > 0 |
      max(Im(z2.t)) < 0) {
    return(NULL)
  } else {
    return(crossing.t * complex(mod = 1, arg = theta) + z1[1])
  }
}


# ----------------------------------------
# Minimal example logic from EG's
# multi-polygon demonstration
# ----------------------------------------

#Example:
# library(sf)
#
# set.seed(4)
# n.sides = 20
# poly1 = 2i + complex(mod = rep(1, n.sides) + rnorm(8, sd = .2), 
#                       arg = seq(2 * pi, 0, length = n.sides))
# 
# poly1[1] = poly1[n.sides]
# poly2 = (poly1 - mean(poly1)) * exp(1i) * 1.2 + 5i + 1 
# poly3 = (poly1 - mean(poly1)) * exp(2i) * .8 + 10i - .1

# make_st = function(z) {
#   st_polygon(list(cbind(Re(z), Im(z))))
# }
#
# polys = lapply(list(poly1, poly2, poly3), make_st) |> st_sfc()
# start = 0 + 0i
# end = 0.5 + 15i
#
# cc.sequences = expand.grid(
#   replicate(length(polys), c(TRUE, FALSE), simplify = FALSE)
# )
#
# allcombinations = apply(cc.sequences, 1, function(cc) {
#   circumnavigateMultiPoly(start, end, polys, cc.sequence = cc, plotme = TRUE)
# })
#
# distances = sapply(allcombinations, function(x) x$distance)
# which.min(distances)
# cc.sequences[which.min(distances), ]


# ========================================
# Appendix S3 Table S1
# Number of movement events by GPS fix interval category
#
# GPS fix interval is defined as the median temporal interval
# (hours) of the six on-land reference steps surrounding each
# focal transit event, using the three steps before and three
# steps after the event.
# ========================================
library(dplyr)
library(tidyr)

# make sure crossing_duration is numeric and in hours
all_event_s3 = all_event %>%
  mutate(
    crossing_duration = as.numeric(crossing_duration),
    crossing_duration_hours = crossing_duration * 24
  )

# calculate event-level median GPS fix interval from reference steps
ref_duration = all_event_s3 %>%
  filter(type != 0) %>%
  group_by(event_id) %>%
  summarise(
    duration_hours_median = median(crossing_duration_hours, na.rm = TRUE),
    .groups = "drop"
  )

# keep only focal event rows and attach median reference-step interval
focal_event_s3 = all_event_s3 %>%
  filter(type == 0) %>%
  left_join(ref_duration, by = "event_id")

# classify into GPS fix interval bins
focal_event_s3 = focal_event_s3 %>%
  mutate(
    fix_interval_cat = case_when(
      duration_hours_median <= 1 ~ "<=1",
      duration_hours_median > 1  & duration_hours_median <= 3  ~ "1-3",
      duration_hours_median > 3  & duration_hours_median <= 6  ~ "3-6",
      duration_hours_median > 6  & duration_hours_median <= 12 ~ "6-12",
      duration_hours_median > 12 & duration_hours_median <= 24 ~ "12-24",
      duration_hours_median > 24 ~ ">=24",
      TRUE ~ NA_character_
    )
  )

# preserve the desired row order
focal_event_s3$fix_interval_cat = factor(
  focal_event_s3$fix_interval_cat,
  levels = c("<=1", "1-3", "3-6", "6-12", "12-24", ">=24")
)

# count events by fix interval category and event type
summary_counts = focal_event_s3 %>%
  count(fix_interval_cat, location) %>%
  pivot_wider(
    names_from = location,
    values_from = n,
    values_fill = 0
  )

# build final table
table_s3_1 = summary_counts %>%
  mutate(
    Total_Event = `Known Crossing Event` +
      `Known Circumnavigating Event` +
      `Unknown Crossing Event`,
    Unknown_Event_Percent = round(
      100 * `Unknown Crossing Event` / Total_Event,
      1
    )
  ) %>%
  rename(
    Crossing_Event = `Known Crossing Event`,
    Circumnavigating_Event = `Known Circumnavigating Event`,
    Unknown_Event = `Unknown Crossing Event`
  )

# optional: add total row
table_s3_1_total = table_s3_1 %>%
  summarise(
    fix_interval_cat = "Total Number",
    Crossing_Event = sum(Crossing_Event, na.rm = TRUE),
    Circumnavigating_Event = sum(Circumnavigating_Event, na.rm = TRUE),
    Unknown_Event = sum(Unknown_Event, na.rm = TRUE),
    Total_Event = sum(Total_Event, na.rm = TRUE),
    Unknown_Event_Percent = NA_real_
  )

table_s3_1_final = bind_rows(table_s3_1, table_s3_1_total)

print(table_s3_1_final)


# ========================================
# Appendix S4
# Timing of albedo threshold transitions
# ========================================

library(raster)
library(ggplot2)
library(dplyr)
library(tidyr)
library(R.matlab)
library(RColorBrewer)

year = 2019
data_filename = paste0("output_", year, ".mat")
climatology_data = readMat(data_filename)
arr_climatology = array(unlist(climatology_data), dim = c(426, 724, 365))

arr_climatology_brick = brick(
  arr_climatology,
  xmn = -111.7125,
  xmx = -108.6958,
  ymn = 64.375,
  ymx = 66.15,
  crs = "+proj=longlat +datum=WGS84",
  transpose = FALSE
)

arr_climatology_brick_crop = crop(arr_climatology_brick, extent(polys1))
albedo_data = mask(arr_climatology_brick_crop, polys1)

# APR from DOY 92 to 280
percentile_rank = calc(albedo_data[[92:280]], fun = function(x) {
  rank(x, na.last = "keep") / sum(!is.na(x))
})

percentiles = c(0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.05)
pal_apr = rev(brewer.pal(11, "RdYlBu"))


# ----------------------------------------
# Appendix S4  Figure S1
# First spring date APR declined below percentile thresholds
# DOY 92-182
# ----------------------------------------

calc_first_decline_date = function(stack, percentile, doy_start) {
  calc(stack, fun = function(x) {
    idx = which(x < percentile)
    if (length(idx) == 0) return(NA)
    return(doy_start + idx[1] - 1)
  })
}

spring_stack = percentile_rank[[1:(182 - 92 + 1)]]

spring_rasters = lapply(percentiles, function(p) {
  calc_first_decline_date(
    stack = spring_stack,
    percentile = p,
    doy_start = 92
  )
})

spring_dfs = lapply(seq_along(spring_rasters), function(i) {
  df = as.data.frame(spring_rasters[[i]], xy = TRUE)
  colnames(df)[3] = "Date"
  df$Percentile = paste0(percentiles[i] * 100, "th")
  df
})

spring_df = bind_rows(spring_dfs) %>%
  drop_na(Date)

spring_df$Percentile = factor(
  spring_df$Percentile,
  levels = paste0(percentiles * 100, "th")
)

p_s4_figS1 = ggplot(spring_df, aes(x = x, y = y, fill = Date)) +
  geom_raster() +
  scale_fill_gradientn(
    colours = pal_apr,
    limits = c(92, 182),
    na.value = "transparent"
  ) +
  facet_wrap(~ Percentile, ncol = 2) +
  theme_classic() +
  labs(
    x = "Longitude",
    y = "Latitude",
    fill = "Day of Year"
  ) +
  theme(legend.position = "top")

print(p_s4_figS1)

# ----------------------------------------
# Appendix S4 Figure S2
# First fall date after seasonal minimum when APR exceeded
# percentile thresholds
# DOY 183-280
# ----------------------------------------

# identify lake-wide minimum APR date in 2019
layer_mean_values = sapply(1:nlayers(percentile_rank), function(i) {
  mean(values(percentile_rank[[i]]), na.rm = TRUE)
})

min_layer_index = which.min(layer_mean_values)   # layer index within 92-280
min_doy = 92 + min_layer_index - 1

calc_first_fall_recovery_date = function(stack, percentile, doy_start, min_layer_index) {
  calc(stack, fun = function(x) {
    threshold = quantile(x, probs = percentile, na.rm = TRUE)
    
    idx = which(x > threshold)
    
    if (length(idx) == 0) return(NA)
    
    idx_after_min = idx[idx > min_layer_index]
    
    if (length(idx_after_min) == 0) return(NA)
    
    return(doy_start + idx_after_min[1] - 1)
  })
}

fall_rasters = lapply(percentiles, function(p) {
  calc_first_fall_recovery_date(
    stack = percentile_rank,
    percentile = p,
    doy_start = 92,
    min_layer_index = min_layer_index
  )
})

fall_dfs = lapply(seq_along(fall_rasters), function(i) {
  df = as.data.frame(fall_rasters[[i]], xy = TRUE)
  colnames(df)[3] = "Date"
  df$Percentile = paste0(percentiles[i] * 100, "th")
  df
})

fall_df = bind_rows(fall_dfs) %>%
  drop_na(Date) %>%
  filter(Date >= 183)

fall_df$Percentile = factor(
  fall_df$Percentile,
  levels = paste0(rev(percentiles * 100), "th")
)

p_s4_figS2 = ggplot(fall_df, aes(x = x, y = y, fill = Date)) +
  geom_raster() +
  scale_fill_gradientn(
    colours = pal_apr,
    limits = c(183, 280),
    na.value = "transparent"
  ) +
  facet_wrap(~ Percentile, ncol = 2) +
  theme_classic() +
  labs(
    x = "Longitude",
    y = "Latitude",
    fill = "Day of Year"
  ) +
  theme(legend.position = "top")

print(p_s4_figS2)

# ----------------------------------------
# Appendix S4 Figure S3
# First date APR declined below percentile thresholds
# across the full observation period (DOY 92-280)
# ----------------------------------------

whole_rasters = lapply(percentiles, function(p) {
  calc_first_decline_date(
    stack = percentile_rank,
    percentile = p,
    doy_start = 92
  )
})

whole_dfs = lapply(seq_along(whole_rasters), function(i) {
  df = as.data.frame(whole_rasters[[i]], xy = TRUE)
  colnames(df)[3] = "Date"
  df$Percentile = paste0(percentiles[i] * 100, "th")
  df
})

whole_df = bind_rows(whole_dfs) %>%
  drop_na(Date)

whole_df$Percentile = factor(
  whole_df$Percentile,
  levels = paste0(percentiles * 100, "th")
)

p_s4_figS3 = ggplot(whole_df, aes(x = x, y = y, fill = Date)) +
  geom_raster() +
  scale_fill_gradientn(
    colours = pal_apr,
    limits = c(92, 280),
    na.value = "transparent"
  ) +
  facet_wrap(~ Percentile, ncol = 2) +
  theme_classic() +
  labs(
    x = "Longitude",
    y = "Latitude",
    fill = "Day of Year"
  ) +
  theme(legend.position = "top")

print(p_s4_figS3)

# ----------------------------------------
# Appendix S4 Figure S4
# Spatial pattern of known crossing events
# ----------------------------------------

library(ggplot2)
library(sf)
library(ggpubr)
library(dplyr)
library(RColorBrewer)

# =========================
# 1. Load data
# =========================
# Replace with your local file/object as needed
# on_lake <- read.csv("on_lake.csv", stringsAsFactors = FALSE)

on_lake <- as.data.frame(on_lake)

# on_lake <- on_lake %>%
#   filter(location == "Known Crossing Event")

# clean variables
on_lake <- on_lake %>%
  mutate(
    Season = as.factor(Season),
    Year = as.numeric(Year),
    Crossing_date = as.numeric(Crossing_date),
    Albedo = as.numeric(Albedo_percentile_rank),
    Location = as.factor(Location)
  )

# lake polygon should already be loaded as an sf object
# lake <- st_read("Contwoyto_pure_water.shp", quiet = TRUE) %>% st_transform(4326)

# =========================
# 2. Panel A
# Season
# =========================
p_s4_A <- ggplot(lake) +
  geom_sf(fill = "grey", alpha = 0.3) +
  geom_point(data = on_lake, aes(x = Lon, y = Lat, color = Season)) +
  scale_colour_manual(values = c("#4575B4", "#D73027")) +
  theme_minimal() +
  theme(legend.position = "top")

# =========================
# 3. Panel B
# Year
# =========================
p_s4_B <- ggplot(lake) +
  geom_sf(fill = "grey", alpha = 0.5) +
  geom_point(data = on_lake, aes(x = Lon, y = Lat, color = Year)) +
  scale_colour_gradientn(colours = brewer.pal(n = 10, name = "RdYlBu")) +
  theme_minimal() +
  theme(legend.position = "top")

# =========================
# 4. Panel C
# Crossing date
# =========================
p_s4_C <- ggplot(lake) +
  geom_sf(fill = "grey", alpha = 0.5) +
  geom_point(data = on_lake, aes(x = Lon, y = Lat, color = Crossing_date)) +
  scale_colour_gradientn(colours = brewer.pal(n = 11, name = "RdYlBu")) +
  theme_minimal() +
  theme(legend.position = "top") +
  labs(color = "Crossing Date (Day of Year)")

# =========================
# 5. Panel D
# Albedo percentile rank
# =========================
p_s4_D <- ggplot(lake) +
  geom_sf(fill = "grey", alpha = 0.5) +
  geom_point(data = on_lake, aes(x = Lon, y = Lat, color = Albedo)) +
  scale_colour_gradientn(colours = brewer.pal(n = 11, name = "RdYlBu")) +
  theme_minimal() +
  theme(legend.position = "top") +
  labs(color = "Albedo Percentile Rank")

# =========================
# 6. Panel E
# Crossing location type
# =========================
p_s4_E <- ggplot(lake) +
  geom_sf(fill = "grey", alpha = 0.5) +
  geom_point(data = on_lake, aes(x = Lon, y = Lat, color = Location)) +
  scale_colour_manual(values = c("#D73027", "#4575B4")) +
  theme_minimal() +
  theme(legend.position = "top")

# =========================
# 7. Combine panels
# =========================
fig_s4_combined <- ggarrange(
  p_s4_A, p_s4_B, p_s4_C,
  p_s4_D, p_s4_E,
  labels = c("A", "B", "C", "D", "E"),
  ncol = 3,
  nrow = 2
)

print(fig_s4_combined)


# ========================================
# Appendix S5 Figure S1
# Model fitting results for fall-only and combined datasets
#
# Panels:
#   A. RF variable importance, fall-only dataset
#   B. BL fitted relationship, fall-only dataset
#   C. RF variable importance, combined dataset
#   D. BL fitted relationship, combined dataset
# ========================================
library(ggplot2)
library(dplyr)
library(randomForest)
library(RColorBrewer)
library(ggpubr)

# ------------------------
# Panel A
# RF variable importance, fall-only dataset
# Source: 09_fall_models.R
# ------------------------
importance_fall <- as.data.frame(randomForest::importance(rf_fall))
importance_fall$Variable <- c(
  "Proportion of Time Spent in Lake Area",
  "APR at the Nearest Pixel",
  "APR Along the Potential Crossing Path",
  "APR of the Entire Lake Area",
  "Lake Width",
  "The Ratio of Circumnavigate Speeds Over Reference Speeds (log)",
  "The Ratio of Direct Speeds Over Reference Speeds (log)",
  "The Ratio of Circumnavigate Speeds Over Direct Speeds (log)",
  "Sex",
  "Migration Season"
)

p_s5_s1_A <- ggplot(
  importance_fall,
  aes(x = reorder(Variable, MeanDecreaseGini), y = MeanDecreaseGini)
) +
  geom_point(size = 3, color = "blue") +
  labs(
    x = "Variables",
    y = "Importance (Gini Decrease)"
  ) +
  coord_flip() +
  theme_classic()

print(p_s5_s1_A)

# ------------------------
# Panel B
# BL fitted relationship, fall-only dataset
# Source: 09_fall_models.R
# ------------------------
threshold_data_fall <- data.frame(
  log_ratio_cir_ref = RF_fall$log_ratio_cir_ref,
  predicted_probs = predicted_probs,
  location = case_when(
    RF_fall$response == 1 ~ "Known Crossing Event",
    RF_fall$response == 0 ~ "Known Circumnavigating Event"
  )
)

log_ratio_cir_ref_grid <- seq(
  min(RF_fall$log_ratio_cir_ref),
  max(RF_fall$log_ratio_cir_ref),
  length = 100
)

new_data_ratio_fall <- data.frame(
  log_ratio_cir_ref = log_ratio_cir_ref_grid
)

predicted_probs_with_se_fall <- predict(
  logistic_model_fall,
  newdata = new_data_ratio_fall,
  type = "response",
  se.fit = TRUE
)

new_data_ratio_fall$lower_bound <- pmax(
  0,
  predicted_probs_with_se_fall$fit - 1.96 * predicted_probs_with_se_fall$se.fit
)
new_data_ratio_fall$upper_bound <- pmin(
  1,
  predicted_probs_with_se_fall$fit + 1.96 * predicted_probs_with_se_fall$se.fit
)
new_data_ratio_fall$predicted_probs <- predicted_probs_with_se_fall$fit

p_s5_s1_B <- ggplot(
  new_data_ratio_fall,
  aes(x = log_ratio_cir_ref, y = predicted_probs)
) +
  geom_line(color = "blue") +
  geom_ribbon(
    aes(ymin = lower_bound, ymax = upper_bound),
    alpha = 0.2
  ) +
  geom_point(
    data = threshold_data_fall,
    aes(x = log_ratio_cir_ref, y = predicted_probs, color = location),
    size = 2
  ) +
  labs(
    x = "The Ratio of Circumnavigate Speeds Over Reference Speeds (log)",
    y = "Predicted Probability of Crossing",
    color = "Event"
  ) +
  theme_classic()

print(p_s5_s1_B)


# ------------------------
# Panel C
# RF variable importance, combined dataset
# Source: 07_combined_models.R
# ------------------------
importance_combined <- as.data.frame(randomForest::importance(rf))
importance_combined$Variable <- c(
  "Proportion of Time Spent in Lake Area",
  "APR at the Nearest Pixel",
  "APR Along the Potential Crossing Path",
  "APR of the Entire Lake Area",
  "Lake Width",
  "The Ratio of Circumnavigate Speeds Over Reference Speeds (log)",
  "The Ratio of Direct Speeds Over Reference Speeds (log)",
  "The Ratio of Circumnavigate Speeds Over Direct Speeds (log)",
  "Sex",
  "Migration Season"
)

p_s5_s1_C <- ggplot(
  importance_combined,
  aes(x = reorder(Variable, MeanDecreaseGini), y = MeanDecreaseGini)
) +
  geom_point(size = 3, color = "blue") +
  labs(
    x = "Variables",
    y = "Importance (Gini Decrease)"
  ) +
  coord_flip() +
  theme_classic()

print(p_s5_s1_C)

# ------------------------
# Panel D
# BL fitted relationship, combined dataset
# Source: 07_combined_models.R
# ------------------------
threshold_data_combined <- data.frame(
  albedo_linear = RF_data$albedo_linear,
  lake_spent_proportion = RF_data$lake_spent_proportion,
  log_ratio_cir_ref = RF_data$log_ratio_cir_ref,
  predicted_probs = predicted_probs,
  location = case_when(
    RF_data$response == 1 ~ "Known Crossing Event",
    RF_data$response == 0 ~ "Known Circumnavigating Event"
  )
)

albedo_grid <- seq(
  min(RF_data$albedo_linear),
  max(RF_data$albedo_linear),
  length = 100
)

new_data_albedo <- data.frame(
  albedo_linear = albedo_grid,
  lake_spent_proportion = median(RF_data$lake_spent_proportion),
  log_ratio_cir_ref = median(RF_data$log_ratio_cir_ref)
)

pred_albedo <- predict(
  model_all_logistic,
  newdata = new_data_albedo,
  type = "response",
  se.fit = TRUE
)

new_data_albedo$lower_bound <- pmax(
  0,
  pred_albedo$fit - 1.96 * pred_albedo$se.fit
)
new_data_albedo$upper_bound <- pmin(
  1,
  pred_albedo$fit + 1.96 * pred_albedo$se.fit
)
new_data_albedo$predicted_probs <- pred_albedo$fit

p_s5_s1_D <- ggplot(
  new_data_albedo,
  aes(x = albedo_linear, y = predicted_probs)
) +
  geom_line(color = "blue") +
  geom_ribbon(
    aes(ymin = lower_bound, ymax = upper_bound),
    alpha = 0.2
  ) +
  geom_point(
    data = threshold_data_combined,
    aes(x = albedo_linear, y = predicted_probs, color = location),
    size = 1.5
  ) +
  labs(
    x = "APR Along the Potential Crossing Path",
    y = "Predicted Probability of Crossing",
    color = "Event"
  ) +
  theme_classic()

print(p_s5_s1_D)


# ------------------------
# Combine panels
# ------------------------
fig_s5_s1 <- ggarrange(
  p_s5_s1_A, p_s5_s1_B,
  p_s5_s1_C, p_s5_s1_D,
  labels = c("A", "B", "C", "D"),
  ncol = 2,
  nrow = 2
)

print(fig_s5_s1)



# ----------------------------------------
# Appendix S5 Figure S2
# Predicted crossing probabilities for unknown fall events
# from BL and RF models
#
# Source: 09_fall_models.R
# ----------------------------------------

library(dplyr)
library(ggplot2)

unknown_data_fall <- Unknown_event_fall %>%
  mutate(
    sex = as.factor(sex),
    season = as.factor(season)
  ) %>%
  dplyr::select(
    lake_spent_proportion,
    albedo_nearest_pixel,
    albedo_linear,
    albedo_whole_lake,
    lake_width,
    log_ratio_cir_ref,
    log_ratio_dir_ref,
    log_ratio_cir_dir,
    sex,
    season
  ) %>%
  filter(rowSums(is.na(.)) == 0)

# RF prediction
unknown_predicted_probs_fall_rf <- predict(
  rf_fall,
  unknown_data_fall,
  type = "prob"
)
unknown_data_fall$rf_predicted_probs <- unknown_predicted_probs_fall_rf[, 2]

# BL prediction
unknown_predicted_probs_fall_logit <- predict(
  logistic_model_fall,
  newdata = unknown_data_fall,
  type = "response",
  se.fit = TRUE
)

unknown_data_fall$predicted_probs <- unknown_predicted_probs_fall_logit$fit
unknown_data_fall$lower_bound <- pmax(
  0,
  unknown_predicted_probs_fall_logit$fit - 1.96 * unknown_predicted_probs_fall_logit$se.fit
)
unknown_data_fall$upper_bound <- pmin(
  1,
  unknown_predicted_probs_fall_logit$fit + 1.96 * unknown_predicted_probs_fall_logit$se.fit
)

p_s5_s2 <- ggplot(
  new_data_ratio_fall,
  aes(x = log_ratio_cir_ref, y = predicted_probs)
) +
  geom_line(color = "blue", linewidth = 1) +
  geom_point(
    data = threshold_data_fall,
    aes(x = log_ratio_cir_ref, y = predicted_probs, color = location),
    size = 2
  ) +
  geom_point(
    data = unknown_data_fall,
    aes(
      x = log_ratio_cir_ref,
      y = predicted_probs,
      color = "Unknown Event (Binomial Logistic Model)"
    ),
    size = 2
  ) +
  geom_point(
    data = unknown_data_fall,
    aes(
      x = log_ratio_cir_ref,
      y = rf_predicted_probs,
      color = "Unknown Event (Random Forest Model)"
    ),
    size = 2
  ) +
  labs(
    x = "The Ratio of Circumnavigate Speeds Over Reference Speeds (log)",
    y = "Predicted Probability of Crossing",
    color = "Event"
  ) +
  scale_color_manual(
    values = c(
      "Known Crossing Event" = "#00BFC4",
      "Known Circumnavigating Event" = "#F8766D",
      "Unknown Event (Binomial Logistic Model)" = "orange",
      "Unknown Event (Random Forest Model)" = "#C77CFF"
    )
  ) +
  theme_classic()

print(p_s5_s2)


# ========================================
# Appendix S6
# Inter-annual 56th-threshold analysis across Contwoyto Lake
#
# Final definitions used in the appendix:
#   Value_56 = annual 56th-percentile albedo value for each lake pixel
#   DOY_56   = first spring day during the melt-down phase when
#              pixel-level APR first falls below 0.56
# ========================================

library(R.matlab)
library(terra)
library(sf)
library(dplyr)
library(ggplot2)
library(viridis)
library(RColorBrewer)
library(svglite)
library(scales)

## -----------------------------
# Load lake polygon
## -----------------------------
lake_shp <- sf::st_read("Contwoyto_pure_water.shp", quiet = TRUE)
lake_vect <- terra::vect(lake_shp)
terra::crs(lake_vect) <- "EPSG:4326"

dir.create("Output_Maps_Value", showWarnings = FALSE)
dir.create("Output_Maps_Date", showWarnings = FALSE)
dir.create("terra_tmp", showWarnings = FALSE)

terra::terraOptions(tempdir = normalizePath("terra_tmp"))

years <- 2005:2021

list_values_rasters <- vector("list", length(years))
list_dates_rasters  <- vector("list", length(years))
full_pixel_list     <- vector("list", length(years))

## -----------------------------
## Loop through years
## -----------------------------

for (i in seq_along(years)) {
  year <- years[i]
  message("Processing year: ", year)
  
  data_filename <- paste0("output_", year, ".mat")
  climatology_data <- readMat(data_filename)
  arr_climatology <- array(unlist(climatology_data), dim = c(426, 724, 365))
  
  r <- rast(arr_climatology)
  ext(r) <- c(-111.7125, -108.6958, 64.375, 66.15)
  crs(r) <- "+proj=longlat +datum=WGS84"
  
  # crop and mask
  r_crop <- crop(r, polys1)
  albedo_data <- mask(r_crop, polys1)
  
  # observation period: DOY 92-280
  albedo_subset <- albedo_data[[92:280]]
  
  ## Pixel-level percentile rank within each year and pixel
  percentile_rank <- app(albedo_subset, fun = function(x) {
    if (all(is.na(x))) return(rep(NA, length(x)))
    rank(x, na.last = "keep", ties.method = "first") / sum(!is.na(x))
  })
  
  # ------------------------
  # Panel A
  # Annual 56th-percentile albedo value for each pixel
  # ------------------------
  
  pixel_threshold_values <- app(albedo_subset, fun = function(x) {
    if (all(is.na(x))) return(NA)
    quantile(x, probs = 0.56, na.rm = TRUE, names = FALSE)
  })
  names(pixel_threshold_values) <- "Value_56"
  
  writeRaster(
    pixel_threshold_values,
    filename = paste0("Output_Maps_Value/Value_56th_", year, ".tif"),
    overwrite = TRUE
  )
  
  list_values_rasters[[i]] <- pixel_threshold_values
  
  # ------------------------
  # Panel B 
  # First spring DOY when APR falls below annual 56th percentile
  # Use melting phase only, defined as period up to annual minimum APR
  # ------------------------
  pixel_threshold_idx <- app(percentile_rank, fun = function(x) {
    
    if (all(is.na(x))) return(NA)
    
    idx_min <- max(which(x == min(x, na.rm = TRUE)))
    if (length(idx_min) == 0 || is.infinite(idx_min)) return(NA)
    
    x_melting_phase <- x[1:idx_min]
    below_threshold_indices <- which(x_melting_phase <= 0.56)
    
    if (length(below_threshold_indices) > 0) {
      return(below_threshold_indices[1])
    } else {
      idx_closest <- which.min(abs(x_melting_phase - 0.56))
      if (length(idx_closest) == 0) return(NA)
      return(idx_closest)
    }
  })
  
  pixel_threshold_doy <- pixel_threshold_idx + 91
  names(pixel_threshold_doy) <- "DOY_56"
  
  writeRaster(
    pixel_threshold_doy,
    filename = paste0("Output_Maps_Date/Date_DOY_56th_", year, ".tif"),
    overwrite = TRUE
  )
  
  list_dates_rasters[[i]] <- pixel_threshold_doy
  
  ## Extract pixel-level values for plotting
  stack_current <- c(pixel_threshold_values, pixel_threshold_doy)
  df_pixels <- as.data.frame(stack_current, na.rm = TRUE)
  df_pixels$Year <- year
  
  full_pixel_list[[i]] <- df_pixels
}

## -----------------------------
## Combine all pixel-level data
## -----------------------------
big_data_df <- bind_rows(full_pixel_list) %>%
  filter(DOY_56 > 92) %>%   
  mutate(
    Year = factor(Year, levels = years),
    Year_num = as.numeric(as.character(Year))
  )

write.csv(big_data_df, "Full_Pixel_Data_Violin_Ready.csv", row.names = FALSE)

save(
  big_data_df,
  list_dates_rasters,
  list_values_rasters,
  file = "interannual_56th_trend_DOY_Albedo.RData"
)

## -----------------------------
## Figure S1A: violin plot for annual 56th-percentile albedo
## -----------------------------
p_s1a <- ggplot(big_data_df, aes(x = Year, y = Value_56, fill = Year)) +
  geom_violin(trim = FALSE, scale = "width", alpha = 0.7, color = NA) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA, alpha = 0.8) +
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 23,
    size = 2,
    fill = "red",
    color = "black"
  ) +
  scale_fill_viridis_d() +
  labs(
    x = "Year",
    y = "Albedo"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

## -----------------------------
## Figure S1B: violin plot for DOY of threshold crossing
## -----------------------------
p_s1b <- ggplot(big_data_df, aes(x = Year, y = DOY_56, fill = Year)) +
  geom_violin(trim = TRUE, scale = "width", alpha = 0.7, color = NA) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA, alpha = 0.8) +
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 23,
    size = 2,
    fill = "red",
    color = "black"
  ) +
  stat_summary(
    fun = min,
    geom = "point",
    shape = 21,
    size = 1.8,
    fill = "white",
    color = "black"
  ) +
  scale_fill_viridis_d() +
  labs(
    x = "Year",
    y = "Day of Year (DOY)"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

## -----------------------------
## Prepare data for spatial maps
## -----------------------------
list_dates_rasters<- lapply(list_dates_rasters, function(r) {
  r[r <= 92] <- NA
  r
})

df_doy <- bind_rows(lapply(seq_along(list_dates_rasters), function(i) {
  r <- list_dates_rasters[[i]]
  names(r) <- "DOY_56"
  
  as.data.frame(r, xy = TRUE, na.rm = TRUE) %>%
    filter(DOY_56 > 92) %>%
    mutate(Year = years[i])
}))

df_value <- bind_rows(lapply(seq_along(list_values_rasters), function(i) {
  r <- list_values_rasters[[i]]
  names(r) <- "Value_56"
  
  as.data.frame(r, xy = TRUE, na.rm = TRUE) %>%
    mutate(Year = years[i])
}))

## -----------------------------
## Figure S2: spatial maps of annual 56th-percentile albedo
## -----------------------------
p_s2 <- ggplot(df_value, aes(x = x, y = y, fill = Value_56)) +
  geom_raster() +
  coord_equal(expand = FALSE) +
  facet_wrap(~ Year, ncol = 3) +
  scale_fill_gradientn(
    # colours = rev(brewer.pal(n = 11, name = "YlOrRd")),
    colours = rev(brewer.pal(n = 11, name = "RdYlBu")),
    oob = scales::squish,
    name = "Albedo at APR = 0.56"
  ) +
  scale_x_continuous(
    breaks = seq(
      floor(min(df_value$x, na.rm = TRUE)),
      ceiling(max(df_value$x, na.rm = TRUE)),
      by = 1
    )
  )+
  theme_classic() +
  labs(
    x = "Longitude",
    y = "Latitude",
    fill = "Albedo at APR = 0.56")+
  theme(
    strip.text = element_text(size = 9),
    legend.position = "top"
  )

## -----------------------------
## Figure S3: spatial maps of DOY threshold crossing
## -----------------------------
p_s3 <- ggplot(df_doy, aes(x = x, y = y, fill = DOY_56)) +
  geom_raster() +
  coord_equal(expand = FALSE) +
  facet_wrap(~ Year, ncol = 3) +
  scale_fill_gradientn(
    # colours = rev(brewer.pal(n = 11, name = "YlOrRd")),
    colours = rev(brewer.pal(n = 11, name = "RdYlBu")),
    oob = scales::squish,
    name = "DOY at APR = 0.56"
  ) +
  scale_x_continuous(
    breaks = seq(
      floor(min(df_doy$x, na.rm = TRUE)),
      ceiling(max(df_doy$x, na.rm = TRUE)),
      by = 1
    )
  )+
  theme_classic() +
  labs(
    x = "Longitude",
    y = "Latitude",
    fill = "DOY at APR = 0.56")+
  theme(
    strip.text = element_text(size = 9),
    legend.position = "top"
  )

print(p_s3)

