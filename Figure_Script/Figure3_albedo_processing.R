# ----------------------------------------
# Figure3_albedo_processing.R
# Purpose:
#   Plot albedo preprocessing results and the derived
#   Albedo Percentile Rank (APR) for an example year (2019).
#
# Panels:
#   A. Albedo of Contwoyto Lake on DOY 92 (April 2, 2019)
#   B. Lake-wide average APR from DOY 92 to DOY 280
#   C. Spatial distribution of APR on 10 selected dates
#
# Upstream logic:
#   - Annual albedo climatology outputs were generated previously
#   - Lake polygon was defined from the Contwoyto pure-water shapefile
#
# Required input files:
#   - output_2019.mat
#   - Contwoyto_pure_water.shp
# ----------------------------------------

library(terra)
library(raster)
library(R.matlab)
library(sf)
library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)

# =========================
# 1. Load lake polygon
# =========================
polys1_sf = st_read("Contwoyto_pure_water.shp", quiet = TRUE)
polys1_sf = st_transform(polys1_sf, 4326)
polys1 = as(polys1_sf, "Spatial")

# =========================
# 2. Load annual albedo data
# =========================
target_year = 2019
data_filename = paste0("output_", target_year, ".mat")

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

# =========================
# 3. Calculate APR for DOY 92-280
# =========================
percentile_rank = calc(albedo_data[[92:280]], fun = function(x) {
  percentile = rank(x, na.last = "keep") / sum(!is.na(x))
  return(percentile)
})

# common palette
color_palette = rev(brewer.pal(11, "RdYlBu"))

# =========================
# 4. Figure 3A
# Albedo on DOY 92 (April 2, 2019)
# =========================
panelA_df = as.data.frame(albedo_data[[92]], xy = TRUE, na.rm = TRUE)
colnames(panelA_df) = c("x", "y", "Albedo")

p_fig3A = ggplot(panelA_df, aes(x = x, y = y, fill = Albedo)) +
  geom_raster() +
  scale_fill_gradientn(
    colours = color_palette,
    limits = c(0, 1),
    na.value = "transparent"
  ) +
  theme_classic() +
  labs(
    x = "Longitude",
    y = "Latitude",
    fill = "Albedo"
  )

print(p_fig3A)

# =========================
# 5. Figure 3B
# Lake-wide average APR from DOY 92 to DOY 280
# =========================
layer_mean_values = sapply(1:nlayers(percentile_rank), function(i) {
  mean(values(percentile_rank[[i]]), na.rm = TRUE)
})

mean_values_df = data.frame(
  Day = 92:280,
  MeanValue = layer_mean_values
)

p_fig3B = ggplot(mean_values_df, aes(x = Day, y = MeanValue)) +
  geom_line(color = "black", size = 0.5) +
  geom_point(color = "darkgrey", size = 0.8) +
  theme_classic() +
  labs(
    x = "Day of the Year",
    y = "Lake-wide Average APR"
  )

print(p_fig3B)

# =========================
# 6. Figure 3C
# Spatial APR on 10 selected dates between DOY 92 and 280
# =========================
start_date = as.Date("2019-04-02")
end_date = as.Date("2019-10-07")
dates = seq.Date(start_date, end_date, by = "day")

evenly_spaced_dates = seq.Date(start_date, end_date, length.out = 10)

layer_indices = sapply(evenly_spaced_dates, function(date_i) {
  which.min(abs(as.numeric(dates - date_i)))
})

raster_list = list()

for (i in seq_along(layer_indices)) {
  layer_index = layer_indices[i]
  date_label = format(evenly_spaced_dates[i], "%Y-%m-%d")
  
  raster_layer = percentile_rank[[layer_index]]
  raster_df = as.data.frame(raster_layer, xy = TRUE)
  colnames(raster_df)[3] = "APR"
  raster_df$Date = date_label
  
  raster_list[[i]] = raster_df
}

combined_raster_df = bind_rows(raster_list) %>%
  drop_na(APR)

p_fig3C = ggplot(combined_raster_df, aes(x = x, y = y, fill = APR)) +
  geom_raster() +
  scale_fill_gradientn(
    colours = color_palette,
    na.value = "transparent"
  ) +
  facet_wrap(~ Date, ncol = 5) +
  theme_classic2() +
  labs(
    x = "Longitude",
    y = "Latitude",
    fill = "Albedo in Percentile Rank"
  ) +
  theme(legend.position = "top")

print(p_fig3C)

# =========================
# 7. Optional: combine panels for review
# =========================
fig3_combined = ggarrange(
  p_fig3A,
  p_fig3B,
  p_fig3C,
  labels = c("A", "B", "C"),
  ncol = 1,
  nrow = 3
)

print(fig3_combined)

# Optional save
# ggsave("Figure3_albedo_processing.png", fig3_combined, width = 10, height = 14, dpi = 300)