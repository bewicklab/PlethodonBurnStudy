
# ==========================================================
# 1. REQUIRED PACKAGES
# ==========================================================

library(sf)
library(dplyr)
library(tmap)
library(grid)


# ==========================================================
# 2. LOAD AND CLEAN METADATA
# ==========================================================

meta <- read.csv(
  "Coweeta_TEWL_Seasons_Sites.csv",
  stringsAsFactors = FALSE
)

meta <- meta %>%
  filter(
    !is.na(QIIME_ID_aligned),
    QIIME_ID_aligned != ""
  ) %>%
  distinct(QIIME_ID_aligned, .keep_all = TRUE)

meta$ID <- tolower(
  gsub("-", ".", meta$QIIME_ID_aligned)
)


# ==========================================================
# 3. REMOVE UNWANTED SAMPLES
# ==========================================================

meta <- meta %>%
  filter(
    !Site %in% c(
      "Blank",
      "Control",
      "2750",
      "2850",
      "EF",
      "PN",
      "Positive_control"
    ),
    
    !title %in% c(
      "TS.31.5",
      "TS.31.6",
      "TS.31.7",
      "TS.31.8",
      "TS.31.9",
      "TS.31.10",
      "TS.31.11"
    ),
    
    ID != "c1.25.27.noburn",
    ID != "ts.32.1",
    
    !Season %in% c(
      "Summer_2024",
      "Summer_2025",
      "Fall_2024"
    )
  )


# ==========================================================
# 4. KEEP THREE STUDY SITES
# ==========================================================

meta <- meta %>%
  filter(
    Site %in% c(
      "31_burn",
      "32_no_burn",
      "2950"
    )
  )


# ==========================================================
# 5. PREPARE POINT DATA
# ==========================================================

points <- meta %>%
  transmute(
    ID = ID,
    
    longitude = as.numeric(
      long_7_Collection_locatio
    ),
    
    latitude = as.numeric(
      lat_7_Collection_locatio
    ),
    
    site = recode(
      Site,
      "31_burn" = "Burn",
      "32_no_burn" = "Control_2",
      "2950" = "Control_1"
    ),
    
    year = sub(
      ".*_",
      "",
      Season
    )
  ) %>%
  filter(
    !is.na(longitude),
    !is.na(latitude),
    is.finite(longitude),
    is.finite(latitude),
    between(longitude, -180, 180),
    between(latitude, -90, 90),
    !is.na(site),
    !is.na(year),
    year %in% c("2024", "2025")
  )


# ==========================================================
# 6. CHECK DATA
# ==========================================================

cat("\nNumber of samples:", nrow(points), "\n")

cat("\nSamples per site:\n")
print(table(points$site))

cat("\nSamples per year:\n")
print(table(points$year))

if (nrow(points) == 0) {
  stop("No valid mapping points remain.")
}

sites <- c(
  "Control_1",
  "Control_2",
  "Burn"
)

if (!all(sites %in% points$site)) {
  stop("One or more study sites have no valid mapping points.")
}


# ==========================================================
# 7. CONVERT TO SPATIAL DATA
# ==========================================================

pts <- st_as_sf(
  points,
  coords = c("longitude", "latitude"),
  crs = 4326,
  remove = FALSE
)

pts$site <- factor(
  pts$site,
  levels = sites
)

pts$year <- factor(
  pts$year,
  levels = c("2024", "2025")
)

stopifnot(!anyNA(st_coordinates(pts)))
stopifnot(!anyNA(pts$site))
stopifnot(!anyNA(pts$year))


# ==========================================================
# 8. SITE COLORS
# ==========================================================

my_colors <- c(
  "Control_1" = "#228A22",
  "Control_2" = "#00FE00",
  "Burn" = "red"
)

tmap_mode("plot")


# ==========================================================
# 9. OVERALL STUDY AREA
# ==========================================================

area <- st_union(pts) %>%
  st_transform(3857) %>%
  st_buffer(100) %>%
  st_transform(4326)


# ==========================================================
# 10. INDIVIDUAL SITE BUFFERS
# ==========================================================

site_areas <- pts %>%
  group_by(site) %>%
  summarise(
    geometry = st_union(geometry),
    .groups = "drop"
  ) %>%
  st_transform(3857) %>%
  st_buffer(30) %>%
  st_transform(4326)


# ==========================================================
# 11. CALCULATE SITE DIMENSIONS
# ==========================================================

pts_3857 <- st_transform(pts, 3857)

site_info <- lapply(sites, function(s) {
  
  p <- pts_3857[pts_3857$site == s, ]
  
  if (nrow(p) == 0) {
    stop(paste("No points found for site:", s))
  }
  
  bb <- st_bbox(p)
  
  width <- unname(bb["xmax"] - bb["xmin"]) + 200
  
  height <- unname(bb["ymax"] - bb["ymin"]) + 200
  
  data.frame(
    site = s,
    width = width,
    height = height,
    aspect = width / height
  )
  
}) %>%
  bind_rows()


# ==========================================================
# 12. CREATE INDIVIDUAL SITE MAPS
# ==========================================================

maps <- lapply(sites, function(s) {
  
  p <- pts[pts$site == s, ]
  
  area_site <- site_areas[
    site_areas$site == s,
  ]
  
  # Site center
  
  p_3857 <- st_transform(p, 3857)
  
  bb <- st_bbox(p_3857)
  
  x_center <- (
    unname(bb["xmin"]) +
      unname(bb["xmax"])
  ) / 2
  
  y_center <- (
    unname(bb["ymin"]) +
      unname(bb["ymax"])
  ) / 2
  
  
  # Site dimensions
  
  info <- site_info[
    site_info$site == s,
  ]
  
  w <- info$width
  h <- info$height
  
  
  # Geographic bounding box
  
  site_bbox_3857 <- st_bbox(
    c(
      xmin = x_center - w / 2,
      xmax = x_center + w / 2,
      ymin = y_center - h / 2,
      ymax = y_center + h / 2
    ),
    crs = st_crs(3857)
  )
  
  site_bbox <- st_bbox(
    st_transform(
      st_as_sfc(site_bbox_3857),
      4326
    )
  )
  
  
  # Map
  
  tm_shape(
    area_site,
    bbox = site_bbox
  ) +
    
    tm_basemap("OpenStreetMap") +
    
    tm_shape(p) +
    
    tm_dots(
      fill = "site",
      
      fill.scale = tm_scale_categorical(
        values = my_colors
      ),
      
      shape = "year",
      
      shape.scale = tm_scale_categorical(
        values = c(21, 24)
      ),
      
      size = 0.85
    ) +
    
    tm_legend_hide() +
    
    tm_layout(
      frame = TRUE,
      frame.lwd = 1,
      frame.col = "black",
      inner.margins = c(0, 0, 0, 0),
      outer.margins = c(0, 0, 0, 0)
    )
})


# ==========================================================
# 13. FULL THREE-SITE MAP
# ==========================================================

full_map <- tm_shape(area) +
  
  tm_basemap("OpenStreetMap") +
  
  tm_shape(pts) +
  
  tm_dots(
    fill = "site",
    
    fill.scale = tm_scale_categorical(
      values = my_colors
    ),
    
    shape = "year",
    
    shape.scale = tm_scale_categorical(
      values = c(21, 24)
    ),
    
    size = 1
  ) +
  
  tm_legend_hide()


# ==========================================================
# 14. OVERALL FIGURE LAYOUT
#
# LEFT: 40%
# RIGHT: 60%
# ==========================================================

graphics.off()

grid.newpage()

pushViewport(
  viewport(
    layout = grid.layout(
      1,
      2,
      widths = unit(
        c(0.4, 0.6),
        "null"
      )
    )
  )
)


# ==========================================================
# 15. FULL MAP ON LEFT
# ==========================================================

print(
  full_map,
  vp = viewport(
    layout.pos.row = 1,
    layout.pos.col = 1
  )
)


# ==========================================================
# 16. RIGHT-HAND MAP ORDER
# ==========================================================

plot_order <- c(
  "Burn",
  "Control_2",
  "Control_1"
)

site_info_right <- site_info[
  match(
    plot_order,
    site_info$site
  ),
]


# ==========================================================
# 17. CALCULATE COMMON MAP SCALE
# ==========================================================

right_width_mm <- convertWidth(
  unit(0.6, "npc"),
  "mm",
  valueOnly = TRUE
)

right_height_mm <- convertHeight(
  unit(1, "npc"),
  "mm",
  valueOnly = TRUE
)

max_width <- max(
  site_info_right$width
)

total_height <- sum(
  site_info_right$height
)

common_scale <- min(
  right_width_mm / max_width,
  right_height_mm / total_height
)


# ==========================================================
# 18. PHYSICAL MAP DIMENSIONS
# ==========================================================

site_info_right <- site_info_right %>%
  mutate(
    panel_width_mm = width * common_scale,
    panel_height_mm = height * common_scale
  )


# ==========================================================
# 19. VERTICAL SPACING
# ==========================================================

total_panel_height_mm <- sum(
  site_info_right$panel_height_mm
)

remaining_height_mm <- (
  right_height_mm -
    total_panel_height_mm
)

gap_mm <- remaining_height_mm / 7


# ==========================================================
# 20. HORIZONTAL POSITION
#
# SHIFT SECOND COLUMN 25 MM TO THE RIGHT
# ==========================================================

right_left_mm <- convertWidth(
  unit(0.4, "npc"),
  "mm",
  valueOnly = TRUE
)

burn_info <- site_info_right[
  site_info_right$site == "Burn",
]

burn_width_mm <- burn_info$panel_width_mm

burn_center_x_mm <- right_left_mm +
  right_width_mm / 9

# Horizontal adjustment
horizontal_shift_mm <- 25

common_left_edge_mm <- burn_center_x_mm -
  burn_width_mm / 2 +
  horizontal_shift_mm


# ==========================================================
# 21. DRAW RIGHT-HAND MAPS
# ==========================================================

current_top_mm <- right_height_mm - gap_mm

for (s in plot_order) {
  
  info <- site_info_right[
    site_info_right$site == s,
  ]
  
  panel_width_mm <- info$panel_width_mm
  
  panel_height_mm <- info$panel_height_mm
  
  panel_center_y_mm <- current_top_mm -
    panel_height_mm / 2
  
  # Adjust middle panel
  if (s == "Control_2") {
    
    panel_center_y_mm <-
      panel_center_y_mm + 2.75
    
  }
  
  # All maps share the same left edge
  panel_center_x_mm <- common_left_edge_mm +
    panel_width_mm / 2
  
  map_vp <- viewport(
    
    x = unit(
      panel_center_x_mm,
      "mm"
    ),
    
    y = unit(
      panel_center_y_mm,
      "mm"
    ),
    
    width = unit(
      panel_width_mm,
      "mm"
    ),
    
    height = unit(
      panel_height_mm,
      "mm"
    ),
    
    just = c(
      "center",
      "center"
    )
  )
  
  print(
    maps[[which(sites == s)]],
    vp = map_vp
  )
  
  current_top_mm <- current_top_mm -
    panel_height_mm -
    gap_mm
}

popViewport()


# ==========================================================
# 22. LEGEND POSITION
#
# LEGEND ALSO MOVES WITH SECOND COLUMN
# ==========================================================

control1_info <- site_info_right[
  site_info_right$site == "Control_1",
]

control1_width_mm <- control1_info$panel_width_mm

control1_height_mm <- control1_info$panel_height_mm

legend_x_mm <- common_left_edge_mm +
  control1_width_mm + 10

control1_center_y_mm <- current_top_mm +
  gap_mm +
  control1_height_mm / 2


# ==========================================================
# 23. DRAW LEGEND
# ==========================================================

pushViewport(
  viewport(
    
    x = unit(
      legend_x_mm,
      "mm"
    ),
    
    y = unit(
      control1_center_y_mm,
      "mm"
    ),
    
    width = unit(
      50,
      "mm"
    ),
    
    height = unit(
      50,
      "mm"
    ),
    
    just = c(
      "left",
      "center"
    )
  )
)


# LEGEND TITLE

grid.text(
  "Site / Year",
  x = unit(0, "mm"),
  y = unit(47, "mm"),
  just = c("left", "center"),
  gp = gpar(
    fontsize = 12,
    fontface = "bold"
  )
)


# 2024

grid.points(
  x = unit(5, "mm"),
  y = unit(37, "mm"),
  pch = 21,
  size = unit(5, "mm"),
  gp = gpar(
    col = "black",
    fill = "white",
    lwd = 1
  )
)

grid.text(
  "2024",
  x = unit(12, "mm"),
  y = unit(37, "mm"),
  just = "left",
  gp = gpar(fontsize = 11)
)


# 2025

grid.points(
  x = unit(5, "mm"),
  y = unit(27, "mm"),
  pch = 24,
  size = unit(5, "mm"),
  gp = gpar(
    col = "black",
    fill = "white",
    lwd = 1
  )
)

grid.text(
  "2025",
  x = unit(12, "mm"),
  y = unit(27, "mm"),
  just = "left",
  gp = gpar(fontsize = 11)
)


# SITE COLORS

legend_sites <- c(
  "Control_1",
  "Control_2",
  "Burn"
)

legend_labels <- c(
  "Control 1",
  "Control 2",
  "Burn"
)

legend_y <- c(
  17,
  9,
  1
)

for (i in seq_along(legend_sites)) {
  
  grid.points(
    x = unit(5, "mm"),
    y = unit(legend_y[i], "mm"),
    pch = 21,
    size = unit(5, "mm"),
    gp = gpar(
      col = "black",
      fill = my_colors[legend_sites[i]],
      lwd = 1
    )
  )
  
  grid.text(
    legend_labels[i],
    x = unit(12, "mm"),
    y = unit(legend_y[i], "mm"),
    just = "left",
    gp = gpar(fontsize = 11)
  )
}

popViewport()


# ==========================================================
# END
# ==========================================================
