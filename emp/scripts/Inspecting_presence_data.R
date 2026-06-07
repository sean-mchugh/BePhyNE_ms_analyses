## ============================================================
## Coordinate sanity check for GBIF occurrence data
##
## Does:
##   1. Loads GBIF presence object
##   2. Handles list entries that are data.frames, matrices, or named vectors
##   3. Runs CoordinateCleaner on all occurrence coordinates
##   4. Summarizes failures overall and by species
##   5. Plots sea-flagged points relative to US + Canada
##   6. Saves CSV summaries and PNG figures
##
## ============================================================


# So we have a large set of occurrence data filtered down 
# and confirmed by eye they fit expert range maps (+/- a modest buffer) 
# we also need to confirm our occurrence arent violating any major buffers using coordinate cleaner
# this isnt as straight forward as throwing out any occurrence that fails a check, 
# many species are coastal and may erroneously violate a sea check (provided they totall off the coast), 
# some species are so widespread I wouldnt be surprised if they were found in a county centroid
# some species are so widespread and successful (ie Red backs, or duskies) I wouldnt be surprised if they were found a country capital
# the goal here is to make sure there isnt an alarming outsize proportion of these points and confirm for the species that do fail a check, it makes sense it does


library(CoordinateCleaner)
library(data.table)
library(sf)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)

## ------------------------------------------------------------
## Load GBIF occurrence object
## ------------------------------------------------------------

gbif_object <- readRDS("../data/scaled_GBIF_clim_pres.RDS")

coord_cols <- c("decimallongitude", "decimallatitude")

if(is.null(names(gbif_object))) {
  stop("gbif_object is not a named list.")
}

## ------------------------------------------------------------
## Standardize GBIF list
## ------------------------------------------------------------

gbif_list <- lapply(gbif_object, function(x) {
  
  if(is.null(dim(x))) {
    
    if(is.null(names(x))) {
      stop("Found an unnamed vector entry in gbif_object.")
    }
    
    x <- as.data.frame(as.list(x), stringsAsFactors = FALSE)
    
  } else {
    
    x <- as.data.frame(x, stringsAsFactors = FALSE)
  }
  
  names(x) <- trimws(names(x))
  row.names(x) <- NULL
  
  x
})

## ------------------------------------------------------------
## Bind into one dataframe
## ------------------------------------------------------------

gbif_df <- rbindlist(
  lapply(names(gbif_list), function(sp) {
    x <- gbif_list[[sp]]
    x$species_gbif <- sp
    x
  }),
  fill = TRUE
)

gbif_df <- as.data.frame(gbif_df)

missing_coord_cols <- setdiff(coord_cols, names(gbif_df))

if(length(missing_coord_cols) > 0) {
  stop(
    "Missing coordinate columns after binding: ",
    paste(missing_coord_cols, collapse = ", ")
  )
}

gbif_df$decimallongitude <- as.numeric(gbif_df$decimallongitude)
gbif_df$decimallatitude  <- as.numeric(gbif_df$decimallatitude)

## ------------------------------------------------------------
## Basic coordinate flags
## ------------------------------------------------------------

gbif_df$bad_missing_coord <- is.na(gbif_df$decimallongitude) |
  is.na(gbif_df$decimallatitude)

gbif_df$bad_lon_range <- gbif_df$decimallongitude < -180 |
  gbif_df$decimallongitude > 180

gbif_df$bad_lat_range <- gbif_df$decimallatitude < -90 |
  gbif_df$decimallatitude > 90

gbif_df$bad_zero_zero <- gbif_df$decimallongitude == 0 &
  gbif_df$decimallatitude == 0

gbif_df$bad_equal_lonlat <- gbif_df$decimallongitude ==
  gbif_df$decimallatitude

gbif_df$record_id <- seq_len(nrow(gbif_df))

## ------------------------------------------------------------
## Run CoordinateCleaner
## ------------------------------------------------------------

cc_input <- gbif_df[
  !gbif_df$bad_missing_coord &
    !gbif_df$bad_lon_range &
    !gbif_df$bad_lat_range,
  ,
  drop = FALSE
]


#filter through coordinates

cc_valid <- clean_coordinates(
  x = cc_input,
  lon = "decimallongitude",
  lat = "decimallatitude",
  species = "species_gbif",
  tests = c(
    "capitals",
    "centroids",
    "equal",
    "gbif",
    "institutions",
    "seas",
    "zeros"
  ),
  value = "spatialvalid",
  verbose = TRUE
)

## In this CoordinateCleaner version, cc_valid is a dataframe with:
## .val, .equ, .zer, .cap, .cen, .sea, .gbf, .inst, .summary

if(!".summary" %in% names(cc_valid)) {
  stop("CoordinateCleaner output does not contain .summary column.")
}

## ------------------------------------------------------------
## Overall failure summary
## ------------------------------------------------------------

# some of these tests we really need to not fail, (ie invalid coordinates, equal lonlat, zero, gbif), 
# but again others (capital, centroid, sea,  institution) could realistically fail for soem points

overall_failures <- data.frame(
  test = c(
    "invalid_coordinate",
    "equal_lonlat",
    "zero",
    "capital",
    "centroid",
    "sea",
    "gbif",
    "institution",
    "any"
  ),
  n_fail = c(
    sum(cc_valid$.val == FALSE, na.rm = TRUE),
    sum(cc_valid$.equ == FALSE, na.rm = TRUE),
    sum(cc_valid$.zer == FALSE, na.rm = TRUE),
    sum(cc_valid$.cap == FALSE, na.rm = TRUE),
    sum(cc_valid$.cen == FALSE, na.rm = TRUE),
    sum(cc_valid$.sea == FALSE, na.rm = TRUE),
    sum(cc_valid$.gbf == FALSE, na.rm = TRUE),
    sum(cc_valid$.inst == FALSE, na.rm = TRUE),
    sum(cc_valid$.summary == FALSE, na.rm = TRUE)
  ),
  stringsAsFactors = FALSE
)

overall_failures$prop_fail <- overall_failures$n_fail / nrow(cc_valid)

cat("\n--- Overall CoordinateCleaner failures ---\n")
print(overall_failures)

#total number is around 1013/60k+ not bad most of them are capital with soem being sea, lets looks at the species that fail capital
# none are failing based on our absolute avoids!

## ------------------------------------------------------------
## Species-level failure summary
## ------------------------------------------------------------

flag_summary <- do.call(
  rbind,
  lapply(split(cc_valid, cc_valid$species_gbif), function(x) {
    
    data.frame(
      species = x$species_gbif[1],
      n_records = nrow(x),
      
      invalid_fail    = sum(x$.val     == FALSE, na.rm = TRUE),
      equal_fail      = sum(x$.equ     == FALSE, na.rm = TRUE),
      zero_fail       = sum(x$.zer     == FALSE, na.rm = TRUE),
      capital_fail    = sum(x$.cap     == FALSE, na.rm = TRUE),
      centroid_fail   = sum(x$.cen     == FALSE, na.rm = TRUE),
      sea_fail        = sum(x$.sea     == FALSE, na.rm = TRUE),
      gbif_fail       = sum(x$.gbf     == FALSE, na.rm = TRUE),
      institution_fail = sum(x$.inst   == FALSE, na.rm = TRUE),
      
      total_fail      = sum(x$.summary == FALSE, na.rm = TRUE),
      
      stringsAsFactors = FALSE
    )
  })
)

row.names(flag_summary) <- NULL

flag_summary$prop_fail <- flag_summary$total_fail /
  flag_summary$n_records

flag_summary <- flag_summary[
  order(-flag_summary$total_fail, -flag_summary$prop_fail),
]

flag_summary_nonzero <- flag_summary[
  flag_summary$total_fail > 0,
]

cat("\n--- Species with at least one CoordinateCleaner failure ---\n")
print(flag_summary_nonzero)

# the fails are in expected species, redbacks, desmogs, that are found around US and/or canadian capitals

## ------------------------------------------------------------
## All bad records dataframe
## ------------------------------------------------------------

bad_points <- cc_valid[
  cc_valid$.summary == FALSE,
]

bad_points_out <- bad_points[
  ,
  c(
    "species_gbif",
    "decimallongitude",
    "decimallatitude",
    "bio1",
    "bio12",
    ".val",
    ".equ",
    ".zer",
    ".cap",
    ".cen",
    ".sea",
    ".gbf",
    ".inst",
    ".summary"
  )
]

bad_points_out$failed_tests <- apply(
  bad_points_out[
    ,
    c(".val", ".equ", ".zer", ".cap", ".cen", ".sea", ".gbf", ".inst")
  ],
  1,
  function(z) {
    labs <- c(
      "invalid_coordinate",
      "equal_lonlat",
      "zero",
      "capital",
      "centroid",
      "sea",
      "gbif",
      "institution"
    )
    
    paste(labs[z == FALSE], collapse = ";")
  }
)

cat("\n--- Total bad records ---\n")
cat("Bad records:", nrow(bad_points_out), "\n")

## ------------------------------------------------------------
## Sea-flagged visual inspection
## ------------------------------------------------------------

#lets make sure none of the sea-points are glaringly wrong, if they are clippingt he land polygon for a species found on the coast, thats reasonable

sea_bad <- cc_valid[
  cc_valid$.sea == FALSE,
]

cat("\n--- Sea-flagged records ---\n")
cat("Sea-flagged records:", nrow(sea_bad), "\n")

if(nrow(sea_bad) > 0) {
  
  sea_bad$decimallongitude <- as.numeric(sea_bad$decimallongitude)
  sea_bad$decimallatitude  <- as.numeric(sea_bad$decimallatitude)
  
  sea_bad_sf <- st_as_sf(
    sea_bad,
    coords = c("decimallongitude", "decimallatitude"),
    crs = 4326,
    remove = FALSE
  )
  
  land <- ne_countries(
    scale = "medium",
    returnclass = "sf"
  )
  
  land <- land[
    land$admin %in% c(
      "United States of America",
      "Canada"
    ),
  ]
  
  if(nrow(land) == 0) {
    stop("US/Canada polygons were not found.")
  }
  
  land_crop <- st_crop(
    land,
    xmin = -105,
    xmax = -50,
    ymin = 20,
    ymax = 60
  )
  
  sf::sf_use_s2(FALSE)
  
  land_crop <- st_make_valid(land_crop)
  sea_bad_sf <- st_make_valid(sea_bad_sf)
  
  land_5070 <- st_transform(land_crop, 5070)
  sea_5070  <- st_transform(sea_bad_sf, 5070)
  
  land_5070 <- st_make_valid(land_5070)
  sea_5070  <- st_make_valid(sea_5070)
  
  dist_mat_m <- st_distance(
    sea_5070,
    land_5070
  )
  
  nearest_land_dist_m <- apply(
    as.matrix(dist_mat_m),
    1,
    min,
    na.rm = TRUE
  )
  
  sea_bad$distance_to_us_canada_land_km <-
    nearest_land_dist_m / 1000
  
  sea_bad$distance_bin <- cut(
    sea_bad$distance_to_us_canada_land_km,
    breaks = c(-Inf, 0.1, 1, 5, 10, 25, 50, Inf),
    labels = c(
      "<0.1 km",
      "0.1-1 km",
      "1-5 km",
      "5-10 km",
      "10-25 km",
      "25-50 km",
      ">50 km"
    ),
    right = TRUE
  )
  
  cat("\nDistance to nearest US/Canada land, km:\n")
  print(summary(sea_bad$distance_to_us_canada_land_km))
  
  cat("\nDistance bins:\n")
  print(table(sea_bad$distance_bin, useNA = "ifany"))
  
  cat("\nTop species by sea flags:\n")
  print(sort(table(sea_bad$species_gbif), decreasing = TRUE))
  
  sea_bad_sf <- st_as_sf(
    sea_bad,
    coords = c("decimallongitude", "decimallatitude"),
    crs = 4326,
    remove = FALSE
  )
  
  p_sea_map <- ggplot() +
    geom_sf(
      data = land_crop,
      fill = "grey95",
      color = "grey40",
      linewidth = 0.3
    ) +
    geom_sf(
      data = sea_bad_sf,
      aes(color = distance_bin),
      size = 2,
      alpha = 0.85
    ) +
    coord_sf(
      xlim = c(-105, -50),
      ylim = c(20, 60),
      expand = FALSE
    ) +
    labs(
      title = "CoordinateCleaner sea-flagged occurrence records",
      subtitle = "US + Canada shown; points colored by distance to nearest land polygon",
      x = "Longitude",
      y = "Latitude",
      color = "Distance to land"
    ) +
    theme_bw()
  
  print(p_sea_map)
  
  ggsave(
    "sea_flagged_occurrence_records_us_canada.png",
    p_sea_map,
    width = 10,
    height = 7,
    dpi = 300
  )
  
  p_sea_zoom <- ggplot() +
    geom_sf(
      data = land_crop,
      fill = "grey95",
      color = "grey40",
      linewidth = 0.3
    ) +
    geom_sf(
      data = sea_bad_sf,
      aes(color = distance_bin),
      size = 2,
      alpha = 0.85
    ) +
    coord_sf(
      xlim = c(-95, -60),
      ylim = c(25, 50),
      expand = FALSE
    ) +
    labs(
      title = "Sea-flagged occurrence records: eastern North America zoom",
      subtitle = "Points colored by distance to nearest US/Canada land polygon",
      x = "Longitude",
      y = "Latitude",
      color = "Distance to land"
    ) +
    theme_bw()
  
  print(p_sea_zoom)
  
  ggsave(
    "sea_flagged_occurrence_records_us_canada_zoom.png",
    p_sea_zoom,
    width = 10,
    height = 7,
    dpi = 300
  )
  
  p_sea_dist <- ggplot(
    sea_bad,
    aes(x = distance_to_us_canada_land_km)
  ) +
    geom_histogram(
      bins = 30,
      color = "black",
      fill = "grey80"
    ) +
    labs(
      title = "Distance of sea-flagged records to nearest US/Canada land",
      x = "Distance to nearest land polygon (km)",
      y = "Number of records"
    ) +
    theme_bw()
  
  print(p_sea_dist)
  
  ggsave(
    "sea_flagged_occurrence_records_distance_histogram.png",
    p_sea_dist,
    width = 8,
    height = 5,
    dpi = 300
  )
  
  write.csv(
    sea_bad,
    "sea_flagged_occurrence_records_with_distance.csv",
    row.names = FALSE
  )
}

# yeah none of the points are particularly surprising, I am choosing to keep them

## ------------------------------------------------------------
## Save outputs
## ------------------------------------------------------------

write.csv(
  overall_failures,
  "coordinatecleaner_failures_overall.csv",
  row.names = FALSE
)

write.csv(
  flag_summary,
  "coordinatecleaner_failures_by_species.csv",
  row.names = FALSE
)

write.csv(
  bad_points_out,
  "coordinatecleaner_all_bad_occurrence_records.csv",
  row.names = FALSE
)

cat("\nDone.\n")



#after reviewing the coordinate cleaner flags I decided not to drop any data from the final occurrence set, it all looks fairly reasonable!
