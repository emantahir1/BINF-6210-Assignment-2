###############################################################
## Author: Eman Tahir | Date: 2025-10
## Course: Software Tools (BINF6210)
## Project: Requiem Sharks on BOLD: BIN Richness and Basin Similarity
###############################################################

## BACKGROUND:
# Requiem sharks (family Carcharhinidae) are common warm-sea sharks with many species spread across tropical–warm temperate oceans. The BOLD database clusters DNA barcodes into BINs (species-like genetic units), which we can use to compare diversity across regions.

## WHY THIS IS INTERESTING:
# Classic “latitudinal diversity gradient” (LDG) ideas expect higher diversity in the tropics, but recent marine work shows variation across groups. Also, not all ocean basins mix freely: some lineages are “basin-linked.” Our BOLD data let us check these patterns objectively.

## MY QUESTION:
# Do major ocean basins share requiem-shark BINs, and do adjacent basins share more than distant ones? Also, do tropical areas hold more BINs than temperate ones?

## PROJECT TYPE:
# Hypothesis-based with exploratory elements.
#   H1 (latitudinal): Tropical > Temperate in BIN richness.
#   H2 (connectivity): Adjacent basins share more BINs than distant basins.
#   H0: No differences (patterns could come from uneven sampling).

# 3 storyboard figures:
#   Fig 1: BIN richness per basin (presence-only; record counts annotated)
#   Fig 2: Basin similarity heatmap (Jaccard; presence/absence of BINs)
#   Fig 3: Global sampling locations by ocean basin

## REPRODUCIBILITY NOTES:
# Keep the .Rproj open at the project root so relative paths work.
# If not using .Rproj, set the working directory once so "../data" and "../figures" resolve.

##========= PACKAGES & BASIC SETUP ===================================
## WHAT: load packages, set a clean plot theme, make randomness reproducible.

# If this is your first time on this machine, UNCOMMENT these once to install:
# install.packages("tidyverse")
# install.packages("vegan")
# install.packages("viridis")
# install.packages("maps")

## ---- Load packages (every run) ----
library(tidyverse)   # pipes, wrangling, ggplot2
library(vegan)       # vegdist()
library(viridis)     # color-blind friendly palettes
library(maps)        # supplies world polygons

theme_set(theme_light(base_size = 12))
set.seed(6210)  # reproducible resampling

##-------------------- Working directory / paths (simple) ---------------------------
# Open the project .Rproj at the project root.
# If not, run setwd("path/to/binf6210_assignment_1/R") once.
# Data is at:    ../data/Carcharhinidae_BOLD_data.tsv
# Figures to:     ../figures/

##========================  1) READ DATA + QUICK CHECKS  ===========================
## WHAT: Read the BOLD export and check basic structure.
## WHY: Confirm the columns we need exist.
dfBOLD <- read_tsv("../data/Carcharhinidae_BOLD_data.tsv")

dim(dfBOLD)                     # rows x columns
names(dfBOLD)[1:30]             # first 30 column names for transparency
stopifnot("bin_uri" %in% names(dfBOLD))  # BINs are required

# Keep only the fields we will actually use:
use_cols <- intersect(c("bin_uri", "country/ocean", "region", "coord"), names(dfBOLD))
df <- dfBOLD %>% select(all_of(use_cols))

##========================  2) BUILD BASIN LABELS (WHAT & WHY)  ====================
## WHAT: Create one "basin" label per record from:
#   (a) coordinates (objective, if present),
#   (b) explicit ocean words in text,
#   (c) a simple country -> basin map (based on this dataset).
## WHY: We need consistent basins to compare BIN richness and overlap.

## GOAL: create one reliable 'basin' label per record so we can compare BIN richness and overlap.

## 2a) STEP: Pull latitude/longitude out of the text column `coord`.
## WHY: many rows don’t have separate lat/lon fields; we need numeric coords to assign basins objectively. coord is a mixed text field; I extract the first two numbers as lat/lon to make an objective basin label later.
df_geo <- df %>%
  mutate(coord = na_if(coord, "")) %>%                                     # turn "" into NA
  mutate(
    nums = stringr::str_extract_all(coord, "[-+]?\\d*\\.?\\d+"),           # grab all numbers in the string
    lat  = purrr::map_dbl(nums, ~ if (length(.x) >= 1)                     # first number = latitude (if present)
      suppressWarnings(as.numeric(.x[1])) else NA_real_),
    lon  = purrr::map_dbl(nums, ~ if (length(.x) >= 2)                     # second number = longitude (if present)
      suppressWarnings(as.numeric(.x[2])) else NA_real_)
  ) %>%
  select(-nums)                                                             # drop the temporary list column

# QUICK CHECK (optional): look at parsed coordinates alongside the raw field
# df_geo %>% select(`country/ocean`, region, coord, lat, lon) %>% head(8)

## 2b) STEP: Define a tiny helper that turns (lat, lon) into a coarse ocean basin.
## WHY: we want a simple, reproducible rule to label Atlantic / Indian / Pacific (plus poles). 
infer_basin_from_lon <- function(lat, lon) {
  if (is.na(lat) || is.na(lon)) return(NA_character_)
  if (lat >= 66)  return("Arctic")     # very high latitudes
  if (lat <= -50) return("Southern")   # Southern Ocean band
  L <- ((lon + 180) %% 360) - 180      # wrap longitude to (-180, 180]
  if (L >= -100 & L <=  20)  return("Atlantic")
  if (L >    20 & L <= 150)  return("Indian")
  return("Pacific")                     # everything else
}

## 2c) STEP: Apply that helper to every row to get `basin_from_coords`.
## WHY: this gives us an objective basin label wherever coords exist.
df_geo <- df_geo %>%
  mutate(basin_from_coords = mapply(infer_basin_from_lon, lat, lon) %>% as.character())

## 2d) STEP: Prints a coverage report (how many usable coords, and their range).
## WHY: shows data quality and justifies later choices if coords are sparse.
df_geo %>%
  summarize(
    lat_min        = min(lat, na.rm = TRUE),
    lat_max        = max(lat, na.rm = TRUE),
    lon_min        = min(lon, na.rm = TRUE),
    lon_max        = max(lon, na.rm = TRUE),
    n_with_coords  = sum(!is.na(lat) & !is.na(lon))
  ) %>% print()

## 2e) STEP: Detect explicit ocean words in the text fields (e.g., “Atlantic Ocean”).
## WHY: when coords are missing, we can still tag a basin from the metadata text.
ocean_regex <- "(?i)\\b(Atlantic|Pacific|Indian|Arctic|Southern)\\b|\\bOcean\\b"

df_geo <- df_geo %>%
  mutate(
    region_search   = coalesce(`country/ocean`, region),
    basin_from_words = case_when(
      str_detect(coalesce(region_search, ""), ocean_regex) &
        str_detect(coalesce(region_search, ""), regex("Atlantic", TRUE)) ~ "Atlantic",
      str_detect(coalesce(region_search, ""), ocean_regex) &
        str_detect(coalesce(region_search, ""), regex("Pacific", TRUE))  ~ "Pacific",
      str_detect(coalesce(region_search, ""), ocean_regex) &
        str_detect(coalesce(region_search, ""), regex("Indian", TRUE))   ~ "Indian",
      str_detect(coalesce(region_search, ""), ocean_regex) &
        str_detect(coalesce(region_search, ""), regex("Arctic", TRUE))   ~ "Arctic",
      str_detect(coalesce(region_search, ""), ocean_regex) &
        str_detect(coalesce(region_search, ""), regex("Southern", TRUE)) ~ "Southern",
      TRUE ~ NA_character_
    )
  )

## 2f) STEP: Map country names to a dominant basin (only for countries in THIS dataset).
## WHY: fills the remaining gaps when neither coords nor ocean words are available.
country_map <- tibble::tibble(
  country_ocean_label = c(
    # Indian
    "Australia","India","Sri Lanka","Oman","United Arab Emirates","Saudi Arabia","Kuwait","Qatar",
    "Pakistan","Reunion","Madagascar","Seychelles","Tanzania","South Africa",
    # Pacific
    "Indonesia","Malaysia","Philippines","Japan","Taiwan","Papua New Guinea","New Zealand","Vietnam",
    "Brunei","Singapore","Fiji","Palau","French Polynesia","South Korea","China","Thailand","Myanmar",
    # Atlantic (+ Mediterranean as Atlantic gateway)
    "Brazil","United States","United States Virgin Islands","Canada","Mexico","Belize","Trinidad and Tobago",
    "Bahamas","Costa Rica","Panama","Nigeria","Sierra Leone","Ghana","Senegal","Portugal","Spain",
    "France","Italy","Malta","Tunisia","Algeria","Libya","Egypt","Greece","Israel",
    # Explicit ocean phrases
    "Indian Ocean","Pacific Ocean","Atlantic Ocean","North Atlantic Ocean",
    # Missing labels
    "Unspecified country","Unrecoverable"
  ),
  basin_map = c(
    rep("Indian", 14),                     # Indian
    rep("Pacific", 17),                    # Pacific
    rep("Atlantic", 25),                   # Atlantic
    "Indian","Pacific","Atlantic","Atlantic",  # ocean phrases
    NA_character_, NA_character_           # missing
  )
)

df_geo <- df_geo %>%
  left_join(country_map, by = c("country/ocean" = "country_ocean_label")) %>%
  mutate(
    ## Priority order: use coords if we have them; else ocean words; else country map.
    basin = coalesce(basin_from_coords, basin_from_words, basin_map, "Unknown/Offshore")
  ) %>%
  select(-basin_map)

## 2g) STEP: Collapse any leftover variations to the five major basins.
## WHY: keep one consistent level per basin (Atlantic / Indian / Pacific / Arctic / Southern).
df_geo <- df_geo %>%
  mutate(basin = case_when(
    str_detect(basin, regex("Atlantic", TRUE)) ~ "Atlantic",
    str_detect(basin, regex("Pacific",  TRUE)) ~ "Pacific",
    str_detect(basin, regex("Indian",   TRUE)) ~ "Indian",
    str_detect(basin, regex("Arctic",   TRUE)) ~ "Arctic",
    str_detect(basin, regex("Southern", TRUE)) ~ "Southern",
    TRUE ~ basin
  ))

## 2h) STEP: Make a latitude band (‘Tropical’ vs ‘Temperate’) for H1.
## WHY: we will compare rarefied BIN richness between these two bands.
df_geo <- df_geo %>%
  mutate(band = case_when(
    !is.na(lat) & abs(lat) <= 23.5 ~ "Tropical",
    !is.na(lat) & abs(lat) >  23.5 ~ "Temperate",
    TRUE ~ NA_character_
  ))

##========================  3) FILTER FOR ANALYSIS  ================================
## GOAL: keep only rows that make a fair presence/absence comparison possible.

## 3a) STEP: basic filters for presence/absence analysis.
## WHY: we need a BIN value, a named major basin, and enough records per basin.
min_n_per_basin <- 10

df_use <- df_geo %>%
  filter(!is.na(bin_uri)) %>%                                   # must have a BIN
  mutate(bin_uri = as.character(bin_uri)) %>%                   # ensure character
  filter(basin %in% c("Atlantic","Indian","Pacific","Arctic","Southern")) %>%  # major basins only
  add_count(basin, name = "n_records_basin") %>%                # how many records per basin?
  filter(n_records_basin >= min_n_per_basin) %>%                # drop tiny basins to avoid noisy plots
  select(-n_records_basin)

## 3b) STEP: quick transparency prints (how many BINs/records survive; per-basin counts).
## WHY: shows coverage after filtering.
cat("# unique BINs (filtered):", n_distinct(df_use$bin_uri), "\n")
cat("# total records (filtered):", nrow(df_use), "\n")
print(df_use %>% count(basin, sort = TRUE), n = 20)

can_analyze <- n_distinct(df_use$basin) >= 2                    # at least 2 basins needed

##========================  FIGURE 1: BIN richness by basin  ======================
## GOAL: show unique BINs per basin using presence only (not specimen counts).

## WHAT: count distinct BINs per basin; annotate bars with record counts for context.
## WHY: presence-only richness avoids overweighting basins with many specimens.

if (can_analyze) {
  bins_by_basin <- df_use %>%
    count(basin, name = "n_records") %>%
    left_join(
      df_use %>% group_by(basin) %>% summarise(bin_richness = n_distinct(bin_uri), .groups = "drop"),
      by = "basin"
    ) %>%
    arrange(desc(bin_richness))
  
  p1 <- ggplot(bins_by_basin,
               aes(x = forcats::fct_reorder(basin, bin_richness),
                   y = bin_richness, fill = basin)) +
    geom_col(show.legend = FALSE, width = 0.75) +
    geom_text(aes(label = paste0("n=", n_records)),
              hjust = -0.15, size = 3.2) +
    coord_flip(ylim = c(0, max(bins_by_basin$bin_richness) * 1.12)) +
    labs(
      title    = "Figure 1: Unique BINs per ocean basin",
      subtitle = "Presence-only richness (records annotated as n)",
      x = "Basin", y = "BIN richness"
    ) +
    scale_fill_viridis(discrete = TRUE) +
    theme(
      plot.title    = element_text(face = "bold"),
      plot.subtitle = element_text(size = 10),
      axis.title.y  = element_text(margin = margin(r = 8))
    )
  
  print(p1)
  ggsave("../figures/Fig1_BIN_richness_by_basin.png", p1,
         width = 9, height = 6, dpi = 300, bg = "white")
}

##=====================  FIGURE 2: Jaccard similarity heatmap  =====================
## GOAL: show how similar basin BIN lists are (shared composition; Jaccard)

# 1) Build a presence/absence table: rows = basins, columns = BINs (0/1)
comm <- df_use %>%
  count(basin, bin_uri) %>%        # tally records per basin x BIN
  mutate(pres = 1L) %>%            # presence only
  select(-n) %>%                   # drop counts (we only want presence)
  tidyr::pivot_wider(
    names_from  = bin_uri,         # one column per BIN
    values_from = pres,
    values_fill = 0
  ) %>%
  as.data.frame()

# 2) Move row names & keep only numeric columns
rownames(comm) <- comm$basin
comm$basin <- NULL
comm <- dplyr::select(comm, where(is.numeric)) %>% as.matrix()

# Quick sanity check:
cat("Matrix dims (basins x BINs):", paste(dim(comm), collapse = " x "), "\n")

if (nrow(comm) >= 2 && ncol(comm) >= 1) {

# 3) Jaccard dissimilarity -> similarity
  dist_jac <- vegan::vegdist(comm, method = "jaccard", binary = TRUE)
  sim_jac  <- 1 - as.matrix(dist_jac)
  
# 4) Nice ordering by hierarchical clustering
  ord <- hclust(as.dist(1 - sim_jac))$order
  sim_plot <- sim_jac[ord, ord]
  
# 5) Long format for ggplot
  sim_df <- as.data.frame(sim_plot) %>%
    tibble::rownames_to_column("row_basin") %>%
    tidyr::pivot_longer(-row_basin, names_to = "col_basin", values_to = "similarity")
  
  p2 <- ggplot(sim_df, aes(col_basin, row_basin, fill = similarity)) +
    geom_tile() +
    geom_text(aes(label = sprintf("%.2f", similarity)), size = 3) +
    scale_fill_viridis(name = "Jaccard similarity", limits = c(0, 1)) +
    coord_fixed() +
    labs(title = "Figure 2: BIN composition similarity among basins (Jaccard)",
         x = "Basin", y = "Basin") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  print(p2)
  ggsave("../figures/Fig2_Jaccard_similarity_heatmap.png", p2,
         width = 7, height = 6, dpi = 300, bg = "white")
} else {
  message("Figure 2 skipped: need ≥2 basins and ≥1 BIN column after filtering.")
}

##=================  FIGURE 3: Global sampling map =================
## WHAT: one world map; the points coloured by basin; ocean labels indicated

# world polygons (uses ggplot2 helper)
world <- ggplot2::map_data("world")

df_pts <- df_geo %>%
  filter(!is.na(lat), !is.na(lon)) %>%
  filter(!is.na(basin) & basin != "Unknown/Offshore")

if (nrow(df_pts) > 0) {
  p3 <- ggplot() +
    # base map
    geom_polygon(
      data = world,
      aes(long, lat, group = group),
      fill = "grey95", color = "grey70", linewidth = 0.2
    ) +

# sampling points
    geom_point(
      data = df_pts,
      aes(x = lon, y = lat, color = basin),
      alpha = 0.75, size = 1.2
    ) +
    coord_quickmap(xlim = c(-180, 180), ylim = c(-65, 80)) +
    scale_color_viridis(discrete = TRUE, name = "Basin") +

  # ocean labels (positions chosen to avoid clutter)
    annotate("text", x = -40,  y =  5,  label = "Atlantic",  fontface = "bold",
             size = 5, alpha = 0.7) +
    annotate("text", x =  80,  y =  0,  label = "Indian",    fontface = "bold",
             size = 5, alpha = 0.7) +
    annotate("text", x = -150, y = 10,  label = "Pacific",   fontface = "bold",
             size = 5, alpha = 0.7) +
    annotate("text", x =  150, y =  8,  label = "Pacific",   fontface = "bold",
             size = 5, alpha = 0.7) +
    annotate("text", x =   10, y = 70,  label = "Arctic",    fontface = "bold",
             size = 4.5, alpha = 0.7) +
    annotate("text", x =   15, y = -55, label = "Southern",  fontface = "bold",
             size = 4.5, alpha = 0.7) +
    labs(
      title = "Figure 3: Global sampling locations by ocean basin (Carcharhinidae)",
      subtitle = "Points = records with coordinates; colours = basin assigned from metadata/coords",
      x = "Longitude", y = "Latitude"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold", size = 13, hjust = 0.5)
    )
  
  print(p3)
  ggsave("../figures/Fig3_Global_sampling_map.png",
         p3, width = 10, height = 5.7, dpi = 300, bg = "white")
} else {
  message("Map skipped: no rows with both coordinates and a final basin label.")
}

## ===================== Hypothesis snapshot (in console) ======================

cat("\n# ===== RESULTS SNAPSHOT =====\n")

## H1: Wwhich basins are richest (from filtered data)?
if (exists("df_use") && nrow(df_use) > 0) {
  h1_tbl <- df_use %>%
    dplyr::summarise(bins = dplyr::n_distinct(bin_uri), .by = basin) %>%
    dplyr::arrange(dplyr::desc(bins))
  cat("# H1 (richness ranking): ", paste(h1_tbl$basin, collapse = " > "), "\n", sep = "")
} else {
  cat("# H1: skipped (no filtered data)\n")
}

## H2: (Adjacent > Distant?) Are Jaccard similarities higher for adjacent basins?

if (exists("sim_jac") && nrow(sim_jac) >= 2) {

# Pair table from similarity matrix
  pairs <- as.data.frame(as.table(sim_jac), stringsAsFactors = FALSE) |>
    dplyr::rename(a = Var1, b = Var2, sim = Freq) |>
    dplyr::filter(a != b) |>
    dplyr::mutate(a = as.character(a), b = as.character(b)) |>
    dplyr::mutate(key_sorted = paste(pmin(a,b), pmax(a,b), sep = "|")) |>
    dplyr::distinct(key_sorted, .keep_all = TRUE)
  
# Define adjacency ONLY for basins present
  present <- sort(unique(c(pairs$a, pairs$b)))
  
# Start empty, then add rules that apply to your data
  adj_keys <- character(0)
  
# Core majors
  if (all(c("Indian","Pacific") %in% present))   adj_keys <- c(adj_keys, "Indian|Pacific")
  if (all(c("North Pacific","South Pacific") %in% present)) adj_keys <- c(adj_keys, "North Pacific|South Pacific")
  if (all(c("North Atlantic","South Atlantic") %in% present)) adj_keys <- c(adj_keys, "North Atlantic|South Atlantic")
  
# Poles to neighbors
  if (all(c("Atlantic","Arctic") %in% present))  adj_keys <- c(adj_keys, "Arctic|Atlantic")
  if (all(c("Pacific","Arctic") %in% present))   adj_keys <- c(adj_keys, "Arctic|Pacific")
  if (all(c("Indian","Southern") %in% present))  adj_keys <- c(adj_keys, "Indian|Southern")
  if (all(c("Pacific","Southern") %in% present)) adj_keys <- c(adj_keys, "Pacific|Southern")
  if (all(c("Atlantic","Southern") %in% present))adj_keys <- c(adj_keys, "Atlantic|Southern")
  
# Label groups
  pairs <- pairs |>
    dplyr::mutate(group = ifelse(key_sorted %in% adj_keys, "Adjacent", "Distant"))
  
  if (any(pairs$group == "Adjacent") && any(pairs$group == "Distant")) {
    obs <- with(pairs, mean(sim[group == "Adjacent"]) - mean(sim[group == "Distant"]))
    set.seed(6210)
    perm <- replicate(3000, {
      lab <- sample(pairs$group)  # shuffle labels
      mean(pairs$sim[lab == "Adjacent"]) - mean(pairs$sim[lab == "Distant"])
    })
    p_one_sided <- mean(perm <= 0)  # H2: Adjacent > Distant
    
    cat("# H2 (Adj − Dist Jaccard): ",
        round(obs,3), " ; permutation p ≈ ", round(p_one_sided,3), "\n", sep = "")
  } else {
    cat("# H2: skipped — need both Adjacent and Distant pairs present.\n")
  }
} else {
  cat("# H2: skipped — run Fig 2 to create `sim_jac`.\n")
}

##========================  FINISH =================

message("Done. Figures saved in ../figures/")
