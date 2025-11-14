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


# Edit 2: Robust latitude-band assignment for H1
# Replaced text-based classification with coordinate-based bands (|lat| < 23.5 = Tropical)
# Added centroid-based imputation for records missing coordinates
                          
suppressPackageStartupMessages({
  library(dplyr)
  library(sf)
  library(rnaturalearth)
  library(rnaturalearthdata)
  library(stringr)
})

# helper: normalize names for fuzzy join
.norm <- function(x) {
  x <- tolower(gsub("[^A-Za-z ]", " ", x))
  str_squish(x)
}

df_geo <- df_geo %>% mutate(.row_id = dplyr::row_number())

# 1) derive band directly from latitude (safe numeric)
df_geo <- df_geo %>%
  mutate(lat_num = suppressWarnings(as.numeric(lat))) %>%
  mutate(
    band_from_coords = dplyr::case_when(
      !is.na(lat_num) & abs(lat_num) < 23.5 ~ "Tropical",
      !is.na(lat_num)                        ~ "Temperate",
      TRUE ~ NA_character_
    )
  )

# 2) impute band using country centroid latitude when coords are missing
# choose a country-like column if available
if ("country" %in% names(df_geo)) {
  country_col <- "country"
} else if ("country/ocean" %in% names(df_geo)) {
  country_col <- "country/ocean"
} else {
  country_col <- NA_character_
}

if (!is.na(country_col)) {
  # Natural Earth countries: fix invalid geometries, use point-on-surface for stable "centroid"
  ne <- rnaturalearth::ne_countries(scale = 110, returnclass = "sf") %>%
    dplyr::select(name_long, geometry) %>%
    sf::st_make_valid()
  rep_pts <- sf::st_point_on_surface(ne)
  ne <- ne %>%
    mutate(lat_cen = as.numeric(sf::st_coordinates(rep_pts)[, 2])) %>%
    select(name_long, lat_cen) %>%
    mutate(.ne_norm = .norm(name_long))
  
  df_geo <- df_geo %>%
    mutate(.country_norm = .norm(.data[[country_col]])) %>%
    left_join(ne %>% select(.ne_norm, lat_cen),
              by = c(".country_norm" = ".ne_norm")) %>%
    mutate(
      band_from_country = dplyr::case_when(
        is.na(band_from_coords) & !is.na(lat_cen) & abs(lat_cen) < 23.5 ~ "Tropical",
        is.na(band_from_coords) & !is.na(lat_cen)                        ~ "Temperate",
        TRUE ~ NA_character_
      )
    )
} else {
  df_geo <- df_geo %>% mutate(band_from_country = NA_character_)
}

# 3) finalize band label (coords > country centroid > existing)
df_geo <- df_geo %>%
  mutate(
    band = dplyr::coalesce(band_from_coords, band_from_country),
    band = factor(band, levels = c("Tropical", "Temperate"))
  ) %>%
  select(-.row_id, -lat_num, -band_from_coords, -band_from_country,
         -lat_cen, -.country_norm)
                          
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

# Edit 1: Rarefaction-based richness test
# Added rarefaction and permutation test to standardize sampling effort and formally test H1 (Tropical > Temperate)
                          
# 0) Keep only valid bands (Tropical/Temperate) and valid basins
valid <- df_use %>%
  filter(!is.na(band), band %in% c("Tropical", "Temperate")) %>%
  filter(!is.na(basin), basin != "Unknown/Offshore")

if (nrow(valid) == 0L || length(unique(valid$band)) < 2L) {
  cat("# H1 rarefaction: skipped (need both Tropical and Temperate)\n")
} else {
  
  # 1) Choose a global draw size based on the smallest (basin, band) cell
  cell_counts <- valid %>% count(basin, band, name = "n_cell")
  n_draw <- max(5L, min(cell_counts$n_cell, na.rm = TRUE))  # lower bound avoids trivial draws
  cat("# H1 rarefaction n_draw =", n_draw, "\n")
  
  # 2) Rarefaction helper: subsample to n_draw and count unique BINs
  rarefy_one <- function(d, n_draw, reps = 499) {
    if (nrow(d) < n_draw) return(NA_real_)
    mean(replicate(reps, {
      s <- d[sample.int(nrow(d), n_draw), , drop = FALSE]
      dplyr::n_distinct(s$bin_uri)
    }), na.rm = TRUE)
  }
  
  # 3) Rarefy per (basin, band), then average within each band
  rich_by_band <- valid %>%
    group_by(band, basin) %>%
    summarise(rich_raref = rarefy_one(cur_data_all(), n_draw, reps = 499),
              .groups = "drop") %>%
    filter(!is.na(rich_raref)) %>%
    group_by(band) %>%
    summarise(mean_rich = mean(rich_raref), .groups = "drop")
  
  # 4) Test only if both bands have finite means
  if (nrow(rich_by_band) == 2L && all(is.finite(rich_by_band$mean_rich))) {
    set.seed(6210)
    
    # Observed effect: Tropical − Temperate (rarefied means)
    obs <- with(rich_by_band,
                mean_rich[band == "Tropical"] - mean_rich[band == "Temperate"])
    
    # Permutation test (one-sided; H1: Tropical > Temperate)
    perm <- replicate(3000, {
      shuffled <- sample(rich_by_band$band)
      mean(rich_by_band$mean_rich[shuffled == "Tropical"]) -
        mean(rich_by_band$mean_rich[shuffled == "Temperate"])
    })
    p_value <- mean(perm <= 0)  # proportion of permuted diffs <= 0
    
    cat("H1 rarefied richness difference (Tropical − Temperate): ",
        round(obs, 2), "   p-value ≈ ", round(p_value, 3), "\n", sep = "")
  } else {
    cat("# H1 rarefaction: insufficient finite means across bands.\n")
    print(rich_by_band)
  }
}                          
                          
## H2: (Adjacent > Distant?) Are Jaccard similarities higher for adjacent basins?


# Edit 3: Graph-based adjacency analysis for H2
# Replaced manual/heuristic adjacency rules with an explicit adjacency graph
# Based on standardized ocean basins (Atlantic / Indian / Pacific / Arctic / Southern)
# Added effect size (Adjacent − Distant), bootstrap 95% CI, and permutation test

suppressPackageStartupMessages(library(dplyr))

if (!exists("sim_jac")) {
  cat("# H2 skipped: sim_jac not found\n")

} else {
  sim_mat <- as.matrix(sim_jac)

  # 1) Standardize basin names to five major oceans
  std <- function(x) case_when(
    grepl("Atlantic", x, TRUE) ~ "Atlantic",
    grepl("Pacific",  x, TRUE) ~ "Pacific",
    grepl("Indian",   x, TRUE) ~ "Indian",
    grepl("Arctic",   x, TRUE) ~ "Arctic",
    grepl("Southern", x, TRUE) ~ "Southern",
    TRUE ~ NA_character_
  )
  rownames(sim_mat) <- std(rownames(sim_mat))
  colnames(sim_mat) <- std(colnames(sim_mat))
  sim_mat <- sim_mat[!is.na(rownames(sim_mat)), !is.na(colnames(sim_mat)), drop = FALSE]
  sim_mat <- sim_mat[!duplicated(rownames(sim_mat)), !duplicated(colnames(sim_mat)), drop = FALSE]

  if (min(dim(sim_mat)) < 2) {
    cat("# H2 skipped: fewer than 2 basins after standardization\n")
  } else {

    # 2) Build unique basin pairs
    pair_df <- as.data.frame(as.table(sim_mat), stringsAsFactors = FALSE) |>
      rename(a = Var1, b = Var2, sim = Freq) |>
      filter(a != b) |>
      mutate(a2 = pmin(a, b), b2 = pmax(a, b)) |>
      distinct(a2, b2, .keep_all = TRUE)

    # 3) Define adjacency edges (Atlantic–Pacific excluded)
    edges_all <- rbind(
      c("Atlantic","Arctic"),
      c("Atlantic","Indian"),
      c("Atlantic","Southern"),
      c("Pacific","Arctic"),
      c("Pacific","Indian"),
      c("Pacific","Southern"),
      c("Indian","Southern"),
      c("Arctic","Southern")
    )
    present <- sort(unique(c(pair_df$a, pair_df$b)))
    keep <- apply(edges_all, 1, function(e) all(e %in% present))
    edges <- edges_all[keep, , drop = FALSE]

    # 4) Label Adjacent vs Distant
    is_adj <- apply(pair_df[, c("a","b")], 1, function(x)
      any(apply(edges, 1, function(e) all(x == e) || all(rev(x) == e))))
    pair_df$group <- ifelse(is_adj, "Adjacent", "Distant")

    if (!any(is_adj) || all(is_adj)) {
      cat("# H2 skipped: need both Adjacent and Distant groups\n")
    } else {

      # 5) Effect size, bootstrap CI, and permutation test
      set.seed(6210)
      obs <- with(pair_df, mean(sim[group == "Adjacent"]) - mean(sim[group == "Distant"]))
      boot <- replicate(1000, {
        mA <- mean(sample(pair_df$sim[pair_df$group == "Adjacent"], replace = TRUE))
        mD <- mean(sample(pair_df$sim[pair_df$group == "Distant"],  replace = TRUE))
        mA - mD
      })
      ci <- quantile(boot, c(.025, .975), na.rm = TRUE)
      perm <- replicate(2000, {
        lab <- sample(pair_df$group)
        mean(pair_df$sim[lab == "Adjacent"]) - mean(pair_df$sim[lab == "Distant"])
      })
      p <- mean(perm >= obs)

      cat(sprintf("H2 (Adj−Dist): effect=%.3f; 95%% CI [%.3f, %.3f]; p≈%.3f\n",
                  obs, ci[1], ci[2], p))
    }
  }
}
                          
##========================  FINISH =================

message("Done. Figures saved in ../figures/")
