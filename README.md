# Dynamic Lake Ice Conditions Shape Caribou Water-Crossing Behavior in the Arctic

<b>Authors:</b> Qianru Liao<sup>1*</sup>, Eliezer Gurarie<sup>2</sup>, William F. Fagan<sup>1</sup>

<sup>1</sup>Department of Biology, University of Maryland, College Park, MD, United States

<sup>2</sup>State University of New York, College of Environmental Science and Forestry, NY, United States 


<sup>*</sup>Correspondence Author

Code, derived products, and figure scripts supporting analyses of how dynamic lake ice conditions shape caribou crossing and circumnavigation behavior at Contwoyto Lake in the Canadian Arctic. This repository contains materials supporting the manuscript.

---

## Abstract
Successful animal migration hinges on navigation and decision-making in dynamic environments. Yet, how individuals navigate transient, fine-scale landscape barriers, such as seasonally ice-covered water bodies, remains poorly understood. Understanding these responses is critical for forecasting migration routes and connectivity under global change. In the Arctic, rising temperatures are causing earlier ice melt and later freeze-up, reshaping landscape permeability and potentially disrupting migration routes for overland migrants, such as barren-ground caribou (Rangifer tarandus), a keystone Arctic species, which relies on frozen lakes and rivers for efficient spring travel to calving grounds. While caribou generally prefer ice to open water, behavioral responses to changing ice conditions have not been quantitatively assessed. We analyzed 20 years (2001-2021) of GPS data for 406 adult caribou and daily MODIS land surface albedo to examine lake-crossing decisions at Contwoyto Lake, a long (>100 km) glacial lake in northern Canada. We classified transit events as crossing or circumnavigation based on GPS trajectories relative to lake boundaries and linked behavioral decisions to spatially and temporally resolved ice conditions. Our models revealed distinct seasonal drivers. Spring crossing decisions were shaped by intermediate-scale ice conditions, with a behavioral threshold at a path-averaged annual albedo percentile rank of 0.56, corresponding to intermediate late-spring melt conditions when lake ice transitions from continuous cover toward fragmented surfaces. In fall, when the lake was ice-free, movement-related factors such as relative speeds along alternative routes better explained behavior. Our findings show how ice acts as a seasonal behavior filter, shaping functional connectivity through perceptual and energetic constraints. Although developed for caribou, this framework is transferable across species and systems. By linking high-resolution, spatiotemporal remote sensing to individual behavior, our framework identifies quantitative behavioral thresholds in response to dynamic, climate-sensitive landscape features, supporting predictive monitoring of climate-driven shifts in migratory behavior and emerging constraints on movement.


---

## Repository structure

```text
Caribou_Crossing_Lake/
├── data/
├── Figure_Script/
├── landcover/
├── 00_MODIS_Albedo_Processing.m
├── 01_Visual_inspection.R
├── 02_Points_within_lake.R
├── 03_Crossing_event.R
├── 03b_unknown_crossing_event.R
├── 04_Circumnavigating_event.R
├── 05_prepare_modeling_dataset.R
├── 06_prepare_model_inputs.R
├── 07_combined_models.R
├── 08_spring_models.R
├── 09_fall_models.R
└── README.md
```

---

## Data availability and restrictions

### Caribou GPS data

Raw caribou GPS tracking data are **not included** in this repository. These data are subject to third-party sharing restrictions from the Government of the Northwest Territories, Department of Environment and Climate Change (GNWT-ECC).

Researchers interested in using these data should request access directly from GNWT-ECC.

### Remote sensing inputs

This workflow uses annual albedo products derived from MODIS datasets, including:

- **MCD43A3.061** shortwave black-sky albedo
- **MCD12Q1.061** land cover

The repository includes code for preprocessing these data, but large raw raster inputs may not be stored here depending on repository size limits.

---

## Core workflow

### Step 0. Albedo preprocessing

#### `00_MODIS_Albedo_Processing.m`

This MATLAB script:

- generates a water mask from MODIS land cover
- compiles annual daily shortwave albedo layers
- applies quality filtering
- fills cloudy-sky gaps using a Kalman filter approach informed by climatology and neighboring pixels
- exports annual daily albedo outputs as `output_<year>.mat`

These processed annual `.mat` files are the environmental inputs used by later R scripts.

---

### Step 1. Visual inspection of tracks

#### `01_Visual_inspection.R`

This script visualizes annual caribou tracks relative to Contwoyto Lake and supports manual identification of candidate on-lake events.

Main purpose:

- plot annual tracks in an interactive leaflet map
- manually record IDs that appear to traverse the lake
- export candidate IDs for downstream checking

---

### Step 2. Identify on-lake points and extract albedo

#### `02_Points_within_lake.R`

This script:

- finds all GPS points within the lake polygon during DOY 92 to 280
- creates a year-by-ID candidate list
- loads annual processed albedo data
- calculates pixel-level albedo percentile rank (APR)
- extracts APR values for points falling within the lake

Main outputs include:

- `on_lake_candidates.csv`
- `on_lake.csv`

---

### Step 3. Define known crossing events

#### `03_Crossing_event.R`

This script identifies confirmed crossing events and computes event-level metrics, including:

- Crossing duration
- Lake width
- Straight distance and speed
- Least-cost circumnavigation distance and speed
- APR along the potential crossing path
- APR at nearest lake pixels
- APR of the entire lake area 
- Reference steps before and after transit events

Main output:

- `crossing_event.csv`

---

### Step 3b. Define unknown crossing events

#### `03b_unknown_crossing_event.R`

This script identifies and processes lake-intersection events that could not be confidently classified *a priori* as either confirmed crossings or confirmed circumnavigations.

The workflow includes three parts:

1. **Broad candidate generation**  
   All Year-ID-season combinations whose seasonal trajectory intersects the lake polygon are identified as candidate unknown crossing events.

2. **Manual review support**  
   The script provides helper functions, including `suggest_crossing_indices()`, to visualize candidate trajectories, label seasonal GPS points by index, and suggest possible `before_index` and `after_index` values based on trajectory-lake intersection geometry.

3. **Event-level metric extraction**  
   After manual review, the checked candidate table is used to calculate focal-event and reference-step variables, including:
   - crossing duration
   - lake width
   - straight-line distance and speed
   - least-cost circumnavigation distance and speed
   - APR along the potential crossing path
   - APR at the nearest lake pixel
   - APR of the entire lake area  

Manual review is required before the final event table can be built.

Main intermediate files:

- `lake_intersection_candidates.csv`
- `unknown_crossing_candidates_checked.csv`

Main output:

- `unknown_crossing_event.csv`

---

### Step 4. Define circumnavigation events

#### `04_Circumnavigating_event.R`

This script identifies candidate circumnavigation events using spatial rules, followed by manual confirmation.

The workflow includes:

- automated candidate screening using a lake buffer and long-axis intersection
- manual trajectory review
- assignment of start and end indices for confirmed events
- extraction of event-level movement and APR variables for confirmed circumnavigation events
- generation of reference rows around focal events

Main output:

- `circumnavigate_event.csv`

---

### Step 5. Build combined event dataset

#### `05_prepare_modeling_dataset.R`

This script combines crossing, unknown, and circumnavigation event tables into a unified event dataset.

It also:

- creates unique `event_id` values
- derives reference-based movement summaries
- estimates time spent in lake area using `ctmm` occurrence distributions
- creates focal-event and full-event output tables

Main outputs:

- `all_event_combined.csv`
- `all_focal_events_with_lake_spent.csv`
- `chapter1_cleaned_events.rda`

---

### Step 6. Prepare model inputs

#### `06_prepare_model_inputs.R`

This script standardizes data types, derives predictor variables, and creates analysis-ready subsets for known and unknown events.

It outputs:

- full focal-event feature table
- known vs unknown subsets
- spring-only and fall-only subsets

Main outputs:

- `feature_df.csv`
- `known_event_all.csv`
- `unknown_event_all.csv`
- `known_event_spring.csv`
- `known_event_fall.csv`
- `unknown_event_spring.csv`
- `unknown_event_fall.csv`

---

### Step 7. Combined models

#### `07_combined_models.R`

This script fits models using all known events across both seasons.

Methods include:

- Random Forest classification
- binomial logistic regression
- variable screening and diagnostics
- prediction for unknown events

Main outputs:

- `combined_model_input.csv`
- `combined_rf_test_predictions.csv`
- `combined_unknown_predictions.csv`

---

### Step 8. Spring-only models

#### `08_spring_models.R`

This script fits spring-specific Random Forest and binomial logistic models and predicts unknown spring events.

Main outputs:

- `spring_model_input.csv`
- `spring_rf_test_predictions.csv`
- `spring_unknown_predictions.csv`

---

### Step 9. Fall-only models

#### `09_fall_models.R`

This script fits fall-specific Random Forest and binomial logistic models and predicts unknown fall events.

Main outputs:

- `fall_model_input.csv`
- `fall_rf_test_predictions.csv`
- `fall_unknown_predictions.csv`

---

## Figure scripts

The repository includes separate scripts for generating manuscript figures and supplementary figures.

### Main-text figure scripts

Examples include:

- `Figure1_maps.R`
- `Figure2_example_paths.R`
- `Figure3_albedo_processing.R`
- `Figure4_known_event_maps.R`
- `Figure5_unknown_event_trajectories.R`
- `Figure6_spring_model_results.R`

These scripts generate figure panels used in the main manuscript. In some cases, final panel assembly was completed manually in external software such as PowerPoint.

### Supplementary figure scripts

Supplementary figures are generated in scripts such as:

- `Supplementary_Figures.R`
- appendix-specific scripts for threshold exploration, event summaries, and interannual 56th-threshold analyses

These include:

- exploratory albedo distributions
- alternative path methods
- fix-interval summary tables
- threshold timing maps across multiple percentile levels
- interannual 56th-percentile albedo and DOY analyses
  
---

## Expected script order

A typical end-to-end run order is:

1. `00_MODIS_Albedo_Processing.m`
2. `01_Visual_inspection.R`
3. `02_Points_within_lake.R`
4. `03_Crossing_event.R`
5. `03b_unknown_crossing_event.R`
6. `04_Circumnavigating_event.R`
7. `05_prepare_modeling_dataset.R`
8. `06_prepare_model_inputs.R`
9. `07_combined_models.R` (option)
10. `08_spring_models.R`
11. `09_fall_models.R`

Several steps require manual confirmation, especially:

- visual identification of candidate lake-use events
- manual confirmation of circumnavigation cases and their `before_index` / `after_index` values
- manual confirmation of unknown crossing candidates and their `before_index` / `after_index` values

Figure scripts can be run after the relevant upstream data objects and modeling results have been generated.

---

## Software requirements

### MATLAB

Used for MODIS albedo preprocessing:

- MATLAB with image and file I/O support

### R

Most analyses were conducted in R.

Commonly used packages include:

- `dplyr`
- `lubridate`
- `sf`
- `terra`
- `raster`
- `R.matlab`
- `leaflet`
- `ggplot2`
- `ggpubr`
- `RColorBrewer`
- `randomForest`
- `caret`
- `ROCR`
- `pROC`
- `MuMIn`
- `car`
- `glmnet`
- `MASS`
- `ctmm`
- `gdistance`
- `fasterize`
- `geosphere`
- `ggspatial`
- `rnaturalearth`
- `rnaturalearthdata`
- `ggmap`
- `svglite`
- `mapview`

---

## Reproducibility note

Because some input data are restricted and several steps require manual inspection, this repository is best understood as a transparent analysis archive rather than a fully automated reproducible pipeline. All scripts used in the manuscript are included, but some intermediate decisions require user input and some upstream data must be obtained separately.

---

## Important notes

- Several scripts include manual verification steps, and these are required parts of the final workflow.
- In particular, `01_Visual_inspection.R`, `03b_unknown_crossing_event.R`, and `04_Circumnavigating_event.R` all require user review before downstream analyses can proceed.
- Some file paths in the scripts were originally local and may need to be updated before reuse.
- Some figures were assembled outside R after panel export.
- Some scripts preserve exploratory or alternative analytical approaches for transparency, even when the final manuscript used a simplified implementation.
- This repository is best understood as a transparent analysis archive rather than a fully automated one-command pipeline.

---

## Contact

**Qianru Liao**  
Department of Biology  
University of Maryland, College Park

For questions about the code, workflow, or manuscript, please contact the corresponding author.
