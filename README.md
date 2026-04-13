# Dynamic Lake Ice Conditions Shape Caribou Water-Crossing Behavior in the Arctic

**Archive DOI:** [![DOI](https://zenodo.org/badge/19364012.svg)](https://doi.org/10.5281/zenodo.19364012)  
**Article DOI:** [10.1111/gcb.70858](https://doi.org/10.1111/gcb.70858)  *(in production)*

**License:** MIT

<b>Authors:</b> Qianru Liao<sup>1*</sup>, Eliezer Gurarie<sup>2</sup>, William F. Fagan<sup>1</sup>

<sup>1</sup>Department of Biology, University of Maryland, College Park, MD, United States

<sup>2</sup>State University of New York, College of Environmental Science and Forestry, NY, United States 


<sup>*</sup>Correspondence Author

This repository contains code, derived products, and figure scripts supporting the analyses in the manuscript:

**Dynamic Lake Ice Conditions Shape Caribou Water-Crossing Behavior in the Arctic**

The workflow links long-term caribou GPS tracking data with MODIS-derived lake albedo dynamics to identify lake-crossing and circumnavigation events at Contwoyto Lake and to quantify how seasonal ice conditions shape movement decisions.

---

## Abstract
Successful animal migration hinges on navigation and decision-making in dynamic environments. Yet, how individuals navigate transient, fine-scale landscape barriers, such as seasonally ice-covered water bodies, remains poorly understood. Understanding these responses is critical for forecasting migration routes and connectivity under global change. In the Arctic, rising temperatures are causing earlier ice melt and later freeze-up, reshaping landscape permeability and potentially disrupting migration routes for overland migrants, such as barren-ground caribou (Rangifer tarandus), a keystone Arctic species, which relies on frozen lakes and rivers for efficient spring travel to calving grounds. While caribou generally prefer ice to open water, behavioral responses to changing ice conditions have not been quantitatively assessed. We analyzed 20 years (2001-2021) of GPS data for 406 adult caribou and daily MODIS land surface albedo to examine lake-crossing decisions at Contwoyto Lake, a long (>100 km) glacial lake in northern Canada. We classified transit events as crossing or circumnavigation based on GPS trajectories relative to lake boundaries and linked behavioral decisions to spatially and temporally resolved ice conditions. Our models revealed distinct seasonal drivers. Spring crossing decisions were shaped by intermediate-scale ice conditions, with a behavioral threshold at a path-averaged annual albedo percentile rank of 0.56, corresponding to intermediate late-spring melt conditions when lake ice transitions from continuous cover toward fragmented surfaces. In fall, when the lake was ice-free, movement-related factors such as relative speeds along alternative routes better explained behavior. Our findings show how ice acts as a seasonal behavior filter, shaping functional connectivity through perceptual and energetic constraints. Although developed for caribou, this framework is transferable across species and systems. By linking high-resolution, spatiotemporal remote sensing to individual behavior, our framework identifies quantitative behavioral thresholds in response to dynamic, climate-sensitive landscape features, supporting predictive monitoring of climate-driven shifts in migratory behavior and emerging constraints on movement.

---

## Repository overview

This repository is intended as a transparent analysis archive for the manuscript. It includes the scripts used to preprocess environmental data, identify lake-use events, build event-level datasets, fit predictive models, and generate manuscript figures.

Because some input data are restricted and several steps require manual review, this repository should not be interpreted as a one-click fully automated pipeline. Instead, it provides the full analytical logic, code base, and derived files needed to understand and reproduce the workflow as closely as possible.

---

## Repository structure

```text
Caribou_Crossing_Lake/
├── data/
│   ├── output_2001.mat
│   ├── output_2002.mat
│   ├── ...
│   ├── output_2021.mat
│   ├── Contwoyto_pure_water.shp
│   ├── Contwoyto_pure_water.shx
│   ├── Contwoyto_pure_water.dbf
│   ├── Contwoyto_pure_water.prj
│   ├── all_focal_events_with_lake_spent.csv
│   ├── feature_df.csv
│   ├── known_event_all.csv
│   ├── known_event_spring.csv
│   ├── known_event_fall.csv
│   ├── unknown_event_all.csv
│   ├── unknown_event_spring.csv
│   └── unknown_event_fall.csv
├── Figure/
├── Figure_Script/
│   ├── Figure1_maps.R
│   ├── Figure2_example_paths.R
│   ├── Figure3_albedo_processing.R
│   ├── Figure4_known_event_maps.R
│   ├── Figure5_unknown_event_trajectories.R
│   ├── Figure6_spring_model_results.R
│   └── Supplementary_Figures.R
├── landcover/
│   ├── MCD12Q1-061-QC-lookup.csv
│   └── MCD12Q1.061_LC_Type1_doy2011001_aid0001.tif
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

### Environmental inputs

This workflow uses MODIS-based environmental products, including:

- **MCD43A3.061** shortwave black-sky albedo
- **MCD12Q1.061** land cover

The repository includes:
	•	the land-cover raster used to generate the lake water mask
	•	annual processed .mat files (output_2001.mat to output_2021.mat) used in downstream analyses

Large raw annual MODIS albedo tiles are not distributed here as a complete raw archive. Users who wish to rerun the full preprocessing workflow from the original raw rasters will need to obtain the relevant MODIS products separately and organize them according to the expected directory structure used in 00_MODIS_Albedo_Processing.m.

### Derived files and intermediate results

The data/ directory also contains selected derived outputs and intermediate files used during the analysis. Some of these are included as reference materials and may not exactly match the final cleaned workflow after later code reorganization and updates. The finalized analytical logic is represented by the scripts in the main project directory.

---

## External methodological reference
The cloudy-sky albedo gap-filling workflow implemented in 00_MODIS_Albedo_Processing.m was adapted from:

Jia, A., et al. (2023). Improved cloudy-sky snow albedo estimates using passive microwave and VIIRS data. ISPRS Journal of Photogrammetry and Remote Sensing, 196, 340-355.

Users who reuse or adapt this preprocessing approach should cite that paper in addition to citing this repository and the associated manuscript.

---

## Core workflow

### Step 0. Albedo preprocessing

#### `00_MODIS_Albedo_Processing.m`

This MATLAB script:

- generates a water mask from the MODIS land-cover raster
- compiles daily annual shortwave black-sky albedo layers
- applies quality filtering
- fills cloudy-sky gaps using a climatology-informed Kalman filtering approach
- exports annual daily albedo outputs as `output_<year>.mat`

These processed annual `.mat` files are used as environmental inputs in the later R scripts.

---

### Step 1. Visual inspection of tracks

#### `01_Visual_inspection.R`

This script visualizes annual caribou tracks relative to Contwoyto Lake and supports manual identification of candidate on-lake events.

Main purpose:

- display annual trajectories in an interactive map
- support manual screening of possible lake-use events
- export candidate IDs for downstream checking

---

### Step 2. Identify on-lake points and extract albedo

#### `02_Points_within_lake.R`

This script:

- identifies GPS points falling within the Contwoyto Lake polygon during DOY 92 to 280
- creates a year-by-ID candidate list
- loads processed annual albedo `.mat` files
- calculates pixel-level albedo percentile rank (APR)
- extracts APR values associated with on-lake points

**Inputs**

- restricted `nwt_lakes.RData`
- `data/Contwoyto_pure_water.shp` and companion shapefile files
- annual processed albedo files `output_<year>.mat`

---

### Step 3. Define known crossing events

#### `03_Crossing_event.R`

This script identifies confirmed crossing events and computes event-level movement and environmental metrics for known crossings.

Derived variables include:

- Crossing duration
- Lake width
- Straight distance and speed
- Least-cost circumnavigation distance and speed
- APR along the potential crossing path
- APR at nearest lake pixels
- APR of the entire lake area 
- Reference steps before and after transit events

**Inputs**

- restricted `nwt_lakes.RData`
- lake shapefile
- annual processed albedo files
- `on_lake.csv` or equivalent checked crossing table

---

### Step 3b. Define unknown crossing events

#### `03b_unknown_crossing_event.R`

This script identifies lake-intersection events that could not be confidently classified in advance as either known crossings or known circumnavigations.

The workflow includes:

1. broad candidate generation from seasonal lake-intersection cases  
2. manual review support for indexing candidate events  
3. event-level metric extraction after manual confirmation  

Manual review is required before generating the final unknown-event dataset.

**Main intermediate outputs**

- candidate tables for manual review
- manually checked unknown-event table with `before_index` and `after_index`

---

### Step 4. Define circumnavigation events

#### `04_Circumnavigating_event.R`

This script identifies and processes circumnavigation events.

The workflow includes:

- automated candidate screening using a lake buffer and long-axis intersection
- manual visual confirmation
- assignment of event start and end indices
- extraction of event-level movement and APR metrics
- generation of reference rows around focal events

Manual review is required before final event extraction.

**Inputs**

- restricted `nwt_lakes.RData`
- lake shapefile
- annual processed albedo files
- manually checked circumnavigation candidate table
  
---

### Step 5. Build combined event dataset

#### `05_prepare_modeling_dataset.R`

This script combines the crossing, unknown-crossing, and circumnavigation event tables into a unified event-level dataset.

It also:

- creates event IDs
- derives reference-based movement summaries
- estimates time spent in the lake area using `ctmm` occurrence distributions
- exports combined and focal-event tables for downstream modeling

**Inputs**

- `crossing_event.csv`
- `unknown_crossing_event.csv`
- `circumnavigate_event.csv`
- restricted `nwt_lakes.RData`
- lake shapefile

---

### Step 6. Prepare model inputs

#### `06_prepare_model_inputs.R`

This script standardizes data types, creates feature-engineered variables, and prepares analysis-ready datasets for modeling.

It outputs:

- a focal-event feature table
- known-event subsets
- unknown-event subsets
- spring-only and fall-only subsets

**Input**

- `all_focal_events_with_lake_spent.csv`

---

### Step 7. Combined models

#### `07_combined_models.R`

This script fits models using all known events across both seasons.

Methods include:

- Random Forest classification
- binomial logistic regression
- variable screening and diagnostics
- prediction for unknown events

**Inputs**

- `known_event_all.csv`
- `unknown_event_all.csv`

---

### Step 8. Spring-only models

#### `08_spring_models.R`

This script fits spring-specific Random Forest and binomial logistic models and predicts unknown spring events.

**Inputs**

- `known_event_spring.csv`
- `unknown_event_spring.csv`

**Main outputs**

- `spring_model_input.csv`
- `spring_rf_test_predictions.csv`
- `spring_unknown_predictions.csv`

---

### Step 9. Fall-only models

#### `09_fall_models.R`

This script fits fall-specific Random Forest and binomial logistic models and predicts unknown fall events.

**Inputs**

- `known_event_fall.csv`
- `unknown_event_fall.csv`

**Main outputs**

- `fall_model_input.csv`
- `fall_rf_test_predictions.csv`
- `fall_unknown_predictions.csv`

---

## Figure scripts

The `Figure_Script/` directory contains scripts used to generate manuscript figures and supplementary figures.

### Main-text figure scripts

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

Additional packages may be required for some figure-generation scripts.

---

## Reproducibility note

This repository is intended as a transparent and complete analysis archive, but not all components are fully automated.

There are three main reasons:

1. the raw caribou GPS data are restricted and cannot be redistributed here  
2. some preprocessing inputs are large and may need to be obtained separately  
3. several event-identification steps require manual review and confirmation  

For these reasons, the repository provides the full analysis logic and the main derived products used in the manuscript, while acknowledging that some steps require access to restricted data and user-guided decisions.

---

## Important notes

- Several scripts include manual verification steps, and these are required parts of the final workflow.
- In particular, `01_Visual_inspection.R`, `03b_unknown_crossing_event.R`, and `04_Circumnavigating_event.R` all require user review before downstream analyses can proceed.
- Some file paths in the scripts were originally local and may need to be updated before reuse.
- Some figures were assembled outside R after panel export.
- Some scripts preserve exploratory or alternative analytical approaches for transparency, even when the final manuscript used a simplified implementation.
- This repository is best understood as a transparent analysis archive rather than a fully automated one-command pipeline.
- Some intermediate files included in `data/` are retained as reference materials and may not perfectly match the final cleaned code structure after later reorganization.


---

## Contact

**Qianru Liao**  
Department of Biology  
University of Maryland, College Park
Email: qianru@terpmail.umd.edu

For questions about the code, workflow, or manuscript, please contact the corresponding author.
