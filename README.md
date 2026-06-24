# Imidacloprid in English Freshwaters

This repository contains the data processing and analysis code supporting the manuscript:

> **Sources and Mitigation of Nationwide Parasiticide Water Pollution**\
> *Submitted for peer review, 2026*

## Structure

### `code/`

Contains scripts used for data processing, statistical analyses, and model fitting.

-   `EAdata.R`: This script processes Environment Agency (EA) monitoring data to extract imidacloprid concentration records and generate non-detect observations from positive-only EA records.

-   `WWTPs.R`: This script (1) processes wastewater treatment plant (WWTP) data obtained from DEFRA; and (2) performs logistic regression models assessing whether upstream WWTP presence is associated with downstream imidacloprid detection – detailed methods, interpretation, and results are provided in the manuscript.

-   `bayesian-models-WP.R`: This script fits Bayesian hurdle models relating upstream wastewater population equivalent to observed imidacloprid concentrations in receiving waters – detailed methods, interpretation, and results are provided in the manuscript.

-   `bayesian-models_mwe-GR.R`: This script conducts sensitivity analyses using alternative substitutions for concentrations below the limit of detection (LOD), revisiting the Bayesian regression models presented in the manuscript – detailed methods, interpretation, and results are provided in the manuscript SI.

### `processed_data/`

Contains processed datasets generated from the original raw data sources.

**Note:** Raw data are not included in this repository due to licensing, privacy, or size constraints. Information on obtaining the original datasets is provided in the manuscript.

-   `Raw_EA_2014OW.rds`: Imidacloprid concentration data for English freshwater monitoring sites between 2014 and 2023, including both positive detections and reconstructed non-detect observations. Generated from EA monitoring data using `code/EAdata.R`.

-   `allSamples1909RBD.rds`: Imidacloprid concentration records between 2019 and 2023 with associated River Basin District (RBD) information. RBD assignments were derived from catchment spatial data and subsequently verified in ArcGIS Pro.

-   `WWpopulation.rds`: Dataset linking upstream wastewater treatment population equivalents with downstream EA site imidacloprid concentrations. Watersheds were delineated in ArcGIS Pro (workflow described below).

### `watershed_delineation_workflow.pdf`

Documents the ArcGIS Pro workflow used to delineate upstream catchments for monitoring locations, including:

1.  Sink filling

2.  Flow direction generation

3.  Flow accumulation calculation

4.  Stream network derivation

5.  Monitoring site snapping

6.  Watershed delineation

7.  Validation

## Software Requirements

The analyses were performed using:

-   R 4.5.2

-   ArcGIS Pro 3.2

Additional R package dependencies are listed within the individual scripts.

```         
git clone https://github.com/JiayiLiu354/Imidacloprid.git 
cd Imidacloprid/processed_data
```
