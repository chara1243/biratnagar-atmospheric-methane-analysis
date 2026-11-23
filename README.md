# Biratnagar Atmospheric Methane Study (2019–2025)

## Overview

This repository contains Google Earth Engine (GEE) JavaScript code for analysing atmospheric methane (XCH4) variability and its environmental drivers in Biratnagar, Koshi Province, Nepal.

The analysis uses Sentinel-5P methane observations and complementary environmental datasets to investigate temporal variability, spatial differences and relationships between methane concentrations and environmental conditions.

## Research Objectives

1. Quantify monthly, seasonal and annual atmospheric methane variability from 2019 to 2025.
2. Compare methane concentrations across urban, peri-urban and regional-background zones.
3. Examine relationships between methane concentrations and rainfall, temperature, soil moisture, wind, vegetation, inundation, cropland and built-up land.

## Study Area

Biratnagar Metropolitan City, Nepal, and its surrounding regional analytical domain.

The analysis considers three spatial zones:

- Urban core: Biratnagar boundary.
- Peri-urban: 0–10 km outside the urban boundary.
- Regional background: 25–50 km surrounding the urban boundary.

## Data Sources

| Dataset | Provider | Application |
|---|---|---|
| Sentinel-5P OFFL CH4 | Copernicus | Atmospheric methane |
| CHIRPS Daily | UCSB Climate Hazards Center | Rainfall |
| ERA5-Land Monthly | ECMWF | Temperature, soil moisture and wind |
| MODIS MOD13Q1 | NASA | Vegetation index |
| Sentinel-1 GRD | ESA/Copernicus | Surface water and inundation screening |
| ESA WorldCover 2021 | ESA | Cropland and built-up land |
| Global Surface Water | JRC | Permanent water masking |

## Repository Structure

```text
biratnagar-atmospheric-methane/
├── README.md
├── gee/
│   └── biratnagar_methane_analysis.js
├── data/
│   └── README.md
├── outputs/
│   └── README.md
├── docs/
│   └── methodology.md
├── .gitignore
└── LICENSE
```

## Software and Platform

- Google Earth Engine
- JavaScript API
- Git
- GitHub

## How to Run

1. Open the Google Earth Engine Code Editor.
2. Access the JavaScript file in the `gee` directory.
3. Copy the script into the Code Editor.
4. Ensure access to the required study-area asset.
5. Run the script.
6. Review the Console, Map layers and Tasks tab.
7. Start the required exports to Google Drive.

Google Earth Engine Code Editor: https://code.earthengine.google.com/

## Important Notes

The study-area asset is referenced using a Google Earth Engine project asset ID. Users must have appropriate permissions to access this asset or replace the asset ID with their own study-area asset.

The exported CSV tables and GeoTIFF files are generated through Google Earth Engine and are not included in this repository.

The methane analysis applies spatial and temporal filtering and physical-range screening. Refer to the source code and methodology documentation for implementation details.

## Reproducibility

The repository provides the analysis scripts and dataset references. Reproduction requires access to Google Earth Engine and the necessary input assets.

## Citation

If you use this code in research, please cite the repository and the associated research publication when available.

## Author

Charan Bhattarai

## License

See the LICENSE file for reuse conditions.
