This repository is intended to provide additional information for the paper "Carbon Farming practices assessment: Modelling spatial changes of Soil Organic Carbon in Flanders, Belgium" - By Spotorno, Gobin, Vanongeval, Del Borghi and Gallo - 2024 

(doi: https://doi.org/10.1016/j.scitotenv.2024.171267)

The repository is divided in two main sections:
1. Scripts
2. Files/Output

The scripts section contains the scripts related only to the RothC model: 
- One script for the Spin Up phase
- One script for the Warm Up phase
- Six different scripts for the Forward phase (one for each scenario studied: BAU, Cover Crops and four different improved crop rotations)

The Flies/Output section is divided into sub-folders containing the output shapefiles (and related shx, cpg, etc...) of each phase:

1.SPIN_UP

2.WARM_UP

3.BAU

4.CoverCrops

5.ImprovedCropRotation

Each FORWARD result is a point-shapefile, with three columns attribute table:
1. ESA_L_C -> refers to the Land Use Class of each point (2 means agricultural arable land)
2. SOC_t0 -> the Soil Organic Carbon content in 2022 (resulting from the Warm Up phase)
3. SOC_***_2 -> the Soil Organic Carbon content at 20 years for the specific scenario
