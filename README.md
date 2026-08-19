# Emission Inventory Workflow

An application for preparing anthropogenic greenhouse-gas (GHG) emission inputs for the Dutch Atmospheric Large-Eddy Simulation (DALES) model.

> **Current status:** CO₂ and CH₄ emissions have been tested and are currently supported.
>
> **Coordinate system:** The workflow currently operates in Dutch RD coordinates (Rijksdriehoekscoördinaten), which are the original coordinates provided by the Dutch Emission Registration (Emissieregistratie) data.

---

## Overview

The **Emission Inventory Workflow** converts annual Dutch emission-inventory data into spatially and temporally resolved emission input suitable for DALES.

The workflow can:

- prepare explicit point-source emissions;
- prepare gridded surface/area emissions;
- subtract explicit point sources from gridded emissions where necessary;
- refine coarse area emissions using proxy/activity data;
- redistribute emissions from coarse grids to the LES grid;
- reproject emissions from RD coordinates to HARMONIE coordinates;
- increase the spatial resolution of area emissions;
- apply temporal disaggregation to hourly emissions;
- distribute emissions vertically where appropriate;
- generate final 3-D NetCDF emission input files for DALES.

The workflow is designed to keep the different processing stages modular. Individual modules can be run independently when required.

---

# Project Structure

The main project folders are organized as follows:


new_emission_prep_workflow/
│
├── main_program/
│   ├── run_script.ipynb
│   ├── emission_preparation_setting.py
│   ├── er_data_converter.py
│   ├── point_source_explicit_input_netcdf.py
│   ├── ruisdael_area_residuals.py
│   ├── ruisdael_area_total_static.py
│   ├── reproject_to_les.py
│   ├── reproject_to_les_snap_3.py
│   ├── reproject_to_les_snap_7.py
│   ├── reproject_to_les_snap_8.py
│   ├── reproject_to_les_snap_10.py
│   ├── create_hourly_emissions_3D.py
│   ├── emission_construction_functions.py
│   └── ...
│
├── emission_raw_data/
│   └── emission_files/
│       └── Raw emission data downloaded from
│           the Dutch Emission Registration
│
├── refinement_area_emiss/
│   └── Proxy and activity data used for
│       spatial refinement/disaggregation
│       of area emissions
│
└── final_input/
    └── Final emission input files prepared
        for DALES


### `main_program/`

Contains the main Python modules and the main workflow entry point.

### `emission_raw_data/emission_files/`

Contains the raw emission data downloaded from the Dutch Emission Registration.

### `refinement_area_emiss/`

Contains proxy and activity datasets used to spatially refine/disaggregate area emissions.

These datasets include information related to:

* households; (not included)
* industrial buildings; (not included)
* road-traffic NOₓ emissions; (not included)
* CAMS TEMPO REG for temporal disaggregation (not included)
* inland waterways; (included)
* agricultural land (50x50m2) (included)

IMPORTANT:
Because of the size limits, not all proxy and activity datasets are added and they need to be downloaded from other resources and placed to refinement_area_emiss/:

Household statistical files for any year can be downloaded from:  
https://www.cbs.nl/nl-nl/dossier/nederland-regionaal/geografische-data/kaart-van-100-meter-bij-100-meter-met-statistieken

industrial buildings data set is originally from: 
https://service.pdok.nl/lv/bag/atom/bag.xml
(use BAG (EPSG:28992) Geopackage)

After, need to extract data to prepare proxy dataset: "BAG_industrial_snap3.gpkg"
Use prepare_BAG_industrial_snap3.py in refinement_area_emiss/ to prepare "BAG_industrial_snap3.gpkg"

CAMS TEMPO REG for temporal disaggregation can be directly downloaded from https://eccad.sedoo.fr/#/metadata/506

road-traffic NOₓ emissions for the Nederlands can be requested from here:  https://doi.org/10.5281/zenodo.14961517

Other proxy/activity data are placed to refinement_area_emiss/
inland waterways from https://www.nationaalgeoregister.nl/geonetwork/srv/api/records/701d4eb8-8aae-4708-bba5-3edf6987676d/formatters/xsl-view?view=advanced&portalLink=

agricultural land generated from https://api.pdok.nl/rvo/gewaspercelen/ogc/v1/collections/brpgewas/items (only "Bouwland", "Grasland")

### `final_input/`

Contains the final emission input files generated for DALES, including hourly 3-D NetCDF files.

---

# Input Data

The workflow requires annual emission data from the Dutch Emission Registration (Emissieregistratie).

Open-access emission data can be downloaded from:

[https://data.emissieregistratie.nl/export](https://data.emissieregistratie.nl/export)

Two main types of emission data are used:

1. Gridded emission data
2. Company/point-source emission data

---

## 1. Gridded Emission Data

Gridded emissions are used as the basis for area-source emissions.

Recommended settings in the Emission Registration export interface are:


Compartiment:
    Lucht

Stof:
    Koolstofdioxide
    or
    Methaan

Gebiedsindeling:
    Combinatie kaart vk5 NCP en vk1 Nederland

Bronniveau:
    Subsector

Jaren:
    Required year(s)


The resulting data can contain emissions on a **1 × 1 km** or **5 × 5 km** grid, depending on the selected export.

---

## 2. Company / Point-Source Emission Data

Company-level emissions are used to prepare explicit point sources when sufficient stack/source parameters are available.

Recommended settings in the Emission Registration export interface are:


Compartiment:
    Lucht

Stof:
    Koolstofdioxide
    or
    Methaan

Gebiedsindeling:
    Nederland

Bronniveau:
    Bedrijf

Jaren:
    Required year(s)


### Explicit Point-Source Requirement

Explicit point-source processing requires the following source parameters to be present in the input data:


HOOGTE
UITSTROOMOPENING_M2
TEMPERATUUR
VOLUMESTROOM


If these parameters are not available, the explicit point-source preparation is skipped.

---

# Running the Workflow

The main workflow can be started using:


main_program/run_script.ipynb


The main configuration is located in:


main_program/emission_preparation_setting.py


This file contains the paths, simulation settings, emission species, dates, model domains, and other parameters required by the workflow.

Before running the workflow, check and update:


main_program/emission_preparation_setting.py


The workflow is modular, so individual processing modules can also be executed separately when required.

# Core Modules

## 1. `er_data_converter.py`

### Description

Converts Dutch Emission Registration Excel exports into CSV files used by the workflow.

The converter reads the `Emissies` worksheet and converts the emission column from the Dutch decimal-comma format to a decimal-point format.

This module is useful because reading large ER Excel files directly can be considerably slower than reading the converted CSV files.

The converter can prepare both:

* point-source emission CSV files;
* gridded emission CSV files.

If the target CSV file already exists, the converter can skip the conversion to avoid unnecessary processing of large Excel files.

---

# 2. `point_source_explicit_input_netcdf.py`

## Description

Prepares explicit point-source emissions from Dutch Emission Registration company-level data and generates hourly point-source input for DALES.

The module performs, where applicable:

* emission-category identification;
* gap filling of missing point-source information;
* assignment of source coordinates;
* processing of stack/source parameters;
* temporal processing;
* preparation of explicit point-source NetCDF input.

## Input

Company-level point-source emission data.

For explicit point-source processing, the following parameters are required:


HOOGTE
UITSTROOMOPENING_M2
TEMPERATUUR
VOLUMESTROOM


If these parameters are not available, explicit point-source processing is skipped.

## Output

Hourly NetCDF files containing explicit point-source emissions.

Example:


pointsources.{yyyy}{mm}{dd}{hh}{mm}.{tracer}.nc


The module supports processing over a simulation period rather than only a single time step.

---

# 3. `ruisdael_area_residuals.py`

## Description

Prepares residual area emissions by subtracting explicit point-source emissions from gridded emissions.

This step is required when point-source emissions are already included in the gridded emission inventory and are subsequently treated as explicit point sources.

Subtracting the point-source contribution prevents **double counting** of emissions.

The module also performs initial area-emission processing and refinement for selected SNAP categories.

## Input

* 1 × 1 km or 5 × 5 km gridded emission data;
* point-source emission data.

## Output

Residual area-emission CSV files separated by SNAP category.

Example:


co2_2022_SNAP_1_residual.csv
co2_2022_SNAP_2_residual.csv
co2_2022_SNAP_3_residual.csv
...


For SNAP 2, higher-resolution proxy/activity data can be used where available.

---

# Area-Emission Spatial Refinement

The workflow uses proxy and activity data to redistribute coarse emission data onto a higher-resolution LES grid.

Different SNAP categories use different refinement approaches depending on the available proxy/activity data.

The relevant proxy and activity data are stored in:


refinement_area_emiss/


---

# 4. `reproject_to_les.py`

## Description

Reprojects and redistributes coarse gridded emission data onto the LES grid.

For categories where suitable proxy/activity data are not available, spatial interpolation/redistribution is used.

---

# 5. `reproject_to_les_snap_3.py`

## Description

Refines SNAP 3 emissions using industrial-related building/activity information.

Proxy data are stored in:


refinement_area_emiss/


The proxy data are based on industrial-related buildings from the Dutch **BAG (Basisregistratie Adressen en Gebouwen)** database.

---

# 6. `reproject_to_les_snap_7.py`

## Description

Refines SNAP 7 emissions using road-traffic proxy/activity data.

The refinement uses information such as:

* road network data;
* road-traffic NOₓ emissions.

Proxy data are stored in:


refinement_area_emiss/


---

# 7. `reproject_to_les_snap_8.py`

## Description

Refines SNAP 8 emissions using inland-waterway/shipping-related spatial information.

Proxy data are stored in:

refinement_area_emiss/


---

# 8. `reproject_to_les_snap_10.py`

## Description

Refines SNAP 10 emissions using agricultural land/activity information, including:

* arable land;
* grassland;
* agricultural fields.

Proxy data are stored in:

refinement_area_emiss/


# 9. `ruisdael_area_total_static.py`

## Description

Generates static annual total area emissions for a selected model domain.

The module can also incorporate point-source emissions into the gridded emission field when point sources are **not** treated as explicit point sources.

Depending on the available information, point-source emissions can be distributed according to:

* company location;
* source coordinates;
* stack height, where available.

The module also prepares 3-D emission fields by applying vertical emission distributions for relevant SNAP categories.

---

# Vertical Distribution

Vertical distributions are assigned based on typical stack heights, plume heights, plume-rise information, and emission-inventory literature.

The current parameterization includes the following categories.

---

## SNAP 1 — Power Generation

Typical stack heights are approximately **70–150 m**.

A representative plume height of approximately **80 m** with a vertical spread of approximately **20 m** is used.

Reference:

**EMEP/EEA Air Pollutant Emission Inventory Guidebook 2019**, Large Combustion Plants.

---

## SNAP 3 — Industrial Combustion

Typical stack heights are approximately **30–80 m**.

A representative plume height of approximately **60 m** with a vertical spread of approximately **15 m** is used.

Reference:

**Dutch National Inventory Report**, including emission factors and stack information.

---

## SNAP 4 — Industrial Processes

Industrial process sources are generally assumed to have lower effective release heights than large combustion sources.

A representative plume height of approximately **60 m** with a vertical spread of approximately **15 m** is used.

---

## SNAP 8 — Mobile Sources / Shipping

Ship exhaust stacks commonly release emissions at approximately **20–40 m**.

A representative plume height of approximately **30 m** with a vertical spread of approximately **10 m** is used.

---

## SNAP 9 — Waste

Waste-related sources such as incinerators and landfill-gas flares are represented using a representative plume height of approximately **20 m** with a vertical spread of approximately **15 m**.

---

## Ground-Level Categories

SNAP categories without an explicitly parameterized vertical component are currently assigned to the first model layer:


SNAP 2
SNAP 5
SNAP 7
SNAP 10

This represents ground-level or diffuse emissions.

> These vertical profiles apply to gridded/mapped emissions for which explicit stack and plume information is unavailable. Explicit point sources with available stack parameters are treated separately.

---

# 10. `create_hourly_emissions_3D.py`

## Description

Creates the final hourly 3-D area-emission input for DALES.

Temporal disaggregation accounts for variations in emission strength associated with:

* month of year;
* day of week;
* hour of day.

Temporal profiles are based on CAMS-REG-TEMPO and literature data, including:

> Denier van der Gon, H. A. C. et al. (2011).

The resulting files contain hourly 3-D emissions for the selected simulation domain.

Example:


{tracer}_emis_{yyyy}{mm}{dd}{hh}{mm}_3d.nc


The final files are date-specific and cover the complete simulation domain.

They are independent of DALES processor decomposition and starting-time settings.

---

# Supporting Module

## `emission_construction_functions.py`

Contains shared functions used by multiple modules throughout the workflow.

These functions provide common operations for:

* emission processing;
* SNAP classification;
* spatial processing;
* temporal processing;
* coordinate transformations;
* emission preparation.

---

# Coordinate Systems

The original Dutch Emission Registration data are provided in Dutch RD coordinates.

The workflow currently performs the initial emission preparation in **RD coordinates**.

For LES applications, emissions can subsequently be reprojected to the **HARMONIE coordinate system** and redistributed onto the required LES grid.

> **Current limitation:** The workflow has currently been tested using RD-coordinate emission data. Other coordinate systems may require additional testing and configuration.

---

# Emission Species

Currently tested and supported:

* **CO₂**
* **CH₄**

Other emission species may require additional testing and/or modifications to the workflow.

> **Important:** Do not assume that the workflow is validated for other species without testing the complete processing chain.

---

# Output

The final workflow produces emission input suitable for DALES, including:

## Explicit Point Sources

Hourly NetCDF files containing explicit point-source emissions.

Example:


pointsources.{yyyy}{mm}{dd}{hh}{mm}.{tracer}.nc

## Area Emissions

Hourly 3-D NetCDF files containing spatially distributed area emissions.

Example:

{tracer}_emis_{yyyy}{mm}{dd}{hh}{mm}_3d.nc

Final files are stored in:

final_input/

---
# Quick Start

## 1. Download Emission Data

Download the required emission data from the Dutch Emission Registration:

[https://data.emissieregistratie.nl/export](https://data.emissieregistratie.nl/export)

Select the required:

* emission species;
* year;
* spatial resolution;
* source level.

See the [Input Data](#input-data) section for recommended export settings.

---

## 2. Store Raw Data

Place the downloaded ER files in:

emission_raw_data/emission_files/

For example:

emission_raw_data/
└── emission_files/
    ├── ER_DataExport-YYYY-MM-DD-XXXXXX.xlsx
    └── ...

---

## 3. Configure the Workflow

Open:

main_program/emission_preparation_setting.py

Check and configure:

* input directories;
* output directories;
* emission species;
* simulation year;
* simulation period;
* model domain;
* grid resolution;
* emission files;
* point-source options;
* area-emission options.

---

## 4. Convert ER Excel Data

If required, use:

main_program/er_data_converter.py

to convert the ER Excel files to CSV.

This can significantly improve reading performance for large emission datasets.

The converter creates semicolon-separated CSV files compatible with the rest of the workflow.

---

## 5. Run the Main Workflow

Start:

main_program/run_script.ipynb

The notebook provides the main entry point for the complete workflow.

---

## 6. Retrieve Final Input

The final emission files are written to the configured output directory, including:

final_input/

The resulting NetCDF files can then be used as input for DALES simulations.

---

# Typical Processing Sequence

A typical complete processing sequence is:


1. Download ER data
        │
        ▼
2. Store raw Excel files
   emission_raw_data/emission_files/
        │
        ▼
3. Convert Excel → CSV
   er_data_converter.py
        │
        ├─────────────────────────┐
        │                         │
        ▼                         ▼
4. Point-source data       5. Gridded data
        │                         │
        ▼                         ▼
6. Explicit point-source   7. Point-source subtraction
   preparation                 from area emissions
        │                         │
        │                         ▼
        │                  8. Area-emission refinement
        │                         │
        │                         ▼
        │                  9. Reprojection /
        │                     spatial redistribution
        │                         │
        │                         ▼
        │                  10. Vertical distribution
        │                         │
        └────────────┬────────────┘
                     │
                     ▼
             11. Temporal
                 disaggregation
                     │
                     ▼
             12. Final hourly
                 3-D NetCDF
                     │
                     ▼
                final_input/
                

Not every simulation requires every step.

For example:

* explicit point sources may be disabled;
* some SNAP categories may use direct spatial interpolation;
* other categories may use proxy-based spatial refinement;
* vertical distribution is only applied where required;
* temporal disaggregation can be configured according to the simulation requirements.

---

# Important Notes and Limitations

* **CO₂ and CH₄** are currently the tested emission species.
* The workflow currently uses **Dutch RD coordinates** for the initial emission preparation.
* Explicit point-source processing requires:

  * `HOOGTE`
  * `UITSTROOMOPENING_M2`
  * `TEMPERATUUR`
  * `VOLUMESTROOM`
* If the required explicit point-source parameters are unavailable, explicit point-source processing is skipped.
* The workflow assumes that point-source emissions may already be included in the gridded inventory. When explicit point sources are used, the corresponding emissions can therefore be subtracted from the gridded inventory to avoid double counting.
* Vertical emission profiles are parameterizations and should be considered assumptions where explicit stack/plume information is unavailable.
* Proxy-based spatial refinement depends on the availability and quality of the corresponding activity/proxy datasets.
* The quality of the final high-resolution emission field depends on both the original emission inventory and the spatial proxy data.
* Users should validate the resulting emissions for their specific:

  * simulation domain;
  * year;
  * emission species;
  * spatial resolution;
  * temporal period;
  * DALES configuration.
* The workflow is under active development, and individual modules may be updated independently.

---

# References

The workflow uses information and parameterizations from Dutch and European emission-inventory sources, including:

### Dutch Emission Registration

Dutch National Emission Registration data:

[https://data.emissieregistratie.nl/export](https://data.emissieregistratie.nl/export)

### EMEP/EEA Air Pollutant Emission Inventory Guidebook

EMEP/EEA Air Pollutant Emission Inventory Guidebook 2019.

[https://www.eea.europa.eu/publications/emep-eea-guidebook-2019](https://www.eea.europa.eu/publications/emep-eea-guidebook-2019)

### Dutch National Inventory Report

Dutch National Inventory Report and associated information on emission factors and source characteristics:

[https://unfccc.int/](https://unfccc.int/)

### Temporal Profiles

Temporal disaggregation uses CAMS-REG-TEMPO data and literature-based temporal profiles, including:

Denier van der Gon, H. A. C. et al. (2011).

---

# Contributing / Feedback

The workflow is under ongoing development.

Suggestions, bug reports, improvements, and contributions are welcome.

If you use the workflow and identify issues or have suggestions for improving the methodology, please get in touch.

---

# Creator

**Dr. Arseni Doyennel**

Vrije Universiteit Amsterdam

### E-mail

* [a.doyennel@vu.nl](mailto:a.doyennel@vu.nl)
* [arsenidoyennel@gmail.com](mailto:arsenidoyennel@gmail.com)