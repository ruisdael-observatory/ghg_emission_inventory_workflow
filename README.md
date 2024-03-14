"Emission Inventory Workflow" is an application designed to prepare area and point source GHG emission inputs for high-resolution models that use the HARMONIE coordinates, such as DALES.

The workflow is used to create explicit point source input, residual area emissions, reproject and recalculate area emissions from RD to HARMONIE coordinates, and from low to high resolution. It refines the area emissions using proxy data (household data, traffic NOx emissions).

This workflow also includes temporal desegregation to get emissions per hour, accounting for month-to-month, data-to-day, hour-to-hour changes. Also, vertical distribution of area emissions is added for some snap categories.  

The input files for this workflow consist of annual raw point source emission data and 1x1 km (or 5x5 km) emission CSV files. Open-access emission data can be downloaded from https://data-preview.emissieregistratie.nl/.

Output: ncdf files of point source emission and 3-D area emission model input.

To run this application use run script: emiss_prep_workflow_run_script.ipynb 

File emission_preparation_setting.py is the setting script, to specify links, experiment and domain setting required in emission preparation.

The workflow modules:

1. point_source_explicit_input_netcdf.py

Module description:
Gapfilling of Emissieregistratie point sources
and creation of explicit point source input for DALES 

INPUT  (CSV) Emissiregistratie raw point source data

This file includes columns like DATASET,
EMISSIEJAAR, STOF, CODE_EMISSIEOORZAAK, EMISSIEOORZAAK, CODE_BEDRIJF,
NAAM_BEDRIJF, XCO_BEDRIJF, YCO_BEDRIJF, EMISSIEPUNT_CODE,
EMISSIEPUNT_NAAM, XCO_EMISSIEPUNT, YCO_EMISSIEPUNT,
TYPE, HOOGTE, LENGTE, BREEDTE, HOEK, UITSTROOMOPENING_M2,
TEMPERATUUR, UITTREEDSNELHEID, VOLUMESTROOM, WARMTEINHOUD,
EMISSIE, EENHEID, XCO_EMISSIEPUNT_HARM, YCO_EMISSIEPUNT_HARM,
XCO_BEDRIJF_HARM, YCO_BEDRIJF_HARM, delimited by a comma (,).

OUTPUT (netCDF) Gapfilled point source input ncdf files per processor block and per hour
Code updated with a loop over the simulation period to prepare point sources for the whole period



2. ruisdael_area_residuals.py

Module Description:
The main script to subtract point source emissions
from the area emissions to create ‘residual’ area emissions.
This data is then saved into appropriate CSV files
for subsequent downscaling and realignment to the HARM grid.
The refinement of area emissions data for the SNAP 2 (consumption) category
includes the incorporation of household data.

Input Data:
The script requires raw emission data provided by RIVM,
which comprises area emissions maps (1x1 and 5x5km)
and point sources obtained from open-source and privately provided data.

Example of Files:

open_source_data_points: This file includes columns such as
DATASET, NIC, BEDRIJF, STRAAT, HUISNUMMER, TOEVOEGING, POSTCODE,
WOONPLAATS, XCOORD, YCOORD, SBI_CODE, SBI, DOELGROEP, SUBDOELGROEP,
CODE_EMISSIEOORZAAK, EMISSIEOORZAAK, COMPARTIMENT, STOFCODE,
STOF, EENHEID, 2015, 2017, 2018, delimited by a semicolon (;).

point_source_rd_file: This file includes columns like DATASET,
EMISSIEJAAR, STOF, CODE_EMISSIEOORZAAK, EMISSIEOORZAAK, CODE_BEDRIJF,
NAAM_BEDRIJF, XCO_BEDRIJF, YCO_BEDRIJF, EMISSIEPUNT_CODE,
EMISSIEPUNT_NAAM, XCO_EMISSIEPUNT, YCO_EMISSIEPUNT,
TYPE, HOOGTE, LENGTE, BREEDTE, HOEK, UITSTROOMOPENING_M2,
TEMPERATUUR, UITTREEDSNELHEID, VOLUMESTROOM, WARMTEINHOUD,
EMISSIE, EENHEID, delimited by a comma (,).

(in general, private_data_points might include all open_source_data_points, but check this)

mapped_data: This file contains columns such as DATASET,
GEBIEDSINDELING, CODE_GEBIED, GEBIED, DOELGROEP, SUBDOELGROEP,
COMPARTIMENT, STOFCODE, STOF, EENHEID, EMISSIE_KG, delimited by a semicolon (;).

Output:
This script generates CSV files with residual area emissions,
segregated per SNAP category, e.g., co2_2017_SNAP_1_residual.csv.
Notably, for SNAP2, the resolution is 100m, while for other categories, it's 1000m.



3. ruisdael_area_csv2gpkg.py

Module Description:
This script converts area emissions data (all SNAP categories except 2)
from CSV-files into GeoFrame GeoPackage (gpkg) files for subsequent downsizing
and reassignment to the HARMONIE grid.

The two reasons for this are that this format allows for the exact definition of spatial extent
and to select a subdomain of the Netherlands that is relevant for a simulation.
Note: Be sure that this subdomain contains the simulated domain (in HARMONIE CRS).

It operates independently of the QGIS application, utilizing standard Python functions
as replacements for QGIS-related operations.



4. consu_ruisdael_area_csv2gpkg.py

Module Description:
This script converts area emissions data from CSV-files into GeoFrame GeoPackage (gpkg) files,
exclusively for SNAP 2 (consumption category), used for subsequent downscaling and coordinate
reassignment onto the HARMONIE grid.

The two reasons for this are that this format allows for the exact definition of spatial extent
and to select a subdomain of the Netherlands that is relevant for a simulation.
Note: Be sure that this subdomain contains the simulated domain (in HARMONIE CRS).

It operates independently of the QGIS application, leveraging standard Python functions
to replace QGIS-related operations.

Each script of this workflow is an independent module that can be run separately.



5. ruisdael_area_RD2HARM.py

Module Description:
This module serves for reprojecting area emissions from RD coordinates and
assigning them to the new HARMONIE coordinates.

The program intersects the Rijksdriehoek gridcells with defined HARMONIE coordinates
and reassigns emissions to the new grid. This method is exact and scaling independent.

However, this is the most demanding operation in the workflow! To prevent running out of RAM,
this operation is done in small blocks of the simulation domain and defined per SNAP category.

TODO: Possible upgrading with implementation of parallel programming for this script
(TODO: Implementing Python MPI) can be done for xmin in xminlist: and for ymin in yminlist loops.

This script no longer relies on QGIS functionalities but utilises standard Python functions.



6. cousu_ruisdael_area_RD2HARM.py

This module serves for reprojecting area emissions from RD coordinates and
assigning them to the new HARMONIE coordinates.

It is used only for SNAP 2 category (consummen).
The program intersects the Rijksdriehoek gridcells with defined HARMONIE coordinates
and reassigns emissions to the new grid.

This method is exact and scaling independent. However, this is the most demanding operation in the workflow!
To prevent running out of RAM, this operation is done in small blocks of the simulation domain.

TODO: Possible upgrading with implementation of parallel programming for this script
(TODO: Implementing Python MPI) can be done for xmin in xminlist: and for ymin in yminlist loops.

This script no longer relies on QGIS functionalities but utilizes standard Python functions.




7. traffic_ruisdael_area_RD2HARM.py

Module Description:
This module serves for reprojecting area emissions from RD coordinates and
assigning them to the new HARMONIE coordinates.

It is used only for SNAP 7 category (traffic).
The program intersects the Rijksdriehoek gridcells with defined HARMONIE coordinates
and reassigns emissions to the new grid.

IMPORTANT: To make this script work, you must have traffic activity and NOx emissions
calculated by Dat.mobility, e.g., high-resolution traffic intensity on road level
('Shape_Length', 'LVtot_N', 'MVtot_N', 'ZVtot_N') for spatially aggregation of traffic emissions.
These data can be obtained from RIVM. If you have no these data, use standard ruisdael_area_RD2HARM for SNAP7.

This method is exact and scaling independent. However, this is the most demanding operation in the workflow
(+ included refinement)! To prevent running out of RAM, this operation is done in small blocks of the simulation domain.

TODO: Possible upgrading with implementation of parallel programming for this script
(TODO: Implementing Python MPI) can be done for xmin in xminlist: and for ymin in yminlist loops.

This script no longer relies on QGIS functionalities but utilizes standard Python functions.




8. ruisdael_area_gpkg2csv.py

Module Description:
This script serves the purpose of translating GeoPackage (GPKG) files,
which have been early reprojected onto the HARMONIE grid emissions, into CSV-files.
This translation enables the subsequent retranslation of the reprojected emissions back to NetCDF format.
Notably, this script operates independently of the QGIS application,
as all QGIS-related functionalities have been meticulously substituted with analogous functions native to Python.




9. ruisdael_area_csv2netc.py

Module Description:
This script serves the purpose of translating reprojected
CSV files with area emissions into NetCDF-files per SNAP category.



10. ruisdael_area_total_static.py

Module Description:

This script is designed to generate static annual total area emissions within a specified domain.
Additionally, it incorporates the following functionalities:

Refinement of area emissions using point source data (if point sources are not utilized as explicit input).
Generation of 3-D emission files per SNAP category to facilitate the final preparation of area emission inputs.
( area emission are devided by number of vertical layers, covered with plume until the plume rise height layer to get fractions of emission!)
On vertical component of area emissions: currently, the vertical distribution of area emissions 
is based on results presented in Brunner et al. (2019) [https://acp.copernicus.org/articles/19/4541/2019/].


11. create_hourly_emissions_3D.py

This script is designed to create a final area emission input for DALES at 
the selected simulation domain.
Temporal aggregation is applied to account for
changes in the emission strength of a source
 for a certain hour in a day, day in the week, or month of the day.

Final input files are rather “DALES-unlike” as they are date specific and
cover the full simulation domain.
They are independent of processors settings or the starting time of the simulation.
The files are labeled as e.g. {tracer}_emis_{yyyy}{mm}{dd}{hh}{mm}_3d.nc

###############################################################
In addition, the emission_construction_functions.py file that collects functions used in workflow

All modules are independent and can be run separately.

If you have any suggestions to improve/update this workflow, do not hesitate to write me.

Creator:
Dr. Arseni Doyennel
E-mail: a.doyennel@vu.nl
