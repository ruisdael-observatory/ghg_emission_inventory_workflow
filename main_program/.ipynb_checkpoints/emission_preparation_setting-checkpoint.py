"""
Module Description:
The general setting file to prepare 
area and point source emission input for DALES model to simulate the Netherlands area.

This is the main setting file, common for all modules in workflow.

In addition, in some modules there are local settings to set up the emission processing.

Also, functions in emission_construction_functions can be updated, i.e., the list of subcategories are listed there.

Currently, workflow supports only RD coordinates, as original RIVM emissions are in RD..

Author: Dr. Arseni Doyennel
Place: VU Amsterdam

"""
import os
from GridDales import GridDales
import numpy as np

author= ''
email=''

snaplist = [1,2,3,4,5,7,8,9,10,11] #Note: snap 11 is used for aviation only for CO2!

year = 2022   #year for which emissions will be used from

year_start = 2022
month_start = 8
day_start = 23
hour_start = 0

year_end = 2022
month_end = 8
day_end = 24
hour_end = 0

# Note currently  workflow supports co2 and ch4 emissions only (other compounds have not been tested!)

#spec_name = 'ch4'
#chemcomp = 'Methaan'
spec_name = 'co2'
chemcomp = 'Koolstofdioxide'

domain='test_excel' #coarse is a mother domain

explicit_point=False  # If True: keep explicit point source emissions for DALES # If False: aggregate all point sources into area emissions
#Note: to use explicit_point you must have stack height, exist temperature, volume and area in your point emission file 
#(HOOGTE, LENGTE, BREEDTE, UITSTROOMOPENING_M2, TEMPERATUUR, VOLUMESTROOM) (usually do not available in data from ER portal, so contact RIVM for this )!

ER_excel_to_csv=True #True: If you have excel files directly from ER portal and need to convert them to csv to be able to use in this workflow (use only once!)
map_ER_excel = 'ER_DataExport-2026-08-19-011609.xlsx'
point_ER_excel= 'ER_DataExport-2026-08-19-010416.xlsx'

point_file_harm_p_only=f'ERI_{spec_name}_{year}_metemk_processed_RD.csv' #Automatically created by the program (name should be other from point_source_harm_file!)
point_source_harm_file=f'ERI_{spec_name}_{year}_metemk_new_RD_example.csv' 
point_ER_direct=True #True if direct ER file from https://data.emissieregistratie.nl/export is used in point_source_harm_file (converted to csv!)

mapped_data = f'ER_map_{year}_{spec_name}_example.csv'
map_ER_direct=True #True if direct ER file from https://data.emissieregistratie.nl/export is used in mapped_data (converted to csv!)

cbs_dataset_name = f'cbs_vk100_{year_start}.gpkg'  #name of your CBS file (household data)
#agro_file_name = "agro_grassland_25m.gpkg" #LGN data (https://lgn.nl/bestanden), here, only "agrarisch_gras", "natuurlijk_beheerde_graslanden" are 
#included. Original data raserized to 25m
agro_file_name = "grassland_processed_A.gpkg" #LGN data (https://lgn.nl/bestanden), here, only "agrarisch_gras" is inc. Original data raserized to 25m
vaarweg_file = "vaarwegvakken.gpkg"
industrial_file_name = "BAG_industrial_snap3.gpkg"
nox_multistring_dataset="NOx_traffic.gpkg" #NOx_traffic.gpkg created from RIVM NOx data from traffic length data
temp_prof_file='CAMS-REG-TEMPO_EUR_0.1x0.1_tmp_weights_v3.1_daily_2020.nc'
#################################################################################

#folders location in point_source_plume_processing:
datadir_point_source_plume_processing = os.path.join(
    os.path.dirname(os.getcwd()),
    'emission_raw_data',
    'emission_files/'
)

# Define the main target directory (you should replace the leading '...' with an actual path)
targetdir_main = os.path.join(
    os.path.dirname(os.getcwd()),  # This is the current directory
    'final_input',
    'ncdfs'
)

# Assuming you have defined variables: year, month_start, day_start, domain
subfolder_name = f'file_set_{year}_{month_start}_{day_start}_{domain}'
targetdir_point_source_plume_processing = os.path.join(
    targetdir_main,
    subfolder_name,
    'point_sources/'
)

rundir_ruisdael_area_residuals=os.path.join(os.getcwd(), 'residuals_inputs/') #current folder

#folders location in ruisdael_area_residuals:
datadir_ruisdael_area_residuals = datadir_point_source_plume_processing
rundir_ruisdael_area_residuals=os.path.join(os.getcwd(), 'residuals_inputs/') #current folder
cbs_loc_consu_preprocessing = os.path.join(os.path.dirname(os.getcwd()),'refinement_area_emiss/')

outputdir_main_ruisdael_area_csv2netc=os.path.join(targetdir_main, subfolder_name,'area/' ) #create a specific folder 
noxdir =  os.path.join(os.path.dirname(os.getcwd()),'refinement_area_emiss/')

#folders location in ruisdael_area_total_static: 
targetdir_ruisdael_area_total_static=os.path.join(outputdir_main_ruisdael_area_csv2netc,'3D_static_emiss/') 

#folders location in create_hourly_emissions_3D:
output_dir_create_hourly_emissions_3D=os.path.join(outputdir_main_ruisdael_area_csv2netc,'DALES_input_area_emissions_final/')

#################################################################################

#X0 and Y0 of the LES domain:
x0 = 57535.0
y0 = 411180.0

#unneeded, can 0,0 is used
x_offset, y_offset = 0,0

# Number of grid points LES
itot = 512
jtot = 512
ktot = 128

xn = itot 
yn = jtot

xres = 100  #from namoptions (experiment setting)
yres = 100  #from namoptions (experiment setting)

zmax=92 #approx. vertical top of initial emission profile in meters

xsize = int(xn*xres)
ysize = int(yn*yres)

# Number of x,y MPI processes (not used within the workflow yet)
nprocx = 16
nprocy = 16

proj_params = {
    'proj4': '+proj=sterea +lat_0=52.15616055555555 +lon_0=5.38763888888889 +k=0.9999079 +x_0=155000 +y_0=463000 +ellps=bessel +towgs84=565.417,50.3319,465.552,-0.398957,0.343988,-1.8774,4.0725 +units=m +no_defs'
}

# Define LES grid (target grid with 100x100m resolution)
grid_params = {
    'xsize': xsize,
    'ysize': ysize,
    'itot': itot,
    'jtot': jtot,
    'kmax': 128,
    'dz0': 20,  #if resolution is 400m we can make tihs also 100m or so?
    'alpha': 0.012,
    'southwest_x': 0,  # Placeholder for RD coordinates of the southwest corner
   'southwest_y': 0   # Placeholder for RD coordinates of the southwest corner
}

les_grid = GridDales(grid_params)

# Find the closest index (number of vertical layers for emissions)
layer_num_emiss_z = np.argmin(np.abs(les_grid.zt - zmax))+1 #IMPORTANT parameter: if SNAP category has vertical component, emiss will be placed at some number of model layers (depends on zmax) from the emission bottom to the emission top
                                    #NOTE: for point sources, if you use explicit input method,
                                    # plume rise height is calculated from stack height interactively; preliminary
                                    # the number of vertical layers, where emissions are allocated, is unknown and
                                    # the emiss will be distributed directly in the model code !

zmax=int(np.floor(les_grid.zt[layer_num_emiss_z-1])) #replace zmax with the value exactly from les_grid.zt (rounded down) for the top of emission profile

print(layer_num_emiss_z, les_grid.zt[layer_num_emiss_z-1])

#(TODO: this in principle allow running this application in parallel using MPI)

dx = xn/nprocx
dy = yn/nprocy