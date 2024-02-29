"""
Module Description:
The general setting file to prepare 
area and point source emission input for DALES model to simulate the Netherlands enviroment

This is the main setting file, common for all scripts needed to prepare emissions

In addition, in file ruisdael_area_residuals there is a setting, if you need to update some of SNAP categories
and large country domain settings

Also, functions in emission_construction_functions can be updated as the list of subcategories are listed there

Author: Dr. Arseni Doyennel
e-mail: a.doyennel@vu.nl

"""
import os


snaplist = [1,2,3,4,5,7,8,9,10]

year = 2017   #year for which emissions will be used from

year_start = 2018   #can be different from the "year", since here it is a simulation period setting
month_start = 6
day_start = 20
hour_start = 0

year_end = 2018
month_end = 6
day_end = 26
hour_end = 6


spec_name = 'co2'
chemcomp = 'Koolstofdioxide'


point_source_rd_file = f'ERI_{spec_name}_metemk.csv'
open_source_data_points=f'ER_bronnen_{spec_name}.csv'
private_data_points=f'ERI_{spec_name}_metemk_harm.csv'       
mapped_data = f'ER_map_{year}_{spec_name}.csv'

#################################################################################
#folders location in point_source_plume_processing:
datadir_point_source_plume_processing=os.path.join(os.path.dirname(os.getcwd()),'emission_raw_data/emission_files/') #one folder back
targetdir_point_source_plume_processing=os.path.join(os.path.dirname(os.getcwd()),'point_sources/point_sources_prepared/') 


#folders location in ruisdael_area_residuals:
datadir_ruisdael_area_residuals = datadir_point_source_plume_processing
rundir_ruisdael_area_residuals=os.path.join(os.getcwd(), 'residuals_inputs/') #current folder
cbs_loc_consu_preprocessing = os.path.join(os.path.dirname(os.getcwd()),'refinement_area_emiss/consu/')

#folders location in ruisdael_area_csv2gpkg:
distin_folder_ruisdael_area_csv2gpkg_new_python_only=os.path.join(os.path.dirname(os.getcwd()),'gpkg_files/') 

#folders location in traffic_ruisdael_area_RD2HARM:
noxdir_script_traffic_new_python_only =  os.path.join(os.path.dirname(os.getcwd()),'refinement_area_emiss/traffic/') 


#folders location in ruisdael_area_gpkg2csv_new_python_only:
outputdir_main_ruisdael_area_gpkg2csv_new_python_only=os.path.join(os.path.dirname(os.getcwd()),'csv_files/') 


#folders location in ruisdael_area_csv2netc:
outputdir_main_ruisdael_area_csv2netc=os.path.join(os.path.dirname(os.getcwd()),'ncdfs/file_set_1/') #create a specific folder for ncdfs and final output

#folders location in ruisdael_area_total_static: 
targetdir_ruisdael_area_total_static=os.path.join(outputdir_main_ruisdael_area_csv2netc,'final_hourly_3D_emiss/') 

#folders location in create_hourly_emissions_3D:
rundir_create_hourly_emissions_3D =os.path.join(os.path.dirname(os.path.dirname(targetdir_ruisdael_area_total_static)),'')
output_dir_create_hourly_emissions_3D =os.path.join(targetdir_ruisdael_area_total_static,'DALES_input_area_emissions_final/')


#################################################################################

#Domain specific parameters (in HARM coordinates):

x0 = 874468     # 
y0 = 926741     # 

# Number of grid points LES
itot = 1536
jtot = 768
ktot = 128

xn = itot #number of gridcells in x direction (itot)
yn = jtot #number of gridcells in y direction (jtot)

xres = 156.25  #from namoptions (experiment setting)
yres = 156.25  #from namoptions (experiment setting)

zmax=200 #the vertical height of emissions in meters (emission plume rise high for area emission) 
#(the rounded height of "kemis" vertical grid point from namoptions) why 200m? ->
#see https://acp.copernicus.org/articles/19/4541/2019/ fig. 2 this could be tuned, I guess

xsize = int(xn*xres)
ysize = int(yn*yres)

# Number of x,y MPI processes

nprocx = 16 # 256 cores = 2 nodes
nprocy = 16

#You can add dx and dy (this is only for codes starting from RD2HARMONIE. Before, both dx and dy either 1000 or 100, correspondingly): (TODO: this in principle allow running this application in parallel using MPI)

dx = xn/nprocx
dy = yn/nprocy

