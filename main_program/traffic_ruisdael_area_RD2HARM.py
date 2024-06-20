"""
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

Creator: Dr. Arseni Doyennel
Contact Email: a.doyennel@vu.nl
"""

import glob
import os
import numpy as np
import time
import os.path
import shutil
import math
import fiona
import pandas as pd
import geopandas as gpd
from shapely.geometry import Polygon, MultiPolygon, LineString, MultiLineString
from pyproj import Proj, transform
import pyproj
from shapely.geometry import box
from fiona.crs import CRS

# Import required settings and libraries
from emission_preparation_setting import (spec_name, year, x0, y0, xn, yn,
                                          xres, yres, dx, dy,
                                          noxdir_script_traffic_new_python_only,
                                          distin_folder_ruisdael_area_csv2gpkg_new_python_only)




class AreaEmissionsReprojection_Traffic:
    def __init__(self):
        self.sourcedir = distin_folder_ruisdael_area_csv2gpkg_new_python_only
        self.noxdir = noxdir_script_traffic_new_python_only
        self.snaplist = ['7']  # always 7, because this script is for SNAP 7 only!
        self.spec_name = spec_name
        self.year = year
        self.x0 = x0
        self.y0 = y0
        self.xn = xn
        self.yn = yn
        self.xres = xres
        self.yres = yres
        self.dx = int(dx)
        self.dy = int(dy)
        self.proj_HARM = '+proj=lcc +lat_1=52.500000 +lat_2=52.500000 +lat_0=52.500000 +lon_0=.000000 +k_0=1.0 +x_0=649536.512574 +y_0=1032883.739533 +a=6371220.000000 +b=6371220.000000'
        self.road_network = self.read_nox_data() #we call read_nox_data(self) method right in the constructor
        self.emis_res = 1000 #resolution of area emission input

    def create_directories(self, rundir):
        if not os.path.exists(rundir):
            os.makedirs(rundir)
            print(f'The {rundir} directory is created!')
        else:
            print(f'The {rundir} directory exists!')
            for f in glob.glob(rundir + '*'):
                os.remove(f)

    #Method to read NOx traffic emission data required for refinement
    def read_nox_data(self):
        columns_to_read = ['Shape_Length', 'LVtot_N', 'MVtot_N', 'ZVtot_N']
        layer_name = f'OD2018_jr{self.year}'

        with fiona.open(os.path.join(self.noxdir, f'wegenbestand_{self.year}.gdb'), layer=layer_name) as src:
            filtered_features = [{'geometry': feature['geometry'],
                                  'properties': {key: feature['properties'][key] for key in columns_to_read}} for
                                 feature in src]

        road_network = gpd.GeoDataFrame.from_features(filtered_features)
        
        # Add a new column 'TVtot_N_sum' containing the sum of 'LVtot_N', 'MVtot_N', and 'ZVtot_N'
        road_network['TVtot_N_sum'] = road_network[['LVtot_N', 'MVtot_N', 'ZVtot_N']].sum(axis=1)
        
        # Multiply the TVtot_N_sum column by 1000 to convert from g/year to kg/year:
        road_network['TVtot_N_sum'] = road_network['TVtot_N_sum'] * 1000
        
        road_network.crs = 'EPSG:28992'
        road_network = road_network.to_crs(self.proj_HARM)

        return road_network


    def reproject_emission_and_refine(self):

        co2_column = self.spec_name + '_sum'

        if not self.xn % self.dx and not self.yn % self.dy:
            print('dx/dy is a multiple of dx_scale/dy_scale!')
        else:
            raise Exception("check dx and dy before proceed!")

        rundir = self.sourcedir + 'snap_7_' + spec_name + '_gpkg_HARM/'

        self.create_directories(rundir)

        inputdir_snap = self.sourcedir + 'snap_7_' + self.spec_name + '_gpkg_RD/'
        source = inputdir_snap + '' + self.spec_name + '_' + str(self.year) + '_SNAP_7_residual.gpkg'
        destination = rundir + '' + self.spec_name + '_' + str(self.year) + '_SNAP_7_residual.gpkg'

        if os.path.isfile(source):
            shutil.copy(source, destination)
            print('Copied', '' + self.spec_name + '_' + str(self.year) + '_SNAP_7_residual.gpkg')

        emis_gdf_full = gpd.read_file(rundir + f'{self.spec_name}_{self.year}_SNAP_7_residual.gpkg')

        xminlist = range(self.x0, self.x0 + int(self.xn * self.xres), int(self.dx * self.xres))
        yminlist = range(self.y0, self.y0 + int(self.yn * self.yres), int(self.dy * self.yres))

        for xmin in xminlist:
            for ymin in yminlist:

                xmax = int(xmin + self.dx * self.xres)
                ymax = int(ymin + self.dy * self.yres)

                extname = f'{xmin}_{xmax}_{ymin}_{ymax}'
                targetfile = rundir + 'HARM_snap_7_' + extname + '.gpkg'

                print(f'extname = {xmin}_{xmax}_{ymin}_{ymax} processing...')

                # Define grid parameters
                grid_xmin, grid_xmax = xmin, xmax
                grid_ymin, grid_ymax = ymin, ymax
                grid_xres, grid_yres = xres, yres

                # Generate x and y coordinates for grid cells
                x_coords = np.arange(grid_xmin, grid_xmax, grid_xres)
                y_coords = np.arange(grid_ymin, grid_ymax, grid_yres)

                # Initialize an empty list to store grid cell properties
                grid_cells = []

                # Initialize id
                id_counter = 1

                # Generate grid cells
                for x in x_coords:
                    for y in y_coords:
                        top = grid_ymax - y + grid_ymin
                        right = x + grid_xres
                        bottom = grid_ymax - (y + grid_yres) + grid_ymin

                        # Create a Polygon geometry
                        grid_cell = Polygon([(x, top), (right, top), (right, bottom), (x, bottom)])

                        # Create a dictionary to represent the properties
                        grid_properties = {
                            'id': id_counter,
                            'left': x,
                            'top': top,
                            'right': right,
                            'bottom': bottom,
                            'geometry': grid_cell
                        }

                        grid_cells.append(grid_properties)

                        # Increment the id counter
                        id_counter += 1

                # Create a GeoDataFrame from the list of grid cell properties
                grid_gdf = gpd.GeoDataFrame(grid_cells, crs=self.proj_HARM)
                grid_buff = grid_gdf.copy()      
                grid_buff['geometry'] = grid_buff['geometry'].buffer(-0.0000001, cap_style=1, join_style=2, mitre_limit=2)

                # Load road network data
                road_gdf_1 = self.road_network

                # Perform the intersection
                intersection_gdf = gpd.overlay(road_gdf_1, grid_gdf, how='intersection')

                # Rename the 'left' column in the intersection result to 'H_left'
                intersection_gdf = intersection_gdf.rename(columns={'left': 'H_left'})

                # Select only the desired columns
                selected_columns = ['Shape_Length', 'TVtot_N_sum', 'geometry']
                intersection_gdf = intersection_gdf[selected_columns]

                # Convert the geometry back to EPSG:28992
                intersection_gdf = intersection_gdf.to_crs('EPSG:28992')
                
                # Ensure that geometry is MultiLineString
                intersection_gdf['geometry'] = intersection_gdf['geometry'].apply(lambda x: x if x.geom_type == 'MultiLineString' else MultiLineString([x]))

                # Select relevant part of gridded 1x1km spec_name emissions
                emis_gdf = emis_gdf_full.copy()

                # Round up the extent coordinates to the nearest self.emis_res:
                xmin_rounded = np.floor(xmin / self.emis_res) * self.emis_res
                ymin_rounded = np.floor(ymin / self.emis_res) * self.emis_res
                xmax_rounded =  np.ceil(xmax / self.emis_res) * self.emis_res
                ymax_rounded =  np.ceil(ymax / self.emis_res) * self.emis_res

                # Define the bounding box extent using Shapely's box function
                bbox = box(xmin_rounded, ymin_rounded, xmax_rounded, ymax_rounded)

                # Create a GeoDataFrame from the extent
                extent_gdf = gpd.GeoDataFrame(geometry=[bbox], crs=CRS.from_string(self.proj_HARM))  

                # Reproject extent_gdf to the CRS of input_gpkg
                extent_gdf = extent_gdf.to_crs(emis_gdf.crs)

                # Clip the input GeoPackage file by the extent
                emis_gdf = gpd.overlay(emis_gdf, extent_gdf, how='intersection')
         

                if not emis_gdf.empty and not intersection_gdf.empty:

                    # Ensure that geometry is MultiLineString
                    intersection_gdf['geometry'] = intersection_gdf['geometry'].apply(
                        lambda x: x if x.geom_type == 'MultiLineString' else MultiLineString([x]))

                    # Save the selected emissions to a GeoPackage
                    grid_gdf.to_file(rundir + 'grid_' + extname + '.gpkg', driver='GPKG')

                    # Save the result to 'road_int_1'
                    intersection_gdf.to_file(rundir + 'road_int_1_' + extname + '.gpkg', driver='GPKG')

                    # Save the selected emissions to a GeoPackage
                    emis_gdf.to_file(rundir + 'emis_' + extname + '.gpkg', driver='GPKG')

                    # Perform the intersection between road_int_1 and emis

                    road_int_2 = gpd.overlay(
                        intersection_gdf[['Shape_Length', 'TVtot_N_sum', 'H_left', 'geometry']],
                        # Select relevant columns from 'road_int_1'
                        emis_gdf[['left', 'right', 'top', 'bottom', 'geometry']],
                        how='intersection',
                        keep_geom_type=True
                    )

                    weight = (road_int_2.geometry.length / road_int_2['Shape_Length']) * (
                                road_int_2['TVtot_N_sum'])

                    # Assign the weight to a new column
                    road_int_2['weight'] = weight

                    # Rename
                    road_int_2 = road_int_2.rename(columns={'left': 'RD_left'})

                    # Ensure that geometry is MultiLineString

                    if not road_int_2.is_empty.all():


                        # Save the final road_int_2 with specified columns
                        selected_columns = ['Shape_Length', 'TVtot_N_sum', 'RD_left', 'weight',
                                        'geometry']
                        road_int_2 = road_int_2[selected_columns]
                        road_int_2.to_file(rundir + 'road_int_2_' + extname + '.gpkg', driver='GPKG')

                        # Delete elements with weight == 0
                        road_int_3 = road_int_2[road_int_2['weight'] != 0]
                        road_int_3.to_file(rundir + 'road_int_3_' + extname + '.gpkg', driver='GPKG')

                        emis_gdf = None
                        # Load the emission and road intersection data
                        emis_gdf = gpd.read_file(f'{rundir}emis_{extname}.gpkg')

                        emis_gdf = emis_gdf.to_crs('EPSG:28992')  # Reproject to  CRS 'EPSG:28992'
                        road_gdf = road_int_3.to_crs('EPSG:28992')  # Reproject to  CRS 'EPSG:28992'

                        # Ensure emission geometry is Polygon
                        emis_gdf['geometry'] = emis_gdf['geometry'].apply(
                            lambda geom: geom if geom.geom_type == 'Polygon' else Polygon(geom))
                        
                        # Add a small buffer to emission grid polygons
                        emis_gdf['geometry'] = emis_gdf['geometry'].buffer(-0.0000001)  # Adjust the buffer size as needed


                        weight_gdf = emis_gdf

                        # Initialize the total_weight column
                        weight_gdf['weight_sum'] = 0.0

                        # Iterate through each emission grid polygon
                        for index, grid_polygon in weight_gdf.iterrows():
                            # Calculate the total weight for the current polygon by summing weights of intersecting road_intersections
                            total_weight = road_gdf[road_gdf.intersects(grid_polygon['geometry'])]['weight'].sum()
                            weight_gdf.at[index, 'weight_sum'] = total_weight

                        # Select the desired columns
                        selected_columns = ['id', 'left', 'top', 'right', 'bottom', self.spec_name, 'layer', 'path',
                                        'weight_sum', 'geometry']
                        weight_gdf = weight_gdf[selected_columns]

                        # Replace 0 with NaN in 'weight_sum' column
                        weight_gdf.loc[weight_gdf['weight_sum'] == 0, 'weight_sum'] = np.nan

                        weight_gdf.to_file(rundir + 'emis_weight_' + extname + '.gpkg', driver='GPKG')

                        # Convert the entire emis_weight and road_int_3 GeoDataFrames to the custom CRS
                        emis_weight = weight_gdf.to_crs('EPSG:28992')
                        road_int_3 = road_int_3.to_crs('EPSG:28992')

                        # Broadcast total weight back to all road segments together with emission data
                        intersections = gpd.overlay(road_int_3, emis_weight, how='intersection', keep_geom_type=False)

                        # Remove features with weight_sum = 0
                        intersections = intersections[intersections['weight_sum'] != 0]

                        # Remove features where spec_name (e.g., 'co2') is equal to 0
                        intersections = intersections[
                            intersections[self.spec_name] != 0]

                        # Calculate the spec_name emission per road segment
                        intersections[f'{self.spec_name}_sum'] = intersections.apply(
                            lambda row: row['weight'] / row['weight_sum'] * row[self.spec_name], axis=1)

                        # Save the result as a GeoPackage file:

                        intersections.to_file(rundir + f'road_{self.spec_name}_cleaned_' + extname + '.gpkg', driver='GPKG')


                        road_data = gpd.read_file(rundir + f'road_{self.spec_name}_cleaned_{extname}.gpkg')
                        road_data = road_data.to_crs(self.proj_HARM)  # Reproject to HARMONIE CRS

                        # Define the spatial join using GeoPandas
                        joined = gpd.sjoin(grid_buff, road_data, how='inner', predicate='intersects')

                        # Group by the grid ID and calculate the sum of spec_name_sum for each group
                        grouped = joined.groupby('id_left')[f'{self.spec_name}_sum'].sum().reset_index()

                        # Replace 0 with NaN
                        grouped[f'{self.spec_name}_sum'].replace(0, np.nan, inplace=True)

                        # Merge the grouped data back into the grid data
                        grid_buff = grid_buff.merge(grouped, left_on='id', right_on='id_left', how='left')

                        grid_buff.to_file(targetfile, driver='GPKG')

                print(f'extname = {xmin}_{xmax}_{ymin}_{ymax} done.')

                # Optionally, you can clear the GeoDataFrame if needed

                grid_gdf = None
                emis_gdf = None
                weight_gdf = None
                intersection_gdf = None
                road_gdf_1 = None
                road_int_1 = None
                road_int_2 = None
                road_int_3 = None
                intersections = None
                grid_buff = None
                road_data = None

                #raise Exception("stop upto now!")

        print(f'SNAP7 Traffic is done.')

if __name__ == "__main__":
    # Instantiate the AreaEmissionsReprojection class
    AER=AreaEmissionsReprojection_Traffic()
    # Run the reprojection and downscaling process
    AER.reproject_emission_and_refine()


