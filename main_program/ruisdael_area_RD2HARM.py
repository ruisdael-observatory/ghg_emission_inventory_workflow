"""
Module Description:
This module serves for reprojecting area emissions from RD coordinates and
assigning them to the new HARMONIE coordinates.

The program intersects the Rijksdriehoek gridcells with defined HARMONIE coordinates
and reassigns emissions to the new grid. This method is exact and scaling independent.

However, this is the most demanding operation in the workflow! To prevent running out of RAM,
this operation is done in small blocks of the simulation domain and defined per SNAP category.

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
from shapely.geometry import Polygon, MultiPolygon

# Import required settings and libraries
from emission_preparation_setting import (spec_name, year, x0, y0, xn, yn,
                                          xres, yres, snaplist, dx, dy,
                                          distin_folder_ruisdael_area_csv2gpkg_new_python_only)

class AreaEmissionsReprojection:
    def __init__(self):
        self.sourcedir=distin_folder_ruisdael_area_csv2gpkg_new_python_only
        self.snaplist = [x for x in snaplist if x not in (2, 7)]  # Exclude 2 and 7
        self.spec_name = spec_name
        self.year = year
        self.dx = int(dx) # must be an equal divisor of xn (nprox)
        self.dy = int(dy) # must be an equal divisor of yn (nproy)
        self.xn = xn
        self.yn = yn
        self.x0 = x0
        self.y0 = y0
        self.xres = xres
        self.yres = yres
        self.proj_HARM = '+proj=lcc +lat_1=52.500000 +lat_2=52.500000 +lat_0=52.500000 +lon_0=.000000 +k_0=1.0 +x_0=649536.512574 +y_0=1032883.739533 +a=6371220.000000 +b=6371220.000000' # Define the custom Proj4 string for HARMONIE CRS
        self.emiss_sc_factor = 1e6 #1e6 because input is in 1000x1000m in all SNAPs used here!
    def create_directories(self, rundir):
        if not os.path.exists(rundir):
            os.makedirs(rundir)
            print(f'The {rundir} directory is created!')
        else:
            print(f'The {rundir} directory exists!')
            for f in glob.glob(rundir + '*'):
                os.remove(f)

    def reproject_emissions(self):

        if not self.xn % self.dx and not self.yn % self.dy:
            print('dx/dy is a multiple of dx_scale/dy_scale!')
        else:
            raise Exception("check dx and dy before proceed!")

        for snap in self.snaplist:
            rundir = self.sourcedir + 'snap_' + str(snap) + '_' + self.spec_name + '_gpkg_HARM/'
            self.create_directories(rundir)

            inputdir_snap = self.sourcedir + 'snap_' + str(snap) + '_' + self.spec_name + '_gpkg_RD/'

            source = inputdir_snap + self.spec_name + '_' + str(self.year) + '_SNAP_' + str(snap) + '_residual.gpkg'
            destination = rundir + self.spec_name + '_' + str(self.year) + '_SNAP_' + str(snap) + '_residual.gpkg'

            if os.path.isfile(source):
                shutil.copy(source, destination)
                print(f'Copied: {self.spec_name}_{self.year}_SNAP_{snap}_residual.gpkg')

            xminlist = range(self.x0, self.x0 + int(self.xn * self.xres), int(self.dx * self.xres))
            yminlist = range(self.y0, self.y0 + int(self.yn * self.yres), int(self.dy * self.yres))

            for xmin in xminlist:
                for ymin in yminlist:
                    xmax = int(xmin + self.dx * self.xres)
                    ymax = int(ymin + self.dy * self.yres)
                    extname = f'{xmin}_{xmax}_{ymin}_{ymax}'
                    print(f'extname = {xmin}_{xmax}_{ymin}_{ymax} is processing...')
                    target_str = os.path.join(rundir + 'HARM_snap_' + str(snap) + '_' + extname + '.gpkg')

                    # Adjust the grid extent and resolution
                    # Generate grid cells using NumPy
                    x_coords = np.arange(xmin, xmax, self.xres)
                    y_coords = np.arange(ymin, ymax, self.yres)
                    grid_cells = [
                        Polygon([(x, y), (x + self.xres, y), (x + self.xres, y + self.yres), (x, y + self.yres)]) for x
                        in x_coords for y in y_coords]
                    # Create a GeoDataFrame from the grid cell polygons
                    grid_gdf = gpd.GeoDataFrame({'geometry': grid_cells}, crs=self.proj_HARM)

                    # Load emissions data (EPSG:28992)
                    input_gdf = gpd.read_file(
                        destination,
                        crs='EPSG:28992')


                    # Convert emissions data to HARMONIE CRS:
                    input_gdf = input_gdf.to_crs(self.proj_HARM)
                    # Add a small buffer to emission grid polygons:
                    input_gdf['geometry'] = input_gdf['geometry'].buffer(-0.0000001)    #Adjust the buffer size as needed

                    # Initialize an empty GeoDataFrame to store the intersection results:
                    intersection_gdf = gpd.GeoDataFrame(columns=['left', 'bottom', self.spec_name, 'emis', 'geometry'], crs=self.proj_HARM)

                    # Initialize an empty GeoDataFrame to store the intersection results:
                    intersection = gpd.overlay(grid_gdf, input_gdf, how='intersection')

                    if not intersection.empty:
                        # Calculate emission for the intersected grid cell:
                        intersection['emis'] = (intersection['geometry'].area / self.emiss_sc_factor) * intersection[self.spec_name]

                        # Calculate the left, top, right, and bottom values based on the intersection geometry:
                        intersection['left'] = intersection.geometry.bounds['minx']
                        intersection['top'] = intersection.geometry.bounds['maxy']
                        intersection['right'] = intersection.geometry.bounds['maxx']
                        intersection['bottom'] = intersection.geometry.bounds['miny']

                        try:
                            intersection[['left', 'bottom', self.spec_name, 'emis', 'geometry']].to_file(
                                rundir + 'inter.gpkg', driver='GPKG', geometry='MULTIPOLYGON')
                        except Exception as e:
                            print(f"Error saving emissions data to {rundir + 'inter.gpkg'}: {e}")

                        # Append the intersection results to the intersection_gdf:
                        intersection_gdf = pd.concat(
                                [intersection_gdf, intersection[['left', 'bottom', self.spec_name, 'emis', 'geometry']]],
                                ignore_index=True)

                        emis_buff_gdf = intersection_gdf[['left', 'bottom', self.spec_name, 'emis', 'geometry']].copy()
                        # Apply a small buffer to the geometry:
                        emis_buff_gdf['geometry'] = emis_buff_gdf['geometry'].buffer(-0.0000001, cap_style=1,
                                                                                     join_style=2, mitre_limit=2)

                        # Save the buffered data to a GeoPackage file:
                        #emis_buff_gdf.to_file(rundir + 'inter_buff.gpkg', driver='GPKG')

                        # ---- Broadcast total CO2 emission to target grid ----

                        # Load the buffered emissions GeoDataFrame:
                        #emis_buff_gdf = gpd.read_file(rundir + 'inter_buff.gpkg')
                        # Spatially join the two datasets and summarize emissions using 'sum' function:
                        joined_data = gpd.sjoin(grid_gdf, emis_buff_gdf, predicate='contains', how='left')
                        summarized_data = joined_data.groupby(joined_data.index)['emis'].sum().reset_index()
                        # Merge the summarized emissions back to the target grid :
                        result = grid_gdf.merge(summarized_data, left_index=True, right_index=True, how='left')
                        # Fill NaN values in the 'emis' column with 0:
                        result['emis'] = result['emis'].fillna(0)
                        harm_resi_gdf = result
                        # Rename the 'emis' column to self.spec_name:
                        harm_resi_gdf.rename(columns={'emis': self.spec_name}, inplace=True)
                        # Save the result to a GeoPackage file:
                        #harm_resi_gdf.to_file(rundir + 'harm_resi.gpkg', driver='GPKG')

                        # ---- Remove cells without emission ----
                        # Load the HARM residual emissions GeoDataFrame:
                        #harm_resi_gdf = gpd.read_file(rundir + 'harm_resi.gpkg')
                        # Filter the DataFrame using boolean indexing and make a copy:
                        harm_resi_filtered_gdf = harm_resi_gdf[harm_resi_gdf[self.spec_name] > 0].copy()
                        # Add columns for 'left', 'top', 'right', and 'bottom' based on the geometry:
                        harm_resi_filtered_gdf.loc[:, 'left'] = harm_resi_filtered_gdf.geometry.bounds['minx']
                        harm_resi_filtered_gdf.loc[:, 'top'] = harm_resi_filtered_gdf.geometry.bounds['maxy']
                        harm_resi_filtered_gdf.loc[:, 'right'] = harm_resi_filtered_gdf.geometry.bounds['maxx']
                        harm_resi_filtered_gdf.loc[:, 'bottom'] = harm_resi_filtered_gdf.geometry.bounds['miny']
                        # Add an 'id' column:
                        harm_resi_filtered_gdf['id'] = range(1, len(harm_resi_filtered_gdf) + 1)
                        # Reorder columns to match your desired format:
                        harm_resi_filtered_gdf = harm_resi_filtered_gdf[
                            ['id', 'left', 'top', 'right', 'bottom', self.spec_name, 'geometry']]

                        # Save the filtered data to the target file:
                        harm_resi_filtered_gdf.to_file(target_str, driver='GPKG')

                        # Optionally, you can clear the GeoDataFrame if needed:
                        grid_gdf = None
                        intersection_gdf = None
                        emis_buff_gdf = None
                        harm_resi_gdf = None
                        harm_resi_filtered_gdf = None

                        # Remove intermediate files:
                        for f in glob.glob(rundir + 'inter.gpkg'):
                            os.remove(f)
                        #for f in glob.glob(rundir + 'inter_buff.gpkg'):
                            #os.remove(f)
                        #for f in glob.glob(rundir + 'harm_resi.gpkg'):
                            #os.remove(f)

            print(f'SNAP{snap} is done.')

if __name__ == "__main__":
    # Instantiate the AreaEmissionsReprojection class
    AER=AreaEmissionsReprojection()
    # Run the reprojection and downscaling process
    AER.reproject_emissions()