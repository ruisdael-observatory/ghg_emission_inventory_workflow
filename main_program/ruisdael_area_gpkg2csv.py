"""
Module Description:
This script serves the purpose of translating GeoPackage (GPKG) files,
which have been early reprojected onto the HARMONIE grid emissions, into CSV-files.
This translation enables the subsequent retranslation of the reprojected emissions back to NetCDF format.
Notably, this script operates independently of the QGIS application,
as all QGIS-related functionalities have been meticulously substituted with analogous functions native to Python.

Creator Information:
- Creator: Dr. Arseni Doyennel
- Contact Email: a.doyennel@vu.nl
"""


import glob
import geopandas as gpd
import os
import numpy as np
import time
import os.path
import shutil
import math
import sys
import csv
import fiona
import pandas as pd

#from qgis.core import QgsVectorLayer, QgsProject, NULL
#sys.path.insert(1,os.path.join(os.getcwd(), ''))
from emission_preparation_setting import (spec_name, year, x0, y0, xn, yn,
                                          xres, yres, snaplist, dx, dy,
                                          distin_folder_ruisdael_area_csv2gpkg_new_python_only,
                                          outputdir_main_ruisdael_area_gpkg2csv_new_python_only ) #setting script


class EmissionPreparation:
    def __init__(self):
        self.rundir_main = distin_folder_ruisdael_area_csv2gpkg_new_python_only
        self.outputdir_main = outputdir_main_ruisdael_area_gpkg2csv_new_python_only
        self.spec_name = spec_name
        self.year = year
        self.x0 = x0
        self.y0 = y0
        self.xn = xn
        self.yn = yn
        self.xres = xres
        self.yres = yres
        self.snaplist = snaplist
        self.dx = int(dx)
        self.dy = int(dy)


    def prepare_emission_data(self):
        xmax_mod = self.x0 + int(self.xn * self.xres)
        ymax_mod = self.y0 + int(self.yn * self.yres)

        if not self.xn % self.dx and not self.yn % self.dy:
            print('dx/dy is a multiple of dx_scale/dy_scale!')
        else:
            raise Exception("check dx and dy before proceed!")

        xminlist = range(self.x0, self.x0 + int(self.xn * self.xres), int(self.dx * self.xres))
        yminlist = range(self.y0, self.y0 + int(self.yn * self.yres), int(self.dy * self.yres))

        # --- Merge parts
        for snap in self.snaplist:
            rundir = self.rundir_main + 'snap_' + str(snap) + '_' + self.spec_name + '_gpkg_HARM/'
            outputdir = self.outputdir_main + 'snap_' + str(snap) + '_csv/'

            # Check whether the specified path exists or not
            if not os.path.exists(outputdir):

                # Create a new directory because it does not exist
                os.makedirs(outputdir)
                print(f'The {outputdir} directory is created!')
            else:
                print(f'The {outputdir} directory exists!')

                # delete all files in the output directory--------------------------
                for f in glob.glob(outputdir + '*'):
                    os.remove(f)

            for xmin in xminlist:
                xmax = int(xmin + self.dx * self.xres)
                extname = str(xmin) + '_' + str(xmax)

                regex_source = 'HARM_snap_' + str(snap) + '_' + extname + '_*.gpkg'
                regex_target = 'HARM_snap_' + str(snap) + '_' + extname + '_aggr'
                regex_target_gpkg = regex_target + '.gpkg'
                regex_target_csv = regex_target + '.csv'

                if os.path.isfile(outputdir + regex_target_csv):
                    print(outputdir + regex_target_csv + ' exists')
                else:
                    print(outputdir + regex_target_csv + ' processing...')
                    filelist = glob.glob(os.path.join(rundir, regex_source))
                    merged_gdf = None       # Initialize the merged GeoDataFrame outside the loop

                    if filelist:
                        # Merge vector layers
                        for input_gpkg in filelist:
                            input_gdf = gpd.read_file(input_gpkg)
                            merged_gdf = gpd.GeoDataFrame(pd.concat([merged_gdf, input_gdf], ignore_index=True), crs=input_gdf.crs)
                        # Save the merged GeoDataFrame as GeoPackage
                        merged_gdf.to_file(os.path.join(outputdir, regex_target_gpkg), driver='GPKG')
                        # Read merged vector layer
                        with fiona.open(outputdir + regex_target_gpkg, 'r') as layer:
                            fieldnames = list(layer.schema['properties'].keys())
                        # Remove 'id' from the fieldnames list
                        fieldnames.remove('id')
                        # Specify the CSV path
                        csv_path = os.path.join(outputdir, regex_target_csv)
                        # Define the desired column order excluding 'id'
                        desired_order = ['left', 'top', 'right', 'bottom', self.spec_name]
                        # Rearrange fieldnames to match the desired order
                        fieldnames = [name for name in desired_order if name in fieldnames] + [name for name in fieldnames if name not in desired_order]
                        # Write to CSV
                        with open(csv_path, 'w', newline='') as output_file:
                            writer = csv.writer(output_file)
                            writer.writerow(fieldnames)
                            # Iterate over the features and write them to the CSV
                            for _, row in merged_gdf.iterrows():
                                writer.writerow([row[name] for name in fieldnames])
                        # Optionally, you can clear the GeoDataFrame if needed
                        merged_gdf = None
                        fieldnames = None
                        print('Output file written at ' + csv_path)
        print('Processing completed.')

if __name__ == "__main__":
    emission_prep = EmissionPreparation()
    emission_prep.prepare_emission_data()