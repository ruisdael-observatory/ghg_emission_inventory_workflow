"""
Module Description:
This script converts area emissions data from CSV-files into GeoFrame GeoPackage (gpkg) files,
exclusively for SNAP 2 (consumption category), used for subsequent downscaling and coordinate
reassignment onto the HARMONIE grid.

The two reasons for this are that this format allows for the exact definition of spatial extent
and to select a subdomain of the Netherlands that is relevant for a simulation.
Note: Be sure that this subdomain contains the simulated domain (in HARMONIE CRS).

It operates independently of the QGIS application, leveraging standard Python functions
to replace QGIS-related operations.

Author: Dr. Arseni Doyennel
Email: a.doyennel@vu.nl
"""

# Import required settings and libraries
import os
import math
import glob
import shutil
import pyproj
from pyproj import Transformer
import pandas as pd
import geopandas as gpd
from shapely.geometry import Polygon

# Set your input and output directories (from setting)

# Import required settings and libraries
from emission_preparation_setting import (spec_name, year, x0, y0, xn, yn,
                                          distin_folder_ruisdael_area_csv2gpkg_new_python_only,
                                          rundir_ruisdael_area_residuals, xres, yres)


class EmissionConverter_consu:
    def __init__(self):
        # Set your input and output directories
        self.sourcedir = rundir_ruisdael_area_residuals
        self.rundir_main = distin_folder_ruisdael_area_csv2gpkg_new_python_only
        self.snaplist = ['2']  # Always 2 in this script
        self.spec_name = spec_name
        self.co2_column = spec_name + '_sum'
        self.xres = 100        #the input csv of snap2 is 100x100m!
        self.yres = 100        #the input csv of snap2 is 100x100m!
        self.dx = 5000  # to fit this check your xmin and xmax (in RD), in principle should be ok for any domains, except too little
        self.dy = 5000  # to fit this check your xmin and xmax (in RD), in principle should be ok for any domains, except too little
        self.xn = xn
        self.yn = yn
        self.x0 = x0
        self.y0 = y0
        self.xres_mod = xres    #the resolution of the current domain in HARM (from settings)
        self.yres_mod = yres    #the resolution of the current domain in HARM  (from settings)



    @staticmethod
    def roundup(x):
        return int(math.ceil(x / 1000.0)) * 1000

    @staticmethod
    def create_grid(xmin, xmax, ymin, ymax, xres, yres):
        grid = []
        for x in range(round(xmin), round(xmax), xres):
            for y in range(round(ymin), round(ymax), yres):
                left, right, bottom, top = x, x + xres, y, y + yres
                grid.append(Polygon([(left, bottom), (right, bottom), (right, top), (left, top)]))
        return gpd.GeoDataFrame({'geometry': grid}, crs="EPSG:28992")

    def prepare_emission_data(self):

        # Copy input file to the rundir ------------------------------------

        # Delete all files in the rundir directory
        rundir = self.rundir_main + 'snap_2_' + self.spec_name + '_gpkg_RD/'
        if not os.path.exists(rundir):
            os.makedirs(rundir)
            print(f'The {rundir} directory is created!')
        else:
            print(f'The {rundir} directory exists!')
            for f in glob.glob(rundir + '*'):
                os.remove(f)

        # Domain settings
        proj_RD = pyproj.Proj('+proj=sterea +lat_0=52.15616055555555 +lon_0=5.38763888888889 +k=0.9999079 '
                              '+x_0=155000 +y_0=463000 +ellps=bessel +towgs84=565.417,50.3319,465.552,-0.398957,'
                              '0.343988,-1.8774,4.0725 +units=m +no_defs')
        proj_HARM = pyproj.Proj('+proj=lcc +lat_1=52.500000 +lat_2=52.500000 +lat_0=52.500000 +lon_0=.000000 '
                                '+k_0=1.0 +x_0=649536.512574 +y_0=1032883.739533 +a=6371220.000000 +b=6371220.000000')

        transformer = Transformer.from_proj(proj_HARM, proj_RD)

        x0_rd, y0_rd = transformer.transform(self.x0,  self.y0)

        x0_rd = self.roundup(x0_rd) - 1000
        y0_rd = self.roundup(y0_rd) - 1000

        xn_mod = self.xn  # from namoptions of your experiment number of gridcells in x direction (itot)
        yn_mod = self.yn  # from namoptions of your experiment number of gridcells in y direction (jtot)

        xres_mod = self.xres_mod
        yres_mod = self.yres_mod

        xmax_mod = self.x0 + int(xn_mod * xres_mod)
        ymax_mod = self.y0 + int(yn_mod * yres_mod)

        transformer = Transformer.from_proj(proj_HARM, proj_RD)
        xmax_mod_rd, ymax_mod_rd = transformer.transform(xmax_mod, ymax_mod)

        xmax_mod_rd = self.roundup(xmax_mod_rd) + 1000
        ymax_mod_rd = self.roundup(ymax_mod_rd) + 1000

        xn = (xmax_mod_rd - x0_rd) / self.dx
        xn = int(xn) + 1

        yn = (ymax_mod_rd - y0_rd) / self.dx
        yn = int(yn) + 1

        source = self.sourcedir + '' + str(self.spec_name) + '_' + str(
            year) + '_SNAP_2_refined_' + str(self.xres) + 'm_residual.csv'
        destination = rundir + '' + str(self.spec_name) + '_' + str(year) + '_SNAP_2_refined_' + str(
            self.xres) + 'm_residual.csv'

        print(source)

        # copy only files
        if os.path.isfile(source):
            shutil.copy(source, destination)
            print('copied',
                  '' + str(self.spec_name) + '_' + str(year) + '_SNAP_2_refined_' + str(self.xres) + 'm_residual.csv')

        # Iterate over snaplist
        for snap in self.snaplist:
            for xmin in range(x0_rd, x0_rd + xn * self.dx, self.dx):
                for ymin in range(y0_rd, y0_rd + yn * self.dy, self.dy):
                    xmax = xmin + self.dx
                    ymax = ymin + self.dy
                    extname = f"{xmin}_{xmax}_{ymin}_{ymax}"

                    # Create target grid
                    grid = self.create_grid(xmin, xmax, ymin, ymax, self.xres, self.yres)
                    grid.to_file(os.path.join(rundir, "grid_consu.gpkg"), driver="GPKG")
                    target_str = os.path.join(rundir, f"grid_consu_{self.spec_name}_{extname}.gpkg")

                    if os.path.isfile(target_str):
                        print(target_str + " exists")
                    else:
                        print(target_str + " processing...")

                    # Read the CSV file into a DataFrame and convert it to a GeoDataFrame
                    csv_file = os.path.join(rundir,
                                            f"{rundir}{self.spec_name}_{year}_SNAP_2_refined_{self.xres}m_residual.csv")
                    snap_df = pd.read_csv(csv_file, delimiter=",", usecols=['x', 'y', self.co2_column])

                    snap_gdf = gpd.GeoDataFrame(
                        {
                            'id': range(1, len(snap_df) + 1),  # Assign unique IDs
                            'left': [xmin] * len(snap_df),
                            'top': [ymax] * len(snap_df),
                            'right': [xmax] * len(snap_df),
                            'bottom': [ymin] * len(snap_df),
                            self.co2_column: snap_df[self.co2_column],
                            'layer': f"grid_consu_{self.spec_name}_{extname}",
                            'path': target_str,
                        },
                        geometry=gpd.points_from_xy(snap_df.x, snap_df.y),
                        crs="EPSG:28992"
                    )

                    # Reproject the grid to EPSG:28992
                    grid = grid.to_crs("EPSG:28992")

                    # Perform a spatial join between the grid and snap data (intersects method)
                    join_gdf = gpd.sjoin(grid, snap_gdf, how="left", predicate="intersects")

                    # Drop unnecessary columns if needed
                    join_gdf = join_gdf.drop(columns=['index_right'])

                    # Remove cells without emission
                    join_gdf = join_gdf[join_gdf[self.co2_column] > 0]

                    # Check if join_gdf is empty before saving
                    if not join_gdf.empty:
                        # Save the GeoDataFrame as the final GeoPackage
                        join_gdf.to_file(target_str, driver="GPKG")
                        print(f"Processing completed for {target_str}")
                    else:
                        print(f"Warning: {target_str} is empty and will not be saved.")

    def merge_and_remove_data(self):
        rundir = self.rundir_main + 'snap_2_' + self.spec_name + '_gpkg_RD/'

        # Merge all grid_consu_{self.spec_name}_{extname}.gpkg into a single Geoframe file
        merged_gdf = gpd.GeoDataFrame()

        for snap in self.snaplist:
            filelist = []

            # Use glob to find grid_consu files with the specified snap value
            for name in glob.glob(os.path.join(rundir, f'grid_consu_{self.spec_name}_*.gpkg')):
                filelist.append(name)

            # Merge the grid_consu files into a single GeoDataFrame
            for file in filelist:
                gdf = gpd.read_file(file)
                merged_gdf = pd.concat([merged_gdf, gdf], ignore_index=True)

        # Define the output path for the merged GeoPackage
        merged_gpkg = os.path.join(rundir, f'' + str(self.spec_name) + '_' + str(year) + '_SNAP_2_refined_' +
                                   str(self.xres) + 'm_residual.gpkg')

        # Save the merged data to a new GeoPackage
        merged_gdf.to_file(merged_gpkg, driver="GPKG")

        # Remove intermediate files
        for f in glob.glob(rundir + '/grid_consu_' + str(self.spec_name) + '_*0.gpkg*'):
            os.remove(f)


if __name__ == "__main__":
    EC = EmissionConverter_consu()
    EC.prepare_emission_data()
    EC.merge_and_remove_data()



# Define the path to your GeoPackage file

#rundir21 = rundir
#gpkg_file1 = f'co2_{year}_SNAP_2_refined_100m_residual.gpkg'

# Read the GeoPackage file into a GeoDataFrame
#gdf1 = gpd.read_file(rundir21+gpkg_file1)

#import matplotlib.pyplot as plt
#import geopandas as gpd
#from matplotlib import colors

# Assuming gdf1 is your GeoDataFrame with a 'co2_sum' column

# Define the logarithmic color scale
#norm = colors.LogNorm(vmin=1e3, vmax=1e6)

# Create a plot of the GeoDataFrame
#ax = gdf1.plot(column='co2_sum', cmap='Spectral_r', norm=norm, legend=True)

# Add color bar
#cbar = ax.get_figure().get_axes()[1]
#cbar.set_ylabel('CO2 Sum')

# Set plot title and labels
#plt.title('CO2 emissions SNAP2 (consu)')
#plt.xlabel('X Coordinate')
#plt.ylabel('Y Coordinate')

