"""
Module Description:
This script converts area emissions data (all SNAP categories except 2)
from CSV-files into GeoFrame GeoPackage (gpkg) files for subsequent downsizing
and reassignment to the HARMONIE grid.

The two reasons for this are that this format allows for the exact definition of spatial extent
and to select a subdomain of the Netherlands that is relevant for a simulation.
Note: Be sure that this subdomain contains the simulated domain (in HARMONIE CRS).

It operates independently of the QGIS application, utilizing standard Python functions
as replacements for QGIS-related operations.

Creator Information:
Author: Dr. Arseni Doyennel
Email: a.doyennel@vu.nl
"""

import os
import math
import glob
import shutil
import pyproj
from pyproj import Transformer
import pandas as pd
import geopandas as gpd
from shapely.geometry import Polygon


# Import required settings and libraries
from emission_preparation_setting import (spec_name, year, x0, y0, xn, yn,
                                          xres, yres, snaplist,rundir_ruisdael_area_residuals,
                                          distin_folder_ruisdael_area_csv2gpkg_new_python_only)


class EmissionConverter:
    def __init__(self):
        self.sourcedir = rundir_ruisdael_area_residuals
        self.distin_folder = distin_folder_ruisdael_area_csv2gpkg_new_python_only
        self.spec_name = spec_name
        self.year = year
        self.snaplist = [x for x in snaplist if x != 2]  # Exclude 2 (if no refinement, include it)
        self.x0 = x0
        self.y0 = y0
        self.xn = xn
        self.yn = yn
        self.dx = 5000
        self.dy = 5000
        self.xres = 1000
        self.yres = 1000
        self.xres_mod = xres
        self.yres_mod = yres

        self.create_directories()

    def create_directories(self):
        if not os.path.exists(self.distin_folder):
            os.makedirs(self.distin_folder)
            print(f'The {self.distin_folder} directory is created!')
        else:
            print(f'The {self.distin_folder} directory exists!')

    @staticmethod
    def roundup(x, base=5000):
        return int(math.ceil(x / base)) * base

    @staticmethod
    def create_grid(xmin, xmax, ymin, ymax, xres, yres):
        grid = []
        for x in range(round(xmin), round(xmax), xres):
            for y in range(round(ymin), round(ymax), yres):
                left, right, bottom, top = x, x + xres, y, y + yres
                grid.append(Polygon([(left, bottom), (right, bottom), (right, top), (left, top)]))
        return gpd.GeoDataFrame({'geometry': grid}, crs="EPSG:28992")

    def convert_emissions(self):

        # Domain settings
        proj_RD = pyproj.Proj('+proj=sterea +lat_0=52.15616055555555 +lon_0=5.38763888888889 +k=0.9999079 '
                              '+x_0=155000 +y_0=463000 +ellps=bessel +towgs84=565.417,50.3319,465.552,-0.398957,'
                              '0.343988,-1.8774,4.0725 +units=m +no_defs')
        proj_HARM = pyproj.Proj('+proj=lcc +lat_1=52.500000 +lat_2=52.500000 +lat_0=52.500000 +lon_0=.000000 '
                                '+k_0=1.0 +x_0=649536.512574 +y_0=1032883.739533 +a=6371220.000000 +b=6371220.000000')

        # Create a transformer for coordinate transformation
        transformer = pyproj.Transformer.from_proj(proj_HARM, proj_RD)

        # Calculate modified coordinates
        x0_rd, y0_rd = transformer.transform(self.x0, self.y0)
        xmax_rd, ymax_rd = transformer.transform(self.x0 + self.xn * self.xres_mod, self.y0 + self.yn * self.xres_mod)

        # Round up the coordinates
        x0_rd = self.roundup(x0_rd) - 5000
        y0_rd = self.roundup(y0_rd) - 5000
        xmax_rd = self.roundup(xmax_rd) + 5000
        ymax_rd = self.roundup(ymax_rd) + 5000

        # Calculate new grid dimensions and resolutions
        xn = int((xmax_rd - x0_rd) / self.dx)
        yn = int((ymax_rd - y0_rd) / self.dy)


        # Generate lists of xmin and ymin
        xminlist = list(range(x0_rd, x0_rd + xn * self.dx, self.dx))
        yminlist = list(range(y0_rd, y0_rd + yn * self.dy, self.dy))


        for snap in self.snaplist:
            rundir = os.path.join(self.distin_folder, f'snap_{snap}_{self.spec_name}_gpkg_RD/')
            if not os.path.exists(rundir):
                os.makedirs(rundir)
                print(f'The new directory SNAP {snap} is created in gpkg folder!')
            else:
                print(f'The directory SNAP {snap} exists in gpkg folder!')
                for f in glob.glob(rundir + '*.gpkg*') + glob.glob(rundir + '*.csv*'):
                    os.remove(f)
            source = os.path.join(self.sourcedir, f'{self.spec_name}_{self.year}_SNAP_{snap}_residual.csv')
            destination = os.path.join(rundir, f'{self.spec_name}_{self.year}_SNAP_{snap}_residual.csv')
            if os.path.isfile(source):
                shutil.copy(source, destination)
                print(f'copied {self.spec_name}_{self.year}_SNAP_{snap}_residual.csv')
            for xmin in xminlist:
                for ymin in yminlist:
                    xmax = xmin + self.dx
                    ymax = ymin + self.dy
                    extname = f"{xmin}_{xmax}_{ymin}_{ymax}"
                    grid = self.create_grid(xmin, xmax, ymin, ymax, self.xres, self.yres)
                    grid.to_file(os.path.join(rundir, "grid.gpkg"), driver="GPKG")
                    target_str = os.path.join(rundir, f"ER_snap_{snap}_{extname}.gpkg")
                    if os.path.isfile(target_str):
                        print(target_str + " exists")
                    else:
                        print(target_str + " processing...")
                    csv_file = os.path.join(rundir, f"{self.spec_name}_{self.year}_SNAP_{snap}_residual.csv")
                    snap_df = pd.read_csv(csv_file, delimiter=",", usecols=['x', 'y', 'co2'])
                    snap_gdf = gpd.GeoDataFrame(
                        {
                            'id': range(1, len(snap_df) + 1),
                            'left': [xmin] * len(snap_df),
                            'top': [ymax] * len(snap_df),
                            'right': [xmax] * len(snap_df),
                            'bottom': [ymin] * len(snap_df),
                            self.spec_name: snap_df[self.spec_name],
                            'layer': f"ER_snap_{snap}_{extname}",
                            'path': target_str,
                        },
                        geometry=gpd.points_from_xy(snap_df.x, snap_df.y),
                        crs="EPSG:28992"
                    )
                    grid = grid.to_crs("EPSG:28992")
                    join_gdf = gpd.sjoin(grid, snap_gdf, how="left", predicate="intersects")
                    snap_gdf = None
                    grid = None
                    join_gdf = join_gdf.drop(columns=['index_right'])
                    join_gdf = join_gdf[join_gdf[self.spec_name] > 0]
                    if not join_gdf.empty:
                        join_gdf['id'] = range(1, len(join_gdf) + 1)
                        join_gdf.to_file(target_str, driver="GPKG")
                        print(f"Processing completed for {target_str}")
                    else:
                        print(f"Warning: {target_str} is empty and will not be saved.")
                    join_gdf = None
                    print(f"Processing completed for {target_str}")

        print(f'SNAP{snap} is done.')
    def merge_and_remove_data(self):
        for snap in self.snaplist:
            merged_gdf = gpd.GeoDataFrame()
            rundir = os.path.join(self.distin_folder, f'snap_{snap}_{self.spec_name}_gpkg_RD')
            filelist = []
            for name in glob.glob(os.path.join(rundir, f'ER_snap_{snap}_*.gpkg')):
                filelist.append(name)
            for file in filelist:
                gdf = gpd.read_file(file)
                merged_gdf = pd.concat([merged_gdf, gdf], ignore_index=True)
            merged_gpkg = os.path.join(rundir, f'{self.spec_name}_{self.year}_SNAP_{snap}_residual.gpkg')
            merged_gdf.to_file(merged_gpkg, driver="GPKG")
            merged_gdf = None

        # Remove intermediate files
        for snap in self.snaplist:
            rundir = os.path.join(self.distin_folder, f'snap_{snap}_{self.spec_name}_gpkg_RD')
            for f in glob.glob(rundir + '/ER_snap_'+str(snap)+'_*0.gpkg*'):
                os.remove(f)




if __name__ == "__main__":
    # Instantiate the EmissionConverter class
    EC = EmissionConverter()
    # Run the conversion process
    EC.convert_emissions()
    # Merge the generated GeoPackages and remove intermediate files
    EC.merge_and_remove_data()
#