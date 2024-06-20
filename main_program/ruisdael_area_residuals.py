"""
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

(in general, point_source_rd_file might include all open_source_data_points, but check this)

mapped_data: This file contains columns such as DATASET,
GEBIEDSINDELING, CODE_GEBIED, GEBIED, DOELGROEP, SUBDOELGROEP,
COMPARTIMENT, STOFCODE, STOF, EENHEID, EMISSIE_KG, delimited by a semicolon (;).

Output:
This script generates CSV files with residual area emissions,
segregated per SNAP category, e.g., co2_2017_SNAP_1_residual.csv.
Notably, for SNAP2, the resolution is 100m, while for other categories is 1000m.

Creator Information:
Creators: Dr. Marko de Bruine, Dr. Arseni Doyennel
Email: a.doyennel@vu.nl
"""
import numpy as np
import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import pandas as pd
import matplotlib.patches as patches
import netCDF4 as netc
import glob
import os
import shutil
from numpy.lib.stride_tricks import as_strided
from emission_construction_functions import *
import sys
import pyproj as pyproj
from pyproj import Proj, transform
from pyproj import Transformer


sys.path.insert(1, os.path.join(os.path.dirname(os.getcwd()), ''))
from emission_preparation_setting import (
    year, spec_name, year_start, snaplist, open_source_data_points,
    mapped_data,
    point_source_rd_file, datadir_ruisdael_area_residuals,
    rundir_ruisdael_area_residuals, cbs_loc_consu_preprocessing
)


class Emission_initial_process_get_residuals:
    def __init__(self):
        self.datadir = datadir_ruisdael_area_residuals
        self.rundir = rundir_ruisdael_area_residuals
        self.cbs_dir = cbs_loc_consu_preprocessing
        self.cbs_name = f'cbs_vk100_{year_start}.gpkg'  #name of your CBS file (household data)
        self.year = year
        self.spec_name = spec_name
        self.snaplist = [x for x in snaplist if
                    x not in [7, 2]]  # Exclude consu (2) and traffic (7), since they are processed differently
        self.DoPlot = False
        self.DoPrivate = False
        self.DoNetCDF = False
        self.DoCSV = True
        self.DO_SNAP_2 = True
        self.DO_SNAP_7 = True
        self.crs = 'RD'
        self.nsnap = 10     #always 10 here as we have 10 SNAP categories anyway in files
        self.x_start, self.x_end, self.dx = 0, 300000, 1000
        self.y_start, self.y_end, self.dy = 300000, 900000, 1000
        self.subdoelgroepen_exclude = ['Smeermiddelengebruik-verkeer',
                                       'Wegverkeer - uitlaatgassen',
                                       'Wegverkeer - niet uitlaatgassen',
                                       'Railverkeer']   #Exclude categories included in SNAP7 as roads have no point sources to extract
        self.subdoelgroepen_exclude_point = []
        self.subdoelgroepen_include = self.subdoelgroepen_exclude

        #################################################################

        
        
        
        
    def read_point(self):
        df_point = pd.read_csv(self.datadir + open_source_data_points, delimiter=';')
        df_point_private = pd.read_csv(self.datadir + point_source_rd_file, delimiter=';', decimal=",",
                                       encoding="ISO-8859-1")
        df_point_private['SUBDOELGROEP'] = 'LEEG'
        year_to_keep = str(year)
        numeric_columns = [col for col in df_point.columns if col.isnumeric()]
        columns_to_keep = [col for col in df_point.columns if not col.isnumeric()]
        df_point = df_point[columns_to_keep + [year_to_keep]]

        return df_point, df_point_private

    def preprocess_point_data(self, df_point, df_point_private):
        # Add the first part of the preprocessing
        for index, row in df_point_private.iterrows():
            df_point, df_point_private = self.process_row(index, row, df_point, df_point_private)

        df_point = df_point_private

        # Replace points without X/YCO_EMISSIEPUNT with X/YCO_BEDRIJF
        df_noloc = df_point[(~np.isfinite(df_point['XCO_EMISSIEPUNT'])) &
                            (df_point['SUBDOELGROEP'] != 'Luchtvaart')]

        for index, row in df_noloc.iterrows():
            df_point.loc[index, 'XCO_EMISSIEPUNT'] = row['XCO_BEDRIJF']
            df_point.loc[index, 'YCO_EMISSIEPUNT'] = row['YCO_BEDRIJF']

        return df_point

    def process_row(self, index, row, df_point, df_point_private):
        df_hit = df_point[df_point['NIC'] == row['CODE_BEDRIJF']]
        ngroep = len(df_hit['SUBDOELGROEP'].unique())

        if ngroep > 1:
            df_hit = df_hit[df_hit['CODE_EMISSIEOORZAAK'] == row['CODE_EMISSIEOORZAAK']]
            ngroep = len(df_hit['SUBDOELGROEP'].unique())

        if ngroep == 1:
            df_point_private.at[index, 'SUBDOELGROEP'] = df_hit['SUBDOELGROEP'].unique()[0]
        else:
            print(f'Nothing found for: {row["EMISSIEOORZAAK"]}')

        return df_point, df_point_private

    
    def save_points_harm(self, df_point):
        
        proj_RD   = pyproj.Proj('+proj=sterea +lat_0=52.15616055555555 +lon_0=5.38763888888889 +k=0.9999079 +x_0=155000 +y_0=463000 +ellps=bessel +towgs84=565.417,50.3319,465.552,-0.398957,0.343988,-1.8774,4.0725 +units=m +no_defs')
        proj_HARM = pyproj.Proj('+proj=lcc +lat_1=52.500000 +lat_2=52.500000 +lat_0=52.500000 +lon_0=.000000 +k_0=1.0 +x_0=649536.512574 +y_0=1032883.739533 +a=6371220.000000 +b=6371220.000000')
        
        df_point_P = df_point[(df_point['TYPE'] == 'P')].copy()  # Create a copy of the DataFrame with P-type emis only
        
        idx_x_emiss=df_point_P.columns.get_loc("XCO_EMISSIEPUNT")
        idx_y_emiss=df_point_P.columns.get_loc("YCO_EMISSIEPUNT")

        
        ff=df_point_P.to_numpy()
        x_emiss=ff[:,idx_x_emiss]
        y_emiss=ff[:,idx_y_emiss]
        
        xemissHARM=np.zeros(np.size(x_emiss))
        yemissHARM=np.zeros(np.size(y_emiss))
        
        for j in range(0,np.size(x_emiss)):    
    
            transformer = Transformer.from_proj(proj_RD, proj_HARM)
            xemissHARM[j], yemissHARM[j] = transformer.transform(x_emiss[j], y_emiss[j])
        
        
        xemissHARM = xemissHARM.tolist()
        yemissHARM = yemissHARM.tolist()
        
        # Using '' as the column name
        # and equating it to the list
        df_point_P['XCO_EMISSIEPUNT_HARM'] = xemissHARM
        df_point_P['YCO_EMISSIEPUNT_HARM'] = yemissHARM
        
        # Assuming self.targetdir contains the directory where you want to save the CSV file
        csv_file_path = os.path.join(self.datadir, f'ERI_{self.spec_name}_{self.year}_metemk_harm_P_only.csv')
        
        for f in glob.glob(csv_file_path):
            os.remove(f)

        # Save df_point_P to CSV
        df_point_P.to_csv(csv_file_path, index=False)
    
    
    def substract_point_from_map(self, df_point):

        print('substract of point source from area emission processing...')

        x = np.arange(self.x_start, self.x_end, self.dx).astype(float)
        y = np.arange(self.y_start, self.y_end, self.dy).astype(float)
        df_map = pd.read_csv(self.datadir + mapped_data, delimiter=';')
        field_source = np.zeros([self.nsnap, len(x), len(y)])

        for i, item in enumerate(df_map['SUBDOELGROEP'].unique()):
            isnap = determine_snap(item)
            print(isnap, item, )

            if item not in self.subdoelgroepen_exclude:
                df_map_subdoelgroep = df_map[df_map['SUBDOELGROEP'] == item]
                # df_map_subdoelgroep.loc[:, 'EMISSIE_KG'] = df_map_subdoelgroep['EMISSIE_KG'].str.replace(',','.')

                # PUT MAP DATA ON 2D FIELD
                field_subdoelgroep_source = np.zeros([len(x), len(y)])
                load_ERmap(df_map_subdoelgroep, field_subdoelgroep_source, x, y, year=year)

                if item not in self.subdoelgroepen_exclude_point:
                    df_point_subdoelgroep = df_point[df_point['SUBDOELGROEP'] == item]

                    if len(df_point_subdoelgroep) > 0:
                        print("points", len(df_point_subdoelgroep))
                        join_map_points(df_point_subdoelgroep, field_subdoelgroep_source,
                                        x, y,
                                        mode='subtract',
                                        dim=2,
                                        xloc='XCO_EMISSIEPUNT',
                                        yloc='YCO_EMISSIEPUNT')

                field_source[isnap - 1] += field_subdoelgroep_source
            print('')

        del field_subdoelgroep_source
        print('Done.')
        return field_source, df_map, x, y

    def populate_road_emissions(self, df_map, x, y):
        # Initialize the emission field
        field_subdoelgroep_source_roads = np.zeros([len(x), len(y)])

        # Iterate over unique subdoelgroepen
        for i, item in enumerate(df_map['SUBDOELGROEP'].unique()):
            if item in self.subdoelgroepen_include:
                print(item)
                # Filter data for the specific subdoelgroep
                df_map_subdoelgroep = df_map[df_map['SUBDOELGROEP'] == item]

                # Split data into 5km and 1km dataframes
                df_5km = df_map_subdoelgroep[df_map_subdoelgroep['CODE_GEBIED'].str.contains('km5')]
                df_1km = df_map_subdoelgroep[df_map_subdoelgroep['CODE_GEBIED'].str.contains('km1')]

                # Iterate over 5km data
                if len(df_5km):
                    for index, row in df_5km.iterrows():
                        mapemission = float(row['EMISSIE_KG'].replace(',', '.'))
                        xidx = np.argmax(x >= 1e3 * float(row['CODE_GEBIED'][1:4]))
                        yidx = np.argmax(y >= 1e3 * float(row['CODE_GEBIED'][5:8]))
                        field_subdoelgroep_source_roads[xidx:xidx + 5, yidx:yidx + 5] += mapemission / 25

                # Iterate over 1km data
                if len(df_1km):
                    for index, row in df_1km.iterrows():
                        mapemission = float(row['EMISSIE_KG'].replace(',', '.'))
                        xidx = np.argmax(x >= 1e3 * float(row['CODE_GEBIED'][0:3]))
                        yidx = np.argmax(y >= 1e3 * float(row['CODE_GEBIED'][3:6]))
                        field_subdoelgroep_source_roads[xidx, yidx] += mapemission

        return field_subdoelgroep_source_roads

    def save_data(self, field_source, field_subdoelgroep_source_roads, x, y):

        if not os.path.exists(self.rundir):
            os.makedirs(self.rundir)
            print(f'The {self.rundir} directory is created!')
        else:
            print(f'The {self.rundir} directory exists!')
            for f in glob.glob(os.path.join(self.rundir, f'{self.spec_name}*.csv')):
                os.remove(f)
            for f in glob.glob(os.path.join(self.rundir, f'{self.spec_name}*.nc')):
                os.remove(f)

        for isnap in self.snaplist:
            fieldsnap = field_source[isnap - 1, :, :].T
            if self.DoNetCDF:
                fobj = netc.Dataset(os.path.join(self.rundir, f'{self.spec_name}_{self.year}_SNAP_{isnap}_residual.nc'),
                                    'w')
                fobj.description = f"emissions of {self.spec_name} are not attributable to point sources (year {self.year}, SNAP category {isnap})"
                dim_x = fobj.createDimension('x', len(x))
                dim_y = fobj.createDimension('y', len(y))
                var_x = fobj.createVariable('x', 'f4', ('x',))
                var_y = fobj.createVariable('y', 'f4', ('y',))
                var_x[:] = np.array(x)
                var_y[:] = np.array(y)
                var_e = fobj.createVariable(self.spec_name, 'f8', ('y', 'x',))
                var_e[:, :] = fieldsnap
                if self.crs == 'RD':
                    var_x.units = 'Rijksdriehoekcoordinaat x in meters'
                    var_y.units = 'Rijksdriehoekcoordinaat y in meters'
                elif self.crs == 'HARM':
                    var_x.units = 'Harmonie coordinate x in meters'
                    var_y.units = 'Harmonie coordinate y in meters'
                var_e.units = 'Kilogram per year'
                fobj.close()

            if self.DoCSV:
                x_csv = x + (x[1] - x[0]) / 2
                y_csv = y + (y[1] - y[0]) / 2
                xv, yv = np.meshgrid(x_csv.astype(int), y_csv.astype(int))
                df = pd.DataFrame(np.array([xv.flatten(), yv.flatten(), fieldsnap.flatten()]).T,
                                  columns=['x', 'y', self.spec_name])
                df = df.loc[df[self.spec_name] > 0]
                df.to_csv(os.path.join(self.rundir, f'{self.spec_name}_{self.year}_SNAP_{isnap}_residual.csv'),
                          index=False)
                print(f'Done for SNAP category {isnap}.\n')

        if self.DoCSV and self.DO_SNAP_7:
            fieldroad = field_subdoelgroep_source_roads[:, :].T
            x_csv = x + (x[1] - x[0]) / 2
            y_csv = y + (y[1] - y[0]) / 2
            xv, yv = np.meshgrid(x_csv.astype(int), y_csv.astype(int))
            df = pd.DataFrame(np.array([xv.flatten(), yv.flatten(), fieldroad.flatten()]).T,
                              columns=['x', 'y', self.spec_name])
            df = df.loc[df[self.spec_name] > 0]
            df.to_csv(os.path.join(self.rundir, f'{self.spec_name}_{self.year}_SNAP_7_residual.csv'), index=False)
            print(f'Done for SNAP category traffic (SNAP 7).\n')

    ##################################################################

    # This method is used to refine residual area emission for
    # the SNAP 2 category (consummenten) category using household statistical data (CBS files)
    # and saving the data into apropriate csv-file for further reprojection and reasignment to HARM grid.
    # Note: cbs_vk100_....gpkg can be downloaded from:
    # https://www.cbs.nl/nl-nl/dossier/nederland-regionaal/geografische-data/kaart-van-100-meter-bij-100-meter-met-statistieken
    # these files are available for different years. Here,cbs_vk100 file is read for the year specified in year_start:

    def refine_and_save_snap2_data(self, field_source):

        print('Refinement of area emission SNAP2 category (consu) processing...')
        print('')

        if not os.path.exists(self.rundir):
            os.makedirs(self.rundir)
            print(f'The {self.rundir} directory is created!')
        else:
            print(f'The {self.rundir} directory exists!')
            for f in glob.glob(os.path.join(self.rundir, '*SNAP_2_refined_100m_residual.csv')):
                os.remove(f)

        data = gpd.read_file(os.path.join(self.cbs_dir, self.cbs_name),
                             include_fields=["crs28992res100m","aantal_inwoners", "gemiddeld_gasverbruik_woning","gemiddeld_elektriciteitsverbruik_woning",  "geometry"])

        # Use 'gemiddeld_gasverbruik_woning' primarily, fall back to 'aantal_inwoners' if missing
        gasverbruik = data["gemiddeld_gasverbruik_woning"].fillna(0)
        aantal_inwoners = data["aantal_inwoners"].fillna(0)
        data2 = data['geometry']

        coordpolygonX = np.zeros(np.size(data2))
        coordpolygonY = np.zeros(np.size(data2))

        for i in range(0, np.size(data2)):
            data3 = data2[i]
            a = np.hstack([np.array(g.boundary.xy) for g in data3.geoms])
            coordpolygonX[i] = np.round(a[0, 3])
            coordpolygonY[i] = np.round(a[1, 3])

        x_start, x_end, dx = self.x_start, self.x_end, 100 #dx is 100 because CRS data are 100x100m
        y_start, y_end, dy = self.y_start, self.y_end, 100 #dy is 100 because CRS data are 100x100m

        x = np.arange(x_start, x_end, dx).astype(float)
        y = np.arange(y_start, y_end, dy).astype(float)

        data_selected = np.zeros([len(x), len(y), 3])

        for i in range(0, np.size(data2)):
            x_shape, y_shape = coordpolygonX[i], coordpolygonY[i]
            xidx = np.argmax(x >= x_shape) - 1
            yidx = np.argmax(y >= y_shape) - 1
            
            # Prioritize gasverbruik but fall back to aantal_inwoners if gasverbruik is zero:
            if gasverbruik[i] > 0:
                data_selected[xidx, yidx, 0] = gasverbruik[i]
            else:
                data_selected[xidx, yidx, 0] = aantal_inwoners[i]

        inwoners_100_100 = data_selected[..., 0].copy()
        inwoners_100_100[inwoners_100_100 < 0] = 0.

        inwoners_1k_1k = shrink(inwoners_100_100, len(x) // 10, len(y) // 10)
        inwoners_1k_100 = tile_array(inwoners_1k_1k, 10, 10)

        # save CSV-file with emissions:
        # spatial aggegation of consu emiss using a number of people settled:

        if self.DoCSV:
            # Check for zero values in the denominator (inwoners_1k_100) before division
            denominator = inwoners_1k_100
            denominator[denominator == 0] = np.nan # Replace zeros with NaN to avoid division by zero
            # Calculate the result while handling division by zero or invalid values
            with np.errstate(divide='ignore', invalid='ignore'):
                # do not sure what is better to use here, gas, elect or number of sitizen (number of sitizen now)
                spec_sum = np.nan_to_num(tile_array(field_source[1], 10, 10) * (inwoners_100_100 / denominator))

            spec_sum = spec_sum.T

            xx = x
            yy = y

            x_csv = xx + (xx[1] - xx[0]) / 2
            y_csv = yy + (yy[1] - yy[0]) / 2

            xv, yv = np.meshgrid(x_csv.astype(int), y_csv.astype(int))
            df = pd.DataFrame(np.array([xv.flatten(), yv.flatten(), spec_sum.flatten()]).T,
                              columns=['x', 'y', f'{self.spec_name}_sum'])
            df.replace(np.nan, 0)
            df = df.loc[df[f'{self.spec_name}_sum'] > 0]

            df.to_csv(os.path.join(self.rundir, f'{self.spec_name}_{self.year}_SNAP_2_refined_100m_residual.csv'), index=False)

            print(f'Done with Consumenten category (SNAP 2).\n')


if __name__ == "__main__":
    # Instantiate the class and call the methods:
    emission_init_prep = Emission_initial_process_get_residuals()
    df_point, df_point_private = emission_init_prep.read_point()
    df_point_2 = emission_init_prep.preprocess_point_data(df_point, df_point_private)
    emission_init_prep.save_points_harm(df_point_2)
    field_source, df_map, x, y = emission_init_prep.substract_point_from_map(df_point_2)
    field_subdoelgroep_source_roads = emission_init_prep.populate_road_emissions(df_map, x, y)
    emission_init_prep.save_data(field_source, field_subdoelgroep_source_roads, x, y)
    emission_init_prep.refine_and_save_snap2_data(field_source)
