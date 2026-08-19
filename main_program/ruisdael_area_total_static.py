"""
Module Description:

This script is designed to generate static annual total area emissions within a specified domain.
Additionally, it incorporates the following functionalities:

Refinement of area emissions using point source data (if point sources are not utilized as explicit input).
Generation of 3-D emission files per SNAP category to facilitate the final preparation of area emission inputs.
( area emission are devided by number of vertical layers, covered with plume until the emission top height layer!)
On vertical destribution of emissions of emissions: 

Vertical distribution of gridded emissions by SNAP category is based on typical stack heights and plume rise values reported in the literature and emission inventories.

- SNAP 1 (Power generation): Typical stack heights range from 70 to 150 m; a plume height of ~80 m with vertical spread (sigma) ~20 m reflects common values.
  Reference: EMEP/EEA Air Pollutant Emission Inventory Guidebook 2019, Chapter on Large Combustion Plants
  (https://www.eea.europa.eu/publications/emep-eea-guidebook-2019/part-b-sectoral-guidance-chapters/1-energy/1-a-combustion/1-a-1-large-combustion-plants)

- SNAP 3 (Industrial combustion): Stack heights usually range between 30 and 80 m; here ~60 m plume height with sigma ~15 m is used.
  Reference: Dutch National Inventory Report 2021, Annex on Emission Factors and Stack Parameters
  (https://unfccc.int/documents/267887)

- SNAP 4 (Industrial processes): Often lower stacks or area sources; plume height ~60 m with sigma ~15 m chosen as proxy.

- SNAP 8 (Mobile sources, including shipping): Ship stack heights commonly around 20-40 m, plume height ~30 m, sigma ~10 m.
  Reference: EMEP/EEA Guidebook 2019, Chapter on Mobile Sources

- SNAP 9 (Waste): Emission sources like waste incinerators or landfill gas flares typically have stack heights around 20-40 m; plume height ~20 m, sigma ~15 m is a reasonable assumption.

- SNAP categories without known vertical components (2, 5, 7, 10) are assigned emissions to the first model layer only, reflecting ground-level or diffuse sources.

This approach applies only to gridded/mapped emissions lacking explicit stack or plume data. Point sources with available stack parameters and plume rise calculations are treated separately.

This methodology is consistent with common practice in air quality modeling and emission inventory processing (e.g., EMEP/EEA Guidebook, Dutch emission inventories).

Author: Dr. Arseni Doyennel
Contact Email: a.doyennel@vu.nl
"""

import matplotlib.pyplot as plt
import numpy as np
import netCDF4 as netc
import pandas as pd
import matplotlib.colors as colors
from emission_construction_functions import create3Dfromscratch, df2numpy, emissieoorzaak_snap, join_map_points
import os
import glob
import shutil

import sys

#sys.path.insert(1, os.path.join(os.path.dirname(os.getcwd()), ''))
from emission_preparation_setting import (spec_name, snaplist, explicit_point, year_start,
                                          outputdir_main_ruisdael_area_csv2netc,
                                          datadir_ruisdael_area_residuals,
                                          targetdir_point_source_plume_processing,
                                          targetdir_ruisdael_area_total_static,
                                          point_file_harm_p_only, grid_params, layer_num_emiss_z )

class Emission_3Dfield_creation:
    def __init__(self):
        self.year = year_start  #if emissions are available only for another year but here you should put a simulation year...
        self.spec_name = spec_name
        self.snaplist = snaplist
        self.datadir = outputdir_main_ruisdael_area_csv2netc
        self.input_points = datadir_ruisdael_area_residuals
        self.input_points_unassign = targetdir_point_source_plume_processing
        self.targetdir = targetdir_ruisdael_area_total_static
        self.crs = 'RD'
        # refinement using point source (alternative method: they can be used explicitly as DALES input (see point_source folder....)
        self.point_source_refinement = not explicit_point
        self.layer_num_emiss_z = layer_num_emiss_z

    def gaussian_weights(self, z, plume_height, sigma):
        weights = np.exp(-0.5 * ((z - plume_height) / sigma) ** 2)
        return weights / weights.sum()

    def check_emission_dataframe(self, df):
        print("--------------------------------------------------")
        print("GLOBAL EMISSION DATA CHECK")
        print("--------------------------------------------------")

        total = len(df)

        print(f"TOTAL ROWS: {total}")

        # SNAP check
        if "SNAP" in df.columns:
            snap_missing = df["SNAP"].isna().sum()
            print(f"MISSING SNAP: {snap_missing}")

            snap_counts = df["SNAP"].value_counts(dropna=False)
            print("\nSNAP DISTRIBUTION:")
            print(snap_counts)

        else:
            print("WARNING: SNAP column missing!")

        # TYPE check (P / O / others)
        if "TYPE" in df.columns:
            type_counts = df["TYPE"].value_counts(dropna=False)

            print("\nTYPE DISTRIBUTION:")
            print(type_counts)

            n_p = type_counts.get("P", 0)
            n_o = type_counts.get("O", 0)
            n_other = total - n_p - n_o

            print("\nSUMMARY:")
            print(f"TYPE P: {n_p}")
            print(f"TYPE O: {n_o}")
            print(f"TYPE OTHER: {n_other}")

        else:
            print("WARNING: TYPE column missing!")

        # coordinate sanity check
        coord_cols = ["XCO_EMISSIEPUNT_HARM", "YCO_EMISSIEPUNT_HARM"]
        missing_coords = 0

        if all(c in df.columns for c in coord_cols):
            missing_coords = df[coord_cols].isna().any(axis=1).sum()
            print(f"\nROWS WITH MISSING COORDINATES: {missing_coords}")

        # emission sanity check
        if "EMISSIE" in df.columns:
            bad_emis = (~df["EMISSIE"].apply(lambda x: pd.notna(x) and (x >= 0))).sum()
            print(f"INVALID EMISSION VALUES: {bad_emis}")

        print("--------------------------------------------------")
                

    def process_emissions(self):
        
        if not self.point_source_refinement:
            #read unasigned points:
            df_selection_rem_full = pd.read_csv(self.input_points_unassign + 'point_source_unassigned.csv')
            self.check_emission_dataframe(df_selection_rem_full)
        else:
            df_all = pd.read_csv(self.input_points + point_file_harm_p_only, delimiter=',', decimal=".", encoding="ISO-8859-1")
            self.check_emission_dataframe(df_all)

        for isnap in self.snaplist:
            
            
            if isnap == 1:
                icatname = 'power'
            elif isnap == 2:
                icatname = 'residential_commercial'
                print(self.datadir + f'SNAP category {isnap} has no (mostly) vertical component (yet): emission from first model level only')
            elif isnap == 3:
                icatname = 'industrial_combustion'
            elif isnap == 4:
                icatname = 'industrial_processes'
            elif isnap == 5:
                icatname = 'fossil_fuels'
                print(self.datadir + f'SNAP category {isnap} has no vertical component (yet): emission from first model level only')
            elif isnap == 7:
                print(self.datadir + f'SNAP category {isnap} has no vertical component (yet): emission from first model level only')
                icatname = 'road'
            elif isnap == 8: 
                icatname = 'mobile'
            elif isnap == 9:
                icatname = 'waste'
            elif isnap == 10:
                print(self.datadir + f'SNAP category {isnap} has no vertical component (yet): emission from first model level only')
                icatname = 'agriculture'
            elif isnap == 11:
                icatname = 'aviation'

            if isnap in [6]:
                if isnap == 6:
                    print('SNAP category 6 (SOLVENT USE) not to be processed.')
                continue
                
            if not os.path.exists(self.targetdir):
                # Create a new directory because it does not exist:
                os.makedirs(self.targetdir)
                print(f'The {self.targetdir} directory is created!')
            else:
                print(f'The {self.targetdir} directory exists!')
                for f in glob.glob(self.targetdir + f'{self.spec_name}_{self.year}_{isnap}_{icatname}.nc'):
                    os.remove(f)

            # --- 1 --- 2D file

            # Load 2D input emission data
            fobj = netc.Dataset(self.datadir + f'HARM_snap_{isnap}_all_{spec_name}.nc')
            x = fobj.variables['x'][:]
            y = fobj.variables['y'][:]
            speci_2d = fobj.variables[self.spec_name][:].T  # Transpose if needed (Y,X)

            # Fill NaNs with zero (DALES does not accept NaNs)
            speci_2d = np.nan_to_num(speci_2d)

            # --- Vertical grid setup ---
            nz = self.layer_num_emiss_z  # e.g. 5
            dz0 = grid_params['dz0']                    # base layer thickness
            alpha = grid_params['alpha']                # stretching factor

            dz = dz0 * (1 + alpha) ** np.arange(nz)
            zh = np.zeros(nz + 1)
            zh[1:] = np.cumsum(dz)
            z = 0.5 * (zh[:-1] + zh[1:])
            zsize = zh[-1]

            # Initialize 3D emission array (shape: Z, Y, X)
            speci_3d = np.zeros((nz, *speci_2d.shape))


            # Inside process_emissions for each SNAP category after loading speci_2d and vertical grid z:
            #Based on esrtimations from:
            #Maas, C. J., & Winther, M. (2007). Stack heights for point source emissions in Europe. Atmospheric Environment.
            #Zheng et al. (2019), "Review of plume rise algorithms and application in atmospheric modeling".

            if isnap in [2, 5, 7, 10]:  # mostly ground-level categories
                speci_3d[0, :, :] = speci_2d

            elif isnap in [1, 3]:  # power and industrial combustion with plume rise
                plume_height = 80  # plume rise height in meters (adjust as needed)
                sigma = 20         # vertical spread of plume
                weights = self.gaussian_weights(z, plume_height, sigma)
                for k in range(nz):
                    speci_3d[k, :, :] = speci_2d * weights[k]

            elif isnap == 4:  # industrial processes - if no info, distribute equally or gaussian
                plume_height = 60  # lower plume rise, for example
                sigma = 15
                weights = self.gaussian_weights(z, plume_height, sigma)
                for k in range(nz):
                    speci_3d[k, :, :] = speci_2d * weights[k]

            elif isnap == 8:  # mobile sources (including shipping) with vertical component
                plume_height = 30  # example average ship stack height in meters
                sigma = 10         # spread around plume height
                weights = self.gaussian_weights(z, plume_height, sigma)
                for k in range(nz):
                    speci_3d[k, :, :] = speci_2d * weights[k]

            elif isnap == 11:      # Aviation (LTO)

                plume_height = 50.0      # m
                sigma = 30.0             # m

                weights = self.gaussian_weights(z, plume_height, sigma)
                for k in range(nz):
                    speci_3d[k, :, :] = speci_2d * weights[k]
            else:
                # fallback: equal distribution
                for k in range(nz):
                    speci_3d[k, :, :] = speci_2d / nz

            #if isnap in [5, 7, 10]:
                # Place all emissions at lowest level only
                #speci_3d[0, :, :] = speci_2d
            #else:
                # Distribute emissions equally across all vertical layers
                #for k in range(nz):
                    #speci_3d[k, :, :] = speci_2d / nz

            # --- 2 --- Point Data (if point source emissions used as a refinement)
                            
            # --- fill all remaining NaNs with zero before continue, since DALES does not accept nans in input----
            speci = np.nan_to_num(speci_3d)

            speci = np.transpose(speci, (1, 2, 0))  # shape: (x, y, z)
            
            df_all_chunk = pd.DataFrame()
            
            x_min, x_max = x.min(), x.max()
            y_min, y_max = y.min(), y.max()
            
            if self.point_source_refinement:
                print('The refinement of area emiss with point sources is processing.... ') 
                
                # Select all rows for current SNAP category
                df_all_chunk = df_all[df_all['SNAP'] == isnap].copy()

                # Apply spatial mask
                mask = (
                    (df_all_chunk['XCO_EMISSIEPUNT_HARM'] >= x_min) &
                    (df_all_chunk['XCO_EMISSIEPUNT_HARM'] <= x_max) &
                    (df_all_chunk['YCO_EMISSIEPUNT_HARM'] >= y_min) &
                    (df_all_chunk['YCO_EMISSIEPUNT_HARM'] <= y_max)
                )

                num_points_within_domain = mask.sum()
                print(f"Total points in df_all_chunk: {len(df_all_chunk)}")
                print(f"Points within LES grid domain: {num_points_within_domain}")

                missing_points = df_all_chunk[~mask]
                print(f"Points outside the LES grid: {len(missing_points)}")
                print(missing_points[['XCO_EMISSIEPUNT_HARM', 'YCO_EMISSIEPUNT_HARM']].head(10))

                # Only map points if there is something to map
                if not df_all_chunk.empty:
                    join_map_points(df_all_chunk, speci,
                    x, y, z,
                    mode='add',
                    dim=3,
                    xloc='XCO_EMISSIEPUNT_HARM',
                    yloc='YCO_EMISSIEPUNT_HARM')

            else:
                
                print('Point sources will be used as an explicit model input using separate nc-files, so here we add only point sourses with unknown plume parameters.... ')
                
                speci = np.array(speci.data)
                
                df_point_subdoelgroep = df_selection_rem_full[df_selection_rem_full['SNAP'] == isnap]
                
                #Add point sources, which were not included to explicit model input files:
                # Only map points if there is something to map
                if not df_point_subdoelgroep.empty:
                    join_map_points(df_point_subdoelgroep, speci,
                    x, y, z,
                    mode='add',
                    dim=3,
                    xloc='XCO_EMISSIEPUNT_HARM',
                    yloc='YCO_EMISSIEPUNT_HARM')
                
                
                df_point_subdoelgroep = pd.DataFrame()
                
                
        
            speci = np.transpose(speci, (2, 1, 0))  # back to (z, y, x)    

            fobj = netc.Dataset(self.targetdir + f'{self.spec_name}_{self.year}_{isnap}_{icatname}.nc', 'w')

            fobj.description = f"{self.spec_name.upper()} emissions (year {self.year}, SNAP category {isnap})"

            dim_x = fobj.createDimension('x', len(x))
            dim_y = fobj.createDimension('y', len(y))
            dim_z = fobj.createDimension('z', len(z))

            var_x = fobj.createVariable('x', 'f4', ('x',))
            var_y = fobj.createVariable('y', 'f4', ('y',))
            var_z = fobj.createVariable('z', 'f4', ('z',))
            var_e = fobj.createVariable(self.spec_name.upper(), 'f8', ('z', 'y', 'x',))

            var_x[:] = np.array(x)
            var_y[:] = np.array(y)
            var_z[:] = np.array(z)

            var_e[:, :, :] = speci

            if self.crs == 'RD':
                var_x.units = 'Rijksdriehoekcoordinaat x in meters'
                var_y.units = 'Rijksdriehoekcoordinaat y in meters'
            elif self.crs == 'HARM':
                var_x.units = 'Harmonie coordinate x in meters'
                var_y.units = 'Harmonie coordinate y in meters'

            var_z.units = 'Altitude z in meters'
            var_e.units = 'Kilogram per year'

            fobj.close()

if __name__ == "__main__":
    E_3D = Emission_3Dfield_creation()
    E_3D.process_emissions()
