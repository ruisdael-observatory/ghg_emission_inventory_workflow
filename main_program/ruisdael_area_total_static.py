"""
Module Description:

This script is designed to generate static annual total area emissions within a specified domain.
Additionally, it incorporates the following functionalities:

Refinement of area emissions using point source data (if point sources are not utilized as explicit input).
Generation of 3-D emission files per SNAP category to facilitate the final preparation of area emission inputs.
( area emission are devided by number of vertical layers, covered with plume until the emission top height layer!)
On vertical allocation of emissions: currently, the vertical allocation of area emissions 
is based on results presented in Brunner et al. (2019) (but applying even destribution, yet) [https://acp.copernicus.org/articles/19/4541/2019/].

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
from emission_preparation_setting import (spec_name, snaplist, year_start,
                                          outputdir_main_ruisdael_area_csv2netc,
                                          datadir_ruisdael_area_residuals,
                                          targetdir_point_source_plume_processing,
                                          targetdir_ruisdael_area_total_static,
                                          point_file_harm_p_only)

class Emission_3Dfield_creation:
    def __init__(self):
        self.year = year_start  #if emissions are available only for another year but here you should put a simulation year...
        self.spec_name = spec_name
        self.snaplist = snaplist
        self.datadir = outputdir_main_ruisdael_area_csv2netc
        self.input_points = datadir_ruisdael_area_residuals
        self.input_points_unassign = targetdir_point_source_plume_processing
        self.targetdir = targetdir_ruisdael_area_total_static
        self.crs = 'HARM'
        # refinement using point source (alrernative method: they can be used explicitly as DALES input (see point_source folder....)
        self.point_source_refinement = False
        self.layer_num_emiss_z = 6  #IMPORTANT parameter: if SNAP category has vertical component, emiss will be placed at some number of model layers (depends on zmax) from the emission bottom to the emission top
                                    #currently, the plume top height is set as ~150m, so emiss should be devided by number of vertical layers from the fround to 150m
                                    # where emiss will be placed, to make equal destribution and observe the mass conservation.
                                    # If you change emiss top height, check the layers and change this number!
                                    #NOTE: for point sources, if you use explicit input method,
                                    # plume rise height is calculated from stack height interactively; preliminary
                                    # the number of vertical layers, where emissions are allocated, is unkown and
                                    # the emiss will be destributed directly in the model code !
                                    #TODO: now, emiss are equally destributed between vertical layers until plume top height, but in the nature, there might be a variability in vertical concentration...

                

    def process_emissions(self):
        
        if not self.point_source_refinement:
            #read unasigned points:
            df_selection_rem_full = pd.read_csv(self.input_points_unassign + 'point_source_unassigned.csv', delimiter=',', decimal=".", encoding="ISO-8859-1")
        else:
            df_all = pd.read_csv(self.input_points + point_file_harm_p_only, delimiter=',', decimal=".", encoding="ISO-8859-1")
            
         
        
        for isnap in self.snaplist:
            
            
            if isnap == 1:
                icatname = 'power'
            elif isnap == 2:
                icatname = 'residential_commercial'
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

            fobj = netc.Dataset(self.datadir + f'HARM_snap_{isnap}_all.nc')
            x = fobj.variables['x'][:]
            y = fobj.variables['y'][:]
            speci = fobj.variables[self.spec_name][:]

            nz = self.layer_num_emiss_z
            dz0 = 25
            alpha = 0.017

            dz = np.zeros(nz)
            z = np.zeros(nz)
            zh = np.zeros(nz + 1)

            dz[:] = dz0 * (1 + alpha) ** np.arange(nz)
            zh[1:] = np.cumsum(dz)
            z[:] = 0.5 * (zh[1:] + zh[:-1])
            zsize = zh[-1]

            # --- 2 --- Point Data (if point source emissions used as a refinement)
                            
            # --- fill all remaining NaNs with zero before continue, since DALES does not accept nans in input----
            speci = np.nan_to_num(speci)
            
            df_all_chunk = pd.DataFrame()
            
            if self.point_source_refinement:
                print('The refinement of area emiss using point sources is processing.... ')
                
                #Add here the name of your point source file:
                speci = np.array(speci.data)

                for i, item in enumerate(df_all['EMISSIEOORZAAK'].unique()):
                    snap_cat = emissieoorzaak_snap(item)
        
                    
                    if snap_cat == isnap:
                        
                        df_all_chunk = pd.concat([df_all_chunk, df_all[df_all['EMISSIEOORZAAK'] == item]], ignore_index=True)
                        
                
                #Add point sources, which were early substracted, to refine area residuals:
                if not df_all_chunk.empty:
                    join_map_points(df_all_chunk, speci,
                                            x, y,
                                            mode='add',
                                            dim=2,
                                            xloc='XCO_EMISSIEPUNT_HARM',
                                            yloc='YCO_EMISSIEPUNT_HARM')
                               


            else:
                
                print('Point sources will be used as an explicit model input using separate nc-files.... ')
                
                speci = np.array(speci.data)
                
                df_point_subdoelgroep = df_selection_rem_full[df_selection_rem_full['SNAP'] == isnap]
                
                #Add point sources, which were not included to explicit model input files:
                join_map_points(df_point_subdoelgroep, speci,
                                        x, y,
                                        mode='add',
                                        dim=2,
                                        xloc='XCO_EMISSIEPUNT_HARM',
                                        yloc='YCO_EMISSIEPUNT_HARM')
               
                  

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

            if isnap in [5, 7, 10]:

                var_e[:, :, :] = speci.T * 0 #do not need to devide by number of layers, since here all emiss are at first layer
                var_e[0, :, :] = speci.T
            else:
                var_e[:, :, :] = speci.T/self.layer_num_emiss_z   #devide by number of vertical layers, where emiss are alocated (check model layers +1)

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
