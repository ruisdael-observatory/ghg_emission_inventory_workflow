"""
Module description:
Gapfilling of Emissieregistratie point sources
and creation of explicit point source input for DALES 

INPUT  (CSV) Emissiregistratie raw point source data

This file includes columns like DATASET,
EMISSIEJAAR, STOF, CODE_EMISSIEOORZAAK, EMISSIEOORZAAK, CODE_BEDRIJF,
NAAM_BEDRIJF, XCO_BEDRIJF, YCO_BEDRIJF, EMISSIEPUNT_CODE,
EMISSIEPUNT_NAAM, XCO_EMISSIEPUNT, YCO_EMISSIEPUNT,
TYPE, HOOGTE, LENGTE, BREEDTE, HOEK, UITSTROOMOPENING_M2,
TEMPERATUUR, UITTREEDSNELHEID, VOLUMESTROOM, WARMTEINHOUD,
EMISSIE, EENHEID, XCO_EMISSIEPUNT_HARM, YCO_EMISSIEPUNT_HARM,
XCO_BEDRIJF_HARM, YCO_BEDRIJF_HARM, delimited by a comma (,).

OUTPUT (netCDF) Gapfilled point source input ncdf files per processor block and per hour
Code updated with a loop over the simulation period to prepare point sources for the whole period

Author: Dr. Arseni Doyennel
Email: a.doyennel@vu.nl
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from emission_construction_functions import (loadsnap, emissieoorzaak_snap, regrescat, gapfill,
                                             write_df_to_csv, df2list, data2netc_old, data2netc)
import datetime as datetime
import os
import glob


from emission_preparation_setting import (year_start, month_start, day_start,
                                            hour_start, year_end, month_end, day_end, hour_end,
                                            x0, y0, xres, yres, xsize, ysize, itot, jtot, ktot, nprocx, nprocy,
                                            spec_name, point_source_harm_file,
                                            datadir_point_source_plume_processing,
                                            targetdir_point_source_plume_processing, point_file_harm_p_only )


class Pointsource_input_preparation:
    def __init__(self):

        self.datadir = datadir_point_source_plume_processing
        self.targetdir = targetdir_point_source_plume_processing
        self.pointsourcefile = point_source_harm_file
        self.year_sim = year_start  # if emissions are used for other year, put the needed year in year_start of setting file!
        self.month_start = month_start
        self.day_start = day_start
        self.hour_start = hour_start
        self.year_end = year_end
        self.month_end = month_end
        self.day_end = day_end
        self.hour_end = hour_end
        self.minute = 0
        self.x0 = x0
        self.y0 = y0
        self.xres = xres
        self.yres = yres
        self.xsize = xsize
        self.ysize = ysize
        self.itot = itot
        self.jtot = jtot
        self.ktot = ktot
        self.nprocx = nprocx
        self.nprocy = nprocy
        self.spec_name = spec_name
        # === Script
        self.DoPlot = True
        self.DoVerbose = True
        self.DoWrite = True
        # === Regression (parameters in function gapfill (see function script)):
        self.t_order = 0
        self.v_order = 1
        self.h_order = 1
        self.h_shape = 'log'
        


    def create_target_directory(self):
        if not os.path.exists(self.targetdir):
            os.makedirs(self.targetdir)
            print(f'The {self.targetdir} directory is created!')
        else:
            print(f'The {self.targetdir} directory exists!')

            for f in glob.glob(self.targetdir + '*'):
                os.remove(f)


    def prepare_time_array(self):

        time_array = pd.Series(
            pd.date_range(
                start=f'{self.year_sim}-{self.month_start}-{self.day_start} {self.hour_start}:00:00',
                end=f'{self.year_end}-{self.month_end}-{self.day_end} {self.hour_end}:00:00',
                freq='H',
            )
        )

        return time_array

    def load_raw_emission_data(self):
        #return pd.read_csv(
            #self.datadir + self.pointsourcefile,
            #delimiter=',',
            #decimal=".",
            #encoding="ISO-8859-1",
        #)
    
        return pd.read_csv(
            self.datadir + self.pointsourcefile,
            delimiter=';', 
            low_memory=False
            
        )



    def identify_emission_categories(self):
        # === Identify number of "Emissieoorzaken" and add numerical label
        # Note: if you've got "Unknown emissieoorzaak", check emissieoorzaak_snap function,
        #maybe some categories are missing and must be added

        df = self.load_raw_emission_data()
        df['EMISSIEOORZAAKLABEL'] = 0
        df['SNAP'] = 0
        df_orig = df.copy()

        emissieoorzaken = sorted(df.EMISSIEOORZAAK.unique().tolist())

        for irow in range(len(df)):
            df.loc[irow, 'EMISSIEOORZAAKLABEL'] = emissieoorzaken.index(
                df.loc[irow, 'EMISSIEOORZAAK']
            )
            df.loc[irow, 'SNAP'] = emissieoorzaak_snap(df.loc[irow, 'EMISSIEOORZAAK'])

        return df, df_orig

    # NOTE: Qualitative grouping of emission causes suitable for regression ===============================

    def prepare_emission_data(self, df):
        # Code for preparing emission data
        
        # Initialize df_selection_non_rem_full and df_gapfilled_full as empty DataFrames:
        df_selection_non_rem_full = pd.DataFrame()
        df_gapfilled_full = pd.DataFrame()
        
        df_selection_rem_full = pd.DataFrame()
        
        #Create a target dir:
        self.create_target_directory()

        # Create an empty DataFrame
        df_empty = pd.DataFrame(columns=[])

        # Save the empty DataFrame to a CSV file
        #Remaining points for which unsufficient cohesive data is available for regression are put together 
        #into point_source_unassigned.csv to be used in the area emissions (see static module!)
        
        df_empty.to_csv(self.targetdir + 'point_source_unassigned.csv', index=False)

        #these groups will be concidered in point sources, all the rest will be put in unassigned and used in area emissions:
        #probably, the reason why PS from only these groups are considered is that point sources from these groups 
        #are originate from specific sources like chimneys, exhaust stacks, or vents, 
        #making them distinct from the more diffuse emissions associated with other industrial processes...
        
        #Commented emissieoorzaakgroepen are for the year 2018 (in other years, they may change)...
        emissieoorzaakgroepen = [
            #['delfstoffen', [2]],
            ['voed', [3, 4, 5, 6, 7, 8, 9, 10]],
            ['papier', [13, 14, 15]],
            ['olie', [16]],
            ['chem', [17, 18, 19, 20, 21, 22, 23]],
            ['verf', [25]],
            ['kunststof', [33]],
            ['glas/keramiek', [34, 35, 36]],
            #['bakstenen', [37]],
            ['metaal', [40, 41, 42, 43, 44, 45, 46]],
            ['gieten', [47]],
            #['metaal_prod', [48]],
            ['elektro', [49]],
            ['scheepsbouw', [53]],
            ['elektriciteit', [55, 56]],
            ['gas', [57]],
            ['avi', [59]],
            ['bouw', [62]],
            ['remaining', []]    
        ]

        remaining_idxs = list(range(len(df.EMISSIEOORZAAK.unique())))

        for grp in emissieoorzaakgroepen:
            for idx in grp[1]:
                remaining_idxs.remove(idx)

        emissieoorzaakgroepen[-1][1] = remaining_idxs

        n_replaced_v, n_replaced_t, n_replaced_h = 0, 0, 0
        # Loop over emissieoorzaakgroepen
        for igroep in emissieoorzaakgroepen:

            df_selection = df[
                (df['EMISSIEOORZAAKLABEL'].isin(list(igroep[1])))
                & (df['EMISSIE'] > 0)
                & (df['TYPE'] == 'P') #We consider only point sources here, I suppose, since in area emissions, to get residuals we substract only emissions with P-type, i.g., point sources...
                ]

            if igroep[0] != 'remaining':
                
                # Concatenate df_selection with df_selection_non_rem_full:
                df_selection_non_rem_full = pd.concat([df_selection_non_rem_full, df_selection], ignore_index=True)
                
                df_gapfilled = df_selection.copy()
                # Apply regression to existing values
                replace_v, replace_t, replace_h = gapfill(df_gapfilled)

                # Use regression to fill missing values
                for replace in replace_v:
                    df_gapfilled.loc[replace[0], 'VOLUMESTROOM'] = replace[1]
                for replace in replace_t:
                    df_gapfilled.loc[replace[0], 'TEMPERATUUR'] = replace[1]
                for replace in replace_h:
                    df_gapfilled.loc[replace[0], 'HOOGTE'] = replace[1]

                n_replaced_v += len(replace_v)
                n_replaced_t += len(replace_t)
                n_replaced_h += len(replace_h)
                
                # Assign df_gapfilled to df_test
                df_test = df_gapfilled.copy()
                
       
                # Concatenate df_selection with df_gapfilled_full:
                df_gapfilled_full = pd.concat([df_gapfilled_full, df_gapfilled], ignore_index=True)

                if self.DoPlot:
                    fig, ax = plt.subplots(1, 3, figsize=(10, 6))
                    plotvars = ['VOLUMESTROOM', 'TEMPERATUUR', 'HOOGTE']
                    for ivar, plotvar in enumerate(plotvars):
                        ax[ivar].scatter(
                            x=df_test['EMISSIE'], y=df_test[plotvar], c='crimson', s=10
                        )
                        ax[ivar].scatter(
                            x=df_selection['EMISSIE'],
                            y=df_selection[plotvar],
                            color='royalblue',
                            s=30,
                        )
                        ax[ivar].set_title(plotvar)
                        ax[ivar].set_xscale('log')

                        if plotvar == 'VOLUMESTROOM':
                            ax[ivar].set_yscale('log')
                        elif plotvar == 'TEMPERATUUR':
                            ax[ivar].set_ylim(bottom=10)

                    fig.suptitle(igroep[0])
                    plt.show()
            else:
                
                #NOTE: in some of remaining emissions (with P-type) there are still values of 
                #VOLUMESTROOM, TEMPERATUUR, HOOGTE, and EMISSIE which all >0, 
                #but since they are in remaining groups they are filtered out above. 
                #In my perspective, they should be added:
                
                #Check if in 'remaining', there are point sources which have all conditions to be included
                #in df_gapfilled_full, since they have no gaps:
                
                df_points_also_to_include = df_selection[
                        (df_selection['VOLUMESTROOM'] > 0) &
                        (df_selection['TEMPERATUUR'] > 0) &
                        (df_selection['HOOGTE'] > 0) &
                        (df_selection['EMISSIE'] > 0)].copy()  # Create a copy of the DataFrame
                
                  
                #Concatenate df_selection with df_gapfilled_full:
                df_gapfilled_full = pd.concat([df_gapfilled_full, df_points_also_to_include], ignore_index=True)
                
                #Remove the selected rows from df_selection:
                df_selection = df_selection.drop(df_points_also_to_include.index)
             
                #Concatenate df_selection with df_selection_rem_full:    
                df_selection_rem_full = pd.concat([df_selection_rem_full, df_selection], ignore_index=True)
                

                
                #save to point_source_unassigned only those which will be not used in point sources to add them to area emiss:
                df_selection_rem_full.to_csv(self.targetdir + 'point_source_unassigned.csv', index=False)
        
        
        print(
            f'Filled values, Volumestroom:{n_replaced_v}, Temperatuur:{n_replaced_t}, Hoogte:{n_replaced_h}'
        )
        


        return df_selection_rem_full, df_gapfilled_full
    

    def plot_emission_data(self, df, df_orig): #use here df_selection_non_rem_full, df_gapfilled_full
        # Code for plotting emission data
        fig, ax = plt.subplots(1, 3, figsize=(12, 6))
        plotvars = ['VOLUMESTROOM', 'TEMPERATUUR', 'HOOGTE']
        for ivar, plotvar in enumerate(plotvars):
            ax[ivar].scatter(x=df['EMISSIE'], y=df[plotvar], color='royalblue', s=20)
            ax[ivar].scatter(x=df_orig['EMISSIE'], y=df_orig[plotvar], c='crimson', s=10)
            ax[ivar].set_title(plotvar)
            ax[ivar].set_xscale('log');
            ax[ivar].set_xlim([1e3, 1e10])

            if plotvar == 'VOLUMESTROOM':
                ax[ivar].set_yscale('log')
                ax[ivar].set_ylim([1e-2, 1e4])
            elif plotvar == 'TEMPERATUUR':
                ax[ivar].set_ylim(bottom=10)

        plt.show()

        
        #Prepare final model input [per hour] for each hour of simulaiton period:

    def prepare_final_input(self, df):
        #A loop over time to create emissions for each hour of simulation period per block:
        #here, all point sourses from df are selected from the full list: data = df2list(df,xmin, xmax, ymin, ymax, dx, dy,1)
        # ============================================================================
        print(f"Point sources input ncdf files are preparing...")

        time_array = self.prepare_time_array()   #get time_array
                
        
        
        for timepoint in range(self.hour_start, time_array.size, 1):
            year = time_array[timepoint].year
            month = time_array[timepoint].month
            day = time_array[timepoint].day
            hour = time_array[timepoint].hour

            tprof_hour, tprof_week, tprof_mnth = loadsnap()

            imonth = month - 1 #-1 because in python counting starts from 0
            ihour = hour
            iweek = datetime.datetime(self.year_sim, month, day, hour).weekday()

            blockx = self.xsize / self.nprocx
            blocky = self.ysize / self.nprocy
            npoints = 0
            pointlist = []

            for mpiidx in range(0, self.nprocx):
                for mpiidy in range(0, self.nprocy):

                    xmin = self.x0 + mpiidx * blockx
                    xmax = xmin + blockx
                    dx = self.xsize / self.itot

                    ymin = self.y0 + mpiidy * blocky
                    ymax = ymin + blocky
                    dy = self.ysize / self.jtot

                    data = df2list(df, xmin, xmax, ymin, ymax, dx, dy, 1)

                    npoints += len(data)
                    pointlist.append(len(data))

                    data[:, 4] = data[:, 4] / (365 * 24) #emissions per hour!

                    if len(data) > 0:

                        #Temporal aggregation is applied to account for changes in the emission strength of a source
                        #for a certain hour in a day, day in the week, or month of the day,
                        #e.g. traffic rush hour or energy needed for heating household in summer vs. winter.
                        #We apply scaling factors from TNO (Denier van der Gon, 2012)
                        #which is based on SNAP (Standard Nomenclature for Air Pollution) emission sectors:

                        for source in data:
                            isnap = int(source[-1] - 1)
                            source[4] = (
                                    source[4]
                                    * tprof_hour[isnap, ihour]
                                    * tprof_week[isnap, iweek]
                                    * tprof_mnth[isnap, imonth]
                            )
                            
                        data

                        data2netc(
                            data,
                            self.targetdir,
                            nprocx=mpiidx,
                            nprocy=mpiidy,
                            tracer=spec_name,
                            minute=0,
                            hour=hour,
                            day=day,
                            month=month,
                            year=self.year_sim,
                        )
                        print(
                            f"Point sources input is prepared for: {self.year_sim}-{month}-{day} {hour}:00:00"
                        )



if __name__ == "__main__":
    PSEP = Pointsource_input_preparation()
    df, df_org = PSEP.identify_emission_categories()
    df_selection_rem_full, df_gapfilled_full = PSEP.prepare_emission_data(df)
    #PSEP.plot_emission_data(df_selection_non_rem_full, df_gapfilled_full)
    PSEP.prepare_final_input(df_gapfilled_full)
