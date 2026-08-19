"""
Module description:
Gapfilling of Emissieregistratie point sources
and creation of explicit point source input for DALES 

INPUT  (CSV) Emissiregistratie raw point source data

Note: explicit point source input can be used only if 
HOOGTE, UITSTROOMOPENING_M2, TEMPERATUUR, VOLUMESTROOM are available in input file!

OUTPUT (netCDF) Gapfilled hourly point source input ncdf files
Code updated with a loop over the simulation period to prepare point sources for the whole period

Author: Dr. Arseni Doyennel
Email: a.doyennel@vu.nl
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from emission_construction_functions import (get_local_time, emissieoorzaak_snap, regrescat, gapfill,
                                             write_df_to_csv, df2list_new, data2netc_point )
import temp_disaggregation_new_cams_temp_prof

import datetime as datetime
import os
import glob
import statsmodels.api as sm
from sklearn.preprocessing import PolynomialFeatures
from netCDF4 import Dataset
from projection import Transform
from GridDales import GridDales
from emission_preparation_setting import (year, year_start, month_start, day_start,
                                            hour_start, year_end, month_end, day_end, hour_end,
                                            x0, y0, xres, yres, xsize, ysize, itot, jtot, ktot, nprocx, nprocy,
                                            spec_name, point_source_harm_file,
                                            datadir_point_source_plume_processing,
                                            targetdir_point_source_plume_processing, point_file_harm_p_only,
                                            proj_params, grid_params, x_offset, y_offset, noxdir, temp_prof_file, author, email )


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
        self.x_offset = x_offset
        self.y_offset = y_offset
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
        self.a_order = 1
        self.h_shape = 'log'
        self.area_shape = 'log'
        self.proj_params = proj_params
        self.grid_params = grid_params
        # Initialize transformer and grid
        self.transformer = Transform(self.proj_params)
        self.les_grid = GridDales(self.grid_params)
        self.lcc_start_x, self.lcc_start_y = x0+self.x_offset, y0+self.y_offset
        # Update grid southwest corner
        self.grid_params['southwest_x'] = self.lcc_start_x
        self.grid_params['southwest_y'] = self.lcc_start_y
        
        


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
                freq='h',
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

        #Fill missing emission location with location of the company (save fallback to not lost emissions):
        missing_xy = (
            df["XCO_EMISSIEPUNT_HARM"].isna() |
            df["YCO_EMISSIEPUNT_HARM"].isna()
        )

        df.loc[missing_xy, "XCO_EMISSIEPUNT_HARM"] = df.loc[missing_xy, "XCO_BEDRIJF"]
        df.loc[missing_xy, "YCO_EMISSIEPUNT_HARM"] = df.loc[missing_xy, "YCO_BEDRIJF"]
        df.loc[missing_xy, "XCO_EMISSIEPUNT"] = df.loc[missing_xy, "XCO_BEDRIJF"]
        df.loc[missing_xy, "YCO_EMISSIEPUNT"] = df.loc[missing_xy, "YCO_BEDRIJF"]
        
        print(df.columns)
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
        
        # Initialize unassigned_points and df_gapfilled_full as empty DataFrames:
        df_gapfilled_full = pd.DataFrame()
        unassigned_points = pd.DataFrame()
        
        #Create a target dir:
        self.create_target_directory()

        # Create an empty DataFrame
        df_empty = pd.DataFrame(columns=[])

        # Before starting the loop: add a helper unique identifier if one doesn't exist
        df = df.reset_index(drop=False).rename(columns={"index": "ORIGINAL_INDEX"})

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
            ['delfstoffen', [2]],
            ['voed', [3, 4, 5, 6, 7, 8, 9, 10]],
            ['papier', [13, 14, 15]],
            ['olie', [16]],
            ['chem', [17, 18, 19, 20, 21, 22, 23]],
            ['kunststof', [33]],
            ['glas/keramiek', [34, 35, 36]],
            ['bakstenen', [37]],
            #['metaal', [40, 41, 42, 43, 44, 45, 46]],
            ['metaal_prod', [48]],
            ['scheepsbouw', [53]],
            ['elektriciteit', [55, 56]],
            ['gas', [57]],
            ['avi', [59]],
            ['bouw', [62]],
            ['remaining', []]
        ]

        # Use EMISSIEOORZAAKLABEL, since you filter on this column
        remaining_idxs = df.EMISSIEOORZAAKLABEL.unique().tolist()

        n_replaced_v = 0
        n_replaced_t = 0
        n_replaced_h = 0
        n_replaced_a = 0

        for grp in emissieoorzaakgroepen[:-1]:  # exclude last group for now
            for idx in grp[1]:
                if idx in remaining_idxs:
                    remaining_idxs.remove(idx)

        emissieoorzaakgroepen[-1][1] = remaining_idxs

        print("Remaining labels assigned to last group:", remaining_idxs)

        # Now print sizes again
        total_count = 0
        for igroep in emissieoorzaakgroepen:
            group_name = igroep[0]
            group_indices = igroep[1]
            df_selection = df[
                (df['EMISSIEOORZAAKLABEL'].isin(group_indices))
                & (df['EMISSIE'] > 0)
                & (df['TYPE'] == 'P')
            ]
            print(f"{group_name}: {df_selection.shape}")
            total_count += df_selection.shape[0]

            print(f"Total points summed across groups: {total_count}")
            
            
            
            n_samples = len(df_selection)

            # Print the count of samples for the current group
            print(f"Emissieoorzaakgroep '{group_name}' has {n_samples} samples in the data.")   
            

            # Define the minimum number of samples required for gap filling
            min_samples_for_gapfill = 5  # Set this to the appropriate threshold....
            
            if n_samples <= min_samples_for_gapfill:
                
                 print(f"Emissieoorzaakgroep '{group_name}' has {n_samples} samples --> no enough for regression!")

            # Check if the group is not 'remaining' and has enough samples for gap filling
            if igroep[0] != 'remaining' and n_samples > min_samples_for_gapfill:
                df_gapfilled = df_selection.copy()
                replace_v, replace_t, replace_h, replace_a  = gapfill(df_gapfilled)

                for replace in replace_v:
                    df_gapfilled.loc[replace[0], 'VOLUMESTROOM'] = replace[1]
                for replace in replace_t:
                    df_gapfilled.loc[replace[0], 'TEMPERATUUR'] = replace[1]
                for replace in replace_h:
                    df_gapfilled.loc[replace[0], 'HOOGTE'] = replace[1]
                for replace in replace_a:
                    df_gapfilled.loc[replace[0], 'UITSTROOMOPENING_M2'] = replace[1]

                n_replaced_v += len(replace_v)
                n_replaced_t += len(replace_t)
                n_replaced_h += len(replace_h)
                n_replaced_a += len(replace_a)

                # ✅ Only keep valid rows
                valid_rows = df_gapfilled[
                    (df_gapfilled['VOLUMESTROOM'] > 0) &
                    (df_gapfilled['TEMPERATUUR'] > 0) &
                    (df_gapfilled['HOOGTE'] > 0) &
                    (df_gapfilled['UITSTROOMOPENING_M2'] > 0) &
                    (df_gapfilled['EMISSIE'] > 0)
                ]
                df_gapfilled_full = pd.concat([df_gapfilled_full, valid_rows], ignore_index=True)


                invalid_rows = df_gapfilled.loc[~df_gapfilled.index.isin(valid_rows.index)]
                unassigned_points = pd.concat([unassigned_points, invalid_rows], ignore_index=True)

                # Assign df_gapfilled to df_test
                df_test = df_gapfilled.copy()

                if self.DoPlot:
                    fig, ax = plt.subplots(1, 4, figsize=(10, 6))
                    plotvars = ['VOLUMESTROOM', 'TEMPERATUUR', 'HOOGTE', 'UITSTROOMOPENING_M2']
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
                # Step 1: Select extra valid point sources (no gaps)
                df_points_also_to_include = df_selection[
                    (df_selection['VOLUMESTROOM'] > 0) &
                    (df_selection['TEMPERATUUR'] > 0) &
                    (df_selection['HOOGTE'] > 0) &
                    (df_selection['UITSTROOMOPENING_M2'] > 0) &
                    (df_selection['EMISSIE'] > 0)
                ].copy()

                # Step 2: Add these to gapfilled only
                df_gapfilled_full = pd.concat([df_gapfilled_full, df_points_also_to_include], ignore_index=True)

                # Step 3: Remove these rows from df_selection
                df_selection = df_selection[~df_selection['ORIGINAL_INDEX'].isin(df_points_also_to_include['ORIGINAL_INDEX'])]

                # Step 4: Add remaining (unusable) df_selection rows to non-rem group
                unassigned_points = pd.concat([unassigned_points, df_selection], ignore_index=True)

        # Final cleanup
        df_gapfilled_full = df_gapfilled_full.drop_duplicates(subset=['ORIGINAL_INDEX'])
        unassigned_points = unassigned_points.drop_duplicates(subset=['ORIGINAL_INDEX'])

        # Save and return the corrected remainder
        unassigned_points.to_csv(self.targetdir + 'point_source_unassigned.csv', index=False)        
        print(
            f'Filled values, Volumestroom:{n_replaced_v}, Temperatuur:{n_replaced_t}, Hoogte:{n_replaced_h}, , exit area:{n_replaced_h}'
        )


        total_point_sources = len(df[df['TYPE'] == 'P'])
        print("Total point sources in original df:", total_point_sources)

        print("Gapfilled + directly accepted:", len(df_gapfilled_full))
        print("Leftover (unassigned):", len(unassigned_points))
        print("Sum:", len(df_gapfilled_full) + len(unassigned_points))
        


        return unassigned_points, df_gapfilled_full
    

    def plot_emission_data(self, df, df_orig, v_order=1, t_order=0, h_order=1, a_order=1, h_shape='log', a_shape='log'):
        fig, ax = plt.subplots(1, 4, figsize=(15, 5))
        plotvars = ['VOLUMESTROOM', 'TEMPERATUUR', 'HOOGTE' ,'UITSTROOMOPENING_M2']
        log_vars = ['VOLUMESTROOM', 'HOOGTE','UITSTROOMOPENING_M2']  # Variables to be plotted on a log scale
        orders = {'VOLUMESTROOM': v_order, 'TEMPERATUUR': t_order, 'HOOGTE': h_order, 'UITSTROOMOPENING_M2': a_order}
    
        emissieoorzaakgroepen = [
            3, 4, 5, 6, 7, 8, 9, 10, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 25, 33, 34, 35, 36,
            40, 41, 42, 43, 44, 45, 46, 47, 49, 53, 55, 56, 57, 59, 62
        ]
    
        df_filtered = df[df['EMISSIEOORZAAKLABEL'].isin(emissieoorzaakgroepen)]
        df_orig_filtered = df_orig[df_orig['EMISSIEOORZAAKLABEL'].isin(emissieoorzaakgroepen)]
    
        # Loop through the emission cause groups
        for ivar, plotvar in enumerate(plotvars):
            for group in emissieoorzaakgroepen:
                group_data = df_filtered[df_filtered['EMISSIEOORZAAKLABEL'] == group]
                group_orig_data = df_orig_filtered[df_orig_filtered['EMISSIEOORZAAKLABEL'] == group]
            
                # Scatter plot for original data (input, before regression)
                ax[ivar].scatter(
                    x=group_orig_data['EMISSIE'], y=group_orig_data[plotvar], c='royalblue', s=10, alpha=0.6, label=f"Group {group} (Original)"
            )
                # Scatter plot for gap-filled data
                ax[ivar].scatter(
                    x=group_data['EMISSIE'], y=group_data[plotvar], c='crimson', s=20, alpha=0.6, label=f"Group {group} (Gap-filled)"
                )

                # Perform regression for each group
                df_valid = group_data[group_data[plotvar] > 0]
            
                if df_valid.empty:
                    continue  # Skip this group if there is no valid data
            
                x = np.log10(df_valid['EMISSIE'])
                y = np.log10(df_valid[plotvar]) if plotvar in log_vars else df_valid[plotvar]

                poly = PolynomialFeatures(orders[plotvar])
                x_poly = poly.fit_transform(x.values.reshape(-1, 1))
                model = sm.OLS(y, x_poly).fit()

                x_range = np.linspace(min(x), max(x), 100)
                x_range_poly = poly.transform(x_range.reshape(-1, 1))
                y_pred = model.predict(x_range_poly)
                y_pred = 10**y_pred if plotvar in log_vars else y_pred

                # Regression line for the current group
                ax[ivar].plot(10**x_range, y_pred, 'k--', label=f'Group {group} Fit')

            ax[ivar].set_title(plotvar)
            ax[ivar].set_xscale('log')
            ax[ivar].set_xlim([1e3, 1e10])

            if plotvar in log_vars:
                ax[ivar].set_yscale('log')
                ax[ivar].set_ylim([1e-2, 1e4])
            elif plotvar == 'TEMPERATUUR':
                ax[ivar].set_ylim(bottom=10)

        fig.suptitle('Emission Cause Group Analysis')
        plt.tight_layout()
        plt.show()

        
        #Prepare final model input [per hour] for each hour of simulaiton period:

    def prepare_final_input(self, df):
        """
        Prepares unified point source NetCDF inputs for all tracers and writes them once per timestep.
        """

        print("Preparing unified point sources input NetCDF files...")

        time_array = self.prepare_time_array()  # array of datetime objects
        
        # Load temporal profiles (hour, week, month)
        tprof_hour, tprof_week, tprof_mnth = temp_disaggregation_new_cams_temp_prof.loadsnap_cams(self.spec_name, self.year_sim, noxdir, temp_prof_file, x0,  y0, xres, yres, itot, jtot )
    
        for timepoint in range(self.hour_start, time_array.size):
            timestamp = time_array[timepoint]
            #year, month, day, hour = timestamp.year, timestamp.month, timestamp.day, timestamp.hour

            #ihour = hour
            #imonth = month - 1  # 0-indexed
            #iweek = datetime.datetime(self.year_sim, month, day, hour).weekday()


            # Original timestamp is UTC
            utc_t = timestamp

            # Convert UTC -> local time ONLY for temporal profiles
            local_t = get_local_time(utc_t)

            year  = local_t.year
            month = local_t.month
            day   = local_t.day
            hour  = local_t.hour

            ihour  = hour
            imonth = month - 1
            iweek  = local_t.weekday()

            blockx = self.xsize / self.nprocx
            blocky = self.ysize / self.nprocy
            dx = self.xsize / self.itot
            dy = self.ysize / self.jtot

            all_data = []
            total = 0.0

            # Collect all emissions from all subdomains
            for mpiidx in range(self.nprocx):
                xmin = (self.lcc_start_x+ self.les_grid.xt[0]) + mpiidx * blockx
                xmax = xmin + blockx

                for mpiidy in range(self.nprocy):
                    ymin = (self.lcc_start_y+ self.les_grid.yt[0]) + mpiidy * blocky
                    ymax = ymin + blocky


                    xt_list=self.les_grid.xt+self.lcc_start_x
                    yt_list=self.les_grid.yt+self.lcc_start_y

                    data = df2list_new(df, xt_list,yt_list, tracer_name=self.spec_name, domain_bounds=(xmin, xmax, ymin, ymax))

                    if data.size == 0:
                        continue
                    total += data[:, 4].sum()
                    
                    # Convert annual emissions to hourly
                    data[:, 4] = data[:, 4] / (365 * 24)
                    

                    for source in data:
                        
                            isnap = int(source[-1] - 1)
                            source[4] = (
                                        source[4]
                                        * tprof_hour[isnap, ihour]
                                        * tprof_week[isnap, iweek]
                                        * tprof_mnth[isnap, imonth]
                                )

                    all_data.append(data)

            if not all_data:
                print(f"No point sources for {utc_t.year}-{utc_t.month:02d}-{utc_t.day:02d} {utc_t.hour:02d}:00")
                continue

            # Concatenate all subdomain data into one array
            all_data = np.vstack(all_data)

            # Write using ORIGINAL UTC timestamp
            data2netc_point(
                all_data,
                author,
                email,
                self.targetdir,
                tracer=self.spec_name,
                minute=utc_t.minute,
                hour=utc_t.hour,
                day=utc_t.day,
                month=utc_t.month,
                year=utc_t.year
            )

            print('total',total)
            print(f"✔ Point sources written for: {utc_t.year}-{utc_t.month:02d}-{utc_t.day:02d} {utc_t.hour:02d}:00")



if __name__ == "__main__":
    PSEP = Pointsource_input_preparation()
    df, df_org = PSEP.identify_emission_categories()
    unassigned_points, df_gapfilled_full = PSEP.prepare_emission_data(df)
    #PSEP.plot_emission_data(df_selection_non_rem_full, df_gapfilled_full)
    PSEP.prepare_final_input(df_gapfilled_full)
