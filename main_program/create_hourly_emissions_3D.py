"""
Module Description:

This script is designed to create a final area emission input for DALES at 
the selected simulation domain.
Temporal aggregation is applied to account for
changes in the emission strength of a source
 for a certain hour in a day, day in the week, or month of the day.

Final input files are rather “DALES-unlike” as they are date specific and
cover the full simulation domain.
They are independent of processors settings or the starting time of the simulation.
The files are labeled as e.g. {tracer}_emis_{yyyy}{mm}{dd}{hh}{mm}_3d.nc

Author: Dr. Arseni Doyennel
Contact Email: a.doyennel@vu.nl
"""

import numpy as np
import matplotlib.pyplot as plt
from emission_construction_functions import reademisoptions, writereademission_3d
import xarray as xr
import glob
import shutil
import os
import pandas as pd

import sys
from emission_preparation_setting import (targetdir_ruisdael_area_total_static,
                                          output_dir_create_hourly_emissions_3D, spec_name, xres, yres,
                                          xn, yn, snaplist, x0, y0, zmax, year_start,
                                          month_start, day_start, hour_start, year_end,
                                          month_end, day_end, hour_end
)


#self.setting = EmissionPreparationSettings()
class Final_3D_input:
    def __init__(self):
        self.targetdir = targetdir_ruisdael_area_total_static
        self.output_dir = output_dir_create_hourly_emissions_3D
        self.spec_name = spec_name
        self.xres = xres
        self.yres = yres
        self.xn = xn
        self.yn = yn
        self.snaplist = snaplist
        self.x0 = x0
        self.y0 = y0
        self.zmax = zmax
        self.year_start = year_start
        self.month_start = month_start
        self.day_start = day_start
        self.hour_start = hour_start
        self.year_end = year_end
        self.month_end = month_end
        self.day_end = day_end
        self.hour_end = hour_end

         # Create the dictionary mapping SNAP numbers to categories:
        self.snap_category_map = {
            1: '1_power',
            2: '2_residential_commercial',
            3: '3_industrial_combustion',
            4: '4_industrial_processes',
            5: '5_fossil_fuels',
            7: '7_road',
            8: '8_mobile',
            9: '9_waste',
            10: '10_agriculture'
        }

    def create_output_directory(self):
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)
            print(f'The {self.output_dir} directory is created!')
        else:
            print(f'The {self.output_dir} directory exists!')
            for f in glob.glob(self.output_dir + '*.nc'):
                os.remove(f)

    def create_settings_file(self):
        xmin = self.x0
        ymin = self.y0
        startyear = self.year_start
        startmonth = self.month_start
        startday = self.day_start
        starthour = self.hour_start

        time_array = pd.Series(
            pd.date_range(start=f'{startyear}-{startmonth}-{startday} {starthour}:00:00',
                          end=f'{self.year_end}-{self.month_end}-{self.day_end} '
                              f'{self.hour_end}:00:00', freq='H'))
        runlength = time_array.size - 1  # should be minus one time point to work correctly

        xmax = xmin + int(self.xn * self.xres)
        ymax = ymin + int(self.yn * self.yres)

        xmax = xmax - self.xres  # should be minus one grid point to work correctly
        ymax = ymax - self.yres  # should be minus one grid point to work correctly

        zmin = 0
        tracer = self.spec_name

        # Create the string tras_and_categories dynamically
        tras_and_categories_list = ['tracer ' + str(self.spec_name)]
        for snap_number in self.snaplist:
            category = self.snap_category_map.get(snap_number)
            if category:
                tras_and_categories_list.append(f'{category} {snap_number}')

        tras_and_categories = ' '.join(tras_and_categories_list)

        lines = ['', f'startyear {str(startyear)}', f'startmonth {str(startmonth)}', f'startday {str(startday)}',
                 f'starthour {str(starthour)}', f'runlength {str(runlength)}']
        with open('settings_hourly.txt', 'w') as f:
            for line in lines:
                f.write(line)
                f.write('\n')

        more_lines = ['', f'xmin {str(xmin)}', f'xmax {str(xmax)}', f'ymin {str(ymin)}', f'ymax {str(ymax)}',
                      f'zmin {str(zmin)}', f'zmax {str(self.zmax)}']

        with open('settings_hourly.txt', 'a') as ff:
            ff.write('\n'.join(more_lines))

        more_lines_2 = ['', '', f'{tras_and_categories}']

        with open('settings_hourly.txt', 'a') as fff:
            fff.write('\n'.join(more_lines_2))

    def create_emissions(self):
        domain_bounds, tstart, tend, tracers, sources, categories = reademisoptions("settings_hourly.txt", show_log=True)
        print(domain_bounds)
        writereademission_3d(self.targetdir, self.output_dir, domain_bounds, tstart, tend, tracers,
                             sources, categories, show_log=True)


if __name__ == "__main__":
    emission_final_prep = Final_3D_input()
    emission_final_prep.create_output_directory()
    emission_final_prep.create_settings_file()
    emission_final_prep.create_emissions()
