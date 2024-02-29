"""
Module Description:
This script serves the purpose of translating reprojected
CSV files with area emissions into NetCDF-files per SNAP category.

Creator Information:
- Creator: Dr. Arseni Doyennel
- Contact Email: a.doyennel@vu.nl
"""

from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import netCDF4 as netc
import os
import xarray as xr
import glob
import pandas as pd
import sys

# from qgis.core import QgsVectorLayer, QgsProject, NULL
# sys.path.insert(1,os.path.join(os.getcwd(), ''))
from emission_preparation_setting import (spec_name, year, x0, y0, xn, yn,
                                          xres, yres, nprocx, nprocy, dx, dy,  # setting script
                                          outputdir_main_ruisdael_area_gpkg2csv_new_python_only,
                                          outputdir_main_ruisdael_area_csv2netc , snaplist)

import numpy as np
import pandas as pd
import os
import glob
import netCDF4 as netc


class CSV_to_NCDF_emiss_convertor:
    def __init__(self):
        self.rundir_main = outputdir_main_ruisdael_area_gpkg2csv_new_python_only
        self.outputdir = outputdir_main_ruisdael_area_csv2netc
        self.snaplist = snaplist
        self.spec_name = spec_name
        self.x0 = x0
        self.y0 = y0
        self.xn = xn
        self.yn = yn
        self.xres = xres
        self.yres = yres
        self.dx = int(dx)
        self.dy = int(dy)

    def csv2netc_1to1_hr(self, snap_number, rundir, netc_target, xmin, xmax, dx, ymin, ymax, dy,
                         fillvalue=0, pos_only=False, getvals=False, writefile=False, verbose=False):

        if snap_number in (2, 7):
            tracer = self.spec_name + '_sum'
        else:
            tracer = self.spec_name


        if (os.path.isfile(netc_target) and not getvals):
            print(f'Target file {netc_target} already exists.')
        else:
            print(f'Creating: {netc_target}')

            xmin, xmax, dx = int(xmin), int(xmax), dx
            ymin, ymax, dy = int(ymin), int(ymax), dy

            x = np.arange(xmin, xmax, dx)
            y = np.arange(ymin, ymax, dy)

            print(x, y)

            data = np.ones([int((xmax - xmin) // dx), int((ymax - ymin) // dy)]) * fillvalue

            # Determine which csv are needed

            filelist = []
            regex = rundir + 'HARM_snap_' + str(snap_number) + '_*_aggr.csv'
            for name in glob.glob(regex):
                filelist.append(name)

            print(regex)
            print(f'Files found: {len(filelist)}')

            # Loop over csv filelist
            for ifile, filename in enumerate(filelist):
                print(f'{ifile + 1}/{len(filelist)} {filename}')
                df = pd.read_csv(filename, delimiter=',')
                df = df.loc[df[str(tracer)] != 'NULL']

                for idx, content in df.iterrows():
                    xloc = content['left']
                    yloc = content['bottom']
                    co2 = content[str(tracer)]
                    if pos_only:
                        co2 = max(0, co2)

                    xidx = int((xloc - xmin) // dx)
                    yidx = int((yloc - ymin) // dy)

                    if (xidx >= 0) and (yidx >= 0) and (xidx < len(x)) and (yidx < len(y)):
                        data[xidx, yidx] = co2
                    else:
                        if verbose: print('Outside domain bounds:', 'xloc', xidx, 'yloc', yloc, 'data', data)

            if writefile and not os.path.isfile(netc_target):
                with netc.Dataset(netc_target, 'w') as fobj:
                    dim_x = fobj.createDimension('x', len(x))
                    dim_y = fobj.createDimension('y', len(y))

                    var_x = fobj.createVariable('x', 'f4', ('x',))
                    var_y = fobj.createVariable('y', 'f4', (    'y',))

                    var_x[:] = np.array(x)
                    var_y[:] = np.array(y)

                    print(xr.DataArray(data[:, :]))

                    var_e = fobj.createVariable(self.spec_name, 'f8', ('x', 'y',))
                    var_e[:, :] = data[:, :]

                    var_x.units = 'HARMONIE x in meters'
                    var_y.units = 'HARMONIE y in meters'
                    var_e.units = 'Kilogram per year'

            if getvals: return np.array(x), np.array(y), data

    def prep_ncdfs(self):
        if not os.path.exists(self.outputdir):
            os.makedirs(self.outputdir)
            print(f'The {self.outputdir} directory is created!')
        else:
            print(f'The {self.outputdir} directory exists!')

        if not self.xn % self.dx and not self.yn % self.dy:
            print('dx/dy is a multiple of dx_scale/dy_scale!')
        else:
            raise Exception("check dx and dy before proceed!")



        xmax = self.x0 + int(self.xn * self.xres)
        ymax = self.y0 + int(self.yn * self.yres)

        xmin, xmax, dx = self.x0, xmax, self.xres
        ymin, ymax, dy = self.y0, ymax, self.yres


        for snap in self.snaplist:

            rundir = self.rundir_main + 'snap_' + str(snap) + '_csv/'

            for f in glob.glob(self.outputdir + '*_' + str(snap) + '_all.nc'):
                os.remove(f)

            self.csv2netc_1to1_hr(snap, rundir, self.outputdir + f'HARM_snap_{snap}_all.nc', xmin, xmax, dx,
                                  ymin, ymax, dy, getvals=False, writefile=True)


if __name__ == "__main__":
    # Instantiate the class
    CSV2NCDF = CSV_to_NCDF_emiss_convertor()
    CSV2NCDF.prep_ncdfs()