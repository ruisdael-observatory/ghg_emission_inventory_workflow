import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
from pyproj import Transformer
from GridDales import GridDales
from scipy.interpolate import LinearNDInterpolator
from shapely.geometry import Polygon
from emission_preparation_setting import spec_name, proj_params, grid_params, x0, y0, x_offset, y_offset, rundir_ruisdael_area_residuals, outputdir_main_ruisdael_area_csv2netc, year_start, snaplist
from projection import Transform
from netCDF4 import Dataset
from shapely.geometry import box
import rasterio.features
from tqdm import tqdm  # optional, progress bar
from scipy.ndimage import gaussian_filter

class EmissionResampler:
    def __init__(self, snap, spec_name, proj_params, grid_params, x0 , y0, x_offset, y_offset, rundir_ruisdael_area_residuals, outputdir_main_ruisdael_area_csv2netc, year_start):
        self.snap = snap
        self.spec_name = spec_name
        self.proj_params = proj_params
        self.grid_params = grid_params
        self.input_folder = rundir_ruisdael_area_residuals
        self.output_folder = outputdir_main_ruisdael_area_csv2netc
        self.year = year_start
        
        # Initialize transformer and grid
        self.transformer = Transform(self.proj_params)
        self.les_grid = GridDales(self.grid_params)
        self.x_offset = x_offset
        self.y_offset = y_offset
        self.lcc_start_x, self.lcc_start_y = x0+self.x_offset, y0+self.y_offset
        # Update grid southwest corner
        self.grid_params['southwest_x'] = self.lcc_start_x
        self.grid_params['southwest_y'] = self.lcc_start_y

        # placeholder for data that should be assigned during processing:
        self.resampled_data = None  # 2D np.array for resampled data
        self.lat = None             # 2D np.array of latitudes
        self.lon = None             # 2D np.array of longitudes

    def create_emission_polygons(df):
        cell_size = 1000 if self.snap != 2 else 100  # in meters
        half = cell_size / 2
        polygons = []
        for i in range(len(df)):
            x = df['lcc_x'].iloc[i]
            y = df['lcc_y'].iloc[i]
            val = df[self.spec_name].iloc[i]
            poly = box(x - half, y - half, x + half, y + half)
            polygons.append((poly, val))
        return polygons

    def smooth_redistributed_data(self):

        cell_size = 1000 if self.snap != 2 else 100  # in meters
        sigma = max((cell_size / (self.les_grid.xt[2]-self.les_grid.xt[1])) / 2, 0.5)
        # Step 1: Compute original total mass
        original_mass = np.nansum(self.resampled_data) * self.cell_area_km2_resampled

        # Step 2: Apply smoothing filter
        smoothed = gaussian_filter(self.resampled_data, sigma=sigma)

        # Step 3: Compute mass after smoothing
        smoothed_mass = np.nansum(smoothed) * self.cell_area_km2_resampled

        # Step 4: Scale smoothed data to preserve mass
        if smoothed_mass > 0:
            correction_factor = original_mass / smoothed_mass
            smoothed *= correction_factor
        else:
            correction_factor = 1.0  # fallback if mass is zero

        # Step 5: Store smoothed and corrected data
        self.resampled_data = smoothed
        print(f"Applied smoothing (σ={sigma}) with correction factor {correction_factor:.4f} to preserve mass.")


    def load_data(self):
        if self.snap == 2:
            filename = os.path.join(self.input_folder, f"{self.spec_name}_{self.year}_SNAP_2_refined_100m_residual.csv")
            print(f"[INFO] SNAP 2 detected: loading 100×100m refined data from {filename}") 
        else:
            filename = os.path.join(self.input_folder, f"{self.spec_name}_{self.year}_SNAP_{self.snap}_residual.csv")
            print(f"Loading data from {filename} ...")

        df = pd.read_csv(filename)
        df['lcc_x'], df['lcc_y'] = df['x'].values, df['y'].values
        self.df = df

    def filter_data(self):
        les_grid = self.les_grid
        lcc_start_x, lcc_start_y = self.lcc_start_x, self.lcc_start_y

        target_x_min = lcc_start_x+les_grid.dx/2
        target_y_min = lcc_start_y+les_grid.dx/2
        target_x_max = target_x_min + les_grid.dx * self.les_grid.itot
        target_y_max = target_y_min + les_grid.dy * self.les_grid.jtot

        buffer_x = 0  # uneeded,as increases the error
        buffer_y = 0  # uneeded,as increases the error

        expanded_target_x_min = target_x_min - buffer_x
        expanded_target_x_max = target_x_max + buffer_x
        expanded_target_y_min = target_y_min - buffer_y
        expanded_target_y_max = target_y_max + buffer_y

        mask = (
            (self.df['lcc_x'] >= expanded_target_x_min) & (self.df['lcc_x'] <= expanded_target_x_max) &
            (self.df['lcc_y'] >= expanded_target_y_min) & (self.df['lcc_y'] <= expanded_target_y_max)
        )
        self.filtered_df = self.df[mask]
        self.expanded_target_bounds = (expanded_target_x_min, expanded_target_x_max, expanded_target_y_min, expanded_target_y_max)

    def create_original_grid(self):
        fdf = self.filtered_df

        x_min, x_max = fdf['lcc_x'].min(), fdf['lcc_x'].max()
        y_min, y_max = fdf['lcc_y'].min(), fdf['lcc_y'].max()

        resolution = 100 if self.snap == 2 else 1000

        nx = int((x_max - x_min) / resolution) + 1
        ny = int((y_max - y_min) / resolution) + 1

        self.original_grid_x, self.original_grid_y = np.meshgrid(
            np.linspace(x_min, x_max, nx),
            np.linspace(y_min, y_max, ny)
        )
        self.original_co2_grid = np.zeros_like(self.original_grid_x)

        for i in range(len(fdf)):
            original_x = fdf['lcc_x'].iloc[i]
            original_y = fdf['lcc_y'].iloc[i]
            original_value = fdf[self.spec_name].iloc[i]

            target_x_idx = np.argmin(np.abs(self.original_grid_x[0, :] - original_x))
            target_y_idx = np.argmin(np.abs(self.original_grid_y[:, 0] - original_y))

            self.original_co2_grid[target_y_idx, target_x_idx] += original_value

    def plot_original(self):
        plt.figure(figsize=(10, 5))
        plt.pcolormesh(self.original_grid_x, self.original_grid_y, self.original_co2_grid, cmap="viridis", shading='auto')
        plt.colorbar(label=f"Original {self.spec_name} concentration")
        plt.title(f"Original 1x1km {self.spec_name} Data in Expanded Domain")
        plt.xlabel("LCC X")
        plt.ylabel("LCC Y")
        plt.show()

    def resample(self):
        les_grid = self.les_grid
        dx, dy = les_grid.dx, les_grid.dy
        nx, ny = len(les_grid.xt), len(les_grid.yt)

        # Coordinates of LES grid centers
        self.les_x = les_grid.xt + self.lcc_start_x
        self.les_y = les_grid.yt + self.lcc_start_y
        les_xc = self.les_x
        les_yc = self.les_y

        # Define LES grid cell polygons
        les_polygons = []
        for j in range(ny):
            row = []
            for i in range(nx):
                x0 = les_xc[i] - dx / 2
                x1 = les_xc[i] + dx / 2
                y0 = les_yc[j] - dy / 2
                y1 = les_yc[j] + dy / 2
                row.append(box(x0, y0, x1, y1))
            les_polygons.append(row)

        # Init mass grid
        mass_grid = np.zeros((ny, nx), dtype='float64')

        # Loop over each original emission cell
        for _, row in self.filtered_df.iterrows():
            center_x, center_y = row['lcc_x'], row['lcc_y']
            half = 50 if self.snap == 2 else 500
            poly = box(center_x - half, center_y - half, center_x + half, center_y + half)
            mass = row[self.spec_name]  # already per cell

            # Determine bounding box in LES index space
            col_min = int((center_x - half - self.lcc_start_x) // dx)
            col_max = int((center_x + half - self.lcc_start_x) // dx)
            row_min = int((center_y - half - self.lcc_start_y) // dy)
            row_max = int((center_y + half - self.lcc_start_y) // dy)

            for j in range(max(0, row_min), min(ny, row_max + 1)):
                for i in range(max(0, col_min), min(nx, col_max + 1)):
                    cell_poly = les_polygons[j][i]
                    inter_area = cell_poly.intersection(poly).area
                    if inter_area > 0:
                        frac = inter_area / poly.area
                        mass_grid[j, i] += mass * frac

        self.resampled_data = mass_grid
        self.cell_area_km2_resampled = 1.0  # Units already in kg/cell/year, no need to multiply with area


    def plot_resampled(self):
        plt.figure(figsize=(10, 5))
        plt.pcolormesh(self.les_x, self.les_y, self.resampled_data * self.cell_area_km2_resampled,
                       cmap="viridis", shading='auto')
        plt.colorbar(label=f"Resampled {self.spec_name} concentration")
        plt.title(f"Resampled 100m {self.spec_name} Data in Target Domain")
        plt.xlabel("LCC X")
        plt.ylabel("LCC Y")
        plt.show()

    def calculate_mass_error(self):
        expanded_target_x_min, expanded_target_x_max, expanded_target_y_min, expanded_target_y_max = self.expanded_target_bounds
        fdf = self.filtered_df

        # Build expanded target polygon
        target_domain = box(expanded_target_x_min, expanded_target_y_min, expanded_target_x_max, expanded_target_y_max)

        adjusted_original_mass = 0
        for _, row in fdf.iterrows():
            center_x, center_y = row['lcc_x'], row['lcc_y']
            half = 500 if self.snap != 2 else 50
            poly = box(center_x - half, center_y - half, center_x + half, center_y + half)
            overlap_area = poly.intersection(target_domain).area
            if overlap_area > 0:
                frac = overlap_area / poly.area
                adjusted_original_mass += row[self.spec_name] * frac

        resampled_mass = np.nansum(self.resampled_data) * self.cell_area_km2_resampled
        error_pct = 100 * (resampled_mass - adjusted_original_mass) / adjusted_original_mass

        print(f"SNAP {self.snap}: Adjusted original mass (target domain): {adjusted_original_mass:.2f}, Resampled mass: {resampled_mass:.2f}")
        print(f"SNAP {self.snap}: Mass conservation error: {error_pct:.2f}%")

    def run(self):
        self.load_data()
        self.filter_data()
        self.create_original_grid()
        self.plot_original()
        self.resample()
        self.smooth_redistributed_data()
        self.plot_resampled()
        self.calculate_mass_error()

    def save_to_netcdf(self, filename=None):
        
        # Convert LES x, y coordinates to lat/lon using your existing transformer
        x_les_2d,y_les2d=np.meshgrid(self.les_x, self.les_y)
        # Convert LES x, y coordinates to lat/lon using your existing transformer
        self.lat, self.lon = self.transformer.xy_to_latlon(x_les_2d,y_les2d)
        
        if filename is None:
            filename = os.path.join(self.output_folder, f'HARM_snap_{self.snap}_all_{self.spec_name}.nc')
    
        output_dir = os.path.dirname(filename)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        ncfile = Dataset(filename, 'w', format='NETCDF4')

        # Define dimensions
        ncfile.createDimension('x', len(self.les_grid.xt))
        ncfile.createDimension('y', len(self.les_grid.yt))

        # Create variables
        x_var = ncfile.createVariable('x', 'f4', ('x',))
        y_var = ncfile.createVariable('y', 'f4', ('y',))
        lat_var = ncfile.createVariable('lat', 'f4', ('y', 'x'))
        lon_var = ncfile.createVariable('lon', 'f4', ('y', 'x'))
        spec_var = ncfile.createVariable(self.spec_name, 'f4', ('y', 'x')) 

        # Assign values
        x_var[:] = self.les_grid.xt + self.lcc_start_x
        y_var[:] = self.les_grid.yt + self.lcc_start_y
        lat_var[:, :] = self.lat
        lon_var[:, :] = self.lon
        spec_var[:, :] = self.resampled_data

        # Add global attribute
        ncfile.title = f'{self.spec_name} Emissions for {self.year} Data'

        ncfile.close()
        print(f'NetCDF file created successfully: {filename}')

if __name__ == "__main__":
    for snap in snaplist:
        if snap == 7:
            continue  # Skip SNAP 7
        print(f"Processing SNAP {snap}")
        resampler = EmissionResampler(snap, spec_name, proj_params, grid_params, x0 , y0, x_offset, y_offset, rundir_ruisdael_area_residuals, outputdir_main_ruisdael_area_csv2netc, year_start)
        resampler.run()
        resampler.save_to_netcdf()