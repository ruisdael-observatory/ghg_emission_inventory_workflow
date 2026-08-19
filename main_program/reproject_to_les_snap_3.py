import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
from shapely.geometry import box
from shapely.ops import unary_union
import rasterio
from rasterio.features import rasterize
from netCDF4 import Dataset
from scipy.ndimage import gaussian_filter
from tqdm import tqdm
from projection import Transform
from GridDales import GridDales
from emission_preparation_setting import spec_name, proj_params, grid_params, x0, y0, x_offset, y_offset, rundir_ruisdael_area_residuals, cbs_loc_consu_preprocessing, outputdir_main_ruisdael_area_csv2netc, year_start, snaplist, industrial_file_name
import geopandas as gpd

class SNAP3Refiner:
    def __init__(self, snap, industrial_file_link, industrial_file_name, input_link, proj_params, grid_params, 
                 x0, y0, x_offset, y_offset, spec_name, year, output_folder):
        """
        industrial_file_name: GeoDataFrame of industrial-related buildings from BAG data base
        https://service.pdok.nl/lv/bag/atom/bag.xml
        
        These includes: 
        'industriefunctie',          # Core industry
        'bedrijf',                   # General business (catch-all for small industries)
        'opslagfunctie',             # Industrial storage/logistics
        'energieopwekking',          # Energy/production (if present)
        'nutsvoorziening'            # Utility buildings (could be indirect)
        
        snap_csv: path to SNAP3 CSV with columns x, y, spec_name
        """
        self.snap=snap
        self.industrial_file_link = industrial_file_link
        self.industrial_file_name=industrial_file_name
        self.gdf_target = gpd.read_file(f"{self.industrial_file_link}{self.industrial_file_name}")
        self.input_link=input_link
        self.spec_name = spec_name
        self.year = year
        self.input_file=f'{self.spec_name}_{self.year}_SNAP_3_residual.csv'
        self.snap_csv = f"{self.input_link}{self.input_file}"
        self.proj_params = proj_params
        self.grid_params = grid_params
        self.transformer = Transform(self.proj_params)
        self.les_grid = GridDales(self.grid_params)
        self.x_offset = x_offset
        self.y_offset = y_offset
        self.lcc_start_x, self.lcc_start_y = x0 + self.x_offset, y0 + self.y_offset
        self.grid_params['southwest_x'] = self.lcc_start_x
        self.grid_params['southwest_y'] = self.lcc_start_y
        self.output_folder = output_folder
        self.les_x = self.les_grid.xt + self.lcc_start_x
        self.les_y = self.les_grid.yt + self.lcc_start_y

    def load_snap10(self):
        df = pd.read_csv(self.snap_csv)
        df['lcc_x'], df['lcc_y'] = df['x'].values, df['y'].values
        self.df_snap = df

    def filter_snap10(self):
        target_x_min = self.les_x[0] - self.les_grid.dx / 2
        target_x_max = self.les_x[-1] + self.les_grid.dx / 2
        target_y_min = self.les_y[0] - self.les_grid.dy / 2
        target_y_max = self.les_y[-1] + self.les_grid.dy / 2

        mask = (
            (self.df_snap['lcc_x'] >= target_x_min) & (self.df_snap['lcc_x'] <= target_x_max) &
            (self.df_snap['lcc_y'] >= target_y_min) & (self.df_snap['lcc_y'] <= target_y_max)
        )
        self.df_snap_filtered = self.df_snap[mask].copy()

    def rasterize_brp_aligned(self):
        """
        Rasterize gdf_target aligned to LES domain and refinement resolution.
        This ensures Netherlands territory appears correctly.
        """
        nx_les = len(self.les_x)
        ny_les = len(self.les_y)
        dx = self.les_grid.dx
        dy = self.les_grid.dy

        # Compute raster extent aligned to LES domain
        xmin = self.les_x[0] - dx/2
        xmax = self.les_x[-1] + dx/2
        ymin = self.les_y[0] - dy/2
        ymax = self.les_y[-1] + dy/2

        width = int(np.ceil((xmax - xmin)/dx))
        height = int(np.ceil((ymax - ymin)/dy))

        transform = rasterio.transform.from_origin(xmin, ymax, dx, dy)

        # Rasterize only the polygons inside LES domain
        shapes = ((geom, 1) for geom in self.gdf_target.geometry)
        brp_raster = rasterize(shapes, out_shape=(height, width),
                           transform=transform, fill=0, all_touched=True, dtype='float32')

        self.brp_raster = brp_raster
        self.brp_transform = transform
        self.raster_xmin = xmin
        self.raster_ymax = ymax
        self.raster_res = dx

    def compute_agri_fraction_aligned(self):
        """
        Compute fraction of industrial buildings area per LES cell using aligned raster.
        """
        nx, ny = len(self.les_x), len(self.les_y)
        dx, dy = self.les_grid.dx, self.les_grid.dy

        fraction_grid = np.zeros((ny, nx), dtype='float32')

        for j in range(ny):
            for i in range(nx):
                # LES cell boundaries
                x0 = self.les_x[i] - dx/2
                x1 = self.les_x[i] + dx/2
                y0 = self.les_y[j] - dy/2
                y1 = self.les_y[j] + dy/2

                # Convert to raster indices
                col_start = int((x0 - self.raster_xmin)/self.raster_res)
                col_end = int((x1 - self.raster_xmin)/self.raster_res)
                row_start = int((self.raster_ymax - y1)/self.raster_res)
                row_end = int((self.raster_ymax - y0)/self.raster_res)

                # Clip
                col_start = max(0, col_start); col_end = min(self.brp_raster.shape[1], col_end)
                row_start = max(0, row_start); row_end = min(self.brp_raster.shape[0], row_end)

                cell_pixels = self.brp_raster[row_start:row_end, col_start:col_end]
                if cell_pixels.size > 0:
                    fraction_grid[j, i] = cell_pixels.mean()

        self.agri_fraction = fraction_grid


    def refine_snap10(self):
        """Refine SNAP3 to LES grid using industrial buildings fraction and area-weighted redistribution.
        Preserves mass by normalizing weights per SNAP cell, sets non-industrial LES cells to NaN for plotting,
        and performs domain-adjusted mass diagnostics.
        """
        nx, ny = len(self.les_x), len(self.les_y)
        dx, dy = self.les_grid.dx, self.les_grid.dy
        mass_grid = np.zeros((ny, nx), dtype='float64')

        snap_cell_size = 1000  # SNAP10 original 1 km²
        total_input_mass = 0.0
        total_unassigned_mass = 0.0
        adjusted_original_mass = 0.0

        # Build LES domain polygon
        xmin = self.les_x[0] - dx/2
        xmax = self.les_x[-1] + dx/2
        ymin = self.les_y[0] - dy/2
        ymax = self.les_y[-1] + dy/2
        les_domain = box(xmin, ymin, xmax, ymax)

        snap_mask = (
            (self.df_snap['lcc_x'] >= self.les_x[0] - self.les_grid.dx/2) &
            (self.df_snap['lcc_x'] <= self.les_x[-1] + self.les_grid.dx/2) &
            (self.df_snap['lcc_y'] >= self.les_y[0] - self.les_grid.dy/2) &
            (self.df_snap['lcc_y'] <= self.les_y[-1] + self.les_grid.dy/2)
        )
        self.df_snap_filtered = self.df_snap[snap_mask].copy()

        for idx, row in tqdm(self.df_snap_filtered.iterrows(),
                         total=len(self.df_snap_filtered), desc="Refining SNAP3 (debug)"):
            snap_x, snap_y, snap_val = row['lcc_x'], row['lcc_y'], row[self.spec_name]
            total_input_mass += snap_val

            snap_poly = box(snap_x - snap_cell_size/2, snap_y - snap_cell_size/2,
                        snap_x + snap_cell_size/2, snap_y + snap_cell_size/2)

            # Clip to LES domain
            snap_poly_clipped = snap_poly.intersection(les_domain)
            if snap_poly_clipped.is_empty:
                total_unassigned_mass += snap_val
                #print(f"[DEBUG] SNAP cell {idx} outside LES domain, mass unassigned: {snap_val}")
                continue
            adjusted_original_mass += snap_val * (snap_poly_clipped.area / snap_poly.area)

            # Determine overlapping LES cells
            col_min = max(0, int(np.floor((snap_poly_clipped.bounds[0] - self.lcc_start_x) / dx)))
            col_max = min(nx-1, int(np.floor((snap_poly_clipped.bounds[2] - self.lcc_start_x) / dx)))
            row_min = max(0, int(np.floor((snap_poly_clipped.bounds[1] - self.lcc_start_y) / dy)))
            row_max = min(ny-1, int(np.floor((snap_poly_clipped.bounds[3] - self.lcc_start_y) / dy)))

            weights = []
            idxs = []
            for j in range(row_min, row_max + 1):
                for i in range(col_min, col_max + 1):
                    les_poly = box(self.les_x[i] - dx/2, self.les_y[j] - dy/2,
                               self.les_x[i] + dx/2, self.les_y[j] + dy/2)
                    inter_area = snap_poly_clipped.intersection(les_poly).area
                    if inter_area <= 0:
                        continue
                    agri_frac = float(self.agri_fraction[j, i])
                    w = inter_area * agri_frac
                    if w > 0:
                        weights.append(w)
                        idxs.append((j, i))

            if len(weights) == 0:
                total_unassigned_mass += snap_val
                #print(f"[DEBUG] SNAP cell {idx} overlaps no agri LES cells, mass unassigned: {snap_val}")
                continue

            weights = np.array(weights, dtype='float64')
            W = weights.sum()
            if W <= 0:
                total_unassigned_mass += snap_val
                print(f"[DEBUG] SNAP cell {idx} has zero total weight, mass unassigned: {snap_val}")
                continue

            frac = weights / W
            for (j, i), f in zip(idxs, frac):
                mass_grid[j, i] += snap_val * f
        
        # For visualization
        mask_no_agri = (self.agri_fraction <= 0)
        self.resampled_data=mass_grid
        self.resampled_data_plot = self.resampled_data.copy()
        self.resampled_data_plot[mask_no_agri] = np.nan
        self.resampled_data[mask_no_agri] = 0.0

        total_assigned_mass = np.nansum(self.resampled_data)
        print(f"[INFO] SNAP3 refinement debug:")
        print(f"  Total SNAP input mass: {total_input_mass:,.3f}")
        print(f"  Adjusted original mass (in domain): {adjusted_original_mass:,.3f}")
        print(f"  Assigned mass: {total_assigned_mass:,.3f}")
        print(f"  Unassigned mass: {total_unassigned_mass:,.3f}")
        rel_err = 100.0 * (total_assigned_mass - adjusted_original_mass) / adjusted_original_mass
        print(f"  Relative mass error: {rel_err:.4f}%")
    
        return self.resampled_data

    def smooth_preserve_mass(self, sigma=None):
        if sigma is None:
            sigma = max(self.les_grid.dx / 2, 0.5)
        original_mass = np.nansum(self.resampled_data)
        smoothed = gaussian_filter(self.resampled_data, sigma=sigma)
        smoothed_mass = np.nansum(smoothed)
        if smoothed_mass > 0:
            smoothed *= original_mass / smoothed_mass
        self.resampled_data = smoothed

    def plot_resampled(self):
        plt.figure(figsize=(10, 5))
        plt.pcolormesh(self.les_x, self.les_y, self.resampled_data, cmap='viridis', shading='auto')
        plt.colorbar(label=f"{self.spec_name} (kg/cell/year)")
        plt.title(f"Refined SNAP3 {self.spec_name} emissions ({self.year})")
        plt.xlabel("X RD")
        plt.ylabel("Y RD")
        plt.show()

    def save_to_netcdf(self, filename=None):
        x_2d, y_2d = np.meshgrid(self.les_x, self.les_y)
        lat, lon = self.transformer.xy_to_latlon(x_2d, y_2d)

        if filename is None:
            filename = os.path.join(self.output_folder, f'HARM_snap_{self.snap}_all_{self.spec_name}.nc')
        os.makedirs(os.path.dirname(filename), exist_ok=True)

        ncfile = Dataset(filename, 'w', format='NETCDF4')
        ncfile.createDimension('x', len(self.les_x))
        ncfile.createDimension('y', len(self.les_y))

        x_var = ncfile.createVariable('x', 'f4', ('x',))
        y_var = ncfile.createVariable('y', 'f4', ('y',))
        lat_var = ncfile.createVariable('lat', 'f4', ('y', 'x'))
        lon_var = ncfile.createVariable('lon', 'f4', ('y', 'x'))
        spec_var = ncfile.createVariable(self.spec_name, 'f4', ('y', 'x'))

        x_var[:] = self.les_x
        y_var[:] = self.les_y
        lat_var[:, :] = lat
        lon_var[:, :] = lon
        spec_var[:, :] = self.resampled_data

        ncfile.title = f'{self.spec_name} Emissions for {self.year} Data'
        ncfile.close()
        print(f"NetCDF saved: {filename}")

    def run(self, smooth=False):
        print("[1] Loading SNAP3 data...")
        self.load_snap10()
        print("[2] Filtering SNAP3 to LES domain...")
        self.filter_snap10()
        print("[3] Rasterizing BAG industrial building data at refinement resolution...")
        self.rasterize_brp_aligned()
        print("[4] Computing industrial building area fraction...")
        self.compute_agri_fraction_aligned()
        print("[5] Refining SNAP3 emissions with mass-conserving redistribution...")
        self.refine_snap10()
        if smooth:
            print("[6] Smoothing while preserving mass...")
            self.smooth_preserve_mass()
        print("[7] Plotting refined SNAP3...")
        self.plot_resampled()
        print("[8] Saving NetCDF...")
        self.save_to_netcdf()
        print("Finished SNAP3 refinement.")

if __name__ == "__main__":
    for snap in snaplist:
        if snap == 3:  
            print(f"Processing SNAP {snap}")
            resampler = SNAP3Refiner(snap, cbs_loc_consu_preprocessing, industrial_file_name,rundir_ruisdael_area_residuals, 
                                      proj_params, grid_params, x0, y0,x_offset, y_offset, 
                                      spec_name, year_start, outputdir_main_ruisdael_area_csv2netc)
            resampler.run()