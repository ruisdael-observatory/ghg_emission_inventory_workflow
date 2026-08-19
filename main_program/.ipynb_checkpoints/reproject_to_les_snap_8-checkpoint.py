import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
from shapely.geometry import box
import geopandas as gpd
from tqdm import tqdm
from shapely.ops import unary_union
from netCDF4 import Dataset
from scipy.ndimage import gaussian_filter
from projection import Transform
from GridDales import GridDales
from emission_preparation_setting import (
    spec_name, proj_params, grid_params, x0, y0, x_offset, y_offset,
    rundir_ruisdael_area_residuals, outputdir_main_ruisdael_area_csv2netc,
    cbs_loc_consu_preprocessing, year_start, snaplist, vaarweg_file
)

class SNAP8Refiner:

    def __init__(self, snap, vaarweg_file_link, vaarweg_file, input_link, proj_params, grid_params,
                 x0, y0, x_offset, y_offset, spec_name, year, output_folder):  
                      
        """
        vaarweg_file: GeoDataFrame of inland water pathways
        https://downloads.rijkswaterstaatdata.nl/nwb-wegen/geogegevens/geopackage/Nederland_totaal/
        snap_csv: path to SNAP8 CSV with columns x, y, spec_name
        
        """
        
        self.snap = snap
        self.vaarweg_file_link = vaarweg_file_link
        self.vaarweg_file = f'{self.vaarweg_file_link}{vaarweg_file}'
        self.input_link = input_link
        self.spec_name = spec_name
        self.year = year
        self.input_file = f'{self.spec_name}_{self.year}_SNAP_8_residual.csv'
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

        # Load NWB vaarwegvakken
        self.gdf_vaar = gpd.read_file(self.vaarweg_file, layer="vaarwegvakken")
        assert self.gdf_vaar.crs is not None, "vaarwegvakken must have CRS defined"
        print(f"[INFO] Loaded vaarwegvakken with {len(self.gdf_vaar)} segments")

    def load_snap8(self):
        df = pd.read_csv(self.snap_csv)
        df['lcc_x'], df['lcc_y'] = df['x'].values, df['y'].values
        self.df_snap = df

    def filter_snap8(self):
        target_x_min = self.les_x[0] - self.les_grid.dx / 2
        target_x_max = self.les_x[-1] + self.les_grid.dx / 2
        target_y_min = self.les_y[0] - self.les_grid.dy / 2
        target_y_max = self.les_y[-1] + self.les_grid.dy / 2

        mask = (
            (self.df_snap['lcc_x'] >= target_x_min) & (self.df_snap['lcc_x'] <= target_x_max) &
            (self.df_snap['lcc_y'] >= target_y_min) & (self.df_snap['lcc_y'] <= target_y_max)
        )
        self.df_snap_filtered = self.df_snap[mask].copy()

    def compute_shipping_fraction(self):
        """Compute fraction of shipping lines per LES cell using NWB vaarwegvakken."""
        dx = self.les_grid.dx
        dy = self.les_grid.dy
        xs = np.asarray(self.les_x)
        ys = np.asarray(self.les_y)
        nx = len(xs)
        ny = len(ys)

        polygons = []
        cell_ids = []
        for j, y in enumerate(ys):
            for i, x in enumerate(xs):
                polygons.append(box(x - dx/2, y - dy/2, x + dx/2, y + dy/2))
                cell_ids.append((j, i))

        grid = gpd.GeoDataFrame({
            "cell_id": cell_ids,
            "row": [c[0] for c in cell_ids],
            "col": [c[1] for c in cell_ids],
            "geometry": polygons
        }, crs=self.gdf_vaar.crs)

        candidates = gpd.sjoin(grid, self.gdf_vaar[["vwk_id", "geometry"]],
                               how="inner", predicate="intersects")

        print(f"[INFO] Candidate cell-line pairs: {len(candidates)}")

        lengths = []
        for _, row in candidates.iterrows():
            cell_geom = row.geometry
            line_geom = self.gdf_vaar.loc[row['index_right'], 'geometry']
            inter = cell_geom.intersection(line_geom)
            l = inter.length if not inter.is_empty else 0.0
            lengths.append((row['row'], row['col'], l))

        fraction_mask = np.zeros((ny, nx), dtype=float)
        cell_diag = np.sqrt(dx**2 + dy**2)
        for r, c, l in lengths:
            fraction_mask[r, c] += l / cell_diag

        self.shipping_fraction = np.clip(fraction_mask, 0, 1)
        print(f"[INFO] Shipping fraction grid computed (mean={self.shipping_fraction.mean():.3f})")

    def refine_snap8(self):
        """Mass-conserving redistribution for SNAP8, using shipping fraction."""
        nx, ny = len(self.les_x), len(self.les_y)
        dx, dy = self.les_grid.dx, self.les_grid.dy
        mass_grid = np.zeros((ny, nx), dtype='float64')

        snap_cell_size = 1000
        total_input_mass = 0.0
        total_assigned_mass = 0.0

        # LES domain polygon
        xmin = self.les_x[0] - dx/2
        xmax = self.les_x[-1] + dx/2
        ymin = self.les_y[0] - dy/2
        ymax = self.les_y[-1] + dy/2
        les_domain = box(xmin, ymin, xmax, ymax)

        for idx, row in tqdm(self.df_snap_filtered.iterrows(),
                         total=len(self.df_snap_filtered), desc="Refining SNAP8"):
            snap_x, snap_y, snap_val = row['lcc_x'], row['lcc_y'], row[self.spec_name]
            total_input_mass += snap_val

            # Original SNAP cell polygon (1 km²)
            snap_poly = box(snap_x - snap_cell_size/2, snap_y - snap_cell_size/2,
                        snap_x + snap_cell_size/2, snap_y + snap_cell_size/2)
            snap_poly_clipped = snap_poly.intersection(les_domain)
            if snap_poly_clipped.is_empty:
                continue  # skip completely outside LES

            # Determine overlapping LES cells
            col_min = max(0, int(np.floor((snap_poly_clipped.bounds[0] - self.lcc_start_x) / dx)))
            col_max = min(nx-1, int(np.floor((snap_poly_clipped.bounds[2] - self.lcc_start_x) / dx)))
            row_min = max(0, int(np.floor((snap_poly_clipped.bounds[1] - self.lcc_start_y) / dy)))
            row_max = min(ny-1, int(np.floor((snap_poly_clipped.bounds[3] - self.lcc_start_y) / dy)))

            weights = []
            idxs = []

            # Collect positive shipping fractions for overlapping cells
            for j in range(row_min, row_max + 1):
                for i in range(col_min, col_max + 1):
                    ship_frac = float(self.shipping_fraction[j, i])
                    if ship_frac > 0:
                        weights.append(ship_frac)
                        idxs.append((j, i))

            if len(weights) > 0:
                # Normalize weights to sum to 1 (redistribute mass)
                weights = np.array(weights, dtype='float64')
                frac = weights / weights.sum()
                for (j, i), f in zip(idxs, frac):
                    mass_grid[j, i] += snap_val * f
                total_assigned_mass += snap_val
            else:
                # No positive shipping fraction → distribute evenly over overlapping LES cells
                n_cells = (row_max - row_min + 1) * (col_max - col_min + 1)
                if n_cells > 0:
                    for j in range(row_min, row_max + 1):
                        for i in range(col_min, col_max + 1):
                            mass_grid[j, i] += snap_val / n_cells
                    total_assigned_mass += snap_val

        # Mask non-shipping cells for plotting (optional)
        mask_no_ship = (self.shipping_fraction <= 0)
        self.resampled_data = mass_grid
        self.resampled_data_plot = mass_grid.copy()
        self.resampled_data_plot[mask_no_ship] = np.nan

        print(f"[INFO] SNAP8 refinement summary:")
        print(f"  Total input mass: {total_input_mass:,.3f}")
        print(f"  Assigned mass: {total_assigned_mass:,.3f}")
        rel_err = 100.0 * (total_assigned_mass - total_input_mass) / total_input_mass
        print(f"  Relative mass error: {rel_err:.3f}%")

    def refine_snap8_old(self):
        """Mass-conserving redistribution for SNAP8, using shipping fraction."""
        nx, ny = len(self.les_x), len(self.les_y)
        dx, dy = self.les_grid.dx, self.les_grid.dy
        mass_grid = np.zeros((ny, nx), dtype='float64')

        snap_cell_size = 1000  # SNAP8 nominal cell size
        total_input_mass = 0.0
        total_unassigned_mass = 0.0
        total_assignable_mass = 0.0

        # LES domain polygon
        xmin = self.les_x[0] - dx/2
        xmax = self.les_x[-1] + dx/2
        ymin = self.les_y[0] - dy/2
        ymax = self.les_y[-1] + dy/2
        les_domain = box(xmin, ymin, xmax, ymax)

        for idx, row in tqdm(self.df_snap_filtered.iterrows(),
                         total=len(self.df_snap_filtered), desc="Refining SNAP8"):
            snap_x, snap_y, snap_val = row['lcc_x'], row['lcc_y'], row[self.spec_name]
            total_input_mass += snap_val

            # SNAP cell as polygon
            snap_poly = box(snap_x - snap_cell_size/2, snap_y - snap_cell_size/2,
                        snap_x + snap_cell_size/2, snap_y + snap_cell_size/2)
            snap_poly_clipped = snap_poly.intersection(les_domain)
            if snap_poly_clipped.is_empty:
                total_unassigned_mass += snap_val
                continue

            # Determine overlapping LES cells
            col_min = max(0, int(np.floor((snap_poly_clipped.bounds[0] - self.lcc_start_x) / dx)))
            col_max = min(nx-1, int(np.floor((snap_poly_clipped.bounds[2] - self.lcc_start_x) / dx)))
            row_min = max(0, int(np.floor((snap_poly_clipped.bounds[1] - self.lcc_start_y) / dy)))
            row_max = min(ny-1, int(np.floor((snap_poly_clipped.bounds[3] - self.lcc_start_y) / dy)))

            # Collect cells with positive shipping fraction
            idxs = []
            weights = []
            for j in range(row_min, row_max + 1):
                for i in range(col_min, col_max + 1):
                    frac = float(self.shipping_fraction[j, i])
                    if frac > 0:
                        idxs.append((j, i))
                        weights.append(frac)

            if len(weights) == 0:
                # No shipping cells: cannot assign mass
                total_unassigned_mass += snap_val
                continue

            weights = np.array(weights, dtype='float64')
            weights /= weights.sum()  # normalize to sum=1

            # Assign mass to LES cells
            for (j, i), w in zip(idxs, weights):
                mass_grid[j, i] += snap_val * w

            total_assignable_mass += snap_val

        # Mask non-shipping cells for plotting
        mask_no_ship = (self.shipping_fraction <= 0)
        self.resampled_data = mass_grid
        self.resampled_data_plot = mass_grid.copy()
        self.resampled_data_plot[mask_no_ship] = np.nan
        self.resampled_data[mask_no_ship] = 0.0

        total_assigned_mass = np.nansum(self.resampled_data)
        rel_err = 100.0 * (total_assigned_mass - total_assignable_mass) / total_assignable_mass

        print(f"[INFO] SNAP8 refinement summary:")
        print(f"  Total input mass: {total_input_mass:,.3f}")
        print(f"  Assignable in-domain mass: {total_assignable_mass:,.3f}")
        print(f"  Assigned mass: {total_assigned_mass:,.3f}")
        print(f"  Unassigned mass (no shipping cell): {total_unassigned_mass:,.3f}")
        print(f"  Relative mass error: {rel_err:.3f}%")
    
    def run(self, smooth=False):
        print("[1] Loading SNAP8 data...")
        self.load_snap8()
        print("[2] Filtering SNAP8 to LES domain...")
        self.filter_snap8()
        print("[3] Computing shipping fraction...")
        self.compute_shipping_fraction()
        print("[4] Refining SNAP8 emissions with mass-conserving redistribution...")
        self.refine_snap8()
        if smooth:
            print("[5] Smoothing while preserving mass...")
            self.smooth_preserve_mass()
        print("[6] Plotting and saving results...")
        self.plot_resampled()
        self.save_to_netcdf()
        print("Finished SNAP8 refinement.")

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
        plt.pcolormesh(self.les_x, self.les_y, self.resampled_data_plot,
                       cmap='viridis', shading='auto')
        plt.colorbar(label=f"{self.spec_name} (kg/cell/year)")
        plt.title(f"Refined SNAP8 {self.spec_name} emissions ({self.year})")
        plt.xlabel("RD X")
        plt.ylabel("RD Y")
        plt.show()

    def save_to_netcdf(self, filename=None):
        x_2d, y_2d = np.meshgrid(self.les_x, self.les_y)
        lat, lon = self.transformer.xy_to_latlon(x_2d, y_2d)
        if filename is None:
            filename = os.path.join(
                self.output_folder, f'HARM_snap_{self.snap}_all_{self.spec_name}.nc'
            )
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
        ncfile.title = f'{self.spec_name} Emissions for {self.year} Data (SNAP8 Shipping)'
        ncfile.close()
        print(f"NetCDF saved: {filename}")


if __name__ == "__main__":
    for snap in snaplist:
        if snap == 8:
            print(f"Processing SNAP {snap} (shipping)")
            resampler = SNAP8Refiner(snap, cbs_loc_consu_preprocessing, vaarweg_file,rundir_ruisdael_area_residuals, 
                                      proj_params, grid_params, x0, y0,x_offset, y_offset, 
                                      spec_name, year_start, outputdir_main_ruisdael_area_csv2netc)
            resampler.run()
