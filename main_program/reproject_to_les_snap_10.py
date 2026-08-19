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
from emission_preparation_setting import spec_name, proj_params, grid_params, x0, y0, x_offset, y_offset, rundir_ruisdael_area_residuals, cbs_loc_consu_preprocessing, outputdir_main_ruisdael_area_csv2netc, year_start, snaplist, agro_file_name
import geopandas as gpd

class SNAP10Refiner:
    def __init__(self, snap, agro_file_link, agro_file_name, input_link, proj_params, grid_params, 
                 x0, y0, x_offset, y_offset, spec_name, year, output_folder):
        """
        agro_file_name: GeoDataFrame of Bouwland + Grasland
        snap_csv: path to SNAP10 CSV with columns x, y, co2/ch4
        refinement_res: optional resolution in meters for LES refinement
        """
        self.snap=snap
        self.agro_file_link = agro_file_link
        self.agro_file_name=agro_file_name
        self.gdf_target = gpd.read_file(f"{self.agro_file_link}{self.agro_file_name}")
        self.input_link=input_link
        self.spec_name = spec_name
        self.year = year
        self.input_file=f'{self.spec_name}_{self.year}_SNAP_10_residual.csv'
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
        self.df_snap_filtered = None
        self.brp_raster = None
        self.agri_fraction = None
        self.brp_transform = None
        self.raster_xmin = None
        self.raster_ymax = None
        self.raster_res = None
        self.resampled_data = None
        self.resampled_data_plot = None

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
        
    def plot_snap_mesh(self):
        df = self.df_snap_filtered

        x = df["lcc_x"].values
        y = df["lcc_y"].values
        v = df[self.spec_name].values

        # --- infer grid
        x_unique = np.sort(np.unique(x))
        y_unique = np.sort(np.unique(y))

        nx = len(x_unique)
        ny = len(y_unique)

        # --- map to grid indices
        xi = np.searchsorted(x_unique, x)
        yi = np.searchsorted(y_unique, y)

        grid = np.full((ny, nx), np.nan)

        grid[yi, xi] = v

        plt.figure(figsize=(8, 8))

        plt.pcolormesh(x_unique, y_unique, grid, vmax=5000, shading="auto")

        plt.colorbar(label="SNAP emissions")
        plt.title("SNAP emissions (structured grid view)")
        plt.xlabel("x")
        plt.ylabel("y")
        plt.axis("equal")
        
    def plot_proxy_field(self):
        """
        Plot agricultural proxy (agri_fraction) on LES grid
        BEFORE any emission redistribution.
        """

        if not hasattr(self, "agri_fraction"):
           print("[WARNING] agri_fraction not found")
           return

        proxy = self.agri_fraction

        nx = len(self.les_x)
        ny = len(self.les_y)

        if proxy.shape != (ny, nx):
           print(f"[WARNING] Shape mismatch: proxy {proxy.shape}, expected {(ny, nx)}")

        x = self.les_x
        y = self.les_y

        plt.figure(figsize=(8, 8))

        plt.pcolormesh(
           x,
           y,
           proxy,
           shading="auto",
           cmap="YlGn"
        )

        plt.colorbar(label="Agricultural fraction (proxy)")
        plt.title("Proxy field (agri_fraction) on LES grid")
        plt.xlabel("x (LES / RD)")
        plt.ylabel("y (LES / RD)")
        plt.axis("equal")

        plt.show()
        

    def rasterize_brp_subset(self, buffer_cells=2):
        """
        Rasterize only the subset of self.gdf_target overlapping the LES domain + buffer.
        Produce a fractional agro fraction raster (0–1) using rasterio in one fast pass.
        """
        nx, ny = len(self.les_x), len(self.les_y)
        dx, dy = self.les_grid.dx, self.les_grid.dy

        xmin = self.les_x[0] - dx/2
        xmax = self.les_x[-1] + dx/2
        ymin = self.les_y[0] - dy/2
        ymax = self.les_y[-1] + dy/2

        # Add buffer
        xmin_buf = xmin - buffer_cells * dx
        xmax_buf = xmax + buffer_cells * dx
        ymin_buf = ymin - buffer_cells * dy
        ymax_buf = ymax + buffer_cells * dy

        les_box_buf = box(xmin_buf, ymin_buf, xmax_buf, ymax_buf)

        # Subset polygons intersecting buffered LES domain
        gdf_subset = self.gdf_target[self.gdf_target.intersects(les_box_buf)].copy()
        if gdf_subset.empty:
            print("[WARNING] No polygons intersect LES domain with buffer!")
            self.brp_raster = np.zeros((ny, nx), dtype='float32')
            self.agri_fraction = np.zeros((ny, nx), dtype='float32')
            return

        # Clip polygons to buffered domain
        gdf_subset['geometry'] = gdf_subset.geometry.intersection(les_box_buf)
        gdf_subset = gdf_subset[~gdf_subset.is_empty]

        transform = rasterio.transform.from_origin(xmin, ymax, dx, dy)

        # Rasterize using rasterio → fast fractional coverage
        # use merge_func=sum to count overlapping areas, then divide by max (cap at 1)
        shapes = ((geom, 1) for geom in gdf_subset.geometry)
        raster = rasterize(
            shapes,
            out_shape=(ny, nx),
            transform=transform,
            all_touched=True,
            dtype='float32',
            merge_alg=rasterio.enums.MergeAlg.add  # sum overlapping pixels
        )

        # Cap values to 1 for fractional coverage
        raster[raster > 1.0] = 1.0
        
        raster = np.flipud(raster) #the y-axis is mirrored here, so fix it...

        self.brp_raster = raster
        self.agri_fraction = raster
        self.brp_transform = transform
        self.raster_xmin = xmin
        self.raster_ymax = ymax
        self.raster_res = dx

        print(f"[INFO] Rasterized subset (fast precise): {len(gdf_subset)} polygons → LES grid {nx}x{ny}")
        
    def refine_snap10_vectorized(self):
        """
        Fully vectorized, mass-conserving redistribution of SNAP emissions
        onto LES grid using agro proxy weighting.

        No loops over LES cells. Only SNAP loop remains.
        """

        nx, ny = len(self.les_x), len(self.les_y)
        dx, dy = self.les_grid.dx, self.les_grid.dy

        mass_grid = np.zeros((ny, nx), dtype=np.float64)

        snap_cell_size = 1000.0

        les_x_min = self.les_x - dx / 2
        les_x_max = self.les_x + dx / 2
        les_y_min = self.les_y - dy / 2
        les_y_max = self.les_y + dy / 2

        Xmin = les_x_min[None, :]
        Xmax = les_x_max[None, :]
        Ymin = les_y_min[:, None]
        Ymax = les_y_max[:, None]

        total_input_mass = 0.0

        for _, row in tqdm(self.df_snap_filtered.iterrows(),
                       total=len(self.df_snap_filtered),
                       desc="Vectorized SNAP refinement"):

           x0 = row["lcc_x"]
           y0 = row["lcc_y"]
           m = row[self.spec_name]
  
           total_input_mass += m

           sxmin = x0 - snap_cell_size / 2
           sxmax = x0 + snap_cell_size / 2
           symin = y0 - snap_cell_size / 2
           symax = y0 + snap_cell_size / 2

           ix = np.maximum(0.0, np.minimum(sxmax, Xmax) - np.maximum(sxmin, Xmin))
           iy = np.maximum(0.0, np.minimum(symax, Ymax) - np.maximum(symin, Ymin))

           inter_area = ix * iy

           if not np.any(inter_area):
              continue

           w = inter_area * self.agri_fraction

           w_sum = np.sum(w)

           if w_sum <= 0:
               w = inter_area
               w_sum = np.sum(w)
               if w_sum <= 0:
                  continue

           w /= w_sum

           mass_grid += m * w

        self.resampled_data = mass_grid
        self.resampled_data_plot = mass_grid.copy()
        self.resampled_data_plot[self.resampled_data_plot == 0] = np.nan

        print(f"[INFO] Total input mass: {total_input_mass:,.3f}")
        print(f"[INFO] Total output mass: {mass_grid.sum():,.3f}")
        print(f"[INFO] Mass conservation error: "
            f"{100*(mass_grid.sum()-total_input_mass)/total_input_mass:.6f}%")

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
        plt.figure(figsize=(8, 8))
        plt.pcolormesh(self.les_x, self.les_y, self.resampled_data, cmap='viridis', vmax=250, shading='auto')
        plt.colorbar(label=f"{self.spec_name} (kg/cell/year)")
        plt.title(f"Refined SNAP10 {self.spec_name} emissions ({self.year})")
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
        print("[1] Loading SNAP10 data...")
        self.load_snap10()
        print("[2] Filtering SNAP10 to LES domain and plot...")
        self.filter_snap10()
        self.plot_snap_mesh()
        print("[3] Rasterizing agro grassland at refinement resolution and computing agricultural fraction...")
        self.rasterize_brp_subset()
        self.plot_proxy_field()
        print("[4] Refining SNAP10 emissions with mass-conserving redistribution...")
        self.refine_snap10_vectorized()
        if smooth:
            print("[5] Smoothing while preserving mass...")
            self.smooth_preserve_mass()
        print("[6] Plotting refined SNAP10...")
        self.plot_resampled()
        print("[7] Saving NetCDF...")
        self.save_to_netcdf()
        print("Finished SNAP10 refinement.")

if __name__ == "__main__":
    for snap in snaplist:
        if snap == 10:  
            print(f"Processing SNAP {snap}")
            resampler = SNAP10Refiner(snap, cbs_loc_consu_preprocessing, agro_file_name,rundir_ruisdael_area_residuals, 
                                      proj_params, grid_params, x0, y0,x_offset, y_offset, 
                                      spec_name, year_start, outputdir_main_ruisdael_area_csv2netc)
            resampler.run()