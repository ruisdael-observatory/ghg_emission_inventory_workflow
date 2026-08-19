import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pyproj import Transformer, CRS
from GridDales import GridDales
from shapely.geometry import box, Point
import fiona
import geopandas as gpd
import os
from emission_preparation_setting import spec_name, proj_params, grid_params, x0, y0, x_offset, y_offset, rundir_ruisdael_area_residuals, noxdir, outputdir_main_ruisdael_area_csv2netc, year_start, snaplist, nox_multistring_dataset
from netCDF4 import Dataset
from projection import Transform

class Traffic_Disaggregator:
    def __init__(self):
        self.nox_path = noxdir
        self.spec_name = spec_name
        self.year = year_start
        self.data_csv_path = f"{rundir_ruisdael_area_residuals}{self.spec_name}_{self.year}_SNAP_7_residual.csv"
        self.nox_multistring_dataset=nox_multistring_dataset
        self.proj_params = proj_params
        self.transformer = Transform(self.proj_params)
        self.x_offset = x_offset
        self.y_offset = y_offset
        self.lcc_start_x = x0+self.x_offset
        self.lcc_start_y = y0+self.y_offset
        self.grid_params = grid_params
        # Step 1: Create LES Grid
        self.grid_params['southwest_x'] = self.lcc_start_x
        self.grid_params['southwest_y'] = self.lcc_start_y
        self.les_grid = GridDales(self.grid_params)
        self.les_grid.xt += self.grid_params['southwest_x']
        self.les_grid.yt += self.grid_params['southwest_y']
        self.output_folder = outputdir_main_ruisdael_area_csv2netc
        self.resampled_data = None  # 2D np.array for resampled data
        self.lat = None             # 2D np.array of latitudes
        self.lon = None             # 2D np.array of longitudes

    def read_nox_data(self):
        cols = ['Shape_Length', 'TVtot_N_sum']
        with fiona.open(os.path.join(self.nox_path, self.nox_multistring_dataset)) as src:
            features = [{'geometry': f['geometry'],
                         'properties': {k: f['properties'][k] for k in cols}} for f in src]
        gdf = gpd.GeoDataFrame.from_features(features)
        gdf['TVtot_N_sum'] = gdf['TVtot_N_sum'] / 1000
        gdf.crs = 'EPSG:28992'
        return gdf

    def calculate_road_lengths(self, gdf):
        gdf['orig_len'] = gdf.geometry.length
        return gdf

    def create_grid(self, bounds, dx, dy, crs):
        xmin, ymin, xmax, ymax = bounds
        cols = np.arange(xmin, xmax, dx)
        rows = np.arange(ymin, ymax, dy)
        polygons = [box(x, y, x + dx, y + dy) for x in cols for y in rows]
        grid = gpd.GeoDataFrame({'geometry': polygons})
        grid.crs = crs
        return grid

    def get_bounds_from_les_grid(self, les_grid):
        dx, dy = les_grid.dx, les_grid.dy
        return (
            les_grid.xt[0] - dx / 2,
            les_grid.yt[0] - dy / 2,
            les_grid.xt[-1] + dx / 2,
            les_grid.yt[-1] + dy / 2
        )

    def mass_conserving_rasterize(self, road_network_lcc, grid, value_col='TVtot_N_sum'):


        # Step 1: Setup
        road_network = road_network_lcc.copy()
        if 'line_id' not in road_network.columns:
            road_network['line_id'] = road_network.index
        if 'orig_len' not in road_network.columns:
            road_network['orig_len'] = road_network.geometry.length

        grid = grid.copy()
        grid['cell_id'] = grid.index

        # Step 2: Spatial join to find which roads intersect which grid cells
        joined = gpd.sjoin(
            road_network[['line_id', value_col, 'orig_len', 'geometry']],
            grid[['cell_id', 'geometry']],
            how='inner',
            predicate='intersects'
        )
        print(f"Number of road segments intersecting grid cells: {len(joined)}")

        if joined.empty:
            print("No intersections found. Returning empty GeoDataFrame.")
            return grid.assign(weighted_val=0.0)

        # Step 3: Rename geometry columns to avoid conflict
        joined = joined.rename(columns={'geometry': 'road_geom'})
        grid = grid.rename(columns={'geometry': 'grid_geom'})

        # Step 4: Create true intersections
        joined['geometry'] = joined.apply(lambda row: row['road_geom'].intersection(grid.loc[row['cell_id'], 'grid_geom']), axis=1)
        intersections = gpd.GeoDataFrame(joined, geometry='geometry', crs=road_network.crs)
        print(f"Number of intersections found: {len(intersections)}")

        # Step 5: Drop empty geometries
        intersections = intersections[~intersections.geometry.is_empty & intersections.geometry.is_valid]

        # Step 6: Compute fragment lengths and mass weights
        intersections['frag_len'] = intersections.geometry.length
        intersections['weight'] = intersections['frag_len'] / intersections['orig_len']
        intersections['weighted_val'] = intersections[value_col] * intersections['weight']

        # Step 7: Aggregate NOx by cell
        cell_vals = intersections.groupby('cell_id')['weighted_val'].sum().reset_index()

        # Step 8: Restore grid geometry and merge results
        grid = grid.rename(columns={'grid_geom': 'geometry'})
        result = grid[['cell_id', 'geometry']].merge(cell_vals, on='cell_id', how='left')
        result['weighted_val'] = result['weighted_val'].fillna(0.0)

        result = gpd.GeoDataFrame(result, geometry='geometry', crs=road_network.crs)
        print("✅ Mass-conserving rasterization complete.")

        return result

    def run(self):

        # Step 2: Load and prepare data
        road_network = self.read_nox_data()
        road_network = self.calculate_road_lengths(road_network)
        road_network_lcc = road_network.copy()

        dx = self.grid_params['xsize'] / self.grid_params['itot']
        xmin, ymin, xmax, ymax = self.get_bounds_from_les_grid(self.les_grid)
        grid = self.create_grid((xmin, ymin, xmax, ymax), dx, dx, crs=road_network.crs)

        

        mass_conserved = self.mass_conserving_rasterize(road_network_lcc, grid)
        clipped_total = road_network_lcc[road_network_lcc.intersects(box(xmin, ymax, xmax, ymin))]['TVtot_N_sum'].sum()
        print("Original vector total NOx in domain:", clipped_total)
        print("Mass-conserved rasterized NOx total:", mass_conserved['weighted_val'].sum())

        # Step 3: Prepare CO2 data
        data_df = pd.read_csv(self.data_csv_path)
        data_df['x_lcc'], data_df['y_lcc'] = data_df['x'], data_df['y']
        geometry = [Point(xy) for xy in zip(data_df['x_lcc'], data_df['y_lcc'])]
        crs_lcc = CRS.from_proj4(self.proj_params['proj4'])
        data_gdf = gpd.GeoDataFrame(data_df, geometry=geometry, crs=crs_lcc)

        # Step 4: Rasterize CO₂ using NOx weights
        les_x_1d, les_y_1d = self.les_grid.xt, self.les_grid.yt
        les_x, les_y = np.meshgrid(les_x_1d, les_y_1d)
        data_raster = np.zeros_like(les_x)

        # Ensure that the number of values in 'weighted_val' matches the number of grid cells
        n_grid_cells = len(les_x[:,1]) * len(les_y[1,:])

        # Verify the length of weighted_val
        weighted_vals = mass_conserved['weighted_val'].values

        # If the size of the data does not match the grid size, we can't reshape directly
        if len(weighted_vals) != n_grid_cells:
            raise ValueError(f"Number of 'weighted_val' data points ({len(weighted_vals)}) does not match grid size ({n_grid_cells}).")

        # Reshape the weighted_vals to match the grid shape, and transpose to correct the orientation
        weighted_vals = weighted_vals.reshape(len(les_x_1d), len(les_y_1d)).T # Transpose to match x, y orientation

        for _, row in data_gdf.iterrows():
            x_c, y_c, val = row['x_lcc'], row['y_lcc'], row[self.spec_name]

            # Find grid indices that fall within the 1000m cell centered at (x_c, y_c)
            x_min, x_max = x_c - 500, x_c + 500
            y_min, y_max = y_c - 500, y_c + 500

            mask_x = (les_x_1d >= x_min) & (les_x_1d < x_max)
            mask_y = (les_y_1d >= y_min) & (les_y_1d < y_max)

            x_indices = np.where(mask_x)[0]
            y_indices = np.where(mask_y)[0]

            if x_indices.size > 0 and y_indices.size > 0:
                sub_weights = weighted_vals[np.ix_(y_indices, x_indices)]
                total_weight = sub_weights.sum()

                if total_weight > 0:
                    data_raster[np.ix_(y_indices, x_indices)] += val * (sub_weights / total_weight)
                    
            self.resampled_data=data_raster

        domain_bounds = box(les_x_1d.min(), les_y_1d.min(), les_x_1d.max(), les_y_1d.max())
        data_clipped = data_gdf[data_gdf.intersects(domain_bounds)]
        print(f"Original {self.spec_name} total (clipped to domain): {data_clipped[self.spec_name].sum():.4f}")
        print(f"Disaggregated {self.spec_name} total (from raster):  {data_raster.sum():.4f}")
        print(f"Relative error: {(data_raster.sum() - data_clipped[self.spec_name].sum()) / data_clipped[self.spec_name].sum():.2%}")

        # Visualization
        fig, ax = plt.subplots(figsize=(10, 10))
        mesh = ax.pcolormesh(les_x, les_y, data_raster, shading='auto', cmap='viridis', vmin=1e0, vmax=1e1)
        plt.colorbar(mesh, ax=ax, label=f'{self.spec_name} (kg/year)')
        ax.set_title(f"Refined {self.spec_name} Emissions Using NOx Weights")
        ax.set_xlabel("X RD")
        ax.set_ylabel("Y RD")
        ax.set_aspect('equal')
        plt.grid(True)
        plt.tight_layout()
        plt.show()
        
    def save_to_netcdf(self, filename=None):

        x_les,y_les=np.meshgrid(self.les_grid.xt, self.les_grid.yt)
        # Convert LES x, y coordinates to lat/lon using your existing transformer
        self.lat, self.lon = self.transformer.xy_to_latlon(x_les,y_les)
        
        if filename is None:
            filename = os.path.join(self.output_folder, f'HARM_snap_7_all_{self.spec_name}.nc')
    
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
        x_var[:] = self.les_grid.xt
        y_var[:] = self.les_grid.yt
        lat_var[:, :] = self.lat
        lon_var[:, :] = self.lon
        spec_var[:, :] = self.resampled_data

        # Add global attribute
        ncfile.title = f'{self.spec_name} Emissions for {self.year} Data'

        ncfile.close()
        print(f'NetCDF file created successfully: {filename}')


if __name__ == '__main__':
    processor = Traffic_Disaggregator()
    processor.run()
    processor.save_to_netcdf()