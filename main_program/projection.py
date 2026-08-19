"""
This is a helper class used to transform the grid

Please note that currently only RD coordinates are supported for the workflow output emissions..

Author: Dr. Arseni Doyennel
Email: a.doyennel@vu.nl
"""

from pyproj import Transformer, CRS

class Transform:
    def __init__(self, parameters):
        self.parameters = parameters
        self.crs_latlon = 'epsg:4326'
        self.crs_rd = 'epsg:28992'
        self.latlon_to_xy_transform = Transformer.from_crs(self.crs_latlon, self.parameters['proj4'])
        self.xy_to_latlon_transform = Transformer.from_crs(self.parameters['proj4'], self.crs_latlon)
        self.rd_to_latlon_transform = Transformer.from_crs(self.crs_rd, self.crs_latlon)
        self.latlon_to_rd_transform = Transformer.from_crs(self.crs_latlon, self.crs_rd)
        self.rd_to_lcc_transform = Transformer.from_crs(self.crs_rd, self.parameters['proj4'])

    def latlon_to_xy(self, lat, lon):
        return self.latlon_to_xy_transform.transform(lat, lon)

    def xy_to_latlon(self, x, y):
        return self.xy_to_latlon_transform.transform(x, y)

    def rd_to_latlon(self, x, y):
        return self.rd_to_latlon_transform.transform(x, y)

    def latlon_to_rd(self, lat, lon):
        return self.latlon_to_rd_transform.transform(lat, lon)

    def rd_to_lcc(self, x, y):
        return self.rd_to_lcc_transform.transform(x, y)