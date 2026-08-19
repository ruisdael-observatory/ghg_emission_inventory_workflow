import fiona
import geopandas as gpd

bag_path = "bag-light.gpkg"

# List layers using fiona
layers = fiona.listlayers(bag_path)
print("Available layers in BAG Geopackage:", layers)

# Load the main layer (likely 'pand')
gdf_bag = gpd.read_file(bag_path, layer=layers[0])

SNAP3_KEYWORDS = [
    'industriefunctie',          # Core industry
    'bedrijf',                   # General business (catch-all for small industries)
    'opslagfunctie',             # Industrial storage/logistics
    'energieopwekking',          # Energy/production (if present)
    'nutsvoorziening'            # Utility buildings (could be indirect)
]

def flag_snap3_buildings(gdf):
    """
    Flag BAG buildings relevant for SNAP3 industrial emissions.

    Parameters
    ----------
    gdf : GeoDataFrame
        BAG building layer with column 'gebruiksdoel'.

    Returns
    -------
    GeoDataFrame
        Same gdf with an extra boolean column 'is_snap3' = True/False.
    """
    import pandas as pd
    gdf = gdf.copy()

    # Normalize: lower-case, fill NAs
    gdf['gebruiksdoel'] = gdf['gebruiksdoel'].fillna('').str.lower()

    # Flag if any keyword in SNAP3_KEYWORDS is present
    gdf['is_snap3'] = gdf['gebruiksdoel'].apply(
        lambda x: any(kw in x for kw in SNAP3_KEYWORDS)
    )

    return gdf


def flag_snap3_buildings(gdf):
    """
    Flag BAG buildings relevant for SNAP3 industrial emissions.

    Parameters
    ----------
    gdf : GeoDataFrame
        BAG building layer with column 'gebruiksdoel'.

    Returns
    -------
    GeoDataFrame
        Same gdf with an extra boolean column 'is_snap3' = True/False.
    """
    import pandas as pd
    gdf = gdf.copy()

    # Normalize: lower-case, fill NAs
    gdf['gebruiksdoel'] = gdf['gebruiksdoel'].fillna('').str.lower()

    # Flag if any keyword in SNAP3_KEYWORDS is present
    gdf['is_snap3'] = gdf['gebruiksdoel'].apply(
        lambda x: any(kw in x for kw in SNAP3_KEYWORDS)
    )

    return gdf

snap3_gdf_flagged = flag_snap3_buildings(gdf_bag)

# Filter only SNAP3-relevant buildings
snap3_gdf = snap3_gdf_flagged[snap3_gdf_flagged['is_snap3']]
print(f"Selected {len(snap3_gdf)} SNAP3-related buildings out of {len(gdf_bag)} total.")

snap3_gdf.to_file('bag_snap3.gpkg', driver='GPKG')