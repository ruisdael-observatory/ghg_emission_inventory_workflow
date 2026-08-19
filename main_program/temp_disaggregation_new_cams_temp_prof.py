"""
CAMS-based temporal profile generator for Dutch emission modelling.

This module generates temporal disaggregation coefficients for atmospheric
emission inventories over the Netherlands using CAMS-REG-TEMPO temporal
profiles. It replaces or supplements the original Van der Gon (2011)
SNAP-based temporal profiles with updated spatially representative profiles.

The module provides the following functionality:

    - Reads CAMS-REG-TEMPO NetCDF temporal profile files.
    - Automatically selects the requested CAMS reference year or, if
      unavailable, uses the latest available year.
    - Converts Dutch RD New (EPSG:28992) emission domains into geographic
      coordinates (WGS84 latitude/longitude).
    - Extracts CAMS temporal profiles representative for the emission domain.
    - Supports species-specific profiles (e.g. CO2, CH4).
    - Calculates monthly temporal factors from CAMS daily profiles.
    - Converts CAMS hourly, weekly and monthly profiles into the same
      SNAP structure used by the previous emission processing system.

Output profiles have the same structure as the original Van der Gon
implementation:

    tprof_hour : numpy.ndarray
        Hourly temporal profiles with dimensions:
        (10 SNAP categories, 24 hours)

    tprof_week : numpy.ndarray
        Weekly temporal profiles with dimensions:
        (10 SNAP categories, 7 days)

    tprof_mnth : numpy.ndarray
        Monthly temporal profiles with dimensions:
        (10 SNAP categories, 12 months)


The module uses CAMS profiles where available and automatically falls back
to the original Van der Gon (2011) coefficients for missing sectors or
species-specific information. This ensures compatibility with existing
emission processing workflows while allowing updated Dutch-specific temporal
variability.

The temporal profiles are intended to be applied during emission file
generation, where the emission timestamp is used to select the corresponding
hour, weekday and month factors.

Example workflow:

    CO2/CH4 annual emissions
              |
              v
    loadsnap_cams()
              |
              v
    tprof_hour, tprof_week, tprof_mnth
              |
              v
    writereademission_3d()
              |
              v
    hourly emission fields


References
----------
Denier van der Gon, H. A. C. et al. (2011).
A European temporal emission profile dataset for atmospheric modelling.

Guevara, M., Jorba, O., Perez Garcia-Pando, C. et al. (2021).
Time-resolved emission inventories for atmospheric chemistry modelling.
Earth System Science Data, 13, 367–404.

Dataset:
CAMS-REG-TEMPO temporal profiles (BSC).
"""
import os
import glob
import warnings

import numpy as np
import pandas as pd
import xarray as xr

from pyproj import Transformer
import calendar

from emission_construction_functions import loadsnap
from emission_preparation_setting import temp_prof_file

# =============================================================================
# Part 1 - Imports and helper functions
#
# Utilities for reading CAMS temporal profiles for the Netherlands.
#
# The functions in this module:
#   - convert Dutch RD coordinates to WGS84 latitude/longitude;
#   - locate the appropriate CAMS temporal-profile NetCDF file;
#   - normalize temporal profiles;
#   - provide automatic fallback to legacy (Van der Gon, 2011) profiles.
#
# =============================================================================


# =============================================================================
# Projection transformer
#
# EPSG:28992 : Dutch RD New
# EPSG:4326  : WGS84 latitude/longitude
# =============================================================================

rd2ll = Transformer.from_crs(
    "EPSG:28992",
    "EPSG:4326",
    always_xy=True
)


# =============================================================================
def normalize_profile(profile):
# =============================================================================
    """
    Normalize a temporal profile so that its mean equals one.

    Parameters
    ----------
    profile : array_like

    Returns
    -------
    ndarray
        Normalized profile.
    """

    profile = np.asarray(profile, dtype=float)

    if np.all(np.isnan(profile)):
        raise ValueError("Profile contains only NaN values.")

    profile = np.where(
        np.isnan(profile),
        np.nanmean(profile),
        profile
    )

    mean = profile.mean()

    if np.isclose(mean, 0.0):
        raise ValueError("Cannot normalize profile with zero mean.")

    return profile / mean


# =============================================================================
def safe_profile(profile, fallback):
# =============================================================================
    """
    Return a normalized CAMS profile if available,
    otherwise return the legacy profile.

    Parameters
    ----------
    profile : ndarray or None
        CAMS profile.

    fallback : ndarray
        Legacy (Van der Gon, 2011) profile.

    Returns
    -------
    ndarray
    """

    if profile is None:
        return fallback.copy()

    profile = np.asarray(profile)

    if np.all(np.isnan(profile)):
        return fallback.copy()

    return normalize_profile(profile)


# =============================================================================
def rd_domain_to_latlon(x0, y0, dx, dy, nx, ny):
# =============================================================================
    """
    Construct an RD grid and convert it to latitude/longitude.

    Parameters
    ----------
    x0, y0 : float
        Lower-left RD coordinates [m].

    dx, dy : float
        Grid spacing [m].

    nx, ny : int
        Number of grid cells.

    Returns
    -------
    lon : ndarray (ny,nx)

    lat : ndarray (ny,nx)
    """

    x = x0 + dx * (np.arange(nx) + 0.5)
    y = y0 + dy * (np.arange(ny) + 0.5)

    xx, yy = np.meshgrid(x, y)

    lon, lat = rd2ll.transform(xx, yy)

    return lon, lat


# =============================================================================
def get_available_year(filename):
# =============================================================================
    """
    Read the year represented by a CAMS NetCDF file.

    Parameters
    ----------
    filename : str

    Returns
    -------
    int
    """

    ds = xr.open_dataset(
        filename,
        decode_times=True
    )

    year = pd.to_datetime(ds.time.values[0]).year

    ds.close()

    return int(year)


# =============================================================================
def select_cams_file(requested_year, path_pattern, cams_file_name):
# =============================================================================
    """
    Select the CAMS temporal-profile file matching the requested year.

    If the requested year is unavailable,
    the newest available file is used.

    Parameters
    ----------
    requested_year : int

    path_pattern : str

        Example
        -------
        "/data/CAMS-REG-TEMPO_*_*.nc"

    Returns
    -------
    filename : str

    year_used : int
    """

    files = sorted(
        glob.glob(
            os.path.join(path_pattern, cams_file_name)
        )
    )

    if len(files) == 0:
        raise FileNotFoundError(
            f"No CAMS files found matching:\n{path_pattern}"
        )

    years = np.array(
        [get_available_year(f) for f in files]
    )

    if requested_year in years:

        idx = np.where(years == requested_year)[0][0]

        return files[idx], requested_year

    idx = np.argmax(years)

    latest = years[idx]

    warnings.warn(
        f"CAMS year {requested_year:d} not available. "
        f"Using latest available year ({latest:d})."
    )

    return files[idx], latest


# =============================================================================
# Part 2 - CAMS data reading and spatial extraction
#
# Functions for extracting Netherlands-specific temporal profiles
# from CAMS-REG-TEMPO gridded files.
#
# =============================================================================


# =============================================================================
def open_cams_dataset(filename):
# =============================================================================
    """
    Open CAMS temporal profile NetCDF file.

    Parameters
    ----------
    filename : str

    Returns
    -------
    xarray.Dataset
    """

    ds = xr.open_dataset(
        filename,
        decode_times=True
    )

    return ds



# =============================================================================
def spatial_mask(ds, lon, lat):
# =============================================================================
    """
    Create spatial mask selecting CAMS grid cells
    corresponding to the emission domain.

    Parameters
    ----------
    ds : xarray.Dataset

    lon : ndarray
        Longitude of model grid.

    lat : ndarray
        Latitude of model grid.

    Returns
    -------
    mask : ndarray
    """

    lon_min = np.nanmin(lon)
    lon_max = np.nanmax(lon)

    lat_min = np.nanmin(lat)
    lat_max = np.nanmax(lat)


    mask = (
        (ds.longitude >= lon_min) &
        (ds.longitude <= lon_max) &
        (ds.latitude >= lat_min) &
        (ds.latitude <= lat_max)
    )

    return mask



# =============================================================================
def spatial_mean(ds, variable, lon=None, lat=None):
# =============================================================================
    """
    Calculate spatial mean of a CAMS temporal profile.

    Parameters
    ----------
    ds : xarray.Dataset

    variable : str
        CAMS variable name.

    lon, lat : ndarray
        Optional emission-domain coordinates.

    Returns
    -------
    ndarray

        Temporal profile.

    """

    if variable not in ds.variables:
        return None


    da = ds[variable]


    # -------------------------------------------------------------
    # Select Dutch emission domain if provided
    # -------------------------------------------------------------

    if lon is not None and lat is not None:

        mask_lon = (
            (ds.longitude >= np.nanmin(lon)) &
            (ds.longitude <= np.nanmax(lon))
        )

        mask_lat = (
            (ds.latitude >= np.nanmin(lat)) &
            (ds.latitude <= np.nanmax(lat))
        )

        da = da.where(
            mask_lon & mask_lat,
            drop=True
        )


    # -------------------------------------------------------------
    # Average spatial dimensions only
    # -------------------------------------------------------------

    spatial_dims = [
        d for d in da.dims
        if d in ["latitude", "longitude"]
    ]


    if len(spatial_dims) > 0:

        da = da.mean(
            dim=spatial_dims,
            skipna=True
        )


    return da.values



# =============================================================================
def extract_cams_profiles(ds, species, lon=None, lat=None):
# =============================================================================
    """
    Extract available CAMS temporal profiles.

    Parameters
    ----------
    ds : xarray.Dataset

    species : str
        Example:
        "CO2"
        "CH4"

    lon, lat : ndarray
        Netherlands emission domain.

    Returns
    -------
    profiles : dict

    """

    species = species.lower()


    profiles = {}



    # =============================================================
    # Energy / Public power
    # =============================================================

    profiles["energy_hour"] = spatial_mean(
        ds,
        f"FH_A_ch4" if species == "ch4"
        else f"FH_A_co2",
        lon,
        lat
    )


    profiles["energy_week"] = spatial_mean(
        ds,
        f"FW_A_ch4" if species == "ch4"
        else f"FW_A_co2",
        lon,
        lat
    )



    # =============================================================
    # Residential / stationary combustion
    # =============================================================

    profiles["residential_hour"] = spatial_mean(
        ds,
        "FH_C_others",
        lon,
        lat
    )



    # =============================================================
    # Road transport
    # =============================================================

    profiles["transport_week"] = spatial_mean(
        ds,
        "FW_F",
        lon,
        lat
    )


    profiles["transport_weekday_hour"] = spatial_mean(
        ds,
        "FH_weekday_F",
        lon,
        lat
    )


    profiles["transport_saturday_hour"] = spatial_mean(
        ds,
        "FH_saturday_F",
        lon,
        lat
    )


    profiles["transport_sunday_hour"] = spatial_mean(
        ds,
        "FH_sunday_F",
        lon,
        lat
    )



    # =============================================================
    # Agriculture
    #
    # CAMS does not provide CH4/CO2 agriculture profiles.
    # Livestock NH3/NOx profile is used only as CH4 approximation.
    # =============================================================

    if species == "ch4":

        profiles["agriculture_day"] = spatial_mean(
            ds,
            "FD_K_nh3_nox",
            lon,
            lat
        )

    else:

        profiles["agriculture_day"] = None
        
    
    # =============================================================
    # Aviation
    # =============================================================

    profiles["aviation_week"] = spatial_mean(
        ds,
        "FW_H",
        lon,
        lat
    )

    profiles["aviation_hour"] = None
    profiles["aviation_day"] = None



    return profiles

# =============================================================================
# Part 3 - Conversion of CAMS profiles to SNAP temporal coefficients
#
# Convert CAMS hourly, weekly and daily profiles into:
#
#     tprof_hour : 10 SNAP categories x 24 hours
#     tprof_week : 10 SNAP categories x 7 days
#     tprof_mnth : 10 SNAP categories x 12 months
#
# =============================================================================

# =============================================================================

def daily_to_monthly(
        daily_profile,
        profile_year,
        target_year
):
# =============================================================================
    """
    Convert CAMS daily weights into monthly temporal factors.

    Parameters
    ----------
    daily_profile : ndarray
        Daily CAMS weights (365 or 366 values).

    year : int
        Reference year.

    Returns
    -------
    ndarray
        12 monthly normalized factors.
    """

    """
    Convert CAMS daily weights to monthly factors.

    Parameters
    ----------
    daily_profile : ndarray
        CAMS daily weights.

    profile_year : int
        Year of CAMS file.

    target_year : int
        Simulation year.

    Returns
    -------
    ndarray
        12 monthly normalized coefficients.
    """

    if daily_profile is None:
        return None


    daily_profile = np.asarray(
        daily_profile,
        dtype=float
    )


    dates = pd.date_range(
        start=f"{profile_year}-01-01",
        periods=len(daily_profile),
        freq="D"
    )


    df = pd.DataFrame(
        {
            "weight": daily_profile
        },
        index=dates
    )


    # ---------------------------------------------------------
    # Remove leap day when CAMS year is leap
    # and simulation year is not
    # ---------------------------------------------------------

    if (
        calendar.isleap(profile_year)
        and not calendar.isleap(target_year)
    ):

        df = df[
            ~(
                (df.index.month == 2)
                &
                (df.index.day == 29)
            )
        ]


    monthly = (
        df.groupby(
            df.index.month
        )["weight"]
        .mean()
        .values
    )


    # Make sure shape is (12,)
    monthly = monthly.flatten()


    return normalize_profile(monthly)


# =============================================================================
def combine_transport_hours(
        weekday,
        saturday,
        sunday
):
# =============================================================================
    """
    Create one representative weekly hourly profile
    for road transport.

    Assumption:
        5 weekdays + Saturday + Sunday.

    Parameters
    ----------
    weekday : ndarray
        24 hourly weights.

    saturday : ndarray

    sunday : ndarray


    Returns
    -------
    ndarray
        24 hourly weights.
    """

    if (
        weekday is None or
        saturday is None or
        sunday is None
    ):
        return None


    profile = (
        5.0 * np.asarray(weekday)
        +
        np.asarray(saturday)
        +
        np.asarray(sunday)
    ) / 7.0


    return normalize_profile(profile)


# =============================================================================
def cams_to_snap_profiles(
        profiles,
        year,
        requested_year
):
# =============================================================================
    """
    Convert extracted CAMS profiles into SNAP arrays.

    Parameters
    ----------
    profiles : dict
        Output from extract_cams_profiles()

    year : int
        CAMS reference year.

    Returns
    -------
    tprof_hour
    tprof_week
    tprof_mnth

    """


    nsnap = 11


    tprof_hour = np.full((nsnap,24), np.nan)
    tprof_week = np.full((nsnap,7), np.nan)
    tprof_mnth = np.full((nsnap,12), np.nan)


    # =============================================================
    # SNAP 1 - Energy production
    # =============================================================

    if profiles["energy_hour"] is not None:

        tprof_hour[0] = normalize_profile(
            profiles["energy_hour"]
        )


    if profiles["energy_week"] is not None:

        tprof_week[0] = normalize_profile(
            profiles["energy_week"]
        )



    # =============================================================
    # SNAP 2 - Residential combustion
    # =============================================================

    if profiles["residential_hour"] is not None:

        tprof_hour[1] = normalize_profile(
            profiles["residential_hour"]
        )



    # =============================================================
    # SNAP 7 - Road transport
    # =============================================================

    transport_hour = combine_transport_hours(
        profiles["transport_weekday_hour"],
        profiles["transport_saturday_hour"],
        profiles["transport_sunday_hour"]
    )


    if transport_hour is not None:

        tprof_hour[6] = transport_hour


    if profiles["transport_week"] is not None:

        tprof_week[6] = normalize_profile(
            profiles["transport_week"]
        )



    # =============================================================
    # SNAP 10 - Agriculture
    # =============================================================

    if profiles["agriculture_day"] is not None:

        tprof_mnth[9] = daily_to_monthly(
            profiles["agriculture_day"],
            year,
            requested_year
        )
        
    
    # =============================================================
    # SNAP 11 - Aviation
    # =============================================================

    if profiles["aviation_week"] is not None:

        tprof_week[10] = normalize_profile(
            profiles["aviation_week"]
        )


    # No CAMS hourly aviation profile available.
    # Use constant hourly distribution.

    tprof_hour[10] = np.ones(24)


    # No CAMS monthly aviation profile available.
    # Assume uniform annual distribution.

    tprof_mnth[10] = np.ones(12)

    return (
        tprof_hour,
        tprof_week,
        tprof_mnth
    )

# =============================================================================
# Part 4 - Merge CAMS profiles with legacy Van der Gon profiles
#
# CAMS profiles replace legacy profiles where available.
# Missing CAMS information automatically falls back to the
# original Denier van der Gon (2011) SNAP profiles.
#
# =============================================================================


# =============================================================================
def merge_cams_with_legacy(
        cams_hour,
        cams_week,
        cams_month
):
# =============================================================================
    """
    Merge CAMS-derived SNAP profiles with the original profiles.

    Parameters
    ----------
    cams_hour : ndarray
        CAMS hourly profiles (10,24)

    cams_week : ndarray
        CAMS weekly profiles (10,7)

    cams_month : ndarray
        CAMS monthly profiles (10,12)


    Returns
    -------
    tprof_hour
    tprof_week
    tprof_mnth

    """


    # -------------------------------------------------------------
    # Load original Van der Gon profiles
    # -------------------------------------------------------------

    old_hour, old_week, old_month = loadsnap()



    # -------------------------------------------------------------
    # Start from legacy profiles
    # -------------------------------------------------------------

    tprof_hour = old_hour.copy()
    tprof_week = old_week.copy()
    tprof_mnth = old_month.copy()



    # -------------------------------------------------------------
    # Replace only valid CAMS profiles
    # -------------------------------------------------------------

    for snap in range(11):

        if cams_hour is not None:

            if not np.all(np.isnan(cams_hour[snap])):

                tprof_hour[snap] = normalize_profile(
                    cams_hour[snap]
                )


        if cams_week is not None:

            if not np.all(np.isnan(cams_week[snap])):

                tprof_week[snap] = normalize_profile(
                    cams_week[snap]
                )


        if cams_month is not None:

            if not np.all(np.isnan(cams_month[snap])):

                tprof_mnth[snap] = normalize_profile(
                    cams_month[snap]
                )


    return (
        tprof_hour,
        tprof_week,
        tprof_mnth
    )




# =============================================================================
# Part 5 - Main CAMS temporal profile loader
#
# Main user function.
#
# Returns the same output structure as the original loadsnap():
#
#     tprof_hour : (10,24)
#     tprof_week : (10,7)
#     tprof_mnth : (10,12)
#
# =============================================================================


# =============================================================================
def loadsnap_cams(
        species,
        year,
        cams_path,
        cams_file_name,
        x0,
        y0,
        dx,
        dy,
        nx,
        ny
):
# =============================================================================
    """
    Generate SNAP temporal profiles using CAMS-REG-TEMPO.

    CAMS profiles are used where available.
    Original Van der Gon (2011) profiles are used as fallback.

    Parameters
    ----------
    species : str

        Species name.

        Examples:
            "CO2"
            "CH4"


    year : int

        Simulation year.


    cams_path : str

        CAMS NetCDF path or wildcard.

        Example:
            "/data/CAMS-REG-TEMPO_*_*.nc"


    x0, y0 : float

        Lower-left corner of Dutch RD grid [m].


    dx, dy : float

        Grid resolution [m].


    nx, ny : int

        Number of grid cells.


    Returns
    -------
    tprof_hour : ndarray
        (10,24)

    tprof_week : ndarray
        (10,7)

    tprof_mnth : ndarray
        (10,12)

    """


    # -------------------------------------------------------------
    # Species formatting
    # -------------------------------------------------------------

    species = species.lower()


    if species not in [
        "co2",
        "ch4"
    ]:
        raise ValueError(
            "Species must be 'CO2' or 'CH4'"
        )


    # -------------------------------------------------------------
    # Select CAMS file
    # -------------------------------------------------------------

    cams_file, cams_year = select_cams_file(
        year,
        cams_path,
        cams_file_name
    )


    print(
        f"Using CAMS temporal profiles: {cams_file}"
    )

    print(
        f"CAMS reference year: {cams_year}"
    )

    print(
        f"Emission year: {year}"
    )



    # -------------------------------------------------------------
    # Convert RD domain to latitude/longitude
    # -------------------------------------------------------------

    lon, lat = rd_domain_to_latlon(
        x0,
        y0,
        dx,
        dy,
        nx,
        ny
    )



    # -------------------------------------------------------------
    # Read CAMS data
    # -------------------------------------------------------------

    ds = open_cams_dataset(
        cams_file
    )



    # -------------------------------------------------------------
    # Extract CAMS temporal profiles
    # -------------------------------------------------------------

    cams_profiles = extract_cams_profiles(
        ds,
        species,
        lon,
        lat
    )


    ds.close()



    # -------------------------------------------------------------
    # Create SNAP arrays from CAMS
    # -------------------------------------------------------------

    (
        cams_hour,
        cams_week,
        cams_month

    ) = cams_to_snap_profiles(
        cams_profiles,
        cams_year,
        year
    )



    # -------------------------------------------------------------
    # Special handling for agriculture
    #
    # CAMS agriculture daily profile is only suitable
    # as CH4 livestock approximation.
    #
    # For CO2 keep original Van der Gon agriculture.
    # -------------------------------------------------------------

    if species == "co2":

        cams_month[9] = np.nan



    # -------------------------------------------------------------
    # Merge CAMS + legacy profiles
    # -------------------------------------------------------------

    (
        tprof_hour,
        tprof_week,
        tprof_mnth

    ) = merge_cams_with_legacy(
        cams_hour,
        cams_week,
        cams_month
    )



    return (
        tprof_hour,
        tprof_week,
        tprof_mnth
    )



#tprof_hour, tprof_week, tprof_mnth = loadsnap_cams(
#    species="co2",
#    year=2022,
#    cams_path='',
#    cams_file=temp_prof_file,
#    x0 = 104000.0,
#    y0 = 425910.0,
#    dx=100,
#    dy=100,
#    nx=384,
#    ny=384
#)

