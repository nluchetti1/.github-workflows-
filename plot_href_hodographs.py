import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from metpy.plots import Hodograph
from metpy.units import units
import numpy as np
import datetime
import requests
import os
import sys

# --- Configuration ---
REGION = [-98, -74, 24, 38]   # Southeast US [West, East, South, North]
GRID_SPACING = 10             # Plot density (Higher = fewer hodographs)
FORECAST_HOUR = '01'          # Mean usually starts at f01
OUTPUT_DIR = "images"

# Levels to look for. 
REQUESTED_LEVELS = [1000, 925, 850, 700, 500, 250] * units.hPa

def get_latest_run_time():
    now = datetime.datetime.utcnow()
    if now.hour >= 15:
        run = '12'
        date = now
    elif now.hour >= 3:
        run = '00'
        date = now
    else:
        run = '12'
        date = now - datetime.timedelta(days=1)
    return date.strftime('%Y%m%d'), run

def download_href_mean(date_str, run, fhr):
    base_url = f"https://nomads.ncep.noaa.gov/pub/data/nccf/com/href/prod/href.{date_str}/ensprod"
    filename = f"href.t{run}z.conus.mean.f{fhr}.grib2"
    url = f"{base_url}/{filename}"
    
    print(f"Downloading: {url}")
    headers = {'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64)'}
    
    try:
        with requests.get(url, stream=True, timeout=60, headers=headers) as r:
            if r.status_code == 404:
                print(f"Error 404: File not found: {url}")
                sys.exit(1)
            r.raise_for_status()
            with open(filename, 'wb') as f:
                for chunk in r.iter_content(chunk_size=8192):
                    f.write(chunk)
        print(f"Success: {filename}")
        return filename
    except Exception as e:
        print(f"Download failed: {e}")
        sys.exit(1)

def plot_hodographs(filename, date_str, run, fhr):
    print("Reading GRIB data (U and V components separately)...")
    
    # FIX 1: Open U and V separately to avoid "DatasetBuildError" from mismatching variables (like w or dpt)
    try:
        # Load U Component
        ds_u = xr.open_dataset(filename, engine='cfgrib', 
                               filter_by_keys={'typeOfLevel': 'isobaricInhPa', 'shortName': 'u'})
        
        # Load V Component
        ds_v = xr.open_dataset(filename, engine='cfgrib', 
                               filter_by_keys={'typeOfLevel': 'isobaricInhPa', 'shortName': 'v'})
        
        # Merge them into one dataset
        ds = xr.merge([ds_u, ds_v])
        
    except Exception as e:
        print(f"Failed to read GRIB: {e}")
        sys.exit(1)

    # FIX 2: Do NOT use ds.sel(latitude=...) because lat/lon are 2D coordinates, not dimensions.
    # We will rely on the loop below to filter points spatially.

    # Check Available Levels
    file_levels = ds.isobaricInhPa.values
    levels_to_use = [l for l in REQUESTED_LEVELS.m if l in file_levels]
    
    print(f"Levels found: {levels_to_use}")
    if len(levels_to_use) < 3:
        print("Error: Not enough vertical levels for hodographs.")
        sys.exit(1)
        
    # Filter vertical levels
    ds = ds.sel(isobaricInhPa=levels_to_use)

    # Extract Data
    u = ds['u'].metpy.convert_units('kts')
    v = ds['v'].metpy.convert_units('kts')

    # Setup Plot
    fig = plt.figure(figsize=(18, 12))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.LambertConformal())
    ax.set_extent(REGION)
    
    ax.add_feature(cfeature.COASTLINE, linewidth=1.5)
    ax.add_feature(cfeature.BORDERS, linewidth=1.5)
    ax.add_feature(cfeature.STATES, linewidth=0.5, edgecolor='gray')
    ax.add_feature(cfeature.LAND, facecolor='#F5F5F5')

    # Get data values
    lons = u.longitude.values
    lats = u.latitude.values
    u_data = u.values
    v_data = v.values

    print(f"Plotting hodographs...")

    # Plot Loop
    for i in range(0, lons.shape[0], GRID_SPACING):
        for j in range(0, lons.shape[1], GRID_SPACING):
            
            # 1. NaN check
            if np.isnan(u_data[:, i, j]).any(): continue

            # 2. Coordinate Check
            curr_lon = lons[i, j]
            curr_lat = lats[i, j]
            
            # Normalize lon to -180 to 180 for check
            check_lon = curr_lon - 360 if curr_lon > 180 else curr_lon
            
            # 3. Spatial Filter (This replaces the failed ds.sel step)
            if not (REGION[0] < check_lon < REGION[1] and REGION[2] < curr_lat < REGION[3]):
                continue

            # 4. Transform to Map Projection
            try:
                proj_pnt = ax.projection.transform_point(curr_lon, curr_lat, ccrs.PlateCarree())
            except:
                continue

            # 5. Draw Hodograph Inset
            try:
                box_size = 1.0 # Size in map units
                bounds = [proj_pnt[0] - box_size/2, proj_pnt[1] - box_size/2, box_size, box_size]
                sub_ax = ax.inset_axes(bounds, transform=ax.projection)
                
                h = Hodograph(sub_ax, component_range=80)
                h.add_grid(increment=20, color='gray', alpha=0.3, linewidth=0.5)
                h.plot(u_data[:, i, j], v_data[:, i, j], linewidth=1.2, color='blue')
                
                sub_ax.set_xticklabels([])
                sub_ax.set_yticklabels([])
                sub_ax.axis('off')
            except Exception:
                continue

    # Save
    plt.title(f"HREF Mean Hodographs | {date_str} {run}Z | f{fhr}", fontsize=18, weight='bold')
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    output_path = f"{OUTPUT_DIR}/href_mean_{date_str}_{run}z_f{fhr}.png"
    plt.savefig(output_path, bbox_inches='tight', dpi=100)
    print(f"Saved: {output_path}")

if __name__ == "__main__":
    date_str, run = get_latest_run_time()
    grib_file = download_href_mean(date_str, run, FORECAST_HOUR)
    
    try:
        plot_hodographs(grib_file, date_str, run, FORECAST_HOUR)
    finally:
        if os.path.exists(grib_file):
            os.remove(grib_file)
