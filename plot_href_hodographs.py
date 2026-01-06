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
REGION = [-98, -74, 24, 38]   # Southeast US
GRID_SPACING = 10             # Plot density
FORECAST_HOUR = '01'          # Changed from '00' to '01' (Mean usually starts at f01)
OUTPUT_DIR = "images"

# Levels to look for. 
REQUESTED_LEVELS = [1000, 925, 850, 700, 500, 250] * units.hPa

def get_latest_run_time():
    """Determines the latest available run."""
    now = datetime.datetime.utcnow()
    # Logic: 12Z is usually ready by 15:00Z; 00Z by 03:00Z
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
    """Downloads the HREF Ensemble Mean file from the ensprod directory."""
    
    # URL Structure
    # .../href.20260106/ensprod/href.t00z.conus.mean.f01.grib2
    base_url = f"https://nomads.ncep.noaa.gov/pub/data/nccf/com/href/prod/href.{date_str}/ensprod"
    filename = f"href.t{run}z.conus.mean.f{fhr}.grib2"
    url = f"{base_url}/{filename}"
    
    print(f"Downloading: {url}")
    
    # NOMADS requires a User-Agent or it may block the request (403 or 404)
    headers = {
        'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36'
    }
    
    try:
        with requests.get(url, stream=True, timeout=60, headers=headers) as r:
            if r.status_code == 404:
                print(f"Error 404: File not found. The mean product might not exist for f{fhr}.")
                print("Try checking f01, f02, etc. instead of f00.")
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
    print("Reading GRIB data...")
    
    try:
        ds = xr.open_dataset(filename, engine='cfgrib', 
                             filter_by_keys={'typeOfLevel': 'isobaricInhPa'})
    except Exception as e:
        print(f"Failed to read GRIB: {e}")
        sys.exit(1)

    ds = ds.sel(latitude=slice(REGION[3], REGION[2]), 
                longitude=slice(360+REGION[0], 360+REGION[1]))

    file_levels = ds.isobaricInhPa.values
    levels_to_use = [l for l in REQUESTED_LEVELS.m if l in file_levels]
    
    if len(levels_to_use) < 3:
        print("Error: Not enough vertical levels in this file for a hodograph.")
        sys.exit(1)
        
    ds = ds.sel(isobaricInhPa=levels_to_use)

    u = ds['u'].metpy.convert_units('kts')
    v = ds['v'].metpy.convert_units('kts')

    fig = plt.figure(figsize=(18, 12))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.LambertConformal())
    ax.set_extent(REGION)
    
    ax.add_feature(cfeature.COASTLINE, linewidth=1.5)
    ax.add_feature(cfeature.BORDERS, linewidth=1.5)
    ax.add_feature(cfeature.STATES, linewidth=0.5, edgecolor='gray')
    ax.add_feature(cfeature.LAND, facecolor='#F5F5F5')

    lons = u.longitude.values
    lats = u.latitude.values
    u_data = u.values
    v_data = v.values

    print(f"Plotting using levels: {levels_to_use} hPa...")

    for i in range(0, lons.shape[0], GRID_SPACING):
        for j in range(0, lons.shape[1], GRID_SPACING):
            
            if np.isnan(u_data[:, i, j]).any(): continue

            curr_lon = lons[i, j]
            curr_lat = lats[i, j]
            
            check_lon = curr_lon - 360 if curr_lon > 180 else curr_lon
            if not (REGION[0] < check_lon < REGION[1] and REGION[2] < curr_lat < REGION[3]):
                continue

            proj_pnt = ax.projection.transform_point(curr_lon, curr_lat, ccrs.PlateCarree())
            
            # Draw Hodograph
            try:
                # Inset axes
                box_size = 1.0
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
