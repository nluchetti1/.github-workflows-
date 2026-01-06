import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from metpy.plots import Hodograph
from metpy.units import units
import numpy as np
import datetime
import requests
import os
import sys
import warnings

# Suppress warnings
warnings.filterwarnings("ignore")

# --- Configuration ---
REGION = [-98, -74, 24, 38]   # Southeast US
GRID_SPACING = 10             # Hodograph density
OUTPUT_DIR = "images"

# Levels for Hodographs
REQUESTED_LEVELS = [1000, 925, 850, 700, 500, 250] * units.hPa

# CAPE Color Levels (J/kg)
CAPE_LEVELS = np.arange(250, 5001, 250) 

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
    """Downloads the HREF Mean file for a specific forecast hour."""
    base_url = f"https://nomads.ncep.noaa.gov/pub/data/nccf/com/href/prod/href.{date_str}/ensprod"
    filename = f"href.t{run}z.conus.mean.f{fhr:02d}.grib2"
    url = f"{base_url}/{filename}"
    
    print(f"\n[f{fhr:02d}] Downloading: {url}")
    headers = {'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64)'}
    
    try:
        with requests.get(url, stream=True, timeout=60, headers=headers) as r:
            if r.status_code == 404:
                print(f"File not found (404). Forecast f{fhr:02d} may not be ready yet.")
                return None
            r.raise_for_status()
            with open(filename, 'wb') as f:
                for chunk in r.iter_content(chunk_size=8192):
                    f.write(chunk)
        return filename
    except Exception as e:
        print(f"Download failed: {e}")
        return None

def process_forecast_hour(date_str, run, fhr):
    grib_file = download_href_mean(date_str, run, fhr)
    if not grib_file:
        return

    try:
        # --- 1. Load Wind Data (Isobaric) ---
        # Load U and V separately to avoid conflicts
        ds_u = xr.open_dataset(grib_file, engine='cfgrib', 
                               filter_by_keys={'typeOfLevel': 'isobaricInhPa', 'shortName': 'u'})
        ds_v = xr.open_dataset(grib_file, engine='cfgrib', 
                               filter_by_keys={'typeOfLevel': 'isobaricInhPa', 'shortName': 'v'})
        ds_wind = xr.merge([ds_u, ds_v])
        
        # --- 2. Load CAPE Data (Surface) ---
        # Improved filter: Explicitly ask for Surface CAPE to avoid ambiguity
        try:
            ds_cape = xr.open_dataset(grib_file, engine='cfgrib', 
                                      filter_by_keys={'shortName': 'cape', 'typeOfLevel': 'surface'})
        except Exception:
            try:
                # Fallback: Try just shortName if 'surface' fails
                ds_cape = xr.open_dataset(grib_file, engine='cfgrib', 
                                          filter_by_keys={'shortName': 'cape'})
            except Exception:
                print("CAPE not found in this file. Plotting winds only.")
                ds_cape = None

        # --- 3. Filter Vertical Levels (Winds) ---
        file_levels = ds_wind.isobaricInhPa.values
        levels_to_use = [l for l in REQUESTED_LEVELS.m if l in file_levels]
        if len(levels_to_use) < 3:
            print(f"Skipping f{fhr:02d}: Not enough vertical levels.")
            return

        ds_wind = ds_wind.sel(isobaricInhPa=levels_to_use)
        u = ds_wind['u'].metpy.convert_units('kts')
        v = ds_wind['v'].metpy.convert_units('kts')

        # --- 4. Plotting Setup ---
        fig = plt.figure(figsize=(18, 12))
        ax = fig.add_subplot(1, 1, 1, projection=ccrs.LambertConformal())
        ax.set_extent(REGION)
        
        ax.add_feature(cfeature.COASTLINE, linewidth=1.5, zorder=10)
        ax.add_feature(cfeature.BORDERS, linewidth=1.5, zorder=10)
        ax.add_feature(cfeature.STATES, linewidth=0.5, edgecolor='gray', zorder=10)

        # --- 5. Plot CAPE (Background) ---
        if ds_cape is not None:
            cape = ds_cape['cape']
            cape_plot = ax.contourf(cape.longitude, cape.latitude, cape.values, 
                                    levels=CAPE_LEVELS, cmap='nipy_spectral', 
                                    extend='max', alpha=0.6, transform=ccrs.PlateCarree())
            
            plt.colorbar(cape_plot, ax=ax, orientation='horizontal', pad=0.02, 
                         aspect=50, shrink=0.8, label='SBCAPE (J/kg)')

        # --- 6. Plot Hodographs (Foreground) ---
        lons = u.longitude.values
        lats = u.latitude.values
        u_data = u.values
        v_data = v.values

        # Set Hodograph Size (in Meters for Lambert Conformal)
        # 0.9 was too small (meters). 75km (75000m) is roughly 0.7 degrees.
        box_size = 75000 

        for i in range(0, lons.shape[0], GRID_SPACING):
            for j in range(0, lons.shape[1], GRID_SPACING):
                
                if np.isnan(u_data[:, i, j]).any(): continue
                
                curr_lon = lons[i, j]
                curr_lat = lats[i, j]
                
                # Check Bounds (Longitude correction)
                check_lon = curr_lon - 360 if curr_lon > 180 else curr_lon
                
                if not (REGION[0] < check_lon < REGION[1] and REGION[2] < curr_lat < REGION[3]):
                    continue

                # Transform Lat/Lon to Project Coordinates (Meters)
                try:
                    proj_pnt = ax.projection.transform_point(curr_lon, curr_lat, ccrs.PlateCarree())
                except: continue

                # Define Bounds in PROJECTED coordinates (Meters)
                bounds = [proj_pnt[0] - box_size/2, proj_pnt[1] - box_size/2, box_size, box_size]
                
                # FIX 1: Use ax.transData instead of ax.projection
                # transData handles the data coordinates (which are in meters here)
                sub_ax = ax.inset_axes(bounds, transform=ax.transData, zorder=20)
                
                h = Hodograph(sub_ax, component_range=80)
                h.add_grid(increment=20, color='black', alpha=0.2, linewidth=0.5)
                h.plot(u_data[:, i, j], v_data[:, i, j], linewidth=1.0, color='navy')
                
                sub_ax.set_xticklabels([])
                sub_ax.set_yticklabels([])
                sub_ax.axis('off')

        # Save and Close
        title_time = f"f{fhr:02d}"
        plt.title(f"HREF Mean CAPE & Hodographs | Run: {date_str} {run}Z | Forecast: +{title_time}", fontsize=16, weight='bold')
        
        os.makedirs(OUTPUT_DIR, exist_ok=True)
        out_path = f"{OUTPUT_DIR}/href_hodo_cape_{date_str}_{run}z_f{fhr:02d}.png"
        plt.savefig(out_path, bbox_inches='tight', dpi=100)
        print(f"Saved: {out_path}")
        
        plt.close(fig)

    except Exception as e:
        # Print full traceback if needed, but simple error for now
        print(f"Error processing f{fhr:02d}: {e}")
        import traceback
        traceback.print_exc() # Useful to see exact line of failure
    finally:
        if os.path.exists(grib_file):
            os.remove(grib_file)

if __name__ == "__main__":
    date_str, run = get_latest_run_time()
    print(f"Starting HREF Hodograph + CAPE generation for {date_str} {run}Z")
    
    for fhr in range(1, 49):
        process_forecast_hour(date_str, run, fhr)
