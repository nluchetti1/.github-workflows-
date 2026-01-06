import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from metpy.plots import Hodograph
from metpy.units import units
import numpy as np
import datetime
from datetime import timedelta
import requests
import os
import sys
import warnings
import traceback

# Suppress warnings
warnings.filterwarnings("ignore")

# --- Configuration ---
# REGION: [West, East, South, North]
# Focused on Mid-Atlantic/Carolinas
REGION = [-86.0, -72.0, 32.0, 39.0]   

OUTPUT_DIR = "images"

# --- TUNING SETTINGS ---
# Grid Spacing: 23 (approx 69km). 
# Slightly increased spacing allows us to make the boxes bigger.
GRID_SPACING = 23             

# Box Size: 65000m (65km). 
# Increased from 55km. This fills the gap created by the spacing
# to maximize chart size without overlap.
BOX_SIZE = 65000              

# Levels for Hodographs
REQUESTED_LEVELS = [1000, 925, 850, 700, 500, 250] * units.hPa

# CAPE Color Levels
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
    return date.strftime('%Y%m%d'), run, date

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

def process_forecast_hour(date_obj, date_str, run, fhr):
    grib_file = download_href_mean(date_str, run, fhr)
    if not grib_file:
        return

    try:
        print(f"[f{fhr:02d}] Loading GRIB data...")
        
        # --- 1. Load Wind Data ---
        ds_u = xr.open_dataset(grib_file, engine='cfgrib', 
                               filter_by_keys={'typeOfLevel': 'isobaricInhPa', 'shortName': 'u'})
        ds_v = xr.open_dataset(grib_file, engine='cfgrib', 
                               filter_by_keys={'typeOfLevel': 'isobaricInhPa', 'shortName': 'v'})
        ds_wind = xr.merge([ds_u, ds_v])
        
        # --- 2. Load CAPE Data ---
        try:
            ds_cape = xr.open_dataset(grib_file, engine='cfgrib', 
                                      filter_by_keys={'shortName': 'cape', 'typeOfLevel': 'surface'})
        except Exception:
            try:
                ds_cape = xr.open_dataset(grib_file, engine='cfgrib', 
                                          filter_by_keys={'shortName': 'cape'})
            except Exception:
                print("CAPE not found. Plotting winds only.")
                ds_cape = None

        # --- 3. Filter Vertical Levels ---
        file_levels = ds_wind.isobaricInhPa.values
        levels_to_use = [l for l in REQUESTED_LEVELS.m if l in file_levels]
        
        if len(levels_to_use) < 3:
            print(f"Skipping f{fhr:02d}: Not enough vertical levels found.")
            return

        ds_wind = ds_wind.sel(isobaricInhPa=levels_to_use)
        u = ds_wind['u'].metpy.convert_units('kts')
        v = ds_wind['v'].metpy.convert_units('kts')

        # --- 4. Plotting Setup ---
        print(f"[f{fhr:02d}] Initializing Map...")
        fig = plt.figure(figsize=(18, 12))
        ax = fig.add_subplot(1, 1, 1, projection=ccrs.LambertConformal())
        ax.set_extent(REGION)
        
        ax.add_feature(cfeature.COASTLINE, linewidth=1.5, zorder=10)
        ax.add_feature(cfeature.BORDERS, linewidth=1.5, zorder=10)
        ax.add_feature(cfeature.STATES, linewidth=0.5, edgecolor='gray', zorder=10)

        # --- 5. Plot CAPE (Background) ---
        if ds_cape is not None:
            cape_data = ds_cape['cape']
            
            # Mask values < 250 as transparent
            cape_masked = np.where(cape_data.values < 250, np.nan, cape_data.values)
            
            cape_plot = ax.contourf(cape_data.longitude, cape_data.latitude, cape_masked, 
                                    levels=CAPE_LEVELS, cmap='nipy_spectral', 
                                    extend='max', alpha=0.6, transform=ccrs.PlateCarree())
            
            plt.colorbar(cape_plot, ax=ax, orientation='horizontal', pad=0.02, 
                         aspect=50, shrink=0.8, label='SBCAPE (J/kg)')

        # --- 6. Plot Hodographs (Foreground) ---
        print(f"[f{fhr:02d}] Generating Hodographs (Spacing: {GRID_SPACING})...")
        
        lons = u.longitude.values
        lats = u.latitude.values
        u_data = u.values
        v_data = v.values

        counter = 0
        
        for i in range(0, lons.shape[0], GRID_SPACING):
            for j in range(0, lons.shape[1], GRID_SPACING):
                
                if np.isnan(u_data[:, i, j]).any(): continue
                
                curr_lon = lons[i, j]
                curr_lat = lats[i, j]
                
                check_lon = curr_lon - 360 if curr_lon > 180 else curr_lon
                
                # Check Bounds
                if not (REGION[0] < check_lon < REGION[1] and REGION[2] < curr_lat < REGION[3]):
                    continue

                try:
                    proj_pnt = ax.projection.transform_point(curr_lon, curr_lat, ccrs.PlateCarree())
                except: continue

                # Inset Axis
                bounds = [proj_pnt[0] - BOX_SIZE/2, proj_pnt[1] - BOX_SIZE/2, BOX_SIZE, BOX_SIZE]
                sub_ax = ax.inset_axes(bounds, transform=ax.transData, zorder=20)
                
                h = Hodograph(sub_ax, component_range=80)
                h.add_grid(increment=20, color='black', alpha=0.2, linewidth=0.5)
                h.plot(u_data[:, i, j], v_data[:, i, j], linewidth=1.2, color='navy')
                
                sub_ax.set_xticklabels([])
                sub_ax.set_yticklabels([])
                sub_ax.axis('off')
                
                counter += 1

        print(f"[f{fhr:02d}] Plotted {counter} hodographs.")

        # --- 7. Title with Valid Time ---
        valid_time = date_obj + timedelta(hours=fhr)
        valid_str = valid_time.strftime("%a %H:%MZ") 
        
        plt.title(f"HREF Mean CAPE & Hodographs | Run: {date_str} {run}Z | Valid: {valid_str} (f{fhr:02d})", 
                  fontsize=16, weight='bold')
        
        os.makedirs(OUTPUT_DIR, exist_ok=True)
        out_path = f"{OUTPUT_DIR}/href_hodo_cape_{date_str}_{run}z_f{fhr:02d}.png"
        
        plt.savefig(out_path, bbox_inches='tight', dpi=80) 
        print(f"Saved: {out_path}")
        
        plt.close(fig)

    except Exception as e:
        print(f"Error processing f{fhr:02d}: {e}")
        traceback.print_exc()
    finally:
        if os.path.exists(grib_file):
            os.remove(grib_file)

if __name__ == "__main__":
    date_str, run, date_obj = get_latest_run_time()
    run_dt = datetime.datetime.strptime(f"{date_str} {run}", "%Y%m%d %H")
    
    print(f"Starting HREF Hodograph + CAPE generation for {date_str} {run}Z")
    print(f"Configuration: Region={REGION}")
    
    for fhr in range(1, 49):
        process_forecast_hour(run_dt, date_str, run, fhr)
