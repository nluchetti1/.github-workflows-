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
import traceback

# Suppress warnings for cleaner logs
warnings.filterwarnings("ignore")

# --- Configuration ---
REGION = [-98, -74, 24, 38]   # Southeast US [West, East, South, North]

# CRITICAL SETTINGS FOR PERFORMANCE
GRID_SPACING = 50             # Skip every 50 points (approx every 150km). 
                              # Lower = more density but MUCH slower.
BOX_SIZE = 75000              # Hodograph width in meters (75km). 
                              # Needs to be smaller than GRID_SPACING * 3000

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
        print(f"[f{fhr:02d}] Loading GRIB data...")
        
        # --- 1. Load Wind Data ---
        # Load U and V separately
        ds_u = xr.open_dataset(grib_file, engine='cfgrib', 
                               filter_by_keys={'typeOfLevel': 'isobaricInhPa', 'shortName': 'u'})
        ds_v = xr.open_dataset(grib_file, engine='cfgrib', 
                               filter_by_keys={'typeOfLevel': 'isobaricInhPa', 'shortName': 'v'})
        ds_wind = xr.merge([ds_u, ds_v])
        
        # --- 2. Load CAPE Data ---
        # Try specific surface level first, fallback to generic name
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
            cape = ds_cape['cape']
            # Plot contours
            cape_plot = ax.contourf(cape.longitude, cape.latitude, cape.values, 
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
        
        # Loop through grid
        for i in range(0, lons.shape[0], GRID_SPACING):
            for j in range(0, lons.shape[1], GRID_SPACING):
                
                # Basic NaN check
                if np.isnan(u_data[:, i, j]).any(): continue
                
                curr_lon = lons[i, j]
                curr_lat = lats[i, j]
                
                # Longitude correction (0-360 to -180/180)
                check_lon = curr_lon - 360 if curr_lon > 180 else curr_lon
                
                # Check if point is inside our region box
                if not (REGION[0] < check_lon < REGION[1] and REGION[2] < curr_lat < REGION[3]):
                    continue

                # Transform Lat/Lon to Projected Meters
                try:
                    proj_pnt = ax.projection.transform_point(curr_lon, curr_lat, ccrs.PlateCarree())
                except: continue

                # Define Bounds in PROJECTED coordinates (Meters)
                # This places the box centered on the grid point
                bounds = [proj_pnt[0] - BOX_SIZE/2, proj_pnt[1] - BOX_SIZE/2, BOX_SIZE, BOX_SIZE]
                
                # Create the inset axis
                # FIX: Use ax.transData to handle the map coordinates correctly
                sub_ax = ax.inset_axes(bounds, transform=ax.transData, zorder=20)
                
                # Create Hodograph
                h = Hodograph(sub_ax, component_range=80)
                h.add_grid(increment=20, color='black', alpha=0.2, linewidth=0.5)
                h.plot(u_data[:, i, j], v_data[:, i, j], linewidth=1.2, color='navy')
                
                # Hide axes labels
                sub_ax.set_xticklabels([])
                sub_ax.set_yticklabels([])
                sub_ax.axis('off')
                
                counter += 1

        print(f"[f{fhr:02d}] Plotted {counter} hodographs. Saving image...")

        # Save and Close
        title_time = f"f{fhr:02d}"
        plt.title(f"HREF Mean CAPE & Hodographs | Run: {date_str} {run}Z | Forecast: +{title_time}", fontsize=16, weight='bold')
        
        os.makedirs(OUTPUT_DIR, exist_ok=True)
        out_path = f"{OUTPUT_DIR}/href_hodo_cape_{date_str}_{run}z_f{fhr:02d}.png"
        
        # Use lower DPI (80) first to ensure speed, then increase if stable
        plt.savefig(out_path, bbox_inches='tight', dpi=80) 
        print(f"Success: {out_path}")
        
        plt.close(fig) # Crucial to free memory

    except Exception as e:
        print(f"Error processing f{fhr:02d}: {e}")
        traceback.print_exc()
    finally:
        # Cleanup GRIB file
        if os.path.exists(grib_file):
            os.remove(grib_file)

if __name__ == "__main__":
    date_str, run = get_latest_run_time()
    print(f"Starting HREF Hodograph + CAPE generation for {date_str} {run}Z")
    print(f"Configuration: Grid Spacing={GRID_SPACING}, Box Size={BOX_SIZE}m")
    
    # Process f01 to f48
    for fhr in range(1, 49):
        process_forecast_hour(date_str, run, fhr)
