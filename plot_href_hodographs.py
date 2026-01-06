import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
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
# ZOOMED DOMAIN: Tight focus on SC/NC
REGION = [-83.5, -75.5, 32.5, 37.5]   
OUTPUT_DIR = "images"

# --- TUNING SETTINGS ---
# Grid Spacing: 25 grid points (approx 75km spacing)
GRID_SPACING = 25             

# Box Size: 100km (100,000m). 
# Large enough to see details clearly in this zoomed view.
BOX_SIZE = 100000              

# Preferred Levels (Surface -> Aloft)
REQUESTED_LEVELS = [1000, 925, 850, 700, 500, 250]

# CAPE Color Levels (Starts at 250)
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

def get_segment_color(pressure_start, pressure_end):
    """
    Determines color based on the pressure level of the segment.
    Uses Standard Atmosphere approximation for heights.
    """
    avg_p = (pressure_start + pressure_end) / 2.0
    
    if avg_p >= 900:      # 0-1 km approx
        return 'magenta'
    elif 700 <= avg_p < 900: # 1-3 km approx
        return 'red'
    elif 500 <= avg_p < 700: # 3-6 km approx
        return 'green'
    else:                 # > 6 km
        return 'gold'

def plot_colored_hodograph(ax, u, v, levels):
    """
    Plots a multi-colored hodograph based on actual pressure levels found.
    """
    # Loop through segments
    for k in range(len(u) - 1):
        p_start = levels[k]
        p_end = levels[k+1]
        
        color = get_segment_color(p_start, p_end)
        
        # INCREASED LINE WIDTH HERE (was 1.5, now 3.0)
        ax.plot([u[k], u[k+1]], [v[k], v[k+1]], color=color, linewidth=3.0)

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
        # Robust loading attempt
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

        # --- 3. ROBUST Level Filtering ---
        file_levels = ds_wind.isobaricInhPa.values
        available_levels = sorted([l for l in REQUESTED_LEVELS if l in file_levels], reverse=True)
        
        if len(available_levels) < 3:
            print(f"Skipping f{fhr:02d}: Not enough vertical levels.")
            return
        
        print(f"       Using levels: {available_levels}")

        ds_wind = ds_wind.sel(isobaricInhPa=available_levels)
        u = ds_wind['u'].metpy.convert_units('kts')
        v = ds_wind['v'].metpy.convert_units('kts')
        
        level_values = available_levels

        # --- 4. Plotting Setup ---
        print(f"[f{fhr:02d}] Initializing Map...")
        fig = plt.figure(figsize=(18, 12))
        ax = fig.add_subplot(1, 1, 1, projection=ccrs.LambertConformal())
        ax.set_extent(REGION)
        
        ax.add_feature(cfeature.COASTLINE, linewidth=2.0, zorder=10) # Thicker coast
        ax.add_feature(cfeature.BORDERS, linewidth=2.0, zorder=10)
        ax.add_feature(cfeature.STATES, linewidth=1.5, edgecolor='black', zorder=10)

        # --- 5. Plot CAPE (Background) ---
        if ds_cape is not None:
            cape_data = ds_cape['cape']
            
            # STRICT MASKING: explicitly mask the array where values < 250
            # This ensures they are not plotted at all (transparent)
            cape_vals = cape_data.values
            cape_masked = np.ma.masked_where(cape_vals < 250, cape_vals)
            
            cape_plot = ax.contourf(cape_data.longitude, cape_data.latitude, cape_masked, 
                                    levels=CAPE_LEVELS, cmap='nipy_spectral', 
                                    extend='max', alpha=0.5, transform=ccrs.PlateCarree())
            
            # Make sure we don't accidentally plot a background color for the masked values
            ax.set_facecolor('white')

            plt.colorbar(cape_plot, ax=ax, orientation='horizontal', pad=0.02, 
                         aspect=50, shrink=0.8, label='SBCAPE (J/kg)')

        # --- 6. ADD LEGEND ---
        legend_elements = [
            mlines.Line2D([], [], color='magenta', lw=3, label='0-1 km (>900mb)'),
            mlines.Line2D([], [], color='red', lw=3, label='1-3 km (900-700mb)'),
            mlines.Line2D([], [], color='green', lw=3, label='3-6 km (700-500mb)'),
            mlines.Line2D([], [], color='gold', lw=3, label='6-9 km (<500mb)')
        ]
        ax.legend(handles=legend_elements, loc='upper left', title="Hodograph Layers", 
                  framealpha=0.9, fontsize=12, title_fontsize=13).set_zorder(100)

        # --- 7. Plot Hodographs (Foreground) ---
        print(f"[f{fhr:02d}] Generating Colored Hodographs...")
        
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
                
                if not (REGION[0] < check_lon < REGION[1] and REGION[2] < curr_lat < REGION[3]):
                    continue

                try:
                    proj_pnt = ax.projection.transform_point(curr_lon, curr_lat, ccrs.PlateCarree())
                except: continue

                bounds = [proj_pnt[0] - BOX_SIZE/2, proj_pnt[1] - BOX_SIZE/2, BOX_SIZE, BOX_SIZE]
                sub_ax = ax.inset_axes(bounds, transform=ax.transData, zorder=20)
                
                h = Hodograph(sub_ax, component_range=80)
                # Made grid slightly darker (alpha 0.3) so it's easier to see against white land
                h.add_grid(increment=20, color='black', alpha=0.3, linewidth=0.5)
                
                u_profile = u_data[:, i, j]
                v_profile = v_data[:, i, j]
                
                plot_colored_hodograph(h.ax, u_profile, v_profile, level_values)
                
                sub_ax.set_xticklabels([])
                sub_ax.set_yticklabels([])
                sub_ax.axis('off')
                
                counter += 1

        print(f"[f{fhr:02d}] Plotted {counter} hodographs.")

        # --- 8. Title ---
        valid_time = date_obj + timedelta(hours=fhr)
        valid_str = valid_time.strftime("%a %H:%MZ") 
        
        plt.title(f"HREF Mean CAPE & Hodographs | Run: {date_str} {run}Z | Valid: {valid_str} (f{fhr:02d})", 
                  fontsize=18, weight='bold')
        
        os.makedirs(OUTPUT_DIR, exist_ok=True)
        out_path = f"{OUTPUT_DIR}/href_hodo_cape_{date_str}_{run}z_f{fhr:02d}.png"
        
        plt.savefig(out_path, bbox_inches='tight', dpi=90) 
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
