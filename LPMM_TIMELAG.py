import matplotlib
matplotlib.use('Agg')  # Required for GitHub Actions
import os
import time
import numpy as np
import requests
import pygrib
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import pandas as pd
from datetime import datetime, timedelta
from metpy.plots import USCOUNTIES
from matplotlib.lines import Line2D
import pytz

# --- 1. DYNAMIC FOLDER SETUP ---
output_folder = "href_data"
if not os.path.exists(output_folder):
    os.makedirs(output_folder) # FIXED INDENTATION HERE

# --- 2. DOWNLOAD & UTILITY FUNCTIONS ---
def download_file(url, save_path, retries=3, delay=60):
    headers = {"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) Chrome/87.0.4280.88 Safari/537.36"}
    for attempt in range(retries):
        try:
            response = requests.get(url, headers=headers, timeout=30)
            if response.status_code == 200:
                with open(save_path, 'wb') as file:
                    file.write(response.content)
                return True
        except Exception as e:
            print(f"Attempt {attempt + 1} failed: {e}")
        time.sleep(delay)
    return False

def unpack_total_precipitation(grib_path):
    try:
        if not os.path.exists(grib_path):
            return None, None, None
        with pygrib.open(grib_path) as grb_file:
            for msg in grb_file:
                if msg.parameterName == "Total precipitation":
                    return msg.data()
    except Exception:
        pass
    return None, None, None

# --- 3. TIMING & RUN SELECTION ---
now_utc = datetime.now(pytz.UTC)
current_hour = now_utc.hour
base_url = "https://nomads.ncep.noaa.gov/pub/data/nccf/com/href/prod/"

# Logic to pick the correct 00Z or 12Z dates
date_now = now_utc.strftime('%Y%m%d')
date_prev = (now_utc - timedelta(days=1)).strftime('%Y%m%d')

if 13 <= current_hour <= 23:
    runs = [{"date": date_now, "hour": "12", "f_range": range(1, 24)},
            {"date": date_now, "hour": "00", "f_range": range(13, 36)},
            {"date": date_prev, "hour": "12", "f_range": range(25, 48)}]
    valid_title = f"Valid: 12Z {date_now} to 12Z {(now_utc+timedelta(days=1)).strftime('%Y%m%d')}"
else:
    runs = [{"date": date_now, "hour": "00", "f_range": range(1, 24)},
            {"date": date_prev, "hour": "12", "f_range": range(13, 36)},
            {"date": date_prev, "hour": "00", "f_range": range(25, 48)}]
    valid_title = f"Valid: 00Z {date_now} to 00Z {(now_utc+timedelta(days=1)).strftime('%Y%m%d')}"

# --- 4. DATA COLLECTION ---
all_results = []
mm_to_inch = 0.0393701
final_lats, final_lons = None, None

for idx, run in enumerate(runs):
    hourly_data = []
    print(f"Starting Run: {run['date']} {run['hour']}Z")
    for f_hour in run['f_range']:
        f_str = f"{f_hour:02d}"
        url = f"{base_url}/href.{run['date']}/ensprod/href.t{run['hour']}z.conus.lpmm.f{f_str}.grib2"
        temp_path = os.path.join(output_folder, f"temp_{idx}_{f_str}.grib2")
        
        if download_file(url, temp_path):
            data, lats, lons = unpack_total_precipitation(temp_path)
            if data is not None:
                hourly_data.append(data * mm_to_inch)
                final_lats, final_lons = lats, lons
            if os.path.exists(temp_path): os.remove(temp_path)
    
    if hourly_data:
        all_results.append({"data": np.sum(hourly_data, axis=0), 
                           "time": pd.to_datetime(run['date'] + ' ' + run['hour'] + 'Z')})

# --- 5. PLOTTING ---
if len(all_results) >= 1:
    clevs = [0.0, 0.01, 0.1, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10., 12., 15., 18.]
    cmap_data = [(1,1,1), (0.31,0.81,0.81), (0,1,1), (0,0.87,0.5), (0,0.75,0), (0.5,0.87,0), (1,1,0), (1,0.62,0), (1,0,0), (1,0.12,0.5), (0.94,0.25,1), (0.5,0.12,1), (0.25,0.25,1), (0.12,0.12,0.5), (0.12,0.12,0.12), (0.5,0.5,0.5), (0.87,0.87,0.87), (0.93,0.83,0.73), (0.85,0.65,0.47), (0.62,0.42,0.23), (0.4,0.2,0)]
    cmap = mcolors.ListedColormap(cmap_data, 'precip')
    norm = mcolors.BoundaryNorm(clevs, cmap.N)

    # Comparison Plot
    fig, axes = plt.subplots(1, len(all_results), figsize=(6*len(all_results), 8), subplot_kw={'projection': ccrs.PlateCarree()})
    if len(all_results) == 1: axes = [axes]
    
    for i, res in enumerate(all_results):
        cs = axes[i].contourf(final_lons, final_lats, res['data'], clevs, cmap=cmap, norm=norm, alpha=0.5)
        axes[i].coastlines()
        axes[i].add_feature(cfeature.STATES, linewidth=1.0)
        axes[i].set_extent([-84.8, -74, 31, 39])
        axes[i].set_title(f'{res["time"].strftime("%Y-%m-%d %H:%M Z")}')
    
    fig.suptitle(f'24hr HREF LPMM Run Comparison\n{valid_title}', fontsize=16, fontweight='bold')
    fig.savefig(os.path.join(output_folder, 'HREF_LPMM_RUN_COMPARE.png'), dpi=300)

    # Threshold Plot
    fig2, ax2 = plt.subplots(2, 2, figsize=(15, 12), subplot_kw={'projection': ccrs.PlateCarree()})
    thresholds = [3, 6, 9, 12]
    colors = ['#00008B', '#4169E1', '#87CEFA']
    for i, thresh in enumerate(thresholds):
        row, col = divmod(i, 2)
        for j, res in enumerate(all_results):
            m_data = np.ma.masked_less(res['data'], thresh)
            ax2[row, col].contourf(final_lons, final_lats, m_data, cmap=mcolors.ListedColormap([colors[j]]), levels=[thresh, 99], alpha=0.6)
        ax2[row, col].coastlines(); ax2[row, col].set_extent([-84.8, -74, 31, 39])
        ax2[row, col].set_title(f'> {thresh} inches', fontweight='bold')
    
    fig2.suptitle(f'24hr HREF LPMM Threshold Compare\n{valid_title}', fontsize=16, fontweight='bold')
    fig2.savefig(os.path.join(output_folder, 'HREF_LPMM_THRESHOLD_COMPARE.png'), dpi=300, bbox_inches='tight')

print("Process Complete.")
