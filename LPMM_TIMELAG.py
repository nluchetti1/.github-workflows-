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
# This folder is created automatically in the GitHub workspace
output_folder = "href_data"
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

# --- 2. DOWNLOAD & UTILITY FUNCTIONS ---
def download_file(url, save_path, retries=5, delay=300):
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
        with pygrib.open(grib_path) as grb_file:
            for msg in grb_file:
                if msg.parameterName == "Total precipitation":
                    return msg.data()
    except Exception:
        return None, None, None

# --- 3. TIMING & RUN SELECTION ---
now_utc = datetime.now(pytz.UTC)
current_hour = now_utc.hour
base_url = "https://nomads.ncep.noaa.gov/pub/data/nccf/com/href/prod/"

if 13 <= current_hour <= 23:
    # 12Z Window
    date_now = datetime.now().strftime('%Y%m%d')
    date_prev = (datetime.now() - timedelta(days=1)).strftime('%Y%m%d')
    runs = [{"date": date_now, "hour": "12", "f_range": range(1, 24)},
            {"date": date_now, "hour": "00", "f_range": range(13, 36)},
            {"date": date_prev, "hour": "12", "f_range": range(25, 48)}]
    valid_title = f"Valid: 12Z {date_now} to 12Z {(datetime.now()+timedelta(days=1)).strftime('%Y%m%d')}"
else:
    # 00Z Window
    date_now = datetime.now().strftime('%Y%m%d')
    date_prev = (datetime.now() - timedelta(days=1)).strftime('%Y%m%d')
    runs = [{"date": date_now, "hour": "00", "f_range": range(1, 24)},
            {"date": date_prev, "hour": "12", "f_range": range(13, 36)},
            {"date": date_prev, "hour": "00", "f_range": range(25, 48)}]
    valid_title = f"Valid: 00Z {date_now} to 00Z {(datetime.now()+timedelta(days=1)).strftime('%Y%m%d')}"

# --- 4. PROCESSING ---
all_results = []
mm_to_inch = 0.0393701
lats, lons = None, None

for idx, run in enumerate(runs):
    hourly_data = []
    for f_hour in run['f_range']:
        f_str = f"{f_hour:02d}"
        url = f"{base_url}/href.{run['date']}/ensprod/href.t{run['hour']}z.conus.lpmm.f{f_str}.grib2"
        temp_path = os.path.join(output_folder, f"temp_{idx}_{f_str}.grib2")
        
        if download_file(url, temp_path):
            data, lats, lons = unpack_total_precipitation(temp_path)
            if data is not None:
                hourly_data.append(data * mm_to_inch)
            if os.path.exists(temp_path): os.remove(temp_path)
    
    if hourly_data:
        all_results.append({"data": np.sum(hourly_data, axis=0), 
                           "time": pd.to_datetime(run['date'] + ' ' + run['hour'] + 'Z')})

# --- 5. SAVING ALL PLOTS (DYNAMICALLY) ---
if len(all_results) >= 3:
    clevs = [0.00, 0.01, 0.1, 0.25, 0.50, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10., 12.0, 15.0, 18.0]
    cmap_data = [(1.0, 1.0, 1.0), (0.31, 0.81, 0.81), (0.0, 1.0, 1.0), (0.0, 0.87, 0.5), (0.0, 0.75, 0.0), (0.5, 0.87, 0.0), (1.0, 1.0, 0.0), (1.0, 0.62, 0.0), (1.0, 0.0, 0.0), (1.0, 0.12, 0.5), (0.94, 0.25, 1.0), (0.5, 0.12, 1.0), (0.25, 0.25, 1.0), (0.12, 0.12, 0.5), (0.12, 0.12, 0.12), (0.5, 0.5, 0.5), (0.87, 0.87, 0.87), (0.93, 0.83, 0.73), (0.85, 0.65, 0.47), (0.62, 0.42, 0.23), (0.4, 0.2, 0.0)]
    cmap = mcolors.ListedColormap(cmap_data, 'precipitation')
    norm = mcolors.BoundaryNorm(clevs, cmap.N)

    # Naming Convention 1: 3-Panel Run Comparison
    fig, axes = plt.subplots(1, 3, figsize=(18, 8), subplot_kw={'projection': ccrs.PlateCarree()})
    for i, res in enumerate(all_results):
        cs = axes[i].contourf(lons, lats, res['data'], clevs, cmap=cmap, norm=norm, alpha=0.5)
        axes[i].coastlines()
        axes[i].add_feature(cfeature.STATES, linewidth=1.0)
        axes[i].set_extent([-84.8, -74, 31, 39])
        axes[i].set_title(f'{res["time"].strftime("%Y-%m-%d %H:%M Z")}')
    
    fig.suptitle(f'24hr HREF LPMM Run Comparison\n{valid_title}', fontsize=16, fontweight='bold')
    # SAVE POINT 1
    fig.savefig(os.path.join(output_folder, 'HREF_LPMM_RUN_COMPARE.png'), dpi=300)

    # Naming Convention 2: 4-Panel Threshold Paintball Plots
    fig2, ax2 = plt.subplots(2, 2, figsize=(15, 12), subplot_kw={'projection': ccrs.PlateCarree()})
    thresholds = [3, 6, 9, 12]
    colors = ['#00008B', '#4169E1', '#87CEFA']
    
    for i, thresh in enumerate(thresholds):
        row, col = divmod(i, 2)
        for j, res in enumerate(all_results):
            m_data = np.ma.masked_less(res['data'], thresh)
            ax2[row, col].contourf(lons, lats, m_data, cmap=mcolors.ListedColormap([colors[j]]), levels=[thresh, 99], alpha=0.6)
        
        ax2[row, col].coastlines()
        ax2[row, col].set_extent([-84.8, -74, 31, 39])
        ax2[row, col].set_title(f'> {thresh} inches', fontweight='bold')
    
    fig2.suptitle(f'24hr HREF LPMM Threshold Compare\n{valid_title}', fontsize=16, fontweight='bold')
    # SAVE POINT 2
    fig2.savefig(os.path.join(output_folder, 'HREF_LPMM_THRESHOLD_COMPARE.png'), dpi=300, bbox_inches='tight')

print("All plots generated and saved in the 'href_data' folder.")