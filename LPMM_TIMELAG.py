import matplotlib
matplotlib.use('Agg')
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
    os.makedirs(output_folder)

# --- 2. DOWNLOAD & UTILITY FUNCTIONS ---
def download_file(url, save_path, retries=2, delay=10):
    headers = {"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) Chrome/87.0.4280.88 Safari/537.36"}
    try:
        response = requests.get(url, headers=headers, timeout=30)
        if response.status_code == 200:
            with open(save_path, 'wb') as file:
                file.write(response.content)
            print(f"Success: {os.path.basename(url)}")
            return True
        else:
            print(f"Missing (404): {os.path.basename(url)}")
    except Exception as e:
        print(f"Error downloading {url}: {e}")
    return False

def unpack_total_precipitation(grib_path):
    try:
        with pygrib.open(grib_path) as grb_file:
            msgs = grb_file.select(parameterCategory=1, parameterNumber=8)
            if msgs:
                msg = msgs[0]
                return msg.data()
            msgs = grb_file.select(name='Total precipitation')
            if msgs:
                return msgs[0].data()
    except Exception as e:
        print(f"Unpack error: {e}")
    return None, None, None

# --- 3. TIMING & RUN SELECTION ---
now_utc = datetime.now(pytz.UTC)
current_hour = now_utc.hour
base_url = "https://nomads.ncep.noaa.gov/pub/data/nccf/com/href/prod"
date_now = now_utc.strftime('%Y%m%d')
date_prev = (now_utc - timedelta(days=1)).strftime('%Y%m%d')

if 13 <= current_hour <= 23:
    runs = [{"date": date_now, "hour": "12", "f_range": range(1, 25)},
            {"date": date_now, "hour": "00", "f_range": range(13, 37)},
            {"date": date_prev, "hour": "12", "f_range": range(25, 49)}]
    valid_range = f"Valid: 12Z {date_now} to 12Z {(now_utc+timedelta(days=1)).strftime('%Y%m%d')}"
else:
    runs = [{"date": date_now, "hour": "00", "f_range": range(1, 25)},
            {"date": date_prev, "hour": "12", "f_range": range(13, 37)},
            {"date": date_prev, "hour": "00", "f_range": range(25, 49)}]
    valid_range = f"Valid: 00Z {date_now} to 00Z {(now_utc+timedelta(days=1)).strftime('%Y%m%d')}"

# --- 4. DATA COLLECTION ---
all_results = []
mm_to_inch = 0.0393701
final_lats, final_lons = None, None

for idx, run in enumerate(runs):
    hourly_data = []
    print(f"--- Processing Run: {run['date']} {run['hour']}Z ---")
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

    # 1. Comparison Plot (3-panel)
    # Reducing figsize height and adjusting subplots_adjust to kill top whitespace
    fig, axes = plt.subplots(1, 3, figsize=(18, 6), subplot_kw={'projection': ccrs.PlateCarree()})
    
    # Manual adjustment gives better control over the 'big gap' than constrained_layout sometimes does
    fig.subplots_adjust(left=0.05, right=0.95, bottom=0.15, top=0.82, wspace=0.05)
    
    cs = None
    for i in range(3):
        if i < len(all_results):
            res = all_results[i]
            cs = axes[i].contourf(final_lons, final_lats, res['data'], clevs, cmap=cmap, norm=norm, alpha=0.5)
            axes[i].set_title(f'{res["time"].strftime("%Y-%m-%d %H:%M Z")}\n24hr HREF LPMM [in]', fontsize=9, fontweight='bold', pad=10)
        
        axes[i].coastlines(resolution='10m')
        axes[i].add_feature(cfeature.STATES, linewidth=0.8, edgecolor='black')
        axes[i].add_feature(USCOUNTIES.with_scale('500k'), edgecolor='gray', linewidth=0.4)
        axes[i].set_extent([-122, -114, 32, 37])

    # Colorbar positioning
    cbar_ax = fig.add_axes([0.15, 0.08, 0.7, 0.03])
    cbar = fig.colorbar(cs, cax=cbar_ax, orientation='horizontal', ticks=clevs)
    cbar.set_label('Precipitation (inches)', fontsize=12, fontweight='bold')
    
    fig.suptitle(f'24hr HREF LPMM [in] dprog/dt\n{valid_range}', fontsize=16, fontweight='bold', y=0.95)
    plt.savefig(os.path.join(output_folder, 'HREF_LPMM_RUN_COMPARE.png'), dpi=300)

    # 2. Threshold Plot (4-panel)
    # Adjusting figsize to (12, 12) to make it more square and less stretched horizontally
    fig2, ax2 = plt.subplots(2, 2, figsize=(12, 12), subplot_kw={'projection': ccrs.PlateCarree()})
    
    thresholds = [3, 6, 9, 12]
    blue_shades = ['#00008B', '#4169E1', '#87CEFA']
    legend_elements = [Line2D([0], [0], marker='o', color='w', label=res['time'].strftime("%Y-%m-%d %H:%M Z"),
                              markerfacecolor=blue_shades[idx], markersize=8) for idx, res in enumerate(all_results)]

    for i, thresh in enumerate(thresholds):
        row, col = divmod(i, 2)
        for j, res in enumerate(all_results):
            m_data = np.ma.masked_less(res['data'], thresh)
            ax2[row, col].contourf(final_lons, final_lats, m_data, cmap=mcolors.ListedColormap([blue_shades[j]]), levels=[thresh, 99], alpha=0.6)
        
        ax2[row, col].coastlines(resolution='10m')
        ax2[row, col].add_feature(cfeature.STATES, linewidth=0.8, edgecolor='black')
        ax2[row, col].set_extent([-122, -114, 32, 37])
        ax2[row, col].set_title(f'> {thresh} inches', fontsize=12, fontweight='bold')
        ax2[row, col].legend(handles=legend_elements, loc='lower right', title='HREF Run', fontsize=8)
    
    # Using subplots_adjust here instead of tight_layout for finer control over the suptitle gap
    fig2.subplots_adjust(left=0.05, right=0.95, bottom=0.05, top=0.90, wspace=0.1, hspace=0.15)
    fig2.suptitle(f'24hr HREF LPMM Threshold Compare\n{valid_range}', fontsize=16, fontweight='bold', y=0.96)
    plt.savefig(os.path.join(output_folder, 'HREF_LPMM_THRESHOLD_COMPARE.png'), dpi=300, bbox_inches='tight')

print("Process Complete.")
