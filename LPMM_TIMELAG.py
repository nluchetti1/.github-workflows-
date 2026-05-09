import matplotlib
matplotlib.use('Agg')
import os
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
import glob

# --- 1. DYNAMIC FOLDER SETUP ---
output_folder = "href_data"
os.makedirs(output_folder, exist_ok=True)

# --- 2. DOWNLOAD & UTILITY FUNCTIONS ---
def download_file(url, save_path):
    headers = {"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64)"}
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
    """
    Extracts data and coordinate grid.
    Crucial for fixing missing geography and titles.
    """
    try:
        with pygrib.open(grib_path) as grb_file:
            msgs = grb_file.select(parameterCategory=1, parameterNumber=8)
            if msgs:
                msg = msgs[0]
                data = msg.values
                lats, lons = msg.latlons() # Fix for blank maps
                return data, lats, lons
    except Exception as e:
        print(f"Unpack error: {e}")
    return None, None, None

def clean_up_grib_files(directory):
    """Deletes temporary GRIB files to keep repo clean."""
    files = glob.glob(os.path.join(directory, "*.grib2"))
    for f in files:
        try:
            os.remove(f)
        except:
            pass

# --- 3. TIMING & RUN SELECTION (NC Domain) ---
now_utc = datetime.now(pytz.UTC)
current_hour = now_utc.hour

base_url_href = "https://nomads.ncep.noaa.gov/pub/data/nccf/com/href/prod"
# Direct AWS S3 bucket path for REFS
base_url_refs = "https://noaa-rrfs-pds.s3.amazonaws.com/rrfs_a" 

date_now = now_utc.strftime('%Y%m%d')
date_prev = (now_utc - timedelta(days=1)).strftime('%Y%m%d')

# Logic for 00Z vs 12Z dprog/dt comparison (HREF - 12hr cadence, 24hr accumulation)
if 13 <= current_hour <= 23:
    target_start = now_utc.replace(hour=12, minute=0, second=0, microsecond=0)
    href_runs = [{"date": date_now, "hour": "12", "f_range": range(1, 25)},
                 {"date": date_now, "hour": "00", "f_range": range(13, 37)},
                 {"date": date_prev, "hour": "12", "f_range": range(25, 49)}]
else:
    target_start = now_utc.replace(hour=0, minute=0, second=0, microsecond=0)
    href_runs = [{"date": date_now, "hour": "00", "f_range": range(1, 25)},
                 {"date": date_prev, "hour": "12", "f_range": range(13, 37)},
                 {"date": date_prev, "hour": "00", "f_range": range(25, 49)}]

href_valid_range = f"Valid: {target_start.strftime('%H')}Z {target_start.strftime('%Y%m%d')} to {(target_start+timedelta(hours=24)).strftime('%H')}Z {(target_start+timedelta(hours=24)).strftime('%Y%m%d')}"
refs_valid_range = f"Valid: {target_start.strftime('%H')}Z {target_start.strftime('%Y%m%d')} to {(target_start+timedelta(hours=48)).strftime('%H')}Z {(target_start+timedelta(hours=48)).strftime('%Y%m%d')}"

# Logic for REFS dprog/dt comparison (6hr cadence, 48hr accumulation)
refs_runs = []
for i in range(3):
    run_time = target_start - timedelta(hours=6*i)
    hour_diff = int((target_start - run_time).total_seconds() / 3600)
    refs_runs.append({
        "date": run_time.strftime('%Y%m%d'),
        "hour": f"{run_time.hour:02d}",
        "f_range": range(hour_diff + 1, hour_diff + 49), # 48-hour accumulation maxing out at f60
        "time": run_time
    })

# --- 4. DATA COLLECTION ---
href_results = []
refs_results = []
mm_to_inch = 0.0393701

href_lats, href_lons = None, None
refs_lats, refs_lons = None, None

# 4A. Fetching HREF Data
for idx, run in enumerate(href_runs):
    hourly_data = []
    print(f"Processing HREF Run: {run['date']} {run['hour']}Z")
    for f_hour in run['f_range']:
        f_str = f"{f_hour:02d}"
        url = f"{base_url_href}/href.{run['date']}/ensprod/href.t{run['hour']}z.conus.lpmm.f{f_str}.grib2"
        temp_path = os.path.join(output_folder, f"temp_href_{idx}_{f_str}.grib2")
        if download_file(url, temp_path):
            data, lats, lons = unpack_total_precipitation(temp_path)
            if data is not None:
                hourly_data.append(data * mm_to_inch)
                href_lats, href_lons = lats, lons
            if os.path.exists(temp_path): os.remove(temp_path)
    if hourly_data:
        href_results.append({"data": np.sum(hourly_data, axis=0), 
                             "time": pd.to_datetime(run['date'] + ' ' + run['hour'] + 'Z')})

# 4B. Fetching REFS Data from AWS S3
for idx, run in enumerate(refs_runs):
    hourly_data = []
    print(f"Processing REFS Run: {run['date']} {run['hour']}Z")
    for f_hour in run['f_range']:
        f_str = f"{f_hour:02d}"
        url = f"{base_url_refs}/refs.{run['date']}/{run['hour']}/enspost/refs.t{run['hour']}z.lpmm.f{f_str}.conus.grib2"
        temp_path = os.path.join(output_folder, f"temp_refs_{idx}_{f_str}.grib2")
        if download_file(url, temp_path):
            data, lats, lons = unpack_total_precipitation(temp_path)
            if data is not None:
                hourly_data.append(data * mm_to_inch)
                refs_lats, refs_lons = lats, lons
            if os.path.exists(temp_path): os.remove(temp_path)
    if hourly_data:
        refs_results.append({"data": np.sum(hourly_data, axis=0), 
                              "time": pd.to_datetime(run['date'] + ' ' + run['hour'] + 'Z')})

# --- 5. PLOTTING ---
clevs = [0.0, 0.01, 0.1, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10., 12., 15., 18.]
cmap_data = [(1,1,1), (0.31,0.81,0.81), (0,1,1), (0,0.87,0.5), (0,0.75,0), (0.5,0.87,0), (1,1,0), (1,0.62,0), (1,0,0), (1,0.12,0.5), (0.94,0.25,1), (0.5,0.12,1), (0.25,0.25,1), (0.12,0.12,0.5), (0.12,0.12,0.12), (0.5,0.5,0.5), (0.87,0.87,0.87), (0.93,0.83,0.73), (0.85,0.65,0.47), (0.62,0.42,0.23), (0.4,0.2,0)]
cmap = mcolors.ListedColormap(cmap_data, 'precip')
norm = mcolors.BoundaryNorm(clevs, cmap.N)

# A: HREF Comparison Plot
if len(href_results) >= 1 and href_lons is not None:
    fig, axes = plt.subplots(1, 3, figsize=(18, 7.5), subplot_kw={'projection': ccrs.PlateCarree()})
    fig.subplots_adjust(left=0.05, right=0.95, bottom=0.22, top=0.85, wspace=0.05)
    cs = None
    for i in range(3):
        if i < len(href_results):
            res = href_results[i]
            cs = axes[i].contourf(href_lons, href_lats, res['data'], clevs, cmap=cmap, norm=norm, alpha=0.5)
            axes[i].set_title(f'{res["time"].strftime("%Y-%m-%d %H:%M Z")}\n24hr HREF LPMM [in]', fontsize=10, fontweight='bold', pad=8)
        axes[i].add_feature(cfeature.STATES, linewidth=0.8, edgecolor='black')
        axes[i].add_feature(USCOUNTIES.with_scale('500k'), edgecolor='gray', linewidth=0.4)
        axes[i].set_extent([-84.8, -74, 31, 39]) # NC Domain
    
    if cs:
        cbar_ax = fig.add_axes([0.15, 0.12, 0.7, 0.03])
        cbar = fig.colorbar(cs, cax=cbar_ax, orientation='horizontal', ticks=clevs)
        cbar.ax.tick_params(labelsize=11)
        cbar.set_label('Precipitation (inches)', fontsize=15, fontweight='bold')
    
    fig.suptitle(f'24hr HREF LPMM [in] dprog/dt\n{href_valid_range}', fontsize=16, fontweight='bold', y=0.96)
    plt.savefig(os.path.join(output_folder, 'latest_compare.png'), dpi=300)

# B: REFS Comparison Plot (Updated to 48 hours)
if len(refs_results) >= 1 and refs_lons is not None:
    fig_refs, axes_refs = plt.subplots(1, 3, figsize=(18, 7.5), subplot_kw={'projection': ccrs.PlateCarree()})
    fig_refs.subplots_adjust(left=0.05, right=0.95, bottom=0.22, top=0.85, wspace=0.05)
    cs_refs = None
    for i in range(3):
        if i < len(refs_results):
            res = refs_results[i]
            cs_refs = axes_refs[i].contourf(refs_lons, refs_lats, res['data'], clevs, cmap=cmap, norm=norm, alpha=0.5)
            axes_refs[i].set_title(f'{res["time"].strftime("%Y-%m-%d %H:%M Z")}\n48hr REFS LPMM [in]', fontsize=10, fontweight='bold', pad=8)
        axes_refs[i].add_feature(cfeature.STATES, linewidth=0.8, edgecolor='black')
        axes_refs[i].add_feature(USCOUNTIES.with_scale('500k'), edgecolor='gray', linewidth=0.4)
        axes_refs[i].set_extent([-84.8, -74, 31, 39]) # NC Domain
    
    if cs_refs:
        cbar_ax_refs = fig_refs.add_axes([0.15, 0.12, 0.7, 0.03])
        cbar_refs = fig_refs.colorbar(cs_refs, cax=cbar_ax_refs, orientation='horizontal', ticks=clevs)
        cbar_refs.ax.tick_params(labelsize=11)
        cbar_refs.set_label('Precipitation (inches)', fontsize=15, fontweight='bold')
    
    fig_refs.suptitle(f'48hr REFS LPMM [in] dprog/dt\n{refs_valid_range}', fontsize=16, fontweight='bold', y=0.96)
    plt.savefig(os.path.join(output_folder, 'latest_refs_compare.png'), dpi=300)

# C: HREF Threshold Plot
if len(href_results) >= 1 and href_lons is not None:
    fig2, ax2 = plt.subplots(2, 2, figsize=(14, 11), subplot_kw={'projection': ccrs.PlateCarree()})
    blue_shades = ['#00008B', '#4169E1', '#87CEFA']
    legend_elements = [Line2D([0], [0], marker='o', color='w', label=res['time'].strftime("%Y-%m-%d %H:%M Z"),
                              markerfacecolor=blue_shades[idx], markersize=8) for idx, res in enumerate(href_results)]
    for i, thresh in enumerate([3, 6, 9, 12]):
        row, col = divmod(i, 2)
        for j, res in enumerate(href_results):
            m_data = np.ma.masked_less(res['data'], thresh)
            ax2[row, col].contourf(href_lons, href_lats, m_data, cmap=mcolors.ListedColormap([blue_shades[j]]), levels=[thresh, 99], alpha=0.6)
        ax2[row, col].add_feature(cfeature.STATES, linewidth=0.8, edgecolor='black')
        ax2[row, col].add_feature(USCOUNTIES.with_scale('500k'), edgecolor='gray', linewidth=0.3, alpha=0.5)
        ax2[row, col].set_extent([-84.8, -74, 31, 39]) # NC Domain
        ax2[row, col].set_title(f'> {thresh} inches', fontsize=12, fontweight='bold')
        ax2[row, col].legend(handles=legend_elements, loc='lower right', title='HREF Run', fontsize=8)
    
    fig2.subplots_adjust(left=0.05, right=0.95, bottom=0.05, top=0.88, wspace=0.1, hspace=0.2)
    fig2.suptitle(f'24hr HREF LPMM Threshold Compare\n{href_valid_range}', fontsize=16, fontweight='bold', y=0.96)
    plt.savefig(os.path.join(output_folder, 'latest_threshold.png'), dpi=300, bbox_inches='tight')

# --- 6. CLEANUP ---
clean_up_grib_files(output_folder)
print("Process Complete.")
