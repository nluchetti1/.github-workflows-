import matplotlib
matplotlib.use('Agg')
import os, numpy as np, requests, pygrib, matplotlib.pyplot as plt, matplotlib.colors as mcolors
import cartopy.crs as ccrs, cartopy.feature as cfeature, pandas as pd, pytz
from datetime import datetime, timedelta
from metpy.plots import USCOUNTIES

output_folder = "href_data" # Keep same folder for shared dashboard
os.makedirs(output_folder, exist_ok=True)

def download_file(url, save_path):
    try:
        response = requests.get(url, timeout=30)
        if response.status_code == 200:
            with open(save_path, 'wb') as f: f.write(response.content)
            return True
    except: pass
    return False

def unpack_precip(grib_path):
    try:
        with pygrib.open(grib_path) as grb:
            # RRFS LPMM identifies similarly to HREF in moisture categories
            msgs = grb.select(parameterCategory=1, parameterNumber=8)
            if msgs:
                msg = msgs[0]
                return msg.values, msg.latlons()
    except: pass
    return None, (None, None)

# TIMING: 00Z and 12Z runs
now_utc = datetime.now(pytz.UTC)
date_now = now_utc.strftime('%Y%m%d')
date_prev = (now_utc - timedelta(days=1)).strftime('%Y%m%d')

# RRFS URL structure
base_url = "https://nomads.ncep.noaa.gov/pub/data/nccf/com/rrfs/prod"
run_hour = "12" if 13 <= now_utc.hour <= 23 else "00"
date_run = date_now if run_hour == "12" or now_utc.hour >= 0 else date_prev

all_data = []
final_lats, final_lons = None, None

# Fetch 24 hours of LPMM data for the current cycle
for f_hour in range(1, 25):
    f_str = f"{f_hour:02d}"
    url = f"{base_url}/rrfs_a.{date_run}/{run_hour}/enspost_timelag/rrfs.t{run_hour}z.conus.lpmm.f{f_str}.grib2"
    temp_path = os.path.join(output_folder, f"rrfs_temp_{f_str}.grib2")
    if download_file(url, temp_path):
        data, (lats, lons) = unpack_precip(temp_path)
        if data is not None:
            all_data.append(data * 0.0393701) # mm to inch
            final_lats, final_lons = lats, lons
        if os.path.exists(temp_path): os.remove(temp_path)

if all_data:
    # Plotting code remains identical to HREF script, saving to 'latest_rrfs.png'
    total_precip = np.sum(all_data, axis=0)
    fig, ax = plt.subplots(figsize=(12, 8), subplot_kw={'projection': ccrs.PlateCarree()})
    # ... [Same plotting logic as HREF script] ...
    plt.savefig(os.path.join(output_folder, 'latest_rrfs.png'), dpi=150)
