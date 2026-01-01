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
REGION = [-98, -74, 24, 38]  # SE US [West, East, South, North]
GRID_SPACING = 8             # Plot a hodograph every Nth grid point
LEVELS = [1000, 925, 850, 700, 500, 300] * units.hPa  # Pressure levels
OUTPUT_DIR = "images"

def get_latest_run_time():
    """Determines the latest available 00Z or 12Z run."""
    now = datetime.datetime.utcnow()
    # Check if 12Z is likely available (after ~14:30Z)
    if now.hour >= 15:
        run = '12'
        date = now
    # Check if 00Z is likely available (after ~02:30Z)
    elif now.hour >= 3:
        run = '00'
        date = now
    else:
        # Fallback to previous day's 12Z if too early for 00Z
        run = '12'
        date = now - datetime.timedelta(days=1)
    
    return date.strftime('%Y%m%d'), run

def download_href_file(date_str, run, fhr='00'):
    """Downloads the HREF Ensemble Mean GRIB2 file from NOMADS."""
    # URL structure for HREF v3 Ensemble Mean
    base_url = f"https://nomads.ncep.noaa.gov/pub/data/nccf/com/href/prod/href.{date_str}/ensprod"
    filename = f"href.t{run}z.conus.mean.f{fhr}.grib2"
    url = f"{base_url}/{filename}"
    
    print(f"Downloading {url}...")
    response = requests.get(url, stream=True)
    if response.status_code == 200:
        with open(filename, 'wb') as f:
            for chunk in response.iter_content(chunk_size=1024):
                f.write(chunk)
        return filename
    else:
        print(f"Failed to download. Status: {response.status_code}")
        sys.exit(1)

def plot_hodographs(filename, date_str, run):
    """Reads data and plots hodographs."""
    # Load dataset (filter for isobaric levels)
    # backend_kwargs helps filter specific GRIB messages if needed
    ds = xr.open_dataset(filename, engine='cfgrib', 
                         filter_by_keys={'typeOfLevel': 'isobaricInhPa'})
    
    # Subset for region and levels
    ds = ds.sel(latitude=slice(REGION[3], REGION[2]), 
                longitude=slice(360+REGION[0], 360+REGION[1]))
    ds = ds.sel(isobaricInhPa=LEVELS.m)

    # Extract U and V components
    u = ds['u'].metpy.convert_units('kts')
    v = ds['v'].metpy.convert_units('kts')
    
    # Create Figure
    fig = plt.figure(figsize=(18, 12))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.LambertConformal())
    ax.set_extent(REGION)
    ax.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=2)
    ax.add_feature(cfeature.BORDERS, linewidth=2)
    ax.add_feature(cfeature.STATES.with_scale('50m'), linewidth=1, edgecolor='gray')

    # Decimate grid for plotting (skip points)
    skip = GRID_SPACING
    lons = u.longitude.values[::skip, ::skip]
    lats = u.latitude.values[::skip, ::skip]
    u_data = u.values[:, ::skip, ::skip]
    v_data = v.values[:, ::skip, ::skip]

    # Plot Hodographs
    # Iterate over the spatial grid
    for i in range(lons.shape[0]):
        for j in range(lons.shape[1]):
            # Create inset axes for each hodograph
            # Transform lat/lon to figure coordinates is complex, 
            # simplified approach: map coordinates
            
            # Simple check to ensure we are plotting valid data
            if np.isnan(u_data[0, i, j]): continue

            # Create an inset ax at the map location
            # Note: This loop can be slow. 
            # For efficiency, we place axes manually or use a custom collection.
            # Here we use the simplified `inset_axes` approach or just standard axes placement
            # But calculating transform is tricky.
            
            # ALTERNATIVE: Use the map coordinates to place the inset
            proj_pnt = ax.projection.transform_point(lons[i, j], lats[i, j], ccrs.PlateCarree())
            
            # Check if point is inside plot extent
            if not ax.get_xlim()[0] < proj_pnt[0] < ax.get_xlim()[1]: continue
            if not ax.get_ylim()[0] < proj_pnt[1] < ax.get_ylim()[1]: continue

            # Define inset size (coord units)
            inset_size = 0.8  # degrees approx equivalent
            bounds = [proj_pnt[0] - inset_size/2, proj_pnt[1] - inset_size/2, inset_size, inset_size]
            
            # Create the inset axes
            sub_ax = ax.inset_axes(bounds, transform=ax.projection)
            
            h = Hodograph(sub_ax, component_range=60)
            h.add_grid(increment=20, color='gray', alpha=0.5, linewidth=0.5)
            h.plot(u_data[:, i, j], v_data[:, i, j], linewidth=1.5, color='red')
            
            # Hide tick labels for cleanliness
            sub_ax.set_xticklabels([])
            sub_ax.set_yticklabels([])
            sub_ax.axis('off') # Turn off box

    # Title and Save
    plt.title(f"HREF Ensemble Mean Hodographs | Run: {date_str} {run}Z", fontsize=20, weight='bold')
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    output_path = f"{OUTPUT_DIR}/href_hodo_{date_str}_{run}z.png"
    plt.savefig(output_path, bbox_inches='tight', dpi=100)
    print(f"Saved {output_path}")

if __name__ == "__main__":
    date_str, run = get_latest_run_time()
    grib_file = download_href_file(date_str, run)
    try:
        plot_hodographs(grib_file, date_str, run)
    finally:
        # Cleanup
        if os.path.exists(grib_file):
            os.remove(grib_file)
