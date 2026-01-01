import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import pygrib, requests, os, numpy as np
import metpy.calc as mpcalc
from metpy.plots import Hodograph, USCOUNTIES
from metpy.units import units
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from datetime import datetime, timedelta
import pytz

# --- CONFIGURATION ---
now = datetime.now(pytz.UTC)
DATE_STR = "20260101" # Fixed for your current test run
CYCLE = "00"
START_TIME = datetime.strptime(f"{DATE_STR}{CYCLE}", "%Y%m%d%H").replace(tzinfo=pytz.UTC)

EXTENT = [-92.0, -74.0, 24.5, 38.5] # Southeast US Domain
OUTPUT_DIR = "hodo_data"
os.makedirs(OUTPUT_DIR, exist_ok=True)

def process_hour(f_hour):
    f_str = f"{f_hour:02d}"
    valid_time = START_TIME + timedelta(hours=f_hour)
    time_label = valid_time.strftime('%m/%d/%Y %H:%M UTC')
    
    url = f"https://nomads.ncep.noaa.gov/pub/data/nccf/com/href/prod/href.{DATE_STR}/ensprod/href.t{CYCLE}z.conus.mean.f{f_str}.grib2"
    
    try:
        r = requests.get(url, timeout=60)
        if r.status_code != 200: return
        with open("temp.grib2", "wb") as f: f.write(r.content)
    except: return

    grbs = pygrib.open("temp.grib2")

    # 1. Background SBCAPE (Index 4)
    cape_msg = grbs.select(shortName='cape', level=0)[0]
    cape = cape_msg.values
    lats, lons = cape_msg.latlons()

    # 2. SURFACE DATA (Required for AGL heights)
    # Surface GH (Index 43) and 10m Wind Speed (Index 52)
    h_sfc = grbs.select(shortName='gh', level=0)[0].values * units('m')
    # Note: Using 925mb as a surface proxy if 10m U/V are missing in mean files
    u_sfc = grbs.select(shortName='u', level=925)[0].values * units('m/s')
    v_sfc = grbs.select(shortName='v', level=925)[0].values * units('m/s')

    # 3. VERTICAL LEVELS (925, 850, 700, 500, 250)
    avail_levels = [925, 850, 700, 500, 250]
    u_list, v_list, h_list, p_list = [u_sfc], [v_sfc], [h_sfc], [np.full(u_sfc.shape, 1013)]
    
    for lev in avail_levels:
        u_list.append(grbs.select(shortName='u', level=lev)[0].values)
        v_list.append(grbs.select(shortName='v', level=lev)[0].values)
        h_list.append(grbs.select(shortName='gh', level=lev)[0].values)
        p_list.append(np.full(u_sfc.shape, lev))

    u_stack = np.array(u_list) * units('m/s')
    v_stack = np.array(v_list) * units('m/s')
    h_stack = np.array(h_list) * units('m')
    p_stack = np.array(p_list) * units('hPa')

    # --- PLOTTING ---
    fig = plt.figure(figsize=(18, 12), facecolor='white')
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    ax.set_extent(EXTENT, crs=ccrs.PlateCarree())
    
    # Border Adjustments (Professional Weight)
    ax.add_feature(cfeature.STATES, edgecolor='black', linewidth=1.0, zorder=10)
    ax.add_feature(USCOUNTIES.with_scale('500k'), edgecolor='black', linewidth=0.3, alpha=0.3, zorder=10)
    
    # Pivotal-Style Palette
    cape_levels = [100, 500, 1000, 1500, 2000, 2500, 3000, 4000, 5000]
    cf = ax.contourf(lons, lats, cape, levels=cape_levels, cmap='PuBuGn', alpha=0.4, zorder=1)
    cbar = plt.colorbar(cf, orientation='horizontal', pad=0.03, aspect=50)
    cbar.set_label('SBCAPE (J/kg)', fontweight='bold')

    # Grid Loop
    skip = 40 
    for i in range(0, lats.shape[0], skip):
        for j in range(0, lats.shape[1], skip):
            lon, lat = lons[i,j], lats[i,j]
            if not (EXTENT[0] <= lon <= EXTENT[1] and EXTENT[2] <= lat <= EXTENT[3]): continue

            try:
                # Calculate Bunkers Motion
                rm, lm, mw = mpcalc.bunkers_storm_motion(p_stack[:,i,j], u_stack[:,i,j], v_stack[:,i,j], h_stack[:,i,j])
                
                # Inset Placement
                ax_ins = inset_axes(ax, width="0.5in", height="0.5in", 
                                    bbox_to_anchor=(lon, lat), 
                                    bbox_transform=ax.transData, loc='center', borderpad=0)
                
                h = Hodograph(ax_ins, component_range=60)
                h.add_grid(increment=20, color='gray', alpha=0.4, linewidth=0.5)
                
                # Plot with AGL heights (Height - Surface Height)
                h.plot_colormapped(u_stack[:,i,j].to('kt'), v_stack[:,i,j].to('kt'), h_stack[:,i,j] - h_sfc[i,j],
                                     intervals=[0, 1000, 3000, 6000, 9000] * units.m,
                                     colors=['#ff00ff', '#ff0000', '#00ff00', '#ffff00'], linewidth=1.5)
                
                ax_ins.plot(rm[0].to('kt'), rm[1].to('kt'), marker='o', color='red', markersize=1.5)
                ax_ins.plot(lm[0].to('kt'), lm[1].to('kt'), marker='o', color='blue', markersize=1.5)
                ax_ins.axis('off')
            except: continue

    plt.title(f"HREF Mean SE Dynamics | Valid: {time_label}", loc='left', fontweight='bold')
    plt.savefig(f"{OUTPUT_DIR}/hodo_f{f_str}.png", dpi=120, bbox_inches='tight')
    plt.close(); grbs.close()

# Quick test hours
test_hours = [1, 6, 12, 18, 24]
for hr in test_hours: process_hour(hr)
