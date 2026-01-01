import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import pygrib, requests, os, numpy as np
import metpy.calc as mpcalc
from metpy.plots import Hodograph, USCOUNTIES
from metpy.units import units
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, Bbox
from datetime import datetime, timedelta
import pytz
import shutil

# --- CONFIGURATION ---
now = datetime.now(pytz.UTC)
DATE_STR = now.strftime('%Y%m%d')
CYCLE = "00" 
START_TIME = datetime.strptime(f"{DATE_STR}{CYCLE}", "%Y%m%d%H").replace(tzinfo=pytz.UTC)
EXTENT = [-92.0, -74.0, 24.5, 38.5] 
OUTPUT_DIR = "hodo_data"

# FORCE CLEANUP: Delete old folder to ensure we aren't diffing against bad files
if os.path.exists(OUTPUT_DIR):
    shutil.rmtree(OUTPUT_DIR)
os.makedirs(OUTPUT_DIR, exist_ok=True)

def process_hour(f_hour):
    f_str = f"{f_hour:02d}"
    valid_time = START_TIME + timedelta(hours=f_hour)
    time_label = valid_time.strftime('%m/%d/%Y %H:%M UTC')
    
    url = f"https://nomads.ncep.noaa.gov/pub/data/nccf/com/href/prod/href.{DATE_STR}/ensprod/href.t{CYCLE}z.conus.mean.f{f_str}.grib2"
    print(f">>> PROCESSING F{f_str}...")

    try:
        r = requests.get(url, timeout=60)
        with open("temp.grib2", "wb") as f: f.write(r.content)
        grbs = pygrib.open("temp.grib2")
    except Exception as e:
        print(f"Download Failed: {e}"); return

    # 1. DATA EXTRACTION
    try:
        cape_msg = grbs.select(shortName='cape', level=0)[0]
        cape = cape_msg.values
        lats, lons = cape_msg.latlons()
        h_sfc = grbs.select(shortName='gh', level=0)[0].values * units('m')
    except Exception as e:
        print(f"Data Read Error: {e}"); return

    # Vertical Profile
    avail_levels = [925, 850, 700, 500, 250]
    u_l, v_l, h_l, p_l = [], [], [], []
    for lev in avail_levels:
        u_l.append(grbs.select(shortName='u', level=lev)[0].values)
        v_l.append(grbs.select(shortName='v', level=lev)[0].values)
        h_l.append(grbs.select(shortName='gh', level=lev)[0].values)
        p_l.append(np.full(u_l[0].shape, lev))

    u_stack = np.array(u_l) * units('m/s')
    v_stack = np.array(v_l) * units('m/s')
    h_stack = np.array(h_l) * units('m')
    p_stack = np.array(p_l) * units('hPa')

    # 2. PLOTTING
    fig = plt.figure(figsize=(18, 12), facecolor='white')
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    ax.set_extent(EXTENT, crs=ccrs.PlateCarree())
    
    ax.add_feature(cfeature.STATES, edgecolor='black', linewidth=0.8, zorder=10)
    ax.add_feature(USCOUNTIES.with_scale('500k'), edgecolor='black', linewidth=0.2, alpha=0.3, zorder=10)
    
    # FIX 1: Add transform=ccrs.PlateCarree() so colors actually appear!
    clevs = [100, 500, 1000, 1500, 2000, 2500, 3000, 4000, 5000]
    cf = ax.contourf(lons, lats, cape, levels=clevs, cmap='magma', alpha=0.4, 
                     transform=ccrs.PlateCarree(), zorder=1)
    plt.colorbar(cf, orientation='horizontal', pad=0.03, aspect=50, label='SBCAPE (J/kg)')

    # 3. HODOGRAPH LOOP
    skip = 42 
    success_count = 0
    
    for i in range(0, lats.shape[0], skip):
        for j in range(0, lats.shape[1], skip):
            lon, lat = lons[i,j], lats[i,j]
            if not (EXTENT[0] <= lon <= EXTENT[1] and EXTENT[2] <= lat <= EXTENT[3]): continue
            if np.isnan(u_stack[:,i,j].magnitude).any(): continue

            try:
                # FIX 2: THE ANCHOR FIX
                # We use a explicit Bbox centered at the lat/lon point
                # This solves the "Using relative units..." error
                ax_ins = inset_axes(ax, width="0.45in", height="0.45in", 
                                    bbox_to_anchor=(lon, lat), 
                                    bbox_transform=ax.transData, 
                                    loc='center', borderpad=0)
                
                h = Hodograph(ax_ins, component_range=60)
                h.add_grid(increment=20, color='gray', alpha=0.5, linewidth=0.5)
                
                h.plot_colormapped(u_stack[:,i,j].to('kt'), v_stack[:,i,j].to('kt'), h_stack[:,i,j] - h_sfc[i,j],
                                     intervals=[0, 1000, 3000, 6000, 9000] * units.m,
                                     colors=['#ff00ff', '#ff0000', '#00ff00', '#ffff00'], linewidth=1.5)
                
                rm, lm, mw = mpcalc.bunkers_storm_motion(p_stack[:,i,j], u_stack[:,i,j], v_stack[:,i,j], h_stack[:,i,j])
                ax_ins.plot(rm[0].to('kt'), rm[1].to('kt'), 'ro', markersize=1.2)
                ax_ins.plot(lm[0].to('kt'), lm[1].to('kt'), 'bo', markersize=1.2)
                ax_ins.axis('off')
                success_count += 1
            except Exception as e:
                if success_count == 0: print(f"Plot Error at {lat:.2f},{lon:.2f}: {e}")
                continue

    print(f"   -> Plotted {success_count} hodographs.")
    
    out_path = f"{OUTPUT_DIR}/hodo_f{f_str}.png"
    plt.title(f"HREF Mean SE Dynamics | Valid: {time_label}", loc='left', fontweight='bold')
    plt.savefig(out_path, dpi=120, bbox_inches='tight')
    plt.close()

for hr in [1, 6, 12, 18, 24]: process_hour(hr)
