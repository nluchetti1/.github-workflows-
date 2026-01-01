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

# --- 1. CONFIGURATION ---
now = datetime.now(pytz.UTC)
DATE_STR = now.strftime('%Y%m%d')
CYCLE = "00" 
START_TIME = datetime.strptime(f"{DATE_STR}{CYCLE}", "%Y%m%d%H").replace(tzinfo=pytz.UTC)
EXTENT = [-92.0, -74.0, 24.5, 38.5] 
OUTPUT_DIR = "hodo_data"
os.makedirs(OUTPUT_DIR, exist_ok=True)

def process_hour(f_hour):
    # Fix for NameError: define time_label immediately
    valid_time = START_TIME + timedelta(hours=f_hour)
    time_label = valid_time.strftime('%m/%d/%Y %H:%M UTC')
    f_str = f"{f_hour:02d}"
    
    url = f"https://nomads.ncep.noaa.gov/pub/data/nccf/com/href/prod/href.{DATE_STR}/ensprod/href.t{CYCLE}z.conus.mean.f{f_str}.grib2"
    
    print(f">>> PROCESSING F{f_str} | Valid: {time_label} <<<")
    try:
        r = requests.get(url, timeout=60)
        if r.status_code != 200: return
        with open("temp.grib2", "wb") as f: f.write(r.content)
        grbs = pygrib.open("temp.grib2")
    except Exception as e:
        print(f"File failed: {e}")
        return

    # 2. DATA EXTRACTION
    try:
        cape = grbs.select(shortName='cape', level=0)[0].values
        lats, lons = grbs.select(shortName='cape', level=0)[0].latlons()
        h_sfc = grbs.select(shortName='gh', level=0)[0].values * units('m')
    except: return

    # Vertical levels available in HREF Mean: 925, 850, 700, 500, 250 mb
    avail_levels = [925, 850, 700, 500, 250]
    u_l, v_l, h_l, p_l = [], [], [], []
    for lev in avail_levels:
        u_l.append(grbs.select(shortName='u', level=lev)[0].values)
        v_l.append(grbs.select(shortName='v', level=lev)[0].values)
        h_l.append(grbs.select(shortName='gh', level=lev)[0].values)
        p_l.append(np.full(u_l[0].shape, lev))

    u_stack = np.array(u_l) * units('m/s')
    v_stack = np.array(v_l) * units('v/s')
    h_stack = np.array(h_l) * units('m')
    p_stack = np.array(p_l) * units('hPa')

    # 3. MAPPING
    fig = plt.figure(figsize=(18, 12), facecolor='white')
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    ax.set_extent(EXTENT, crs=ccrs.PlateCarree())
    
    ax.add_feature(cfeature.STATES, edgecolor='black', linewidth=0.8, zorder=10)
    ax.add_feature(USCOUNTIES.with_scale('500k'), edgecolor='black', linewidth=0.2, alpha=0.3, zorder=10)
    
    clevs = [100, 500, 1000, 1500, 2000, 2500, 3000, 4000, 5000]
    cf = ax.contourf(lons, lats, cape, levels=clevs, cmap='magma', alpha=0.4, zorder=1)
    plt.colorbar(cf, orientation='horizontal', pad=0.03, aspect=50, label='SBCAPE (J/kg)')

    # 4. HODOGRAPH GRID LOOP
    skip = 42 
    for i in range(0, lats.shape[0], skip):
        for j in range(0, lats.shape[1], skip):
            lon, lat = lons[i,j], lats[i,j]
            if not (EXTENT[0] <= lon <= EXTENT[1] and EXTENT[2] <= lat <= EXTENT[3]): continue

            try:
                # COORDINATE LOCK: bbox_transform=ax.transData pins inset to Lat/Lon
                ax_ins = inset_axes(ax, width="0.45in", height="0.45in", 
                                    bbox_to_anchor=(lon, lat), 
                                    bbox_transform=ax.transData, loc='center', borderpad=0)
                
                h = Hodograph(ax_ins, component_range=60)
                h.add_grid(increment=20, color='gray', alpha=0.5, linewidth=0.5)
                
                # Plot AGL segments (Height minus Ground GH)
                h.plot_colormapped(u_stack[:,i,j].to('kt'), v_stack[:,i,j].to('kt'), h_stack[:,i,j] - h_sfc[i,j],
                                     intervals=[0, 1000, 3000, 6000, 9000] * units.m,
                                     colors=['#ff00ff', '#ff0000', '#00ff00', '#ffff00'], linewidth=1.5)
                
                # Bunkers Storm Motion
                rm, lm, mw = mpcalc.bunkers_storm_motion(p_stack[:,i,j], u_stack[:,i,j], v_stack[:,i,j], h_stack[:,i,j])
                ax_ins.plot(rm[0].to('kt'), rm[1].to('kt'), 'ro', markersize=1.2)
                ax_ins.plot(lm[0].to('kt'), lm[1].to('kt'), 'bo', markersize=1.2)
                ax_ins.axis('off')
            except: continue

    plt.title(f"HREF Mean SE Dynamics | Valid: {time_label}", loc='left', fontweight='bold')
    plt.savefig(f"{OUTPUT_DIR}/hodo_f{f_str}.png", dpi=120, bbox_inches='tight')
    plt.close(); grbs.close()

for hr in [1, 6, 12, 18, 24]: process_hour(hr)
