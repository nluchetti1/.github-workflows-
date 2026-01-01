import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import pygrib, requests, os, numpy as np
import metpy.calc as mpcalc
from metpy.plots import Hodograph, USCOUNTIES
from metpy.units import units
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from datetime import datetime, timedelta
import pytz

# --- CONFIGURATION ---
now = datetime.now(pytz.UTC)
DATE = (now - timedelta(hours=6)).strftime('%Y%m%d')
CYCLE = "12" if now.hour >= 16 else "00"
EXTENT = [-92.0, -74.0, 24.5, 38.5] # Southeast US
OUTPUT_DIR = "hodo_data"
os.makedirs(OUTPUT_DIR, exist_ok=True)

def process_hour(f_hour):
    f_str = f"{f_hour:02d}"
    url = f"https://nomads.ncep.noaa.gov/pub/data/nccf/com/href/prod/href.{DATE}/ensprod/href.t{CYCLE}z.conus.mean.f{f_str}.grib2"
    
    r = requests.get(url, timeout=45)
    if r.status_code != 200: return
    with open("temp.grib2", "wb") as f: f.write(r.content)
    grbs = pygrib.open("temp.grib2")

    # 1. Grab MUCAPE for Masking and Background
    try:
        # Using SFC CAPE as proxy for MUCAPE background
        cape_grb = grbs.select(shortName='cape', level=0)[0]
        cape = cape_grb.values
        lats, lons = cape_grb.latlons()
    except: return

    # 2. Extract Vertical Profile (Pressure, Height, Wind)
    levels = [1000, 925, 850, 700, 500, 400, 300, 250, 200]
    u_l, v_l, h_l, p_l = [], [], [], []
    for lev in levels:
        try:
            u_l.append(grbs.select(shortName='u', level=lev)[0].values)
            v_l.append(grbs.select(shortName='v', level=lev)[0].values)
            h_l.append(grbs.select(shortName='gh', level=lev)[0].values)
            p_l.append(np.full(u_l[0].shape, lev))
        except: continue 

    u_stack = np.array(u_l) * units('m/s')
    v_stack = np.array(v_l) * units('m/s')
    h_stack = np.array(h_l) * units('m')
    p_stack = np.array(p_l) * units('hPa')

    # --- PLOTTING ---
    fig = plt.figure(figsize=(16, 10), facecolor='white')
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    ax.set_extent(EXTENT)
    
    # SPC-Style CAPE Color Levels
    # Blue -> Green -> Yellow -> Orange -> Red -> Pink
    cape_levels = [100, 250, 500, 750, 1000, 1500, 2000, 2500, 3000, 4000, 5000, 8000]
    cape_colors = ['#CCFFFF', '#99FFFF', '#66CCFF', '#0099FF', '#33FF00', '#33CC00', 
                   '#FFFF00', '#FFCC00', '#FF9900', '#FF0000', '#FF00FF', '#9900CC']
    
    ax.contourf(lons, lats, cape, levels=cape_levels, colors=cape_colors, alpha=0.4)
    ax.add_feature(USCOUNTIES.with_scale('500k'), edgecolor='black', linewidth=0.3, alpha=0.2)

    # 3. Grid Loop with SPC Logic
    skip = 45 
    for i in range(0, lats.shape[0], skip):
        for j in range(0, lats.shape[1], skip):
            # Only plot if CAPE >= 100 J/kg
            if cape[i,j] < 100 or not (EXTENT[0] <= lons[i,j] <= EXTENT[1]): continue

            try:
                # Calculate Bunkers Storm Motion
                rm, lm, mw = mpcalc.bunkers_storm_motion(p_stack[:,i,j], u_stack[:,i,j], v_stack[:,i,j], h_stack[:,i,j])
                
                # Create Hodograph Inset
                ax_ins = inset_axes(ax, width="0.5in", height="0.5in", bbox_to_anchor=(lons[i,j], lats[i,j]), bbox_transform=ax.transData, loc='center')
                hodo = Hodograph(ax_ins, component_range=60)
                hodo.add_grid(increment=20, color='gray', alpha=0.5, linewidth=0.5)
                
                # SPC Intervals: 0-500m (magenta), 0.5-3km (red), 3-6km (green), 6-9km (yellow)
                hodo.plot_colormapped(u_stack[:,i,j].to('kt'), v_stack[:,i,j].to('kt'), h_stack[:,i,j] - h_stack[0,i,j],
                                     intervals=[0, 500, 3000, 6000, 9000] * units.m,
                                     colors=['#ff00ff', '#ff0000', '#00ff00', '#ffff00'], linewidth=1.2)
                
                # Plot Storm Motion: Right (Red Circle), Left (Blue Circle)
                ax_ins.plot(rm[0].to('kt'), rm[1].to('kt'), marker='o', color='red', markersize=1.5)
                ax_ins.plot(lm[0].to('kt'), lm[1].to('kt'), marker='o', color='blue', markersize=1.5)
                
                ax_ins.axis('off')
            except: continue

    plt.title(f"SPC-Style Mesoanalysis | HREF Mean F{f_str}", loc='left', fontweight='bold', pad=10)
    plt.savefig(f"{OUTPUT_DIR}/hodo_f{f_str}.png", dpi=150, bbox_inches='tight')
    plt.close(); grbs.close()

for h in range(1, 25): process_hour(h)
