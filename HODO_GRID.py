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
DATE_STR = now.strftime('%Y%m%d')
# Target the most recent 00Z or 12Z cycle
CYCLE = "12" if now.hour >= 16 else "00"
START_TIME = datetime.strptime(f"{DATE_STR}{CYCLE}", "%Y%m%d%H").replace(tzinfo=pytz.UTC)

EXTENT = [-92.0, -74.0, 24.5, 38.5] # Southeast US
OUTPUT_DIR = "hodo_data"
os.makedirs(OUTPUT_DIR, exist_ok=True)

def process_hour(f_hour):
    f_str = f"{f_hour:02d}"
    valid_time = START_TIME + timedelta(hours=f_hour)
    time_label = valid_time.strftime('%m/%d/%Y %H:%M UTC')
    
    url = f"https://nomads.ncep.noaa.gov/pub/data/nccf/com/href/prod/href.{DATE_STR}/ensprod/href.t{CYCLE}z.conus.mean.f{f_str}.grib2"
    
    print(f"Processing F{f_str} (Valid: {time_label})...")
    try:
        r = requests.get(url, timeout=60)
        if r.status_code != 200: return
        with open("temp.grib2", "wb") as f: f.write(r.content)
    except: return

    grbs = pygrib.open("temp.grib2")

    # 1. Background SBCAPE
    try:
        cape_grb = grbs.select(shortName='cape', level=0)[0]
        cape = cape_grb.values
        lats, lons = cape_grb.latlons()
    except: return

    # 2. Extract Vertical Profile
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
    fig = plt.figure(figsize=(18, 12), facecolor='white')
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    ax.set_extent(EXTENT)
    
    # Bold Borders
    ax.add_feature(cfeature.STATES, edgecolor='black', linewidth=2.0, zorder=5)
    ax.add_feature(USCOUNTIES.with_scale('500k'), edgecolor='black', linewidth=0.6, alpha=0.5, zorder=5)
    
    # SPC-Style CAPE Palette
    cape_levels = [100, 250, 500, 750, 1000, 1500, 2000, 2500, 3000, 4000, 5000]
    cape_colors = ['#CCFFFF', '#99FFFF', '#66CCFF', '#0099FF', '#33FF00', '#33CC00', 
                   '#FFFF00', '#FFCC00', '#FF9900', '#FF0000', '#FF00FF']
    
    cf = ax.contourf(lons, lats, cape, levels=cape_levels, colors=cape_colors, alpha=0.4, zorder=1)
    cbar = plt.colorbar(cf, orientation='horizontal', pad=0.03, aspect=50)
    cbar.set_label('SBCAPE (J/kg)', fontweight='bold')
    
    # 3. Grid Loop
    skip = 45 
    for i in range(0, lats.shape[0], skip):
        for j in range(0, lats.shape[1], skip):
            if not (EXTENT[0] <= lons[i,j] <= EXTENT[1]): continue

            try:
                # Bunkers Storm Motion
                rm, lm, mw = mpcalc.bunkers_storm_motion(p_stack[:,i,j], u_stack[:,i,j], v_stack[:,i,j], h_stack[:,i,j])
                
                ax_ins = inset_axes(ax, width="0.55in", height="0.55in", bbox_to_anchor=(lons[i,j], lats[i,j]), bbox_transform=ax.transData, loc='center')
                h = Hodograph(ax_ins, component_range=60)
                h.add_grid(increment=20, color='gray', alpha=0.5, linewidth=0.7)
                
                # Height-Based color coding
                h.plot_colormapped(u_stack[:,i,j].to('kt'), v_stack[:,i,j].to('kt'), h_stack[:,i,j] - h_stack[0,i,j],
                                     intervals=[0, 500, 3000, 6000, 9000] * units.m,
                                     colors=['#ff00ff', '#ff0000', '#00ff00', '#ffff00'], linewidth=1.8)
                
                ax_ins.plot(rm[0].to('kt'), rm[1].to('kt'), marker='o', color='red', markersize=1.8, zorder=10)
                ax_ins.plot(lm[0].to('kt'), lm[1].to('kt'), marker='o', color='blue', markersize=1.8, zorder=10)
                ax_ins.axis('off')
            except: continue

    plt.title(f"HREF Mean SE Dynamics | Valid: {time_label}", loc='left', fontweight='bold', fontsize=16)
    plt.title(f"F{f_str} | Init: {DATE_STR} {CYCLE}Z", loc='right', fontsize=11)
    
    plt.savefig(f"{OUTPUT_DIR}/hodo_f{f_str}.png", dpi=150, bbox_inches='tight')
    plt.close(); grbs.close()

# LOOP THROUGH ALL 48 HOURS
for hr in range(1, 49): 
    process_hour(hr)
