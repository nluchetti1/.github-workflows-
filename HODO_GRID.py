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
DATE_STR = "20260101" 
CYCLE = "00"
START_TIME = datetime.strptime(f"{DATE_STR}{CYCLE}", "%Y%m%d%H").replace(tzinfo=pytz.UTC)

EXTENT = [-92.0, -74.0, 24.5, 38.5] # SE US Domain
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

    # 1. Background SBCAPE (Index 4) - Reverted to Magma
    cape_msg = grbs.select(shortName='cape', level=0)[0]
    cape = cape_msg.values
    lats, lons = cape_msg.latlons()

    # 2. FORCE GROUND LAYER (Critical Fix)
    # Using Index 43 (GH at 0) and Index 52 (10m Wind Speed) as Ground
    h_sfc = grbs.select(shortName='gh', level=0)[0].values * units('m')
    u_sfc = grbs.select(shortName='u', level=925)[0].values * units('m/s') # Proxy for surface
    v_sfc = grbs.select(shortName='v', level=925)[0].values * units('m/s')

    # 3. BUILD PROFILE FROM AVAILABLE LEVELS
    avail_levels = [925, 850, 700, 500, 250]
    u_l, v_l, h_l, p_l = [u_sfc], [v_sfc], [h_sfc], [np.full(u_sfc.shape, 1013)]
    
    for lev in avail_levels:
        u_l.append(grbs.select(shortName='u', level=lev)[0].values)
        v_l.append(grbs.select(shortName='v', level=lev)[0].values)
        h_l.append(grbs.select(shortName='gh', level=lev)[0].values)
        p_l.append(np.full(u_sfc.shape, lev))

    u_stack = np.array(u_l) * units('m/s')
    v_stack = np.array(v_l) * units('m/s')
    h_stack = np.array(h_l) * units('m')
    p_stack = np.array(p_l) * units('hPa')

    # --- PLOTTING ---
    fig = plt.figure(figsize=(18, 12), facecolor='white')
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    ax.set_extent(EXTENT, crs=ccrs.PlateCarree())
    
    # Geography Borders
    ax.add_feature(cfeature.STATES, edgecolor='black', linewidth=0.8, zorder=10)
    ax.add_feature(USCOUNTIES.with_scale('500k'), edgecolor='black', linewidth=0.2, alpha=0.3, zorder=10)
    
    # CAPE Shading - Original Palette
    clevs = [100, 500, 1000, 1500, 2000, 2500, 3000, 4000, 5000]
    cf = ax.contourf(lons, lats, cape, levels=clevs, cmap='magma', alpha=0.4, zorder=1)
    plt.colorbar(cf, orientation='horizontal', pad=0.03, aspect=50, label='SBCAPE (J/kg)')

    # 4. GRID LOOP
    skip = 40 
    for i in range(0, lats.shape[0], skip):
        for j in range(0, lats.shape[1], skip):
            lon, lat = lons[i,j], lats[i,j]
            if not (EXTENT[0] <= lon <= EXTENT[1] and EXTENT[2] <= lat <= EXTENT[3]): continue

            try:
                # Calculate Bunkers Motion
                rm, lm, mw = mpcalc.bunkers_storm_motion(p_stack[:,i,j], u_stack[:,i,j], v_stack[:,i,j], h_stack[:,i,j])
                
                # Inset Placement with projection fix
                ax_ins = inset_axes(ax, width="0.45in", height="0.45in", 
                                    bbox_to_anchor=(lon, lat), 
                                    bbox_transform=ax.transData, loc='center', borderpad=0)
                
                h = Hodograph(ax_ins, component_range=60)
                h.add_grid(increment=20, color='gray', alpha=0.5, linewidth=0.5)
                
                # Plot Line AGL (Height minus Ground)
                h.plot_colormapped(u_stack[:,i,j].to('kt'), v_stack[:,i,j].to('kt'), h_stack[:,i,j] - h_sfc[i,j],
                                     intervals=[0, 1000, 3000, 6000, 9000] * units.m,
                                     colors=['#ff00ff', '#ff0000', '#00ff00', '#ffff00'], linewidth=1.5)
                
                # Bunkers Markers
                ax_ins.plot(rm[0].to('kt'), rm[1].to('kt'), marker='o', color='red', markersize=1.2)
                ax_ins.plot(lm[0].to('kt'), lm[1].to('kt'), marker='o', color='blue', markersize=1.2)
                ax_ins.axis('off')
            except: continue

    plt.title(f"HREF Mean SE Dynamics | {time_label}", loc='left', fontweight='bold')
    plt.savefig(f"{OUTPUT_DIR}/hodo_f{f_str}.png", dpi=120, bbox_inches='tight')
    plt.close(); grbs.close()

test_hours = [1, 6, 12, 18, 24]
for hr in test_hours: process_hour(hr)
