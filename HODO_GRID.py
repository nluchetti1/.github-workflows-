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
DATE = (datetime.now(pytz.UTC) - timedelta(hours=6)).strftime('%Y%m%d')
CYCLE = "12" 
EXTENT = [-84.5, -75.0, 33.5, 37.5] # North Carolina
OUTPUT_DIR = "hodo_data"
os.makedirs(OUTPUT_DIR, exist_ok=True)

def process_hour(f_hour):
    f_str = f"{f_hour:02d}"
    url = f"https://nomads.ncep.noaa.gov/pub/data/nccf/com/href/prod/href.{DATE}/ensprod/href.t{CYCLE}z.conus.mean.f{f_str}.grib2"
    
    print(f"--- Processing Hour F{f_str} ---")
    r = requests.get(url, stream=True)
    if r.status_code != 200:
        print(f"Error: Could not download F{f_str} from NOMADS.")
        return
    
    with open("temp.grib2", "wb") as f: f.write(r.content)
    grbs = pygrib.open("temp.grib2")

    # SBCAPE Background (shortName 'cape')
    try:
        cape_grb = grbs.select(shortName='cape', level=0)[0]
        cape = cape_grb.values
        lats, lons = cape_grb.latlons()
    except ValueError:
        print("SBCAPE not found; check GRIB inventory.")
        return

    # Vertical Profile Data
    # We use shortNames 'u', 'v', and 'gh' (Geopotential Height)
    levels = [1000, 925, 850, 700, 500, 400, 300]
    u_list, v_list, h_list = [], [], []
    valid_levels = []

    for lev in levels:
        try:
            u_list.append(grbs.select(shortName='u', level=lev)[0].values)
            v_list.append(grbs.select(shortName='v', level=lev)[0].values)
            h_list.append(grbs.select(shortName='gh', level=lev)[0].values)
            valid_levels.append(lev)
        except (ValueError, IndexError):
            continue # Skip missing levels gracefully

    u_all = np.array(u_list) * units('m/s')
    v_all = np.array(v_list) * units('m/s')
    h_all = np.array(h_list) * units('m')

    # --- PLOTTING ---
    fig = plt.figure(figsize=(16, 10), facecolor='black')
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    ax.set_extent(EXTENT)
    
    # SBCAPE Shading
    clevs = [100, 500, 1000, 2000, 3000, 4000, 5000]
    cf = ax.contourf(lons, lats, cape, levels=clevs, cmap='magma', alpha=0.5)
    ax.add_feature(USCOUNTIES.with_scale('500k'), edgecolor='white', alpha=0.15)

    # GRID LOOPING
    skip = 35 # Density control
    for i in range(0, lats.shape[0], skip):
        for j in range(0, lats.shape[1], skip):
            
            # 0-6km Bulk Shear Calculation
            try:
                u_shr, v_shr = mpcalc.bulk_shear(h_all[:,i,j], u_all[:,i,j], v_all[:,i,j], depth=6000*units.m)
                ax.quiver(lons[i,j], lats[i,j], u_shr.to('kt').m, v_shr.to('kt').m, 
                          color='gold', scale=450, width=0.0025, headwidth=4, pivot='middle')
            except: continue

            # Inset Hodograph
            ax_ins = inset_axes(ax, width="0.45in", height="0.45in",
                                bbox_to_anchor=(lons[i,j], lats[i,j]),
                                bbox_transform=ax.transData, loc='center')
            
            hodo = Hodograph(ax_ins, component_range=60)
            hodo.add_grid(increment=20, color='white', alpha=0.15)
            hodo.plot_colormapped(u_all[:,i,j].to('kt'), v_all[:,i,j].to('kt'), h_all[:,i,j], 
                                 intervals=[0, 1000, 3000, 6000] * units.m,
                                 colors=['#ff00ff', '#ff0000', '#00ff00'])
            ax_ins.axis('off')

    plt.title(f"HREF Mean | F{f_str} SBCAPE, 0-6km Shear & Hodographs", color='white', loc='left', pad=10)
    plt.savefig(f"{OUTPUT_DIR}/hodo_f{f_str}.png", dpi=150, facecolor='black', bbox_inches='tight')
    plt.close()
    grbs.close()

# Execute first 24 forecast hours
for hour in range(1, 25):
    process_hour(hour)
