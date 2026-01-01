import pygrib, requests, os
from datetime import datetime, timedelta
import pytz

# --- CONFIGURATION ---
now = datetime.now(pytz.UTC)
DATE_STR = (now - timedelta(hours=6)).strftime('%Y%m%d')
CYCLE = "12" if now.hour >= 16 else "00"

# Target F01 for a quick check
url = f"https://nomads.ncep.noaa.gov/pub/data/nccf/com/href/prod/href.{DATE_STR}/ensprod/href.t{CYCLE}z.conus.mean.f01.grib2"

print(f"--- DIAGNOSING GRIB AT: {url} ---")

try:
    r = requests.get(url, timeout=30)
    with open("diag.grib2", "wb") as f: f.write(r.content)
    grbs = pygrib.open("diag.grib2")
    
    # This loop prints every message in the file to your console
    print(f"{'INDEX':<5} | {'NAME':<40} | {'SHORT':<10} | {'LEVEL':<10}")
    print("-" * 75)
    for i, g in enumerate(grbs):
        print(f"{i:<5} | {g.name:<40} | {g.shortName:<10} | {g.level:<10}")
        
    grbs.close()
except Exception as e:
    print(f"FAILED TO READ GRIB: {e}")
