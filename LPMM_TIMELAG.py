#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  8 10:31:43 2024

@author: nicholasluchetti
"""



#%%
import matplotlib
matplotlib.use('Agg')
import pytz  # Required for time zone handling if running locally
import os
import numpy as np
import requests
import pygrib
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import cartopy.crs as ccrs
import cartopy.feature as cfeature
#from cartopy.feature import USCOUNTIES
import pandas as pd

#from datetime import datetime, timedelta
from urllib.request import urlopen
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
from metpy.plots import USCOUNTIES
from metpy.units import masked_array, units
from netCDF4 import Dataset
import pandas as pd
import matplotlib.patheffects as PathEffects

import os
#import datetime

from datetime import datetime, timedelta

# Get current UTC time
now_utc = datetime.now(pytz.UTC)  # current time in UTC
current_hour = now_utc.hour

# Define which code block to run based on the time range
if 13 <= current_hour <= 23:

    # Define the base URL
    base_url = "https://nomads.ncep.noaa.gov/pub/data/nccf/com/href/prod/"

    current_date = datetime.now()
    date = current_date.strftime('%Y%m%d')

    # Specify the date and forecast hour
    #date = "20241108"  # Example date
    hour_str = "12"  # Example hour (can be 00Z or 12Z)
    forecast_hour = "01"  # First forecast hour

    # Define the output directory and file path
    output_dir = "href_data"
    os.makedirs(output_dir, exist_ok=True)

    # Conversion factor from mm to inches
    mm_to_inch = 0.0393701  # 1 mm = 0.0393701 inches

    # Create the plot projection
    pc_proj = ccrs.PlateCarree()

    # Initialize lists to hold hourly precipitation data
    all_hourly_precip = []

    # Function to unpack "Total precipitation" variable based on parameterName
    def unpack_total_precipitation(grib_path):
        try:
            with pygrib.open(grib_path) as grb_file:
                for grb_message in grb_file:
                    if grb_message.parameterName == "Total precipitation":
                        data, lats, lons = grb_message.data()
                        print(f"Extracted Total precipitation from {grib_path}")
                        return data, lats, lons
                print(f"Total precipitation not found in {grib_path}")
        except Exception as e:
            print(f"Error reading {grib_path}: {e}")
        return None, None, None

    # Function to download files with User-Agent header
    def download_file(url, save_path):
        headers = {
            "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/87.0.4280.88 Safari/537.36"
        }
        try:
            response = requests.get(url, headers=headers)
            response.raise_for_status()  # Check if request was successful
            with open(save_path, 'wb') as file:
                file.write(response.content)
            print(f"Downloaded {url}")
        except requests.exceptions.RequestException as e:
            print(f"Failed to download {url}: {e}")

    # Loop through forecast hours (from 01 to 24)
    for i in range(1, 24):
        forecast_hour_str = f"{i:02d}"  # e.g., '01', '02', ..., '24'
        file_url = f"{base_url}/href.{date}/ensprod/href.t{hour_str}z.conus.lpmm.f{forecast_hour_str}.grib2"
        save_path = os.path.join(output_dir, f"href_{date}_t{hour_str}_f{forecast_hour_str}.grib2")
        
        # Download the file for each forecast hour
        download_file(file_url, save_path)

        # Unpack the "Total precipitation" variable from the downloaded file
        hourly_precip, lats, lons = unpack_total_precipitation(save_path)
        
        # If precipitation data is successfully extracted, store it
        if hourly_precip is not None:
            hourly_precip_inch = hourly_precip * mm_to_inch
            all_hourly_precip.append(hourly_precip_inch)

    # Sum the hourly precipitation data for the 24 hours
    total_precip_1 = np.sum(all_hourly_precip, axis=0)

    # Create the figure and axes for plotting
    fig = plt.figure(figsize=(20, 20))
    ax = fig.add_subplot(111, projection=pc_proj)
    ax.set_extent([-74, -84.8, 31, 39])  # Area of interest

    # Draw coastlines, state and country boundaries, and edges of the map
    ax.coastlines(resolution='10m')
    ax.add_feature(cfeature.BORDERS.with_scale('10m'), linewidth=1.5)
    ax.add_feature(cfeature.STATES.with_scale('10m'), linewidth=2.0)
    ax.add_feature(USCOUNTIES.with_scale('500k'), edgecolor='black')

    # Define the contour levels for precipitation
    clevs = [0.00, 0.01, 0.1, 0.25, 0.50, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0,
             6.0, 8.0, 10., 12.0, 15.0, 18.0]

    # Define the color map for the precipitation levels
    cmap_data = [(1.0, 1.0, 1.0),
                 (0.3137255012989044, 0.8156862854957581, 0.8156862854957581),
                 (0.0, 1.0, 1.0),
                 (0.0, 0.8784313797950745, 0.501960813999176),
                 (0.0, 0.7529411911964417, 0.0),
                 (0.501960813999176, 0.8784313797950745, 0.0),
                 (1.0, 1.0, 0.0),
                 (1.0, 0.6274510025978088, 0.0),
                 (1.0, 0.0, 0.0),
                 (1.0, 0.125490203499794, 0.501960813999176),
                 (0.9411764740943909, 0.250980406999588, 1.0),
                 (0.501960813999176, 0.125490203499794, 1.0),
                 (0.250980406999588, 0.250980406999588, 1.0),
                 (0.125490203499794, 0.125490203499794, 0.501960813999176),
                 (0.125490203499794, 0.125490203499794, 0.125490203499794),
                 (0.501960813999176, 0.501960813999176, 0.501960813999176),
                 (0.8784313797950745, 0.8784313797950745, 0.8784313797950745),
                 (0.9333333373069763, 0.8313725590705872, 0.7372549176216125),
                 (0.8549019694328308, 0.6509804129600525, 0.47058823704719543),
                 (0.6274510025978088, 0.42352941632270813, 0.23529411852359772),
                 (0.4000000059604645, 0.20000000298023224, 0.0)]
    cmap = mcolors.ListedColormap(cmap_data, 'precipitation')
    norm = mcolors.BoundaryNorm(clevs, cmap.N)

    # Plot the contour fill
    cs = ax.contourf(lons, lats, total_precip_1, clevs, alpha=0.5, cmap=cmap, norm=norm, transform=pc_proj)

    # Add the color bar
    cbar = plt.colorbar(cs, orientation='vertical', label='Precipitation (inches)', ticks=clevs)
    cbar.set_label('Precipitation (inches)', fontsize=18, fontweight='bold')
    cbar.ax.tick_params(labelsize=14, which='both')
    cbar.ax.set_yticklabels(["{:4.2f}".format(i) for i in clevs])

    # Convert time to datetime and format the title
    time_1 = pd.to_datetime(date + ' ' + hour_str + 'Z')
    ax.set_title(f'{time_1.strftime("%Y-%m-%d %H:%M Z")}' + '\n 24hr HREF LPMM [in]', fontsize=16, fontweight='bold')

    # # Add city labels
    # locs = {"Raleigh": [-78.6382, 35.7796], "Fayetteville":[-78.8784, 35.0527], "Goldsboro":[-77.9928, 35.3849], 
    #         "Greensboro": [-79.7920, 36.0726], "Rocky Mount": [-77.7905, 35.9382], "Asheboro": [-79.8136, 35.7079], 
    #         "Rockingham": [-79.7739, 34.9393], "Roxboro": [-78.9828, 36.3938], "Henderson": [-78.3992, 36.3294], 
    #         "Roanoke Rapids": [-77.6541, 36.4615]}

    # for key in locs.keys():
    #     x, y = locs[key]
    #     ax.scatter(x, y, c='k', zorder=1, s=50)
    #     ax.annotate(key, xy=(x, y), color='k', fontsize=8, xytext=(x+0.045, y-0.05), fontweight='bold', 
    #                 bbox=dict(boxstyle="round", fc="black", ec="b", alpha=0.2))

    # Display the plot
    #plt.show()


    # Define the base URL
    base_url = "https://nomads.ncep.noaa.gov/pub/data/nccf/com/href/prod/"

    current_date = datetime.now()
    date = current_date.strftime('%Y%m%d')

    # Specify the date and forecast hour
    #date = "20241108"  # Example date
    hour_str = "00"  # Example hour (can be 00Z or 12Z)
    forecast_hour = "01"  # First forecast hour

    # Define the output directory and file path
    output_dir = "href_data"
    os.makedirs(output_dir, exist_ok=True)

    # Conversion factor from mm to inches
    mm_to_inch = 0.0393701  # 1 mm = 0.0393701 inches

    # Create the plot projection
    pc_proj = ccrs.PlateCarree()

    # Initialize lists to hold hourly precipitation data
    all_hourly_precip = []

    # Function to unpack "Total precipitation" variable based on parameterName
    def unpack_total_precipitation(grib_path):
        try:
            with pygrib.open(grib_path) as grb_file:
                for grb_message in grb_file:
                    if grb_message.parameterName == "Total precipitation":
                        data, lats, lons = grb_message.data()
                        print(f"Extracted Total precipitation from {grib_path}")
                        return data, lats, lons
                print(f"Total precipitation not found in {grib_path}")
        except Exception as e:
            print(f"Error reading {grib_path}: {e}")
        return None, None, None

    # Function to download files with User-Agent header
    def download_file(url, save_path):
        headers = {
            "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/87.0.4280.88 Safari/537.36"
        }
        try:
            response = requests.get(url, headers=headers)
            response.raise_for_status()  # Check if request was successful
            with open(save_path, 'wb') as file:
                file.write(response.content)
            print(f"Downloaded {url}")
        except requests.exceptions.RequestException as e:
            print(f"Failed to download {url}: {e}")

    # Loop through forecast hours (from 01 to 24)
    for i in range(13, 36):
        forecast_hour_str = f"{i:02d}"  # e.g., '01', '02', ..., '24'
        file_url = f"{base_url}/href.{date}/ensprod/href.t{hour_str}z.conus.lpmm.f{forecast_hour_str}.grib2"
        save_path = os.path.join(output_dir, f"href_{date}_t{hour_str}_f{forecast_hour_str}.grib2")
        
        # Download the file for each forecast hour
        download_file(file_url, save_path)

        # Unpack the "Total precipitation" variable from the downloaded file
        hourly_precip, lats, lons = unpack_total_precipitation(save_path)
        
        # If precipitation data is successfully extracted, store it
        if hourly_precip is not None:
            hourly_precip_inch = hourly_precip * mm_to_inch
            all_hourly_precip.append(hourly_precip_inch)

    # Sum the hourly precipitation data for the 24 hours
    total_precip_2 = np.sum(all_hourly_precip, axis=0)

    # Create the figure and axes for plotting
    fig = plt.figure(figsize=(20, 20))
    ax = fig.add_subplot(111, projection=pc_proj)
    ax.set_extent([-74, -84.8, 31, 39])  # Area of interest

    # Draw coastlines, state and country boundaries, and edges of the map
    ax.coastlines(resolution='10m')
    ax.add_feature(cfeature.BORDERS.with_scale('10m'), linewidth=1.5)
    ax.add_feature(cfeature.STATES.with_scale('10m'), linewidth=2.0)
    ax.add_feature(USCOUNTIES.with_scale('500k'), edgecolor='black')

    # Define the contour levels for precipitation
    clevs = [0.00, 0.01, 0.1, 0.25, 0.50, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0,
             6.0, 8.0, 10., 12.0, 15.0, 18.0]

    # Define the color map for the precipitation levels
    cmap_data = [(1.0, 1.0, 1.0),
                 (0.3137255012989044, 0.8156862854957581, 0.8156862854957581),
                 (0.0, 1.0, 1.0),
                 (0.0, 0.8784313797950745, 0.501960813999176),
                 (0.0, 0.7529411911964417, 0.0),
                 (0.501960813999176, 0.8784313797950745, 0.0),
                 (1.0, 1.0, 0.0),
                 (1.0, 0.6274510025978088, 0.0),
                 (1.0, 0.0, 0.0),
                 (1.0, 0.125490203499794, 0.501960813999176),
                 (0.9411764740943909, 0.250980406999588, 1.0),
                 (0.501960813999176, 0.125490203499794, 1.0),
                 (0.250980406999588, 0.250980406999588, 1.0),
                 (0.125490203499794, 0.125490203499794, 0.501960813999176),
                 (0.125490203499794, 0.125490203499794, 0.125490203499794),
                 (0.501960813999176, 0.501960813999176, 0.501960813999176),
                 (0.8784313797950745, 0.8784313797950745, 0.8784313797950745),
                 (0.9333333373069763, 0.8313725590705872, 0.7372549176216125),
                 (0.8549019694328308, 0.6509804129600525, 0.47058823704719543),
                 (0.6274510025978088, 0.42352941632270813, 0.23529411852359772),
                 (0.4000000059604645, 0.20000000298023224, 0.0)]
    cmap = mcolors.ListedColormap(cmap_data, 'precipitation')
    norm = mcolors.BoundaryNorm(clevs, cmap.N)

    # Plot the contour fill
    cs = ax.contourf(lons, lats, total_precip_2, clevs, alpha=0.5, cmap=cmap, norm=norm, transform=pc_proj)

    # Add the color bar
    cbar = plt.colorbar(cs, orientation='vertical', label='Precipitation (inches)', ticks=clevs)
    cbar.set_label('Precipitation (inches)', fontsize=18, fontweight='bold')
    cbar.ax.tick_params(labelsize=14, which='both')
    cbar.ax.set_yticklabels(["{:4.2f}".format(i) for i in clevs])

    # Convert time to datetime and format the title
    time_2 = pd.to_datetime(date + ' ' + hour_str + 'Z')
    ax.set_title(f'{time_2.strftime("%Y-%m-%d %H:%M Z")}' + '\n 24hr HREF LPMM [in]', fontsize=16, fontweight='bold')

    # # Add city labels
    # locs = {"Raleigh": [-78.6382, 35.7796], "Fayetteville":[-78.8784, 35.0527], "Goldsboro":[-77.9928, 35.3849], 
    #         "Greensboro": [-79.7920, 36.0726], "Rocky Mount": [-77.7905, 35.9382], "Asheboro": [-79.8136, 35.7079], 
    #         "Rockingham": [-79.7739, 34.9393], "Roxboro": [-78.9828, 36.3938], "Henderson": [-78.3992, 36.3294], 
    #         "Roanoke Rapids": [-77.6541, 36.4615]}

    # for key in locs.keys():
    #     x, y = locs[key]
    #     ax.scatter(x, y, c='k', zorder=1, s=50)
    #     ax.annotate(key, xy=(x, y), color='k', fontsize=8, xytext=(x+0.045, y-0.05), fontweight='bold', 
    #                 bbox=dict(boxstyle="round", fc="black", ec="b", alpha=0.2))

    # Display the plot
    #plt.show()



    # Define the base URL
    base_url = "https://nomads.ncep.noaa.gov/pub/data/nccf/com/href/prod/"

    current_date = datetime.now()
    #date = current_date.strftime('%Y%m%d')

    previous_date = current_date - timedelta(days=1)
    date = previous_date.strftime('%Y%m%d')

    # Specify the date and forecast hour
    #date = "20241108"  # Example date
    hour_str = "12"  # Example hour (can be 00Z or 12Z)
    forecast_hour = "01"  # First forecast hour

    # Define the output directory and file path
    output_dir = "href_data"
    os.makedirs(output_dir, exist_ok=True)

    # Conversion factor from mm to inches
    mm_to_inch = 0.0393701  # 1 mm = 0.0393701 inches

    # Create the plot projection
    pc_proj = ccrs.PlateCarree()

    # Initialize lists to hold hourly precipitation data
    all_hourly_precip = []

    # Function to unpack "Total precipitation" variable based on parameterName
    def unpack_total_precipitation(grib_path):
        try:
            with pygrib.open(grib_path) as grb_file:
                for grb_message in grb_file:
                    if grb_message.parameterName == "Total precipitation":
                        data, lats, lons = grb_message.data()
                        print(f"Extracted Total precipitation from {grib_path}")
                        return data, lats, lons
                print(f"Total precipitation not found in {grib_path}")
        except Exception as e:
            print(f"Error reading {grib_path}: {e}")
        return None, None, None

    # Function to download files with User-Agent header
    def download_file(url, save_path):
        headers = {
            "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/87.0.4280.88 Safari/537.36"
        }
        try:
            response = requests.get(url, headers=headers)
            response.raise_for_status()  # Check if request was successful
            with open(save_path, 'wb') as file:
                file.write(response.content)
            print(f"Downloaded {url}")
        except requests.exceptions.RequestException as e:
            print(f"Failed to download {url}: {e}")

    # Loop through forecast hours (from 01 to 24)
    for i in range(25, 48):
        forecast_hour_str = f"{i:02d}"  # e.g., '01', '02', ..., '24'
        file_url = f"{base_url}/href.{date}/ensprod/href.t{hour_str}z.conus.lpmm.f{forecast_hour_str}.grib2"
        save_path = os.path.join(output_dir, f"href_{date}_t{hour_str}_f{forecast_hour_str}.grib2")
        
        # Download the file for each forecast hour
        download_file(file_url, save_path)

        # Unpack the "Total precipitation" variable from the downloaded file
        hourly_precip, lats, lons = unpack_total_precipitation(save_path)
        
        # If precipitation data is successfully extracted, store it
        if hourly_precip is not None:
            hourly_precip_inch = hourly_precip * mm_to_inch
            all_hourly_precip.append(hourly_precip_inch)

    # Sum the hourly precipitation data for the 24 hours
    total_precip_3 = np.sum(all_hourly_precip, axis=0)

    # Create the figure and axes for plotting
    fig = plt.figure(figsize=(20, 20))
    ax = fig.add_subplot(111, projection=pc_proj)
    ax.set_extent([-74, -84.8, 31, 39])  # Area of interest

    # Draw coastlines, state and country boundaries, and edges of the map
    ax.coastlines(resolution='10m')
    ax.add_feature(cfeature.BORDERS.with_scale('10m'), linewidth=1.5)
    ax.add_feature(cfeature.STATES.with_scale('10m'), linewidth=2.0)
    ax.add_feature(USCOUNTIES.with_scale('500k'), edgecolor='black')

    # Define the contour levels for precipitation
    clevs = [0.00, 0.01, 0.1, 0.25, 0.50, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0,
             6.0, 8.0, 10., 12.0, 15.0, 18.0]

    # Define the color map for the precipitation levels
    cmap_data = [(1.0, 1.0, 1.0),
                 (0.3137255012989044, 0.8156862854957581, 0.8156862854957581),
                 (0.0, 1.0, 1.0),
                 (0.0, 0.8784313797950745, 0.501960813999176),
                 (0.0, 0.7529411911964417, 0.0),
                 (0.501960813999176, 0.8784313797950745, 0.0),
                 (1.0, 1.0, 0.0),
                 (1.0, 0.6274510025978088, 0.0),
                 (1.0, 0.0, 0.0),
                 (1.0, 0.125490203499794, 0.501960813999176),
                 (0.9411764740943909, 0.250980406999588, 1.0),
                 (0.501960813999176, 0.125490203499794, 1.0),
                 (0.250980406999588, 0.250980406999588, 1.0),
                 (0.125490203499794, 0.125490203499794, 0.501960813999176),
                 (0.125490203499794, 0.125490203499794, 0.125490203499794),
                 (0.501960813999176, 0.501960813999176, 0.501960813999176),
                 (0.8784313797950745, 0.8784313797950745, 0.8784313797950745),
                 (0.9333333373069763, 0.8313725590705872, 0.7372549176216125),
                 (0.8549019694328308, 0.6509804129600525, 0.47058823704719543),
                 (0.6274510025978088, 0.42352941632270813, 0.23529411852359772),
                 (0.4000000059604645, 0.20000000298023224, 0.0)]
    cmap = mcolors.ListedColormap(cmap_data, 'precipitation')
    norm = mcolors.BoundaryNorm(clevs, cmap.N)

    # Plot the contour fill
    cs = ax.contourf(lons, lats, total_precip_3, clevs, alpha=0.5, cmap=cmap, norm=norm, transform=pc_proj)

    # Add the color bar
    cbar = plt.colorbar(cs, orientation='vertical', label='Precipitation (inches)', ticks=clevs)
    cbar.set_label('Precipitation (inches)', fontsize=18, fontweight='bold')
    cbar.ax.tick_params(labelsize=14, which='both')
    cbar.ax.set_yticklabels(["{:4.2f}".format(i) for i in clevs])

    # Convert time to datetime and format the title
    time_3 = pd.to_datetime(date + ' ' + hour_str + 'Z')
    ax.set_title(f'{time_3.strftime("%Y-%m-%d %H:%M Z")}' + '\n 24hr HREF LPMM [in]', fontsize=16, fontweight='bold')

    # # Add city labels
    # locs = {"Raleigh": [-78.6382, 35.7796], "Fayetteville":[-78.8784, 35.0527], "Goldsboro":[-77.9928, 35.3849], 
    #         "Greensboro": [-79.7920, 36.0726], "Rocky Mount": [-77.7905, 35.9382], "Asheboro": [-79.8136, 35.7079], 
    #         "Rockingham": [-79.7739, 34.9393], "Roxboro": [-78.9828, 36.3938], "Henderson": [-78.3992, 36.3294], 
    #         "Roanoke Rapids": [-77.6541, 36.4615]}

    # for key in locs.keys():
    #     x, y = locs[key]
    #     ax.scatter(x, y, c='k', zorder=1, s=50)
    #     ax.annotate(key, xy=(x, y), color='k', fontsize=8, xytext=(x+0.045, y-0.05), fontweight='bold', 
    #                 bbox=dict(boxstyle="round", fc="black", ec="b", alpha=0.2))

    # Display the plot
    #plt.show()


    #If you want to remove the grib files use these lines:

    # import os
    # import glob

    # # Specify the directory where the .grib2 files are stored
    # grib_directory = "href_data"

    # # Function to delete all .grib2 files in the specified directory
    # def clean_up_grib_files(directory):
    #     grib_files = glob.glob(os.path.join(directory, "*.grib2"))
    #     for file_path in grib_files:
    #         try:
    #             os.remove(file_path)
    #             print(f"Deleted {file_path}")
    #         except OSError as e:
    #             print(f"Error deleting {file_path}: {e}")

    # # Run the cleanup function at the end of the script
    # clean_up_grib_files(grib_directory)




    '''now we want to combine the three plots side by side'''


    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    import matplotlib.colors as mcolors
    import numpy as np

    # Create the figure and axes for plotting
    fig, ax = plt.subplots(1, 3, figsize=(18, 8), subplot_kw={'projection': ccrs.PlateCarree()})
    #fig.subplots_adjust(wspace=0.05)
    fig.subplots_adjust(top=1.5, bottom=0.2, wspace=0.04, hspace=0.9,left=0.05, right=0.95)
    #fig.tight_layout(pad=8.5, w_pad=0.5, h_pad=2.0)

    # Set up the contour levels and color map
    clevs = [0.00, 0.01, 0.1, 0.25, 0.50, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0,
             6.0, 8.0, 10., 12.0, 15.0, 18.0]
    cmap_data = [
        (1.0, 1.0, 1.0), (0.3137255012989044, 0.8156862854957581, 0.8156862854957581),
        (0.0, 1.0, 1.0), (0.0, 0.8784313797950745, 0.501960813999176),
        (0.0, 0.7529411911964417, 0.0), (0.501960813999176, 0.8784313797950745, 0.0),
        (1.0, 1.0, 0.0), (1.0, 0.6274510025978088, 0.0), (1.0, 0.0, 0.0),
        (1.0, 0.125490203499794, 0.501960813999176), (0.9411764740943909, 0.250980406999588, 1.0),
        (0.501960813999176, 0.125490203499794, 1.0), (0.250980406999588, 0.250980406999588, 1.0),
        (0.125490203499794, 0.125490203499794, 0.501960813999176), (0.125490203499794, 0.125490203499794, 0.125490203499794),
        (0.501960813999176, 0.501960813999176, 0.501960813999176), (0.8784313797950745, 0.8784313797950745, 0.8784313797950745),
        (0.9333333373069763, 0.8313725590705872, 0.7372549176216125), (0.8549019694328308, 0.6509804129600525, 0.47058823704719543),
        (0.6274510025978088, 0.42352941632270813, 0.23529411852359772), (0.4000000059604645, 0.20000000298023224, 0.0)
    ]
    cmap = mcolors.ListedColormap(cmap_data, 'precipitation')
    norm = mcolors.BoundaryNorm(clevs, cmap.N)

    # Sample data arrays (replace with your own)
    #lons, lats = np.meshgrid(np.linspace(-84.8, -74, 1473), np.linspace(31, 39, 1473))
    total_precip_1 = total_precip_1
    total_precip_2 = total_precip_2
    total_precip_3 = total_precip_3

    # Plot each panel
    # Inside the loop for each subplot, comment out or modify the state feature layer
    for i, data in enumerate([total_precip_1, total_precip_2, total_precip_3]):
        cs = ax[i].contourf(lons, lats, data, clevs, alpha=0.5, cmap=cmap, norm=norm, transform=ccrs.PlateCarree())
        ax[i].coastlines(resolution='10m')
        ax[i].add_feature(cfeature.BORDERS.with_scale('10m'), linewidth=1.5)
        ax[i].add_feature(cfeature.STATES.with_scale('10m'), linewidth=1.0)  # Reduced width

        # Optional: Comment out or remove the following line if it’s adding unwanted artifacts
        # ax[i].add_feature(cfeature.NaturalEarthFeature('cultural', 'admin_1_states_provinces_lines', '10m'),
        #                   edgecolor='black', linewidth=0.5)
        ax[i].add_feature(USCOUNTIES.with_scale('500k'), edgecolor='gray', linewidth =0.7)
        ax[i].set_extent([-84.8, -74, 31, 39])
        
        #ax[i].set_title(f'Panel {i + 1}', fontsize=14, fontweight='bold')
        for i, time in enumerate([time_1, time_2, time_3]):
            ax[i].set_title(f'{time.strftime("%Y-%m-%d %H:%M Z")}' + '\n 24hr HREF LPMM [in]')
    # Shared colorbar code remains the same
    # Add overall title

    current_date = datetime.now()
    date_1 = current_date.strftime('%Y%m%d')

    previous_date = current_date + timedelta(days=1)
    date_2 = previous_date.strftime('%Y%m%d')



    # Set the suptitle with proper string concatenation
    fig.suptitle('24hr HREF LPMM [in] dprog/dt\nValid: 12Z ' + date_1 + ' to 12Z ' + date_2, 
                 fontsize=18, fontweight='bold', y=0.88)

    #if you wanna mess with a smaller font subtitle:
        
    # # Set the main title (suptitle)
    # fig.suptitle('24hr HREF LPMM [in] dprog/dt', fontsize=18, fontweight='bold', y=0.9)

    # # Add the subtitle with different styling
    # fig.text(0.5, 0.87, 'Valid: 12Z ' + date_1 + ' to 12Z ' + date_2,
    #          ha='center', fontsize=12, fontweight='normal')  # Adjust fontsize and fontweight

    # Add a single color bar below the plots
    cbar = fig.colorbar(cs, ax=ax, orientation='horizontal', pad=0.05, fraction=0.02, aspect=50)
    cbar.set_label('Precipitation (inches)', fontsize=16, fontweight='bold')
    cbar.ax.tick_params(labelsize=12)

    #plt.show()

    #save = '/Users/nicholasluchetti/Documents/HREF_Plots/'
    #fig.savefig(save + 'HREF_LPMM_RUN_COMPARE', dpi=300, bbox_inches='tight')
   # fig.savefig(save + 'HREF_LPMM_RUN_COMPARE', dpi=300)

    # Create a folder for the output if it doesn't exist
    output_folder = "href_data"
    if not os.path.exists(output_folder):
    os.makedirs(output_folder)

   # Save to the local folder instead of your Documents folder
   fig.savefig(f'{output_folder}/HREF_LPMM_RUN_COMPARE.png', dpi=300)


    '''Now we want to make a 4 panel plot where we mask 
    different thresholds (>3, >6, >9, >12 inches) and plot them 
    like a paintball plot'''

    #First start by masking out the different thresholds: 
        
    ## >3 inches:
        
    three_inch_run_1 = np.ma.masked_less(total_precip_1, 3)
    three_inch_run_2 = np.ma.masked_less(total_precip_2, 3)
    three_inch_run_3 = np.ma.masked_less(total_precip_3, 3)


    #> 6 inches
    six_inch_run_1 = np.ma.masked_less(total_precip_1, 6)
    six_inch_run_2 = np.ma.masked_less(total_precip_2, 6)
    six_inch_run_3 = np.ma.masked_less(total_precip_3, 6)

    #> 9 inches
    nine_inch_run_1 = np.ma.masked_less(total_precip_1, 9)
    nine_inch_run_2 = np.ma.masked_less(total_precip_2, 9)
    nine_inch_run_3 = np.ma.masked_less(total_precip_3, 9)

    #> 12 inches
    twelve_inch_run_1 = np.ma.masked_less(total_precip_1, 12)
    twelve_inch_run_2 = np.ma.masked_less(total_precip_2, 12)
    twelve_inch_run_3 = np.ma.masked_less(total_precip_3, 12)


    #%%

    import numpy as np
    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    import matplotlib.colors as mcolors
    from matplotlib.lines import Line2D  # For custom legends

    # Create a 2x2 subplot grid
    fig, ax = plt.subplots(2, 2, figsize=(15, 12), subplot_kw={'projection': ccrs.PlateCarree()})
    #fig.subplots_adjust(top=1.5, bottom=0.5, wspace=0.04, hspace=0.9,left=0.05, right=0.95)

    # Define the shades of blue for the three simulations (darkest to lightest)
    blue_shades = ['#00008B', '#4169E1', '#87CEFA']  # Dark blue, medium blue, light blue
    alpha_values = [1.0, 0.6, 0.4]  # Transparency for each simulation: simulation 1 = opaque, others faded

    # Define a list of all thresholds and corresponding masked arrays
    thresholds = [3, 6, 9, 12]
    thresholded_data = [
        [three_inch_run_1, three_inch_run_2, three_inch_run_3],  # > 3 inches
        [six_inch_run_1, six_inch_run_2, six_inch_run_3],        # > 6 inches
        [nine_inch_run_1, nine_inch_run_2, nine_inch_run_3],     # > 9 inches
        [twelve_inch_run_1, twelve_inch_run_2, twelve_inch_run_3] # > 12 inches
    ]

    # Titles for the subplots
    titles = ['> 3 inches', '> 6 inches', '> 9 inches', '> 12 inches']

    # Loop through each subplot and plot the corresponding masked data
    for i, (thresh, data) in enumerate(zip(thresholds, thresholded_data)):
        row, col = divmod(i, 2)  # Determine the subplot position (2x2 grid)
        
        # For each simulation, apply the mask and plot with different shades of blue
        for j, (precip_data, shade, alpha) in enumerate(zip(data, blue_shades, alpha_values)):
            # Apply the mask for values less than the threshold
            masked_data = np.ma.masked_less(precip_data, thresh)

            # Plot only the masked data (values above the threshold), applying transparency
            cs = ax[row, col].contourf(lons, lats, masked_data, cmap=mcolors.ListedColormap([shade]),
                                       levels=[0, np.nanmax(masked_data)], transform=ccrs.PlateCarree(), alpha=alpha)

        # Add coastlines and features
        ax[row, col].coastlines(resolution='10m')
        ax[row, col].add_feature(cfeature.BORDERS.with_scale('10m'), linewidth=1.5)
        ax[row, col].add_feature(cfeature.STATES.with_scale('10m'), linewidth=1.0)
        
        # Set extent (optional based on your data)
        ax[row, col].set_extent([-84.8, -74, 31, 39])
       
        
        
        # Set the title for each subplot
        ax[row, col].set_title(titles[i], fontsize=14, fontweight='bold')

        # Create a custom legend with patches for each simulation
        legend_labels = [f'{time_1.strftime("%Y-%m-%d %H:%M Z")}', f'{time_2.strftime("%Y-%m-%d %H:%M Z")}', f'{time_3.strftime("%Y-%m-%d %H:%M Z")}']
        legend_colors = blue_shades  # Corresponding to each simulation
        legend_alphas = alpha_values  # Alpha values for transparency
        legend_patches = [Line2D([0], [0], marker='o', color='w', markerfacecolor=shade, markersize=10, alpha=alpha,
                                 label=label) for shade, alpha, label in zip(legend_colors, legend_alphas, legend_labels)]
        
        # Add the legend to the upper right of the current subplot
        ax[row, col].legend(handles=legend_patches, loc='lower right', fontsize=8, title='HREF Run', title_fontsize=10)

    # Set the suptitle with proper string concatenation
    fig.suptitle('24hr HREF LPMM [in] \nValid: 12Z ' + date_1 + ' to 12Z ' + date_2, 
                 fontsize=14, fontweight='bold', y=0.98)

    # Add a color bar (optional, but useful for understanding the plot)
    #cbar = fig.colorbar(cs, ax=ax.ravel().tolist(), orientation='horizontal', pad=0.05, fraction=0.02, aspect=50)
    #cbar.set_label('Precipitation (inches)', fontsize=12)
    plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1, wspace=0.05, hspace=0.25)

    # Show the plot
    plt.tight_layout()
    #plt.show()

    #save = '/Users/nicholasluchetti/Documents/HREF_Plots/'
    #fig.savefig(save + 'HREF_LPMM_THRESHOLD_COMPARE', dpi=300, bbox_inches='tight')
# Create a folder for the output if it doesn't exist
    output_folder = "href_data"
    if not os.path.exists(output_folder):
    os.makedirs(output_folder)

   # Save to the local folder instead of your Documents folder
   fig.savefig(f'{output_folder}/HREF_LPMM_THRESHOLD_COMPARE.png', dpi=300)


elif 00 <= current_hour <= 13:

  
  # Define the base URL
  base_url = "https://nomads.ncep.noaa.gov/pub/data/nccf/com/href/prod/"

  current_date = datetime.now()
  date = current_date.strftime('%Y%m%d')

  # Specify the date and forecast hour
  #date = "20241108"  # Example date
  hour_str = "00"  # Example hour (can be 00Z or 12Z)
  forecast_hour = "01"  # First forecast hour

  # Define the output directory and file path
  output_dir = "href_data"
  os.makedirs(output_dir, exist_ok=True)

  # Conversion factor from mm to inches
  mm_to_inch = 0.0393701  # 1 mm = 0.0393701 inches

  # Create the plot projection
  pc_proj = ccrs.PlateCarree()

  # Initialize lists to hold hourly precipitation data
  all_hourly_precip = []

  # Function to unpack "Total precipitation" variable based on parameterName
  def unpack_total_precipitation(grib_path):
      try:
          with pygrib.open(grib_path) as grb_file:
              for grb_message in grb_file:
                  if grb_message.parameterName == "Total precipitation":
                      data, lats, lons = grb_message.data()
                      print(f"Extracted Total precipitation from {grib_path}")
                      return data, lats, lons
              print(f"Total precipitation not found in {grib_path}")
      except Exception as e:
          print(f"Error reading {grib_path}: {e}")
      return None, None, None

  # Function to download files with User-Agent header
  def download_file(url, save_path):
      headers = {
          "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/87.0.4280.88 Safari/537.36"
      }
      try:
          response = requests.get(url, headers=headers)
          response.raise_for_status()  # Check if request was successful
          with open(save_path, 'wb') as file:
              file.write(response.content)
          print(f"Downloaded {url}")
      except requests.exceptions.RequestException as e:
          print(f"Failed to download {url}: {e}")

  # Loop through forecast hours (from 01 to 24)
  for i in range(1, 24):
      forecast_hour_str = f"{i:02d}"  # e.g., '01', '02', ..., '24'
      file_url = f"{base_url}/href.{date}/ensprod/href.t{hour_str}z.conus.lpmm.f{forecast_hour_str}.grib2"
      save_path = os.path.join(output_dir, f"href_{date}_t{hour_str}_f{forecast_hour_str}.grib2")
      
      # Download the file for each forecast hour
      download_file(file_url, save_path)

      # Unpack the "Total precipitation" variable from the downloaded file
      hourly_precip, lats, lons = unpack_total_precipitation(save_path)
      
      # If precipitation data is successfully extracted, store it
      if hourly_precip is not None:
          hourly_precip_inch = hourly_precip * mm_to_inch
          all_hourly_precip.append(hourly_precip_inch)

  # Sum the hourly precipitation data for the 24 hours
  total_precip_1 = np.sum(all_hourly_precip, axis=0)

  # Create the figure and axes for plotting
  fig = plt.figure(figsize=(20, 20))
  ax = fig.add_subplot(111, projection=pc_proj)
  ax.set_extent([-74, -84.8, 31, 39])  # Area of interest

  # Draw coastlines, state and country boundaries, and edges of the map
  ax.coastlines(resolution='10m')
  ax.add_feature(cfeature.BORDERS.with_scale('10m'), linewidth=1.5)
  ax.add_feature(cfeature.STATES.with_scale('10m'), linewidth=2.0)
  ax.add_feature(USCOUNTIES.with_scale('500k'), edgecolor='black')

  # Define the contour levels for precipitation
  clevs = [0.00, 0.01, 0.1, 0.25, 0.50, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0,
           6.0, 8.0, 10., 12.0, 15.0, 18.0]

  # Define the color map for the precipitation levels
  cmap_data = [(1.0, 1.0, 1.0),
               (0.3137255012989044, 0.8156862854957581, 0.8156862854957581),
               (0.0, 1.0, 1.0),
               (0.0, 0.8784313797950745, 0.501960813999176),
               (0.0, 0.7529411911964417, 0.0),
               (0.501960813999176, 0.8784313797950745, 0.0),
               (1.0, 1.0, 0.0),
               (1.0, 0.6274510025978088, 0.0),
               (1.0, 0.0, 0.0),
               (1.0, 0.125490203499794, 0.501960813999176),
               (0.9411764740943909, 0.250980406999588, 1.0),
               (0.501960813999176, 0.125490203499794, 1.0),
               (0.250980406999588, 0.250980406999588, 1.0),
               (0.125490203499794, 0.125490203499794, 0.501960813999176),
               (0.125490203499794, 0.125490203499794, 0.125490203499794),
               (0.501960813999176, 0.501960813999176, 0.501960813999176),
               (0.8784313797950745, 0.8784313797950745, 0.8784313797950745),
               (0.9333333373069763, 0.8313725590705872, 0.7372549176216125),
               (0.8549019694328308, 0.6509804129600525, 0.47058823704719543),
               (0.6274510025978088, 0.42352941632270813, 0.23529411852359772),
               (0.4000000059604645, 0.20000000298023224, 0.0)]
  cmap = mcolors.ListedColormap(cmap_data, 'precipitation')
  norm = mcolors.BoundaryNorm(clevs, cmap.N)

  # Plot the contour fill
  cs = ax.contourf(lons, lats, total_precip_1, clevs, alpha=0.5, cmap=cmap, norm=norm, transform=pc_proj)

  # Add the color bar
  cbar = plt.colorbar(cs, orientation='vertical', label='Precipitation (inches)', ticks=clevs)
  cbar.set_label('Precipitation (inches)', fontsize=18, fontweight='bold')
  cbar.ax.tick_params(labelsize=14, which='both')
  cbar.ax.set_yticklabels(["{:4.2f}".format(i) for i in clevs])

  # Convert time to datetime and format the title
  time_1 = pd.to_datetime(date + ' ' + hour_str + 'Z')
  ax.set_title(f'{time_1.strftime("%Y-%m-%d %H:%M Z")}' + '\n 24hr HREF LPMM [in]', fontsize=16, fontweight='bold')

  # # Add city labels
  # locs = {"Raleigh": [-78.6382, 35.7796], "Fayetteville":[-78.8784, 35.0527], "Goldsboro":[-77.9928, 35.3849], 
  #         "Greensboro": [-79.7920, 36.0726], "Rocky Mount": [-77.7905, 35.9382], "Asheboro": [-79.8136, 35.7079], 
  #         "Rockingham": [-79.7739, 34.9393], "Roxboro": [-78.9828, 36.3938], "Henderson": [-78.3992, 36.3294], 
  #         "Roanoke Rapids": [-77.6541, 36.4615]}

  # for key in locs.keys():
  #     x, y = locs[key]
  #     ax.scatter(x, y, c='k', zorder=1, s=50)
  #     ax.annotate(key, xy=(x, y), color='k', fontsize=8, xytext=(x+0.045, y-0.05), fontweight='bold', 
  #                 bbox=dict(boxstyle="round", fc="black", ec="b", alpha=0.2))

  # Display the plot
  #plt.show()

  # Define the base URL
  base_url = "https://nomads.ncep.noaa.gov/pub/data/nccf/com/href/prod/"

  current_date = datetime.now()
  #date = current_date.strftime('%Y%m%d')
  
  previous_date = current_date - timedelta(days=1)
  date = previous_date.strftime('%Y%m%d')
 
  # Specify the date and forecast hour
  #date = "20241108"  # Example date
  hour_str = "12"  # Example hour (can be 00Z or 12Z)
  forecast_hour = "01"  # First forecast hour

  # Define the output directory and file path
  output_dir = "href_data"
  os.makedirs(output_dir, exist_ok=True)

  # Conversion factor from mm to inches
  mm_to_inch = 0.0393701  # 1 mm = 0.0393701 inches

  # Create the plot projection
  pc_proj = ccrs.PlateCarree()

  # Initialize lists to hold hourly precipitation data
  all_hourly_precip = []

  # Function to unpack "Total precipitation" variable based on parameterName
  def unpack_total_precipitation(grib_path):
      try:
          with pygrib.open(grib_path) as grb_file:
              for grb_message in grb_file:
                  if grb_message.parameterName == "Total precipitation":
                      data, lats, lons = grb_message.data()
                      print(f"Extracted Total precipitation from {grib_path}")
                      return data, lats, lons
              print(f"Total precipitation not found in {grib_path}")
      except Exception as e:
          print(f"Error reading {grib_path}: {e}")
      return None, None, None

  # Function to download files with User-Agent header
  def download_file(url, save_path):
      headers = {
          "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/87.0.4280.88 Safari/537.36"
      }
      try:
          response = requests.get(url, headers=headers)
          response.raise_for_status()  # Check if request was successful
          with open(save_path, 'wb') as file:
              file.write(response.content)
          print(f"Downloaded {url}")
      except requests.exceptions.RequestException as e:
          print(f"Failed to download {url}: {e}")

  # Loop through forecast hours (from 01 to 24)
  for i in range(13, 36):
      forecast_hour_str = f"{i:02d}"  # e.g., '01', '02', ..., '24'
      file_url = f"{base_url}/href.{date}/ensprod/href.t{hour_str}z.conus.lpmm.f{forecast_hour_str}.grib2"
      save_path = os.path.join(output_dir, f"href_{date}_t{hour_str}_f{forecast_hour_str}.grib2")
      
      # Download the file for each forecast hour
      download_file(file_url, save_path)

      # Unpack the "Total precipitation" variable from the downloaded file
      hourly_precip, lats, lons = unpack_total_precipitation(save_path)
      
      # If precipitation data is successfully extracted, store it
      if hourly_precip is not None:
          hourly_precip_inch = hourly_precip * mm_to_inch
          all_hourly_precip.append(hourly_precip_inch)

  # Sum the hourly precipitation data for the 24 hours
  total_precip_2 = np.sum(all_hourly_precip, axis=0)

  # Create the figure and axes for plotting
  fig = plt.figure(figsize=(20, 20))
  ax = fig.add_subplot(111, projection=pc_proj)
  ax.set_extent([-74, -84.8, 31, 39])  # Area of interest

  # Draw coastlines, state and country boundaries, and edges of the map
  ax.coastlines(resolution='10m')
  ax.add_feature(cfeature.BORDERS.with_scale('10m'), linewidth=1.5)
  ax.add_feature(cfeature.STATES.with_scale('10m'), linewidth=2.0)
  ax.add_feature(USCOUNTIES.with_scale('500k'), edgecolor='black')

  # Define the contour levels for precipitation
  clevs = [0.00, 0.01, 0.1, 0.25, 0.50, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0,
           6.0, 8.0, 10., 12.0, 15.0, 18.0]

  # Define the color map for the precipitation levels
  cmap_data = [(1.0, 1.0, 1.0),
               (0.3137255012989044, 0.8156862854957581, 0.8156862854957581),
               (0.0, 1.0, 1.0),
               (0.0, 0.8784313797950745, 0.501960813999176),
               (0.0, 0.7529411911964417, 0.0),
               (0.501960813999176, 0.8784313797950745, 0.0),
               (1.0, 1.0, 0.0),
               (1.0, 0.6274510025978088, 0.0),
               (1.0, 0.0, 0.0),
               (1.0, 0.125490203499794, 0.501960813999176),
               (0.9411764740943909, 0.250980406999588, 1.0),
               (0.501960813999176, 0.125490203499794, 1.0),
               (0.250980406999588, 0.250980406999588, 1.0),
               (0.125490203499794, 0.125490203499794, 0.501960813999176),
               (0.125490203499794, 0.125490203499794, 0.125490203499794),
               (0.501960813999176, 0.501960813999176, 0.501960813999176),
               (0.8784313797950745, 0.8784313797950745, 0.8784313797950745),
               (0.9333333373069763, 0.8313725590705872, 0.7372549176216125),
               (0.8549019694328308, 0.6509804129600525, 0.47058823704719543),
               (0.6274510025978088, 0.42352941632270813, 0.23529411852359772),
               (0.4000000059604645, 0.20000000298023224, 0.0)]
  cmap = mcolors.ListedColormap(cmap_data, 'precipitation')
  norm = mcolors.BoundaryNorm(clevs, cmap.N)

  # Plot the contour fill
  cs = ax.contourf(lons, lats, total_precip_2, clevs, alpha=0.5, cmap=cmap, norm=norm, transform=pc_proj)

  # Add the color bar
  cbar = plt.colorbar(cs, orientation='vertical', label='Precipitation (inches)', ticks=clevs)
  cbar.set_label('Precipitation (inches)', fontsize=18, fontweight='bold')
  cbar.ax.tick_params(labelsize=14, which='both')
  cbar.ax.set_yticklabels(["{:4.2f}".format(i) for i in clevs])

  # Convert time to datetime and format the title
  time_2 = pd.to_datetime(date + ' ' + hour_str + 'Z')
  ax.set_title(f'{time_2.strftime("%Y-%m-%d %H:%M Z")}' + '\n 24hr HREF LPMM [in]', fontsize=16, fontweight='bold')

  # # Add city labels
  # locs = {"Raleigh": [-78.6382, 35.7796], "Fayetteville":[-78.8784, 35.0527], "Goldsboro":[-77.9928, 35.3849], 
  #         "Greensboro": [-79.7920, 36.0726], "Rocky Mount": [-77.7905, 35.9382], "Asheboro": [-79.8136, 35.7079], 
  #         "Rockingham": [-79.7739, 34.9393], "Roxboro": [-78.9828, 36.3938], "Henderson": [-78.3992, 36.3294], 
  #         "Roanoke Rapids": [-77.6541, 36.4615]}

  # for key in locs.keys():
  #     x, y = locs[key]
  #     ax.scatter(x, y, c='k', zorder=1, s=50)
  #     ax.annotate(key, xy=(x, y), color='k', fontsize=8, xytext=(x+0.045, y-0.05), fontweight='bold', 
  #                 bbox=dict(boxstyle="round", fc="black", ec="b", alpha=0.2))

  # Display the plot
  #plt.show()


  # Define the base URL
  base_url = "https://nomads.ncep.noaa.gov/pub/data/nccf/com/href/prod/"

  current_date = datetime.now()
  #date = current_date.strftime('%Y%m%d')

  previous_date = current_date - timedelta(days=1)
  date = previous_date.strftime('%Y%m%d')

  # Specify the date and forecast hour
  #date = "20241108"  # Example date
  hour_str = "00"  # Example hour (can be 00Z or 12Z)
  forecast_hour = "01"  # First forecast hour

  # Define the output directory and file path
  output_dir = "href_data"
  os.makedirs(output_dir, exist_ok=True)

  # Conversion factor from mm to inches
  mm_to_inch = 0.0393701  # 1 mm = 0.0393701 inches

  # Create the plot projection
  pc_proj = ccrs.PlateCarree()

  # Initialize lists to hold hourly precipitation data
  all_hourly_precip = []

  # Function to unpack "Total precipitation" variable based on parameterName
  def unpack_total_precipitation(grib_path):
      try:
          with pygrib.open(grib_path) as grb_file:
              for grb_message in grb_file:
                  if grb_message.parameterName == "Total precipitation":
                      data, lats, lons = grb_message.data()
                      print(f"Extracted Total precipitation from {grib_path}")
                      return data, lats, lons
              print(f"Total precipitation not found in {grib_path}")
      except Exception as e:
          print(f"Error reading {grib_path}: {e}")
      return None, None, None

  # Function to download files with User-Agent header
  def download_file(url, save_path):
      headers = {
          "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/87.0.4280.88 Safari/537.36"
      }
      try:
          response = requests.get(url, headers=headers)
          response.raise_for_status()  # Check if request was successful
          with open(save_path, 'wb') as file:
              file.write(response.content)
          print(f"Downloaded {url}")
      except requests.exceptions.RequestException as e:
          print(f"Failed to download {url}: {e}")

  # Loop through forecast hours (from 01 to 24)
  for i in range(25, 48):
      forecast_hour_str = f"{i:02d}"  # e.g., '01', '02', ..., '24'
      file_url = f"{base_url}/href.{date}/ensprod/href.t{hour_str}z.conus.lpmm.f{forecast_hour_str}.grib2"
      save_path = os.path.join(output_dir, f"href_{date}_t{hour_str}_f{forecast_hour_str}.grib2")
      
      # Download the file for each forecast hour
      download_file(file_url, save_path)

      # Unpack the "Total precipitation" variable from the downloaded file
      hourly_precip, lats, lons = unpack_total_precipitation(save_path)
      
      # If precipitation data is successfully extracted, store it
      if hourly_precip is not None:
          hourly_precip_inch = hourly_precip * mm_to_inch
          all_hourly_precip.append(hourly_precip_inch)

  # Sum the hourly precipitation data for the 24 hours
  total_precip_3 = np.sum(all_hourly_precip, axis=0)

  # Create the figure and axes for plotting
  fig = plt.figure(figsize=(20, 20))
  ax = fig.add_subplot(111, projection=pc_proj)
  ax.set_extent([-74, -84.8, 31, 39])  # Area of interest

  # Draw coastlines, state and country boundaries, and edges of the map
  ax.coastlines(resolution='10m')
  ax.add_feature(cfeature.BORDERS.with_scale('10m'), linewidth=1.5)
  ax.add_feature(cfeature.STATES.with_scale('10m'), linewidth=2.0)
  ax.add_feature(USCOUNTIES.with_scale('500k'), edgecolor='black')

  # Define the contour levels for precipitation
  clevs = [0.00, 0.01, 0.1, 0.25, 0.50, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0,
           6.0, 8.0, 10., 12.0, 15.0, 18.0]

  # Define the color map for the precipitation levels
  cmap_data = [(1.0, 1.0, 1.0),
               (0.3137255012989044, 0.8156862854957581, 0.8156862854957581),
               (0.0, 1.0, 1.0),
               (0.0, 0.8784313797950745, 0.501960813999176),
               (0.0, 0.7529411911964417, 0.0),
               (0.501960813999176, 0.8784313797950745, 0.0),
               (1.0, 1.0, 0.0),
               (1.0, 0.6274510025978088, 0.0),
               (1.0, 0.0, 0.0),
               (1.0, 0.125490203499794, 0.501960813999176),
               (0.9411764740943909, 0.250980406999588, 1.0),
               (0.501960813999176, 0.125490203499794, 1.0),
               (0.250980406999588, 0.250980406999588, 1.0),
               (0.125490203499794, 0.125490203499794, 0.501960813999176),
               (0.125490203499794, 0.125490203499794, 0.125490203499794),
               (0.501960813999176, 0.501960813999176, 0.501960813999176),
               (0.8784313797950745, 0.8784313797950745, 0.8784313797950745),
               (0.9333333373069763, 0.8313725590705872, 0.7372549176216125),
               (0.8549019694328308, 0.6509804129600525, 0.47058823704719543),
               (0.6274510025978088, 0.42352941632270813, 0.23529411852359772),
               (0.4000000059604645, 0.20000000298023224, 0.0)]
  cmap = mcolors.ListedColormap(cmap_data, 'precipitation')
  norm = mcolors.BoundaryNorm(clevs, cmap.N)

  # Plot the contour fill
  cs = ax.contourf(lons, lats, total_precip_3, clevs, alpha=0.5, cmap=cmap, norm=norm, transform=pc_proj)

  # Add the color bar
  cbar = plt.colorbar(cs, orientation='vertical', label='Precipitation (inches)', ticks=clevs)
  cbar.set_label('Precipitation (inches)', fontsize=18, fontweight='bold')
  cbar.ax.tick_params(labelsize=14, which='both')
  cbar.ax.set_yticklabels(["{:4.2f}".format(i) for i in clevs])

  # Convert time to datetime and format the title
  time_3 = pd.to_datetime(date + ' ' + hour_str + 'Z')
  ax.set_title(f'{time_3.strftime("%Y-%m-%d %H:%M Z")}' + '\n 24hr HREF LPMM [in]', fontsize=16, fontweight='bold')

  # # Add city labels
  # locs = {"Raleigh": [-78.6382, 35.7796], "Fayetteville":[-78.8784, 35.0527], "Goldsboro":[-77.9928, 35.3849], 
  #         "Greensboro": [-79.7920, 36.0726], "Rocky Mount": [-77.7905, 35.9382], "Asheboro": [-79.8136, 35.7079], 
  #         "Rockingham": [-79.7739, 34.9393], "Roxboro": [-78.9828, 36.3938], "Henderson": [-78.3992, 36.3294], 
  #         "Roanoke Rapids": [-77.6541, 36.4615]}

  # for key in locs.keys():
  #     x, y = locs[key]
  #     ax.scatter(x, y, c='k', zorder=1, s=50)
  #     ax.annotate(key, xy=(x, y), color='k', fontsize=8, xytext=(x+0.045, y-0.05), fontweight='bold', 
  #                 bbox=dict(boxstyle="round", fc="black", ec="b", alpha=0.2))

  # Display the plot
  #plt.show()


  #If you want to remove the grib files use these lines:

  # import os
  # import glob

  # # Specify the directory where the .grib2 files are stored
  # grib_directory = "href_data"

  # # Function to delete all .grib2 files in the specified directory
  # def clean_up_grib_files(directory):
  #     grib_files = glob.glob(os.path.join(directory, "*.grib2"))
  #     for file_path in grib_files:
  #         try:
  #             os.remove(file_path)
  #             print(f"Deleted {file_path}")
  #         except OSError as e:
  #             print(f"Error deleting {file_path}: {e}")

  # # Run the cleanup function at the end of the script
  # clean_up_grib_files(grib_directory)



  '''now we want to combine the three plots side by side'''


  import matplotlib.pyplot as plt
  import cartopy.crs as ccrs
  import cartopy.feature as cfeature
  import matplotlib.colors as mcolors
  import numpy as np

  # Create the figure and axes for plotting
  fig, ax = plt.subplots(1, 3, figsize=(18, 8), subplot_kw={'projection': ccrs.PlateCarree()})
  #fig.subplots_adjust(wspace=0.05)
  fig.subplots_adjust(top=1.5, bottom=0.2, wspace=0.04, hspace=0.9,left=0.05, right=0.95)
  #fig.tight_layout(pad=8.5, w_pad=0.5, h_pad=2.0)

  # Set up the contour levels and color map
  clevs = [0.00, 0.01, 0.1, 0.25, 0.50, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0,
           6.0, 8.0, 10., 12.0, 15.0, 18.0]
  cmap_data = [
      (1.0, 1.0, 1.0), (0.3137255012989044, 0.8156862854957581, 0.8156862854957581),
      (0.0, 1.0, 1.0), (0.0, 0.8784313797950745, 0.501960813999176),
      (0.0, 0.7529411911964417, 0.0), (0.501960813999176, 0.8784313797950745, 0.0),
      (1.0, 1.0, 0.0), (1.0, 0.6274510025978088, 0.0), (1.0, 0.0, 0.0),
      (1.0, 0.125490203499794, 0.501960813999176), (0.9411764740943909, 0.250980406999588, 1.0),
      (0.501960813999176, 0.125490203499794, 1.0), (0.250980406999588, 0.250980406999588, 1.0),
      (0.125490203499794, 0.125490203499794, 0.501960813999176), (0.125490203499794, 0.125490203499794, 0.125490203499794),
      (0.501960813999176, 0.501960813999176, 0.501960813999176), (0.8784313797950745, 0.8784313797950745, 0.8784313797950745),
      (0.9333333373069763, 0.8313725590705872, 0.7372549176216125), (0.8549019694328308, 0.6509804129600525, 0.47058823704719543),
      (0.6274510025978088, 0.42352941632270813, 0.23529411852359772), (0.4000000059604645, 0.20000000298023224, 0.0)
  ]
  cmap = mcolors.ListedColormap(cmap_data, 'precipitation')
  norm = mcolors.BoundaryNorm(clevs, cmap.N)

  # Sample data arrays (replace with your own)
  #lons, lats = np.meshgrid(np.linspace(-84.8, -74, 1473), np.linspace(31, 39, 1473))
  total_precip_1 = total_precip_1
  total_precip_2 = total_precip_2
  total_precip_3 = total_precip_3

  # Plot each panel
  # Inside the loop for each subplot, comment out or modify the state feature layer
  for i, data in enumerate([total_precip_1, total_precip_2, total_precip_3]):
      cs = ax[i].contourf(lons, lats, data, clevs, alpha=0.5, cmap=cmap, norm=norm, transform=ccrs.PlateCarree())
      ax[i].coastlines(resolution='10m')
      ax[i].add_feature(cfeature.BORDERS.with_scale('10m'), linewidth=1.5)
      ax[i].add_feature(cfeature.STATES.with_scale('10m'), linewidth=1.0)  # Reduced width

      # Optional: Comment out or remove the following line if it’s adding unwanted artifacts
      # ax[i].add_feature(cfeature.NaturalEarthFeature('cultural', 'admin_1_states_provinces_lines', '10m'),
      #                   edgecolor='black', linewidth=0.5)
      ax[i].add_feature(USCOUNTIES.with_scale('500k'), edgecolor='gray', linewidth =0.7)
      ax[i].set_extent([-84.8, -74, 31, 39])
      
      #ax[i].set_title(f'Panel {i + 1}', fontsize=14, fontweight='bold')
      for i, time in enumerate([time_1, time_2, time_3]):
          ax[i].set_title(f'{time.strftime("%Y-%m-%d %H:%M Z")}' + '\n 24hr HREF LPMM [in]')
  # Shared colorbar code remains the same
  # Add overall title

  current_date = datetime.now()
  date_1 = current_date.strftime('%Y%m%d')

  previous_date = current_date + timedelta(days=1)
  date_2 = previous_date.strftime('%Y%m%d')



  # Set the suptitle with proper string concatenation
  fig.suptitle('24hr HREF LPMM [in] dprog/dt\nValid: 00Z ' + date_1 + ' to 00Z ' + date_2, 
               fontsize=18, fontweight='bold', y=0.88)

  #if you wanna mess with a smaller font subtitle:
      
  # # Set the main title (suptitle)
  # fig.suptitle('24hr HREF LPMM [in] dprog/dt', fontsize=18, fontweight='bold', y=0.9)

  # # Add the subtitle with different styling
  # fig.text(0.5, 0.87, 'Valid: 12Z ' + date_1 + ' to 12Z ' + date_2,
  #          ha='center', fontsize=12, fontweight='normal')  # Adjust fontsize and fontweight

  # Add a single color bar below the plots
  cbar = fig.colorbar(cs, ax=ax, orientation='horizontal', pad=0.05, fraction=0.02, aspect=50)
  cbar.set_label('Precipitation (inches)', fontsize=16, fontweight='bold')
  cbar.ax.tick_params(labelsize=12)

  #plt.show()

  #save = '/Users/nicholasluchetti/Documents/HREF_Plots/'
  #fig.savefig(save + 'HREF_LPMM_RUN_COMPARE', dpi=300, bbox_inches='tight')
  #fig.savefig(save + 'HREF_LPMM_RUN_COMPARE', dpi=300)
   #save = '/Users/nicholasluchetti/Documents/HREF_Plots/'
    #fig.savefig(save + 'HREF_LPMM_THRESHOLD_COMPARE', dpi=300, bbox_inches='tight')
# Create a folder for the output if it doesn't exist
    output_folder = "href_data"
    if not os.path.exists(output_folder):
    os.makedirs(output_folder)

   # Save to the local folder instead of your Documents folder
   fig.savefig(f'{output_folder}/HREF_LPMM_RUN_COMPARE.png', dpi=300)


  #%%

  '''Now we want to make a 4 panel plot where we mask 
  different thresholds (>3, >6, >9, >12 inches) and plot them 
  like a paintball plot'''

  #First start by masking out the different thresholds: 
      
  ## >3 inches:
      
  three_inch_run_1 = np.ma.masked_less(total_precip_1, 3)
  three_inch_run_2 = np.ma.masked_less(total_precip_2, 3)
  three_inch_run_3 = np.ma.masked_less(total_precip_3, 3)


  #> 6 inches
  six_inch_run_1 = np.ma.masked_less(total_precip_1, 6)
  six_inch_run_2 = np.ma.masked_less(total_precip_2, 6)
  six_inch_run_3 = np.ma.masked_less(total_precip_3, 6)

  #> 9 inches
  nine_inch_run_1 = np.ma.masked_less(total_precip_1, 9)
  nine_inch_run_2 = np.ma.masked_less(total_precip_2, 9)
  nine_inch_run_3 = np.ma.masked_less(total_precip_3, 9)

  #> 12 inches
  twelve_inch_run_1 = np.ma.masked_less(total_precip_1, 12)
  twelve_inch_run_2 = np.ma.masked_less(total_precip_2, 12)
  twelve_inch_run_3 = np.ma.masked_less(total_precip_3, 12)


  import numpy as np
  import matplotlib.pyplot as plt
  import cartopy.crs as ccrs
  import cartopy.feature as cfeature
  import matplotlib.colors as mcolors
  from matplotlib.lines import Line2D  # For custom legends

  # Create a 2x2 subplot grid
  fig, ax = plt.subplots(2, 2, figsize=(15, 12), subplot_kw={'projection': ccrs.PlateCarree()})
  #fig.subplots_adjust(top=1.5, bottom=0.5, wspace=0.04, hspace=0.9,left=0.05, right=0.95)

  # Define the shades of blue for the three simulations (darkest to lightest)
  blue_shades = ['#00008B', '#4169E1', '#87CEFA']  # Dark blue, medium blue, light blue
  alpha_values = [1.0, 0.6, 0.4]  # Transparency for each simulation: simulation 1 = opaque, others faded

  # Define a list of all thresholds and corresponding masked arrays
  thresholds = [3, 6, 9, 12]
  thresholded_data = [
      [three_inch_run_1, three_inch_run_2, three_inch_run_3],  # > 3 inches
      [six_inch_run_1, six_inch_run_2, six_inch_run_3],        # > 6 inches
      [nine_inch_run_1, nine_inch_run_2, nine_inch_run_3],     # > 9 inches
      [twelve_inch_run_1, twelve_inch_run_2, twelve_inch_run_3] # > 12 inches
  ]

  # Titles for the subplots
  titles = ['> 3 inches', '> 6 inches', '> 9 inches', '> 12 inches']

  # Loop through each subplot and plot the corresponding masked data
  for i, (thresh, data) in enumerate(zip(thresholds, thresholded_data)):
      row, col = divmod(i, 2)  # Determine the subplot position (2x2 grid)
      
      # For each simulation, apply the mask and plot with different shades of blue
      for j, (precip_data, shade, alpha) in enumerate(zip(data, blue_shades, alpha_values)):
          # Apply the mask for values less than the threshold
          masked_data = np.ma.masked_less(precip_data, thresh)

          # Plot only the masked data (values above the threshold), applying transparency
          cs = ax[row, col].contourf(lons, lats, masked_data, cmap=mcolors.ListedColormap([shade]),
                                     levels=[0, np.nanmax(masked_data)], transform=ccrs.PlateCarree(), alpha=alpha)

      # Add coastlines and features
      ax[row, col].coastlines(resolution='10m')
      ax[row, col].add_feature(cfeature.BORDERS.with_scale('10m'), linewidth=1.5)
      ax[row, col].add_feature(cfeature.STATES.with_scale('10m'), linewidth=1.0)
      
      # Set extent (optional based on your data)
      ax[row, col].set_extent([-84.8, -74, 31, 39])
     
      
      
      # Set the title for each subplot
      ax[row, col].set_title(titles[i], fontsize=14, fontweight='bold')

      # Create a custom legend with patches for each simulation
      legend_labels = [f'{time_1.strftime("%Y-%m-%d %H:%M Z")}', f'{time_2.strftime("%Y-%m-%d %H:%M Z")}', f'{time_3.strftime("%Y-%m-%d %H:%M Z")}']
      legend_colors = blue_shades  # Corresponding to each simulation
      legend_alphas = alpha_values  # Alpha values for transparency
      legend_patches = [Line2D([0], [0], marker='o', color='w', markerfacecolor=shade, markersize=10, alpha=alpha,
                               label=label) for shade, alpha, label in zip(legend_colors, legend_alphas, legend_labels)]
      
      # Add the legend to the upper right of the current subplot
      ax[row, col].legend(handles=legend_patches, loc='lower right', fontsize=8, title='HREF Run', title_fontsize=10)

  # Set the suptitle with proper string concatenation
  fig.suptitle('24hr HREF LPMM [in] \nValid: 00Z ' + date_1 + ' to 00Z ' + date_2, 
               fontsize=14, fontweight='bold', y=0.98)

  # Add a color bar (optional, but useful for understanding the plot)
  #cbar = fig.colorbar(cs, ax=ax.ravel().tolist(), orientation='horizontal', pad=0.05, fraction=0.02, aspect=50)
  #cbar.set_label('Precipitation (inches)', fontsize=12)
  plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1, wspace=0.05, hspace=0.25)

  # Show the plot
  plt.tight_layout()
  #plt.show()

 # save = '/Users/nicholasluchetti/Documents/HREF_Plots/'
  #fig.savefig(save + 'HREF_LPMM_THRESHOLD_COMPARE', dpi=300, bbox_inches='tight')
# Create a folder for the output if it doesn't exist
    output_folder = "href_data"
    if not os.path.exists(output_folder):
    os.makedirs(output_folder)

   # Save to the local folder instead of your Documents folder
   fig.savefig(f'{output_folder}/HREF_LPMM_THRESHOLD_COMPARE.png', dpi=300)


else:
    print("Time is out of expected range. Double-check your setup.")



# #%%
# import os
# import numpy as np
# import requests
# import pygrib
# import matplotlib.pyplot as plt
# import matplotlib.colors as mcolors
# import cartopy.crs as ccrs
# import cartopy.feature as cfeature
# #from cartopy.feature import USCOUNTIES
# import pandas as pd

# from datetime import datetime, timedelta
# from urllib.request import urlopen
# import cartopy.crs as ccrs
# import cartopy.feature as cfeature
# import matplotlib.colors as mcolors
# import matplotlib.pyplot as plt
# from metpy.plots import USCOUNTIES
# from metpy.units import masked_array, units
# from netCDF4 import Dataset
# import pandas as pd
# import matplotlib.patheffects as PathEffects

# import os
# import datetime

# # Define the base URL
# base_url = "https://nomads.ncep.noaa.gov/pub/data/nccf/com/href/prod/"

# current_date = datetime.datetime.now()
# date = current_date.strftime('%Y%m%d')

# # Specify the date and forecast hour
# #date = "20241108"  # Example date
# hour_str = "12"  # Example hour (can be 00Z or 12Z)
# forecast_hour = "01"  # First forecast hour

# # Define the output directory and file path
# output_dir = "href_data"
# os.makedirs(output_dir, exist_ok=True)

# # Conversion factor from mm to inches
# mm_to_inch = 0.0393701  # 1 mm = 0.0393701 inches

# # Create the plot projection
# pc_proj = ccrs.PlateCarree()

# # Initialize lists to hold hourly precipitation data
# all_hourly_precip = []

# # Function to unpack "Total precipitation" variable based on parameterName
# def unpack_total_precipitation(grib_path):
#     try:
#         with pygrib.open(grib_path) as grb_file:
#             for grb_message in grb_file:
#                 if grb_message.parameterName == "Total precipitation":
#                     data, lats, lons = grb_message.data()
#                     print(f"Extracted Total precipitation from {grib_path}")
#                     return data, lats, lons
#             print(f"Total precipitation not found in {grib_path}")
#     except Exception as e:
#         print(f"Error reading {grib_path}: {e}")
#     return None, None, None

# # Function to download files with User-Agent header
# def download_file(url, save_path):
#     headers = {
#         "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/87.0.4280.88 Safari/537.36"
#     }
#     try:
#         response = requests.get(url, headers=headers)
#         response.raise_for_status()  # Check if request was successful
#         with open(save_path, 'wb') as file:
#             file.write(response.content)
#         print(f"Downloaded {url}")
#     except requests.exceptions.RequestException as e:
#         print(f"Failed to download {url}: {e}")

# # Loop through forecast hours (from 01 to 24)
# for i in range(1, 24):
#     forecast_hour_str = f"{i:02d}"  # e.g., '01', '02', ..., '24'
#     file_url = f"{base_url}/href.{date}/ensprod/href.t{hour_str}z.conus.lpmm.f{forecast_hour_str}.grib2"
#     save_path = os.path.join(output_dir, f"href_{date}_t{hour_str}_f{forecast_hour_str}.grib2")
    
#     # Download the file for each forecast hour
#     download_file(file_url, save_path)

#     # Unpack the "Total precipitation" variable from the downloaded file
#     hourly_precip, lats, lons = unpack_total_precipitation(save_path)
    
#     # If precipitation data is successfully extracted, store it
#     if hourly_precip is not None:
#         hourly_precip_inch = hourly_precip * mm_to_inch
#         all_hourly_precip.append(hourly_precip_inch)

# # Sum the hourly precipitation data for the 24 hours
# total_precip_1 = np.sum(all_hourly_precip, axis=0)

# # Create the figure and axes for plotting
# fig = plt.figure(figsize=(20, 20))
# ax = fig.add_subplot(111, projection=pc_proj)
# ax.set_extent([-74, -84.8, 31, 39])  # Area of interest

# # Draw coastlines, state and country boundaries, and edges of the map
# ax.coastlines(resolution='10m')
# ax.add_feature(cfeature.BORDERS.with_scale('10m'), linewidth=1.5)
# ax.add_feature(cfeature.STATES.with_scale('10m'), linewidth=2.0)
# ax.add_feature(USCOUNTIES.with_scale('500k'), edgecolor='black')

# # Define the contour levels for precipitation
# clevs = [0.00, 0.01, 0.1, 0.25, 0.50, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0,
#          6.0, 8.0, 10., 12.0, 15.0, 18.0]

# # Define the color map for the precipitation levels
# cmap_data = [(1.0, 1.0, 1.0),
#              (0.3137255012989044, 0.8156862854957581, 0.8156862854957581),
#              (0.0, 1.0, 1.0),
#              (0.0, 0.8784313797950745, 0.501960813999176),
#              (0.0, 0.7529411911964417, 0.0),
#              (0.501960813999176, 0.8784313797950745, 0.0),
#              (1.0, 1.0, 0.0),
#              (1.0, 0.6274510025978088, 0.0),
#              (1.0, 0.0, 0.0),
#              (1.0, 0.125490203499794, 0.501960813999176),
#              (0.9411764740943909, 0.250980406999588, 1.0),
#              (0.501960813999176, 0.125490203499794, 1.0),
#              (0.250980406999588, 0.250980406999588, 1.0),
#              (0.125490203499794, 0.125490203499794, 0.501960813999176),
#              (0.125490203499794, 0.125490203499794, 0.125490203499794),
#              (0.501960813999176, 0.501960813999176, 0.501960813999176),
#              (0.8784313797950745, 0.8784313797950745, 0.8784313797950745),
#              (0.9333333373069763, 0.8313725590705872, 0.7372549176216125),
#              (0.8549019694328308, 0.6509804129600525, 0.47058823704719543),
#              (0.6274510025978088, 0.42352941632270813, 0.23529411852359772),
#              (0.4000000059604645, 0.20000000298023224, 0.0)]
# cmap = mcolors.ListedColormap(cmap_data, 'precipitation')
# norm = mcolors.BoundaryNorm(clevs, cmap.N)

# # Plot the contour fill
# cs = ax.contourf(lons, lats, total_precip_1, clevs, alpha=0.5, cmap=cmap, norm=norm, transform=pc_proj)

# # Add the color bar
# cbar = plt.colorbar(cs, orientation='vertical', label='Precipitation (inches)', ticks=clevs)
# cbar.set_label('Precipitation (inches)', fontsize=18, fontweight='bold')
# cbar.ax.tick_params(labelsize=14, which='both')
# cbar.ax.set_yticklabels(["{:4.2f}".format(i) for i in clevs])

# # Convert time to datetime and format the title
# time_1 = pd.to_datetime(date + ' ' + hour_str + 'Z')
# ax.set_title(f'{time_1.strftime("%Y-%m-%d %H:%M Z")}' + '\n 24hr HREF LPMM [in]', fontsize=16, fontweight='bold')

# # # Add city labels
# # locs = {"Raleigh": [-78.6382, 35.7796], "Fayetteville":[-78.8784, 35.0527], "Goldsboro":[-77.9928, 35.3849], 
# #         "Greensboro": [-79.7920, 36.0726], "Rocky Mount": [-77.7905, 35.9382], "Asheboro": [-79.8136, 35.7079], 
# #         "Rockingham": [-79.7739, 34.9393], "Roxboro": [-78.9828, 36.3938], "Henderson": [-78.3992, 36.3294], 
# #         "Roanoke Rapids": [-77.6541, 36.4615]}

# # for key in locs.keys():
# #     x, y = locs[key]
# #     ax.scatter(x, y, c='k', zorder=1, s=50)
# #     ax.annotate(key, xy=(x, y), color='k', fontsize=8, xytext=(x+0.045, y-0.05), fontweight='bold', 
# #                 bbox=dict(boxstyle="round", fc="black", ec="b", alpha=0.2))

# # Display the plot
# #plt.show()


# #%%


# import os
# import numpy as np
# import requests
# import pygrib
# import matplotlib.pyplot as plt
# import matplotlib.colors as mcolors
# import cartopy.crs as ccrs
# import cartopy.feature as cfeature
# #from cartopy.feature import USCOUNTIES
# import pandas as pd

# from datetime import datetime, timedelta
# from urllib.request import urlopen
# import cartopy.crs as ccrs
# import cartopy.feature as cfeature
# import matplotlib.colors as mcolors
# import matplotlib.pyplot as plt
# from metpy.plots import USCOUNTIES
# from metpy.units import masked_array, units
# from netCDF4 import Dataset
# import pandas as pd
# import matplotlib.patheffects as PathEffects

# import os
# import datetime

# # Define the base URL
# base_url = "https://nomads.ncep.noaa.gov/pub/data/nccf/com/href/prod/"

# current_date = datetime.datetime.now()
# date = current_date.strftime('%Y%m%d')

# # Specify the date and forecast hour
# #date = "20241108"  # Example date
# hour_str = "00"  # Example hour (can be 00Z or 12Z)
# forecast_hour = "01"  # First forecast hour

# # Define the output directory and file path
# output_dir = "href_data"
# os.makedirs(output_dir, exist_ok=True)

# # Conversion factor from mm to inches
# mm_to_inch = 0.0393701  # 1 mm = 0.0393701 inches

# # Create the plot projection
# pc_proj = ccrs.PlateCarree()

# # Initialize lists to hold hourly precipitation data
# all_hourly_precip = []

# # Function to unpack "Total precipitation" variable based on parameterName
# def unpack_total_precipitation(grib_path):
#     try:
#         with pygrib.open(grib_path) as grb_file:
#             for grb_message in grb_file:
#                 if grb_message.parameterName == "Total precipitation":
#                     data, lats, lons = grb_message.data()
#                     print(f"Extracted Total precipitation from {grib_path}")
#                     return data, lats, lons
#             print(f"Total precipitation not found in {grib_path}")
#     except Exception as e:
#         print(f"Error reading {grib_path}: {e}")
#     return None, None, None

# # Function to download files with User-Agent header
# def download_file(url, save_path):
#     headers = {
#         "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/87.0.4280.88 Safari/537.36"
#     }
#     try:
#         response = requests.get(url, headers=headers)
#         response.raise_for_status()  # Check if request was successful
#         with open(save_path, 'wb') as file:
#             file.write(response.content)
#         print(f"Downloaded {url}")
#     except requests.exceptions.RequestException as e:
#         print(f"Failed to download {url}: {e}")

# # Loop through forecast hours (from 01 to 24)
# for i in range(13, 36):
#     forecast_hour_str = f"{i:02d}"  # e.g., '01', '02', ..., '24'
#     file_url = f"{base_url}/href.{date}/ensprod/href.t{hour_str}z.conus.lpmm.f{forecast_hour_str}.grib2"
#     save_path = os.path.join(output_dir, f"href_{date}_t{hour_str}_f{forecast_hour_str}.grib2")
    
#     # Download the file for each forecast hour
#     download_file(file_url, save_path)

#     # Unpack the "Total precipitation" variable from the downloaded file
#     hourly_precip, lats, lons = unpack_total_precipitation(save_path)
    
#     # If precipitation data is successfully extracted, store it
#     if hourly_precip is not None:
#         hourly_precip_inch = hourly_precip * mm_to_inch
#         all_hourly_precip.append(hourly_precip_inch)

# # Sum the hourly precipitation data for the 24 hours
# total_precip_2 = np.sum(all_hourly_precip, axis=0)

# # Create the figure and axes for plotting
# fig = plt.figure(figsize=(20, 20))
# ax = fig.add_subplot(111, projection=pc_proj)
# ax.set_extent([-74, -84.8, 31, 39])  # Area of interest

# # Draw coastlines, state and country boundaries, and edges of the map
# ax.coastlines(resolution='10m')
# ax.add_feature(cfeature.BORDERS.with_scale('10m'), linewidth=1.5)
# ax.add_feature(cfeature.STATES.with_scale('10m'), linewidth=2.0)
# ax.add_feature(USCOUNTIES.with_scale('500k'), edgecolor='black')

# # Define the contour levels for precipitation
# clevs = [0.00, 0.01, 0.1, 0.25, 0.50, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0,
#          6.0, 8.0, 10., 12.0, 15.0, 18.0]

# # Define the color map for the precipitation levels
# cmap_data = [(1.0, 1.0, 1.0),
#              (0.3137255012989044, 0.8156862854957581, 0.8156862854957581),
#              (0.0, 1.0, 1.0),
#              (0.0, 0.8784313797950745, 0.501960813999176),
#              (0.0, 0.7529411911964417, 0.0),
#              (0.501960813999176, 0.8784313797950745, 0.0),
#              (1.0, 1.0, 0.0),
#              (1.0, 0.6274510025978088, 0.0),
#              (1.0, 0.0, 0.0),
#              (1.0, 0.125490203499794, 0.501960813999176),
#              (0.9411764740943909, 0.250980406999588, 1.0),
#              (0.501960813999176, 0.125490203499794, 1.0),
#              (0.250980406999588, 0.250980406999588, 1.0),
#              (0.125490203499794, 0.125490203499794, 0.501960813999176),
#              (0.125490203499794, 0.125490203499794, 0.125490203499794),
#              (0.501960813999176, 0.501960813999176, 0.501960813999176),
#              (0.8784313797950745, 0.8784313797950745, 0.8784313797950745),
#              (0.9333333373069763, 0.8313725590705872, 0.7372549176216125),
#              (0.8549019694328308, 0.6509804129600525, 0.47058823704719543),
#              (0.6274510025978088, 0.42352941632270813, 0.23529411852359772),
#              (0.4000000059604645, 0.20000000298023224, 0.0)]
# cmap = mcolors.ListedColormap(cmap_data, 'precipitation')
# norm = mcolors.BoundaryNorm(clevs, cmap.N)

# # Plot the contour fill
# cs = ax.contourf(lons, lats, total_precip_2, clevs, alpha=0.5, cmap=cmap, norm=norm, transform=pc_proj)

# # Add the color bar
# cbar = plt.colorbar(cs, orientation='vertical', label='Precipitation (inches)', ticks=clevs)
# cbar.set_label('Precipitation (inches)', fontsize=18, fontweight='bold')
# cbar.ax.tick_params(labelsize=14, which='both')
# cbar.ax.set_yticklabels(["{:4.2f}".format(i) for i in clevs])

# # Convert time to datetime and format the title
# time_2 = pd.to_datetime(date + ' ' + hour_str + 'Z')
# ax.set_title(f'{time_2.strftime("%Y-%m-%d %H:%M Z")}' + '\n 24hr HREF LPMM [in]', fontsize=16, fontweight='bold')

# # # Add city labels
# # locs = {"Raleigh": [-78.6382, 35.7796], "Fayetteville":[-78.8784, 35.0527], "Goldsboro":[-77.9928, 35.3849], 
# #         "Greensboro": [-79.7920, 36.0726], "Rocky Mount": [-77.7905, 35.9382], "Asheboro": [-79.8136, 35.7079], 
# #         "Rockingham": [-79.7739, 34.9393], "Roxboro": [-78.9828, 36.3938], "Henderson": [-78.3992, 36.3294], 
# #         "Roanoke Rapids": [-77.6541, 36.4615]}

# # for key in locs.keys():
# #     x, y = locs[key]
# #     ax.scatter(x, y, c='k', zorder=1, s=50)
# #     ax.annotate(key, xy=(x, y), color='k', fontsize=8, xytext=(x+0.045, y-0.05), fontweight='bold', 
# #                 bbox=dict(boxstyle="round", fc="black", ec="b", alpha=0.2))

# # Display the plot
# #plt.show()

# #%%


# import os
# import numpy as np
# import requests
# import pygrib
# import matplotlib.pyplot as plt
# import matplotlib.colors as mcolors
# import cartopy.crs as ccrs
# import cartopy.feature as cfeature
# #from cartopy.feature import USCOUNTIES
# import pandas as pd

# from datetime import datetime, timedelta
# from urllib.request import urlopen
# import cartopy.crs as ccrs
# import cartopy.feature as cfeature
# import matplotlib.colors as mcolors
# import matplotlib.pyplot as plt
# from metpy.plots import USCOUNTIES
# from metpy.units import masked_array, units
# from netCDF4 import Dataset
# import pandas as pd
# import matplotlib.patheffects as PathEffects

# import os
# import datetime

# # Define the base URL
# base_url = "https://nomads.ncep.noaa.gov/pub/data/nccf/com/href/prod/"

# current_date = datetime.datetime.now()
# #date = current_date.strftime('%Y%m%d')

# previous_date = current_date - datetime.timedelta(days=1)
# date = previous_date.strftime('%Y%m%d')

# # Specify the date and forecast hour
# #date = "20241108"  # Example date
# hour_str = "12"  # Example hour (can be 00Z or 12Z)
# forecast_hour = "01"  # First forecast hour

# # Define the output directory and file path
# output_dir = "href_data"
# os.makedirs(output_dir, exist_ok=True)

# # Conversion factor from mm to inches
# mm_to_inch = 0.0393701  # 1 mm = 0.0393701 inches

# # Create the plot projection
# pc_proj = ccrs.PlateCarree()

# # Initialize lists to hold hourly precipitation data
# all_hourly_precip = []

# # Function to unpack "Total precipitation" variable based on parameterName
# def unpack_total_precipitation(grib_path):
#     try:
#         with pygrib.open(grib_path) as grb_file:
#             for grb_message in grb_file:
#                 if grb_message.parameterName == "Total precipitation":
#                     data, lats, lons = grb_message.data()
#                     print(f"Extracted Total precipitation from {grib_path}")
#                     return data, lats, lons
#             print(f"Total precipitation not found in {grib_path}")
#     except Exception as e:
#         print(f"Error reading {grib_path}: {e}")
#     return None, None, None

# # Function to download files with User-Agent header
# def download_file(url, save_path):
#     headers = {
#         "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/87.0.4280.88 Safari/537.36"
#     }
#     try:
#         response = requests.get(url, headers=headers)
#         response.raise_for_status()  # Check if request was successful
#         with open(save_path, 'wb') as file:
#             file.write(response.content)
#         print(f"Downloaded {url}")
#     except requests.exceptions.RequestException as e:
#         print(f"Failed to download {url}: {e}")

# # Loop through forecast hours (from 01 to 24)
# for i in range(25, 48):
#     forecast_hour_str = f"{i:02d}"  # e.g., '01', '02', ..., '24'
#     file_url = f"{base_url}/href.{date}/ensprod/href.t{hour_str}z.conus.lpmm.f{forecast_hour_str}.grib2"
#     save_path = os.path.join(output_dir, f"href_{date}_t{hour_str}_f{forecast_hour_str}.grib2")
    
#     # Download the file for each forecast hour
#     download_file(file_url, save_path)

#     # Unpack the "Total precipitation" variable from the downloaded file
#     hourly_precip, lats, lons = unpack_total_precipitation(save_path)
    
#     # If precipitation data is successfully extracted, store it
#     if hourly_precip is not None:
#         hourly_precip_inch = hourly_precip * mm_to_inch
#         all_hourly_precip.append(hourly_precip_inch)

# # Sum the hourly precipitation data for the 24 hours
# total_precip_3 = np.sum(all_hourly_precip, axis=0)

# # Create the figure and axes for plotting
# fig = plt.figure(figsize=(20, 20))
# ax = fig.add_subplot(111, projection=pc_proj)
# ax.set_extent([-74, -84.8, 31, 39])  # Area of interest

# # Draw coastlines, state and country boundaries, and edges of the map
# ax.coastlines(resolution='10m')
# ax.add_feature(cfeature.BORDERS.with_scale('10m'), linewidth=1.5)
# ax.add_feature(cfeature.STATES.with_scale('10m'), linewidth=2.0)
# ax.add_feature(USCOUNTIES.with_scale('500k'), edgecolor='black')

# # Define the contour levels for precipitation
# clevs = [0.00, 0.01, 0.1, 0.25, 0.50, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0,
#          6.0, 8.0, 10., 12.0, 15.0, 18.0]

# # Define the color map for the precipitation levels
# cmap_data = [(1.0, 1.0, 1.0),
#              (0.3137255012989044, 0.8156862854957581, 0.8156862854957581),
#              (0.0, 1.0, 1.0),
#              (0.0, 0.8784313797950745, 0.501960813999176),
#              (0.0, 0.7529411911964417, 0.0),
#              (0.501960813999176, 0.8784313797950745, 0.0),
#              (1.0, 1.0, 0.0),
#              (1.0, 0.6274510025978088, 0.0),
#              (1.0, 0.0, 0.0),
#              (1.0, 0.125490203499794, 0.501960813999176),
#              (0.9411764740943909, 0.250980406999588, 1.0),
#              (0.501960813999176, 0.125490203499794, 1.0),
#              (0.250980406999588, 0.250980406999588, 1.0),
#              (0.125490203499794, 0.125490203499794, 0.501960813999176),
#              (0.125490203499794, 0.125490203499794, 0.125490203499794),
#              (0.501960813999176, 0.501960813999176, 0.501960813999176),
#              (0.8784313797950745, 0.8784313797950745, 0.8784313797950745),
#              (0.9333333373069763, 0.8313725590705872, 0.7372549176216125),
#              (0.8549019694328308, 0.6509804129600525, 0.47058823704719543),
#              (0.6274510025978088, 0.42352941632270813, 0.23529411852359772),
#              (0.4000000059604645, 0.20000000298023224, 0.0)]
# cmap = mcolors.ListedColormap(cmap_data, 'precipitation')
# norm = mcolors.BoundaryNorm(clevs, cmap.N)

# # Plot the contour fill
# cs = ax.contourf(lons, lats, total_precip_3, clevs, alpha=0.5, cmap=cmap, norm=norm, transform=pc_proj)

# # Add the color bar
# cbar = plt.colorbar(cs, orientation='vertical', label='Precipitation (inches)', ticks=clevs)
# cbar.set_label('Precipitation (inches)', fontsize=18, fontweight='bold')
# cbar.ax.tick_params(labelsize=14, which='both')
# cbar.ax.set_yticklabels(["{:4.2f}".format(i) for i in clevs])

# # Convert time to datetime and format the title
# time_3 = pd.to_datetime(date + ' ' + hour_str + 'Z')
# ax.set_title(f'{time_3.strftime("%Y-%m-%d %H:%M Z")}' + '\n 24hr HREF LPMM [in]', fontsize=16, fontweight='bold')

# # # Add city labels
# # locs = {"Raleigh": [-78.6382, 35.7796], "Fayetteville":[-78.8784, 35.0527], "Goldsboro":[-77.9928, 35.3849], 
# #         "Greensboro": [-79.7920, 36.0726], "Rocky Mount": [-77.7905, 35.9382], "Asheboro": [-79.8136, 35.7079], 
# #         "Rockingham": [-79.7739, 34.9393], "Roxboro": [-78.9828, 36.3938], "Henderson": [-78.3992, 36.3294], 
# #         "Roanoke Rapids": [-77.6541, 36.4615]}

# # for key in locs.keys():
# #     x, y = locs[key]
# #     ax.scatter(x, y, c='k', zorder=1, s=50)
# #     ax.annotate(key, xy=(x, y), color='k', fontsize=8, xytext=(x+0.045, y-0.05), fontweight='bold', 
# #                 bbox=dict(boxstyle="round", fc="black", ec="b", alpha=0.2))

# # Display the plot
# #plt.show()


# #If you want to remove the grib files use these lines:

# # import os
# # import glob

# # # Specify the directory where the .grib2 files are stored
# # grib_directory = "href_data"

# # # Function to delete all .grib2 files in the specified directory
# # def clean_up_grib_files(directory):
# #     grib_files = glob.glob(os.path.join(directory, "*.grib2"))
# #     for file_path in grib_files:
# #         try:
# #             os.remove(file_path)
# #             print(f"Deleted {file_path}")
# #         except OSError as e:
# #             print(f"Error deleting {file_path}: {e}")

# # # Run the cleanup function at the end of the script
# # clean_up_grib_files(grib_directory)


# #%%



# '''now we want to combine the three plots side by side'''


# import matplotlib.pyplot as plt
# import cartopy.crs as ccrs
# import cartopy.feature as cfeature
# import matplotlib.colors as mcolors
# import numpy as np

# # Create the figure and axes for plotting
# fig, ax = plt.subplots(1, 3, figsize=(18, 8), subplot_kw={'projection': ccrs.PlateCarree()})
# #fig.subplots_adjust(wspace=0.05)
# fig.subplots_adjust(top=1.5, bottom=0.2, wspace=0.04, hspace=0.9,left=0.05, right=0.95)
# #fig.tight_layout(pad=8.5, w_pad=0.5, h_pad=2.0)

# # Set up the contour levels and color map
# clevs = [0.00, 0.01, 0.1, 0.25, 0.50, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0,
#          6.0, 8.0, 10., 12.0, 15.0, 18.0]
# cmap_data = [
#     (1.0, 1.0, 1.0), (0.3137255012989044, 0.8156862854957581, 0.8156862854957581),
#     (0.0, 1.0, 1.0), (0.0, 0.8784313797950745, 0.501960813999176),
#     (0.0, 0.7529411911964417, 0.0), (0.501960813999176, 0.8784313797950745, 0.0),
#     (1.0, 1.0, 0.0), (1.0, 0.6274510025978088, 0.0), (1.0, 0.0, 0.0),
#     (1.0, 0.125490203499794, 0.501960813999176), (0.9411764740943909, 0.250980406999588, 1.0),
#     (0.501960813999176, 0.125490203499794, 1.0), (0.250980406999588, 0.250980406999588, 1.0),
#     (0.125490203499794, 0.125490203499794, 0.501960813999176), (0.125490203499794, 0.125490203499794, 0.125490203499794),
#     (0.501960813999176, 0.501960813999176, 0.501960813999176), (0.8784313797950745, 0.8784313797950745, 0.8784313797950745),
#     (0.9333333373069763, 0.8313725590705872, 0.7372549176216125), (0.8549019694328308, 0.6509804129600525, 0.47058823704719543),
#     (0.6274510025978088, 0.42352941632270813, 0.23529411852359772), (0.4000000059604645, 0.20000000298023224, 0.0)
# ]
# cmap = mcolors.ListedColormap(cmap_data, 'precipitation')
# norm = mcolors.BoundaryNorm(clevs, cmap.N)

# # Sample data arrays (replace with your own)
# #lons, lats = np.meshgrid(np.linspace(-84.8, -74, 1473), np.linspace(31, 39, 1473))
# total_precip_1 = total_precip_1
# total_precip_2 = total_precip_2
# total_precip_3 = total_precip_3

# # Plot each panel
# # Inside the loop for each subplot, comment out or modify the state feature layer
# for i, data in enumerate([total_precip_1, total_precip_2, total_precip_3]):
#     cs = ax[i].contourf(lons, lats, data, clevs, alpha=0.5, cmap=cmap, norm=norm, transform=ccrs.PlateCarree())
#     ax[i].coastlines(resolution='10m')
#     ax[i].add_feature(cfeature.BORDERS.with_scale('10m'), linewidth=1.5)
#     ax[i].add_feature(cfeature.STATES.with_scale('10m'), linewidth=1.0)  # Reduced width

#     # Optional: Comment out or remove the following line if it’s adding unwanted artifacts
#     # ax[i].add_feature(cfeature.NaturalEarthFeature('cultural', 'admin_1_states_provinces_lines', '10m'),
#     #                   edgecolor='black', linewidth=0.5)
#     ax[i].add_feature(USCOUNTIES.with_scale('500k'), edgecolor='gray', linewidth =0.7)
#     ax[i].set_extent([-84.8, -74, 31, 39])
    
#     #ax[i].set_title(f'Panel {i + 1}', fontsize=14, fontweight='bold')
#     for i, time in enumerate([time_1, time_2, time_3]):
#         ax[i].set_title(f'{time.strftime("%Y-%m-%d %H:%M Z")}' + '\n 24hr HREF LPMM [in]')
# # Shared colorbar code remains the same
# # Add overall title

# current_date = datetime.datetime.now()
# date_1 = current_date.strftime('%Y%m%d')

# previous_date = current_date + datetime.timedelta(days=1)
# date_2 = previous_date.strftime('%Y%m%d')



# # Set the suptitle with proper string concatenation
# fig.suptitle('24hr HREF LPMM [in] dprog/dt\nValid: 12Z ' + date_1 + ' to 12Z ' + date_2, 
#              fontsize=18, fontweight='bold', y=0.88)

# #if you wanna mess with a smaller font subtitle:
    
# # # Set the main title (suptitle)
# # fig.suptitle('24hr HREF LPMM [in] dprog/dt', fontsize=18, fontweight='bold', y=0.9)

# # # Add the subtitle with different styling
# # fig.text(0.5, 0.87, 'Valid: 12Z ' + date_1 + ' to 12Z ' + date_2,
# #          ha='center', fontsize=12, fontweight='normal')  # Adjust fontsize and fontweight

# # Add a single color bar below the plots
# cbar = fig.colorbar(cs, ax=ax, orientation='horizontal', pad=0.05, fraction=0.02, aspect=50)
# cbar.set_label('Precipitation (inches)', fontsize=16, fontweight='bold')
# cbar.ax.tick_params(labelsize=12)

# #plt.show()

# save = '/Users/nicholasluchetti/Documents/HREF_Plots/'
# #fig.savefig(save + 'HREF_LPMM_RUN_COMPARE', dpi=300, bbox_inches='tight')
# fig.savefig(save + 'HREF_LPMM_RUN_COMPARE', dpi=300)



# #%%

# '''Now we want to make a 4 panel plot where we mask 
# different thresholds (>3, >6, >9, >12 inches) and plot them 
# like a paintball plot'''

# #First start by masking out the different thresholds: 
    
# ## >3 inches:
    
# three_inch_run_1 = np.ma.masked_less(total_precip_1, 3)
# three_inch_run_2 = np.ma.masked_less(total_precip_2, 3)
# three_inch_run_3 = np.ma.masked_less(total_precip_3, 3)


# #> 6 inches
# six_inch_run_1 = np.ma.masked_less(total_precip_1, 6)
# six_inch_run_2 = np.ma.masked_less(total_precip_2, 6)
# six_inch_run_3 = np.ma.masked_less(total_precip_3, 6)

# #> 9 inches
# nine_inch_run_1 = np.ma.masked_less(total_precip_1, 9)
# nine_inch_run_2 = np.ma.masked_less(total_precip_2, 9)
# nine_inch_run_3 = np.ma.masked_less(total_precip_3, 9)

# #> 12 inches
# twelve_inch_run_1 = np.ma.masked_less(total_precip_1, 12)
# twelve_inch_run_2 = np.ma.masked_less(total_precip_2, 12)
# twelve_inch_run_3 = np.ma.masked_less(total_precip_3, 12)


# #%%

# import numpy as np
# import matplotlib.pyplot as plt
# import cartopy.crs as ccrs
# import cartopy.feature as cfeature
# import matplotlib.colors as mcolors
# from matplotlib.lines import Line2D  # For custom legends

# # Create a 2x2 subplot grid
# fig, ax = plt.subplots(2, 2, figsize=(15, 12), subplot_kw={'projection': ccrs.PlateCarree()})
# #fig.subplots_adjust(top=1.5, bottom=0.5, wspace=0.04, hspace=0.9,left=0.05, right=0.95)

# # Define the shades of blue for the three simulations (darkest to lightest)
# blue_shades = ['#00008B', '#4169E1', '#87CEFA']  # Dark blue, medium blue, light blue
# alpha_values = [1.0, 0.6, 0.4]  # Transparency for each simulation: simulation 1 = opaque, others faded

# # Define a list of all thresholds and corresponding masked arrays
# thresholds = [3, 6, 9, 12]
# thresholded_data = [
#     [three_inch_run_1, three_inch_run_2, three_inch_run_3],  # > 3 inches
#     [six_inch_run_1, six_inch_run_2, six_inch_run_3],        # > 6 inches
#     [nine_inch_run_1, nine_inch_run_2, nine_inch_run_3],     # > 9 inches
#     [twelve_inch_run_1, twelve_inch_run_2, twelve_inch_run_3] # > 12 inches
# ]

# # Titles for the subplots
# titles = ['> 3 inches', '> 6 inches', '> 9 inches', '> 12 inches']

# # Loop through each subplot and plot the corresponding masked data
# for i, (thresh, data) in enumerate(zip(thresholds, thresholded_data)):
#     row, col = divmod(i, 2)  # Determine the subplot position (2x2 grid)
    
#     # For each simulation, apply the mask and plot with different shades of blue
#     for j, (precip_data, shade, alpha) in enumerate(zip(data, blue_shades, alpha_values)):
#         # Apply the mask for values less than the threshold
#         masked_data = np.ma.masked_less(precip_data, thresh)

#         # Plot only the masked data (values above the threshold), applying transparency
#         cs = ax[row, col].contourf(lons, lats, masked_data, cmap=mcolors.ListedColormap([shade]),
#                                    levels=[0, np.nanmax(masked_data)], transform=ccrs.PlateCarree(), alpha=alpha)

#     # Add coastlines and features
#     ax[row, col].coastlines(resolution='10m')
#     ax[row, col].add_feature(cfeature.BORDERS.with_scale('10m'), linewidth=1.5)
#     ax[row, col].add_feature(cfeature.STATES.with_scale('10m'), linewidth=1.0)
    
#     # Set extent (optional based on your data)
#     ax[row, col].set_extent([-84.8, -74, 31, 39])
   
    
    
#     # Set the title for each subplot
#     ax[row, col].set_title(titles[i], fontsize=14, fontweight='bold')

#     # Create a custom legend with patches for each simulation
#     legend_labels = [f'{time_1.strftime("%Y-%m-%d %H:%M Z")}', f'{time_2.strftime("%Y-%m-%d %H:%M Z")}', f'{time_3.strftime("%Y-%m-%d %H:%M Z")}']
#     legend_colors = blue_shades  # Corresponding to each simulation
#     legend_alphas = alpha_values  # Alpha values for transparency
#     legend_patches = [Line2D([0], [0], marker='o', color='w', markerfacecolor=shade, markersize=10, alpha=alpha,
#                              label=label) for shade, alpha, label in zip(legend_colors, legend_alphas, legend_labels)]
    
#     # Add the legend to the upper right of the current subplot
#     ax[row, col].legend(handles=legend_patches, loc='lower right', fontsize=8, title='HREF Run', title_fontsize=10)

# # Set the suptitle with proper string concatenation
# fig.suptitle('24hr HREF LPMM [in] \nValid: 12Z ' + date_1 + ' to 12Z ' + date_2, 
#              fontsize=14, fontweight='bold', y=0.98)

# # Add a color bar (optional, but useful for understanding the plot)
# #cbar = fig.colorbar(cs, ax=ax.ravel().tolist(), orientation='horizontal', pad=0.05, fraction=0.02, aspect=50)
# #cbar.set_label('Precipitation (inches)', fontsize=12)
# plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1, wspace=0.05, hspace=0.25)

# # Show the plot
# plt.tight_layout()
# #plt.show()

# save = '/Users/nicholasluchetti/Documents/HREF_Plots/'
# fig.savefig(save + 'HREF_LPMM_THRESHOLD_COMPARE', dpi=300, bbox_inches='tight')




#%%










