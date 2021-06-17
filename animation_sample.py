# -*- coding: utf-8 -*-
"""
Created on Mon N  9 15:35:26 2020

@author: arevill
"""
import geopandas as gpd
import matplotlib as plt
import rasterio
from rasterio import features
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib_scalebar.scalebar import ScaleBar

import glob
import salem
from salem import DataLevels, GoogleVisibleMap, Map
from PIL import Image
from pyproj import Proj, transform
import shapefile

import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import datetime as dt
import glob
from datetime import date, timedelta
import matplotlib.dates as mdates
from numpy import savetxt

# Extract point name from file name.
def get_pt(year_files):
       pt = year_files.split("_")
       #breakpoint()
       pt2= pt[7].split('.')
       return pt2

dalec_shape = 'C:/Users/arevi/OneDrive/BBSRC_IAA_project/SCRIPTS/WF_DALEC.shp'

mm_background = r'C:\Users\arevi\OneDrive\BBSRC_IAA_project\farmer_visit\WestFortune_MM.shp'

data = gpd.read_file(dalec_shape)

min_x = np.min(data.field_2)
min_y = np.min(data.field_3)
max_x = np.max(data.field_2)
max_y = np.max(data.field_3)

# 55.997983, -2.764435
# 56.010688, -2.734454

background = gpd.read_file(mm_background)

g = GoogleVisibleMap(x=[min_x, max_x], y=[min_y, max_y],
                      crs=background.crs,
                      scale=3,  # scale is for more details
                      maptype='satellite')  # try out also: 'terrain'

# g = GoogleVisibleMap(x=[-2.764435, -2.734454], y=[55.997983, 56.010688],
#                       scale=3,  # scale is for more details
#                       maptype='satellite')  # try out also: 'terrain'



data['date'] = pd.to_datetime(data['date'], format='%Y-%m-%d')

date_array = pd.date_range(start='1-01-2021', end='6-07-2021', freq='D')



for tt in date_array[range(len(date_array))]:
       print (tt)
       data_3 = data[data.date==tt]
       
       la_array = np.array(data_3['lai'])
       
       mean_la = np.mean(la_array[np.isfinite(la_array)])       
       
       mean_la_text = str(round(mean_la,2))     

       # data_2.to_file('doy1.shp')
       minx, miny, maxx, maxy = data_3.to_crs(32630).total_bounds
            
       f, ax = plt.subplots(figsize=(8,4))    
      
       background.to_crs(32630).plot(ax=ax,column='OBJECTID',legend=False,color='white', edgecolor='k', alpha=0.5)
       
       data_3.to_crs(32630).plot(ax=ax, column='lai', cmap='Greens', legend=True,norm=plt.Normalize(vmin=0, vmax=5, clip=False))
       
       ax.set_title('Field average Leaf area: ' + mean_la_text + ' m$^2$/m$^2$')   
       ax.get_xaxis().set_visible(False)
       ax.get_yaxis().set_visible(False)

       minx, miny, maxx, maxy = data_3.to_crs(32630).total_bounds
       offset= 170
       ax.set_xlim(minx-offset,maxx+offset)
       ax.set_ylim(miny-offset,maxy+offset)
       
       ax.add_artist(ScaleBar(2, units='m',))

       fig = ax.get_figure()
       
       plt.show(block=False)
       
       # add text box
       out_date = tt.strftime('%d %B, %Y')
       
       
       # NB
       #ax.text(517148,6210330, out_date, fontsize=13, bbox=dict(facecolor='White', alpha=0.8)) 
       
       # WF
       ax.text(515010,6206161, out_date, fontsize=13, bbox=dict(facecolor='White', alpha=0.8))
       
       f.tight_layout()
       
       
       out_dir = 'C:/Users/arevi/OneDrive/BBSRC_IAA_project/farmer_visit/IMAGES/WF/'
       plt.savefig(out_dir+'lai_'+str(tt)[0:10]+'.png')   
 
 
       plt.close(fig)


fp_in = out_dir+'lai_*.png'
fp_out = out_dir+"LAI_animation.gif"
 
img, *imgs = [Image.open(f) for f in sorted(glob.glob(fp_in))]

img.save(fp=fp_out, format='GIF', append_images=imgs,
         save_all=True, duration=180, loop=0)

# ##  Leaf N
for tt in date_array[range(len(date_array))]:
       print (tt)
       data_3 = data[data.date==tt]

       ln_array = np.array(data_3['leaf_n'])
       
       mean_ln = np.mean(ln_array[np.isfinite(ln_array)])       
       
       mean_ln_text = str(round(mean_ln,2))     
       
       # data_2.to_file('doy1.shp')
       minx, miny, maxx, maxy = data_3.to_crs(32630).total_bounds
            
       f, ax = plt.subplots(figsize=(8,4))
       
       background.to_crs(32630).plot(ax=ax,column='OBJECTID',legend=False,color='white', edgecolor='k', alpha=0.5)
       
       data_3.to_crs(32630).plot(ax=ax, column='leaf_n', cmap='Blues', legend=True,norm=plt.Normalize(vmin=1, vmax=10, clip=False))
       
       ax.set_title('Field averaged canopy nitrogen: ' + mean_ln_text + ' g N/m$^2$')   
       ax.get_xaxis().set_visible(False)
       ax.get_yaxis().set_visible(False)

       minx, miny, maxx, maxy = data_3.to_crs(32630).total_bounds
       offset= 170
       ax.set_xlim(minx-offset,maxx+offset)
       ax.set_ylim(miny-offset,maxy+offset)
       
       ax.add_artist(ScaleBar(2, units='m',))

       fig = ax.get_figure()
       
       plt.show(block=False)
       
       # add text box
       out_date = tt.strftime('%d %B, %Y')
       # NB
       #ax.text(517148,6210330, out_date, fontsize=13, bbox=dict(facecolor='White', alpha=0.8)) 
       
       # WF
       ax.text(515010,6206161, out_date, fontsize=13, bbox=dict(facecolor='White', alpha=0.8))
       
       f.tight_layout()
             
       out_dir = 'C:/Users/arevi/OneDrive/BBSRC_IAA_project/farmer_visit/IMAGES/WF/'
       plt.savefig(out_dir+'leaf_n_'+str(tt)[0:10]+'.png')   
 
       plt.close(fig)


fp_in = out_dir+'leaf_n_*.png'
fp_out = out_dir+"leaf_n_animation_v1.gif"
 
img, *imgs = [Image.open(f) for f in sorted(glob.glob(fp_in))]

img.save(fp=fp_out, format='GIF', append_images=imgs,
         save_all=True, duration=180, loop=0)


# ##  GS
for tt in date_array[range(len(date_array))]:
       print (tt)
       data_3 = data[data.date==tt]

       mean_gs = str(round(np.nanmean(np.array(data_3['gs'])),2))       

       # data_2.to_file('doy1.shp')
       minx, miny, maxx, maxy = data_3.to_crs(32630).total_bounds
            
       f, ax = plt.subplots(figsize=(8,4))
       
       background.to_crs(32630).plot(ax=ax,column='OBJECTID',legend=False,color='white', edgecolor='k', alpha=0.5)
       
       data_3.to_crs(32630).plot(ax=ax, column='gs', cmap='Oranges', legend=True,norm=plt.Normalize(vmin=0, vmax=80, clip=False))
       
       ax.set_title('Field average GS: ' + mean_gs)       
       ax.get_xaxis().set_visible(False)
       ax.get_yaxis().set_visible(False)

       minx, miny, maxx, maxy = data_3.to_crs(32630).total_bounds
       offset= 170
       ax.set_xlim(minx-offset,maxx+offset)
       ax.set_ylim(miny-offset,maxy+offset)
       
       ax.add_artist(ScaleBar(2, units='m',))

       fig = ax.get_figure()
       
       plt.show(block=False)
       
       # add text box
       out_date = tt.strftime('%d %B, %Y')
       # NB
       #ax.text(517148,6210330, out_date, fontsize=13, bbox=dict(facecolor='White', alpha=0.8)) 
       
       # WF
       ax.text(515010,6206161, out_date, fontsize=13, bbox=dict(facecolor='White', alpha=0.8))
       
       f.tight_layout()
             
       out_dir = 'C:/Users/arevi/OneDrive/BBSRC_IAA_project/farmer_visit/IMAGES/WF/'
       plt.savefig(out_dir+'gs_'+str(tt)[0:10]+'.png')   
 
       plt.close(fig)


fp_in = out_dir+'gs_*.png'
fp_out = out_dir+"gs_animation_v1.gif"
 
img, *imgs = [Image.open(f) for f in sorted(glob.glob(fp_in))]

img.save(fp=fp_out, format='GIF', append_images=imgs,
         save_all=True, duration=180, loop=0)


# for tt in date_array[range(len(date_array))]:
#        print (tt)
#        data_3 = data[data.date==tt]

#        # data_2.to_file('doy1.shp')
#        f, ax = plt.subplots(figsize=(8,8))
       
#        background.to_crs(32630).plot(ax=ax,column='OBJECTID',legend=False, alpha=1.0,color='white', edgecolor='k')
#        data_3.to_crs(32630).plot(ax=ax, column='lai_float', cmap='RdYlGn',legend=True,norm=plt.Normalize(vmin=0, vmax=5, clip=False))

#        ax.legend(['A simple line'])
#        ax.get_xaxis().set_visible(False)
#        ax.get_yaxis().set_visible(False)

#        minx, miny, maxx, maxy = data_3.to_crs(32630).total_bounds
#        offset= 170
#        ax.set_xlim(minx-offset,maxx+offset)
#        ax.set_ylim(miny-offset,maxy+offset)
#        # ax.set_xlim(-2.260,-2.230)
#        # ax.set_ylim(55.714, 55.734)

#        ax.add_artist(ScaleBar(2, units='m',))

#        fig = ax.get_figure()
       
#        plt.show(block=False)
       
#        # add text box
#        out_date = tt.strftime('%d %B, %Y')
#        ax.text(546691, 6176440, out_date, fontsize=13, bbox=dict(facecolor='White', alpha=0.8)) 
       
#        f.tight_layout()
       
#        out_dir = 'C:/Users/arevi/OneDrive/ATEC/PAPERS/PAPER-5_N_field/vect_animate/LAI_IMAGES/'
#        plt.savefig(out_dir+'lai_'+str(tt)[0:10]+'.png')
 
#        plt.close(fig)


# fp_in = out_dir+'lai_*.png'
# fp_out = out_dir+"LAI_animation_v1.gif"
 
# img, *imgs = [Image.open(f) for f in sorted(glob.glob(fp_in))]

# img.save(fp=fp_out, format='GIF', append_images=imgs,
#          save_all=True, duration=200, loop=0)
