"""
=================================
Creating a mask for cor2 C2 data
=================================

In this example, we will manually create a mask to block the occulter in an
unprocessed cor2 C2 coronagraph image.
"""
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np

import astropy.units as u
from astropy.coordinates import SkyCoord

from sunpy.map import Map
from sunpy.map.maputils import all_coordinates_from_map
from sunpy.net import Fido
from sunpy.net import attrs as a
from sunpy.coordinates import Helioprojective
import os
from datetime import datetime as dtdt
from datetime import timedelta as dttd
from skimage import filters

import warnings

warnings.filterwarnings("error")
###############################################################################
# First, download some unprocessed cor2 C2 data with `~sunpy.net.Fido`.
start=1
diff=2
shift=2.53*diff/2.
dir_cor2='./cor2/'
dir_hi1='./hi_1/'
dir_hi2='./hi_2/'
files_cor2=os.listdir(dir_cor2)
files_cor2=np.sort(np.array([i for i in files_cor2 if '.fts' in i and '._' not in i]))
files_hi1=os.listdir(dir_hi1)
files_hi1=np.sort(np.array([i for i in files_hi1 if '.fts' in i and '._' not in i]))
files_hi2=os.listdir(dir_hi2)
files_hi2=np.sort(np.array([i for i in files_hi2 if '.fts' in i and '._' not in i]))

times_cor2=np.array([dtdt(int(i[0:4]),int(i[4:6]),int(i[6:8]),int(i[9:11]),int(i[11:13]),int(i[13:15])) for i in files_cor2])
times_hi1=np.array([dtdt(int(i[0:4]),int(i[4:6]),int(i[6:8]),int(i[9:11]),int(i[11:13]),int(i[13:15])) for i in files_hi1])
times_hi2=np.array([dtdt(int(i[0:4]),int(i[4:6]),int(i[6:8]),int(i[9:11]),int(i[11:13]),int(i[13:15])) for i in files_hi2])

time_exp=dtdt(2023,3,19,16,0,0)

def removenan(im,key=0):
    im2 = np.copy(im)
    arr = np.isnan(im2)
    im2[arr] = key
    arr2 = np.isinf(im2)
    im2[arr2] = key
def disk(M, N, r0):
    X, Y = np.meshgrid(np.arange(int(-(N / 2)), int(N / 2)),
                    np.linspace(-int(M / 2), int(M / 2) - 1, M))
    r = (X) ** 2 + (Y) ** 2
    r = (r ** 0.5)
    im = r < r0
    return im
def movdiff(map1,map2,shift):
    if shift==0:
        return map2-map1
    elif shift>0:
        shift_int=int(np.ceil(shift))
        shift_res=float(shift_int)-float(shift)
        dmap=np.zeros(map2.shape)
        if -shift_int+1==0:
            shift_end=len(map1[0])
        else:
            shift_end=-shift_int+1
        dmap[:,shift_int:]=map2[:,shift_int:]-map1[:,:-shift_int]*(1-shift_res)-map1[:,1:shift_end]*shift_res
    else:
        shift_int=abs(int(np.floor(shift)))
        shift_res=shift-shift_int
        dmap=np.zeros(map2.shape)
        dmap[:,:-shift_int]=map2[:,:-shift_int]-map1[:,shift_int:]*(1-shift_res)-map1[:,shift_int-1:-1]*shift_res
    return dmap


while time_exp.day < 22:
    time_exp=time_exp+dttd(minutes=15)
    i_cor2=np.where(np.abs(times_cor2-time_exp)==np.min(np.abs(times_cor2-time_exp)))[0][0]-6
    i_hi1=np.where(np.abs(times_hi1-time_exp)==np.min(np.abs(times_hi1-time_exp)))[0][0]-6
    i_hi2=np.where(np.abs(times_hi2-time_exp)==np.min(np.abs(times_hi2-time_exp)))[0][0]-diff
    cor2_map = Map(dir_cor2+files_cor2[i_cor2])
    cor2_map2=Map(dir_cor2+files_cor2[i_cor2+6])
    expt1=cor2_map.fits_header['EXPTIME']/50
    expt2=cor2_map2.fits_header['EXPTIME']/50
    masked_cor2=Map(cor2_map2.data.astype(float)/expt2-cor2_map.data.astype(float)/expt1,cor2_map2.meta)
    masked_cor2=masked_cor2.rotate(order=3)
    #masked_cor2=Map(cor2_map2.data.astype(float)[:]-cor2_map.data.astype(float)[:],cor2_map2.meta)
    # Before plotting the map, we need to create a new colormap to ensure we mask
    # the bad values correctly.
    
    
    
    
    
    
    
    
    
    
    
    hi1_map = Map(dir_hi1+files_hi1[i_hi1])
    hi1_map2=Map(dir_hi1+files_hi1[i_hi1+6])
    expt1=hi1_map.fits_header['EXPTIME']/1200.
    expt2=hi1_map2.fits_header['EXPTIME']/1200.
    masked_hi1=Map(hi1_map2.data.astype(float)[:]/expt2-hi1_map.data.astype(float)[:]/expt1,hi1_map2.meta)
    masked_hi1=masked_hi1.rotate(order=3)
    
    
    
    
    #map_comp=Map(np.array(masked_cor2,masked_hi1,masked_hi2),composite=True)
    
    fits_cor2=masked_cor2.wcs.to_header()
    minx_cor2=fits_cor2['CRVAL1']-fits_cor2['CDELT1']*(fits_cor2['CRPIX1'])
    maxx_cor2=fits_cor2['CRVAL1']+fits_cor2['CDELT1']*(len(masked_cor2.data[0])-fits_cor2['CRPIX1'])
    maxy_cor2=fits_cor2['CRVAL2']+fits_cor2['CDELT2']*(fits_cor2['CRPIX2'])
    miny_cor2=fits_cor2['CRVAL2']-fits_cor2['CDELT2']*(len(masked_cor2.data)-fits_cor2['CRPIX2'])
    
    fits_hi1=masked_hi1.wcs.to_header()
    minx_hi1=fits_hi1['CRVAL1']-fits_hi1['CDELT1']*(fits_hi1['CRPIX1'])
    maxx_hi1=fits_hi1['CRVAL1']+fits_hi1['CDELT1']*(len(masked_hi1.data[0])-fits_hi1['CRPIX1'])
    maxy_hi1=fits_hi1['CRVAL2']+fits_hi1['CDELT2']*(fits_hi1['CRPIX2'])
    miny_hi1=fits_hi1['CRVAL2']-fits_hi1['CDELT2']*(len(masked_hi1.data)-fits_hi1['CRPIX2'])
    
    
    minx=-6
    miny=-12
    
    print('generating_figure...')
    fullx=32
    fully=24
    
    fig = plt.figure(figsize=(fullx/10,fully/10),facecolor='black')
    #ax1 = fig.add_subplot(1, 3, 1, projection=masked_cor2)
    ax1=plt.axes([(minx_cor2-minx)/fullx, (miny_cor2-miny)/fully, (maxx_cor2-minx_cor2)/fullx, (maxy_cor2-miny_cor2)/fully],projection=masked_cor2)
    #cor2_map.plot(clip_interval=(2, 98)*u.percent, cmap=occult_colormap, axes=ax1)
    #cor2_map.draw_limb()
    #ax1.set_title("Level 1 cor2 C2")
    #map_comp.set_plot_settings(index=0,plot_settings={'norm':colors.Normalize(vmin=-100, vmax=100),'cmap':'afmhot'})
    #map_comp.set_plot_settings(index=1,plot_settings={'norm':colors.Normalize(vmin=-200, vmax=200),'cmap':masked_hi1.cmap})
    #map_comp.set_plot_settings(index=2,plot_settings={'norm':colors.Normalize(vmin=-100, vmax=100),'cmap':masked_hi2.cmap})
    #with Helioprojective.assume_spherical_screen(masked_cor2.observer_coordinate):
        #masked_hi1 = masked_hi1.reproject_to(masked_cor2.wcs)
        #masked_hi2 = masked_hi2.reproject_to(masked_cor2.wcs)
    masked_cor2.plot(norm=colors.Normalize(vmin=-100, vmax=100),cmap='afmhot',zorder=1)
        #map_comp.draw_limb()
    #with Helioprojective.assume_spherical_screen(masked_hi1.observer_coordinate):
    plt.axis('off')
    #ax1.set_position([(minx_cor2+minx)/fullx, (miny_cor2+miny)/fully, (maxx_cor2-minx_cor2)/fullx, (maxy_cor2-miny_cor2)/fully])
    plt.grid(True, alpha=0)
    plt.title(None)
    #ax2 = fig.add_subplot(1,3,2,projection=masked_hi1)
    ax2=plt.axes([(minx_hi1-minx)/fullx, (miny_hi1-miny)/fully, (maxx_hi1-minx_hi1)/fullx, (maxy_hi1-miny_hi1)/fully],projection=masked_hi1,zorder=3)
    masked_hi1.plot(norm=colors.Normalize(vmin=-150, vmax=150), autoalign=False, axes=ax2)
    #with Helioprojective.assume_spherical_screen(masked_hi2.observer_coordinate):
    plt.axis('off')
    plt.grid(True, alpha=0)
    plt.title(None)
    #ax3 = fig.add_subplot(1,3,3,projection=masked_hi2)
    #plt.xlabel('Solar_X')
    #plt.ylabel('Solar_Y')
    #plt.subplots_adjust(left=0.08,bottom=0.07,top=0.95,right=1,wspace=0,hspace=0)
    plt.subplots_adjust(left=0,bottom=0,top=1,right=1,wspace=0,hspace=0)
    ax4=plt.axes([0,0,1,1],facecolor='none')
    plt.text(0.01,0.02,str(time_exp),fontsize=20/100*fullx,color='white',fontweight='bold')
    
    plt.savefig('./composite_small/'+str(time_exp)+'.png',dpi=300)
    plt.close()
    #break
