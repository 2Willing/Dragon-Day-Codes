"""
=================================
Creating a mask for LASCO C2 data
=================================

In this example, we will manually create a mask to block the occulter in an
unprocessed LASCO C2 coronagraph image.
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
import os
from skimage import filters

from sunpy.coordinates import get_body_heliographic_stonyhurst

###############################################################################
# First, download some unprocessed LASCO C2 data with `~sunpy.net.Fido`.
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

files_all=os.listdir()
files=np.sort(np.array([i for i in files_all if '.fts' in i and '._' not in i]))
replace=1
start=1
diff=2
shift=2.53*diff/2.
for i in range(33,len(files)-diff):
    lasco_map = Map(files[i])
    lasco_map2=Map(files[i+diff])
    map01=Map(files[i-start])
    map02=Map(files[i+diff-start])
    ldata01=map01.data.astype(float)/map01.fits_header['EXPTIME']*4950.
    ldata02=map02.data.astype(float)/map02.fits_header['EXPTIME']*4950.
    if replace==0 and 'HI2A_'+str(lasco_map2.date)+'.png' in files_all:
        continue
    #lasco_map.data[np.isnan(lasco_map.data)]=0
    #lasco_map = lasco_map.rotate(angle=float(lasco_map.fits_header['CROTA2'])*u.deg)

    if shift!=0:
        ldata1=lasco_map.data.astype(float)/lasco_map.fits_header['EXPTIME']*4950.-ldata01
        ldata2=lasco_map2.data.astype(float)/lasco_map2.fits_header['EXPTIME']*4950.-ldata02
        dmap=movdiff(ldata1,ldata2,shift)
        masked_lasco=Map(dmap,lasco_map2.meta)
    else:
        masked_lasco=Map(lasco_map2.data.astype(float)-lasco_map.data.astype(float),lasco_map2.meta)
    # Before plotting the map, we need to create a new colormap to ensure we mask
    # the bad values correctly.
    occult_colormap = lasco_map.cmap.copy()
    occult_colormap.set_bad('black')

    #lasco_map.plot(clip_interval=(2, 98)*u.percent, cmap=occult_colormap, axes=ax1)
    #lasco_map.draw_limb()
    #ax1.set_title("Level 1 LASCO C2")
    min_ldata=masked_lasco.data.min()
    ldata=masked_lasco.data-min_ldata
    f=np.fft.fft2(ldata)
    f=np.fft.fftshift(f)
    filt = disk(1024,1024, 25)
    filt = filters.gaussian(filt*1.0, 5)
    z=np.fft.ifft2(np.fft.fftshift(f*filt))
    z=(np.abs(z))+min_ldata
    masked_lasco.data[:]=z[:]
    masked_lasco=masked_lasco.rotate(order=3)
    
    fig = plt.figure(figsize=(1,8),facecolor='white')
    plt.rcParams.update({"font.size":20,'font.family':"Arial",'font.weight':'bold'})
    #ax1 = fig.add_subplot(1, 2, 1, projection=lasco_map)
    ax2 = fig.add_subplot(1, 1, 1, projection=masked_lasco)

    masked_lasco.plot(norm=colors.Normalize(vmin=-1000, vmax=1000), cmap=occult_colormap, axes=ax2)
    masked_lasco.draw_limb()
    earth = get_body_heliographic_stonyhurst('earth', masked_lasco.date, observer=masked_lasco.observer_coordinate)
    ax2.plot_coord(earth, '*', color='green', markersize=25, label='Earth')
    ax2.set_title(masked_lasco.name.replace('white-light ','Rdiff.'),fontsize=25,fontweight='bold',y=0.02,x=0.35,color='cyan')
    ax2.grid(False)
    ax2.set_xlabel('Solar-X',fontsize=20,fontweight='bold')
    ax2.set_ylabel('Solar-Y',fontsize=20,fontweight='bold')
    ax2.legend()
    plt.subplots_adjust(left=0.1,bottom=0.1,top=0.99,right=0.99,wspace=0,hspace=0)

    plt.savefig('./earth/HI2A_'+str(masked_lasco.date)+'.pdf',dpi=100)
    plt.show()
    #plt.show()
    plt.close()
    #03231200:92.55
    #03231800:92.54
    break
