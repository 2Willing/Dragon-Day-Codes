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

###############################################################################
# First, download some unprocessed LASCO C2 data with `~sunpy.net.Fido`.
replace=1
files_all=os.listdir()
files=np.sort(np.array([i for i in files_all if '.fts' in i and '._' not in i]))
for i in range(339-6,len(files)-6):
    lasco_map = Map(files[i])
    lasco_map2=Map(files[i+6])
    if replace==0 and 'COR2A_'+str(lasco_map2.date)+'.png' in files_all:
        continue
    #lasco_map.data[np.isnan(lasco_map.data)]=0
    #lasco_map = lasco_map.rotate(angle=float(lasco_map.fits_header['CROTA2'])*u.deg)

    ###############################################################################
    # The LASCO C2 coronagraph has a field of view extending from 2-6 solar
    # radii. So, our mask will have two parts: an inner component which masks the
    # occulter and an outer component which masks data outside the field of view.
    #
    # We will follow a process similar to
    # :ref:`sphx_glr_generated_gallery_computer_vision_techniques_finding_masking_bright_pixels.py`
    # to express the coordinates relative to the occulter center.

    pixel_coords = all_coordinates_from_map(lasco_map)
    solar_center = SkyCoord(0*u.deg, 0*u.deg, frame=lasco_map.coordinate_frame)
    pixel_radii = np.sqrt((pixel_coords.Tx-solar_center.Tx)**2 +
                        (pixel_coords.Ty-solar_center.Ty)**2)
    # Note that the inner mask extends just beyond 2 solar radii to mask the
    # Fresnel diffraction caused by the occulter edge.
    #mask_inner = pixel_radii < lasco_map.rsun_obs*2.4
    #mask_outer = pixel_radii > lasco_map.rsun_obs*6
    #final_mask = mask_inner + mask_outer

    ###############################################################################
    # To apply the final mask, we must create a new map.

    #masked_lasco = Map(lasco_map.data, lasco_map.meta, mask=final_mask)
    #masked_lasco=lasco_map2
    expt1=lasco_map.fits_header['EXPTIME']/50
    expt2=lasco_map2.fits_header['EXPTIME']/50
    masked_lasco=Map(lasco_map2.data.astype(float)/expt2-lasco_map.data.astype(float)/expt1,lasco_map2.meta)

    #masked_lasco=Map(lasco_map2.data.astype(float)[:]-lasco_map.data.astype(float)[:],lasco_map2.meta)
    # Before plotting the map, we need to create a new colormap to ensure we mask
    # the bad values correctly.
    occult_colormap = lasco_map.cmap.copy()
    occult_colormap.set_bad('black')

    fig = plt.figure(figsize=(9.285714285714285,9.285714285714285))
    masked_lasco=masked_lasco.rotate(order=3)
    #ax1 = fig.add_subplot(1, 2, 1, projection=lasco_map)
    ax2 = fig.add_subplot(1, 1, 1, projection=masked_lasco)

    #lasco_map.plot(clip_interval=(2, 98)*u.percent, cmap=occult_colormap, axes=ax1)
    #lasco_map.draw_limb()
    #ax1.set_title("Level 1 LASCO C2")

    masked_lasco.plot(norm=colors.Normalize(vmin=-100, vmax=100), cmap='afmhot', axes=ax2)
    masked_lasco.draw_limb()
    ax2.set_title(masked_lasco.name.replace('Orange white-light ',''))
    ax2.grid(False)
    ax2.set_xlabel('Solar X')
    ax2.set_ylabel('Solar_Y')
    plt.subplots_adjust(left=0.08,bottom=0.07,top=0.95,right=1,wspace=0,hspace=0)

    plt.savefig('COR2A_'+str(masked_lasco.date)+'.png',dpi=100)
    plt.close()
    
    
    fig = plt.figure(figsize=(12,8),facecolor='white')
    plt.rcParams.update({"font.size":20,'font.family':"Arial",'font.weight':'bold'})
    #ax1 = fig.add_subplot(1, 2, 1, projection=lasco_map)
    ax2 = fig.add_subplot(1, 1, 1, projection=masked_lasco)

    masked_lasco.plot(norm=colors.Normalize(vmin=-50, vmax=50), cmap='afmhot', axes=ax2)
    masked_lasco.draw_limb()
    #earth = get_body_heliographic_stonyhurst('earth', masked_lasco.date, observer=masked_lasco.observer_coordinate)
    #ax2.plot_coord(earth, '*', color='green', markersize=25, label='Earth')
    ax2.set_title(masked_lasco.name.replace('white-light','Rdiff.'),fontsize=25,fontweight='bold',y=0.02,x=0.42,color='cyan')
    ax2.grid(False)
    ax2.set_xlabel('Solar-X',fontsize=20,fontweight='bold')
    ax2.set_ylabel('Solar-Y',fontsize=20,fontweight='bold')
    ax2.legend()
    plt.subplots_adjust(left=0.1,bottom=0.1,top=0.99,right=0.99,wspace=0,hspace=0)

    plt.savefig('COR2A_'+str(masked_lasco.date)+'.pdf',dpi=100)
    plt.show()
    #plt.show()
    plt.close()
    break
