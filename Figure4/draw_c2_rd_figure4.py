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

ftype='.pdf'
files=os.listdir()
files=np.sort(np.array([i for i in files if '.fts' in i and '._' not in i]))
for i in [80-5]:#range(0,len(files)-5):
    lasco_map = Map(files[i])
    lasco_map2=Map(files[i+5])
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

    #masked_lasco = lasco_map#Map(lasco_map.data, lasco_map.meta, mask=final_mask)
    expt1=lasco_map.fits_header['EXPTIME']/50
    expt2=lasco_map2.fits_header['EXPTIME']/50
    masked_lasco=Map(lasco_map2.data.astype(float)/expt2-lasco_map.data.astype(float)/expt1,lasco_map2.meta)

    # Before plotting the map, we need to create a new colormap to ensure we mask
    # the bad values correctly.
    occult_colormap = lasco_map.cmap.copy()
    occult_colormap.set_bad('black')
    
    fig=plt.figure(figsize=(7*3.5/4,(7-0.5)*3.5/4),facecolor='white')
    plt.rcParams.update({"font.size":10,'font.family':"Arial",'font.weight':'bold'})
    masked_lasco=masked_lasco.rotate(order=3)
    masked_lasco.plot(norm=colors.Normalize(vmin=-100, vmax=100))
    plt.xlabel('Solar-X',fontsize=10,fontweight='bold')
    plt.ylabel('Solar-Y',fontsize=10,fontweight='bold')
    plt.subplots_adjust(left=0.12,bottom=0.1,top=0.985,right=0.99,wspace=0,hspace=0)
    plt.gca().grid(False)
    plt.gca().set_facecolor('black')
    plt.title(masked_lasco.name.replace('Orange white-light 2023-','').replace('2023-',''),x=0.4,y=0.01,color='white',fontsize=18,fontweight='bold')#.replace('.0 Angstrom',r' ${\rm \AA}$')
    plt.savefig('/Volumes/金乌3号/JSOC_20231128_2155/figure4/'+masked_lasco.name+ftype,dpi=100)
    plt.close()
