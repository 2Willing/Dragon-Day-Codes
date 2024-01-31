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
#from sunpy.map.maputils import all_coordinates_from_map
from sunpy.coordinates import propagate_with_solar_surface
from sunpy.net import Fido
from sunpy.net import attrs as a
import os

###############################################################################
# First, download some unprocessed LASCO C2 data with `~sunpy.net.Fido`.

aia_map=Map('/Volumes/金乌3号/JSOC_20231128_2155/aia.lev1_euv_12s.2023-03-20T045810Z.171.image_lev1.fits')
x0=-500
y0=-600
dx=1000
dy=1300

files=os.listdir()
files=np.sort(np.array([i for i in files if '.fits' in i and '._' not in i]))

so_map=Map(files[0])
so_map=so_map.rotate(order=3)
xx0=-600
dxx=2200
yy0=-1300
dyy=2860
center=SkyCoord((xx0+dxx/2)*u.arcsec,(yy0+dyy/2)*u.arcsec,frame=so_map.coordinate_frame)

for i in range(len(files)):
    lasco_map = Map(files[i])
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

    '''pixel_coords = all_coordinates_from_map(lasco_map)
    solar_center = SkyCoord(0*u.deg, 0*u.deg, frame=lasco_map.coordinate_frame)
    pixel_radii = np.sqrt((pixel_coords.Tx-solar_center.Tx)**2 +
                        (pixel_coords.Ty-solar_center.Ty)**2)'''
    # Note that the inner mask extends just beyond 2 solar radii to mask the
    # Fresnel diffraction caused by the occulter edge.
    #mask_inner = pixel_radii < lasco_map.rsun_obs*2.4
    #mask_outer = pixel_radii > lasco_map.rsun_obs*6
    #final_mask = mask_inner + mask_outer

    ###############################################################################
    # To apply the final mask, we must create a new map.

    lasco_map=lasco_map.rotate(order=3)
    with propagate_with_solar_surface():
        center_i=center.transform_to(lasco_map.coordinate_frame)
    xx0_i=center_i.Tx/u.arcsec-dxx/2
    yy0_i=center_i.Ty/u.arcsec-dyy/2
    bottom_left = np.array([xx0_i, yy0_i]) * u.arcsec
    top_right = np.array([xx0_i+dxx, yy0_i+dyy]) * u.arcsec

    lasco_map = lasco_map.submap(SkyCoord(*bottom_left, frame=lasco_map.coordinate_frame),
                          top_right=SkyCoord(*top_right, frame=lasco_map.coordinate_frame))
    masked_lasco = Map(lasco_map.data/lasco_map.fits_header['XPOSURE']*10.+1, lasco_map.meta)

    # Before plotting the map, we need to create a new colormap to ensure we mask
    # the bad values correctly.

    fig=plt.figure(figsize=(5+0.6,5/dx*dy),facecolor='white') 
    #ax1 = fig.add_subplot(1, 2, 1, projection=lasco_map)
    #ax2 = fig.add_subplot(1, 1, 1, projection=masked_lasco)

    #lasco_map.plot(clip_interval=(2, 98)*u.percent, cmap=occult_colormap, axes=ax1)
    #lasco_map.draw_limb()
    #ax1.set_title("Level 1 LASCO C2")
    #masked_lasco.plot(norm=colors.LogNorm(vmin=1, vmax=5000), cmap=occult_colormap, axes=ax2)
    plt.rcParams.update({"font.size":15,'font.family':"Arial",'font.weight':'bold'})
    masked_lasco.data[:]=masked_lasco.data[:]
    masked_lasco.data[masked_lasco.data<0]=0
    masked_lasco.data[:]=masked_lasco.data[:]+1
    masked_lasco.plot(norm=colors.LogNorm(vmin=100, vmax=6000))
    plt.xlabel('Solar-X',fontsize=15,fontweight='bold')
    plt.ylabel('Solar-Y',fontsize=15,fontweight='bold')
    #plt.gca().tick_params(axis="x", direction="in", which="major", length=7,width=1.5)
    #plt.gca().tick_params(axis="y", direction="in", which="major", length=7,width=1.5)
    plt.subplots_adjust(left=0.15,bottom=0.1,top=0.98,right=1,wspace=0,hspace=0)
    plt.title(masked_lasco.name.replace('.0 Angstrom',r' ${\rm \AA}$').replace('2023-',''),y=0.92,color='white',fontsize=22,fontweight='bold')
    if i==0:
        with propagate_with_solar_surface():
            coords = SkyCoord(
                Tx=(x0, x0+dx) * u.arcsec,
                Ty=(y0, y0+dy) * u.arcsec,
                frame=aia_map.coordinate_frame,
            )
            masked_lasco.draw_quadrangle(
                coords,
                axes=plt.gca(),
                edgecolor="white",
                linestyle="--",
                linewidth=2,
                label='2-element SkyCoord'
            )
    plt.grid(True,alpha=0.)

    plt.savefig('EUI_fsi174_'+str(masked_lasco.date)+'.pdf',dpi=200)
    plt.close()
    print(i)
    #break
