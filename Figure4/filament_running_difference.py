"""
=============================
Differentially rotating a map
=============================

How to apply differential rotation to a Map.

.. note::
   This example requires `reproject` 0.6 or later to be installed.

The example uses the :func:`~sunpy.coordinates.propagate_with_solar_surface`
context manager to apply differential rotation during coordinate
transformations.
"""
import matplotlib.pyplot as plt
import matplotlib.colors as colors

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS

import sunpy.map
from sunpy.coordinates import Helioprojective, propagate_with_solar_surface
import os
import numpy as np
from sunpy.physics.differential_rotation import diff_rot, solar_rotate_coordinate
from skimage import filters
from datetime import datetime as dtdt
##############################################################################
# First, load an AIA observation.

files=os.listdir()
files=np.sort(np.array([i for i in files if '.304.image' in i and '._' not in i]))[150:]

times=np.array([dtdt(int(i[17:21]),int(i[22:24]),int(i[25:27]),int(i[28:30]),int(i[30:32]),int(i[32:34])) for i in files])
time_exp=dtdt(2023,3,20,0,54,0)
i=np.where(np.abs(times-time_exp)==np.min(np.abs(times-time_exp)))[0][0]
ftype='.pdf'
filter=0
triple=0

x0=-250
y0=-600
dx=750
dy=400

map0=sunpy.map.Map(files[0])
aia_bottom_left = SkyCoord(x0 * u.arcsec,
                           y0 * u.arcsec,
                           frame=map0.coordinate_frame)
aia_top_right = SkyCoord((x0+dx) * u.arcsec,
                         (y0+dy) * u.arcsec,
                         frame=map0.coordinate_frame)
center=SkyCoord((x0+dx/2)*u.arcsec,(y0+dy/2)*u.arcsec,frame=map0.coordinate_frame)

#map0c=map0.submap(aia_bottom_left,top_right=aia_top_right)
#wcs0=map0.wcs
di=20
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
#for i in [29+30+50+di]:#range(di+30+50,len(files),di):
map0m1=sunpy.map.Map(files[i-di-1])
map0=sunpy.map.Map(files[i-di])
map0p1=sunpy.map.Map(files[i-di+1])
map_im1=sunpy.map.Map(files[i-1])
map_i=sunpy.map.Map(files[i])
map_ip1=sunpy.map.Map(files[i+1])
bottom_left_i = SkyCoord((x0-100) * u.arcsec,
                        (y0-100) * u.arcsec,
                        frame=map0.coordinate_frame)
top_right_i = SkyCoord((x0+dx+100) * u.arcsec,
                        (y0+dy+100) * u.arcsec,
                        frame=map0.coordinate_frame)
map0m1=map0m1.submap(bottom_left_i,top_right=top_right_i)
map0=map0.submap(bottom_left_i,top_right=top_right_i)
map0p1=map0p1.submap(bottom_left_i,top_right=top_right_i)
map0data=(map0.data[:]/map0.fits_header['EXPTIME']*3+map0m1.data[:]/map0m1.fits_header['EXPTIME']*3+map0p1.data[:]/map0p1.fits_header['EXPTIME']*3).astype(float)/3
map0=sunpy.map.Map(map0data,map0.meta)
wcsi=map_i.wcs
with propagate_with_solar_surface():
    map0_i = map0.reproject_to(wcsi)
    center_i=center.transform_to(map_i.coordinate_frame)
x0_i=center_i.Tx.arcsec-dx/2
y0_i=center_i.Ty.arcsec-dy/2
bottom_left_i = SkyCoord(x0_i * u.arcsec,
                        y0_i * u.arcsec,
                        frame=map_i.coordinate_frame)
top_right_i = SkyCoord((x0_i+dx) * u.arcsec,
                        (y0_i+dy) * u.arcsec,
                        frame=map_i.coordinate_frame)
#
map_i=map_i.submap(bottom_left_i,top_right=top_right_i)
#
map_idata=map_i.data.astype(float)/map_i.fits_header['EXPTIME']*3
if triple==1:
    map_idata*=0
    map_ip1=map_ip1.submap(bottom_left_i,top_right=top_right_i)
    map_im1=map_im1.submap(bottom_left_i,top_right=top_right_i)
    map_idata[:670,:1252]=(map_im1.data[:670,:1252]/map_im1.fits_header['EXPTIME']*3+map_i.data[:670,:1252]/map_i.fits_header['EXPTIME']*3+map_ip1.data[:670,:1252]/map_ip1.fits_header['EXPTIME']*3).astype(float)/3
map_i=sunpy.map.Map(map_idata,map_i.meta)
map0c=map0_i.submap(bottom_left_i,top_right=top_right_i)
map_i.data[:670,:1252]=map_i.data[:670,:1252]-map0c.data[:670,:1252]
if filter==1:
    map_i.data[:670,:1252]=map_i.data[:670,:1252]-map0c.data[:670,:1252]
    min_ldata=map_i.data[:670,:1252].min()
    ldata=(map_i.data[:670,:1252]-min_ldata)
    f=np.fft.fft2(ldata)
    f=np.fft.fftshift(f)
    filt = disk(670,1252, 40)
    filt = filters.gaussian(filt*1.0, 10)
    z=np.fft.ifft2(np.fft.fftshift(f*filt))
    z=(np.abs(z))+min_ldata
    map_i.data[:670,:1252]=z[:]










fig=plt.figure(figsize=(7,(7+0.5)/dx*dy),facecolor='white')
plt.rcParams.update({"font.size":10,'font.family':"Arial",'font.weight':'bold'})
ax=plt.subplot(projection=map_i)
map_i.plot(norm=colors.Normalize(vmin=-8, vmax=8))
plt.xlabel('Solar-X',fontsize=10,fontweight='bold')
plt.ylabel('Solar-Y',fontsize=10,fontweight='bold')
plt.subplots_adjust(left=0.1,bottom=0.08,top=0.98,right=0.98,wspace=0,hspace=0)
plt.gca().grid(False)

small_x0=x0_i+dx/2+140-20
small_y0=y0_i+dy/2-30-20
dxy=150
bottom_left_i = SkyCoord((small_x0) * u.arcsec,
                        (small_y0) * u.arcsec,
                        frame=map_i.coordinate_frame)
top_right_i = SkyCoord((small_x0+dxy) * u.arcsec,
                        (small_y0+dxy) * u.arcsec,
                        frame=map_i.coordinate_frame)

coords = SkyCoord(
    Tx=(small_x0, small_x0+dxy) * u.arcsec,
    Ty=(small_y0, small_y0+dxy) * u.arcsec,
    frame=map_i.coordinate_frame,
)
map_i.draw_quadrangle(
    coords,
    axes=plt.gca(),
    edgecolor="cyan",
    linestyle="-",
    linewidth=1,
    label='2-element SkyCoord'
)

if filter==1:
    plt.title(map_i.name.replace('AIA 304.0 Angstrom',r'304 ${\rm \AA}$ Filtered Rdiff').replace('2023-',''),x=0.35,y=0.0,color='cyan',fontsize=18,fontweight='bold')#.replace('.0 Angstrom',r' ${\rm \AA}$')
    plt.savefig('./304rundiff_filter/filter_'+map_i.name[-19:]+ftype,dpi=100)
else:
    plt.title(map_i.name.replace('AIA 304.0 Angstrom',r'304 ${\rm \AA}$ Rdiff.').replace('2023-',''),x=0.28,y=0.0,color='cyan',fontsize=18,fontweight='bold')#.replace('.0 Angstrom',r' ${\rm \AA}$')
    plt.savefig('./304rundiff_filter/'+map_i.name[-19:]+ftype,dpi=100)
plt.close()

map_ic=map_i.submap(bottom_left_i,top_right=top_right_i)
fig=plt.figure(figsize=(map_ic.data.shape[1]/100,map_ic.data.shape[0]/100))
plt.imshow(map_ic.data[::-1],cmap='sdoaia304',vmin=-8,vmax=8)
plt.gca().axis(None)
plt.subplots_adjust(left=0,bottom=0,top=1,right=1,wspace=0,hspace=0)
plt.savefig('./304rundiff_filter/ampview_'+map_i.name[-19:]+'.png',dpi=100)
plt.close()
#break