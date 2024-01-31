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
from datetime import datetime as dtdt
##############################################################################
# First, load an AIA observation.

wvl='211'
files=os.listdir()
files=np.sort(np.array([i for i in files if '.'+wvl+'.image' in i and '._' not in i]))[100:]

times=np.array([dtdt(int(i[17:21]),int(i[22:24]),int(i[25:27]),int(i[28:30]),int(i[30:32]),int(i[32:34])) for i in files])
time_exp=dtdt(2023,3,20,3,30,0)
i=np.where(np.abs(times-time_exp)==np.min(np.abs(times-time_exp)))[0][0]
ftype='.png'

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
di=1
#for i in range(di,len(files),di):
map_i=sunpy.map.Map(files[i])
wcsi=map_i.wcs
with propagate_with_solar_surface():
    center_i=center.transform_to(map_i.coordinate_frame)
x0_i=center_i.Tx/u.arcsec-dx/2
y0_i=center_i.Ty/u.arcsec-dy/2
bottom_left_i = SkyCoord(x0_i * u.arcsec,
                        y0_i * u.arcsec,
                        frame=map_i.coordinate_frame)
top_right_i = SkyCoord((x0_i+dx) * u.arcsec,
                        (y0_i+dy) * u.arcsec,
                        frame=map_i.coordinate_frame)
map_i=map_i.submap(bottom_left_i,top_right=top_right_i)
map_i.data[map_i.data<1]=1
map_i.data[np.isnan(map_i.data)]=1

fig=plt.figure(figsize=(7,(7+0.5)/dx*dy),facecolor='white')
plt.rcParams.update({"font.size":10,'font.family':"Arial",'font.weight':'bold'})
map_i.plot(norm=colors.LogNorm(vmin=10, vmax=1000))
plt.xlabel('Solar-X',fontsize=10,fontweight='bold')
plt.ylabel('Solar-Y',fontsize=10,fontweight='bold')
plt.subplots_adjust(left=0.1,bottom=0.08,top=0.98,right=0.98,wspace=0,hspace=0)
plt.gca().grid(False)
plt.title(map_i.name.replace('.0 Angstrom',r'${\rm \AA}$').replace('2023-',''),x=0.25,y=0.0,color='white',fontsize=18,fontweight='bold')#.replace('.0 Angstrom',r' ${\rm \AA}$')
plt.savefig('./'+wvl+'origin/'+map_i.name[-19:]+ftype,dpi=100)
plt.close()
print(i)
#break