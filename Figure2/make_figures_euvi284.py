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
#from proplot import rc
##############################################################################
# First, load an AIA observation.


wvl='304'
typef='.pdf'
dpi=200
files=os.listdir()
files=np.sort(np.array([i for i in files if '.'+wvl+'.image' in i and '._' not in i]))

dir_284='/Volumes/金乌3号/JSOC_20231128_2155/euvi284/'
files_284=os.listdir(dir_284)
files_284=np.sort(np.array([i for i in files_284 if '.fts' in i and '._' not in i]))
rd=1
times=np.array([dtdt(int(i[0:4]),int(i[4:6]),int(i[6:8]),int(i[9:11]),int(i[11:13]),int(i[13:15])) for i in files_284])
time_exp=dtdt(2023,3,19,20,30,0)
i=np.where(np.abs(times-time_exp)==np.min(np.abs(times-time_exp)))[0][0]

x0=-500
y0=-600
dx=1000
dy=1300

map0=sunpy.map.Map(files[0])
center=SkyCoord((x0+dx/2)*u.arcsec,(y0+dy/2)*u.arcsec,frame=map0.coordinate_frame)

#map0c=map0.submap(aia_bottom_left,top_right=aia_top_right)
#wcs0=map0.wcs
if rd==0:
    map_i=sunpy.map.Map(dir_284+files_284[i])
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
    fig=plt.figure(figsize=(5+0.6,5/dx*dy),facecolor='white')
    plt.rcParams.update({"font.size":15,'font.family':"Arial",'font.weight':'bold'})
    map_i.data[:]=map_i.data[:]-730
    map_i.data[map_i.data<0]=0
    map_i.data[:]=map_i.data[:]+1
    map_i.plot(norm=colors.LogNorm(vmin=5, vmax=2400))
    plt.xlabel('Solar-X',fontsize=15,fontweight='bold')
    plt.ylabel('Solar-Y',fontsize=15,fontweight='bold')
    #plt.gca().tick_params(axis="x", direction="in", which="major", length=7,width=1.5)
    #plt.gca().tick_params(axis="y", direction="in", which="major", length=7,width=1.5)
    plt.subplots_adjust(left=0.15,bottom=0.1,top=0.98,right=1,wspace=0,hspace=0)
    plt.title(map_i.name.replace('.0 Angstrom',r' ${\rm \AA}$').replace('2023-',''),y=0.92,color='white',fontsize=22,fontweight='bold')#,y=0.92,color='white',fontsize=18,fontweight='bold')
    plt.gca().grid(False)
    plt.savefig('./Figure2/'+map_i.name+typef,dpi=dpi)
    plt.close()
else:
    map_i=sunpy.map.Map(dir_284+files_284[i])
    map284_0=sunpy.map.Map(dir_284+files_284[0])
    wcsi=map_i.wcs
    with propagate_with_solar_surface():
        map0_i = map284_0.reproject_to(wcsi)
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
    map0c=map0_i.submap(bottom_left_i,top_right=top_right_i)
    map_i.data[map_i.data<0]=0
    #map_i.data[:]=map_i.data[:]-map0c.data[:]
    
    map_i=sunpy.map.Map(map_i.data[:]-map0c.data[:],map_i.meta)

    fig=plt.figure(figsize=(5+0.6,5/dx*dy),facecolor='white')
    plt.rcParams.update({"font.size":15,'font.family':"Arial",'font.weight':'bold'})
    
    map_i.plot(norm=colors.Normalize(vmin=-100, vmax=100))
    plt.xlabel('Solar-X',fontsize=15,fontweight='bold')
    plt.ylabel('Solar-Y',fontsize=15,fontweight='bold')
    plt.subplots_adjust(left=0.15,bottom=0.1,top=0.98,right=1,wspace=0,hspace=0)
    plt.title('Base Diff. '+map_i.name.replace('EUVI-A 284.0 Angstrom','').replace('2023-',''),y=0.92,color='cyan',fontsize=22,fontweight='bold')#,y=0.92,color='cyan',fontsize=18,fontweight='bold')#.replace('.0 Angstrom',r' ${\rm \AA}$')
    plt.gca().grid(False)
    plt.savefig('./Figure2/'+map_i.name+typef,dpi=dpi)
    plt.close()
    #break