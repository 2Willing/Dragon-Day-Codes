import sunpy
import sunpy.map
import numpy as np
from math import *
import astropy.units as u
from astropy.io import fits
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from sunpy.coordinates import frames
import matplotlib.gridspec as gridspec
from scipy.interpolate import interp1d
from astropy.coordinates import SkyCoord
from matplotlib.patches import ConnectionPatch
from sunpy.coordinates import propagate_with_solar_surface
import matplotlib.colors as colors
import os
filename='/Volumes/金乌3号/风暴/CHASE/RSM20230320T000022_0000_HA.fits'
dir_aia='/Volumes/金乌3号/JSOC_20231128_2155/'
files=os.listdir(dir_aia)
files=np.sort(np.array([dir_aia+i for i in files if '.304.image' in i and '._' not in i]))[100:]
map0=sunpy.map.Map(files[0])

x0=-250
y0=-600
dx=750
dy=400

center=SkyCoord((x0+dx/2)*u.arcsec,(y0+dy/2)*u.arcsec,frame=map0.coordinate_frame)


def rot(array,arc,dxi,dyi):
    #global box2,i1,i2,i11,i22,max11,max12,max21,max22,med1,med2,med11,med22
    box2=np.zeros([int(np.shape(array)[0]*1.5*2)]*2)
    i1=np.arange(-int(np.shape(array)[0]/2),np.shape(array)[0]-int(np.shape(array)[0]/2),0.5)
    i2=np.arange(-int(np.shape(array)[1]/2),np.shape(array)[1]-int(np.shape(array)[1]/2),0.5)
    i1,i2=np.meshgrid(i1,i2)
    i11=(i1+dyi)*np.cos(arc)+(i2+dxi)*np.sin(arc)
    i22=-(i1+dyi)*np.sin(arc)+(i2+dxi)*np.cos(arc)
    #i11=i11*2
    #i22=i22*2
    i11=i11.astype(int)
    i22=i22.astype(int)
    i1=i1.astype(int)
    i2=i2.astype(int)
    i1-=np.min(i1)
    i2-=np.min(i2)
    i11+=int(np.shape(array)[0]*1.5)
    i22+=int(np.shape(array)[0]*1.5)
    box2[i11,i22]=array[i1,i2]
    len1=np.shape(array)[0]
    len2=np.shape(box2)[0]
    rmin=int((len2-len1)/2)
    rmax=int(len1+(len2-len1)/2)
    box2=box2[rmin:rmax,rmin:rmax]
    return box2

def removenan(im,key=0):
    im2 = np.copy(im)
    arr = np.isnan(im2)
    im2[arr] = key
    arr2 = np.isinf(im2)
    im2[arr2] = key
    return im2

rsm = fits.open(filename)
print(rsm[1].data.shape)
hacore = rsm[1].data[71, :, :]
hablue = np.mean(rsm[1].data[71-10:71-4, :, :],axis=0)
hared = np.mean(rsm[1].data[71+4:71+10, :, :],axis=0)
# 0.5218" per pixel in bin = 1 mode and 1.0436" per pixel in bin = 2 mode
time=filename[-28:-24]+'-'+filename[-24:-22]+'-'+filename[-22:-20]+' '+filename[-19:-17]+':'+filename[-17:-15]+':'+filename[-15:-13]
coord_HIS = SkyCoord(0 * u.arcsec, 0 * u.arcsec, obstime = time, observer = 'earth', \
                     frame = frames.Helioprojective)
'''headerwingb = sunpy.map.make_fitswcs_header(hablue, coord_HIS,
                                       reference_pixel = \
                                       [rsm[1].header['CRPIX1']+29, rsm[1].header['CRPIX2']-25] * u.pixel,
                                       scale = [0.5218 * 2, 0.5218 * 2] * u.arcsec / u.pixel,
                                       telescope = 'CHASE', instrument = 'RSM')
headerwingr = sunpy.map.make_fitswcs_header(hared, coord_HIS,
                                       reference_pixel = \
                                       [rsm[1].header['CRPIX1']+29, rsm[1].header['CRPIX2']-25] * u.pixel,
                                       scale = [0.5218 * 2, 0.5218 * 2] * u.arcsec / u.pixel,
                                       telescope = 'CHASE', instrument = 'RSM')'''
hacore=removenan(hacore,0)
hacore=rot(hacore,-rsm[1].header['INST_ROT']/180.*np.pi,0,0)
headercore = sunpy.map.make_fitswcs_header(hacore, coord_HIS,
                                       reference_pixel = \
                                       [rsm[1].header['CRPIX1']+32, rsm[1].header['CRPIX2']-30] * u.pixel,
                                       scale = [0.5218 * 2, 0.5218 * 2] * u.arcsec / u.pixel,
                                       telescope = 'CHASE', instrument = 'RSM')
#hawing_map = sunpy.map.Map(hared-hablue, headerwingr)
hacore_map = sunpy.map.Map(hacore, headercore)
#hacore_map = sunpy.map.Map(hacore_map.data,hacore_map.meta)
#hacore_map=hacore_map.rotate(order=3)
with propagate_with_solar_surface():
    center_i=center.transform_to(hacore_map.coordinate_frame)
    
x0_i=center_i.Tx.arcsec-dx/2
y0_i=center_i.Ty.arcsec-dy/2

'''fig = plt.figure()
fig = plt.figure(figsize = (20, 20))
hawing_map.plot(cmap = 'bwr' , vmin = -100, vmax = 100)
hawing_map.draw_grid(color = 'cyan',linewidth=2)
hawing_map.draw_limb(color='black')
plt.colorbar()
plt.show()'''
fig = plt.figure()
fig = plt.figure(figsize = (20, 20))
hacore_map.plot(cmap = 'afmhot' , vmin = 0, vmax = 4 * hacore.mean())
hacore_map.draw_grid(color = 'white')
hacore_map.draw_limb()
plt.show()
ang_res = 0.5218 * 2

'''fig = plt.figure(figsize = (26, 12))
gs = gridspec.GridSpec(6, 13, wspace = 0, hspace = 1)


ax1 = fig.add_subplot(gs[0:6, 0:5], projection = hacore_map)
hacore_map.plot(axes = ax1, title = '', cmap = 'afmhot', vmin = 0, vmax = 4 * hacore_map.mean())
plt.grid(False)'''

'''left, right, bottom, top = x0_i-100, x0_i+dx+100, y0_i-100, y0_i+dy+100
left_corner = SkyCoord(Tx = left * u.arcsec, Ty = bottom * u.arcsec, frame = hacore_map.coordinate_frame)
right_corner = SkyCoord(Tx = right * u.arcsec, Ty = top * u.arcsec, frame = hacore_map.coordinate_frame)

hacore_map_medium = hacore_map.submap(left_corner, top_right = right_corner)
hacore_map_medium=hacore_map_medium.rotate(order=3)'''



left, right, bottom, top = x0_i, x0_i+dx, y0_i, y0_i+dy
left_corner = SkyCoord(Tx = left * u.arcsec, Ty = bottom * u.arcsec, frame = hacore_map.coordinate_frame)
right_corner = SkyCoord(Tx = right * u.arcsec, Ty = top * u.arcsec, frame = hacore_map.coordinate_frame)

#left2, right2, bottom2, top2 = -700, -400, -100, 200
#left_corner2 = SkyCoord(Tx = left2 * u.arcsec, Ty = bottom2 * u.arcsec, frame = hacore_map.coordinate_frame)
#right_corner2 = SkyCoord(Tx = right2 * u.arcsec, Ty = top2 * u.arcsec, frame = hacore_map.coordinate_frame)

'''len_arc = (right - left) / 100 + 1
tick_arcsecx = list(np.linspace(left, right, int(len_arc)))
tick_arcsecx = [int(i) for i in tick_arcsecx]
tick_arcsecy = list(np.linspace(bottom, top, int(len_arc)))
tick_arcsecy = [int(i) for i in tick_arcsecy]
tick_pixelx = list(np.linspace(0, (right - left) / ang_res, int(len_arc)))
tick_pixely = list(np.linspace(0, (top - bottom) / ang_res, int(len_arc)))'''

'''len_arc2 = (right2 - left2) / 100 + 1
tick_arcsec2x = list(np.linspace(left2, right2, int(len_arc2)))
tick_arcsec2x = [int(i) for i in tick_arcsec2x]
tick_arcsec2y = list(np.linspace(bottom2, top2, int(len_arc2)))
tick_arcsec2y = [int(i) for i in tick_arcsec2y]
tick_pixel2x = list(np.linspace(0, (right2 - left2) / ang_res, int(len_arc2)))
tick_pixel2y = list(np.linspace(0, (top2 - bottom2) / ang_res, int(len_arc2)))'''

'''for coord in ax1.coords:
    coord.frame.set_linewidth(0)
    coord.set_ticks_visible(False)
    coord.set_ticklabel_visible(False)

hacore_map.draw_quadrangle(left_corner, top_right = right_corner, edgecolor = 'black', lw = 1)
hacore_map.draw_quadrangle(left_corner2, top_right = right_corner2, edgecolor = 'black', lw = 1)'''





hacore_map_small = hacore_map.submap(left_corner, top_right = right_corner)
fig=plt.figure(figsize=(7,(7+0.5)/dx*dy),facecolor='white')
plt.rcParams.update({"font.size":10,'font.family':"Arial",'font.weight':'bold'})
hacore_map_small.plot(norm=colors.Normalize(vmin=250, vmax=450),cmap='afmhot')
plt.xlabel('Solar-X',fontsize=10,fontweight='bold')
plt.ylabel('Solar-Y',fontsize=10,fontweight='bold')
plt.subplots_adjust(left=0.1,bottom=0.08,top=0.98,right=0.98,wspace=0,hspace=0)
plt.gca().grid(False)
plt.title(r'CHASE H$\alpha$'+hacore_map_small.name.replace('2023-',''),x=0.28,y=0.0,color='white',fontsize=18,fontweight='bold')#.replace('.0 Angstrom',r' ${\rm \AA}$')
small_x0=x0_i+dx/2+140-20
small_y0=y0_i+dy/2-30-20
dxy=150
bottom_left_i = SkyCoord((small_x0) * u.arcsec,
                        (small_y0) * u.arcsec,
                        frame=hacore_map_small.coordinate_frame)
top_right_i = SkyCoord((small_x0+dxy) * u.arcsec,
                        (small_y0+dxy) * u.arcsec,
                        frame=hacore_map_small.coordinate_frame)

coords = SkyCoord(
    Tx=(small_x0, small_x0+dxy) * u.arcsec,
    Ty=(small_y0, small_y0+dxy) * u.arcsec,
    frame=hacore_map_small.coordinate_frame,
)
hacore_map_small.draw_quadrangle(
    coords,
    axes=plt.gca(),
    edgecolor="cyan",
    linestyle="-",
    linewidth=1,
    label='2-element SkyCoord'
)
plt.savefig(time+'.pdf',dpi=100)
plt.close()
'''hacore_map_small = hacore_map_small.data
ax2 = fig.add_subplot(gs[0:3, 6:9])
im2 = ax2.imshow(hacore_map_small, origin = 'lower',cmap = 'afmhot', vmin = 0, vmax = 4 * hacore_map.mean())
plt.xlabel('Solar X (arcsec)', fontsize = 18)
plt.ylabel('Solar Y (arcsec)', fontsize = 18)
plt.xticks(ticks = tick_pixelx, labels = tick_arcsecx, fontsize = 13)
plt.yticks(ticks = tick_pixely, labels = tick_arcsecy, fontsize = 13)

plt.savefig(time+'.pdf',dpi=100)
plt.close()'''
'''wing_small = hawing_map.submap(left_corner, top_right = right_corner)
wing_small = wing_small.data


ax3 = fig.add_subplot(gs[0:3, 10:13])
im3 = ax3.imshow(wing_small, origin = 'lower', cmap = 'afmhot', vmin = 0, vmax = 4 * hawing_map.mean())
plt.xlabel('Solar X (arcsec)', fontsize = 18)
plt.ylabel('Solar Y (arcsec)', fontsize = 18)
plt.xticks(ticks = tick_pixelx, labels = tick_arcsecx, fontsize = 13)
plt.yticks(ticks = tick_pixely, labels = tick_arcsecy, fontsize = 13)


xpix, ypix = hacore_map.world_to_pixel(right_corner)
con1 = ConnectionPatch(
    (0, 1), (xpix.value, ypix.value), 'axes fraction', 'data', axesA = ax2, axesB = ax1,
    arrowstyle = '-', color = 'black', lw = 1
)
xpix, ypix = hacore_map.world_to_pixel(
    SkyCoord(right_corner.Tx, left_corner.Ty, frame = hacore_map.coordinate_frame))
con2 = ConnectionPatch(
    (0, 0), (xpix.value, ypix.value), 'axes fraction', 'data', axesA = ax2, axesB = ax1,
    arrowstyle = '-', color = 'black', lw = 1
)
ax2.add_artist(con1)
ax2.add_artist(con2)

# ax2.contour(hmi_small_dim, levels, origin = 'lower', cmap = 'Accent', linewidths = 2)
# ax3.contour(hmi_small_dim, levels, origin = 'lower', cmap = 'Accent', linewidths = 2)
'''




'''core_small2 = hacore_map.submap(left_corner2, top_right = right_corner2)
core_small2 = core_small2.data
ax4 = fig.add_subplot(gs[3:6, 6:9])
im4 = ax4.imshow(core_small2, origin = 'lower', cmap = 'afmhot', vmin = 0, vmax = 4 * hacore_map.mean())
plt.xlabel('Solar X (arcsec)', fontsize = 18)
plt.ylabel('Solar Y (arcsec)', fontsize = 18)
plt.xticks(ticks = tick_pixel2x, labels = tick_arcsec2x, fontsize = 13)
plt.yticks(ticks = tick_pixel2y, labels = tick_arcsec2y, fontsize = 13)'''


'''wing_small2 = hawing_map.submap(left_corner2, top_right = right_corner2)
wing_small2 = wing_small2.data


ax5 = fig.add_subplot(gs[3:6, 10:13])
im5 = ax5.imshow(wing_small2, origin = 'lower', cmap = 'afmhot', vmin = 0, vmax = 4 * hawing_map.mean())
plt.xlabel('Solar X (arcsec)', fontsize = 18)
plt.ylabel('Solar Y (arcsec)', fontsize = 18)
plt.xticks(ticks = tick_pixel2x, labels = tick_arcsec2x, fontsize = 13)
plt.yticks(ticks = tick_pixel2y, labels = tick_arcsec2y, fontsize = 13)


xpix, ypix = hacore_map.world_to_pixel(right_corner2)
con1 = ConnectionPatch(
    (0, 1), (xpix.value, ypix.value), 'axes fraction', 'data', axesA = ax4, axesB = ax1,
    arrowstyle = '-', color = 'black', lw = 1
)
xpix, ypix = hacore_map.world_to_pixel(
    SkyCoord(right_corner2.Tx, left_corner2.Ty, frame = hacore_map.coordinate_frame))
con2 = ConnectionPatch(
    (0, 0), (xpix.value, ypix.value), 'axes fraction', 'data', axesA = ax4, axesB = ax1,
    arrowstyle = '-', color = 'black', lw = 1
)
ax4.add_artist(con1)
ax4.add_artist(con2)

# ax4.contour(hmi_small2_dim, levels, origin = 'lower', cmap = 'Accent', linewidths = 2)
# ax5.contour(hmi_small2_dim, levels, origin = 'lower', cmap = 'Accent', linewidths = 2)
'''



'''position = fig.add_axes([0.08, 0.12, 0.006, 0.76]) # 位置[左, 下, 宽, 高]
cb = fig.colorbar(im4, cax = position)
cb.ax.tick_params(labelsize = 15)
cb.set_label('Count Rate', fontsize = 25)'''

'''position = fig.add_axes([0.92, 0.12, 0.006, 0.76]) # 位置[左, 下, 宽, 高]
cb = fig.colorbar(im5, cax = position)
cb.ax.tick_params(labelsize = 15)
cb.set_label('Count Rate', fontsize = 25)
'''

