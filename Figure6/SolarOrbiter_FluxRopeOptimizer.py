#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  7 19:07:54 2023

@author: twinwilling
"""

import numpy as np
from scipy import interpolate
from scipy.optimize import minimize
import matplotlib.pyplot as plt
import spacepy.coordinates as spc
from spacepy.time import Ticktock
from datetime import datetime as dtdt
from datetime import timedelta as dttd
import cdflib

vsw_so=550#km/s
r_sun=696300#km
dz=0.002*r_sun
t_cad=60#s
longitude_fr=0
latitude_fr=0
Au2Rs=215.035
N_low=2
#N_high=170
dlon_fr=10
dlon_gcs=36
inclin_so=5.2*np.pi/180.#angle between solar orbiter RTN and ecliptic cs

nxyz=np.loadtxt('nxyz.dat')
nx=int(nxyz[0])
ny=int(nxyz[1])
nz=int(nxyz[2])
#datacube=np.loadtxt('flux_rope.dat')
br=np.loadtxt('bz.dat')#datacube[:,:,:,2]#bz
br=br.reshape([nz,ny,nx])
br=br.transpose((2,1,0))[:,:,N_low:]
bt=np.loadtxt('bx.dat')#datacube[:,:,:,0]#bx
bt=bt.reshape([nz,ny,nx])
bt=bt.transpose((2,1,0))[:,:,N_low:]
bn=np.loadtxt('by.dat')#datacube[:,:,:,1]#by
bn=bn.reshape([nz,ny,nx])
bn=bn.transpose((2,1,0))[:,:,N_low:]
longitudes=np.loadtxt('lon_rope.dat')
latitudes=np.loadtxt('lat_rope.dat')
longitudes_index=np.arange(len(longitudes))
latitudes_index=np.arange(len(latitudes))
so_br=np.loadtxt('br.dat')
so_bt=np.loadtxt('bt.dat')
so_bn=np.loadtxt('bn.dat')
so_bt_y=so_bt*np.cos(inclin_so)+so_bn*np.sin(inclin_so)
so_bn_z=so_bn*np.cos(inclin_so)-so_bt*np.sin(inclin_so)
so_bt=so_bt_y
so_bn=so_bn_z
solar_orbiter_locations=np.loadtxt('Loc_SO.dat')
solar_orbiter_lon=solar_orbiter_locations[0]
solar_orbiter_lat=solar_orbiter_locations[1]
solar_orbiter_dis=solar_orbiter_locations[2]*Au2Rs
expansion_so=solar_orbiter_dis**2#(2*0.78)
nz=nz-N_low

t_range=len(so_br)*60
dt=dz*np.sqrt(expansion_so)/vsw_so

#deflections are longitudinal
def cut(deflection_lon,deflection_lat):
    global where_lon,where_lat
    lon=deflection_lon#solar_orbiter_lon-deflection_lon
    lat=deflection_lat#solar_orbiter_lat-deflection_lat
    where_lon=np.where(np.abs(longitudes_index-lon)==np.min(np.abs(longitudes_index-lon)))[0][0]
    where_lat=np.where(np.abs(latitudes_index-lat)==np.min(np.abs(latitudes_index-lat)))[0][0]
    br_cut=br[where_lon,where_lat,::-1]
    bt_cut=bt[where_lon,where_lat,::-1]
    bn_cut=bn[where_lon,where_lat,::-1]
    return br_cut,bt_cut,bn_cut,dt

def b_smooth(deflection_lon,deflection_lat):
    dlon_floor=int(np.floor(deflection_lon))
    dlon_ceil=int(np.ceil(deflection_lon))
    dlat_floor=int(np.floor(deflection_lat))
    dlat_ceil=int(np.ceil(deflection_lat))
    dlon_res=deflection_lon-dlon_floor
    dlat_res=deflection_lat-dlat_floor
    lonclatc_br,lonclatc_bt,lonclatc_bn,dt=cut(dlon_ceil,dlat_ceil)
    lonflatc_br,lonflatc_bt,lonflatc_bn,dt=cut(dlon_floor,dlat_ceil)
    lonclatf_br,lonclatf_bt,lonclatf_bn,dt=cut(dlon_ceil,dlat_floor)
    lonflatf_br,lonflatf_bt,lonflatf_bn,dt=cut(dlon_floor,dlat_floor)
    br_cut=(lonclatc_br*dlon_res+lonflatc_br*(1-dlon_res))*dlat_res+(lonclatf_br*dlon_res+lonflatf_br*(1-dlon_res))*(1-dlat_res)
    bt_cut=(lonclatc_bt*dlon_res+lonflatc_bt*(1-dlon_res))*dlat_res+(lonclatf_bt*dlon_res+lonflatf_bt*(1-dlon_res))*(1-dlat_res)
    bn_cut=(lonclatc_bn*dlon_res+lonflatc_bn*(1-dlon_res))*dlat_res+(lonclatf_bn*dlon_res+lonflatf_bn*(1-dlon_res))*(1-dlat_res)
    return br_cut,bt_cut,bn_cut,dt

def fill_array(br_cut,bt_cut,bn_cut):#,N_points):#,dt,time_amp):
    N_points=int(len(br_cut)*dt/t_cad)
    f=interpolate.interp1d(np.arange(len(br_cut)),br_cut,kind='cubic')
    br_exp=f(np.linspace(0,len(br_cut)-1,N_points))
    f=interpolate.interp1d(np.arange(len(bt_cut)),bt_cut,kind='cubic')
    bt_exp=f(np.linspace(0,len(bt_cut)-1,N_points))
    f=interpolate.interp1d(np.arange(len(bn_cut)),bn_cut,kind='cubic')
    bn_exp=f(np.linspace(0,len(bn_cut)-1,N_points))
    return br_exp,bt_exp,bn_exp


t_amp=1

def loss(x):
    global br_exp,bt_exp,bn_exp,so_br_i,so_bn_i,so_bt_i,so_br_f,so_br_c,mag_amp#,time_amp
    deflection_lon,deflection_lat,t_start_so,mag_amp=x
    #deflection_lon,t_start_so,mag_amp=x
    #deflection_lat=solar_orbiter_lat_index
    t_start_floor=int(np.floor(t_start_so))
    t_start_res=t_start_so-t_start_floor
    br_cut,bt_cut,bn_cut,dt=b_smooth(deflection_lon,deflection_lat)
    #N_points=int(np.floor(len(br_cut)*dt*time_amp/t_cad))
    #N_res=len(br_cut)*dt*time_amp/t_cad-np.floor(len(br_cut)*dt*time_amp/t_cad)
    br_exp,bt_exp,bn_exp=fill_array(br_cut,bt_cut,bn_cut)#,N_points)
    #br_exp_1,bt_exp_1,bn_exp_1=fill_array(br_cut,bt_cut,bn_cut,N_points+1)
    
    br_exp=br_exp/expansion_so*mag_amp*1e5#Gauss2nT
    bt_exp=bt_exp/expansion_so*mag_amp*1e5
    bn_exp=bn_exp/expansion_so*mag_amp*1e5
    if t_start_floor+len(br_exp)>len(so_br):
        so_br_f=np.append(so_br[int(t_start_floor):],np.zeros(int(t_start_floor+len(br_exp)-len(so_br))))
        so_bt_f=np.append(so_bt[int(t_start_floor):],np.zeros(int(t_start_floor+len(br_exp)-len(so_br))))
        so_bn_f=np.append(so_bn[int(t_start_floor):],np.zeros(int(t_start_floor+len(br_exp)-len(so_br))))
        so_br_c=np.append(so_br[int(t_start_floor+1):],np.zeros(int(t_start_floor+len(br_exp)-len(so_br)+1)))
        so_bt_c=np.append(so_bt[int(t_start_floor+1):],np.zeros(int(t_start_floor+len(br_exp)-len(so_br)+1)))
        so_bn_c=np.append(so_bn[int(t_start_floor+1):],np.zeros(int(t_start_floor+len(br_exp)-len(so_br)+1)))
    else:
        so_br_f=so_br[int(t_start_floor):int(t_start_floor+len(br_exp))]
        so_bt_f=so_bt[int(t_start_floor):int(t_start_floor+len(br_exp))]
        so_bn_f=so_bn[int(t_start_floor):int(t_start_floor+len(br_exp))]
        so_br_c=so_br[int(t_start_floor+1):int(t_start_floor+1+len(br_exp))]
        so_bt_c=so_bt[int(t_start_floor+1):int(t_start_floor+1+len(br_exp))]
        so_bn_c=so_bn[int(t_start_floor+1):int(t_start_floor+1+len(br_exp))]
    so_br_i=so_br_f*(1-t_start_res)+so_br_c*t_start_res
    so_bt_i=so_bt_f*(1-t_start_res)+so_bt_c*t_start_res
    so_bn_i=so_bn_f*(1-t_start_res)+so_bn_c*t_start_res
    loss_br=np.mean((br_exp-so_br_i)**2)
    loss_bt=np.mean((bt_exp-so_bt_i)**2)
    loss_bn=np.mean((bn_exp-so_bn_i)**2)

    '''br_exp_1=br_exp_1/expansion_so*mag_amp*1e5#Gauss2nT
    bt_exp_1=bt_exp_1/expansion_so*mag_amp*1e5
    bn_exp_1=bn_exp_1/expansion_so*mag_amp*1e5
    if t_start_floor+len(br_exp_1)>len(so_br):
        so_br_f=np.append(so_br[int(t_start_floor):],np.zeros(int(t_start_floor+len(br_exp_1)-len(so_br))))
        so_bt_f=np.append(so_bt[int(t_start_floor):],np.zeros(int(t_start_floor+len(br_exp_1)-len(so_br))))
        so_bn_f=np.append(so_bn[int(t_start_floor):],np.zeros(int(t_start_floor+len(br_exp_1)-len(so_br))))
        so_br_c=np.append(so_br[int(t_start_floor+1):],np.zeros(int(t_start_floor+len(br_exp_1)-len(so_br)+1)))
        so_bt_c=np.append(so_bt[int(t_start_floor+1):],np.zeros(int(t_start_floor+len(br_exp_1)-len(so_br)+1)))
        so_bn_c=np.append(so_bn[int(t_start_floor+1):],np.zeros(int(t_start_floor+len(br_exp_1)-len(so_br)+1)))
    else:
        so_br_f=so_br[int(t_start_floor):int(t_start_floor+len(br_exp_1))]
        so_bt_f=so_bt[int(t_start_floor):int(t_start_floor+len(br_exp_1))]
        so_bn_f=so_bn[int(t_start_floor):int(t_start_floor+len(br_exp_1))]
        so_br_c=so_br[int(t_start_floor+1):int(t_start_floor+1+len(br_exp_1))]
        so_bt_c=so_bt[int(t_start_floor+1):int(t_start_floor+1+len(br_exp_1))]
        so_bn_c=so_bn[int(t_start_floor+1):int(t_start_floor+1+len(br_exp_1))]
    so_br_i=so_br_f*(1-t_start_res)+so_br_c*t_start_res
    so_bt_i=so_bt_f*(1-t_start_res)+so_bt_c*t_start_res
    so_bn_i=so_bn_f*(1-t_start_res)+so_bn_c*t_start_res
    loss_br_1=np.mean((br_exp_1-so_br_i)**2)
    loss_bt_1=np.mean((bt_exp_1-so_bt_i)**2)
    loss_bn_1=np.mean((bn_exp_1-so_bn_i)**2)'''

    loss_tot=loss_br+loss_bt+loss_bn#loss_br*(1-N_res)+loss_br_1*N_res+loss_bt*(1-N_res)+loss_bt_1*N_res+loss_bn*(1-N_res)+loss_bn_1*N_res
    return loss_tot

'''solar_orbiter_lon_index=200#np.where(np.abs(solar_orbiter_lon-longitudes)==np.min(np.abs(solar_orbiter_lon-longitudes)))[0][0]
solar_orbiter_lat_index=250#np.where(np.abs(solar_orbiter_lat-latitudes)==np.min(np.abs(solar_orbiter_lat-latitudes)))[0][0]
x_0=np.array([solar_orbiter_lon_index,solar_orbiter_lat_index,500.,2.])#,1.])
results=minimize(loss,x_0,bounds=[(np.min(longitudes_index),np.max(longitudes_index)),(np.min(latitudes_index),np.max(latitudes_index)),(0,800),(1,4)])#,(1.0,2.0)])'''

solar_orbiter_lon_index=np.where(np.abs(solar_orbiter_lon-longitudes)==np.min(np.abs(solar_orbiter_lon-longitudes)))[0][0]
solar_orbiter_lat_index=np.where(np.abs(latitudes-solar_orbiter_lat)==np.min(np.abs(latitudes-solar_orbiter_lat)))[0][0]
results_all=[]
loss_all=[]
its=np.arange(200.,400.,50.)
#its=np.arange(400.,800.,50.)
iams=np.arange(1.,4.,0.25)
isols=np.arange(50,200,50)
istarts=[]
iamps=[]
isolons=[]
for itstart in its:
    for iamp in iams:
        for solar_orbiter_lon_index in isols:
            x_0=np.array([solar_orbiter_lon_index,solar_orbiter_lat_index,itstart,iamp])
            results=minimize(loss,x_0,bounds=[(np.min(longitudes_index),np.max(longitudes_index)),(np.min(latitudes_index),np.max(latitudes_index)),(200,800),(1,4)])#,(1.0,2.0)])
            results_all.append(results)
            loss_all.append(results.fun)
            print(itstart,results.fun)
            istarts.append(itstart)
            iamps.append(iamp)
            isolons.append(solar_orbiter_lon_index)

iii=np.where(loss_all==np.min(loss_all))[0][0]
#tstart0=its[int(np.floor(iii/len(iams)/len(isols)))]
#so_lon0=isols[int(np.floor((iii-int(np.floor(iii/len(iams)/len(isols)))*len(iams)*len(isols))/len(isols)))]
#amp0=iams[int(iii-int(np.floor(iii/len(iams)/len(isols)))*len(iams)*len(isols)-int(np.floor((iii-its[int(np.floor(iii/len(iams)/len(isols)))]*len(iams)*len(isols))/len(isols)))*len(isols))]
tstart0=istarts[iii]
so_lon0=isolons[iii]
amp0=iamps[iii]
print(so_lon0,tstart0,amp0)

loss_all=np.array(loss_all)
results=results_all[iii]

###shiyan
x_0=np.array([188.,solar_orbiter_lat_index,328.,1.79])#,1.])
results=minimize(loss,x_0,bounds=[(np.min(longitudes_index),np.max(longitudes_index)),(np.min(latitudes_index),np.max(latitudes_index)),(0,800),(1,4)])


print('lon_index:',longitudes[int(results.x[0])],'lat_index:',latitudes[int(results.x[1])])#latitudes[int(results.x[1])])
paras=results.x#here, the best-fit results.x=(188.53,408.91,328.5,1.78)
fun=loss(paras)#results.x)


cdf_21 = cdflib.CDF("solo_L2_mag-rtn-normal-1-minute_20230321_V01.cdf")
cdf_22 = cdflib.CDF("solo_L2_mag-rtn-normal-1-minute_20230322_V01.cdf")
epoch_time21=cdf_21.varget("EPOCH")
epoch_time22=cdf_22.varget("EPOCH")
epoch_time=np.append(epoch_time21,epoch_time22)
time0=dtdt(2000,1,1,12,0,0)#J2000.0
time_so=np.array([time0+dttd(seconds=i/1e9) for i in epoch_time])
time_so=time_so[int(paras[2]):int(paras[2])+len(so_bn_i)]
fig=plt.figure(figsize=(8,5.5))
plt.rcParams.update({"font.size":15,'font.family':"Arial",'font.weight':'bold'})
f1,=plt.plot(time_so,-br_exp,color='blue',linewidth=3)
m1,=plt.plot(time_so,-so_br_i,'b:',label=r'$B_{x}$_SO',linewidth=3)
f2,=plt.plot(time_so,-bt_exp,color='green',linewidth=3)
m2,=plt.plot(time_so,-so_bt_i,'g:',label=r'$B_{y}$_SO',linewidth=3)
f3,=plt.plot(time_so,bn_exp,color='red',linewidth=3)
m3,=plt.plot(time_so,so_bn_i,'r:',label=r'$B_{z}$_SO',linewidth=3)
plt.grid(True)
plt.ylabel('B (nT)',fontweight='bold',fontsize=15)
plt.xlabel('Time (UT)',fontweight='bold',fontsize=15)
plt.gca().set_xticklabels(['06:00','08:00','10:00','12:00','14:00','16:00','18:00','20:00'])
a=plt.legend([f1,f2,f3],[r'$B_{x}$_f',r'$B_{y}$_f',r'$B_{z}$_f'],bbox_to_anchor=(0.21,1),fontsize=15)
plt.legend([m1,m2,m3],[r'$B_{x}$_m',r'$B_{y}$_m',r'$B_{z}$_m'],loc="upper right",fontsize=15)
plt.gca().add_artist(a)
plt.legend(loc='upper right',)
#plt.xlim([-10,len(bx_gsm)+10])
plt.ylim([-60,60])
plt.subplots_adjust(left=0.12, bottom=0.14, right=0.99, top=0.98, wspace=None, hspace=0.03)
plt.savefig('Solar_Orbiter_fitting.pdf')
plt.show()
plt.close()

#plt.plot(bn_exp,color='red')
#plt.plot(br_exp,color='blue')
#plt.plot(bt_exp,color='green')
#plt.plot(so_bn_i,color='red')
#plt.plot(so_br_i,color='blue')
#plt.plot(so_bt_i,color='green')
#plt.grid(True)
#plt.show()
#print(results)

t_start_e=int(1900-260)
file=np.loadtxt('OMNI_HRO2_1MIN_2074319.csv',dtype=str)
data=np.array([i.split(',') for i in file])
print(data[0])
data=data[1:]
times_fact=data[:,0][:7200][t_start_e:]
times_fact=np.array([dtdt(int(i[:4]),int(i[5:7]),int(i[8:10]),int(i[11:13]),int(i[14:16]),int(i[17:19])) for i in times_fact])
bx_gsm_fact=data[:,2].astype(float)[:7200][t_start_e:]
by_gsm_fact=data[:,3+2].astype(float)[:7200][t_start_e:]
bz_gsm_fact=data[:,4+2].astype(float)[:7200][t_start_e:]

expansion_e=Au2Rs**2
amp_e=paras[3]*Au2Rs/solar_orbiter_dis
dt_e=dz*np.sqrt(expansion_e)/vsw_so
e_angle=-23.5*np.pi/180
loc_E=np.loadtxt('Loc_E.dat')
lon_E=loc_E[0]
lat_E=loc_E[1]
so_lon_exp=longitudes[int(paras[0])]
so_lat_exp=latitudes[int(paras[1])]
lon_E_df=so_lon_exp+(lon_E-solar_orbiter_lon)*dlon_fr/dlon_gcs
lat_E_df=so_lat_exp+(lat_E-solar_orbiter_lat)
where_e_lon=np.where(np.abs(longitudes-lon_E_df)==np.min(np.abs(longitudes-lon_E_df)))[0][0]
where_e_lat=np.where(np.abs(latitudes-lat_E_df)==np.min(np.abs(latitudes-lat_E_df)))[0][0]
br_exp_e=br[where_e_lon,where_e_lat,::-1]/expansion_e*amp_e*1e5
bt_exp_e=bt[where_e_lon,where_e_lat,::-1]/expansion_e*amp_e*1e5
bn_exp_e=bn[where_e_lon,where_e_lat,::-1]/expansion_e*amp_e*1e5
bx_gse=-br_exp_e[:185]
by_gse=-bt_exp_e[:185]#*np.cos(e_angle)+bn_exp_e*np.sin(e_angle)
bz_gse=bn_exp_e[:185]#*np.cos(e_angle)+bt_exp_e*np.sin(e_angle)
times=Ticktock(['2023-03-23T16:00:00']*len(bz_gse),'ISO','UTC')
b_gse=np.vstack([bx_gse, by_gse, bz_gse]).T
gse_coords = spc.Coords(b_gse, 'GSE', 'car', units=['nT', 'nT', 'nT'],ticks=times)
gsm_coords = gse_coords.convert('GSM', 'car')
bx_gsm = gsm_coords.data[:,0]
by_gsm = gsm_coords.data[:,1]
bz_gsm = gsm_coords.data[:,2]
bxs=[]
bys=[]
bzs=[]
for i in range(len(bz_gsm)):
    time_i=int(i*dt_e/60.)
    if bx_gsm_fact[time_i]>900:
        bxs.append(np.nan)
        bys.append(np.nan)
        bzs.append(np.nan)
    else:
        bxs.append(bx_gsm_fact[time_i])
        bys.append(by_gsm_fact[time_i])
        bzs.append(bz_gsm_fact[time_i])
plt.plot(bx_gsm,color='blue')
plt.plot(bxs,'b:')
plt.plot(by_gsm,color='green')
plt.plot(bys,'g:')
plt.plot(bz_gsm,color='red')
plt.plot(bzs,'r:')
plt.grid(True)
plt.show()

gamma=0.75
tau=7.7*60.*60.

#time=data[:,0][:7200][t_start_e:]
#n=data[:,6+2].astype(float)[:7200][t_start_e:]
#pressure_kinetic=data[:,8+2].astype(float)[:7200][t_start_e:]#n*vsw_so**2*1e-2
sym_h_fact=data[:,12+2].astype(float)[:7200][t_start_e:]
bz_gsm_fact=data[:,4+2].astype(float)[:7200][t_start_e:]
#bz_gsm=bz_gsm_fact
#dt_e=60
Ey=-(vsw_so*1e3)*(bz_gsm*1e-9)*1e3#mV
sym_h=np.zeros(len(bz_gsm))
sym_h_stars=np.zeros(len(bz_gsm))
n_so=60
n_e=n_so/(Au2Rs/solar_orbiter_dis)**(2)
pressure_kinetic=((1.67e-27*n_e)*(vsw_so*1e3)**2)/1.6e-19#eV/cm^3
sym_h_star_0=0-0.2*gamma*pressure_kinetic**(1/2)+20*gamma
sym_h_stars[0]=sym_h_star_0
sym_h_fact_array=np.zeros(len(bz_gsm))
bz_gsm_fact_array=np.zeros(len(bz_gsm))
#bz_gsm_last=bz_gsm[0]
time_is=[times_fact[0]]
for i in range(1,len(bz_gsm)):
    time_i=int(i*dt_e/60.)
    #pressure_i=pressure_kinetic[time_i]
    #if bz_gsm[i]>900:
    #    bz_gsm[i]=bz_gsm_last
    #else:
    #    bz_gsm_last=bz_gsm[i]
    Ey_i=Ey[i]
    if Ey_i>=0.5:
        Qe=-1.5e-3*gamma*(Ey_i-0.5)
        sym_h_stars[i]=sym_h_stars[i-1]+(Qe-sym_h_stars[i-1]/tau)*dt_e
    else:
        Qe=0
        sym_h_stars[i]=sym_h_stars[i-1]+(Qe-sym_h_stars[i-1]/tau)*dt_e
    sym_h[i]=sym_h_stars[i]+0.2*gamma*pressure_kinetic**(1/2)-20*gamma#pressure_i**(1/2)-20
    sym_h_fact_array[i]=sym_h_fact[time_i]
    if bz_gsm_fact[time_i]<999:
        bz_gsm_fact_array[i]=bz_gsm_fact[time_i]
    else:
        bz_gsm_fact_array[i]=bz_gsm_fact_array[i-1]
    time_is.append(times_fact[time_i])
plt.plot(sym_h)
plt.plot(sym_h_fact_array)
#plt.plot(bz_gsm_fact_array)
#plt.plot(bz_gsm)
plt.show()
plt.close()




fig=plt.figure(figsize=(8,8))
plt.rcParams.update({"font.size":15,'font.family':"Arial",'font.weight':'bold'})
ax1=fig.add_subplot(211)
p1,=plt.plot(time_is,bx_gsm,color='blue',linewidth=3)
m1,=plt.plot(time_is,bxs,'b:',label=r'$B_{x}$_m',linewidth=3)
p2,=plt.plot(time_is,by_gsm,color='green',linewidth=3)
m2,=plt.plot(time_is,bys,'g:',label=r'$B_{y}$_m',linewidth=3)
p3,=plt.plot(time_is,bz_gsm,color='red',linewidth=3)
m3,=plt.plot(time_is,bzs,'r:',label=r'$B_{z}$_m',linewidth=3)
plt.grid(True)
plt.ylabel('B (nT)',fontweight='bold',fontsize=15)
plt.gca().set_xticklabels(['']*9)
a=plt.legend([p1,p2,p3],[r'$B_{x}$_p',r'$B_{y}$_p',r'$B_{z}$_p'],bbox_to_anchor=(0.21,1),fontsize=15)
plt.legend([m1,m2,m3],[r'$B_{x}$_m',r'$B_{y}$_m',r'$B_{z}$_m'],loc="upper right",fontsize=15)
plt.gca().add_artist(a)
plt.legend(loc='upper right',)
#plt.xlim([-10,len(bx_gsm)+10])
plt.ylim([-25,25])
plt.grid(True)
ax2=fig.add_subplot(212)
plt.plot(time_is,sym_h,linewidth=3,label='SYM-H_p')
plt.plot(time_is,sym_h_fact_array,":",linewidth=3,label='SYM-H_m')
plt.ylim([-200,50])
#plt.xlim([-10,len(bx_gsm)+10])
plt.legend(loc='upper right',fontsize=15)
plt.ylabel(r'SYM-H (nT)',fontweight='bold',fontsize=15)
plt.xlabel('Time (UT)',fontweight='bold',fontsize=15,labelpad=15)
plt.gca().set_xticklabels(['03:00','06:00','09:00','12:00','15:00','18:00','21:00','00:00','03:00','06:00'])
plt.grid(True)
plt.subplots_adjust(left=0.12, bottom=0.12, right=0.99, top=0.98, wspace=None, hspace=0.03)
plt.savefig('storm_prediction.pdf')
plt.show()
