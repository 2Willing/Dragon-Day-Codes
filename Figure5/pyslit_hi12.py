#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 13 10:55:26 2022

@author: twinwilling
"""

import matplotlib.pyplot as plt
import numpy as np
import os
import scipy.interpolate as spi
from datetime import datetime as dtdt
from datetime import timedelta as dttd

def interpcurve(N,pX,pY):
    global cumarc,tbins,chordlen
    N=np.transpose(np.linspace(0,1,N))
    nt=N.size
    
    #number of points on the curve
    n=pX.size
    pxy=np.array((pX,pY)).T
    p1=pxy[0,:]
    pend=pxy[-1,:]
    last_segment= np.linalg.norm(np.subtract(p1,pend))
    epsilon= 10*np.finfo(float).eps
    
    #IF the two end points are not close enough lets close the curve
    if last_segment > epsilon*np.linalg.norm(np.amax(abs(pxy),axis=0)):
        pxy=np.vstack((pxy,pend))
        nt = nt + 1
    else:
        print('Contour already closed')
    
    pt=np.zeros((nt,2))
    
    #Compute the chordal arclength of each segment.
    chordlen = (np.sum(np.diff(pxy,axis=0)**2,axis=1))**(1/2)
    #Normalize the arclengths to a unit total
    chordlen = chordlen/np.sum(chordlen)
    #cumulative arclength
    cumarc = np.append(0,np.cumsum(chordlen))
    tbins= np.digitize(N,cumarc) # bin index in which each N is in
    
    #catch any problems at the ends
    tbins[np.where(tbins<=0 | (N<=0))]=1
    tbins[np.where(tbins >= n | (N >= 1))] = n - 1      
    
    s = np.divide((N - cumarc[tbins]),chordlen[tbins-1])
    pt = pxy[tbins,:] + np.multiply((pxy[tbins,:] - pxy[tbins-1,:]),(np.vstack([s]*2)).T)
    return pt 

time0=dtdt(2023,3,19,16,0,0)
l_num=0
u_num=207
ipng_slit=60
times=np.zeros(u_num-l_num).astype(str)
for w in range(u_num-l_num):
    timei=time0+dttd(minutes=(w+l_num)*30)
    times[w]=str(timei)[11:16]
gogo
pngdir='/Volumes/金乌3号/secchi/L0/a/img/composite/'
pngs=os.listdir(pngdir)
pngs=[i for i in pngs if '.png' in i and '._' not in i]
pngs=np.sort(pngs)[l_num:u_num]
png1=plt.imread(pngdir+pngs[ipng_slit])
png1[:,:,:3]=1-(1-png1[:,:,:3])*1.5
N_points=800
plt.figure(figsize=np.array([10,8])/2)
plt.rcParams.update({"font.size":10,'font.family':"Arial",'font.weight':'bold'})
plt.title(None)
plt.imshow(png1)
xu=-6
xl=-6+100
yu=-40
yl=40
plt.xticks(np.arange(0+6*10,int((xl-xu)*10),100),np.arange(xu+6,xl,10).astype(int))
plt.yticks(np.arange(int((yl-yu)*10),0,-100),np.arange(yu,yl,10).astype(int))
plt.xlabel(r'Helioprojective X ($^\circ$)',fontsize=12,fontweight='bold')
plt.ylabel(r'Helioprojective Y ($^\circ$)',fontsize=12,fontweight='bold')
#pos=plt.ginput(2)
#sun-earth
pos=np.array([(6*10, 40*10),
 ((6+80)*10, 40*10)])#plt.ginput(2)
xs=np.linspace(pos[0][0],pos[1][0],N_points+2)
ys=np.linspace(pos[0][1],pos[1][1],N_points+2)
plt.plot(xs,ys,'--',color='cyan')

pos2=np.array([(6*10, 40*10),
 ((6+20)*10, (40+10)*10)])
plt.plot(pos2[:,0],pos2[:,1],':',color='yellow',linewidth=2)

plt.subplots_adjust(left=0.12,bottom=0.12,top=0.98,right=0.99,wspace=0,hspace=0)
plt.savefig(pngdir+'slit_location_eruption.pdf',dpi=200)#,dpi=500)
plt.show()
pool=np.zeros([N_points,len(pngs),3]).astype(int)
base=0
running=0
slit_width=5
if base==0:
    running=running
else:
    print('using base difference instead of running difference')
    running='disabled'
    png2=plt.imread(pngdir+pngs[0])
for i in range(len(pngs)):
    print(pngs[i])
    pici=plt.imread(pngdir+pngs[i])
    if running==1:
        if i==0:
            continue
        png2=plt.imread(pngdir+pngs[i-1])
        pici[:,:,:-1]=pici[:,:,:-1]-png2[:,:,:-1]
    elif base==1:
        if i==0:
            continue
        pici[:,:,:-1]=pici[:,:,:-1]-png2[:,:,:-1]
    for k in range(1,N_points+1):
        x_i=int(xs[k])
        y_i=int(ys[k])
        dyds=(ys[k+1]-ys[k-1])/((xs[k+1]-xs[k-1])**2+(ys[k+1]-ys[k-1])**2)**(1/2)
        dxds=(xs[k+1]-xs[k-1])/((xs[k+1]-xs[k-1])**2+(ys[k+1]-ys[k-1])**2)**(1/2)
        positions_x=(xs[k]-np.arange(-slit_width,slit_width)*dyds).astype(int)
        positions_y=(ys[k]+np.arange(-slit_width,slit_width)*dxds).astype(int)
        positions=np.array([positions_y,positions_x])
        positions=positions.tolist()
        positions=tuple(positions)
        pici_sub=pici[:,:,:-1]
        #for h in range(len(positions_x)):
        #    pool[k-1,i]=pool[k-1,i]+np.sum(pici[positions_x[h],positions_y[h],:-1])
        #pool[k-1,i]=np.sum(pici[y_i-1:y_i+1,x_i-1:x_i+1,:-1])
        pool[k-1,i,:]=(np.sum(pici_sub[positions],axis=0)/len(positions[0])*255).astype(int)
#times=os.listdir('/Volumes/金乌3号/BaiduNetdiskDownload/BBSO数据/20150803/20150803/Hα/')
#times=[i[-10:-8]+':'+i[-8:-6]+':'+i[-6:-4] for i in times if '.fts' in i]#dtdt(2015,8,3,int(i[-10:-8]),int(i[-8:-6]),int(i[-6:-4]))
#times=np.sort(times)[:370]
#timeslices=times
#if (running==0) and (base==0):
#    vmax=np.max(pool)/1.5
#    pool[pool>vmax]=vmax
plt.figure(figsize=(16,8))
plt.rcParams.update({"font.size":20,'font.family':"Arial",'font.weight':'bold'})
plt.imshow(pool[::-1],aspect='auto',cmap='Blues_r')
plt.xticks(np.arange(len(pool[0]))[4:-1:12],times[4:-1:12],fontsize=18)
plt.tick_params(pad=5)
plt.yticks(np.arange(80*10,-1,-100),np.arange(0,80+1,10))
plt.xlabel(r'Time (UT)',fontsize=24,fontweight='bold',labelpad=20)
plt.ylabel(r'Helioprojective X ($^\circ$)',fontsize=24,fontweight='bold')
#plt.yticks([0,len(pool)],[r'P$_{end}$',r'Sun'],fontsize=20)
#plt.ylabel(r'${\rm s}$',fontsize=25)
#plt.xlabel(r'${\rm time}$',fontsize=25)
plt.subplots_adjust(left=0.12, bottom=0.12, right=0.99, top=0.98, wspace=None, hspace=None)
pos=np.array([[24.97844828,759.27688953],
    [44.80603448,697.61046512],
    [60.75071839,635.9440407],
    [78.26508621,549.26199128],
    [89.33548851,496.9037064],
    [100.57112069,435.23728198],
    [114.69827586,353.20930233],
    [129.48635057,277.5806686]])#plt.ginput(8))
print(pos)
func=np.polyfit(pos[:,0],pos[:,1],2)
func=np.poly1d(func)
xs=np.linspace(pos[0,0],pos[-1,0]+60,100)
ys=func(xs)
plt.plot([0,len(pool[0])],np.array([80-91.95,80-92.55])*10,color='green',linewidth=5,label='Earth')
plt.plot(xs,ys,'--',color='cyan',linewidth=6,alpha=0.7,label='CME Front')
plt.grid(True,linestyle=':',color='white',linewidth=1.5)
plt.plot([186,186],[-15*10,len(pool)],linestyle=':',color='limegreen',linewidth=7)
plt.plot([83,83],[len(pool)-35*10,len(pool)],linestyle=':',color='deepskyblue',linewidth=7)#27degrees
plt.xlim([0,len(pool[0])])
plt.ylim([80*10,-15*10])
plt.legend(loc='upper left',fontsize=24)
plt.savefig(pngdir+'slit_map_eruption.pdf')#png',dpi=100)
plt.show()
