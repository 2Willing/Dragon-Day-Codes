import cdflib
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime as dtdt
from datetime import timedelta as dttd


inclin_so=5.2/180.*np.pi
cdf_21 = cdflib.CDF("/Users/twinwilling/Downloads/solar orbiter/solo_L2_mag-rtn-normal-1-minute_20230321_V01.cdf")
cdf_22 = cdflib.CDF("/Users/twinwilling/Downloads/solar orbiter/solo_L2_mag-rtn-normal-1-minute_20230322_V01.cdf")#"solo_L2_mag-rtn-normal_20230320_V01.cdf")
#print(cdf_21.cdf_info())
epoch_time21=cdf_21.varget("EPOCH")
B_RTN21=cdf_21.varget("B_RTN")
epoch_time22=cdf_22.varget("EPOCH")
B_RTN22=cdf_22.varget("B_RTN")
epoch_time=np.append(epoch_time21,epoch_time22)
time0=dtdt(2000,1,1,12,0,0)#J2000.0
time=np.array([time0+dttd(seconds=i/1e9) for i in epoch_time])
times=dtdt(2023,3,21,6,0,0)
timee=dtdt(2023,3,22,6,0,0)
time_fr=dtdt(2023,3,21,9,30,0)
time_frstt=dtdt(2023,3,21,15,40,0)
time_frend=dtdt(2023,3,21,19,50,0)
start=np.where(np.abs(time-times)==np.min(np.abs(time-times)))[0][0]
end=np.where(np.abs(time-timee)==np.min(np.abs(time-timee)))[0][0]
fr=np.where(np.abs(time-time_fr)==np.min(np.abs(time-time_fr)))[0][0]
fr_end=np.where(np.abs(time-time_frend)==np.min(np.abs(time-time_frend)))[0][0]
fr_stt=np.where(np.abs(time-time_frstt)==np.min(np.abs(time-time_frstt)))[0][0]
B_RTN=np.vstack([B_RTN21,B_RTN22])
B_total=(B_RTN[:,0]**2+B_RTN[:,1]**2+B_RTN[:,2]**2)**(1/2)
br=B_RTN[:,0]
bt=B_RTN[:,1]
bn=B_RTN[:,2]
by=bt*np.cos(inclin_so)+bn*np.sin(inclin_so)
bz=bn*np.cos(inclin_so)-bt*np.sin(inclin_so)
br=-br
by=-by

bxy=np.sqrt((-br)**2+(-by)**2)
phi=np.ones(B_total.shape)*999
sign=((by>0)-0.5)*2
phi[(by<999)*(bxy<999)*(br<999)]=sign[(by<999)*(bxy<999)*(br<999)]*np.arccos((br/bxy)[(by<999)*(bxy<999)*(br<999)])/np.pi*180
phi[phi<0]=360+phi[phi<0]
phi=phi-180
theta=np.ones(B_total.shape)*999
theta[(bz<999)*(br<999)*(by<999)]=np.arctan(bz[(bz<999)*(br<999)*(by<999)]/np.sqrt(br[(bz<999)*(br<999)*(by<999)]**2+by[(bz<999)*(br<999)*(by<999)]**2))/np.pi*180.

fig=plt.figure(figsize=(14,8))
plt.rcParams.update({"font.size":20,'font.family':"Arial",'font.weight':'bold'})
ax1=fig.add_subplot(211)
plt.plot([time[fr],time[fr]],[-65,65],color='deepskyblue',linewidth=5)
plt.plot([time[fr_end],time[fr_end]],[-65,65],color='gray',linewidth=5)
plt.plot(time[start:end],br[start:end],color='blue',linewidth=5,label=r'$B_{x}$')
plt.plot(time[start:end],by[start:end],color='green',linewidth=5,label=r'$B_{y}$')
plt.plot(time[start:end],bz[start:end],color='red',linewidth=5,label=r'$B_{z}$')
plt.plot(time[start:end],B_total[start:end],color='black',linewidth=5,label=r'$\boldsymbol{\rm |B|}$')
plt.ylabel('B (nT)',fontweight='bold',fontsize=24)
plt.gca().set_xticklabels(['']*9)
plt.legend(loc='upper right',fontsize=24)
plt.ylim([-65,65])
plt.grid(True)
ax2=fig.add_subplot(212)
plt.plot([time[fr],time[fr]],[-180,180],color='deepskyblue',linewidth=5)
plt.plot([time[fr_end],time[fr_end]],[-180,180],color='gray',linewidth=5)
plt.plot(time[start:end],phi[start:end],color='lime',linewidth=5,label=r'$\phi-180^{\circ}$')
plt.plot(time[start:end],theta[start:end],color='blue',linewidth=5,label=r'$\theta$')
plt.ylim([-180,180])
plt.legend(loc='upper right',fontsize=24)
plt.ylabel(r'$\phi$ & $\theta$ ($^{\circ}$)',fontweight='bold',fontsize=24)
plt.xlabel('Time (UT)',fontweight='bold',fontsize=24)
plt.grid(True)
plt.subplots_adjust(left=0.12, bottom=0.12, right=0.99, top=0.98, wspace=None, hspace=0.03)
plt.savefig('solar_orbiter_mag.pdf')
plt.show()


#MVA
# These functions are written by *** Harry Wheeler ***
# https://github.com/harrywheeleriv/minvarly/tree/master

def magneticVariance(B_vector):  
    '''     
        the matrix M (magnetic variance matrix) is constructed and returned, along with its eigenvalues and eigenvectors. The eigenvalue are normalized by the minimum eigenvalue.
    '''
    M = np.zeros([3,3])
    
    for i in range(3):
        for j in range(3):
            M[i,j] = np.mean(np.multiply(B_vector[i,:],B_vector[j,:])) - np.mean(B_vector[i,:])*np.mean(B_vector[j,:])
        # M[i,0] = np.mean(np.multiply(B_vector[i,:],B_vector[0,:])) - np.mean(B_vector[i,:])*np.mean(B_vector[0,:])
        # M[i,1] = np.mean(np.multiply(B_vector[i,:],B_vector[1,:])) - np.mean(B_vector[i,:])*np.mean(B_vector[1,:])
        # M[i,2] = np.mean(np.multiply(B_vector[i,:],B_vector[2,:])) - np.mean(B_vector[i,:])*np.mean(B_vector[2,:])

    eigen_values, eigen_vectors = np.linalg.eig(M)     
    #sorting the vectors so max = 0, inter = 1, min = 2
    group = []
    for index, row in enumerate(eigen_values):
        group.append([row,eigen_vectors[:,index]])

    group.sort(key = lambda row: row[:][0])

    del eigen_vectors
    for index, row in enumerate(group):
        eigen_values[index] = (group[2-index][0])
        try:
            eigen_vectors = np.concatenate((eigen_vectors, group[2-index][1].reshape(3,1)), axis = 1)
        except:
            eigen_vectors = group[2-index][1].reshape(3,1)


    return M, eigen_values/eigen_values[2], eigen_vectors, group

def minvarRotate(B_vector):   
    ''' 
        This function calls the magneticVariance function which generates the variance matrix then it rotates in inputted data
        into minimum variance coordinates. It also returns the angle between the minimum variance direction and B0.
    '''

    M, eigenValues, eigenVectors, grouping = magneticVariance(B_vector)
    
    b0Direction = np.mean(B_vector,axis=1)/np.linalg.norm(np.mean(B_vector,axis=1))

    print(eigenVectors[:,1])
    angle = 180/np.pi*np.arccos(np.dot(eigenVectors[:,2],b0Direction))
    b = np.zeros([3,np.shape(B_vector[1,:])[0]])
    for i in range(np.shape(B_vector[1,:])[0]):
        b[0,i]=eigenVectors[0,0]*B_vector[0,i] + eigenVectors[1,0]*B_vector[1,i] + eigenVectors[2,0]*B_vector[2,i]
        b[1,i]=eigenVectors[0,1]*B_vector[0,i] + eigenVectors[1,1]*B_vector[1,i] + eigenVectors[2,1]*B_vector[2,i]
        b[2,i]=eigenVectors[0,2]*B_vector[0,i] + eigenVectors[1,2]*B_vector[1,i] + eigenVectors[2,2]*B_vector[2,i]
    
    RMS = np.sqrt( sum( np.power(b[2,:],2) ) / len(b[2,:]) );

    return b, angle, RMS, eigenValues, eigenVectors


bx=br[fr_stt:fr_end]#2870]#3380]
by=by[fr_stt:fr_end]#2870]#3380]
bz=bz[fr_stt:fr_end]#2870]#3380]
bx=bx[bx<999]
by=by[by<999]
bz=bz[bz<999]
#moving average
order=31
bx = np.convolve(bx, np.ones((order,))/order, mode='valid')
by = np.convolve(by, np.ones((order,))/order, mode='valid')
bz = np.convolve(bz, np.ones((order,))/order, mode='valid')
plt.plot(bx,color='blue')
plt.plot(by,color='green')
plt.plot(bz,color='red')
plt.show()
plt.close()
bxyz=np.vstack((bx,by,bz))
b,angle,RMS,eigenValues,eigenVectors=minvarRotate(bxyz)
exmax,eymax,ezmax=eigenVectors[:,1]

theta=np.arctan(ezmax/eymax)/np.pi*180.
print('flux rope tilt angle:',theta)



cdf = cdflib.CDF("solo_L2_swa-pas-grnd-mom_20230321_V02.cdf")
cdf_2=cdflib.CDF("solo_L2_swa-pas-grnd-mom_20230322_V02.cdf")
#print(cdf.cdf_info())
epoch_time=cdf.varget("Epoch")
V_RTN=cdf.varget("V_RTN")
epoch_time2=cdf_2.varget("Epoch")
V_RTN2=cdf_2.varget("V_RTN")
epoch_time=np.append(epoch_time,epoch_time2)
V_RTN=np.vstack([V_RTN,V_RTN2])
T=cdf.varget("T")
T2=cdf_2.varget("T")
T=np.append(T,T2)
time0=dtdt(2000,1,1,12,0,0)#J2000.0
time=np.array([time0+dttd(seconds=i/1e9) for i in epoch_time])
start=np.where(np.abs(time-times)==np.min(np.abs(time-times)))[0][0]
end=np.where(np.abs(time-timee)==np.min(np.abs(time-timee)))[0][0]
fr=np.where(np.abs(time-time_fr)==np.min(np.abs(time-time_fr)))[0][0]
fr_end=np.where(np.abs(time-time_frend)==np.min(np.abs(time-time_frend)))[0][0]


fig=plt.figure(figsize=(14,8))
plt.rcParams.update({"font.size":20,'font.family':"Arial",'font.weight':'bold'})
ax1=fig.add_subplot(211)
plt.plot([time[fr],time[fr]],[450,1000],color='deepskyblue',linewidth=5)
plt.plot([time[fr_end],time[fr_end]],[450,1000],color='gray',linewidth=5)
plt.plot(time[start:end],V_RTN[:,0][start:end],color='green',linewidth=5)
#plt.plot(time[start:end],V_RTN[:,1][start:end],color='green')
#plt.plot(time[start:end],V_RTN[:,2][start:end],color='red')
plt.ylabel(r'$V_{sw}$ (km s$^{-1}$)',fontweight='bold',fontsize=24)
plt.gca().set_xticklabels(['']*9)
plt.legend(loc='upper right',fontsize=24)
plt.ylim([450,1000])
plt.grid(True)
ax2=fig.add_subplot(212)
plt.plot([time[fr],time[fr]],[0,1],color='deepskyblue',linewidth=5)
plt.plot([time[fr_end],time[fr_end]],[0,1],color='gray',linewidth=5)
plt.plot(time[start:end],T[start:end]*11605/1e6,color='red',linewidth=5)
plt.ylim([0,1])
plt.legend(loc='upper right',fontsize=24)
plt.ylabel(r'$T_{p}$ (MK)',fontweight='bold',fontsize=24)
plt.xlabel('Time (UT)',fontweight='bold',fontsize=24)
plt.grid(True)
plt.subplots_adjust(left=0.12, bottom=0.12, right=0.99, top=0.98, wspace=None, hspace=0.03)
plt.savefig('solar_orbiter_swa.pdf')
plt.show()


