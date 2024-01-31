import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime as dtdt
from datetime import timedelta as dttd
#file=np.loadtxt('OMNI_HRO_1MIN_124141.csv',dtype=str)
file=np.loadtxt('OMNI_HRO2_1MIN_2074319.csv',dtype=str)
data=np.array([i.split(',') for i in file])
print(data[0])
data=data[1:]
time=data[:,0][:7200]
b=data[:,1].astype(float)[:7200]
bx_gse=data[:,2].astype(float)[:7200]

by_gse=data[:,3].astype(float)[:7200]
bz_gse=data[:,4].astype(float)[:7200]

by_gsm=data[:,3+2].astype(float)[:7200]
bz_gsm=data[:,4+2].astype(float)[:7200]
v_sw=data[:,5+2].astype(float)[:7200]
density_p=data[:,6+2].astype(float)[:7200]
tem=data[:,7+2].astype(float)[:7200]
pressure_kinetic=data[:,8+2].astype(float)[:7200]
beta=data[:,9+2].astype(float)[:7200]
m_a=data[:,10+2].astype(float)[:7200]
sym_d=data[:,11+2].astype(float)[:7200]
sym_h=data[:,12+2].astype(float)[:7200]
asy_d=data[:,13+2].astype(float)[:7200]
asy_h=data[:,14+2].astype(float)[:7200]
bxy=np.sqrt(bx_gse**2+by_gsm**2)
phi=np.ones(b.shape)*999
#sign=((bx_gse<0)-0.5)*2
#phi[(by_gsm<999)*(bxy<999)*(bx_gse<999)]=360*(-sign[(by_gsm<999)*(bxy<999)*(bx_gse<999)]+1)/2.+sign[(by_gsm<999)*(bxy<999)*(bx_gse<999)]*np.arccos((by_gsm/bxy)[(by_gsm<999)*(bxy<999)*(bx_gse<999)])/np.pi*180
#phi[phi>180]=phi[phi>180]-360
sign=((by_gsm>0)-0.5)*2
phi[(by_gsm<999)*(bxy<999)*(bx_gse<999)]=sign[(by_gsm<999)*(bxy<999)*(bx_gse<999)]*np.arccos((bx_gse/bxy)[(by_gsm<999)*(bxy<999)*(bx_gse<999)])/np.pi*180
phi[phi<0]=360+phi[phi<0]
phi=phi-180
theta=np.ones(b.shape)*999
theta[(bz_gsm<999)*(bx_gse<999)*(by_gsm<999)]=np.arctan(bz_gsm[(bz_gsm<999)*(bx_gse<999)*(by_gsm<999)]/np.sqrt(bx_gse[(bz_gsm<999)*(bx_gse<999)*(by_gsm<999)]**2+by_gsm[(bz_gsm<999)*(bx_gse<999)*(by_gsm<999)]**2))/np.pi*180.#360
labels=['']*10

time=np.array([dtdt(int(i[:4]),int(i[5:7]),int(i[8:10]),int(i[11:13]),int(i[14:16]),int(i[17:19])) for i in time])
fig=plt.figure(figsize=(10,10))
plt.title('Time profiles of solar wind properties',fontsize=15)
plt.rcParams['axes.axisbelow'] = True
plt.axis('off')
ax1=plt.subplot(10,1,1)
plt.plot([dtdt(2023,3,23,7,43,21),dtdt(2023,3,23,7,43,21)],[-1000,1000],color='blue')
plt.plot([dtdt(2023,3,23,14,20,00),dtdt(2023,3,23,14,20,00)],[-1000,1000],color='orange')
plt.plot([dtdt(2023,3,23,17,40,00),dtdt(2023,3,23,17,40,00)],[-1000,1000],color='red')
plt.plot([dtdt(2023,3,24,8,20,00),dtdt(2023,3,24,8,20,00)],[-1000,1000],color='gray')
ax1.plot(time[b<999],b[b<999],linewidth=0.5,color='blue')
ax1.set_ylabel(r'$|B|$ (nT)',labelpad=10,fontsize=7)
ax1.set_xticklabels(labels)
ax1.set_ylim([0,25])
ax1.set_yticks([0,5,10,15,20,25])
plt.yticks(size = 6)
plt.grid(True, axis='y')
#ax2=plt.subplot(10,1,2)
#ax2.plot(time[by_gsm<999],by_gsm[by_gsm<999],linewidth=0.5,color='green')
#ax2.set_ylabel(r'$B_{y}$',labelpad=14,fontsize=7)
#ax2.set_xticklabels(labels)
#ax2.set_ylim([-25,25])
#plt.grid(True, axis='y')
ax3=plt.subplot(10,1,2)
plt.plot([dtdt(2023,3,23,7,43,21),dtdt(2023,3,23,7,43,21)],[-1000,1000],color='blue')
plt.plot([dtdt(2023,3,23,14,20,00),dtdt(2023,3,23,14,20,00)],[-1000,1000],color='orange')
plt.plot([dtdt(2023,3,23,17,40,00),dtdt(2023,3,23,17,40,00)],[-1000,1000],color='red')
plt.plot([dtdt(2023,3,24,8,20,00),dtdt(2023,3,24,8,20,00)],[-1000,1000],color='gray')
ax3.plot(time[bz_gsm<999],bz_gsm[bz_gsm<999],linewidth=0.5,color='red')
plt.fill_between(time[bz_gsm<999],bz_gsm[bz_gsm<999],0,where=(bz_gsm[bz_gsm<999]<0),color='purple',alpha=0.7)
ax3.set_ylabel(r'$B_{z}$ (nT)',labelpad=5,fontsize=7)
ax3.set_xticklabels(labels)
ax3.set_ylim([-22,12])
ax3.set_yticks([-20,-15,-10,-5,0,5,10])
plt.yticks(size = 6)
plt.grid(True, axis='y')

ax35=plt.subplot(10,1,3)
plt.plot([dtdt(2023,3,23,7,43,21),dtdt(2023,3,23,7,43,21)],[-1000,1000],color='blue')
plt.plot([dtdt(2023,3,23,14,20,00),dtdt(2023,3,23,14,20,00)],[-1000,1000],color='orange')
plt.plot([dtdt(2023,3,23,17,40,00),dtdt(2023,3,23,17,40,00)],[-1000,1000],color='red')
plt.plot([dtdt(2023,3,24,8,20,00),dtdt(2023,3,24,8,20,00)],[-1000,1000],color='gray')
ax35.plot(time[bx_gse<999],bx_gse[bx_gse<999],linewidth=0.5,color='blue',label=r'$B_{x}$')
ax35.plot(time[by_gsm<999],by_gsm[by_gsm<999],linewidth=0.5,color='green',label=r'$B_{y}$')
ax35.set_ylabel(r'$B_{x}$ & $B_{y}$ (nT)',labelpad=5,fontsize=7)
ax35.set_xticklabels(labels)
ax35.set_ylim([-22,22])
plt.legend(loc='upper right',fontsize=8)
ax35.set_yticks([-20,-10,0,10,20])
plt.yticks(size = 6)
plt.grid(True, axis='y')

ax3=plt.subplot(10,1,4)
plt.plot([dtdt(2023,3,23,7,43,21),dtdt(2023,3,23,7,43,21)],[-1000,1000],color='blue')
plt.plot([dtdt(2023,3,23,14,20,00),dtdt(2023,3,23,14,20,00)],[-1000,1000],color='orange')
plt.plot([dtdt(2023,3,23,17,40,00),dtdt(2023,3,23,17,40,00)],[-1000,1000],color='red')
plt.plot([dtdt(2023,3,24,8,20,00),dtdt(2023,3,24,8,20,00)],[-1000,1000],color='gray')
ax3.plot(time[phi<999],phi[phi<999],linewidth=0.5,color='lime',label=r'$\phi-180^{\circ}$')
ax3.plot(time[theta<999],theta[theta<999],linewidth=0.5,color='blue',label=r'$\theta$')
ax3.set_ylabel(r'$\phi$ & $\theta$ (Degree)',labelpad=2,fontsize=7)
ax3.set_xticklabels(labels)
ax3.set_ylim([-180,180])
plt.legend(loc='upper right',fontsize=8)
ax3.set_yticks([-180,-90,0,90,180])
plt.yticks(size = 6)
plt.grid(True, axis='y')

ax4=plt.subplot(10,1,5)
plt.plot([dtdt(2023,3,23,7,43,21),dtdt(2023,3,23,7,43,21)],[-1000,1000],color='blue')
plt.plot([dtdt(2023,3,23,14,20,00),dtdt(2023,3,23,14,20,00)],[-1000,1000],color='orange')
plt.plot([dtdt(2023,3,23,17,40,00),dtdt(2023,3,23,17,40,00)],[-1000,1000],color='red')
plt.plot([dtdt(2023,3,24,8,20,00),dtdt(2023,3,24,8,20,00)],[-1000,1000],color='gray')
ax4.plot(time[v_sw<999],v_sw[v_sw<999],linewidth=0.5,color='green')
ax4.set_ylabel(r'$v_{sw}$ (km s$^{-1}$)',labelpad=6,fontsize=7)
ax4.set_xticklabels(labels)
ax4.set_ylim([300,700])
ax4.set_yticks([400,450,500,550,600])
plt.yticks(size = 6)
plt.grid(True, axis='y')

ax5=plt.subplot(10,1,6)
plt.plot([dtdt(2023,3,23,7,43,21),dtdt(2023,3,23,7,43,21)],[-1000,1000],color='blue')
plt.plot([dtdt(2023,3,23,14,20,00),dtdt(2023,3,23,14,20,00)],[-1000,1000],color='orange')
plt.plot([dtdt(2023,3,23,17,40,00),dtdt(2023,3,23,17,40,00)],[-1000,1000],color='red')
plt.plot([dtdt(2023,3,24,8,20,00),dtdt(2023,3,24,8,20,00)],[-1000,1000],color='gray')
ax5.plot(time[density_p<99],density_p[density_p<99],linewidth=0.5,color='black')
#ax5.plot(time[pressure_kinetic<99],pressure_kinetic[pressure_kinetic<99],linewidth=0.5,color='red')
ax5.set_ylabel(r'$n_{p}$',labelpad=10,fontsize=7)
ax5.set_xticklabels(labels)
plt.ylim(-2,50)
plt.yticks(size = 6)
plt.grid(True, axis='y')
#ax5.set_ylim([-2,22])
#ax6=plt.subplot(10,1,6)
#plt.plot([dtdt(2023,3,23,7,43,21),dtdt(2023,3,23,7,43,21)],[-1000,1000],color='blue')
#plt.plot([dtdt(2023,3,23,14,20,00),dtdt(2023,3,23,14,20,00)],[-1000,1000],color='orange')
#ax6.plot(time[tem<999999],tem[tem<999999]/1e6,linewidth=0.5,color='black')#plot(time[density_p<99],density_p[density_p<99],linewidth=0.5,color='red')
#ax6.set_ylabel(r'$T_{p} (MK)$',labelpad=6,fontsize=7)
#ax6.set_xticklabels(labels)
#plt.yticks(size = 6)
#plt.ylim(-0.2,1)
#plt.grid(True, axis='y')
ax7=plt.subplot(10,1,7)
plt.plot([dtdt(2023,3,23,7,43,21),dtdt(2023,3,23,7,43,21)],[-1000,1000],color='blue')
plt.plot([dtdt(2023,3,23,14,20,00),dtdt(2023,3,23,14,20,00)],[-1000,1000],color='orange')
plt.plot([dtdt(2023,3,23,17,40,00),dtdt(2023,3,23,17,40,00)],[-1000,1000],color='red')
plt.plot([dtdt(2023,3,24,8,20,00),dtdt(2023,3,24,8,20,00)],[-1000,1000],color='gray')
ax7.plot(time[pressure_kinetic<99],pressure_kinetic[pressure_kinetic<99],linewidth=0.5,color='red')#plot(time[density_p<99],density_p[density_p<99],linewidth=0.5,color='red')
ax7.set_ylabel(r'$p_{k}$ (nPa)',labelpad=10,fontsize=7)
ax7.set_xticklabels(labels)
plt.yticks(size = 6)
plt.ylim(-2,25)
plt.grid(True, axis='y')

ax8=plt.subplot(10,1,8)
plt.plot([dtdt(2023,3,23,7,43,21),dtdt(2023,3,23,7,43,21)],[-1000,1000],color='blue')
plt.plot([dtdt(2023,3,23,14,20,00),dtdt(2023,3,23,14,20,00)],[-1000,1000],color='orange')
plt.plot([dtdt(2023,3,23,17,40,00),dtdt(2023,3,23,17,40,00)],[-1000,1000],color='red')
plt.plot([dtdt(2023,3,24,8,20,00),dtdt(2023,3,24,8,20,00)],[-1000,1000],color='gray')
ax8.plot(time[tem<999999],tem[tem<999999]/1e6,linewidth=0.5,color='black')#plot(time[density_p<99],density_p[density_p<99],linewidth=0.5,color='red')
ax8.set_ylabel(r'$T_{p} (MK)$',labelpad=6,fontsize=7)
ax8.set_xticklabels(labels)
plt.yticks(size = 6)
ax8.set_yticks([0,0.1,0.2,0.3,0.4,0.5])
plt.ylim(-0.1,0.5)
plt.grid(True, axis='y')
'''ax8.plot(time[beta<9],beta[beta<9],linewidth=0.5,color='blue')#plot(time[density_p<99],density_p[density_p<99],linewidth=0.5,color='red')
ax8.set_ylabel(r'$\beta$',labelpad=12,fontsize=7)
ax8.set_xticklabels(labels)
plt.yticks(size = 6)
plt.ylim(-1,8)
plt.grid(True, axis='y')'''
#ax9=plt.subplot(10,1,9)
#ax9.plot(time[m_a<99],m_a[m_a<99],linewidth=0.5,color='red')#plot(time[density_p<99],density_p[density_p<99],linewidth=0.5,color='red')
#ax9.set_ylabel(r'$M_A$',labelpad=18,fontsize=7)
#ax9.set_xticklabels(labels)
#plt.grid(True, axis='y')
'''ax11=plt.subplot(10,1,8)
plt.plot([dtdt(2023,3,23,7,43,21),dtdt(2023,3,23,7,43,21)],[-1000,1000],color='blue')
plt.plot([dtdt(2023,3,23,14,20,00),dtdt(2023,3,23,14,20,00)],[-1000,1000],color='orange')
plt.plot([dtdt(2023,3,23,17,40,00),dtdt(2023,3,23,17,40,00)],[-1000,1000],color='red')
ax11.plot(time[asy_d<999],asy_d[asy_d<999],linewidth=0.5,color='lime',label='ASY-D')
ax11.plot(time[asy_h<999],asy_h[asy_h<999],linewidth=0.5,color='cyan',label='ASY-H')#plot(time[density_p<99],density_p[density_p<99],linewidth=0.5,color='red')
ax11.set_ylabel(r'$ASY-D&H$ (nT)',labelpad=6,fontsize=7)
ax11.set_xticklabels(labels)
plt.legend(loc='best',fontsize=8)
plt.yticks(size = 6)
plt.ylim(-2,300)
plt.grid(True, axis='y')'''

Dst=np.loadtxt('DST.txt',dtype=str)
Dst=Dst[:5,2:-1].reshape(24*5).astype(int)
times_dst=np.array([dtdt(2023,3,22,1,0,0)+dttd(hours=i) for i in range(0,5*24,1)])
ax10=plt.subplot(10,1,9)
plt.plot([dtdt(2023,3,23,7,43,21),dtdt(2023,3,23,7,43,21)],[-1000,1000],color='blue')
plt.plot([dtdt(2023,3,23,14,20,00),dtdt(2023,3,23,14,20,00)],[-1000,1000],color='orange')
plt.plot([dtdt(2023,3,23,17,40,00),dtdt(2023,3,23,17,40,00)],[-1000,1000],color='red')
plt.plot([dtdt(2023,3,24,8,20,00),dtdt(2023,3,24,8,20,00)],[-1000,1000],color='gray')
ax10.plot(times_dst,Dst,linewidth=1.5,color='red',label=r'$Dst$')
ax10.plot(time[sym_h<999],sym_h[sym_h<999],linewidth=0.75,color='blue',label=r'$SYM-H$')#plot(time[density_p<99],density_p[density_p<99],linewidth=0.5,color='red')
ax10.set_ylabel(r'$Dst$ & $SYM-H$ (nT)',labelpad=0,fontsize=7)
ax10.set_xticklabels(labels)
ax10.set_ylim(-210,40)
ax10.set_yticks([-170,-150,-100,-50,0,25])
plt.legend(loc='upper right',fontsize=8)
plt.yticks(size = 6)
plt.grid(True, axis='y')

kp=np.loadtxt('kp_ap.txt')
kp=kp[:5,7:15].reshape(8*5)
times_new=np.array([dtdt(2023,3,22,1,30,0)+dttd(hours=i) for i in range(0,5*24,3)])
colors=np.zeros(kp.shape).astype(str)
colors[kp<=3]='lime'
colors[(kp>3)*(kp<=5)]='yellow'
colors[kp>5]='red'
ax12=plt.subplot(10,1,10)
plt.plot([dtdt(2023,3,23,7,43,21),dtdt(2023,3,23,7,43,21)],[-1000,1000],color='blue')
plt.plot([dtdt(2023,3,23,14,20,00),dtdt(2023,3,23,14,20,00)],[-1000,1000],color='orange')
plt.plot([dtdt(2023,3,23,17,40,00),dtdt(2023,3,23,17,40,00)],[-1000,1000],color='red')
plt.plot([dtdt(2023,3,24,8,20,00),dtdt(2023,3,24,8,20,00)],[-1000,1000],color='gray')
ax12.plot(time[sym_h<999],sym_h[sym_h<999],alpha=0)
ax12.bar(times_new,kp,width=0.1,color=colors)#plot(time[density_p<99],density_p[density_p<99],linewidth=0.5,color='red')
ax12.set_ylabel(r'$K_p$',labelpad=8,fontsize=7)
ax12.set_ylim(0,10)
ax12.set_yticks([0,1,2,3,4,5,6,7,8,9,10])
plt.yticks(size = 6)
plt.grid(True, axis='y')


fig.subplots_adjust(0.1,0.05,0.97,0.95,0,0.08)
plt.savefig('IMF_solar_wind_L1.pdf')#,dpi=700)
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


bx=bx_gse[1900:3380]
by=by_gse[1900:3380]
bz=bz_gse[1900:3380]
bx=bx[bx<999]
by=by[by<999]
bz=bz[bz<999]
#moving average
order=151
bx = np.convolve(bx, np.ones((order,))/order, mode='valid')
by = np.convolve(by, np.ones((order,))/order, mode='valid')
bz = np.convolve(bz, np.ones((order,))/order, mode='valid')
plt.plot(bx,color='blue')
plt.plot(by,color='green')
plt.plot(bz,color='red')
bxyz=np.vstack((bx,by,bz))
b,angle,RMS,eigenValues,eigenVectors=minvarRotate(bxyz)
exmax,eymax,ezmax=eigenVectors[:,1]

theta=np.arctan(ezmax/eymax)/np.pi*180.
print('flux rope tilt angle:',theta)