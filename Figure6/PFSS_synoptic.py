"""
HMI PFSS solutions
------------------
Calculating a PFSS solution from a HMI synoptic map.

This example shows how to calcualte a PFSS solution from a HMI synoptic map.
There are a couple of important things that this example shows:

- HMI maps have non-standard metadata, so this needs to be fixed
- HMI synoptic maps are very big (1440 x 3600), so need to be downsampled
  in order to calculate the PFSS solution in a reasonable time.
"""
import os

import astropy.units as u
import matplotlib.pyplot as plt
import sunpy.map
from sunpy.net import Fido
from sunpy.net import attrs as a

import pfsspy
import pfsspy.utils
from scipy.interpolate import interp1d
from datetime import datetime as dtdt

###############################################################################
# Set up the search.
#
# Note that for SunPy versions earlier than 2.0, a time attribute is needed to
# do the search, even if (in this case) it isn't used, as the synoptic maps are
# labelled by Carrington rotation number instead of time
'''time = a.Time('2023/03/19', '2023/03/20')
series = a.jsoc.Series('hmi.synoptic_mr_polfil_720s')
crot = a.jsoc.PrimeKey('CAR_ROT', 2268)

###############################################################################
# Do the search.
#
# If you use this code, please replace this email address
# with your own one, registered here:
# http://jsoc.stanford.edu/ajax/register_email.html
result = Fido.search(time, series, crot,
                     a.jsoc.Notify(os.environ[r"TwinWilling@protonmail.com"]))
files = Fido.fetch(result,path="~/Downloads/pfss.fits")
gogo'''
###############################################################################
# Read in a file. This will read in the first file downloaded to a sunpy Map
# object
hmi_map = sunpy.map.Map('hmi.synoptic_mr_polfil_720s.2268.Mr_polfil.fits')
print('Data shape: ', hmi_map.data.shape)

###############################################################################
# Since this map is far to big to calculate a PFSS solution quickly, lets
# resample it down to a smaller size.
ntheta=180
nphi=360
hmi_map = hmi_map.resample([nphi, ntheta] * u.pix)
print('New shape: ', hmi_map.data.shape)

###############################################################################
# Now calculate the PFSS solution
nrho = 192
rss = 2.5
pfss_in = pfsspy.Input(hmi_map, nrho, rss)
pfss_out = pfsspy.pfss(pfss_in)

###############################################################################
# Using the Output object we can plot the source surface field, and the
# polarity inversion line.
ss_br = pfss_out.source_surface_br
# Create the figure and axes
fig = plt.figure()
ax = plt.subplot(projection=ss_br)

# Plot the source surface map
ss_br.plot()
# Plot the polarity inversion line
ax.plot_coord(pfss_out.source_surface_pils[0])
# Plot formatting
plt.colorbar()
ax.set_title('Source surface magnetic field')

plt.show()

import numpy as np
#tstart=hmi_map.fits_header['T_START']
#tstop=hmi_map.fits_header['T_STOP']
#t_research=dtdt(2023,3,19,16,30,0)
#t_start=dtdt(int(tstart[:4]),int(tstart[5:7]),int(tstart[8:10]),int(tstart[11:13]),int(tstart[14:16]),int(tstart[17:19]))
#t_stop=dtdt(int(tstop[:4]),int(tstop[5:7]),int(tstop[8:10]),int(tstop[11:13]),int(tstop[14:16]),int(tstop[17:19]))
#cut=int((t_research-t_start).total_seconds()/(t_stop-t_start).total_seconds()*nphi)
cut=int(61.6)
def deriv(x,y,h):
    try:
        delta=(y(x+h)-y(x))/h
    except:
        delta=(y(x)-y(x-h))/h
    return delta
R_sun=696300e3*1e2#cm
rs=np.linspace(0,(rss-1)*R_sun,nrho+1)
log_r=np.log(rs)
pfss_dat=np.array(pfss_out.bg)
bs=np.sqrt(pfss_dat[:,:,:,0]**2+pfss_dat[:,:,:,1]**2)[:-1,:,:]#Gauss
log_b=np.log(bs)
f_log_b=interp1d(log_r,log_b)
n=-deriv(log_r,f_log_b,1e-6)
#d_logb=log_b[:,:,1:]-log_b[:,:,:-1]
#n=-d_logb/d_logr
plt.imshow(n[:,:,50].T,vmax=2,vmin=-2)
plt.colorbar()
plt.show()
plt.imshow(n[cut,:,::-1].T,vmax=2,vmin=-2)
plt.colorbar()
plt.show()
n2=n.reshape(nphi,(ntheta+1)*(nrho+1))
np.savetxt('decay_index.dat',n2)
bx=pfss_dat[:,:,:,0][:-1,:,:].reshape(nphi,(ntheta+1)*(nrho+1))
by=pfss_dat[:,:,:,1][:-1,:,:].reshape(nphi,(ntheta+1)*(nrho+1))
bz=pfss_dat[:,:,:,2][:-1,:,:].reshape(nphi,(ntheta+1)*(nrho+1))
np.savetxt('bx_synoptic.dat',bx)
np.savetxt('by_synoptic.dat',by)
np.savetxt('bz_synoptic.dat',bz)