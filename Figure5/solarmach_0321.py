from solarmach import SolarMACH, print_body_list
import matplotlib.pyplot as plt
from pandas import *
import numpy as np

ftype='.pdf'
# optional: get list of available bodies/spacecraft
print(print_body_list().index)

# necessary options
body_list = ['STEREO-A', 'Earth', 'Solar Orbiter']
date = '2023-03-21 09:30:00'

# Previously you needed to define position-sensitive solar wind speed per
# body in body_list, e.g., vsw_list = [400, 400, 400, 400, 400, 400, 400]
# Now you can skip this parameter or provide an empty list. Then solarmach
# will try to automatically obtain measured solar wind speeds from each
# spacecraft
vsw_list = [500]*3

# optional parameters
coord_sys = 'Stonyhurst'                         # 'Carrington' (default) or 'Stonyhurst'
reference_long = 0.00262814-360./365.*(1*24+11.)/24  #GCS Model Shell Axis Longitude at Lat=0.                           # longitude of reference (None to omit)
reference_lat = -7.0                                # latitude of reference (None to omit)
plot_spirals = False                             # plot Parker spirals for each body
plot_sun_body_line = True                        # plot straight line between Sun and body
long_offset = 270+360./365.*(1*24+11.)/24                               # longitudinal offset for polar plot; defines where Earth's longitude is (by default 270, i.e., at "6 o'clock")
reference_vsw = 3e5                             # define solar wind speed at reference
return_plot_object = False                       # figure and axis object of matplotib are returned, allowing further adjustments to the figure
transparent = False                              # make output figure background transparent
markers = 'letters'                              # use 'numbers' or 'letters' for the body markers (use False for colored squares)
filename = 'Solar-MACH_'+date.replace(' ', '_')  # define filename of output figure

# initialize
sm = SolarMACH(date, body_list, vsw_list, reference_long, reference_lat, coord_sys)
#sm = SolarMACH(date, body_list, vsw_list, coord_sys)

coord=np.vstack([np.array(sm.coord_table.columns),sm.coord_table.values])
np.savetxt('positions.txt',coord,fmt='%.100s',delimiter=' ; ')
# make plot
#fig=plt.figure(figsize=(2,1))
plt.rcParams.update({"font.size":20,'font.family':"Arial",'font.weight':'bold'})
sm.plot(
   long_sector=[reference_long-19.27,reference_long+19.27],
   long_sector_vsw=None,
   long_sector_color='cyan',
   plot_spirals=plot_spirals,
   plot_sun_body_line=plot_sun_body_line,
   reference_vsw=reference_vsw,
   transparent=transparent,
   markers=markers,
   long_offset=long_offset,
   return_plot_object=return_plot_object,
   outfile=filename+ftype,
   figsize=(10,6)
)

# obtain data as Pandas DataFrame
display(sm.coord_table)
sm.coord_table.to_csv('Loc_Obs.dat',index=True)