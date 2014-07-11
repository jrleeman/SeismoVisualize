from obspy import read
from scipy.integrate import cumtrapz
from obspy.taup import getTravelTimes
from math import degrees,radians,cos,sin,atan2,sqrt,floor
from obspy.fdsn import Client
from obspy import UTCDateTime
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation

def step(ind):
    # Slice the data to plot
    x_data =  st[0].data[ind:ind+trail]
    y_data =  st[1].data[ind:ind+trail]
    z_data =  st[2].data[ind:ind+trail]
    t      =  time[ind:ind+trail]
    cur_time = t[-1]
    
    # Set 3D Data
    s3d.set_data(x_data, y_data)
    s3d.set_3d_properties(z_data)
    
    # Set Marker Position
    marker_line.set_xdata(cur_time)

    return marker_line,s3d

labelsize = 14
ticksize = 12
trail = 10

#
# Read data and combine into stream object
#
st1 = read('example_data/2014-07-07T11-23-58.IU.ANMO.10.BH1.sac')
st2 = read('example_data/2014-07-07T11-23-58.IU.ANMO.10.BH2.sac')
st3 = read('example_data/2014-07-07T11-23-58.IU.ANMO.10.BHZ.sac')
st = st1 + st2 + st3

# 
# Process the traces
#
st.detrend(type = 'demean')
st.detrend(type = 'linear')
st.integrate(type = 'cumtrapz') # Reduces the length by 1 data point
st.decimate(factor=int(st[0].stats.sampling_rate),no_filter=True)

# Scale to mm
for tr in st:
    tr.data = tr.data * 1000

#
# Make a time array
#
time = np.arange(st[0].stats.npts) + 1
time = time / st[0].stats.sampling_rate

#
# Add zero padding for plotting
#
time = np.concatenate((np.zeros([trail-1]),time))

for tr in st:
    tr.data = np.concatenate((np.zeros([trail-1]),tr.data))


##################
##################
# PLOT
##################
##################


#
# Determine the limits and offsets for the 3D plot
#
ax_lims = 0
offset  = 0
for tr in st:
    tr_max = max(abs(tr.data))
    tr_amp = abs(max(tr.data) - min(tr.data))
    
    if tr_max * 1.2 > ax_lims:
        ax_lims = tr_max * 1.2
        
    if tr_amp > offset:
        offset = tr_amp 

#
# Setup Axes
#
fig = plt.figure(figsize=(9,9))
ax1 = plt.subplot(2, 2, 1,projection='3d') # 3D Motion Plot
ax2 = plt.subplot(2, 1, 2)                 # Seismogram Plot
ax3 = plt.subplot(2, 2, 2)                 # Text Area Plot
ax3.axis('off') # Turn off the border for text area

figtitle_text = plt.suptitle('',fontsize=labelsize+5) # Give the plot a title

#
# 2D Seismogram Plot
#

# Setup Labels
ax2.set_xlabel(r'Seconds Since Earthquake',fontsize=labelsize)
ax2.set_ylabel(r'Motion',fontsize=labelsize)
ax2.tick_params(axis='both', which='major', labelsize=ticksize)
ax2.axes.get_yaxis().set_visible(False)
ax2.text(0,-0.2,'www.johnrleeman.com',fontsize=labelsize-4,transform = ax2.transAxes)

# Plot 2D static seismogram data
for (i,tr) in enumerate(st):
    ax2.plot(time,tr.data + offset * i, color = 'k')

# Set limits
ax2.set_xlim(min(time),max(time))

# Place the vertical time scroll bar
marker_line = ax2.axvline(x=0,color='r',linewidth=2)

#
# Plot points/bars in 3D scatter
#

# Main 3D scatter
s3d, = ax1.plot([], [], [], marker='o', linestyle='None')

# Points on plane projections
xz_points = ax1.plot([], [], [],color='k',marker='o', linestyle='None')
xy_points = ax1.plot([], [], [],color='k',marker='o', linestyle='None')
yz_points = ax1.plot([], [], [],color='k',marker='o', linestyle='None')

# Maker bars to current point
x_marker = ax1.plot([], [], [],color='k')
y_marker = ax1.plot([], [], [],color='k')
z_marker = ax1.plot([], [], [],color='k')

ax1.set_xlim3d(-1*ax_lims,ax_lims)
ax1.set_ylim3d(-1*ax_lims,ax_lims)
ax1.set_zlim3d(-1*ax_lims,ax_lims)
    

#
# Do the animation
#
inds = np.arange(0,len(time))
anim = FuncAnimation(fig, step, frames=inds, interval=50, repeat_delay=2000,blit=True)
anim.save('test_scipy.mp4')

#plt.show()