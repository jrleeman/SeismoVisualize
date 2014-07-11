from obspy import read
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D


def step(ind):
    # Slice the data to plot
    x_data = st[0].data[ind:ind+trail]
    y_data = st[1].data[ind:ind+trail]
    z_data = st[2].data[ind:ind+trail]
    t = time[ind:ind+trail]
    cur_time = t[-1]

    # Set 3D Data
    s3d.set_data(x_data, y_data)
    s3d.set_3d_properties(z_data)

    # Set Marker Position
    marker_line.set_xdata(cur_time)

    # Set points projection data
    xz_points.set_data(x_data, np.ones([trail])*ax_lims)
    xz_points.set_3d_properties(z_data)

    yz_points.set_data(np.ones([trail])*-ax_lims, y_data)
    yz_points.set_3d_properties(z_data)

    xy_points.set_data(x_data, y_data)
    xy_points.set_3d_properties(np.ones([trail])*-ax_lims)

    # Set the marker bars
    x = [-ax_lims, x_data[-1]]
    y = [y_data[-1], y_data[-1]]
    z = [z_data[-1], z_data[-1]]
    x_marker.set_data(x, y)
    x_marker.set_3d_properties(z)

    x = [x_data[-1], x_data[-1]]
    y = [ax_lims, y_data[-1]]
    z = [z_data[-1], z_data[-1]]
    y_marker.set_data(x, y)
    y_marker.set_3d_properties(z)

    x = [x_data[-1], x_data[-1]]
    y = [y_data[-1], y_data[-1]]
    z = [-ax_lims, z_data[-1]]
    z_marker.set_data(x, y)
    z_marker.set_3d_properties(z)

    #
    # Set figure text
    #
    figtitle_text.set_text('%.1f Seconds After Earthquake' % t[-1])

    return marker_line, s3d, xz_points, yz_points, xy_points, x_marker, y_marker, z_marker

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
st.detrend(type='demean')
st.detrend(type='linear')
st.integrate(type='cumtrapz', initial=0)
st.decimate(factor=int(st[0].stats.sampling_rate), no_filter=True)

# Scale to mm
for tr in st:
    tr.data = tr.data * 1000

#
# Make a time array
#
time = np.arange(st[0].stats.npts)
time = time / st[0].stats.sampling_rate

#
# Add zero padding for plotting
#
time = np.concatenate((np.zeros([trail-1]), time))

for tr in st:
    tr.data = np.concatenate((np.zeros([trail-1]), tr.data))

##################
##################
# PLOT
##################
##################


#
# Determine the limits and offsets for the 3D plot
#
ax_lims = 0
offset = 0
for tr in st:
    tr_max = max(abs(tr.data))
    tr_amp = abs(max(tr.data) - min(tr.data))

    if tr_max * 1.3 > ax_lims:
        ax_lims = tr_max * 1.3

    if tr_amp > offset:
        offset = tr_amp

#
# Setup Axes
#
fig = plt.figure(figsize=(9, 9))
ax1 = plt.subplot(2, 2, 1, projection='3d')  # 3D Motion Plot
ax2 = plt.subplot(2, 1, 2)                   # Seismogram Plot
ax3 = plt.subplot(2, 2, 2)                   # Text Area Plot
ax3.axis('off')  # Turn off the border for text area

figtitle_text = plt.suptitle('', fontsize=labelsize+5)

#
# 2D Seismogram Plot
#

# Setup Labels
ax2.set_xlabel(r'Seconds Since Earthquake', fontsize=labelsize)
ax2.set_ylabel(r'Motion', fontsize=labelsize)
ax2.tick_params(axis='both', which='major', labelsize=ticksize)
ax2.axes.get_yaxis().set_visible(False)
ax2.text(0, -0.2, 'www.johnrleeman.com',
         fontsize=labelsize-4, transform=ax2.transAxes)

# Plot 2D static seismogram data
for (i, tr) in enumerate(st):
    ax2.plot(time, tr.data+offset * i, color='k')

# Set limits
ax2.set_xlim(min(time), max(time))

# Place the vertical time scroll bar
marker_line = ax2.axvline(x=0, color='r', linewidth=2)

#
# Plot points/bars in 3D scatter
#

# Main 3D scatter
s3d, = ax1.plot([], [], [], marker='o', linestyle='None')

# Points on plane projections
xz_points, = ax1.plot([], [], [], color='k', marker='o',
                      linestyle='None', markersize=3)

xy_points, = ax1.plot([], [], [], color='k', marker='o',
                      linestyle='None', markersize=3)

yz_points, = ax1.plot([], [], [], color='k', marker='o',
                      linestyle='None', markersize=3)

# Maker bars to current point
x_marker, = ax1.plot([], [], [], color='k')
y_marker, = ax1.plot([], [], [], color='k')
z_marker, = ax1.plot([], [], [], color='k')

ax1.set_xlim3d(-1*ax_lims, ax_lims)
ax1.set_ylim3d(-1*ax_lims, ax_lims)
ax1.set_zlim3d(-1*ax_lims, ax_lims)

#
# Do the animation
#
inds = np.arange(0, len(time))

# Reduce size for testing
# inds = np.arange(500,1000)

anim = FuncAnimation(fig, step, frames=inds, interval=50,
                     repeat_delay=2000, blit=True)

anim.save('test_scipy.mp4')
