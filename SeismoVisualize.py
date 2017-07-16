from obspy import read
from obspy.clients.fdsn import Client
from obspy import UTCDateTime
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D
from obspy.taup import getTravelTimes
from math import degrees, radians, cos, sin, atan2, sqrt, floor
import argparse
from sys import platform as _platform

# TODO
# - Shade surface wave arrivals as an option
# - Add automatic scale bar
# - Add coordinate rotation and relabel of plots


def GetData(t0, net, st0, loc, ch, duration):
    """
    Download data from the IRIS datacenter and output
    with the instrument response removed and calibrated.
    A filter is also placed. Return a station object.
    """
    client = Client('IRIS')
    st = client.get_waveforms(net, st0, loc, ch, t0,
                              t0+duration*60, attach_response=True)
    st.detrend(type='demean')
    st.detrend(type='linear')
    st.taper(max_percentage=0.05)
    st.remove_response(output='DISP')
    st.filter('highpass', freq=0.01, corners=4, zerophase=True)
    return st


def GetStationLocation(t0, net, st0, loc, duration):
    """
    Given a time, duration, loc code, and station network/name, get
    station information from IRIS.  Return a list containing the
    lat, lon, and elevation.
    """
    client = Client('IRIS')
    st0 = client.get_stations(starttime=t0, endtime=t0+duration*60,
                              network=net, station=st0, level='station')
    slat = st0[0][0].latitude
    slon = st0[0][0].longitude
    selev = st0[0][0].elevation
    return [slat, slon, selev]


def MarkPhase(ax, phase, travel_times, yloc, fontsize=16):
    """
    Mark a phase with a letter on the plot.  The letter
    is partially transparent until the phase arrives.
    """
    tx = None
    if phase in travel_times.keys():
        tarr = travel_times[phase]
        tx = ax.text(tarr, yloc, phase, fontsize=fontsize, alpha=0.3, ha='center')
    else:
        tarr = None
    return tx


def UpdatePhaseMarker(marker, travel_time, cur_time):
    if cur_time >= travel_time:
        marker.set_alpha(1.0)
        return marker
    else:
        return marker


def GetTravelTimes(station, earthquake):
    """
    Calculate travel times for phases using obspy and reformat
    the output to be a dictionary with phase name as key and
    arrival time as the value.
    """
    dist = DegreesDistance(station[0], station[1],
                           earthquake[0], earthquake[1])

    tt = getTravelTimes(dist, earthquake[2])
    travel_times = {}
    for item in tt:
        travel_times[item['phase_name']] = item['time']
    return travel_times


def DegreesDistance(lat1, lon1, lat2, lon2):
    """
    Calcaulate the distance in degrees between two
    lat/lon pairs.
    """
    lat1 = radians(lat1)
    lat2 = radians(lat2)
    lon1 = radians(lon1)
    lon2 = radians(lon2)
    dLon = abs(lon2-lon1)

    arg1 = (cos(lat2)*sin(dLon))**2
    arg2 = (cos(lat1)*sin(lat2)-sin(lat1)*cos(lat2)*cos(dLon))**2
    arg3 = sin(lat1)*sin(lat2)+cos(lat1)*cos(lat2)*cos(dLon)

    dist = atan2(sqrt(arg1+arg2), arg3)
    return degrees(dist)


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
    xz_points.set_data(x_data, np.ones_like(x_data)*ax_lims)
    xz_points.set_3d_properties(z_data)

    yz_points.set_data(np.ones_like(y_data)*-ax_lims, y_data)
    yz_points.set_3d_properties(z_data)

    xy_points.set_data(x_data, y_data)
    xy_points.set_3d_properties(np.ones_like(z_data)*-ax_lims)

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

    # Plot phase arrival text markers
    for key in phase_markers:
        phase = phase_markers[key]
        UpdatePhaseMarker(phase[0], phase[1], cur_time)

    #
    # Set figure text
    #
    figtitle_text.set_text('%.1f Seconds After Earthquake' % t[-1])

    #
    # Set text for current position
    #
    ch1value_text.set_text('E/W Position: %8.3f mm' % x_data[-1])
    ch2value_text.set_text('N/S Position: %8.3f mm' % y_data[-1])
    ch3value_text.set_text('U/D Position: %8.3f mm' % z_data[-1])
    return (marker_line, s3d, xz_points, yz_points, xy_points,
            x_marker, y_marker, z_marker)

#
# Set parameters and parse arguments
#
parser = argparse.ArgumentParser(description='Produce a visualization of seismometer data.',
                                 epilog='Example: python SeismoVisualize.py -n IU -s ADK -l 10 -t 2014-07-07T11:23:58 -d 60 -c BH1 BH2 BHZ -p P S -e 14.782 -92.371 92.')

parser.add_argument('-s', required=True, type=str, help='Station code of the seismometer')
parser.add_argument('-n', required=True, type=str, help='Network code of the seismometer')
parser.add_argument('-l', required=True, type=str, help='Location code of the seismometer')
parser.add_argument('-c', required=True, type=str, nargs=3, help='Channels listed in x,y,z order (E/W,N/S,U/D)')
parser.add_argument('-e', required=True, type=float, nargs=3, help='Earthquake latitude, longitude, depth')
parser.add_argument('-d', required=True, type=int, help='Duration of data to plot (in minutes)')
parser.add_argument('-t', required=True, type=str, help='Start time in format: YYYY-MM-DDTTT:MM:SS)')
parser.add_argument('-p', required=True, type=str, nargs='+', help='Phases to plot')

parser.add_argument('-label', default=14, help='Label size for the plot')
parser.add_argument('-tick', default=12, help='Tick size for the plot')
parser.add_argument('-trail', default=15, help='Number of points in the tail of the plot')

args = parser.parse_args()
args.t = UTCDateTime(args.t)
trail = args.trail
labelsize = args.label
ticksize = args.tick

#
# Get Station Information
#
print ('Fetching Station Information...')
station_info = GetStationLocation(args.t, args.n, args.s, args.l, args.d)
print ('Complete\n')


#
# Calculate phase arrival times
#
print ('Calculating travel times...')
travel_times = GetTravelTimes(station_info, args.e)
print ('Complete\n')

#
# Read data and combine into stream object, channels should be in
# order for x,y,z axes. Generall that is E-W,N-S,U-D or
# T,R,Z for rotated systems
#
print ('Downloading station data:')
print ('Ch.1')
st1 = GetData(args.t, args.n, args.s, args.l, args.c[0], args.d)
print ('Ch.2')
st2 = GetData(args.t, args.n, args.s, args.l, args.c[1], args.d)
print ('Ch.3')
st3 = GetData(args.t, args.n, args.s, args.l, args.c[2], args.d)
st = st1 + st2 + st3
print ('Complete\n')


#
# Process the traces
#
st.detrend(type='demean')
st.detrend(type='linear')
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

#
#
# Make the Plot
#
#

print ('Beginning plotting...')

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

ax1.set_xticklabels("")
ax1.set_yticklabels("")
ax1.set_zticklabels("")

ax1.set_xlabel('East - West',fontsize=labelsize-2)
ax1.set_ylabel('North - South',fontsize=labelsize-2)
ax1.set_zlabel('Up - Down',fontsize=labelsize-2)

#
# Setup Text on Figure
#
figtitle_text = plt.suptitle('', fontsize=labelsize+5)

# Static labels on the 2D plot
ch1label_text = ax2.text(0.75*max(time), 0.1*offset,
                         'East - West', fontsize=labelsize-2)

ch2label_text = ax2.text(0.75*max(time), 1.1*offset,
                         'North - South', fontsize=labelsize-2)

ch3label_text = ax2.text(0.75*max(time), 2.1*offset,
                         'Up - Down', fontsize=labelsize-2)

#
# Make the phase marker dictionary
#
phase_markers = {}
for phase in args.p:
    phase_text = MarkPhase(ax2, phase, travel_times, max(st[2])+2.1*offset)
    if phase_text != None:
        print ('Adding phase {0:s} to plot at time {1:.2f}'.format(phase,travel_times[phase]))
        phase_markers[phase] = [phase_text, travel_times[phase]]
    else:
        print ('Phase {0:s} not found in predictions'.format(phase))

# Lables in the text axes
x_text_loc = 0.2
y_text_loc = 0.8
y_text_offset = 0.1
station_text = ax3.text(x_text_loc, y_text_loc-0*y_text_offset,
                        '', transform=ax3.transAxes)

ch3value_text = ax3.text(x_text_loc, y_text_loc-1*y_text_offset,
                         '', transform=ax3.transAxes)

ch2value_text = ax3.text(x_text_loc, y_text_loc-2*y_text_offset,
                         '', transform=ax3.transAxes)

ch1value_text = ax3.text(x_text_loc, y_text_loc-3*y_text_offset,
                         '', transform=ax3.transAxes)

ch3maxvalue_text = ax3.text(x_text_loc, y_text_loc-4*y_text_offset,
                         '', transform=ax3.transAxes)

ch2maxvalue_text = ax3.text(x_text_loc, y_text_loc-5*y_text_offset,
                         '', transform=ax3.transAxes)

ch1maxvalue_text = ax3.text(x_text_loc, y_text_loc-6*y_text_offset,
                         '', transform=ax3.transAxes)

ch1maxvalue_text.set_text('E-W Maximum Dispacement: %5.2f mm' % max(abs(st[0].data)))
ch2maxvalue_text.set_text('N-S Maximum Dispacement: %5.2f mm' % max(abs(st[1].data)))
ch3maxvalue_text.set_text('U-D Maximum Dispacement: %5.2f mm' % max(abs(st[2].data)))

station_text.set_text('Station %s' % args.s)


#
# 2D Seismogram Plot
#

# Setup Labels
ax2.set_xlabel('Seconds Since Earthquake', fontsize=labelsize)
ax2.set_ylabel('Motion', fontsize=labelsize)
ax2.tick_params(axis='both', which='major', labelsize=ticksize)
ax2.spines['right'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax2.yaxis.set_ticks_position('left')
ax2.xaxis.set_ticks_position('bottom')
ax2.axes.get_yaxis().set_visible(False)
ax2.text(0, -0.2, 'www.leemangeophysical.com',
         fontsize=labelsize-4, transform=ax2.transAxes)

# Plot 2D static seismogram data
for (i, tr) in enumerate(st):
    ax2.plot(time, tr.data+offset * i, color='k', linewidth=1)

# Set limits
ax2.set_xlim(min(time), max(time))
ax2.set_ylim(-offset*.75, 2.75*offset)

# Place the vertical time scroll bar
marker_line = ax2.axvline(x=0, color='r', linewidth=2)

#
# Plot points/bars in 3D scatter
#

# Maker bars to current point
x_marker, = ax1.plot([], [], [], color='k')
y_marker, = ax1.plot([], [], [], color='k')
z_marker, = ax1.plot([], [], [], color='k')

# Main 3D scatter
s3d, = ax1.plot([], [], [], marker='o', linestyle='None', color='r')

# Points on plane projections
xz_points, = ax1.plot([], [], [], color='k', marker='o',
                      linestyle='None', markersize=3)

xy_points, = ax1.plot([], [], [], color='k', marker='o',
                      linestyle='None', markersize=3)

yz_points, = ax1.plot([], [], [], color='k', marker='o',
                      linestyle='None', markersize=3)

ax1.set_xlim3d(-1*ax_lims, ax_lims)
ax1.set_ylim3d(-1*ax_lims, ax_lims)
ax1.set_zlim3d(-1*ax_lims, ax_lims)

#
# Do the animation
#
inds = np.arange(0, len(time))

anim = FuncAnimation(fig, step, frames=inds, interval=50,
                     repeat_delay=2000, blit=True)

fname = '%s_%s_%s.mp4' % (args.n, args.s, str(args.t))
if (_platform == 'win32'):
    fname = fname.replace(':', '-')
print('Writing to {0}'.format(fname));
anim.save(fname, bitrate=2500)
plt.close('all')
