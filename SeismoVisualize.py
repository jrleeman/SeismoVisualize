from obspy import read
from scipy.integrate import cumtrapz
from scipy.signal import detrend
from scipy import signal
import os
from obspy.taup import getTravelTimes
from math import degrees,radians,cos,sin,atan2,sqrt,floor
from obspy.fdsn import Client
from obspy import UTCDateTime
import sys

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation

########## TODO
# Get event data given time,window,mag cutoff
# x = e/w, y = n/s, z = u/d
def GetData(t0,net,st0,loc,ch,duration):
    """
    Download data from the IRIS datacenter and output
    with the instrument response removed and calibrated.
    Return a station object.
    """
    client = Client("IRIS")
    st = client.get_waveforms(net, st0, loc, ch, t0, t0+duration*60,attach_response=True)
    #st.remove_response(output="VEL")
    return st
    
def GetStationLocation(t0,net,st0,loc,duration):
    """
    Given a time, duration, loc code, and station network/name, get
    station information from IRIS.  Return a list containing the
    lat, lon, and elevation.
    """    
    client = Client("IRIS")
    st0 = client.get_stations(starttime=t0,endtime=t0+duration*60,network=net,station=st0,level='station')
    slat = st0[0][0].latitude
    slon = st0[0][0].longitude
    selev = st0[0][0].elevation
    return [slat,slon,selev]
    
def DegreesDistance(lat1,lon1,lat2,lon2):
    """
    Calcaulate the distance in degrees between two
    lat/lon pairs.
    """
    lat1 = radians(lat1)
    lat2 = radians(lat2)
    lon1 = radians(lon1)
    lon2 = radians(lon2)
    dLat = abs(lat2-lat1)
    dLon = abs(lon2-lon1)
    
    arg1 = (cos(lat2)*sin(dLon))**2
    arg2 = (cos(lat1)*sin(lat2)-sin(lat1)*cos(lat2)*cos(dLon))**2
    arg3 = sin(lat1)*sin(lat2)+cos(lat1)*cos(lat2)*cos(dLon)
    
    dist = atan2(sqrt(arg1+arg2),arg3)
    return degrees(dist)

def ProcessTrace(data,Fs):
    """
    Take a data array and detrend, remove the mean, and
    integrate.  Multiply by 1000 to get units in mm.
    """
    data = detrend(data)
    data = data - np.mean(data)
    data = cumtrapz(data,dx=(1./Fs),initial=0) * 1000.
    return data

def MarkPhase(ax,phase,t,travel_times,yloc,fontsize=20,alpha=0.3):   
    """ 
    Mark a phase with a letter on the plot.  The letter
    is partially transparent until the phase arrives.
    """ 
    
    if phase in travel_times.keys():
        tarr = travel_times[phase]
        if t >= tarr:
            ax.text(tarr,yloc,phase,fontsize=fontsize)
        else:
            ax.text(tarr,yloc,phase,fontsize=fontsize,alpha=alpha)
    else:
        tarr = None
        
        
def GetTravelTimes(station,earthquake):
    """
    Calculate travel times for phases using obspy and reformat
    the output to be a dictionary with phase name as key and
    arrival time as the value.
    """
    dist = DegreesDistance(station[0],station[1],earthquake[0],earthquake[1])
    tt = getTravelTimes(dist,earthquake[2])
    travel_times={}
    for item in tt:
        travel_times[item['phase_name']] = item['time']
    return travel_times

class Seismogram:
    def __init__(self):
        self.channels = 
        self.time = None
        self.ndecimate = 40
        self.Fs = 40.
        
    def CleanData(self):
        chs = [self.ch1,self.ch2,self.ch3]
        print self.ch1[1000]
        for ch in chs:
            ch = self.todisp(ch,self.Fs)
            ch = self.decimate(ch,self.ndecimate)
            ch = self.scale(ch,1000)
            ch = self.remove_mean(ch)
            ch = self.remove_trend(ch)
        print self.ch1[1000]
        return self
    
    class channel:
        def __init__(self):
            self.label = None
            self.data = None
            self.Fs = None
            
        def decimate(self,data,dec):
            return data[::dec]
        
        def remove_mean(self,data):
            return data - np.mean(data)
            
        def scale(self,data,factor):
            return data * factor
        
        def remove_trend(self,data):
            return detrend(data)
            
        def get_scatterdata(self,i):
            time = self.time[i-10+1:i+1]
            ch1_scatter = self.ch1[i-10+1:i+1]
            ch2_scatter = self.ch2[i-10+1:i+1]
            ch3_scatter = self.ch3[i-10+1:i+1]
            return time,ch1_scatter,ch2_scatter,ch3_scatter
    
        def todisp(self,data,Fs):
            return cumtrapz(data,dx=(1./Fs),initial=0)
    

        
    
    
###############################################################################
# Change these parameters to modify the plotting behavior

trail = 10
labelsize = 14
ticksize  = 12
###############################################################################

### Parse call
# python SeismoVisualize.py NETWORK STATION LOCATION CHN CHE CHZ TIME DURATION
try:
    network = sys.argv[1]
    stationcd = sys.argv[2]
    location = sys.argv[3]
    chn = sys.argv[4]
    che = sys.argv[5]
    chz = sys.argv[6]
    time_str = sys.argv[7]
    duration = int(sys.argv[8])
    time = UTCDateTime(time_str)
    print network,stationcd,location,chn,che,chz,time_str,duration
    
except:
    print "\n\nINVALID CALL!"
    print "python SeismoVisualize.py NETWORK STATION LOCATION CHN CHE CHZ TIME DURATION\n"
    print "python SeismoVisualize.py IU ANMO 10 BH1 BH2 BHZ 2014-07-07T11:23:58 60\n\n"
    sys.exit()
    

seismo = Seismogram()

station = GetStationLocation(time,network,stationcd,location,duration)
earthquake = [14.782,-92.371,92.]

#st0 = GetData(time,network,stationcd,location,'BH1',duration)
#st1 = GetData(time,network,stationcd,location,'BH2',duration)
#st2 = GetData(time,network,stationcd,location,'BHZ',duration)

# Temporary to keep us from downloading data all the time
st1 = read('example_data/2014-07-07T11-23-58.IU.ANMO.10.BH1.sac')
st2 = read('example_data/2014-07-07T11-23-58.IU.ANMO.10.BH2.sac')
st3 = read('example_data/2014-07-07T11-23-58.IU.ANMO.10.BHZ.sac')

# Calculate travel times
travel_times = GetTravelTimes(station,earthquake)

# Get the sampling rate from a channel, assume that they
# are all identical.  Force to be an integer. 
Fs = int(st1[0].stats.sampling_rate)

# Get the length of the trace. Again assume that they are all the same
length = len(st1[0].data)

# Make an array of zeros that is longer than the data by the length of 
# the trail.  This means that the first 10 seconds will have some artifical zeros.
data = np.zeros([length+trail-1,4])

# Make a time column

seismo.time = (1/seismo.Fs) * np.ones(length+trail-1)
seismo.time = np.cumsum(seismo.time) - ((1./seismo.Fs) * trail) 

seismo.ch1 = st1[0].data
seismo.ch2 = st2[0].data
seismo.ch3 = st3[0].data
seismo.CleanData()


# Calculate y-axis Offsetsto plot the traces on a single plot
offset1 = max(seismo.ch1) + abs(min(seismo.ch2))*1.1
offset2 = max(seismo.ch2) + abs(min(seismo.ch3))*1.1
offset  = max(offset1,offset2)

# Setup figure and axes


# Determine the limits for the 3D plot
# This is done by taking the max amplitude with some fudge factor
# so that all axis have the same scale
x_lims = max(abs(min(seismo.ch1))*0.8,max(seismo.ch1)*1.2)
y_lims = max(abs(min(seismo.ch2))*0.8,max(seismo.ch2)*1.2)
z_lims = max(abs(min(seismo.ch3))*0.8,max(seismo.ch3)*1.2)
ax_lims = max([x_lims,y_lims,z_lims])




##################
##################
# PLOT
##################
##################

#
# Setup Axes
#
fig = plt.figure(figsize=(9,9))
ax1 = plt.subplot(2, 2, 1,projection='3d') # 3D Motion Plot
ax2 = plt.subplot(2, 1, 2)                 # Seismogram Plot
ax3 = plt.subplot(2, 2, 2)                 # Text Area Plot

ax3.axis('off') # Turn off the border for text area
figtitle_text = plt.suptitle('',fontsize=labelsize+5) # Give the plot a title

# Add website tagline
ax2.text(0,-0.2,'www.johnrleeman.com',fontsize=labelsize-4,transform = ax2.transAxes)

#
# Plot points in 3D scatter
#
s3d, = ax1.plot([], [], [], marker='o', linestyle='None')

#
# Plot projections of points on axes
#
#s = ax1.scatter(-1*ax_lims*np.ones([trail]),scatter[:,1],scatter[:,2],c='k',s=5)
#s = ax1.scatter(scatter[:,0],ax_lims*np.ones([trail]),scatter[:,2],c='k',s=5)
#s = ax1.scatter(scatter[:,0],scatter[:,1],-1*ax_lims*np.ones([trail]),c='k',s=5)

#
# Plot marker bars to current position
#
#ax1.plot([-1*ax_lims,scatter[-1,0]],[scatter[-1,1],scatter[-1,1]],[scatter[-1,2],scatter[-1,2]],color='k')

#ax1.plot([scatter[-1,0],scatter[-1,0]],[scatter[-1,1],scatter[-1,1]],[-1*ax_lims,scatter[-1,2]],color='k')

#ax1.plot([scatter[-1,0],scatter[-1,0]],[ax_lims,scatter[-1,1]],[scatter[-1,2],scatter[-1,2]],color='k')

#ax1.tick_params(axis='both', which='major', labelsize=ticksize)
#ax1.set_xlabel('West - East',fontsize=labelsize)
#ax1.set_ylabel('South - North',fontsize=labelsize)
#ax1.set_zlabel('Down - Up',fontsize=labelsize)

ax1.set_xlim3d(-1*ax_lims,ax_lims)
ax1.set_ylim3d(-1*ax_lims,ax_lims)
ax1.set_zlim3d(-1*ax_lims,ax_lims)

#
# Plot the seismograms on the 2D plot
#

ax2.plot(seismo.time[9::40],seismo.ch1[::40],color='k')
ax2.plot(seismo.time[9::40],seismo.ch2[::40]+offset,color='k')
ax2.plot(seismo.time[9::40],seismo.ch3[::40]+2*offset,color='k')

# Make and label the scale bar
#ax2.plot([0.03*max(data[:,0]),0.03*max(data[:,0])],[min(data[:,1]),min(data[:,1])+order_of_mag],color='k',linewidth=2)
#ax2.text(0.04*max(data[:,0]),min(data[:,1])+(0.35*order_of_mag),'%.2f mm'%order_of_mag,fontsize=labelsize-4)

#ax2.text(3000,0.2*order_of_mag,'North - South',fontsize=labelsize-2)
#ax2.text(3000,0.2*order_of_mag+offset,'East - West',fontsize=labelsize-2)
#a#x2.text(3000,0.2*order_of_mag+2*offset,'Up - Down',fontsize=labelsize-2)
#a#x2.set_xlabel(r'Seconds Since Earthquake',fontsize=labelsize)
#a#x2.set_ylabel(r'Motion',fontsize=labelsize)
#a#x2.tick_params(axis='both', which='major', labelsize=ticksize)
#a#x2.axes.get_yaxis().set_visible(False)

# Set limits
ax2.set_xlim(min(seismo.time),max(seismo.time))
ax2.set_ylim(min(seismo.ch1)*1.3,max(seismo.ch3)+2.3*offset)

# Time position marker
marker_line = ax2.axvline(x=0,color='r',linewidth=2)



# Place phase markers and color them 
#yloc = max(data[:,3])+2.05*offset
#MarkPhase(ax2,'P',i,travel_times,yloc)
#MarkPhase(ax2,'S',i,travel_times,yloc)


#### Make the text plot
#t_offset = 0.08
#s = 0.9
#stn_text    = ax3.text(0,s,'Station: %s' %stationcd)
#ch1val_text = ax3.text(0,s-t_offset,  'North-South Disp:  %6.2f' %scatter[-1,0])
#ch2val_text = ax3.text(0,s-2*t_offset,'East-West Disp:    %6.2f' %scatter[-1,1])
#ch3val_text = ax3.text(0,s-3*t_offset,'Up-Down Disp:      %6.2f' %scatter[-1,2])
#data = data[::40]

inds = np.arange(10,15)


#TEMP###
ax2.set_xlim(0,3600)


def step(ind):
    cur_time = data[ind,0]
    #print ind,cur_time
    st,x,y,z = seismo.get_scatterdata(ind)
    s3d.set_data(x, y)
    s3d.set_3d_properties(z)
    marker_line.set_xdata(st[-1])
    #print "Moving: ", st[-1]
    return marker_line,s3d


anim = FuncAnimation(fig, step, frames=inds, interval=50, repeat_delay=2000,blit=True)
anim.save('test_seismo_all_dec.mp4')
######### OLD
    
# Pull out the data that will plot in this frame
#scatter = np.zeros([trail,3])
#scatter[:,0] = data[i*Fs:(i+trail)*Fs:Fs,1]
#scatter[:,1] = data[i*Fs:(i+trail)*Fs:Fs,2]
#scatter[:,2] = data[i*Fs:(i+trail)*Fs:Fs,3]
    
   
    
