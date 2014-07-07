import numpy as np
import matplotlib.pylab as plt
from obspy import read
from scipy.integrate import cumtrapz
from scipy.signal import detrend
from scipy import signal
import mpl_toolkits.mplot3d.axes3d as p3
import os
from obspy.taup import getTravelTimes
from math import degrees,radians,cos,sin,atan2,sqrt
from obspy.fdsn import Client
from obspy import UTCDateTime

def GetData(t0,net,stn,loc,ch,duration):
    client = Client("IRIS")
    st = client.get_waveforms(net, stn, loc, ch, t0, t0+duration*60)
    return st
    
def GetStationLocation(t0,net,stn,loc,duration):    
    client = Client("IRIS")
    stn = client.get_stations(starttime=t0,endtime=t0+duration*60,network=net,station=stn,level='station')
    slat = stn[0][0].latitude
    slon = stn[0][0].longitude
    selev = stn[0][0].elevation
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

### LOCATION DATA

network = 'IU'
stationcd = 'ANMO'
location = '10'
duration = 60
time = UTCDateTime('2014-07-07T11:23:58')

station = GetStationLocation(time,network,stationcd,location,duration)
earthquake = [14.782,-92.371,92.]

stN = GetData(time,network,stationcd,location,'BH1',duration)
stE = GetData(time,network,stationcd,location,'BH2',duration)
stZ = GetData(time,network,stationcd,location,'BHZ',duration)

dist = DegreesDistance(station[0],station[1],earthquake[0],earthquake[1])

# Calculate Travel Times from EQ to Station
travel_times = getTravelTimes(dist,earthquake[2])

for item in travel_times:
    if item['phase_name'] == 'P':
        ptt = item['time']
    if item['phase_name'] == 'S':
        stt = item['time']
        

# Settings for Visualization
# Lag
trail = 10
Fs = stN[0].stats.sampling_rate
labelsize = 14
ticksize  = 12
order_of_mag = 0.1 #mm #### MAKE THIS AUTOMATIC

#stN = read('example_data/2014-07-07T11-23-58.IU.ANMO.10.BH1.sac')
#stE = read('example_data/2014-07-07T11-23-58.IU.ANMO.10.BH2.sac')
#stZ = read('example_data/2014-07-07T11-23-58.IU.ANMO.10.BHZ.sac')

length = len(stN[0].data)
data = np.zeros([length+trail,7])

# Make a time column
data[:,0] = (1/40.)
data[:,0] = np.cumsum(data[:,0])

# Detrend each trace
data[trail:,1] = detrend(stN[0].data)
data[trail:,2] = detrend(stE[0].data)
data[trail:,3] = detrend(stZ[0].data)

# Remove the mean of each trace
data[:,1] = data[:,1] - np.mean(data[:,1])
data[:,2] = data[:,2] - np.mean(data[:,2])
data[:,3] = data[:,3] - np.mean(data[:,3])

# Integrate each trace to obtain ground displacement in mm
data[:,4] = cumtrapz(data[:,1],dx=(1./Fs),initial=0) * 1000.
data[:,5] = cumtrapz(data[:,2],dx=(1./Fs),initial=0) * 1000.
data[:,6] = cumtrapz(data[:,3],dx=(1./Fs),initial=0) * 1000.



# Calculate y-axis Offsetsto plot the traces on a single plot
offset1 = max(data[:,4]) + abs(min(data[:,5]))*1.1
offset2 = max(data[:,5]) + abs(min(data[:,6]))*1.1
offset  = max(offset1,offset2)

# Setup figure and axes
fig = plt.figure(figsize=(9,9))

# Determine the limits for the 3D plot
# This is done by taking the max amplitude with some fudge factor
# so that all axis have the same scale
x_lims = max(abs(min(data[:,4]))*0.8,max(data[:,4])*1.2)
y_lims = max(abs(min(data[:,5]))*0.8,max(data[:,5])*1.2)
z_lims = max(abs(min(data[:,6]))*0.8,max(data[:,6])*1.2)
ax_lims = max([x_lims,y_lims,z_lims])


### Fix this to use all of the data
for i in range(0,int(length/Fs-trail)):
    
    print "Time: %f sec" %i
    
    # Pull out the data that will plot in this frame
    scatter = np.zeros([trail,3])
    scatter[:,0] = data[i*Fs:(i+trail)*Fs:Fs,4]
    scatter[:,1] = data[i*Fs:(i+trail)*Fs:Fs,5]
    scatter[:,2] = data[i*Fs:(i+trail)*Fs:Fs,6]
    
    # Set axes on plot
    ax1 = plt.subplot(221,projection='3d')
    ax2 = plt.subplot(212)
    
    # Make the 3D motion plot
    ax1.tick_params(axis='both', which='major', labelsize=ticksize)
    s = ax1.scatter(scatter[:,0],scatter[:,1],scatter[:,2],c=np.linspace(0,1,10))
    s.set_edgecolors = s.set_facecolors = lambda *args:None
    
    s = ax1.scatter(-1*ax_lims*np.ones([trail]),scatter[:,1],scatter[:,2],c='k',s=5)
    s.set_edgecolors = s.set_facecolors = lambda *args:None
    
    s = ax1.scatter(scatter[:,0],ax_lims*np.ones([trail]),scatter[:,2],c='k',s=5)
    s.set_edgecolors = s.set_facecolors = lambda *args:None
    
    s = ax1.scatter(scatter[:,0],scatter[:,1],-1*ax_lims*np.ones([trail]),c='k',s=5)
    s.set_edgecolors = s.set_facecolors = lambda *args:None
    
    ax1.plot([-1*ax_lims,scatter[-1,0]],[scatter[-1,1],scatter[-1,1]],[scatter[-1,2],scatter[-1,2]],color='k')
    
    ax1.plot([scatter[-1,0],scatter[-1,0]],[scatter[-1,1],scatter[-1,1]],[-1*ax_lims,scatter[-1,2]],color='k')
    
    ax1.plot([scatter[-1,0],scatter[-1,0]],[ax_lims,scatter[-1,1]],[scatter[-1,2],scatter[-1,2]],color='k')
    
    #ax1.set_xlabel('West - East',fontsize=labelsize)
    #ax1.set_ylabel('South - North',fontsize=labelsize)
    #ax1.set_zlabel('Down - Up',fontsize=labelsize)
    
    ax1.set_xlim3d(-1*ax_lims,ax_lims)
    ax1.set_ylim3d(-1*ax_lims,ax_lims)
    ax1.set_zlim3d(-1*ax_lims,ax_lims)
    
    # Make the 2D seismogram plot
    ax2.plot(data[:,0],data[:,4],color='k')
    ax2.plot(data[:,0],data[:,5]+offset,color='k')
    ax2.plot(data[:,0],data[:,6]+2*offset,color='k')
    
    # Make and label the scale bar
    ax2.plot([0.03*max(data[:,0]),0.03*max(data[:,0])],[min(data[:,4]),min(data[:,4])+order_of_mag],color='k',linewidth=2)
    ax2.text(0.04*max(data[:,0]),min(data[:,4])+(0.35*order_of_mag),'%.2f mm'%order_of_mag,fontsize=labelsize-4)
    
    ax2.text(3000,0.2*order_of_mag,'North - South',fontsize=labelsize-2)
    ax2.text(3000,0.2*order_of_mag+offset,'East - West',fontsize=labelsize-2)
    ax2.text(3000,0.2*order_of_mag+2*offset,'Up - Down',fontsize=labelsize-2)
    ax2.set_xlabel(r'Seconds Since Earthquake',fontsize=labelsize)
    ax2.set_ylabel(r'Motion',fontsize=labelsize)
    ax2.tick_params(axis='both', which='major', labelsize=ticksize)
    ax2.axes.get_yaxis().set_visible(False)
    
    # Set limits
    ax2.set_xlim(min(data[:,0]),max(data[:,0]))
    ax2.set_ylim(min(data[:,4])*1.3,max(data[:,6])+2.3*offset)
    
    # Time position marker
    ax2.axvline(x=i,color='r')
    
    # Add website tagline
    ax2.text(0,-0.2,'www.johnrleeman.com',fontsize=labelsize-4,transform = ax2.transAxes)
    
    # Place phase markers and color them (MAKE FUNCTION)
    if i>ptt:
        ax2.text(ptt,max(data[:,6])+2*offset,'P',fontsize=20)
    else:
        ax2.text(ptt,max(data[:,6])+2*offset,'P',fontsize=20,alpha=0.3)
        
    if i>stt:
        ax2.text(stt,max(data[:,6])+2*offset,'S',fontsize=20)
    else:
        ax2.text(stt,max(data[:,6])+2*offset,'S',fontsize=20,alpha=0.3)
    
    plt.suptitle('%d Seconds After Earthquake'%data[i*Fs,0],fontsize=labelsize+5)
    plt.savefig('frame%06d.png'%i)
    plt.clf()
    
os.system('ffmpeg -r 25 -i frame%06d.png output_ANMO_2014.mp4')