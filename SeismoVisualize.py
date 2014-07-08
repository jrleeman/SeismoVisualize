import numpy as np
import matplotlib.pylab as plt
from obspy import read
from scipy.integrate import cumtrapz
from scipy.signal import detrend
from scipy import signal
import mpl_toolkits.mplot3d.axes3d as p3
import os
from obspy.taup import getTravelTimes
from math import degrees,radians,cos,sin,atan2,sqrt,floor
from obspy.fdsn import Client
from obspy import UTCDateTime
import sys

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

def ProcessTrace(data,Fs):

    # Detrend, remove the mean, and integrate
    # Multiply by 1000 to make units mm
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
    # Calculate Travel Times from EQ to Station
    dist = DegreesDistance(station[0],station[1],earthquake[0],earthquake[1])
    tt = getTravelTimes(dist,earthquake[2])
    travel_times={}
    for item in tt:
        travel_times[item['phase_name']] = item['time']
    return travel_times
    
###############################################################################
# Change these parameters to modify the plotting behavior

network = 'IU'
stationcd = 'ANMO'
location = '10'
duration = 60
time_str = '2014-07-07T11:23:58'
time = UTCDateTime(time_str)

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
    
except:
    print "\n\nINVALID CALL!"
    print "python SeismoVisualize.py NETWORK STATION LOCATION CHN CHE CHZ TIME DURATION\n"
    print "python SeismoVisualize.py IU ANMO 10 BH1 BH2 BHZ 2014-07-07T11:23:58 60\n\n"
    sys.exit()
    



station = GetStationLocation(time,network,stationcd,location,duration)
earthquake = [14.782,-92.371,92.]

#stN = GetData(time,network,stationcd,location,'BH1',duration)
#stE = GetData(time,network,stationcd,location,'BH2',duration)
#stZ = GetData(time,network,stationcd,location,'BHZ',duration)

stN = read('example_data/2014-07-07T11-23-58.IU.ANMO.10.BH1.sac')
stE = read('example_data/2014-07-07T11-23-58.IU.ANMO.10.BH2.sac')
stZ = read('example_data/2014-07-07T11-23-58.IU.ANMO.10.BHZ.sac')

#stN = read('example_data/2014-06-23T20-53-09.IU.SSPA.00.BH1.sac')
#stE = read('example_data/2014-06-23T20-53-09.IU.SSPA.00.BH2.sac')
#stZ = read('example_data/2014-06-23T20-53-09.IU.SSPA.00.BHZ.sac')



# Calculate travel times
travel_times = GetTravelTimes(station,earthquake)

# Settings for Visualization
# Lag

Fs = int(stN[0].stats.sampling_rate)


length = len(stN[0].data)
data = np.zeros([length+trail,4])

# Make a time column
data[:,0] = (1/40.)
data[:,0] = np.cumsum(data[:,0])

for i,trace in zip(range(1,4),[stN,stE,stZ]):
    data[trail:,i] = ProcessTrace(trace[0].data,Fs)

# Calculate y-axis Offsetsto plot the traces on a single plot
offset1 = max(data[:,1]) + abs(min(data[:,2]))*1.1
offset2 = max(data[:,2]) + abs(min(data[:,3]))*1.1
offset  = max(offset1,offset2)

order_of_mag = 10**floor(np.log10(offset)) #mm #### MAKE THIS AUTOMATIC

# Setup figure and axes
fig = plt.figure(figsize=(9,9))

# Determine the limits for the 3D plot
# This is done by taking the max amplitude with some fudge factor
# so that all axis have the same scale
x_lims = max(abs(min(data[:,1]))*0.8,max(data[:,1])*1.2)
y_lims = max(abs(min(data[:,2]))*0.8,max(data[:,2])*1.2)
z_lims = max(abs(min(data[:,3]))*0.8,max(data[:,3])*1.2)
ax_lims = max([x_lims,y_lims,z_lims])


### Fix this to use all of the data
for i in range(0,int(length/Fs-trail)):
    
    print "Time: %f sec" %i
    
    # Pull out the data that will plot in this frame
    scatter = np.zeros([trail,3])
    scatter[:,0] = data[i*Fs:(i+trail)*Fs:Fs,1]
    scatter[:,1] = data[i*Fs:(i+trail)*Fs:Fs,2]
    scatter[:,2] = data[i*Fs:(i+trail)*Fs:Fs,3]
    
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
    ax2.plot(data[:,0],data[:,1],color='k')
    ax2.plot(data[:,0],data[:,2]+offset,color='k')
    ax2.plot(data[:,0],data[:,3]+2*offset,color='k')
    
    # Make and label the scale bar
    ax2.plot([0.03*max(data[:,0]),0.03*max(data[:,0])],[min(data[:,1]),min(data[:,1])+order_of_mag],color='k',linewidth=2)
    ax2.text(0.04*max(data[:,0]),min(data[:,1])+(0.35*order_of_mag),'%.2f mm'%order_of_mag,fontsize=labelsize-4)
    
    ax2.text(3000,0.2*order_of_mag,'North - South',fontsize=labelsize-2)
    ax2.text(3000,0.2*order_of_mag+offset,'East - West',fontsize=labelsize-2)
    ax2.text(3000,0.2*order_of_mag+2*offset,'Up - Down',fontsize=labelsize-2)
    ax2.set_xlabel(r'Seconds Since Earthquake',fontsize=labelsize)
    ax2.set_ylabel(r'Motion',fontsize=labelsize)
    ax2.tick_params(axis='both', which='major', labelsize=ticksize)
    ax2.axes.get_yaxis().set_visible(False)
    
    # Set limits
    ax2.set_xlim(min(data[:,0]),max(data[:,0]))
    ax2.set_ylim(min(data[:,1])*1.3,max(data[:,3])+2.3*offset)
    
    # Time position marker
    ax2.axvline(x=i,color='r')
    
    # Add website tagline
    ax2.text(0,-0.2,'www.johnrleeman.com',fontsize=labelsize-4,transform = ax2.transAxes)
    
    # Place phase markers and color them 
    yloc = max(data[:,3])+2.05*offset
    MarkPhase(ax2,'P',i,travel_times,yloc)
    MarkPhase(ax2,'S',i,travel_times,yloc)
    
    plt.suptitle('%d Seconds After Earthquake'%data[i*Fs,0],fontsize=labelsize+5)
    plt.savefig('frame%06d.png'%i)
    plt.clf()
    
os.system('ffmpeg -r 25 -i frame%06d.png %s_%s_%s.mp4' %(network,station,time_str))
# Option to cleanup
#os.system('rm frame*.png')