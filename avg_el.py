import pandas as pd
import numpy as np
from datetime import timedelta
from sunang import sunang

def solar_geo(t, lat, lon, zone):
    '''
    Gives the elevation angle, azimuth, normalized Earth-Sun distance,
    and solar hour angle for a specified location. Calls sunang,
    which is based on Whiscombe's algorithm
    SYNTAX:
 	[EL, AZ, SOLDIST, HA] = SolarGeometry(t,lat,long,zone)
    INPUTS:
	t	= datetime array
	lat	= 1x1 value of latitude  (+N, -S)
	long	= 1x1 value of longitude (+E, -W)
	zone	= 1x1 value of time zone from GMTdef solar_geo():   
    '''
    YYYY = t.year
    # Find integar Julian Day    
    JD = np.floor(t.dayofyear)
    # Find decimal hour
    HHMM = t.hour + t.minute / 60
    # Call sunang
    [ela, az, soldist, dec, ha] = sunang(YYYY, JD, HHMM+zone, lat, lon)

    # Force night time
    ela[ela < 0] = 0
    
    return ela

def avg_el(dates, lat, lon, tz, ref):
    ''' 
    Calculates the average cosSZA during the time interval and returns 
    the effective elevation angle
    
    SYNTAX:
    	EL = AVG_EL(TIME,lat,lon,tz,REF)
    
     INPUTS:
    	TIME	= Nx7 matrix - time_builder format time
    	lat		= 1x1 scalar - degrees north
    	lon		= 1x1 scalar - degrees west
    	tz		= 1x1 sclar - # of time zones West of UTC
    	REF		= string - argument describing how the data is referenced to the time stamp
    			The default assumption is that the time stamp is for the beginning of the 
    			averaging interval. This argument must be specified for accuracy if the
    			default value is not true.
    			'END' - averaged data referenced to the end of the interval
    			'MID' - averaged data referenced to the middle of the interval
    			'BEG' - averaged data referenced to the beginning of the
    			interval
    
     OUTPUTS:
    	EL		= Nx1 vector - Average elevation angle [degrees above horizon]
    '''
    rm_lp = len(dates[(dates.dt.month==2) & (dates.dt.day==29)]) == 0    
    dt = abs(np.diff(dates)).mean()/ np.timedelta64(1, 's')

    if ref == "BEG":
        dates -= timedelta(seconds=0)
    elif ref == "MID":
        dates -= timedelta(seconds=dt/2)
    elif ref == "END":
        dates -= timedelta(seconds=dt)
        
    else:
        print("Unrecognzed REF option")
        raise RuntimeError from None
    
    # Instantaneous Elevation Angle
    t_fine = pd.date_range(dates[0], dates.iloc[-1], freq='5min')    
    
    # Calculate solar geometry for fine time step
    ela_fine = solar_geo(t_fine, lat, lon, tz)
    mew_fine = pd.DataFrame(data = np.sin(np.deg2rad(ela_fine)),
                            index = t_fine, columns = ['ela'])
    
    # Average Elevation angle
    cosSZA = mew_fine[mew_fine>0].resample(timedelta(seconds=dt)).mean()
    cosSZA = cosSZA.fillna(0)
    pd.DataFrame(cosSZA).to_csv('cossza.csv', header=False, index=False)

    
    # Convert from avg(cos(SZA)) to elevation angle
    ela = np.rad2deg(np.arcsin(cosSZA))

    if len(ela) != len(dates):
        print("Remove leap:", rm_lp)
        print("Dims do not match", len(ela), len(dates))
        raise RuntimeError from None

    return ela.ela.values
