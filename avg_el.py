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
    YYYY = t.dt.year
    # Find integar Julian Day    
    JD = np.floor(t.dt.dayofyear)
    # Find decimal hour
    HHMM = t.dt.hour + t.dt.minute / 60
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
    dt = abs(np.diff(dates)).mean()/ np.timedelta64(1, 'D')
    
    if ref == "BEG":
        dt = -timedelta(dt)
    elif ref == "MID":
        dt = timedelta(dt)/2
    elif ref == "END":
        dt = timedelta(0)
    else:
        print("Unrecognized REF option")
        raise RuntimeError from None
    
    # Instantaneous Elevation Angle
    t_fine =  pd.date_range(dates[0], dates.iloc[-1]-dt, freq='5min')

    # Calculate solar geometry for fine time step
    ela_fine = solar_geo(t_fine, lat, lon, tz)
    mew_fine = np.sin(np.deg2rad(ela_fine))

    # Average Elevation angle
    cosSZA = mew_fine.resample().mean()
    cosSZA = cosSZA.fillna(0)
    
    # Convert from avg(cos(SZA)) to elevation angle
    ela = np.arcsin(np.deg2rad(cosSZa))

    if len(ela) != len(dates):
        print(rm_lp, "Dims do not match")
        raise RuntimeError from None

    return ela
