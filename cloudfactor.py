import numpy as np

def greg2julian(date):
    year = date.dt.year
    month = date.dt.month
    day = date.dt.day
    hour = date.dt.hour
    minute = date.dt.minute
    sec = date.dt.second
    JD = ( 367 * year ) \
        - np.floor ( 7 * ( year + np.floor( ( month + 9 ) / 12 ) ) / 4 ) \
        - np.floor( 3 * ( np.floor( ( year + ( month - 9 ) / 7 ) / 100 ) + 1 ) / 4 ) \
        + np.floor( ( 275 * month ) / 9 ) + day + 1721028.5 \
        + ( hour + ( minute / 60 ) + ( sec / 3600 ) ) / 24
    return JD

def cloudfactor_Jessica(dates, lat, isw, timestep=1):
    ''' 
    Trying to calculate a cloud factor presuming you have measured SW rad    
    for timestep in hours, Incoming shortwave in W/m^2
    '''
    ## STEP 1 Potential Solar Radiation, J/m^2
    # Note, this is extraterrestrial radiation and so has not been attenuated
    # at all, real would depend on atmospheric transmissivity as well as cloudiness
    N = greg2julian(dates)   # day of year, julian day
    S0 = 1360                # solar constant W/m^2
    phi = lat * np.pi / 180  # Latitude in radians    
    ds = 23.45 * np.sin(2 * np.pi * (284 + N) / 365) # solar declination, in degrees    
    delta = ds * np.pi / 180 #converted to radians.
    hs = np.arccos(-np.tan(phi) * np.tan(delta))    
    Q0 = 86400 * S0 *(hs * np.sin(phi) * np.sin(delta) + \
                    np.cos(phi) * np.cos(delta) * np.sin(hs)) / np.pi
            
    ## STEP 2: Convert observed W per m^2 to daily J per m^2 and calculate
    # clearness index
    ist = dates[dates.dt.hour==0].index.min()
    ifi = dates[dates.dt.hour==23].index.max()+1
    xs = int(24/timestep)
    ys = int(len(isw[ist:ifi])/xs)

    S1 = isw[ist:ifi].values.reshape(ys, xs)
    Q1 = Q0[ist:ifi].values.reshape(ys, xs)
    
    S2 = S1.sum(axis=1) * 3600 * timestep
    Q2 = Q1.mean(axis=1)
    
    # so the clearness index, k, is the ratio of the observed solar flux to the
    # potential solar flux at the outside of the atmosphere, Q0!
    k = S2 / Q2  # clearness index, one value per day.
    
    ## STEP 3:  Convert clearness index to cloud factor
    # Campbell [1985] suggests that c can be linearly
    # interpolated between c = 1.0 at a clearness index of 0.4 for
    # complete cloud cover (kcld) to c = 0.0 at a clearness index of
    # 0.7 (kclr). Others have used values of 0.35 and 0.6 for kcld
    # and kclr, respectively [Flerchinger, 2000; Xiao et al., 2006];
    # the optimum values are unknown and may perhaps depend
    # on the location. (quoted from Flerchinger et al. 2009)    
    lowlimit = 0.2
    highlimit = 0.8

    c2 = np.ones(len(k))        
    sunny = k >= highlimit
    c2[sunny] = 0
    sunbreak = (k>= lowlimit) & (k <highlimit)
    c2[sunbreak] = 1 - ((k[sunbreak]- lowlimit) / (highlimit - lowlimit))
    
    ## STEP 4: Assume we have the same cloud factor every day, but distribute
    # it according to the original timeseries
    # note that right now only certain timesteps are hardcoded in.
    if timestep==2:
        c1 = c2.repeat(12)
    elif (timestep==1):
        c1 = c2.repeat(24)
    elif timestep==3:
        c1 = c2.repeat(8)
    elif timestep==0.5:
        c1 = c2.repeat(48)
    else:
        error('only 1/2 hour, 1 hour, 2 hour and 3 hour timesteps are programmed at this time')

    c = np.zeros(len(isw))
    c[ist:ifi] = c1
            
    return c
