import numpy as np

def sunang(YYYY, DAY, HH, lat, lon, refraction_flag=0):
    '''
    sunang - a converted fortran program from Dr. Wiscombe at NASA.
    Calculates the elevation angle, azimuthal angle, and Earth-Sun distance.
    
    SYNTAX:
       [EL,AZ,SOLDIST,DEC,HA] = sunang(YYYY,DAY,HH,lat,lon)
    
    INPUTS:
	YYYY	= 1xn vector - Year between 1950 and 2050 
				(approximations become invalid outside this range)
	DAY 	= 1xn vector - Day of year. Must be between 1 and 366
	HH 	= 1xn vector - (local hour) + (time zone number) + (Daylight
       	  savings Time correction; -1 or 0) where (local hour) range is 0 to 
       	  24, (time zone number) range is -12 to +12, and (Daylight Time 
       	  correction) is -1 if on Daylight Time (summer half of year), 0
       	  otherwise.
	lat	= 1x1 scalar  - 90 to -90 with + in North
	long	= 1x1 scalar - -180 to 180 with + in East 
    
    OUTPUTS:
	EL	= 1xn vector - Solar elevation angle (0 @ horizon, 90 @ zenith)
	AZ	= 1xn vector - Azimuth angle
	SOLDST	= 1xn vector - Normalized solar distance (D/D_0)
	DEC	= 1xn vector - Declination angle
	HA	= 1xn vector - Solar hour angle (- before noon, + after noon)
    The solar angle/position calculations are based on the sunang program
    from Warren Wiscombe's ftp site:
    ftp://climate1.gsfc.nasa.gov/wiscombe/Solar_Rad/SunAngles/sunae.f
    
    Authors:  Dr. Joe Michalsky (joe@asrc.albany.edu)
              Dr. Lee Harrison (lee@asrc.albany.edu)
              Atmospheric Sciences Research Center
              State University of New York
              Albany, New York
    
    Modified by:  Dr. Warren Wiscombe (wiscombe@climate.gsfc.nasa.gov)
                  NASA Goddard Space Flight Center
                  Code 913
                  Greenbelt, MD 20771
    
    With reference to: 
      Michalsky, J., 1988: The Astronomical Almanac's algorithm for
        approximate solar position (1950-2050), Solar Energy 40,
         227-235 (but the version of this program in the Appendix
         contains errors and should not be used)
    
    The below commands can be used to test the code to make sure it is doing
    what you think it is doing.
    -------------------------------------------------------------------------
    NASA Goddard:
    [EL,AZ,SOLDST,DEC,HA] = sunang(1995,352,9.525+5,39,-76.8);
    [EL;AZ;SOLDST;DEC*180/pi;HA*180/pi]
    These values in sunang.m
     lat = 39.0; long = -76.8; year = 1995; day = 352; zone = 5; lst =
     9.525
    Should give
    solar azimuth = 143.20 deg; elev = 18.08 deg; 
     hour angle = -38.35; declination = -23.38; 
     air mass = 3.19; solar distance = 0.9840
    -------------------------------------------------------------------------
    Atmospheric Sci Res Ctr (Albany, NY):
    [EL,AZ,SOLDST,DEC,HA] = sunang(1996,12,14.6+5,42.7,-73.83);
    [EL;AZ;SOLDST;DEC*180/pi;HA*180/pi]
    These values in sunang.m
     lat = 42.7; long = -73.83; year = 1996; day = 12; zone = 5; lst = 14.5;
    Should give
     solar azimuth = 216.68 deg; elev = 16.79 deg; 
     hour angle = 37.99; declination = -21.66; 
     air mass = 3.43; solar distance = 0.9835
    '''
    
    YYYY = np.array([YYYY]).flatten()
    DAY = np.array([DAY]).flatten()
    HH = np.array([HH]).flatten()
    
    # Initial Checks
    if ((YYYY < 1950) | (YYYY > 2100)).any():
        print("Bad input variable YEAR")
        return
    elif ((DAY < 1) | (DAY > 366)).any():
        print("Bad input variable DAY")
        return
    elif ((HH < -13.0) | (HH > 36)).any():
        print("Bad input variable HOUR")
        return
    elif (lat < -90.0) | (lat > 90):
        print("Bad input variable lat")
        return
    elif (lon < -180) | (lon > 180):
        print("Bad input variable lon")
        return

    RPD = np.pi / 180

    # current Julian date (actually add 2,400,000 
    # for true JD);  LEAP = leap days since 1949;
    # 32916.5 is midnite 0 jan 1949 minus 2.4e6
    DELTA = YYYY - 1949
    LEAP = DELTA / 4
    JD = 32916.5 + (DELTA*365 + LEAP + DAY) + HH / 24

    # last yr of century not leap yr unless divisible
    # by 400 (not executed for the allowed YEAR range,
    # but left in so our successors can adapt this for 
    # the following 100 years)
    JD[(YYYY % 100 == 0) & (YYYY % 400 != 0)] -= 1
        
    # ecliptic coordinates
    # 51545.0 + 2.4e6 = noon 1 jan 2000
    TIME  = JD - 51545.0

    # force mean longitude between 0 and 360 degs
    MNLONG = 280.460 + 0.9856474 * TIME
    MNLONG = MNLONG % 360
    MNLONG[MNLONG < 0] += 360

    # mean anomaly in radians between 0 and 2*pi
    MNANOM = 357.528 + 0.9856003 * TIME
    MNANOM = MNANOM % 360
    MNANOM[MNANOM < 0] += 360

    # Convert to radians
    MNANOM = MNANOM * RPD

    # ecliptic longitude and obliquity 
    # of ecliptic in radians
    ECLONG = MNLONG + 1.915 * np.sin(MNANOM) + 0.020 * np.sin(2.* MNANOM)
    ECLONG = ECLONG % 360
    ECLONG[ECLONG < 0] += 360
    OBLQEC = 23.439 - 0.0000004 * TIME
    ECLONG = ECLONG * RPD
    OBLQEC = OBLQEC * RPD

    # right ascension
    NUM = np.cos(OBLQEC) * np.sin(ECLONG)
    DEN = np.cos(ECLONG)
    RA = np.arctan(NUM/DEN)

    # Force right ascension between 0 and 2*pi
    RA[DEN < 0] += np.pi
    RA[(DEN >= 0) & (NUM < 0)] += 2 * np.pi
 
    # declination
    DEC = np.arcsin(np.sin(OBLQEC) * np.sin(ECLONG))

    # Greenwich mean sidereal time in hours
    GMST = 6.697375 + 0.0657098242 * TIME + HH

    # Hour not changed to sidereal time since 
    # 'time' includes the fractional day
    GMST = GMST % 24
    GMST[GMST < 0] += 24
    
    # local mean sidereal time in radians
    LMST = GMST + lon / 15
    LMST = LMST % 24
    LMST[LMST < 0] += 24    
    LMST = LMST * 15. * RPD

    # hour angle in radians between -pi and pi
    HA  = LMST - RA;

    HA[HA < -np.pi] += 2 * np.pi
    HA[HA > np.pi] -= 2 * np.pi
    

    # solar azimuth and elevation
    EL  = np.arcsin(np.sin(DEC)*np.sin(lat*RPD) + np.cos(DEC)*np.cos(lat*RPD)*np.cos(HA))
    AZ  = np.arcsin(-np.cos(DEC)*np.sin(HA)/np.cos(EL));

    # Put azimuth between 0 and 2*pi radians
    src = np.sin(DEC)-np.sin(EL)*np.sin(lat*RPD) >= 0
    AZ[src & (np.sin(AZ) < 0)] += 2 * np.pi
    AZ[~src] *= -1
    AZ[~src] += np.pi
    
    
    # Convert elevation and azimuth to degrees
    EL = EL / RPD
    AZ = AZ / RPD

    #   ======== Refraction correction for U.S. Standard Atmos. ==========
    #      (assumes elevation in degs) (3.51823=1013.25 mb/288 K)
    if refraction_flag:    
        if EL >= 19.225:
            REFRAC = 0.00452 * 3.51823 / np.tan(EL*RPD)
        elif EL > -0.766 & (EL < 19.225):
            REFRAC = 3.51823 * ( 0.1594 + EL*(0.0196 + 0.00002*EL) ) / \
                ( 1. + EL*(0.505 + 0.0845*EL) )
        elif EL <= -0.766:
            REFRAC = 0
        EL = EL + REFRAC
    
    #  ===================================================================

    # distance to sun in A.U.
    SOLDST = 1.00014 - 0.01671 * np.cos(MNANOM) - 0.00014 * np.cos(2.*MNANOM);

    # checks
    if ((EL < -90.0) | (EL > 90.0)).any():
        print('output argument EL out of range')
        return
    if ((AZ < 0.0) | (AZ > 360.0)).any():
        print('output argument AZ out of range')
        return

    return [EL, AZ, SOLDST, DEC, HA]
