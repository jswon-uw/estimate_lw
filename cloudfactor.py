

def cloudfactor_Jessica(dates, lat, isw, timestep=1):
    ''' 
    Trying to calculate a cloud factor presuming you have measured SW rad    
    for timestep in hours, Sin in W/m^2
    '''
    #dvc = time_in(:,1:6)
    #t_met = time_in(:,7)
    #H = round(dvc(:,4))

    ## STEP 1 Potential Solar Radiation
    N = greg2julian(Y, M, D, H, M, S);%day of year, julian day

    S0 = 1360                # solar constant W/m^2
    phi = lat * np.pi / 180  # Latitude in radians    
    ds = 23.45 * np.sin(np.pi * (284 + N) / 180) # solar declination, in degrees    
    delta = ds * np.pi / 180 #converted to radians.    
    hs = np.arccos(-np.tan(phi) * np.tan(delta))    
    Q0=86400 * S0 *(hs * np.sin(phi) * np.sin(delta) + \
                    np.cos(phi) * np.cos(delta) *np.sin(hs)) / np.pi
    # the above is in J per m^2
    # Note, this is extraterrestrial radiation and so has not been attenuated
    # at all, real would depend on atmospheric transmissivity as well as cloudiness
        
    ## STEP 2: Convert observed W per m^2 to daily J per m^2 and calculate
    # clearness index
    L = length(isw)
    if rem(L,(24/timestep))>0
        Sin(L+1:L+rem(L,(24/timestep)))=0
        t_met(L+1:L+rem(L,(24/timestep)))=0
        Q0(L+1:L+rem(L,(24/timestep)))=0
    end
    
    S1=reshape(Sin,24/timestep,length(Sin)./(24/timestep))
    S2=sum(S1)*3600*timestep
    t1=reshape(t_met,24/timestep,length(Sin)./(24/timestep))
    Q1=reshape(Q0,24/timestep,length(Sin)./(24/timestep))
    Q2=mean(Q1)

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
    lowlimit = 0.2;
    highlimit = 0.8;
    
    c2=ones(1,length(k)) #first, make everything cloudy
    sunny=find(k>=highlimit)
    c2(sunny)=0
    sunbreak=find(k>=lowlimit & k<highlimit)
    c2(sunbreak)=1-((k(sunbreak)-lowlimit)./(highlimit-lowlimit))
    
    ## STEP 4: Check everything with a plot
    # figure
    # subplot(3,1,1)
    # plot(t_met,Q0/10^6,'r')
    # hold on
    # plot(t1(1,:),S2/10^6,'b')
    # ylabel('SW (MJ/day)')
    # legend('Incident at top of atmosphere','Observed')
    # axis([datenum(1997,12,1) datenum(1998,4,1) -inf inf])
    # datetick('keeplimits')
    # subplot(3,1,2)
    # plot(t1(1,:),S2./Q2,'k')
    # ylabel('clearness index')
    # axis([datenum(1997,12,1) datenum(1998,4,1) -inf inf])
    # datetick('keeplimits')
    # subplot(3,1,3)
    # plot(t1(1,:),c2,'k')
    # ylabel('cloud factor')
    # axis([datenum(1997,12,1) datenum(1998,4,1) -inf inf])
    # datetick('keeplimits')
    
    ## STEP 5: Assume we have the same cloud factor every day, but distribute
    # it according to the original timeseries
    # note that right now only certain timesteps are hardcoded in.
    if timestep==2:
        c1=[c2; c2; c2; c2; c2; c2; c2; c2; c2; c2; c2; c2]
    elif (timestep==1):
        c1=[c2; c2; c2; c2; c2; c2; c2; c2; c2; c2; c2; c2; ...
            c2; c2; c2; c2; c2; c2; c2; c2; c2; c2; c2; c2]
    elif timestep==3:
        c1=[c2; c2; c2; c2; c2; c2; c2; c2]
    elif timestep==0.5:
        c1=[c2; c2; c2; c2; c2; c2; c2; c2; c2; c2; c2; c2; ...
            c2; c2; c2; c2; c2; c2; c2; c2; c2; c2; c2; c2; ...
            c2; c2; c2; c2; c2; c2; c2; c2; c2; c2; c2; c2; ...
            c2; c2; c2; c2; c2; c2; c2; c2; c2; c2; c2; c2]
    else:
        error('only 1/2 hour, 1 hour, 2 hour and 3 hour timesteps are programmed at this time')
        
    c=reshape(c1,length(isw),1);
    c=c(1:L); #throw away the last ones    
    return
