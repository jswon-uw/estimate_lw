import numpy as np
import pandas as pd

## Constants
stefan = 5.67 * (10**-8)      # Stefan-Boltzmann constant (J/s/m^2/K^4)
T_C2K = 273.15                # conversion from C to K
E_L2W = 41840 / 86400         # conversion from langleys/day to W/m2

def inc_LW_e_clr(df, elev, m_clr):
    # Prata (1996) approximation for precipitable water
    # NOTE: corrected to be 465 instead of 4650 (error in Flerchinger Table 1)
    w = 465 * (df.eo / df.TEMP) 
    pd.DataFrame(w).to_csv('w.csv', index=False, header=False)
    
    if m_clr == 1:
        # 1 == Angstrom (1918)
        P = [0.83, 0.18, 0.67]
        e_clr = P[0] - P[1] * 10 ** (-1 * P[2] * df.eo)
    elif m_clr == 2:
        # 2 = Brunt (1932)
        P = [0.52, 0.205]
        e_clr = P[0] + P[1] * np.sqrt(df.eo)
    elif m_clr == 3:
        # 3 = Brutsaert (1975)
        P = [1.723, 1/7]
        e_clr = P[0] * (df.eo / df.TEMP) ** (P[1])
    elif m_clr == 7:
        # 7 = Iziomon et al. (2003)
        P = [0.35, 100, 212,     # Low Land sites 
             0.43, 115, 1489]    # Higher elevation sites
        Mxz = (P[3] - P[0]) / (P[5] - P[2])
        Myz = (P[4] - P[1]) / (P[5] - P[2])
        X = Mxz * (elev - P[2]) + P[0]
        Y = Myz * (elev - P[2]) + P[1]
        e_clr = 1 - X * np.exp(-Y * df.eo/df.TEMP)
    elif m_clr == 13:
        # 13 = Dilley and O'Brien (1998)
        P = [59.38, 113.7, 96.96]
        # Note: corrected to be 2.5 instead of 25 (error in Flerchinger Table 1)
        Lclr = P[0] + P[1] * ((df.TEMP / T_C2K) ** 6) + P[2] * np.sqrt(w/2.5)
        print(Lclr)
        e_clr = Lclr / ((stefan) * (df.TEMP ** 4)) # Effective emissivity
        print(df.TEMP)
        print(e_clr)
    else:
        print("Error: Invalid clear-sky option specified")
        raise RuntimeError from None

    return e_clr
        

def inc_LW_e_all(df, elev, m_all, m_clr):
    '''
    Calculates all sky atmospheric emissivity based on a specified method for
    clear sky and a specified method for cloudy skies
    
    RELEASE NOTES
    Written by Mark Raleigh (raleigh@ucar.edu), Oct 2013)
    Version 2.0 Overhauled by Mark Raleigh (Feb 2015) to have structure inputs
    
    SYNTAX
    [LWdwn,e_all,e_clr] = inc_LW_e_all(M_INPUTS, M_PARAMS, M_OPTIONS)
    
    INPUTS
    M_INPUTS = structure with the following variables (must have exact name)
       TIME = time matrix (per time_builder.m format), same number of rows as other inputs
       Ta = air temperature, K or C (will assume C if mean is less than 40 C)
       eo = vapor pressure (kPa) - optional if RH is provided and option Method_vpr>0
       RH = fractional relative humidity
       c = time series of cloud fraction (0-1 range) for Method_all from 1-8
       s = solar index (0-1 range) for Method_all from 9-10
       tau = transmissiivty (0-1 range) for method 16
    
    M_PARAMS = structure with the coefficients used for the method.  If you
    want to use the defaults, set M_PARAMS.P1 = []; or simply exclude the
    variable in the structure.  Variables include:
       STA_Elev = elevation (m)
       P1_clr = 1xn array of parameter values for PARAMETER 1 (method-specific) for clear sky method
       P2_clr = 1xn array of parameter values for PARAMETER 2 (method-specific) for clear sky method
       ...
       Pn_clr = 1xn array of parameter values for PARAMETER n (method-specific) for clear sky method
    
       P1_all = 1xn array of parameter values for PARAMETER 1 (method-specific) for all sky method (cloud correction method)
       P2_all = 1xn array of parameter values for PARAMETER 2 (method-specific) for all sky method (cloud correction method)
       ...
       Pn_all = 1xn array of parameter values for PARAMETER n (method-specific) for all sky method (cloud correction method)
    
    M_OPTIONS = structure with options for the
       Method_vpr = if vapor pressure (eo) is not supplied in M_INPUTS,
       need to specify how to calculate it based on Ta and RH. Methods:
           0 = do not calculate (eo is supplied)
           1 = calculate with Dozier and Marks approach
           2 = calculate with Clausius-Clapeyron e_sat in mb (hPa) from Murray 1967
       Method_clr = enter code for clear-sky longwave method, where:
           1 = Angstrom (1918)
           2 = Brunt (1932)
           3 = Brutsaert (1975)
           4 = Garratt (1992)
           5 = Idso and Jackson (1969) (Idso-1)
           6 = Idso (1981) (Idso-2)
           7 = Iziomon et al. (2003)
           8 = Keding (1989)
           9 = Niemela et al. (2001)
           10 = Prata (1996)
           11 = Satterlund (1979)
           12 = Swinbank (1963)
           13 = Dilley and O'Brien (1998)
    
           %%% additional methods from Juszak and Pellicciotti (2013)
           14 = Maykut and Church (1973)
           15 = Konzelmann et al. (1994)
           16 = Dilley and O'Brien (A) (1998)
    
           %%% other methods
           17 = Campbell and Norman (1998) as cited by Walter et al (2005)
           18 = Long and Turner (2008) - based on Brutsaert (1975)
           19 = Ohmura (1982) as cited by Howard and Stull 2013
           20 = Efimova (1961) as cited by Key et al (1996)
    
       Method_all = enter code for all-sky longwave method, where:
           %%% Cloud cover based methods
           1 = Brutsaert (1982)
           2 = Iziomon et al. (2003)
           3 = Jacobs (1978)
           4 = Keding (1989)
           5 = Maykut and Church (1973)
           6 = Sugita and Brutsaert (1993)
           7 = Unsworth and Monteith (1975)
           8 = Kimball et al. (1982)
           %%% Clearness/solar index based methods
           9 = Crawford and Duchon (1999)
           10 = Lhomme et al. (2007)
    
           %%% additional methods from Juszak and Pellicciotti (2013)
           11 = Bolz (1949)
           12 = Konzelmann et al (1994)
    
           %%% other methods
           13 = Sicart et al (2006) (eqn 9, daily average transmissivity)
           14 = Pirazzini et al (2000) as cited by Gubler (eqn 21)
           15 = Pirazzini et al (2000) as cited by Gubler (eqn 22)
           16 = Gubler et al (2012) eqn 23
           17 = Zillman (1972) as cited by Key et al (1996)
           18 = Duarte et al (2006) eqn 22
           19 = Kruk et al (2010) eqn 18
           20 = TVA (1972) as cited by Bras (1990)
    
    OUTPUTS
    LWdwn = incoming longwave (W m^-2)
    e_all = all-sky emissivity of the atmosphere
    e_clr = clear sky emissivity (prior to all-sky correction)
    
    NOTES
    The list of LW methods is populated based on Flerchinger et al. (2009).
    Additional empirical LW models may be added to this list.
    '''
    
    
    #%%%%%%% at this point, eo should be in kPa and Ta should be in K
    e_clr = inc_LW_e_clr(df, elev, m_clr)
    
    #%% Method-specific default parameters
    if m_all == 1:
        # Brusaert (1982)
        P = [0.22]
        e_all = (1 + P[0] * df.c) * e_clr
    elif m_all == 2:
        # 2 = Iziomon et al. (200)
        P = [0.35, 212, 0.50, 1489, 2]
        Mzz = (P[2] - P[0]) / (P[3] - P[2])
        Z = Mzz * (elev - P[1]) + P[0]
        e_all = (1 + Z * df.c ** P[4]) * e_clr
    elif m_all == 3:
        # 3 = Jacobs (1978)
        P = [0.26]
        e_all = (1 + P[0] * df.c) * e_clr
    elif m_all == 4:
        # 4 = Keding (1989)
        P = [153, 2.183]
        e_all = (1 + P[0] * df.c ** P[1]) * e_clr
    elif m_all == 5:
        # 5 = Maykut and Church (1973)
        P = [0.22, 2.75]
        e_all = (1 + P[0] * df.c ** P[1]) * e_clr
    elif m_all == 6:
        # 6 = Sugita and Brutsaert (1993)
        P = [0.0496, 2.45]
        e_all = (1 + P[0] * df.c ** P[1]) * e_clr        
    elif m_all == 7:
        # 7 = Unsworth and Monteith (1975)
        # Note that this is the same as Campbell and Norman (1998) as cited by Walkter et al (2005)
        P = [0.84, 0.84]
        e_all = (1 - P[0] * df.c) * e_clr + P[1] * df.c
    elif m_all == 8:
        # 8 = Kimball et al. (1982)
        P = [11, -0.6732, 0.6240*10**-2, 0.9140*10**-5,
             0.24, 2.98*10**-6, 3000, 1.4, 0.4]
        
    elif m_all == 9:
        # 9 = Crawford and Duchon (1999)
        P = []
        e_all = (1 - df.s) + df.s * e_clr
    elif m_all == 10:
        # 10 = Lhomme et al. (2007)
        P = [1.37, 0.34]
        e_all = (P[0] - P[1]*df.s) * e_clr
    elif m_all == 11:
        # 11 = Bolz (1949)
        P = [0.22, 2.5]
        e_all = e_clr * (1 + P[0] * df.c ** P[1])        
    elif m_all == 12:
        # 12 = Konzelmann et al (1994)
        P = [4, 0.952]
        e_all = e_clr * (1 - df.c**P[0]) + (P[1] * df.c**P[0])
    elif m_all == 13:
        #    # 13 = Sicart et al (2006) (eqn 9, daily average transmissivity)
        P = [0.44, 0.18]
        #    if len(df.tau) > 1:
        #        tauD = df.tau.resample('D').mean()
        #    else:
        #        tauD = df.tau
        #    e_all = e_clr * (1 + P[0] * df.RH - P[1] * tauD)
        
        ## TODO Handle daily Tau
    elif m_all == 14:
        # 14 = Pirazzini et al (2000) as cited by Gubler (eqn 21)
        P = [0.40, 2.00]
        e_all = e_clr * (1 + P[0] * df.c ** P[1])
    elif m_all == 15:
        # 15 = Pirazzini et al (2000) as cited by Gubler (eqn 22)
        P = [0.979, 6, 4]
        e_all = e_clr * (1 - df.c ** P[1]) + P[0] * df.c ** P[2]
    elif m_all == 16:
        # 16 = Gubler et al (2012) eqn 23
        P = [(0.968+0.985+0.940+0.928+0.987+0.926+0.828)/7,   # e_oc value (average across 7 sites)
             (3.77+2.05+4.08+3.28+2.05+5.02+0.76)/7,          # p1
             (2.97+1.61+2.94+2.57+1.78+3.74+1.24)/7]          # p2

        ## TODO: Handle daily Tau
    elif m_all == 17:
        # 17 = Zillman (1972) as cited by Key et al (1996)
        P = [0.96, 9.2 * 10**-6]
        e_all = e_clr + P[0] * (1 - P[1]) * df.c
    elif m_all == 18:
        # 18 = Duarte et al (2006) eqn 22 - same form as Konzelmann
        P = [0.671, 0.990]    
        e_all = e_clr  * (1 - df.s ** P[0]) + (P[1] * df.s **P[1])
    elif m_all == 19:
        # 19 = Kruk et al (2010) eqn 18 - same form as Bolz, Duarte et al eqn 21, etc
        P = [0.1007, 0.9061]
        e_all = e_clr * (1 + P[0] * df.s ** P[1])
    elif m_all == 20:
        # 20 = TVA (1972) as cited by Bras (1990) - same form as Bolz
        P = [0.17, 2]
        e_all = e_clr * (1 + P[0] * df.c ** P[1])        
    else:
        print("Error: Invalid all-sky option specified")
        raise RuntimeError from None

    
    if m_all != 8:
        e_all[e_all < 0] = 0
        e_all[e_all > 1] = 1
        LWdwn = e_all * stefan * df.TEMP ** 4
    return [LWdwn, e_all, e_clr]    

    
def calc_inc_LW(df, elev, m_all, m_clr, m_vpr):
    ## Common variables
    # Check to see if temp is in celsius, and if so, convert to K
    if df.TEMP.mean() < 40:
        # assume inputs was in celsius
        df.TEMP += T_C2K
    
    # Check RH is fractional
    if df.RH.mean() > 1:
        print("Assuming RH input was in percent. Converting to fractional")
        df.RH /= 100
        

    if m_vpr == 0:
        if 'eo' not in df.columns:
            print("eo must be provided in inputs, or include RH and specify method for computing eo")
            raise RuntimeError from None        

    elif m_vpr == 1:
        print("This has not been coded yet")
        raise RuntimeError from None
        
    elif m_vpr == 2:
        # Clausius-Clapeyron e_sat in mb (hPa) from Murray 1967
        sat_vap_pressure = 6.1078 * np.exp((17.2693882 * (df.TEMP - T_C2K))/ (df.TEMP - 35.86))
        df['eo'] = df.RH * sat_vap_pressure

        # Convert eo from hPa to kPa
        df['eo'] /= 10
        
    else:
        print("Invalid vapor pressure option specified")
        raise RuntimeError from None

    
    return inc_LW_e_all(df, elev, m_all, m_clr)
    


    
'''    
%%% First, calculate clear-sky emissivity for the selected method
e_clr = inc_LW_e_clr(M_INPUTS, M_PARAMS, M_OPTIONS);

%%% double check to make sure we have e_clr within realistic limits [0 1]
e_clr(e_clr<0) = 0;
e_clr(e_clr>1) = 1;

%%% Now, calculate all-sky emissivity
if M_OPTIONS.Method_all==1
    % 1 = Brutsaert (1982)
    e_all = (1+ M_PARAMS.P1_all .*c).*e_clr;
elseif M_OPTIONS.Method_all==2
    % 2 = Iziomon et al. (2003)
    Mzz = (M_PARAMS.P3_all-M_PARAMS.P1_all)/(M_PARAMS.P4_all-M_PARAMS.P2_all);
    Z = Mzz .*(M_PARAMS.STA_Elev - M_PARAMS.P2_all) + M_PARAMS.P1_all;
    e_all = (1+Z.*c.^M_PARAMS.P5_all).*e_clr;
elseif M_OPTIONS.Method_all==3
    % 3 = Jacobs (1978)
    e_all = (1+ M_PARAMS.P1_all .*c).*e_clr;
elseif M_OPTIONS.Method_all==4
    % 4 = Keding (1989)
    e_all = (1+ M_PARAMS.P1_all .*c.^M_PARAMS.P2_all).*e_clr;
elseif M_OPTIONS.Method_all==5
    % 5 = Maykut and Church (1973)
    e_all = (1+ M_PARAMS.P1_all .*c.^M_PARAMS.P2_all).*e_clr;
elseif M_OPTIONS.Method_all==6
    % 6 = Sugita and Brutsaert (1993)
    e_all = (1+ M_PARAMS.P1_all .*c.^M_PARAMS.P2_all).*e_clr;
elseif M_OPTIONS.Method_all==7
    % 7 = Unsworth and Monteith (1975)
    e_all = (1-M_PARAMS.P1_all.*c).*e_clr + M_PARAMS.P2_all.*c;
elseif M_OPTIONS.Method_all==8
    % 8 = Kimball et al. (1982)
    Tc = Ta - M_PARAMS.P1_all; %%% no seasonal variation in Tc yet... could add in later
    f8 = M_PARAMS.P2_all + M_PARAMS.P3_all.*Tc - M_PARAMS.P4_all.* Tc.^2;
    e8z = M_PARAMS.P5_all + M_PARAMS.P6_all .* eo.^2 .* exp(M_PARAMS.P7_all./Ta);
    tau8 = 1- e8z.*(M_PARAMS.P8_all-M_PARAMS.P9_all.*e8z);
    
    Lclr = e_clr.*stefan.*Ta.^4;
    
    Ld = Lclr + tau8.*c.*f8.*stefan.*(Tc.^4);
    
    e_all = Ld./(stefan.*(Tc.^4));  % effective emissivity
    
    %%% make sure e_all is in realistic limits [0 1]
    e_all(e_all<0) = 0;
    e_all(e_all>1) = 1;
    
    %%% recalculate Ld
    LWdwn = e_all.*stefan.*(Tc.^4);        % note it is Tc not Ta for this method
    
%     % 8 = Kimball et al. (1982)
%     nparams = 9;
%     P(1) = 11; P(2) = -0.6732; P(3) = 0.6240 .* (10^-2); P(4) = 0.9140 .* (10^-5);
%     P(5) = 0.24; P(6) = 2.98 .* (10^-6); P(7) = 3000; P(8) = 1.4; P(9) = 0.4;
elseif M_OPTIONS.Method_all==9
    % 9 = Crawford and Duchon (1999)
    e_all = (1-s) + s.*e_clr;
elseif M_OPTIONS.Method_all==10
    % 10 = Lhomme et al. (2007)
    e_all = (M_PARAMS.P1_all - M_PARAMS.P2_all.*s).*e_clr;
elseif M_OPTIONS.Method_all==11
    % 11 = Bolz (1949)
    e_all = e_clr.*(1+M_PARAMS.P1_all.*c.^M_PARAMS.P2_all);
elseif M_OPTIONS.Method_all==12
    % 12 = Konzelmann et al (1994)
    e_all = e_clr .*(1-c.^M_PARAMS.P1_all)+(M_PARAMS.P2_all.*c.^M_PARAMS.P1_all);
elseif M_OPTIONS.Method_all==13
    % 13 = Sicart et al (2006) (eqn 9, daily average transmissivity)
    if numel(tau)>1
        tauD = TimeAverage(M_INPUTS.TIME,tau,'DAY');
    else
        tauD = tau;
    end
    e_all = e_clr.*(1+ M_PARAMS.P1_all.*RH - M_PARAMS.P2_all .* tauD);
elseif M_OPTIONS.Method_all==14
    % 14 = Pirazzini et al (2000) as cited by Gubler (eqn 21)
    e_all = e_clr.*(1+ M_PARAMS.P1_all.*c.^M_PARAMS.P2_all);
elseif M_OPTIONS.Method_all==15
    % 15 = Pirazzini et al (2000) as cited by Gubler (eqn 22)
    e_all = e_clr.*(1-c.^M_PARAMS.P2_all) + M_PARAMS.P1_all.*c.^M_PARAMS.P3_all;
elseif M_OPTIONS.Method_all==16
    % 16 = Gubler et al (2012) eqn 23
    if numel(tau)>1
        tauD = TimeAverage(M_INPUTS.TIME,tau,'DAY');
    else
        tauD = tau;
    end
    e_all= (e_clr.*tauD.^M_PARAMS.P3_all) + M_PARAMS.P1_all.*(1-tauD.^M_PARAMS.P2_all);
elseif M_OPTIONS.Method_all==17
    % 17 = Zillman (1972) as cited by Key et al (1996)
    e_all = e_clr + M_PARAMS.P1_all.*(1-M_PARAMS.P2_all).*c;
elseif M_OPTIONS.Method_all==18
    % 18 = Duarte et al (2006) eqn 22 - same form as Konzelmann
    e_all = e_clr .*(1-s.^M_PARAMS.P1_all)+(M_PARAMS.P2_all.*s.^M_PARAMS.P1_all);
elseif M_OPTIONS.Method_all==19
    % 19 = Kruk et al (2010) eqn 18 - same form as Bolz, Duarte et al eqn 21, etc
    e_all = e_clr.*(1+M_PARAMS.P1_all.*s.^M_PARAMS.P2_all);
elseif M_OPTIONS.Method_all==20
    % 20 = TVA (1972) as cited by Bras (1990) - same form as Bolz
    e_all = e_clr.*(1+M_PARAMS.P1_all.*c.^M_PARAMS.P2_all);
end

if M_OPTIONS.Method_all ~= 8
    %%% make sure e_all is in realistic limits [0 1]
    e_all(e_all<0) = 0;
    e_all(e_all>1) = 1;
    LWdwn = e_all.*stefan.*Ta.^4;
end
'''
