import pandas as pd
import sys
import os

#ifile = sys.args[1]
#elev = sys.args[2]

# estimating longwave from empirical approaches developed in ~/matlab/LongWaveRadiation


scrdir  = '/home/disk/margaret/mauger/2020_12_SnohoCounty_Flooding/scripts';

indir = ['/home/disk/rocinante/DATA/temp/kcp3/data/interp_pnnl/test/bogachiel/' gcm];
outdir= ['/home/disk/rocinante/DATA/temp/kcp3/data/interp_pnnl/test/lwfix/' gcm];


def datenum(d):    
    return 366 + d.toordinal() + (d - dt.fromordinal(d.toordinal())).total_seconds()/(24*60*60)


# get LW (Dilley and O'Brien + Unsworth and Monteith (1975))
# (needs to be done with the temp correction applied)    
def fix_lw(ifile, elev):
    lat = os.path.basename(ifile).split('_')[1]
    lon = os.path.basename(ifile).split('_')[2]
    df = pd.read_csv(ifile, header=None, sep='[ \t]+', engine='python')
    dates = pd.to_datetime(df[0], format='%m/%d/%Y-%H:%M:%S')
    df.columns = ['date', 'TEMP', 'WIND', 'RH', 'ISW', 'GLW', 'PREC']

    dnm = dates.apply(lambda x: datenum(x))
    
    lw = get_LW(lat, lon, elev, start, end, [], df.copy(), 13, 7)
    df.GLW = lw
    

fix_lw("data_47.81651_-124.35223", 185.30138)

    % round up:
    i60=find(dvc(:,5)>55);
    dvc(i60,4) = dvc(i60,4) + 1;
    dvc(i60,5:6) = 0;
    clear i60,

    % round down:
    dvc(find(dvc(:,5)>0 & dvc(:,5)<5),5) = 0;
    dvc(find(dvc(:,6)>0),6) = 0;

    stdvc = datevec(dtstr{1});
    fidvc = datevec(dtstr{end});
    fidvc(:,5) = 1;	# to ensure final hour is counted. If zero it's sometimes dropped by time_builder

    % ------------------------------------------------
    [LWdob, jnk] = get_LW(frclat(fi),frclon(fi),demelv(fi),stdvc,fidvc,[dvc dnm],prcp,temp,rhum,isw,13,7);
        
   
def sind(v):
    return np.sin(np.deg2rad(v))

def get_LW(lat, lon, elev, dates, temp, rh, isw, prec, LWclr_i, LWall_i):
    # time array

    # Pressure. Assume constant, based on elevation [kPa]
    prsr = (101.325*((1-(.0000225577*(3.281*elev)))**5.25588)) * np.ones(len(temp))

    # Convert precip from mm/hr to cm/hr
    df['PREC'] /= 10


    M_inputs -> time, elev, temp, rh, cloudfactor


    # Calculate average elevation angle [degrees above horizon]
    ela = avg_el(dates, lat, lon, 8, 'MID')

    # Generate a cloud fraction
    cld = cloudfactor_Jessica(dates, lat, df.SWDOWN, 1)

    #  Optical air mass at 101.3 kPa
    sind_ela = sind(ela)
    m = 35 * sind_ela * (1224 * (sind_ela**2 + 1)) ** -0.5

    # Empirical expression for the product of the first two transmission coeffs
    tau_r_tau_pg = 1.021 - (0.083 * np.sqrt( m * (0.00949 * pressure + 0.051 )))

    # Transmission coefficient for absorption by water vapor
    tau_w = 1 - (0.077 * (df.PREC * m) ** 0.3)

    # Transmission coefficient for absorption by aerosols
    tau_a = 0.945 ** m

    # As defined in Appendix A of Flechinger
    s_clr = 1360 * sind_ela * tau_r_tau_pg * tau_w * tau_a
    s = df.ISW / s_clr

    # Atmos. Transmissitivity: a common index of cloud-cover from Sincart et al. 2006
    tau = df.ISW / 1360

    

% Empirical expression for the product of the first two transmission coeffs (Rayleigh and permanent gases)
%whos m pressure,
tau_r_tau_pg = 1.021-(0.084.*sqrt(m.*(0.00949.*pressure+0.051)));

% transmission coeffecient for absorption by water vapor
tau_w = 1-(0.077.*((prcp.*m).^0.3));

% transmission coeffecient for absorption by aerosols
tau_a = 0.945.^m;

% as defined in Appendix A of Flechinger
s_clr = 1360.*sind(El).*tau_r_tau_pg.*tau_w.*tau_a

    
    function [LWhly, LWdly] = get_LW(lat,lon,elv,stdvc,fidvc,time_in,prcp,temp,rhum,isw,LWclr_i,LWall_i)
% PREP INPUTS:
path(path,'/home/disk/margaret/mauger/matlab/LongWaveRadiation/lwscript/lwscript'),

M_PARAMS.STA_Elev = elv;

% time arrays:
%M_INPUTS.TIME = time_builder(stdvc(:,1),stdvc(:,2),stdvc(:,3),stdvc(:,4),stdvc(:,5),...
%			     fidvc(:,1),fidvc(:,2),fidvc(:,3),fidvc(:,4),fidvc(:,5),1);

M_INPUTS.TIME = time_in;	% SWITCH TO USING TIME AS AN INPUT VAR

if size(M_INPUTS.TIME,1) ~= length(temp)
  whos, error('TIME and data time series not of same length'),
end
%M_INPUTS, whos temp prcp isw,

ind = TimeID(M_INPUTS.TIME,'DAY');                                   % Daily indices
ind(:,2) = [];

M_INPUTS.Ta = temp;
M_INPUTS.RH = rhum;
disp(['    ' datestr(now) ' -- done w MINPUTS (met)']),

% --------------------------------------------
% RADIATION:

% Average elevation angle [degrees above horizon]
%M_INPUTS,
El=AVG_EL(M_INPUTS.TIME,lat,lon,8,'BEG');
%whos El,
disp(['    ' datestr(now) ' -- done w AVG_EL']),

% generate a cloud fraction (0-1)
M_INPUTS.c = cloudfactor_Jessica(M_INPUTS.TIME,lat,isw,1);
disp(['    ' datestr(now) ' -- done w cloudfactor_Jessica']),

% Optical [depth?] air mass at 101.3 kPa
m = 35.*sind(El).*((1224.*(sind(El).*sind(El))+1).^-0.5);

% Empirical expression for the product of the first two transmission coeffs (Rayleigh and permanent gases)
%whos m pressure,
tau_r_tau_pg = 1.021-(0.084.*sqrt(m.*(0.00949.*pressure+0.051)));

% transmission coeffecient for absorption by water vapor
tau_w = 1-(0.077.*((prcp.*m).^0.3));

% transmission coeffecient for absorption by aerosols 
tau_a = 0.945.^m;

% as defined in Appendix A of Flechinger
s_clr = 1360.*sind(El).*tau_r_tau_pg.*tau_w.*tau_a;

M_INPUTS.s = isw./s_clr;
M_INPUTS.tau = isw/1360; % Atmos. Transmissivity: a common index of cloud-cover from Sicart et al. 2006
disp(['    ' datestr(now) ' -- done w MINPUTS (radiation)']),

% ----------------------------------------------------------------------------------------
% calc LW:

% initialize LW arrays:
LWhly = NaN(size(M_INPUTS.TIME,1),length(LWclr_i),length(LWall_i));
LWdly  = NaN(length(ind),length(LWclr_i),length(LWall_i));

M_PARAMS.P1 = [];
M_OPTIONS.Method_vpr = 2;

for iclr = 1:length(LWclr_i)
for iall = 1:length(LWall_i)
  M_OPTIONS.Method_clr = LWclr_i(iclr);
  M_OPTIONS.Method_all = LWall_i(iall);

  [jnk,e_all,e_clr] = inc_LW_e_all(M_INPUTS, M_PARAMS, M_OPTIONS);
  LWhly(:,iclr,iall) = jnk;
  jnkdly	      = TimeAverage(M_INPUTS.TIME,LWhly,'day');
  LWdly(:,iclr,iall) = jnkdly(ind);
  clear jnk jnkdly e_all e_clr,

  disp(['    ' datestr(now) ' -- ' ...
	num2str([iclr length(LWclr_i) iall length(LWall_i)],'%d of %d and %d of %d')]),
end
end

disp(['    ' datestr(now) ' -- DONE.']),
