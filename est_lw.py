import pandas as pd
import numpy as np
import sys
import os
import argparse
from avg_el import avg_el
from cloudfactor import cloudfactor_Jessica
from inc_LW import calc_inc_LW

#ifile = sys.args[1]
#elev = sys.args[2]
# estimating longwave from empirical approaches developed in ~/matlab/LongWaveRadiation
#scrdir  = '/home/disk/margaret/mauger/2020_12_SnohoCounty_Flooding/scripts';
#indir = ['/home/disk/rocinante/DATA/temp/kcp3/data/interp_pnnl/test/bogachiel/' gcm];
#outdir= ['/home/disk/rocinante/DATA/temp/kcp3/data/interp_pnnl/test/lwfix/' gcm];


class m_inputs():
    def __init__(self, t, rh, elev):
        self.temp = t
        self.rh = rh
        self.elev = elev
        

def sind(v):
    return np.sin(np.deg2rad(v))

def get_LW(lat, lon, elev, dates, df, LWall_i, LWclr_i):
    # Pressure. Assume constant, based on elevation [kPa]
    prsr = (101.325*((1-(.0000225577*(3.281*elev)))**5.25588))
    print(prsr)

    # Convert precip from mm/hr to cm/hr
    df['PREC'] /= 10

    # Calculate average elevation angle [degrees above horizon]
    ela = avg_el(dates, lat, lon, 8, 'BEG')

    # Generate a cloud fraction
    cld = cloudfactor_Jessica(dates, lat, df.ISW, 1)
    
    #  Optical air mass at 101.3 kPa
    sind_ela = sind(ela)
    m = 35 * sind_ela * (1224 * sind_ela**2 + 1) ** -0.5

    # Empirical expression for the product of the first two transmission coeffs
    tau_r_tau_pg = 1.021 - (0.084 * np.sqrt( m * (0.00949 * prsr + 0.051 )))
    
    # Transmission coefficient for absorption by water vapor
    tau_w = 1 - (0.077 * (df.PREC * m) ** 0.3)

    # Transmission coefficient for absorption by aerosols
    tau_a = 0.945 ** m

    # As defined in Appendix A of Flechinger
    s_clr = 1360 * sind_ela * tau_r_tau_pg * tau_w * tau_a
    df['s'] = df.ISW / s_clr

    # Atmos. Transmissitivity: a common index of cloud-cover from Sincart et al. 2006
    df['tau'] = df.ISW / 1360
    df['c'] = cld

    [LWdwn, e_all, e_clr] = calc_inc_LW(df, elev, LWall_i, LWclr_i, 2)
    return LWdwn
    

# get LW (Dilley and O'Brien + Unsworth and Monteith (1975))
# (needs to be done with the temp correction applied)    
def est_lw(ifile, elev):
    lat = os.path.basename(ifile).split('_')[1]
    lon = os.path.basename(ifile).split('_')[2]
    print("Read file: ", ifile)
    df = pd.read_csv(ifile, header=None, sep='[ \t]+', engine='python')    
    dates = pd.to_datetime(df[0], format='%m/%d/%Y-%H:%M:%S')
    df.columns = ['date', 'TEMP', 'WIND', 'RH', 'ISW', 'GLW', 'PREC']

    #dnm = dates.apply(lambda x: datenum(x))
    
    lw = get_LW(float(lat), float(lon), elev, dates, df.copy(), 7, 13)
    print(lw)
    df.GLW = lw

    

    #fix_lw("data_47.81651_-124.35223", 185.30138)
    #% round up:
    #i60=find(dvc(:,5)>55);
    #dvc(i60,4) = dvc(i60,4) + 1;
    #dvc(i60,5:6) = 0;
    #clear i60,

    #% round down:
    #dvc(find(dvc(:,5)>0 & dvc(:,5)<5),5) = 0;
    #dvc(find(dvc(:,6)>0),6) = 0;

    #stdvc = datevec(dtstr{1});
    #fidvc = datevec(dtstr{end});
    #fidvc(:,5) = 1;	# to ensure final hour is counted. If zero it's sometimes dropped by time_builder

    #% ------------------------------------------------
    #[LWdob, jnk] = get_LW(frclat(fi),frclon(fi),demelv(fi),stdvc,fidvc,[dvc dnm],prcp,temp,rhum,isw,13,7);
        
    
def main():
    #get_LW()
    ifile = "/home/disk/rocinante/DATA/temp/kcp3/data/interp_pnnl/test/bogachiel/access1.0_RCP45/data_47.81651_-124.35223"
    elev = 185.30138
    est_lw(ifile, elev)
    
        
if __name__ == "__main__":
    main()


#M_inputs -> time, elev, temp, rh, cloudfactor


#% initialize LW arrays:
#LWhly = NaN(size(M_INPUTS.TIME,1),length(LWclr_i),length(LWall_i));
#LWdly  = NaN(length(ind),length(LWclr_i),length(LWall_i));
#
#M_PARAMS.P1 = [];
#M_OPTIONS.Method_vpr = 2;

#for iclr = 1:length(LWclr_i)
#for iall = 1:length(LWall_i)
#  M_OPTIONS.Method_clr = LWclr_i(iclr);
#  M_OPTIONS.Method_all = LWall_i(iall);
#
#  clear jnk jnkdly e_all e_clr,
#
#  disp(['    ' datestr(now) ' -- ' ...
#	num2str([iclr length(LWclr_i) iall length(LWall_i)],'%d of %d and %d of %d')]),
#end
#end
#
#disp(['    ' datestr(now) ' -- DONE.']),
