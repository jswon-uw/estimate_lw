import pandas as pd
import numpy as np
import sys
import os
import argparse
from avg_el import avg_el
from cloudfactor import cloudfactor_Jessica
from inc_LW import calc_inc_LW

def est_LW(lat, lon, elev, dates, df, LWall_i, LWclr_i):
    # Pressure. Assume constant, based on elevation [kPa]
    prsr = (101.325*((1-(.0000225577*(3.281*elev)))**5.25588))

    # Convert precip from mm/hr to cm/hr
    df['PREC'] /= 10

    # Calculate average elevation angle [degrees above horizon]
    ela = avg_el(dates, lat, lon, 8, 'BEG')

    # Generate a cloud fraction
    cld = cloudfactor_Jessica(dates, lat, df.ISW, 1)
    
    #  Optical air mass at 101.3 kPa
    sind_ela = np.sin(np.deg2rad(ela))
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
def get_lw(ifile, elev, outfile):
    lat = os.path.basename(ifile).split('_')[1]
    lon = os.path.basename(ifile).split('_')[2]
    print("Read file: ", ifile)
    df = pd.read_csv(ifile, header=None, sep='[ \t]+', engine='python')    
    dates = pd.to_datetime(df[0], format='%m/%d/%Y-%H:%M:%S')
    df.columns = ['date', 'TEMP', 'WIND', 'RH', 'ISW', 'GLW', 'PREC']
    
    lw = est_LW(float(lat), float(lon), elev, dates, df.copy(), 7, 13)
    df['GLW'] = lw
    df.to_csv(outfile, sep=' ', header=False, index=False)
    
    
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_file')
    parser.add_argument('elev', type=float)
    parser.add_argument('output_file', nargs='?', default=None)
    args = parser.parse_args()
    ifile = args.input_file
    elev = args.elev
    if args.output_file is not None:
        ofile = args.output_file
    else:
        ofile = ifile

    get_lw(ifile, elev, ofile)
    
        
if __name__ == "__main__":
    main()
