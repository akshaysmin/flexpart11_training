'''
Translated by Marie D. Mulder
from a fortran script of the VAST project.
'''

import os,sys
import numpy as np

import source_parameters

def make_apriori(ph_file):

    vol_alt = source_parameters.vol_alt
    partdens = source_parameters.partdens
    fineashfraction = source_parameters.fineashfraction
    unc_apriori = source_parameters.unc_apriori
    nhlevels = source_parameters.nhlevels
    hlevel_start = source_parameters.hlevel_start
    hstep = source_parameters.hstep
    scaling_total_ash = source_parameters.scaling_total_ash
    
    startdate, starthour, enddate, endhour, plumeheight,ntime = read_ph(ph_file)
    print('first startdate: '+ str(startdate[0]))
    print('first plume_height: '+ str(plumeheight[0]))

    hlevels_low, hlevels_high =get_lows_highs(hlevel_start, hstep, nhlevels)
    print('first hlevels_low: ' + str(hlevels_low[0]))

    totalmass_scaled, dMdt_ash = calculate_eruption_mass(ntime,plumeheight,vol_alt,\
            partdens,startdate,starthour,enddate,endhour,scaling_total_ash)
    print('totalmass_scaled: '+str(totalmass_scaled))

    maxheight_indx, maxheight = find_height_box(ntime,nhlevels, hlevels_high, plumeheight)

    a_priori_2d, a_priori_2d_unc, total_a_priori = create_apriori_matrix(\
    ntime, nhlevels, dMdt_ash, maxheight_indx, unc_apriori)


    current_dir = os.getcwd()
    with open(current_dir+'a_priori_2d.dat', 'w') as f:
        # Write lower and upper boundaries of height levels
        f.write(" ".join(map(str, hlevels_low)) + "\n")
        f.write(" ".join(map(str, hlevels_high)) + "\n")

        # Write data for each time interval
        for i in range(ntime):
            # Write emissions data
            line1 = (
                f"{enddate[i]} {endhour[i]} {total_a_priori[i]} "
                f"{float(plumeheight[i]) * 1000} "
                + " ".join(map(str, a_priori_2d[i]))
            )
            f.write(line1 + "\n")
            # Write uncertainties data
            line2 = (
                f"{enddate[i]} {endhour[i]} {total_a_priori[i]} "
                f"{float(plumeheight[i]) * 1000} "
                + " ".join(map(str, a_priori_2d_unc[i]))
            )
            f.write(line2 + "\n")
    if os.path.isfile(current_dir+'/a_priori_2d.dat'):
        print('Ready writing: '+current_dir+'/a_priori_2d.dat')
    else:
        print('Writing a_priori_2d.dat failed')


def get_lows_highs(hlevel_start: float, hstep: float, nhlevels: int) -> tuple[list[float], list[float]]:
    # Initialize lists
    hlevels_low = [0.0] * nhlevels
    hlevels_high = [0.0] * nhlevels

    # Calculate first level
    hlevels_low[0] = hlevel_start
    hlevels_high[0] = hlevel_start + hstep

    # Calculate remaining levels
    for i in range(2, nhlevels):
        hlevels_low[i] = hlevels_low[i-1] + hstep
        hlevels_high[i] = hlevels_low[i] + hstep

    return hlevels_low, hlevels_high


def read_ph(ph_file):
    startdate = []
    starthour = []
    enddate = []
    endhour = []
    plumeheight = []

    try:
        with open(ph_file, "r") as f:
            # Read and discard the header line
            f_lines = f.readlines()[1:]
            ntime=0
            for line in f_lines:
                parts = line.split()
                startdate.append(parts[0])
                starthour.append(parts[1])
                enddate.append(parts[2])
                endhour.append(parts[3])
                plumeheight.append(parts[4])
                ntime+=1

    except FileNotFoundError:
        print(f"Error: File '{ph_file}' not found.")
        return None, None, None, None, None

    return startdate, starthour, enddate, endhour, plumeheight, ntime

def calculate_eruption_mass(ntime,plumeheight,vol_alt,partdens,startdate,starthour,enddate,endhour,scaling_total_ash):
    """
    Calculates the mass eruption rate using Mastin's (2009) formula and
    total SO2 mass released over given time intervals.
    """
    totalmass = 0.0
    dMdt = np.zeros(ntime)
    dMdt_ash = np.zeros(ntime)
    massrelease = np.zeros(ntime)

    # Reference one second difference in julian date format
    onesec = juldate(20111010, 1) - juldate(20111010, 0)
    for i in range(ntime):
        # H is plume height in km above vent (plumeheight in km - vol_alt in km)
        # Assuming plumeheight is in km and vol_alt is in meters (hence vol_alt / 1000)
        h_km = float(1*plumeheight[i]) - (vol_alt / 1000.0)

        # Mastin (2009) formula solved for dMdt (kg/s)
        # dMdt = ((H / 2.0) ** (1 / 0.241)) * rho
        dMdt[i] = ((h_km / 2.0) ** (1.0 / 0.241)) * partdens

        # Calculate time interval in seconds
        startd = juldate(startdate[i], starthour[i])
        endd = juldate(enddate[i], endhour[i])
        seconds = (endd - startd) / onesec

        # Calculate mass released during the time interval
        massrelease[i] = dMdt[i] * seconds
        totalmass += massrelease[i]

    totalmass_scaled=0
    for i in range(1,ntime):
        dMdt_ash[i] = dMdt[i] * scaling_total_ash/ totalmass
        startd = juldate(startdate[i], starthour[i])
        endd = juldate(enddate[i], endhour[i])
        onesec = juldate(20111010, 1) - juldate(20111010, 0)
        seconds = (endd - startd) / onesec
        massrelease[i] = dMdt_ash[i] * seconds
        totalmass_scaled = totalmass_scaled + massrelease[i]

    return totalmass_scaled, dMdt_ash

def find_height_box(ntime,nhlevels, hlevels_high, plumeheight):
    maxheight_indx = np.zeros(ntime)
    maxheight = np.zeros(ntime)

    for i in range(ntime):
        for j in range(nhlevels):
            print('i,j : '+str(i)+','+str(j))
            plume_height_m = float(plumeheight[i]) * 1000
            print('plume_height_m: '+str(plume_height_m))
            upper_bound = plume_height_m + 100

            if hlevels_high[j] <= upper_bound:
                maxheight_indx[i] = j
                maxheight[i] = hlevels_high[j]
                break  # Exit inner loop once the first matching box is found

    return maxheight_indx, maxheight

def create_apriori_matrix(
    ntime: int,
    nhlevels: int,
    dMdt_ash: list[float],
    maxheight_indx: list[int],
    unc_apriori: float
) -> tuple[list[list[float]], list[list[float]], list[float]]:

    # Initialize matrices with zeros
    a_priori_2d = [[0.0 for _ in range(nhlevels)] for _ in range(ntime)]
    a_priori_2d_unc = [[0.0 for _ in range(nhlevels)] for _ in range(ntime)]
    total_a_priori = [0.0] * ntime

    for i in range(ntime):
        total_apri = 0.0
        max_idx = maxheight_indx[i]

        for j in range(nhlevels):
            if j <= max_idx - 1:
                # Distribute the mass release rate over the number of height levels up to maximum height
                a_priori_2d[i][j] = dMdt_ash[i] / max_idx
                a_priori_2d_unc[i][j] = unc_apriori
            else:
                a_priori_2d[i][j] = 0.0
                a_priori_2d_unc[i][j] = unc_apriori

            total_apri += a_priori_2d[i][j]

        total_a_priori[i] = total_apri

    return a_priori_2d, a_priori_2d_unc, total_a_priori


def juldate(yyyymmdd, time_hhmmss=None, format_type=1):
    from datetime import datetime, timedelta
    date_str = str(yyyymmdd)
    year, month, day = int(date_str[:4]), int(date_str[4:6]), int(date_str[6:8])

    # Handle time component if provided
    if time_hhmmss:
        time_str = str(time_hhmmss).zfill(6)  # Ensure 6 digits (HHMMSS)
        hour, minute, second = int(time_str[:2]), int(time_str[2:4]), int(time_str[4:6])
    else:
        hour, minute, second = 0, 0, 0

    dt = datetime(year, month, day, hour, minute, second)

    # Calculate Julian Date
    a = (14 - month) // 12
    y = year + 4800 - a
    m = month + 12 * a - 3
    jdn = day + (153 * m + 2) // 5 + 365 * y + y // 4 - y // 100 + y // 400 - 32045
    jd = jdn + (dt.hour - 12) / 24 + dt.minute / 1440 + dt.second / 86400

    return jd if format_type == 1 else (jd - 2400000.5)


#---------------

file_plume_heights = sys.argv[1]

make_apriori(file_plume_heights)




