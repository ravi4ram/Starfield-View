# Functions to read, filter sky catalogue data
# Author: Ravi Ram

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from datetime import datetime, timedelta
from sgp4.io import twoline2rv
from sgp4.earth_gravity import wgs72

from star_data import *
from plot import *


"""
TLE format for OBJECT A (23057A)
https://celestrak.org/NORAD/elements/gp.php?FORMAT=TLE&CATNR=23057A
"""

#OBJECT A                
line1 = '1 56308U 23057A   23125.81960475  .00000660  00000+0  00000+0 0  9992'
line2 = '2 56308   9.8856 212.7161 0023120  32.8725 327.3113 14.88280559  2007'

#RISAT-2B                
#line1 = '1 44233U 19028A   23134.10174731  .00006152  00000+0  50662-3 0  9998'
#line2 = '2 44233  36.9998 199.2877 0014452 250.9853 108.9323 14.98747358218208'

# get satellite object from TLE (2 lines data)
def get_satellite(line1, line2):
    global mu, r, a
    
    # create satellite object from TLE
    satellite = twoline2rv(line1, line2, wgs72)
    # constants
    mu = satellite.mu           # Earth’s gravitational parameter (km³/s²)
    r = satellite.radiusearthkm # Radius of the earth (km).
    # orbital parameters
    a = satellite.a * r
    #e = satellite.ecco
    apo, peri = satellite.alta * r, satellite.altp * r
    # perigee and apogee
    print('Perigee : %5.2f km, Apogee : %5.2f km' % (peri, apo))
    # return
    return satellite

# get ra and dec from state vectors
def get_ra_dec_from_sv(r, v):
    # norm
    r_n = np.linalg.norm(r)
    # direction cosines
    l = r[0]/r_n; m = r[1]/r_n; n = r[2]/r_n;   
    # declination
    delta = np.arcsin(n)*180/np.pi                    
    # right ascension
    np.cosd = lambda x : np.cos( np.deg2rad(x) )
    if m >0:  alfa = np.arccos(l/np.cosd(delta))*180/np.pi
    else: alfa = 360 - np.arccos(l/np.cosd(delta))*180/np.pi
    # return
    return alfa, delta,

# returns a list of state vectors, ra, dec for a
# given sgp4 satellite object
def propagate(sat, time_start, time_end, dt):
    # time
    end = np.arange(0.0, time_end, dt)
    # list of datetime
    time_arr = time_start + end.astype('timedelta64[s]')    
    # state vectors, ra, dec for each time step
    position = []; velocity = []
    right_ascension = []; declination = []
    for j in time_arr.tolist():
        p, v = sat.propagate(j.year, j.month, j.day, j.hour, j.minute, j.second)
        ra, dec = get_ra_dec_from_sv(p, v)
        # list packing
        position.append(p); velocity.append(v)
        right_ascension.append(ra); declination.append(dec)

    # slice into columns
    pos, vel   = list(zip(*position)), list(zip(*velocity))
    X, Y, Z    = np.array(pos[0]), np.array(pos[1]), np.array(pos[2])
    VX, VY, VZ = np.array(vel[0]), np.array(vel[1]), np.array(vel[2])
    state_vectors = [X, Y, Z, VX, VY, VZ]
    celestial_coordinates = [np.array(right_ascension), np.array(declination)]
    # return
    return time_arr, state_vectors, celestial_coordinates

# get list of star data in view along with satellite state_vectors
def get_simulation_data(sat, df, start_time, sim_secs, time_step, roll=False):
    # state_vectors, celestial_coordinates
    tr, sc, cc = propagate(sat, start_time, sim_secs, time_step)
    # parse celestial_coordinates
    ra, dec = cc

    # [TESTING] Roll about velocity direction 
    if roll:
        roll_rate_hrs = 0.1        # deg per hr
        roll_rate_sec = 0.1/3600.0 # deg per sec
        # modify ra with roll rate
        ra = [ri + roll_rate_sec * i for i, ri in enumerate (ra)]

    # Create an empty all frames data
    frame_row_list = []
    for frame, (r, d) in enumerate(zip(ra, dec)):
        #print(frame, r,d)
        tdf_values = filter_by_fov(df, r, d).values.tolist()
        #print(tdf_values)
        frame_row_list.append([frame, tdf_values])   
    return tr, sc, frame_row_list,


def main():
    global data
    
    # create satellite object
    satellite = get_satellite(line1, line2)
    # read star data
    df = read_hipparcos_data()
    
    # time period for one revolution
    t_period = 2 * np.pi * (a**3/mu)**0.5
    # each time slice
    t_slice = 50
    # no of frames
    t_step = t_period / t_slice
    # simulation starts from current time to one full orbit
    start = np.datetime64(datetime.datetime.now())
    # times, state_vectors, celestial_coordinates  
    time_arr, state_vectors, celestial_coordinates = get_simulation_data(satellite, df, start, t_period, t_slice)
    # animate
    animate(time_arr, state_vectors, celestial_coordinates, r)
    return

# main
if __name__ == '__main__':
    main()
    
    