# Functions to read and filter sky catalogue data
# Author: Ravi Ram

import datetime
import numpy as np
import pandas as pd

# Hipparcos Catalogue [hip_main.dat]
# http://cdsarc.u-strasbg.fr/ftp/cats/I/239/ 
FILENAME = r'hip_main.dat'

# read hipparcos catalogue 'hip_main.dat'
def read_hipparcos_data(fname=r'hip_main.dat', threshold=6.5):
    # Field H1: Hipparcos Catalogue (HIP) identifier
    # Field H5: V magnitude
    # Fields H8–9:  The right ascension, α , and declination, δ (in degrees)    
    try:
        df = pd.read_csv(FILENAME, header=None,
                         sep = '|', skipinitialspace=True).iloc[:, [1, 5, 8, 9]]
        df.columns = ['hip', 'mag', 'ra_deg', 'de_deg']

        df['mar_size'] = 2*(threshold - df['mag'])
        # filter data above
        q = 'mag <= @threshold'
        df = df.query(q) 

        # convert all columns of DataFrame
        df = df.apply(pd.to_numeric)
        
        return df    
    except FileNotFoundError:
        print("File not found.")
    
# camera fov : 9.31◦ × 7◦
def filter_by_fov(mdf, ra, de): 
    # frame field of view
    width, height = 9.31, 7.0
    # get valid boundaries  
    xmin, ymin, xmax, ymax = get_frame_boundaries( width, height, ra, de)
    # extract useful columns
    mdf = mdf[['ra_deg', 'de_deg', 'mar_size']]     
    # filter data within the boundaries    
    q = 'ra_deg >= @xmin & ra_deg <= @xmax & de_deg >= @ymin & de_deg <= @ymax' 
    mdf = mdf.query(q)
    # return filtered data
    return mdf

# get valid frame boundaries
def get_frame_boundaries(w, h, x, y):
    # set x boundaries
    xmin = float(x) - float(w) / 2.0
    xmin = 0 if xmin<0 else xmin   
    if xmin==0: xmax = w
    else: xmax = float(x) + float(w) / 2.0
    xmax = 360 if xmax>360 else xmax
    if xmax==360: xmin = 360-w
    
    # set y boundaries
    ymin = float(y) - float(h) / 2.0
    ymin = -90 if ymin<-90 else ymin
    if ymin==-90: ymax = (-90 + float(h))
    else: ymax = float(y) + float(h) / 2.0
    if ymax == 90: ymin = (90 - float(h))
    
    # return limits 
    return  xmin, ymin, xmax, ymax
