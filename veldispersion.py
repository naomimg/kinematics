#!/usr/bin/env python
# N. McClure-Griffiths 31 May 2017

from numpy import *
from pylab import *
from astropy.io import fits
from astropy.wcs import WCS
import os
import scipy 
from scipy.interpolate import interp1d

import aplpy

# Generic function to read the terminal velocities from the data
def vterm(fname):
    csv = np.genfromtxt (fname,comments='%')
    l = csv[:,0]
    v = csv[:,1]
    for i in range(len(l)): 
        if l[i]>180.0: l[i]-=360.
    return np.array(l),np.array(v)


# Now do a spline fit to the dv and dr values and return those new arrays
#dv=map(lambda x: 50. if x>50. else x, dv) # Replace all values greater than 50 with 50.
#dv=map(lambda x: 2. if x<2. else x, dv) # Replace all values less than 2 with 2
#dv_int = interp1d(l, dv,bounds_error=False,fill_value=2.)
#dr_int = interp1d(l, dr,kind='cubic',bounds_error=False,fill_value=1.)

# open the FITS file
input_cube='/Users/naomi/Data/GASS/GASSIII/gass_-20_0.fits'
hdulist = fits.open(input_cube)
im = hdulist[0].data
hdr = hdulist[0].header

data_cube=1.*im
dimensions = data_cube.shape

idimensions=0
if dimensions[0]==1:
    idimensions=1        ## this means that the header most likely has a degenerate polarization axis

nch=dimensions[0+idimensions]       ## number of velocity channels
nypix=dimensions[1+idimensions]
nxpix=dimensions[2+idimensions]

w = WCS(input_cube)
x = zeros(nxpix)
y=zeros(nypix)
z = range(nch)
lon,lat,vel=w.wcs_pix2world(x,y,z,0)

# Read the terminal velocity values
l1,v1=vterm('term_vel_q1.dat')   # First quadrant values
l4,v4=vterm('term_vel_q4.dat')  # Fourth quadrant values



