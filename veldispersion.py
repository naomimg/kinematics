#!/usr/bin/env python
# N. McClure-Griffiths 31 May 2017

from numpy import *
from pylab import *
from astropy.io import fits
from astropy.wcs import WCS
import os
import math
import scipy 
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit

import aplpy

# Generic function to read the terminal velocities from the data
def vterm(fname):
    csv = np.genfromtxt (fname,comments='%')
    l = csv[:,0]
    v = csv[:,1]
    for i in range(len(l)): 
        if l[i]>180.0: l[i]-=360.
    return np.array(l),np.array(v)



# Function to fit error functions to the spectra
def errfunc2(v,v_o1,a1,dv1, v_o2,a2,dv2):
	return a1* scipy.special.erfc((v_o1-v)/dv1) + a2* scipy.special.erfc((v_o2-v)/dv2) 



# Function to fit error functions to the spectra
def errfunc3(v,v_o1,a1,dv1, v_o2,a2,dv2,v_o3,a3,dv3):
	return a1* scipy.special.erfc((v_o1-v)/dv1) + a2* scipy.special.erfc((v_o2-v)/dv2) + a3* scipy.special.erfc((v_o3-v)/dv3) 





# Now do a spline fit to the dv and dr values and return those new arrays
#dv=map(lambda x: 50. if x>50. else x, dv) # Replace all values greater than 50 with 50.
#dv=map(lambda x: 2. if x<2. else x, dv) # Replace all values less than 2 with 2
#dv_int = interp1d(l, dv,bounds_error=False,fill_value=2.)
#dr_int = interp1d(l, dr,kind='cubic',bounds_error=False,fill_value=1.)

# open the FITS file
input_cube='/Users/naomi/Data/HI4PI/CAR_E16.fits'
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

# Get the coordinates of each axis, som fiddling here to get the correct number of pixels.
w = WCS(input_cube)
x = zeros(nch)
y=zeros(nch)
z = range(nch)
for i in range(nypix+1):
    y[i]=i
for i in range(nxpix+1):
    x[i]=i

l,b,vel=w.wcs_pix2world(x,y,z,0)
lon =l[0:nxpix]
lat = b[0:nypix]

# Read the terminal velocity values
l1,v1=vterm('term_vel_q1.dat')   # First quadrant values
l4,v4=vterm('term_vel_q4.dat')  # Fourth quadrant values

lon4 = map(lambda x: x-360., lon)  # Write the longitude values as negatives
vel = vel/1000.0
# Do the spline fit the terminal velocity values to extract the values at the longitudes measured.
f = interp1d(l4, v4,bounds_error=False,fill_value="extrapolate")
v4_int = f(lon4)  # These are now sampled on the same longitude points as the data


# Find the y_index that is at zero degrees latitude
y_indx0 = min(range(len(lat)), key=lambda i: abs(lat[i]-0.0))

# Find the channel that corresponds to the terminal velocity as measured
# pack into an array with length equal to lon
indices=[]
for j in range(nxpix):
	indices.append(min(range(len(vel)), key=lambda i: abs(vel[i]-v4_int[j])))

# Check that these look ok:
#for j in range(nxpix):
#	print lon[j],v4_int[j],vel[indices[j]]

# Now let's extract the spectrum from 60 channels around that index
nchannels = 50
k=0
spectra=np.empty([nxpix, nchannels])
for j in range(nxpix):
	imin = indices[j]-nchannels+8
	imax= indices[j]+8
	k=0
	for i in range(imin,imax):
		print k,j,i
		spectra[j,k]= data_cube[i,y_indx0,j]
		k+=1


# Bin the spectra by some amount, say 5 samples
bin = 5
newspec = np.empty([int(nxpix/bin)+1,nchannels])
newlon = np.empty([int(nxpix/bin)+1])
k = 0
newx = int(nxpix/bin)+1
for j in range(0,nxpix,bin):
	print j,k
	for i in range(nchannels):
		newspec[k,i]=np.mean(spectra[j:j+bin-1,i])
	newlon[k] = np.mean(lon[j:j+bin-1])
	k+=1


xdata = arange(0.0,50.)
a1 = []
dv1 = []
vo1 = []
a2 = []
dv2 = []
vo2 = []

for i in range(newx):
#	 Fit a spectrum
	popt, pcov = curve_fit(errfunc2, xdata, newspec[i,:],bounds=(0,[50,100,10,50,100,30]))
	tb = errfunc2(xdata,popt[0],popt[1],popt[2],popt[3],popt[4],popt[5])
	plot(newspec[i],linestyle='-')
	plot(tb,linestyle='--')
	if popt[2]<popt[5]:
		vo1.append(popt[0])
		a1.append(popt[1])
		dv1.append(popt[2])
		vo2.append(popt[3])
		a2.append(popt[4])
		dv2.append(popt[5])
	else:
		vo2.append(popt[0])
		a2.append(popt[1])
		dv2.append(popt[2])
		vo1.append(popt[3])
		a1.append(popt[4])
		dv1.append(popt[5])

a1 = []
dv1 = []
vo1 = []
a2 = []
dv2 = []
vo2 = []
a3 = []
vo3 = []
dv3=[]

for i in range(nxpix):
#	 Fit a spectrum
	popt, pcov = curve_fit(errfunc3, xdata, spectra[i,:],bounds=(0,[50,100,15,50,100,20,50,10,50]))
	tb = errfunc3(xdata,popt[0],popt[1],popt[2],popt[3],popt[4],popt[5],popt[6],popt[7],popt[8])
	plot(spectra[i],linestyle='-')
	plot(tb,linestyle='--')
	if popt[2]<popt[5] and popt[5]<popt[8]:
		vo1.append(popt[0])
		a1.append(popt[1])
		dv1.append(popt[2])
		vo2.append(popt[3])
		a2.append(popt[4])
		dv2.append(popt[5])
		vo3.append(popt[6])
		a3.append(popt[7])
		dv3.append(popt[8])
	elif (popt[2]> popt[5] and popt[5]<popt[8]):
		vo2.append(popt[0])
		a2.append(popt[1])
		dv2.append(popt[2])
		vo1.append(popt[3])
		a1.append(popt[4])
		dv1.append(popt[5])
		vo3.append(popt[6])
		a3.append(popt[7])
		dv3.append(popt[8])
	else:
		vo1.append(popt[0])
		a1.append(popt[1])
		dv1.append(popt[2])
		vo3.append(popt[3])
		a3.append(popt[4])
		dv3.append(popt[5])
		vo2.append(popt[6])
		a2.append(popt[7])
		dv2.append(popt[8])





