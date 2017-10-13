import sys
import time

base = '../'
sys.path.append(base+"utils/Correlation")
sys.path.append(base+"utils/GLOBALutils")
sys.path.append(base+"utils/OptExtract")

# ceres modules
import expresutils
#import correlation
import GLOBALutils
#import Marsh

# other useful modules
#import argparse
#import ephem
#import jplephem
from math import radians as rad
import pickle
import os
import numpy as np
from scipy import interpolate

from astropy.io import fits
from glob import glob

#import statsmodels.api as sm
#lowess = sm.nonparametric.lowess

#######################################################
# CHANGE THE PATH BEFORE USE
path = '/home/joel/scratch/EXPRES/Test'
files = sorted(glob(path + '/*.fits'))

flat_list = files[:-2]

print('loading science data')
science_data = expresutils.OverscanTrim(fits.getdata(files[-1]))

print('creating master flat')
# this already does overscan subtraction
master_extended_flat = expresutils.MedianCombine(flat_list)
master_extended_flat[master_extended_flat < 200] = 200
print('creating normalized master flat')
normalized_master_extended_flat = master_extended_flat/np.max(master_extended_flat)

# generate science flat
print('generating science flat')
science_flat = expresutils.OverscanTrim(fits.getdata(files[-2]))

# flat-field the science data
print('reducing science data')
reduced_science_data = science_data #/normalized_master_extended_flat

print('writing science data')
hdu = fits.PrimaryHDU(reduced_science_data)
hdulist = fits.HDUList([hdu])
hdulist.writeto('reduced_sun_spectrum_noflat.fits', overwrite=True)

# set up simple extraction parameters to get traces from science flat
print('starting simple extraction')
start = time.time()
ext_aperture = 16  # height of orders
trace_degree = 4  # degree of polynomial fit
max_orders = 18
startfrom = 0  # first column
endat = 5786  # last column
nsigmas = 5
#c_all, nord = GLOBALutils.get_them(oscan_subtracted_science_flat,ext_aperture,trace_degree,max_orders,startfrom,endat,nsigmas,mode=1)  # fit polynomial to orders
print('computing order tracing polynomials')
c_all, nord = GLOBALutils.get_them(science_flat,ext_aperture,trace_degree,mode=2)  # fit polynomial to orders
npools = 2
print('performing simple extraction')
s = GLOBALutils.simple_extraction(reduced_science_data,c_all,ext_aperture,startfrom,endat,npools)
print('writing simple extraction science data')
hdu = fits.PrimaryHDU(s)
hdulist = fits.HDUList([hdu])
hdulist.writeto('extracted_solar_spectrum_sliced_noflat.fits', overwrite=True)
end = time.time()
print('simple extraction performed in ', end-start, ' seconds') 

# optimal_extraction next
print('starting optimal extraction')
start = time.time()
ronoise = 1     # obtained from header
gain = 1        # obtained from header
NSigma_Marsh = 5
S_Marsh = 0.4
N_Marsh = 3
Marsh_alg = 0
NCosmic_Marsh = 5 
print('computing weights for all traces')
P = GLOBALutils.obtain_P(reduced_science_data,c_all,ext_aperture,ronoise,gain,NSigma_Marsh,S_Marsh,N_Marsh,Marsh_alg,startfrom,endat,npools)
print('performing optimal extraction')
s_optimal = GLOBALutils.optimal_extraction(reduced_science_data,P,c_all,ext_aperture,ronoise,gain,S_Marsh,NCosmic_Marsh,startfrom,endat,npools)
print('writing optimally extracted data')
hdu = fits.PrimaryHDU(s_optimal)
hdulist = fits.HDUList([hdu])
hdulist.writeto('optimally_extracted_solar_spectrum_marsh_sliced_noflat.fits', overwrite=True)
end = time.time()
print('optimal extraction performed in ', end-start, ' seconds')


