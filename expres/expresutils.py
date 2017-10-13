"""
EXPRES Utils for ceres
"""

import numpy as np
import astropy.io.fits as pyfits

def OverscanTrim(d):
    """
    Overscan correct and trim an EXPRES image. d is the input data.
    You need to set the polynomial order and the first overscan pixel of the CCD.
    """

    # decide on the polynomial order
    poly_order = 1

    # find the coefficients for the polynomial fit
    oscan_fit_coeffs = np.polyfit(range(len(d[:, 0])),
                                  np.mean(d[:, 5786:], axis=1), poly_order)

    # define a polynomial with the previously obtained coefficients
    oscan_function = np.poly1d(oscan_fit_coeffs)

    # calculate the values of the fit at each pixel so that it can be subtracted
    oscan_fit = oscan_function(range(len(d[:, 0])))

    # this trims the image of the overscan region and subtracts the overscan fit.
    # (with appropriate array broadcasting)
    oscan_corrected_data = d[:, :5786] - oscan_fit[:, np.newaxis]

    return oscan_corrected_data

def MedianCombine(ImgList):
    """
    Median combine a list of images.
    """

    # find number of images in the list, raise error if empty

    assert ImgList != [], "empty list provided!"

    # Read images and subtract overscan from images in list

    overscan_subtracted = [OverscanTrim(pyfits.getdata(image))
                           for image in ImgList]
    median_combined = np.median(overscan_subtracted, axis=0)

    return median_combined
