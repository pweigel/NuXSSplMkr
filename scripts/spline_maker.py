import numpy as np
import operator
import re
import warnings

from glob import glob

# Fast sparse-matrix implementation
import photospline
from photospline import glam_fit, ndsparse, bspline
from photospline import SplineTable

import os

"""
Script to make structure functions spline tables.
C.A. Arguelles Delgado - aug.03.14
"""

#warnings.filterwarnings('error')

def ModLog10(x):
    if x <= 0. :
        # print(x)
        return -50
    else:
        return np.log10(x)

ModLog10 = np.vectorize(ModLog10)

def SplineFitMaker1D(filename, outfile='out.fits', scale = 'lin', prefix = '', skip_header = 0, column = 1, N = 50, outname = "", oscale = 'lin'):
    """
    Creates a spline table from a table of numbers. Then
    saves the table using the same filename as given.

    Options:
    scale : linear or log. For logbicubic splines or simple bicubic splines
    prefix : prefix for outputfile
    skip_header : skipe lines in input datafile
    column : z = f(x), asummes x to be the first column.
    """
    if(column < 1):
        print("Error: column < 1.")
        exit()

    datas = np.loadtxt(filename, skiprows = skip_header, delimiter=',')

    x = np.log10(datas[:, 0])
    z = datas[:, 1]

    knots = [np.linspace(x.min()-1,x.max()+1,N,endpoint = True)]
    order = [2]
    smooth = [1.0e-15]
    penaltyorder = [2]

    weight = np.ones(z.shape)
    zs, w = ndsparse.from_data(z, weight)
    result = photospline.glam_fit(zs, w, [x], knots, order, smooth, penaltyorder)

    result.write(outfile)
    print("Done. Generated :"  + outfile)

if __name__ == "__main__":
    infile = '/n/holylfs05/LABS/arguelles_delgado_lab/Everyone/pweigel/sandbox/src/NuXSSplMkr/data/CT18A_NNLO/replica_0/cross_sections/total_neutrino_proton_light.out'
    outfile = 'test_spline.fits'
    SplineFitMaker1D(infile, outfile)