import numpy as np
import operator
import re
import warnings

from glob import glob

import photospline
from photospline import glam_fit, ndsparse, bspline
from photospline import SplineTable

import os

"""
Script to make structure functions spline tables.
C.A. Arguelles Delgado - aug.03.14

modified nearly 10 years later by P.L.R. Weigel - apr.15.24
"""

def ModLog10(x):
    if x <= 0. :
        return -50.0
    else:
        return np.log10(x)

ModLog10 = np.vectorize(ModLog10)

# TODO: reimplement the features...

def SplineFitMaker1D(filename, outfile='out.fits', scale = 'lin', prefix = '', skip_header = 0, column = 1, N = 200, outname = "", oscale = 'lin'):
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
    
    x = ModLog10(datas[:, 0]) - 9.0  # Convert to GeV
    z = ModLog10(datas[:, 1])

    knots = [np.linspace(x.min()-1,x.max()+1,N,endpoint = True)]
    order = [2]
    smooth = [1.0e-15]
    penaltyorder = [2]

    weight = np.ones(z.shape)
    zs, w = ndsparse.from_data(z, weight)
    result = photospline.glam_fit(zs, w, [x], knots, order, smooth, penaltyorder)
    result.write(outfile)
    print("Done. Generated: "  + outfile)

def SplineFitMaker2D(filename, outfile='out.fits', scale = 'lin', prefix = '', skip_header = 0, column = 2, N = [100, 100]):
    """
    Creates a spline table from a table of numbers. Then
    saves the table using the sanem filename as given.

    Options:
    scale : linear or log. For logbicubic splines or simple bicubic splines
    prefix : prefix for outputfile
    skip_header : skipe lines in input datafile
    column : z = f(x,y), asummes x y to be the first two columns.
    """
    if(column < 2):
        print("Error: column < 2.")
        exit()

    datas = np.loadtxt(filename, skiprows=skip_header, delimiter=',')
    print(datas.shape)
    
    # Let's get the E, y vectors
    with open(filename, 'r') as f:
        energies = np.array([float(x) for x in f.readline().rstrip('\n').split(',')[1:]])
        y_values = np.array([float(x) for x in f.readline().rstrip('\n').split(',')[1:]])
    print(energies)
    print(y_values)
    # exit()

    f = lambda x : x
    if scale == "log":
        f = lambda x : ModLog10(x)
    elif scale == "lin":
        pass
    else:
        print("Error: unknown scale.")
        exit()

    x = f(energies) - 9.0 # convert to GeV
    y = f(y_values)
    z = f(datas)
    num_data_points = np.sum(datas > 0)

    nd_data = photospline.ndsparse(int(num_data_points), 2)
    for i in range(datas.shape[0]):
        for j in range(datas.shape[1]):
            if datas[i, j] != 0:
                nd_data.insert(np.log10(datas[i, j]), [i, j])

    knots = [np.linspace(x.min()-1,x.max()+1, N[0],endpoint = True),
             np.linspace(y.min()-1,y.max()+1, N[1],endpoint = True)]
    
    order = [2, 2]
    smooth = [1.0e-15, 1.0e-15]
    penaltyorder = [2, 2]

    w = np.ones(num_data_points)
    result = photospline.glam_fit(nd_data, w, [x, y], knots, order, smooth, penaltyorder)
    result.write(outfile)
    
    print("Done. Generated: "  + outfile)

if __name__ == "__main__":
    infile = '/n/holylfs05/LABS/arguelles_delgado_lab/Everyone/pweigel/sandbox/src/NuXSSplMkr/data/CT18A_NNLO/replica_0/cross_sections/total_neutrino_proton_bottom.out'
    outfile = 'test_spline_bottom.fits'
    # SplineFitMaker1D(infile, outfile)
    
    infile = '/n/holylfs05/LABS/arguelles_delgado_lab/Everyone/pweigel/sandbox/src/NuXSSplMkr/data/CT18A_NNLO/replica_0/cross_sections/dsdy_neutrino_proton_light.out'
    outfile = '2d_spline_light.fits'
    SplineFitMaker2D(infile, outfile, scale='log', skip_header=2, N=[200, 100])
