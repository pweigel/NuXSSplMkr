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

    weights = np.ones(z.shape)
    zs, weights = ndsparse.from_data(z, weights)
    result = photospline.glam_fit(zs, weights, [x], knots, order, smooth, penaltyorder)
    result.write(outfile)
    print("Done. Generated: "  + outfile)

def SplineFitMaker2D(filename, outfile='out.fits', scale = 'lin', prefix = '', skip_header = 2, column = 2, N = [100, 100]):
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

    weights = np.ones(num_data_points)
    result = photospline.glam_fit(nd_data, weights, [x, y], knots, order, smooth, penaltyorder)
    result.write(outfile)
    
    print("Done. Generated: "  + outfile)

def SplineFitMaker3D(filename, outfile, scale = 'lin', prefix = '', skip_header = 3,  N = 50, outname = "", oscale = 'lin'):
    """
    Creates a spline table from a table of numbers. Then
    saves the table using the sanem filename as given.

    Options:
    scale : linear or log. For logbicubic splines or simple bicubic splines
    prefix : prefix for outputfile
    skip_header : skipe lines in input datafile
    column : z = f(x,y,w), asummes x/y/w to be the first/second/third column.
    """


    datas = np.loadtxt(filename, skiprows=skip_header, delimiter=',')
    print(datas.shape)
    
    # Let's get the E, y, x vectors
    with open(filename, 'r') as f:
        energies = np.array([float(x) for x in f.readline().rstrip('\n').split(',')[1:]])
        y_values = np.array([float(x) for x in f.readline().rstrip('\n').split(',')[1:]])
        x_values = np.array([float(x) for x in f.readline().rstrip('\n').split(',')[1:]])
    
    datas = datas.reshape(len(energies), len(y_values), len(x_values))

    f = lambda x : x
    if scale == "log":
        f = lambda x : ModLog10(x)
    elif scale == "lin":
        pass
    else:
        print("Error: unknown scale.")
        exit()

    of = lambda x : x
    if oscale == "log":
        of = lambda x : ModLog10(x)
    elif oscale == "lin":
        pass
    else:
        print("Error: unknown scale.")
        exit()

    x = f(energies) - 9.0 # convert to GeV
    y = f(y_values)
    w = f(x_values)
    z = datas

    knots = [np.linspace(x.min()-1,x.max()+1,N[0],endpoint = True),
             np.linspace(y.min()-1,y.max()+1,N[1],endpoint = True),
             np.linspace(w.min()-1,w.max()+1,N[2],endpoint = True)]
    
    order = [2, 2, 2]
    smooth = [1.0e-15, 1.0e-15, 1.0e-5]
    penaltyorder = [2, 2, 2]
    
    num_data_points = np.sum(datas > 0.0)
    print(num_data_points)

    nd_data = photospline.ndsparse(int(num_data_points), 3)
    maxes = [0, 0, 0]
    for i in range(z.shape[0]):
        for j in range(z.shape[1]):
            for k in range(z.shape[2]):
                if z[i, j, k] > 0:
                    nd_data.insert(np.log10(z[i, j, k]), [i, j, k])
                    if i > maxes[0]: maxes[0] = i
                    if j > maxes[1]: maxes[1] = j
                    if k > maxes[2]: maxes[2] = k
      
    print(maxes, z.shape)

    weights = np.ones(num_data_points)
    result = photospline.glam_fit(nd_data, weights, [x, y, w], knots, order, smooth, penaltyorder)
    result.write(outfile)

if __name__ == "__main__":
    infile = '/n/holylfs05/LABS/arguelles_delgado_lab/Everyone/pweigel/sandbox/src/NuXSSplMkr/data/CT18A_NNLO/replica_0/cross_sections/total_neutrino_proton_bottom.out'
    outfile = 'test_spline_bottom.fits'
    # SplineFitMaker1D(infile, outfile)
    
    infile = '/n/holylfs05/LABS/arguelles_delgado_lab/Everyone/pweigel/sandbox/src/NuXSSplMkr/data/CT18A_NNLO/replica_0/cross_sections/dsdy_neutrino_proton_light.out'
    outfile = '2d_spline_light.fits'
    # SplineFitMaker2D(infile, outfile, scale='log', skip_header=2, N=[200, 100])
    
    infile = '/n/holylfs05/LABS/arguelles_delgado_lab/Everyone/pweigel/sandbox/src/NuXSSplMkr/data/CT18A_NNLO/replica_0/cross_sections/dsdxdy_neutrino_proton_light.out'
    outfile = '3d_spline_light.fits'
    SplineFitMaker3D(infile, outfile, scale='log', oscale='log', skip_header=3, N=[50, 50, 50])
