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
    # print(energies)
    # print(y_values)
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

def SplineFitMaker3D(filenames, outfile, scale = 'lin', prefix = '', skip_header = 3,  N = 50, outname = "", oscale = 'lin', mod_knots=False, factor=1.0):
    """
    Creates a spline table from a table of numbers. Then
    saves the table using the sanem filename as given.

    Options:
    scale : linear or log. For logbicubic splines or simple bicubic splines
    prefix : prefix for outputfile
    skip_header : skipe lines in input datafile
    column : z = f(x,y,w), asummes x/y/w to be the first/second/third column.
    """
    # print('Infile: {}'.format(filename))
    print('Outfile: {}'.format(outfile))

    datas = None
    for infile in filenames:
        if datas is None:
            datas = np.loadtxt(infile, skiprows=skip_header, delimiter=',') * factor
        else:
            datas += np.loadtxt(infile, skiprows=skip_header, delimiter=',') * factor
        
    # Let's get the E, y, x vectors
    with open(filenames[0], 'r') as f:
        energies = np.array([float(x) for x in f.readline().rstrip('\n').split(',')[1:]])
        y_values = np.array([float(x) for x in f.readline().rstrip('\n').split(',')[1:]])
        x_values = np.array([float(x) for x in f.readline().rstrip('\n').split(',')[1:]])
    
    datas = datas.reshape(len(energies), len(y_values), len(x_values))

    x = np.log10(energies) - 9.0 # convert to GeV
    y = np.log10(y_values) # x
    w = np.log10(x_values) # y
    z = datas       # dsigma_dxdy

    print(x.shape, y.shape, z.shape)

    # knots = [np.linspace(x.min()-1,x.max()+1,N[0],endpoint = True),
    #          np.linspace(y.min()-2,y.max()+2,N[1],endpoint = True),
    #          np.linspace(w.min()-2,w.max()+2,N[2],endpoint = True)]
    
    # E_knots = np.linspace(x.min()-1,x.max()+1,N[0],endpoint = True)
    # y_knots = np.concatenate([np.linspace(y.min()-1,-1,N[1]-20,endpoint=False), np.log10(np.linspace(0.1, 1.0, 15)), np.log10([1.1, 1.2, 1.3, 1.5, 2.0])])
    # x_knots = np.concatenate([np.linspace(w.min()-1,-1,N[2]-20,endpoint=False), np.log10(np.linspace(0.1, 1.0, 15)), np.log10([1.1, 1.2, 1.3, 1.5, 2.0])])

    E_knots = np.linspace(0.0, 14, N[0], endpoint=True)
    y_knots = np.linspace(-14, 1, N[1], endpoint=True)
    x_knots = np.linspace(-14, 1, N[2], endpoint=True)

    knots = [E_knots, y_knots, x_knots]
    order = [2, 2, 2]
    # smooth = [1.0e-10, 1.0e-2, 1.0e-2]
    # smooth = [1e-15, 1e-15, 1e-15]
    smooth = [1.0e-10, 1.0e-10, 1.0e-10]
    penaltyorder = [2, 2, 2]

    # num_data_points = len(datas.flatten())
    num_data_points = np.sum(datas>0)

    # Q2_array = np.zeros(z.shape)
    # for i in range(z.shape[0]):
    #     for j in range(z.shape[1]):
    #         for k in range(z.shape[2]):
    #             s = 2 * 0.938 * (energies[i] / 1e9)
    #             Q2 = s * x_values[j] * y_values[k]
    #             Q2_array[i, j, k] = Q2
    
    nd_data = photospline.ndsparse(int(num_data_points), 3)
    for i in range(z.shape[0]):
        for j in range(z.shape[1]):
            for k in range(z.shape[2]):
                # s = 2 * 0.938 * (energies[i] / 1e9)
                # Q2 = Q2_array[i, j, k]
                if (z[i, j, k] > 0.0):
                    nd_data.insert(np.log10(z[i, j, k]), [i, j, k])
                # else:
                    # nd_data.insert(-44, [i, j, k])
    
    weights = np.ones(num_data_points)
    result = photospline.glam_fit(nd_data, weights, [x, y, w], knots, order, smooth, penaltyorder)
    result.write(outfile)

if __name__ == "__main__":
    # infile = '/n/holylfs05/LABS/arguelles_delgado_lab/Everyone/pweigel/sandbox/src/NuXSSplMkr/data/CT18A_NNLO/replica_0/cross_sections/total_neutrino_proton_bottom.out'
    # outfile = 'test_spline_bottom.fits'
    # SplineFitMaker1D(infile, outfile)
    
    # infile = '/n/holylfs05/LABS/arguelles_delgado_lab/Everyone/pweigel/sandbox/src/NuXSSplMkr/data/CT18A_NNLO/replica_0/cross_sections/dsdy_neutrino_proton_light.out'
    # outfile = '2d_spline_light.fits'
    # SplineFitMaker2D(infile, outfile, scale='log', skip_header=2, N=[200, 100])
    
    # infile = '/n/holylfs05/LABS/arguelles_delgado_lab/Everyone/pweigel/sandbox/src/NuXSSplMkr/data/CT18A_NNLO/replica_0/cross_sections/rc_dsdxdy_neutrino_proton_light.out'
    # outfile = 'rc_3d_spline_light.fits'
    # SplineFitMaker3D(infile, outfile, scale='log', oscale='log', skip_header=3, N=[50, 50, 50])
    
    for current in ['NC']:
        for projectile in ['neutrino', 'antineutrino']:
            for target in ['proton', 'neutron', 'isoscalar']:
                    
                # current = 'CC'
                # projectile = 'neutrino'
                # # target = 'proton'
                # target = 'isoscalar'
                
                
                suff = '4'
                if current == 'NC':
                    suff = '1'
                
                infiles = []
                base_path = '/n/holylfs05/LABS/arguelles_delgado_lab/Everyone/pweigel/sandbox/src/NuXSSplMkr/data/CT18ANNLO_nf6_sx_v2/replica_0/cross_sections'
                flavors = ['light', 'charm', 'bottom', 'top']
                for flavor in flavors:
                    if target == 'isoscalar':
                        infiles.append(base_path + '/' + f'dsdxdy_{current}_{projectile}_proton_{flavor}.{suff}.out')
                        infiles.append(base_path + '/' + f'dsdxdy_{current}_{projectile}_neutron_{flavor}.{suff}.out')
                    else:
                        infiles.append(base_path + '/' + f'dsdxdy_{current}_{projectile}_{target}_{flavor}.{suff}.out')
                print(infiles)
                outfile = f'wcg24_dsdxdy_{current}_{projectile}_{target}.fits'
                
                factor = 1.0
                if target=='isoscalar':
                    factor = 0.5

                SplineFitMaker3D(infiles, outfile, scale='log', oscale='log', skip_header=3, N=[12*6, 70, 70], factor=factor)
