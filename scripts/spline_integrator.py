import numpy as np
import operator
import re
import warnings

from glob import glob

import photospline
from photospline import glam_fit, ndsparse, bspline
from photospline import SplineTable
import scipy.integrate
import os

import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib import rcParams

rcParams['font.family']           = 'serif'
rcParams['font.serif']            = 'Computer Modern Roman'
rcParams['font.weight']           = 200
rcParams['font.size']             = 30
rcParams['text.usetex']           = True

rcParams['grid.color']            = 'black'
rcParams['grid.alpha']            = 0.10
rcParams['grid.linestyle']        = '-'

rcParams['axes.grid']             = False
rcParams['axes.linewidth']        = 1.5
rcParams['axes.labelpad']         = 14.0
rcParams['axes.labelsize']        = 30.0

rcParams['xtick.labelsize']       = 30
rcParams['ytick.labelsize']       = 30
rcParams['xtick.direction']       = 'in'
rcParams['ytick.direction']       = 'in'
rcParams['xtick.major.pad']       = 10
rcParams['ytick.major.pad']       = 10
rcParams['xtick.major.width']     = 1.5
rcParams['ytick.major.width']     = 1.5
rcParams['xtick.top']             = True
rcParams['ytick.right']           = True
rcParams['xtick.minor.visible']   = True
rcParams['ytick.minor.visible']   = True

rcParams['legend.frameon']        = False
rcParams['legend.title_fontsize'] = 20
rcParams['legend.fontsize']       = 20


def ModLog10(x):
    if x <= 0. :
        return -50.0
    else:
        return np.log10(x)
ModLog10 = np.vectorize(ModLog10)

def spline_integrator(in_path, dsigma_dy_path='dsdy.fits', sigma_path='sigma.fits'):
    data = np.loadtxt(in_path, skiprows=3, delimiter=',')
    # print(datas.shape)
        
    # Let's get the E, y, x vectors
    with open(in_path, 'r') as f:
        energies = np.array([float(x) for x in f.readline().rstrip('\n').split(',')[1:]])
        y_values = np.array([float(x) for x in f.readline().rstrip('\n').split(',')[1:]])
        x_values = np.array([float(x) for x in f.readline().rstrip('\n').split(',')[1:]])
    
    data = data.reshape(len(energies), len(y_values), len(x_values))
    
    dsdy = np.zeros( (len(energies), len(y_values)) )
    for i, E in enumerate(energies):
        for j, y in enumerate(y_values):
            CUT = 1.0 / (2 * 0.938 * (E/1e9) * y)
            _x_cut = x_values[x_values > CUT]
            M = len(x_values) - len(_x_cut)
            
            if len(_x_cut) == 0:
                dsdy[i, j] = 0.0
            else:
                dsdy[i, j] = scipy.integrate.simpson(data[i, M:, j], _x_cut)
    
    sigma = np.zeros((len(energies)))
    for i, E in enumerate(energies):
        sigma[i] = scipy.integrate.simpson(dsdy[i, :], y_values)
      
    # Make total xs spline
    x = ModLog10(energies) - 9.0  # Convert to GeV
    z = ModLog10(sigma)

    knots = [np.linspace(x.min()-2,x.max()+2,80,endpoint = True)]
    order = [2]
    smooth = [1.0e-15]
    penaltyorder = [2]

    weights = np.ones(z.shape)
    zs, weights = ndsparse.from_data(z, weights)
    sigma_spline = photospline.glam_fit(zs, weights, [x], knots, order, smooth, penaltyorder)
    sigma_spline.write(sigma_path)
    
    # Make dsdy spline
    x = ModLog10(energies) - 9.0 # convert to GeV
    y = ModLog10(y_values)
    z = dsdy
    num_data_points = np.sum(z > 0)

    nd_data = photospline.ndsparse(int(num_data_points), 2)
    for i in range(z.shape[0]):
        for j in range(z.shape[1]):
            if z[i, j] != 0:
                nd_data.insert(np.log10(z[i, j]), [i, j])

    knots = [np.linspace(x.min()-2,x.max()+2, 80,endpoint = True),
             np.linspace(y.min()-2,y.max()+2, 90,endpoint = True)]
    
    order = [2, 2]
    smooth = [1.0e-15, 1.0e-15]
    penaltyorder = [2, 2]

    weights = np.ones(num_data_points)
    dsdy_spline = photospline.glam_fit(nd_data, weights, [x, y], knots, order, smooth, penaltyorder)
    dsdy_spline.write(dsigma_dy_path)
    
    fig = plt.figure(figsize=(12,9))
    ax = fig.add_subplot(111)
    # print(energies[10]/1e9)
    ax.scatter(y_values, dsdy[20, :], s=5)
    selected_energy = np.log10(energies[20])-9.0
    print(selected_energy)
    spline_y = np.linspace(-3, 0, 501)
    ax.plot(10**spline_y, 10**dsdy_spline.evaluate_simple([selected_energy, spline_y]), linewidth=3)
    ax.plot()

    ax.set_xlim([0, 1])
    plt.savefig('test_integrator_dsdy.pdf')

    fig = plt.figure(figsize=(12, 9))
    ax = fig.add_subplot(111)
    ax.scatter(energies/1e9, sigma, s=5)
    spline_energies = np.linspace(2, 9, 101)
    ax.plot(10**spline_energies, 10**sigma_spline.evaluate_simple([spline_energies]))
    ax.set_xscale('log')
    ax.set_yscale('log')
    
    plt.savefig('test_integrator_sigma.pdf')

if __name__ == '__main__':
    in_path = '/n/holylfs05/LABS/arguelles_delgado_lab/Everyone/pweigel/sandbox/src/NuXSSplMkr/data/CSMS/cross_sections/FINAL/dsdxdy_nubar_CC_iso.out'
    dsigma_dy_path = '/n/holylfs05/LABS/arguelles_delgado_lab/Everyone/pweigel/sandbox/src/NuXSSplMkr/data/CSMS/cross_sections/FINAL/dsdy_nubar_CC_iso.fits'
    sigma_path = '/n/holylfs05/LABS/arguelles_delgado_lab/Everyone/pweigel/sandbox/src/NuXSSplMkr/data/CSMS/cross_sections/FINAL/sigma_nubar_CC_iso.fits'
    spline_integrator(in_path, dsigma_dy_path, sigma_path)