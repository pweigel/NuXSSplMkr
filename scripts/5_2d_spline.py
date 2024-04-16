import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib import rcParams
import photospline

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

def load_file(fname):
    energies = []
    dsdy = []
    with open(fname, 'r') as f:
        n = 0
        for line in f.readlines():
            if n == 0:
                energies = np.array([float(x) for x in line.rstrip('\n').split(',')[1:]])
            elif n == 1:
                y_values = np.array([float(x) for x in line.rstrip('\n').split(',')[1:]])
            else:
                d = np.array([float(x) for x in line.rstrip('\n').split(',')])
                dsdy.append(d)
            n += 1
    return energies / 1e9, y_values, np.array(dsdy)


if __name__ == '__main__':
    base_path = '/n/holylfs05/LABS/arguelles_delgado_lab/Everyone/pweigel/sandbox/src/NuXSSplMkr/data/CT18A_NNLO/replica_0/cross_sections'
    
    fig = plt.figure(figsize=(12, 9))
    ax = fig.add_subplot(111)
    
    energies, y_values, dsdy = load_file(base_path + '/dsdy_neutrino_proton_light.out')
    spline_path = '2d_spline_light.fits'
    spline = photospline.SplineTable(spline_path)
    E = np.log10(energies[100])
    _y_values = np.linspace(-6, 0, 1001)
    
    dxs = spline.evaluate_simple([E, _y_values])
    ax.scatter(y_values, dsdy[100, :])
    ax.plot(10**_y_values, dxs)

    ax.set_xscale('log')
    ax.set_xlim([1e-6, 1])
    ax.set_yscale('linear')
    
    plt.savefig('5_2d_spline_y.pdf')
    
    fig = plt.figure(figsize=(12, 9))
    ax = fig.add_subplot(111)
    
    E = np.linspace(1, 9, 1001)
    m = 90
    y = y_values[m]
    ax.scatter(energies, dsdy[:, m])
    print(y)
    dxs = 10**spline.evaluate_simple([E, y])
    print(dxs)

    ax.plot(10**E, dxs)

    ax.set_xscale('log')
    ax.set_yscale('log')
    
    plt.savefig('5_2d_spline_E.pdf')
