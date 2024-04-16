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
    sigma = []
    with open(fname, 'r') as f:
        for line in f.readlines():
            d = [float(x) for x in line.rstrip('\n').split(',')]
            energies.append(d[0])
            sigma.append(np.array(d[1:]))
    return np.array(energies), np.array(sigma)


if __name__ == '__main__':
    base_path = '/n/holylfs05/LABS/arguelles_delgado_lab/Everyone/pweigel/sandbox/src/NuXSSplMkr/data/CT18A_NNLO/replica_0/cross_sections'
    
    fig = plt.figure(figsize=(12, 9))
    ax = fig.add_subplot(111)
    
    energies, xs = load_file(base_path + '/total_neutrino_proton_light.out')
    spline_path = 'test_spline_light.fits'
    spline = photospline.SplineTable(spline_path)
    ax.scatter(energies / 1e9, xs, s=5)
    ax.plot(energies / 1e9, 10**spline.evaluate_simple([np.log10(energies / 1e9)]), color='k')
    
    
    energies, xs = load_file(base_path + '/total_neutrino_proton_charm.out')
    spline_path = 'test_spline_charm.fits'
    spline = photospline.SplineTable(spline_path)
    ax.scatter(energies / 1e9, xs, s=5)
    ax.plot(energies / 1e9, 10**spline.evaluate_simple([np.log10(energies / 1e9)]), color='k')
    
    energies, xs = load_file(base_path + '/total_neutrino_proton_bottom.out')
    spline_path = 'test_spline_bottom.fits'
    spline = photospline.SplineTable(spline_path)
    ax.scatter(energies / 1e9, xs, s=5)
    ax.plot(energies / 1e9, 10**spline.evaluate_simple([np.log10(energies / 1e9)]), color='k')
    
    energies, xs = load_file(base_path + '/total_neutrino_proton_top.out')
    spline_path = 'test_spline_top.fits'
    spline = photospline.SplineTable(spline_path)
    
    ax.scatter(energies / 1e9, xs, s=5)
    top_energies = np.logspace(3, 9, 1001)
    top_xs = 10**spline.evaluate_simple([np.log10(top_energies)])
    ax.plot(top_energies, top_xs, color='k')
    
    ax.set_xscale('log')
    ax.set_yscale('log')
    
    print(xs)

    plt.savefig('5_1d_spline.pdf')
    
