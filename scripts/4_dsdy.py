import numpy as np
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

def load_file(fname):
    energies = []
    dsdy = []
    with open(fname, 'r') as f:
        for line in f.readlines():
            d = [float(x) for x in line.rstrip('\n').split(',')]
            energies.append(d[0])
            dsdy.append(np.array(d[1:]))
    return np.array(energies), 10**np.array(dsdy)


if __name__ == '__main__':
    base_path = '/n/holylfs05/LABS/arguelles_delgado_lab/Everyone/pweigel/sandbox/src/NuXSSplMkr/data/CT18A_NNLO/replica_0/cross_sections'
    
    energies, dsdy = load_file(base_path + '/dsdy_neutrino_proton_light.out')
    print(energies.shape)
    print(dsdy.shape)
    
    y_values = np.logspace(-6, 0, 100)
    fig = plt.figure(figsize=(12,9))
    ax = fig.add_subplot(111)
    
    print(energies[50])
    ax.plot(y_values, dsdy[50, :])
    ax.plot(y_values, dsdy[60, :])
    ax.plot(y_values, dsdy[70, :])
    
    plt.savefig('4_plot.png')