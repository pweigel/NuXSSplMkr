import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib import rcParams
import scipy.integrate

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
    # energies = []
    dsdy = []
    with open(fname, 'r') as f:
        lines = f.readlines()
        energies = np.array([float(x) for x in lines[0].rstrip('\n').split(',')[1:]]) / 1e9
        yvals = np.array([float(x) for x in lines[1].rstrip('\n').split(',')[1:]])
        for line in lines[2:]:
            d = [float(x) for x in line.rstrip('\n').split(',')]
            # energies.append(d[0])
            dsdy.append(np.array(d))
    return np.array(energies), np.array(yvals), np.array(dsdy)


if __name__ == '__main__':
    base_path = '/work/submit/pweigel/icecube/src/NuXSSplMkr/data/CT18A_NNLO/replica_0/cross_sections'
    
    energies, y_values, dsdy = load_file(base_path + '/dsdy_neutrino_proton_total.out')
    # print(energies.shape)
    # print(dsdy.shape)
    
    fig = plt.figure(figsize=(12,9))
    ax = fig.add_subplot(111)
    
    print(energies[-1])
    energies, y_values, dsdy = load_file(base_path + '/dsdy_antineutrino_proton_total.out')
    proton = scipy.integrate.simpson(dsdy[0, :], x=y_values)
    # print(proton)

    ax.plot(y_values, dsdy[0, :])
    ax.plot(y_values, dsdy[10, :])
    ax.plot(y_values, dsdy[20, :])
    energies, y_values, dsdy = load_file(base_path + '/dsdy_antineutrino_neutron_total.out')
    neutron = scipy.integrate.simpson(dsdy[0, :], x=y_values)
    
    ax.plot(y_values, dsdy[0, :], linestyle='--')
    ax.plot(y_values, dsdy[10, :], linestyle='--')
    ax.plot(y_values, dsdy[20, :], linestyle='--')
    
    print(0.5 * (proton + neutron))
    
    plt.savefig('4_plot.png')