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
    dsdxdy = []
    with open(fname, 'r') as f:
        n = 0
        for line in f.readlines():
            if n == 0:
                energies = np.array([float(x) for x in line.rstrip('\n').split(',')[1:]])
            elif n == 1:
                y_values = np.array([float(x) for x in line.rstrip('\n').split(',')[1:]])
            elif n == 2:
                x_values = np.array([float(x) for x in line.rstrip('\n').split(',')[1:]])
            else:
                d = np.array([float(x) for x in line.rstrip('\n').split(',')])
                dsdxdy.append(d)
            n += 1
    return energies / 1e9, y_values, x_values, np.array(dsdxdy)
  
if __name__ == '__main__':
    base_path = '/n/holylfs05/LABS/arguelles_delgado_lab/Everyone/pweigel/sandbox/src/NuXSSplMkr/data/CT18A_NNLO/replica_0/cross_sections'
    filename = base_path + '/rc_dsdxdy_neutrino_proton_light.out'
    xs_rc = np.loadtxt(filename, skiprows=3, delimiter=',')

    filename = base_path + '/no_rc_dsdxdy_neutrino_proton_light.out'
    xs_no_rc = np.loadtxt(filename, skiprows=3, delimiter=',')

    # fig = plt.figure(figsize=(12, 9))
    # ax = fig.add_subplot(111)
    
    print(xs_rc.shape)
    with open(filename, 'r') as f:
        energies = np.array([float(x) for x in f.readline().rstrip('\n').split(',')[1:]])
        y_values = np.array([float(x) for x in f.readline().rstrip('\n').split(',')[1:]])
        x_values = np.array([float(x) for x in f.readline().rstrip('\n').split(',')[1:]])
        
    fig = plt.figure(figsize=(12,9))
    ax = fig.add_subplot(111)
    _x, _y = np.meshgrid(x_values, y_values)
    ax_img = ax.pcolormesh(_x, _y, xs_rc.T / xs_no_rc.T, vmin=-1.0, vmax=1.0)
    
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    
    plt.colorbar(ax_img, ax=ax)
    plt.tight_layout()
    plt.savefig('6_qed_corrections_1.pdf')
    
    fig = plt.figure(figsize=(12,9))
    ax = fig.add_subplot(111)
    _x, _y = np.meshgrid(x_values, y_values)
    ax_img = ax.pcolormesh(_x, _y, xs_no_rc.T, norm=mpl.colors.LogNorm(vmin=1e-42, vmax=1e-32)) #
    
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    
    plt.colorbar(ax_img, ax=ax)
    plt.tight_layout()
    plt.savefig('6_qed_corrections_2.pdf')
    
    fig = plt.figure(figsize=(12,9))
    ax = fig.add_subplot(111)
    _x, _y = np.meshgrid(x_values, y_values)
    ax_img = ax.pcolormesh(_x, _y, xs_rc.T, norm=mpl.colors.LogNorm(vmin=1e-42))
    
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    
    plt.colorbar(ax_img, ax=ax)
    plt.tight_layout()
    plt.savefig('6_qed_corrections_3.pdf')
        