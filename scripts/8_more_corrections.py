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
    spline_path = '/n/home06/pweigel/utils/xs_iso/dsdxdy_nu_CC_iso.fits'
    spline = photospline.SplineTable(spline_path)
    
    spline_corr_path = '/n/holylfs05/LABS/arguelles_delgado_lab/Everyone/pweigel/sandbox/src/NuXSSplMkr/data/CSMS/cross_sections/dsdxdy_nu_CC_iso_corrected.fits'
    spline_corr = photospline.SplineTable(spline_corr_path)
    
    rc_path = '/n/holylfs05/LABS/arguelles_delgado_lab/Everyone/pweigel/sandbox/src/NuXSSplMkr/bin/dsdxdy_nu_CC_iso_rc.out'
    xs_rc = np.loadtxt(rc_path, skiprows=3, delimiter=',')
    
    corr_path = '/n/holylfs05/LABS/arguelles_delgado_lab/Everyone/pweigel/sandbox/src/NuXSSplMkr/data/CSMS/cross_sections/dsdxdy_nu_CC_iso_corrected.out'
    xs_corr = np.loadtxt(corr_path, skiprows=3, delimiter=',')

    with open(rc_path, 'r') as f:
        energies = np.array([float(x) for x in f.readline().rstrip('\n').split(',')[1:]])
        y_values = np.array([float(x) for x in f.readline().rstrip('\n').split(',')[1:]])
        x_values = np.array([float(x) for x in f.readline().rstrip('\n').split(',')[1:]])
    
    xs_rc = xs_rc.reshape(len(energies), len(y_values), len(x_values))
    xs_corr = xs_corr.reshape(len(energies), len(y_values), len(x_values))
    
    selected_energy = np.log10(energies[30]) - 9.0
    print('ENERGY = {} GeV'.format(10**selected_energy))
    xs_rc = xs_rc[30, :, :]
    xs_corr = xs_corr[30, :, :]
    
    xs_spline = np.zeros((len(y_values), len(x_values)))
    for i, y in enumerate(y_values):
        for j, x in enumerate(x_values):
            xs_spline[i, j] = 10**spline.evaluate_simple([selected_energy, np.log10(x), np.log10(y)])
    
    xs_corr_spline = np.zeros((len(y_values), len(x_values)))
    for i, y in enumerate(y_values):
        for j, x in enumerate(x_values):
            v = 10**spline_corr.evaluate_simple([selected_energy, np.log10(x), np.log10(y)])
            if v < 1e-50 or v > 1e-28:
                v = 0
            xs_corr_spline[i, j] = v
            

    fig = plt.figure(figsize=(12*1.1,9*1.1))
    ax = fig.add_subplot(111)
    _x, _y = np.meshgrid(x_values, y_values)
    ax_img = ax.pcolormesh(_x, _y, 100 * (xs_corr / xs_spline - 1.0), vmin=-100.0, vmax=100.0, cmap='RdBu', rasterized=True)
    
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim([1e-5, 1])
    ax.set_ylim([1e-5, 1])
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    
    plt.colorbar(ax_img, ax=ax, label=r'$\frac{d^2 \sigma^{(1)}}{dxdy} / \frac{d^2 \sigma^{(0)}}{dxdy}$ $[\%]$')
    plt.tight_layout()
    plt.savefig('8_more_1.pdf')
    
    fig = plt.figure(figsize=(12, 9))
    ax = fig.add_subplot(111)
    _x, _y = np.meshgrid(x_values, y_values)
    ax_img = ax.pcolormesh(_x, _y, xs_spline, norm=mpl.colors.LogNorm(vmin=1e-42)) #
    
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim([1e-5, 1])
    ax.set_ylim([1e-5, 1])
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    # ax.tick_params(axis='x', which='both')
    # ax.set_aspect(1)
    
    plt.colorbar(ax_img, ax=ax)
    plt.tight_layout()
    plt.savefig('8_more_2.pdf')
    
    fig = plt.figure(figsize=(12,9))
    ax = fig.add_subplot(111)
    _x, _y = np.meshgrid(x_values, y_values)
    ax_img = ax.pcolormesh(_x, _y, xs_corr, norm=mpl.colors.LogNorm(vmin=1e-42))
    
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim([1e-5, 1])
    ax.set_ylim([1e-5, 1])
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    
    plt.colorbar(ax_img, ax=ax)
    plt.tight_layout()
    plt.savefig('8_more_3.pdf')
    
    fig = plt.figure(figsize=(12,9))
    ax = fig.add_subplot(111)
    _x, _y = np.meshgrid(x_values, y_values)
    ax_img = ax.pcolormesh(_x, _y, 100.0 * (xs_corr_spline / xs_spline - 1.0), vmin=-100, vmax=100, cmap='RdBu')
    
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim([1e-5, 1])
    ax.set_ylim([1e-5, 1])
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    
    plt.colorbar(ax_img, ax=ax)
    plt.tight_layout()
    plt.savefig('8_more_4.pdf')
    
    
        