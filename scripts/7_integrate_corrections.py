import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib import rcParams
import photospline
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

if __name__ == '__main__':
    base_path = '/n/holylfs05/LABS/arguelles_delgado_lab/Everyone/pweigel/sandbox/src/NuXSSplMkr/data/CT18A_NNLO/replica_0/cross_sections'
    spline_path = '/n/home06/pweigel/utils/xs_iso/dsdxdy_nu_CC_iso.fits'
    spline = photospline.SplineTable(spline_path)
    
    rc_path = '/n/holylfs05/LABS/arguelles_delgado_lab/Everyone/pweigel/sandbox/src/NuXSSplMkr/bin/dsdxdy_nu_CC_iso_rc.out'
    xs_rc = np.loadtxt(rc_path, skiprows=3, delimiter=',')
    
    print(xs_rc.shape)
    with open(rc_path, 'r') as f:
        energies = np.array([float(x) for x in f.readline().rstrip('\n').split(',')[1:]])
        y_values = np.array([float(x) for x in f.readline().rstrip('\n').split(',')[1:]])
        x_values = np.array([float(x) for x in f.readline().rstrip('\n').split(',')[1:]])
    
    xs_rc = xs_rc.reshape(len(energies), len(y_values), len(x_values))
    
    selected_energy = np.log10(energies[30]) - 9.0
    print('ENERGY = {} GeV'.format(10**selected_energy))
    xs_rc = xs_rc[30, :, :]
    
    xs_spline = np.zeros((len(y_values), len(x_values)))
    for i, y in enumerate(y_values):
        for j, x in enumerate(x_values):
            xs_spline[i, j] = 10**spline.evaluate_simple([selected_energy, np.log10(x), np.log10(y)])
    
    dsdy = np.zeros(len(y_values))
    for i, y in enumerate(y_values):
        dsdy[i] = scipy.integrate.simpson(xs_spline[i, :], x_values)
    
    dsdy_rc = np.zeros(len(y_values))
    for i, y in enumerate(y_values):
        dsdy_rc[i] = scipy.integrate.simpson(xs_rc[i, :], x_values)
    
    dsdy_tot = dsdy + dsdy_rc
    total_xs = scipy.integrate.simpson(dsdy_tot, y_values)
    
    fig = plt.figure(figsize=(12, 9))
    ax = fig.add_subplot(111)
    ax.plot(y_values, dsdy / total_xs, color='r', linewidth=3, label=r'$\frac{d\sigma^{(0)}}{dy}$')
    ax.plot(y_values, dsdy_rc / total_xs, color='b', linewidth=3, label=r'$\frac{d\sigma^{(1)}}{dy}$')
    ax.plot(y_values, (dsdy + dsdy_rc) / total_xs, color='k', linewidth=3, label=r'$\frac{d\sigma^{(0)}}{dy}+\frac{d\sigma^{(1)}}{dy}$')
    ax.plot([0, 1], [0, 0], color='k', alpha=0.33, linestyle='--')
    ax.set_xlim([0, 1])
    ax.set_xlabel(r'$y$')
    ax.set_ylabel(r'$\frac{1}{\sigma}\frac{d\sigma}{dy}$')
    ax.legend()
    plt.tight_layout()
    plt.savefig('7_dsdy_rc.pdf')