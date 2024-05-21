import numpy as np
import scipy.integrate
from matplotlib import pyplot as plt
import matplotlib as mpl
from matplotlib import rcParams

from xsecs import *
from daemonflux import Flux

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

xs = CrossSections(path='/n/holylfs05/LABS/arguelles_delgado_lab/Everyone/pweigel/sandbox/src/NuXSSplMkr/old_data/CT18A_NNLO',
                   pdfname='CT18ANNLO',
                   nreps=59, A=1, Z=0.5, target_name=r'$\textrm{Isoscalar}$')

fig = plt.figure(figsize=(12,9))
ax = fig.add_subplot(111)

nu_i_cen, nu_i_minus, nu_i_plus = xs.total_xs['nu']['total']
nubar_i_cen, nubar_i_minus, nubar_i_plus = xs.total_xs['nubar']['total']
E = xs.total_xs['E']

# ax.fill_between(E, nu_i_cen - nu_i_minus, nu_i_cen + nu_i_plus, color='b', alpha=0.5)
ax.plot(E, nu_i_cen, color='k', label=r'$\nu$')
# ax.fill_between(E, nubar_i_cen - nubar_i_minus, nubar_i_cen + nubar_i_plus, color='r', alpha=0.5)
ax.plot(E, nubar_i_cen, color='k', linestyle='--', label=r'$\bar{\nu}$')

csms_nu_data = np.loadtxt('data/csms_total_xs_nu.txt', skiprows=1)
csms_nubar_data = np.loadtxt('data/csms_total_xs_nubar.txt', skiprows=1)
cteq_nu_data = np.loadtxt('data/cteq_total_xs_nu.txt', skiprows=1)
cteq_nubar_data = np.loadtxt('data/cteq_total_xs_nubar.txt', skiprows=1)

ax.plot(csms_nu_data[:, 0], csms_nu_data[:, 1] * 1e-36, color='magenta', linewidth=1)
ax.plot(csms_nubar_data[:, 0], csms_nubar_data[:, 1] * 1e-36, color='magenta', linewidth=1, linestyle='--')
ax.plot(cteq_nu_data[:, 0], cteq_nu_data[:, 1] * 1e-36, color='orange', linewidth=1)
ax.plot(cteq_nubar_data[:, 0], cteq_nubar_data[:, 1] * 1e-36, color='orange', linewidth=1, linestyle='--')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_ylim([1e-38, 1e-32])
ax.set_xlim([1e2, 1e9])
ax.legend()
plt.tight_layout()
plt.savefig('1_cross_section.png', dpi=300)