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

# Get my stuff
nu_i_meany_cen, nu_i_meany_minus, nu_i_meany_plus = xs.mean_y['nu']['total']
nubar_i_meany_cen, nubar_i_meany_minus, nubar_i_meany_plus = xs.mean_y['nubar']['total']
E = np.logspace(1, 9, 100)

ax.fill_between(E, nu_i_meany_cen - nu_i_meany_minus, nu_i_meany_cen + nu_i_meany_plus, color='b', alpha=0.5)
ax.plot(E, nu_i_meany_cen, color='b', label=r'$\nu$')
ax.fill_between(E, nubar_i_meany_cen - nubar_i_meany_minus, nubar_i_meany_cen + nubar_i_meany_plus, color='r', alpha=0.5)
ax.plot(E, nubar_i_meany_cen, color='r', label=r'$\bar{\nu}$')

fl = Flux(location="generic", use_calibration=True, debug=1)
numu_flux = fl.flux(E, "average", "numu")
numubar_flux = fl.flux(E, "average", "antinumu")
frac_nu = numu_flux / (numu_flux + numubar_flux)
frac_nubar = numubar_flux / (numu_flux + numubar_flux)
flux_avg = frac_nu * nu_i_meany_cen + frac_nubar * nubar_i_meany_cen
flux_avg_lower = frac_nu * (nu_i_meany_cen - nu_i_meany_minus) + frac_nubar * (nubar_i_meany_cen - nubar_i_meany_minus)
flux_avg_upper = frac_nu * (nu_i_meany_cen + nu_i_meany_plus) + frac_nubar * (nubar_i_meany_cen + nubar_i_meany_plus)
ax.fill_between(E, flux_avg_lower, flux_avg_upper, color='purple', alpha=0.5)
ax.plot(E, flux_avg, color='purple', label=r'$\textrm{daemonflux average}$')

# r = n / b 
# Load IceCube Data
# log10(Emin) log10(Emax) log10(Enu) Enu_+err Enu_-err <y> <y>_+err <y>_-err
ic_data = np.loadtxt('data/icecube_inelasticity.txt', skiprows=1)
ic_energy = 10**ic_data[:, 2]

xhi = 10**(ic_data[:, 2] + ic_data[:, 3]) - 10**ic_data[:, 2]
xlo = np.abs(10**(ic_data[:, 2] - ic_data[:, 4]) - 10**ic_data[:, 2])
xerr = np.zeros((2, 5))
xerr[0, :] = xlo
xerr[1, :] = xhi
yhi = ic_data[:, 6]
ylo = ic_data[:, 7]
yerr = np.zeros((2, 5))
yerr[0, :] = ylo
yerr[1, :] = yhi
ax.errorbar(ic_energy, ic_data[:, 5], 
            xerr=xerr,
            yerr=yerr, elinewidth=3, linestyle='none', color='k', label='IceCube (2018)')


ax.set_xscale('log')
ax.set_xlim([1e2, 1e9])
ax.set_ylim([0.0, 0.75])
ax.set_ylabel(r'$\langle y \rangle$')
ax.set_xlabel(r'$E_{\nu}~[\textrm{GeV}]$')
ax.legend()
plt.tight_layout()
plt.savefig('3_mean_inelasticity.png')
plt.savefig('3_mean_inelasticity.pdf')