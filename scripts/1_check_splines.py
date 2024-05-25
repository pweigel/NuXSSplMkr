"""
Check that the splines are not bumpy!
"""

import numpy as np
import photospline
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib import rcParams

rcParams['font.family']           = 'serif'
rcParams['font.serif']            = 'Computer Modern Roman'
rcParams['font.weight']           = 200
rcParams['font.size']             = 14
rcParams['text.usetex']           = True

rcParams['grid.color']            = 'black'
rcParams['grid.alpha']            = 0.10
rcParams['grid.linestyle']        = '-'

rcParams['axes.grid']             = False
rcParams['axes.linewidth']        = 1.5
# rcParams['axes.labelpad']         = 8.0
rcParams['axes.labelsize']        = 14

rcParams['xtick.labelsize']       = 14
rcParams['ytick.labelsize']       = 14
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
rcParams['legend.title_fontsize'] = 14
rcParams['legend.fontsize']       = 14

from matplotlib.ticker import AutoMinorLocator, MultipleLocator

# formats = {
#   'APFEL': {'linestyle': '-'}, 
#   'TMC': {'linestyle': '--'},
#   'CKMT': {'linestyle': '-.'},
#   'PCAC': {'linestyle': ':'},
# }

formats = {
  'TMC': {'linestyle': ':'},
  'CKMT': {'linestyle': '--'},
  'PCAC': {'linestyle': '-'},
}

q2_vals = np.log10([0.1, 1.0, 2.0, 4.0])
x_values = np.linspace(np.log10(5e-3), 1, 250)

base_dir = '/work/submit/pweigel/icecube/src/NuXSSplMkr/data/nCTEQ15_H/replica_0'
sf_types = ['APFEL', 'TMC', 'CKMT', 'PCAC']
# sf_types = ['APFEL', 'TMC', 'CKMT', 'PCAC']
F2s = {}
for sft in sf_types:
    F2s[sft] = {}
    if sft != 'APFEL':
        fn_p = base_dir + '/F2_neutrino_proton_total_'+sft+'.fits'
        fn_n = base_dir + '/F2_neutrino_neutron_total_'+sft+'.fits'
    else:
        fn_p = base_dir + '/F2_neutrino_proton_total.fits'
        fn_n = base_dir + '/F2_neutrino_neutron_total.fits'
    spline_p = photospline.SplineTable(fn_p)
    spline_n = photospline.SplineTable(fn_n)

    for _q2 in q2_vals:
        if 10**_q2 < 2.0 and sft in ['APFEL', 'TMC']:
            continue
        F2s[sft][_q2] = 0.5*(spline_p.evaluate_simple([_q2, x_values]) +  spline_n.evaluate_simple([_q2, x_values]))

F3s = {}
for sft in sf_types:
    F3s[sft] = {}
    if sft != 'APFEL':
        fn_p = base_dir + '/F3_neutrino_proton_total_'+sft+'.fits'
        fn_n = base_dir + '/F3_neutrino_neutron_total_'+sft+'.fits'
    else:
        fn_p = base_dir + '/F3_neutrino_proton_total.fits'
        fn_n = base_dir + '/F3_neutrino_neutron_total.fits'
    spline_p = photospline.SplineTable(fn_p)
    spline_n = photospline.SplineTable(fn_n)

    for _q2 in q2_vals:
        if 10**_q2 < 2.0 and sft in ['APFEL', 'TMC']:
            continue
        F3s[sft][_q2] = 0.5*(spline_p.evaluate_simple([_q2, x_values]) +  spline_n.evaluate_simple([_q2, x_values]))

fig = plt.figure()
ax = fig.add_subplot(111)

ax.plot(10**x_values, F2s['TMC'][np.log10(2)], color='k')
ax.plot(10**x_values, 10**x_values * F3s['TMC'][np.log10(2)], color='k', linestyle='--')
ax.plot(10**x_values, F2s['APFEL'][np.log10(2)], color='gray', linewidth=1)
ax.plot(10**x_values, 10**x_values * F3s['APFEL'][np.log10(2)], color='gray', linestyle='--', linewidth=1)
ax.plot(10**x_values, F2s['CKMT'][np.log10(2)], color='b')
ax.plot(10**x_values, 10**x_values * F3s['CKMT'][np.log10(2)], color='orange', linestyle='--')

ax.set_xlabel(r'$x$')
ax.set_ylim([0, 2.5])
ax.set_ylabel(r'$F_{i}(x,Q^2)$')
ax.set_xscale('log')
ax.set_xlim([5e-3, 1])
# ax.yaxis.tick_right()
ax.yaxis.set_minor_locator(MultipleLocator(0.25/2))
ax.yaxis.set_ticks_position('both')
# ax.legend()
plt.tight_layout()
plt.savefig('ckmt_checks.png', dpi=300, bbox_inches='tight')

fig = plt.figure()
ax = fig.add_subplot(111)

ax.plot(10**x_values, F2s['TMC'][np.log10(4)] / F2s['APFEL'][np.log10(4)], color='k')
ax.plot(10**x_values, F3s['TMC'][np.log10(4)] / F3s['APFEL'][np.log10(4)], color='k', linestyle='--')

ax.set_xlabel(r'$x$')
ax.set_ylim([0.5, 2.5])
ax.set_ylabel(r'Ratio TMC/APFEL')
# ax.set_xscale('log')
ax.set_xlim([5e-3, 1])
# ax.yaxis.tick_right()
ax.yaxis.set_minor_locator(MultipleLocator(0.25/2))
ax.yaxis.set_ticks_position('both')
# ax.legend()
plt.tight_layout()
plt.savefig('tmc_ratio.png', dpi=300, bbox_inches='tight')
