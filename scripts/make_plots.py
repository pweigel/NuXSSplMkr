import numpy as np
from matplotlib import pyplot as plt
import matplotlib
import numpy as np
from matplotlib import pyplot as plt
import matplotlib as mpl
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

def load_gary_data(fname):
    E_nu, E_nubar, frac_nu, frac_nubar = [], [], [], []
    with open(fname, 'r') as f:
        for line in f.readlines():
            z = [float(x) for x in line.rstrip('\n').split(',')]
            E_nu.append(z[0])
            frac_nu.append(z[1])
            E_nubar.append(z[2])
            frac_nubar.append(z[3])
    return np.array(E_nu), np.array(frac_nu), np.array(E_nubar), np.array(frac_nubar)

def load_dsdy(fname):
    NE = 100
    Ny = 100
    energy_vals = np.logspace(2, 9, NE)
    # y_vals = np.linspace(1e-4, 1.0, Ny)
    y_vals = np.logspace(-6, 0, Ny)
    
    xs_values = np.zeros((NE, Ny))
    with open(fname, 'r') as f:
        n = 0
        for line in f.readlines():
            _xs = float(line.rstrip('\n'))
            xs_values[n // Ny, n % Ny] = 10**_xs
            n += 1
    return energy_vals, y_vals, xs_values

def load_xs(fname):
    E = []
    xs = []
    with open(fname, 'r') as f:
        for l in f.readlines():
            d = [float(x) for x in l.rstrip('\n').split(',') ]
            E.append(d[0] / 1e9)
            xs.append(10**d[1])
    return np.array(E), np.array(xs)

def collect_cross_sections(name):
    base_path = '../data/{}/cross_sections/'.format(name)
    
    xs = {}
    xs['nu'] = {'p': {}, 'n': {}, 'i': {}}
    xs['nubar'] = {'p': {}, 'n': {}, 'i': {}}
    
    E, nu_p_total = load_xs(base_path + 'total_neutrino_proton_total.out')
    _, nu_n_total = load_xs(base_path + 'total_neutrino_neutron_total.out')
    _, nubar_p_total = load_xs(base_path + 'total_antineutrino_proton_total.out')
    _, nubar_n_total = load_xs(base_path + 'total_antineutrino_neutron_total.out')
    _, nu_p_charm = load_xs(base_path + 'total_neutrino_proton_charm.out')
    _, nu_n_charm = load_xs(base_path + 'total_neutrino_neutron_charm.out')
    _, nubar_p_charm = load_xs(base_path + 'total_antineutrino_proton_charm.out')
    _, nubar_n_charm = load_xs(base_path + 'total_antineutrino_neutron_charm.out')
    
    xs['E'] = E
    
    xs['nu']['p']['total'] = nu_p_total
    xs['nu']['p']['charm'] = nu_p_charm
    xs['nu']['n']['total'] = nu_n_total
    xs['nu']['n']['charm'] = nu_n_charm
    xs['nubar']['p']['total'] = nubar_p_total
    xs['nubar']['p']['charm'] = nubar_p_charm
    xs['nubar']['n']['total'] = nubar_n_total
    xs['nubar']['n']['charm'] = nubar_n_charm
    xs['nu']['i']['total'] = 0.5 * (nu_p_total + nu_n_total)
    xs['nu']['i']['charm'] = 0.5 * (nu_p_charm + nu_n_charm)
    xs['nubar']['i']['total'] = 0.5 * (nubar_p_total + nubar_n_total)
    xs['nubar']['i']['charm'] = 0.5 * (nubar_p_charm + nubar_n_charm)
    
    dxs = {}
    dxs['nu'] = {'p': {}, 'n': {}, 'i': {}}
    dxs['nubar'] = {'p': {}, 'n': {}, 'i': {}}
    
    E, y, nu_p_total = load_dsdy(base_path + 'dsdy_neutrino_proton_total.out')
    _, _, nu_n_total = load_dsdy(base_path + 'dsdy_neutrino_neutron_total.out')
    _, _, nubar_p_total = load_dsdy(base_path + 'dsdy_antineutrino_proton_total.out')
    _, _, nubar_n_total = load_dsdy(base_path + 'dsdy_antineutrino_neutron_total.out')
    _, _, nu_p_charm = load_dsdy(base_path + 'dsdy_neutrino_proton_charm.out')
    _, _, nu_n_charm = load_dsdy(base_path + 'dsdy_neutrino_neutron_charm.out')
    _, _, nubar_p_charm = load_dsdy(base_path + 'dsdy_antineutrino_proton_charm.out')
    _, _, nubar_n_charm = load_dsdy(base_path + 'dsdy_antineutrino_neutron_charm.out')
    
    dxs['E'] = E
    dxs['y'] = y
    
    dxs['nu']['p']['total'] = nu_p_total
    dxs['nu']['p']['charm'] = nu_p_charm
    dxs['nu']['n']['total'] = nu_n_total
    dxs['nu']['n']['charm'] = nu_n_charm
    dxs['nubar']['p']['total'] = nubar_p_total
    dxs['nubar']['p']['charm'] = nubar_p_charm
    dxs['nubar']['n']['total'] = nubar_n_total
    dxs['nubar']['n']['charm'] = nubar_n_charm
    dxs['nu']['i']['total'] = 0.5 * (nu_p_total + nu_n_total)
    dxs['nu']['i']['charm'] = 0.5 * (nu_p_charm + nu_n_charm)
    dxs['nubar']['i']['total'] = 0.5 * (nubar_p_total + nubar_n_total)
    dxs['nubar']['i']['charm'] = 0.5 * (nubar_p_charm + nubar_n_charm)
    
    return xs, dxs

CSMS_E =  np.array([50,   100,  200, 500, 1000, 2000, 5000, 1e4, 2e4, 5e4, 1e5, 2e5, 5e5, 1e6,  2e6,  5e6,  1e7,  2e7,  5e7,  1e8,  2e8,   5e8,  1e9])
CSMS_nu_i = 1e-36 * np.array([0.32, 0.65, 1.3, 3.2,  6.2,   12,    27, 47,  77, 140, 210, 310, 490, 690,  950, 1400, 1900, 2600, 3700, 4800, 6200, 8700, 11000])
CSMS_nubar_i = 1e-36 * np.array([0.15, 0.33, 0.69, 1.8, 3.6, 7, 17, 31, 55, 110, 180, 270, 460, 660, 920, 1400, 1900, 2600, 3700, 4800, 6200, 8700, 11000])

# Load gary data
gary_nu_E, gary_nu_charm_fraction, gary_nubar_E, gary_nubar_charm_fraction = load_gary_data('gary_data/gary_charm_massive.csv')
# gary_nu_E, gary_nu_charm_fraction, gary_nubar_E, gary_nubar_charm_fraction = load_gary_data('gary_data/gary_charm_massless.csv')

xs, dxs = collect_cross_sections('CT18ANNLO_FONLL-C_pto2')
xs_nocharm, dxs_nocharm = collect_cross_sections('CT18ANNLO_nocharm_FONLL-C_pto2')

plot_path = 'plots/'

# Total Cross Section Plots
fig = plt.figure(figsize=(12, 9 ))
ax = fig.add_subplot(111)
ax.plot(xs['E'], xs['nu']['i']['total'], color='k', linewidth=2, label=r'$\nu~\textrm{Total}$')
ax.plot(xs['E'], xs['nubar']['i']['total'], color='k', linewidth=2, linestyle='--', label=r'$\bar{\nu}~\textrm{Total}$')
ax.plot(xs['E'], xs_nocharm['nu']['i']['charm'], color='m', linewidth=2, label=r'$\nu~\textrm{Charm}$')
ax.plot(xs['E'], xs_nocharm['nubar']['i']['charm'], color='m', linewidth=2, linestyle='--', label=r'$\bar{\nu}~\textrm{Charm}$')
ax.plot(CSMS_E, CSMS_nu_i, color='b', linewidth=1, linestyle='-', label=r'$\textrm{CSMS}~\nu$')
ax.plot(CSMS_E, CSMS_nubar_i, color='b', linewidth=1, linestyle='--', label=r'$\textrm{CSMS}~\bar{\nu}$')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim([xs['E'][0], xs['E'][-1]])
ax.set_xlabel(r'$\textrm{Energy}~[\textrm{GeV}]$')
ax.set_ylabel(r'$\sigma~[\textrm{cm}^2]$')
ax.set_title(r'$\textrm{CT18A NNLO}$')
ax.legend()
ax.grid(True, which='Both')
plt.tight_layout()
plt.savefig('plots/total_cross_section.pdf')
plt.close()

fig = plt.figure(figsize=(12, 9 ))
ax = fig.add_subplot(111)
ax.plot(xs['E'], xs['nu']['i']['total'] / xs['E'], color='k', linewidth=2, label=r'$\nu~\textrm{Total}$')
ax.plot(xs['E'], xs['nubar']['i']['total'] / xs['E'], color='k', linewidth=2, linestyle='--', label=r'$\bar{\nu}~\textrm{Total}$')
ax.plot(xs['E'], xs_nocharm['nu']['i']['charm'] / xs['E'], color='m', linewidth=2, label=r'$\nu~\textrm{Charm}$')
ax.plot(xs['E'], xs_nocharm['nubar']['i']['charm'] / xs['E'], color='m', linewidth=2, linestyle='--', label=r'$\bar{\nu}~\textrm{Charm}$')
ax.plot(CSMS_E, CSMS_nu_i / CSMS_E, color='b', linewidth=1, linestyle='-', label=r'$\textrm{CSMS}~\nu$')
ax.plot(CSMS_E, CSMS_nubar_i / CSMS_E, color='b', linewidth=1, linestyle='--', label=r'$\textrm{CSMS}~\bar{\nu}$')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim([xs['E'][0], xs['E'][-1]])
ax.set_xlabel(r'$\textrm{Energy}~[\textrm{GeV}]$')
ax.set_ylabel(r'$\sigma/E~[\textrm{cm}^2 / \textrm{GeV}]$')
ax.set_title(r'$\textrm{CT18A NNLO}$')
ax.legend()
ax.grid(True, which='Both')
plt.tight_layout()
plt.savefig('plots/total_cross_section_over_E.pdf')
plt.close()

# Charm Fraction
fig = plt.figure(figsize=(12, 9 ))
ax = fig.add_subplot(111)
ax.plot(xs['E'], xs['nu']['i']['charm'] / xs['nu']['i']['total'], linewidth=2, color='g', label=r'$\nu~\textrm{w/ c pdf}$')
ax.plot(xs['E'], xs['nubar']['i']['charm'] / xs['nubar']['i']['total'], linewidth=2, color='g', linestyle='--', label=r'$\bar{\nu}~\textrm{w/ c pdf}$')
ax.plot(xs['E'], xs_nocharm['nu']['i']['charm'] / xs['nu']['i']['total'], linewidth=2, color='orange', label=r'$\nu~\textrm{w/o c pdf}$')
ax.plot(xs['E'], xs_nocharm['nubar']['i']['charm'] / xs['nubar']['i']['total'], linewidth=2, color='orange', linestyle='--', label=r'$\bar{\nu}~\textrm{w/o c pdf}$')
ax.plot(gary_nu_E, gary_nu_charm_fraction, color='k', linestyle='-', label=r'$\textrm{Gary}~\nu$')
ax.plot(gary_nubar_E, gary_nubar_charm_fraction, color='k', linestyle='--', label=r'$\textrm{Gary}~\bar{\nu}$')

ax.set_xscale('log')
ax.set_yscale('linear')
ax.set_xlim([xs['E'][0], xs['E'][-1]])
ax.set_ylim([0, 0.5])
ax.set_xlabel(r'$\textrm{Energy}~[\textrm{GeV}]$')
ax.set_ylabel(r'$\textrm{Charm Fraction}$')
ax.set_title(r'$\textrm{CT18A NNLO}$')
ax.legend()
ax.grid(True, which='Both')
plt.tight_layout()
plt.savefig('plots/charm_fraction.pdf')
plt.close()

# Single differential cross sections


Ei = 0
scale = 1e35
fig = plt.figure(figsize=(12, 9))
ax = fig.add_subplot(111)
ax.plot(dxs['y'], scale*dxs['nu']['i']['total'][Ei, :], linewidth=2, color='k', label=r'$\nu~\textrm{Total}$')
ax.plot(dxs['y'], scale*dxs['nubar']['i']['total'][Ei, :], linewidth=2, color='k', linestyle='--', label=r'$\bar{\nu}~\textrm{Total}$')
ax.plot(dxs['y'], scale*dxs_nocharm['nu']['i']['charm'][Ei, :], linewidth=2, color='m', label=r'$\nu~\textrm{Charm}$')
ax.plot(dxs['y'], scale*dxs_nocharm['nubar']['i']['charm'][Ei, :], linewidth=2, color='m', linestyle='--', label=r'$\bar{\nu}~\textrm{Charm}$')
ax.set_xscale('linear')
ax.set_yscale('linear')
ax.set_xlim([0, 1])
ax.set_ylim([0, None])
ax.set_xlabel(r'$y$')
ax.set_ylabel(r'$\frac{d\sigma}{dy}~[10^{-35}~\textrm{cm}^2]$')
ax.set_title(r'$\textrm{CT18A NNLO},~E = 500~\textrm{GeV}$')
ax.legend()
ax.grid(True, which='both')
plt.tight_layout()
plt.savefig('plots/dsdy_500GeV.pdf')
plt.close()

scale = 1e35
fig = plt.figure(figsize=(12, 9))
ax = fig.add_subplot(111)
ax.plot(dxs['y'], dxs['nu']['i']['total'][Ei, :] / xs['nu']['i']['total'][Ei], linewidth=2, color='k', label=r'$\nu~\textrm{Total}$')
ax.plot(dxs['y'], dxs['nubar']['i']['total'][Ei, :] / xs['nubar']['i']['total'][Ei], linewidth=2, color='k', linestyle='--', label=r'$\bar{\nu}~\textrm{Total}$')
ax.plot(dxs['y'], dxs_nocharm['nu']['i']['charm'][Ei, :] / xs['nu']['i']['total'][Ei], linewidth=2, color='m', label=r'$\nu~\textrm{Charm}$')
ax.plot(dxs['y'], dxs_nocharm['nubar']['i']['charm'][Ei, :] / xs['nubar']['i']['total'][Ei], linewidth=2, color='m', linestyle='--', label=r'$\bar{\nu}~\textrm{Charm}$')
ax.set_xscale('linear')
ax.set_yscale('linear')
ax.set_xlim([0, 1])
ax.set_ylim([0, None])
ax.set_xlabel(r'$y$')
ax.set_ylabel(r'$\frac{1}{\sigma}\frac{d\sigma}{dy}$')
ax.set_title(r'$\textrm{CT18A NNLO},~E = 500~\textrm{GeV}$')
ax.legend()
ax.grid(True, which='both')
plt.tight_layout()
plt.savefig('plots/dsdy_500GeV_normalized.pdf')
plt.close()

Ei = 25
scale = 1e35
fig = plt.figure(figsize=(12, 9))
ax = fig.add_subplot(111)
ax.plot(dxs['y'], scale*dxs['nu']['i']['total'][Ei, :], linewidth=2, color='k', label=r'$\nu~\textrm{Total}$')
ax.plot(dxs['y'], scale*dxs['nubar']['i']['total'][Ei, :], linewidth=2, color='k', linestyle='--', label=r'$\bar{\nu}~\textrm{Total}$')
ax.plot(dxs['y'], scale*dxs_nocharm['nu']['i']['charm'][Ei, :], linewidth=2, color='m', label=r'$\nu~\textrm{Charm}$')
ax.plot(dxs['y'], scale*dxs_nocharm['nubar']['i']['charm'][Ei, :], linewidth=2, color='m', linestyle='--', label=r'$\bar{\nu}~\textrm{Charm}$')
ax.set_xscale('linear')
ax.set_yscale('linear')
ax.set_xlim([0, 1])
ax.set_ylim([0, None])
ax.set_xlabel(r'$y$')
ax.set_ylabel(r'$\frac{d\sigma}{dy}~[10^{-35}~\textrm{cm}^2]$')
ax.set_title(r'$\textrm{CT18A NNLO},~E = 5857~\textrm{GeV}$')
ax.legend()
ax.grid(True, which='both')
plt.tight_layout()
plt.savefig('plots/dsdy_5857GeV.pdf')
plt.close()

scale = 1e35
fig = plt.figure(figsize=(12, 9))
ax = fig.add_subplot(111)
ax.plot(dxs['y'], dxs['nu']['i']['total'][Ei, :] / xs['nu']['i']['total'][Ei], linewidth=2, color='k', label=r'$\nu~\textrm{Total}$')
ax.plot(dxs['y'], dxs['nubar']['i']['total'][Ei, :] / xs['nubar']['i']['total'][Ei], linewidth=2, color='k', linestyle='--', label=r'$\bar{\nu}~\textrm{Total}$')
ax.plot(dxs['y'], dxs_nocharm['nu']['i']['charm'][Ei, :] / xs['nu']['i']['total'][Ei], linewidth=2, color='m', label=r'$\nu~\textrm{Charm}$')
ax.plot(dxs['y'], dxs_nocharm['nubar']['i']['charm'][Ei, :] / xs['nubar']['i']['total'][Ei], linewidth=2, color='m', linestyle='--', label=r'$\bar{\nu}~\textrm{Charm}$')
ax.set_xscale('linear')
ax.set_yscale('linear')
ax.set_xlim([0, 1])
ax.set_ylim([0, None])
ax.set_xlabel(r'$y$')
ax.set_ylabel(r'$\frac{1}{\sigma}\frac{d\sigma}{dy}$')
ax.set_title(r'$\textrm{CT18A NNLO},~E = 5857~\textrm{GeV}$')
ax.legend()
ax.grid(True, which='both')
plt.tight_layout()
plt.savefig('plots/dsdy_5857GeV_normalized.pdf')
plt.close()

exit()
scale = 1e36
_x, _y = np.meshgrid(E, y)
fig = plt.figure(figsize=(12, 9))
ax = fig.add_subplot(111)
ax_img = ax.pcolormesh(_x, _y, nu_i_total.T, norm=mpl.colors.LogNorm(vmin=1e-38, vmax=1e-30), cmap='Spectral_r')
ax.set_xscale('log')
ax.set_yscale('log')
plt.colorbar(ax_img, ax=ax)
plt.savefig('test_dsdy_nocharm.png')


Ei = 0
fig = plt.figure(figsize=(12, 9))
ax = fig.add_subplot(111)
ax.plot(y, scale * nu_i_total[Ei, :], color='r', label=r'$\nu$ Total')
ax.plot(y, scale * nubar_i_total[Ei, :], color='b', label=r'$\bar{\nu}$ Total')
ax.plot(y, scale * (nu_i_total[Ei, :] + nubar_i_total[Ei, :]), color='k', label=r'$\nu+\bar{\nu}$ Total')
ax.plot(y, scale * nu_i_charm[Ei, :], color='r', linestyle='--', label=r'$\nu$ Charm')
ax.plot(y, scale * nubar_i_charm[Ei, :], color='b', linestyle='--', label=r'$\bar{\nu}$ Charm')
ax.plot(y, scale * (nu_i_charm[Ei, :] + nubar_i_charm[Ei, :]), color='k', linestyle='--', label=r'$\nu+\bar{\nu}$ Charm')

ax.set_title(r'$E = 100~\textrm{GeV}$')
plt.savefig('ok.pdf')
print(E[Ei])

# fig = plt.figure(figsize=(8, 6))
# gs = fig.add_gridspec(3, 1, height_ratios=(3, 1, 1), wspace=0.05, hspace=0.125,
#                             top=0.95, bottom=0.1, left=0.1, right=0.9)

# ax = fig.add_subplot(gs[0])

# ax.plot(y, csms_dsdy_iso_total * 1e35, color='g', linestyle='--', label='CSMS Total')
# ax.plot(y, csms_dsdy_iso_charm * 1e35, color='m', linestyle='--', label='CSMS Charm')

# ax.plot(y, ct18_dsdy_iso_total * 1e35, color='g', label='CT18A NNLO Total')
# ax.plot(y, ct18_dsdy_iso_charm * 1e35, color='m', label='CT18A NNLO Charm')
# ax.tick_params(right=True, top=True, which='both', direction='in')
# ax.set_xlim([0, 1])
# ax.set_ylim([0, 0.4])

# ax.set_xlabel(r'$y$', fontsize=16)
# ax.set_ylabel(r'$\frac{d\sigma}{dy}~[10^{-35}~\textrm{cm}^2]$', fontsize=16)
# ax.legend(loc='upper right', prop={'size': 8})

# ax = fig.add_subplot(gs[1, 0])
# ax.plot(y, csms_dsdy_iso_charm / csms_dsdy_iso_total, color='r', label='CSMS')
# ax.plot(y, ct18_dsdy_iso_charm / ct18_dsdy_iso_total, color='b', label='CT18A NNLO')
# ax.set_xlim([0, 1])
# ax.set_ylim([0, 0.3])
# ax.set_ylabel(r'$\textrm{Charm Fraction}$', fontsize=10)
# ax.legend(prop={'size': 8})
# ax.tick_params(right=True, top=True, which='both', direction='in')

# ax = fig.add_subplot(gs[2, 0])
# ax.plot(y, csms_dsdy_iso_total / ct18_dsdy_iso_total, color='g', label='Total')
# ax.plot(y, csms_dsdy_iso_charm / ct18_dsdy_iso_charm, color='m', label='Charm')
# ax.set_xlim([0, 1])
# ax.set_ylim([1, 1.5])
# ax.set_ylabel(r'\textrm{CSMS/CT18 NNLO}', fontsize=10)
# ax.set_xlabel(r'$y$')
# ax.legend(prop={'size': 8})
# ax.tick_params(right=True, top=True, which='both', direction='in')

# plt.savefig('test_dsdy.pdf')


# base_path = '../data/CT18ANNLO_FONLL-C_pto2/cross_sections/'
base_path = '../data/CT18ANNLO_FONLL-C_pto2/cross_sections/'
E, nu_p_total = load_xs(base_path + 'total_neutrino_proton_total.out')
_, nu_n_total = load_xs(base_path + 'total_neutrino_neutron_total.out')
_, nubar_p_total = load_xs(base_path + 'total_antineutrino_proton_total.out')
_, nubar_n_total = load_xs(base_path + 'total_antineutrino_neutron_total.out')
_, nu_p_charm = load_xs(base_path + 'total_neutrino_proton_charm.out')
_, nu_n_charm = load_xs(base_path + 'total_neutrino_neutron_charm.out')
_, nubar_p_charm = load_xs(base_path + 'total_antineutrino_proton_charm.out')
_, nubar_n_charm = load_xs(base_path + 'total_antineutrino_neutron_charm.out')

base_path = '../data/CT18ANNLO_nocharm_FONLL-C_pto2/cross_sections/'
E, no_charm_nu_p_total = load_xs(base_path + 'total_neutrino_proton_total.out')
_, no_charm_nu_n_total = load_xs(base_path + 'total_neutrino_neutron_total.out')
_, no_charm_nubar_p_total = load_xs(base_path + 'total_antineutrino_proton_total.out')
_, no_charm_nubar_n_total = load_xs(base_path + 'total_antineutrino_neutron_total.out')
_, no_charm_nu_p_charm = load_xs(base_path + 'total_neutrino_proton_charm.out')
_, no_charm_nu_n_charm = load_xs(base_path + 'total_neutrino_neutron_charm.out')
_, no_charm_nubar_p_charm = load_xs(base_path + 'total_antineutrino_proton_charm.out')
_, no_charm_nubar_n_charm = load_xs(base_path + 'total_antineutrino_neutron_charm.out')