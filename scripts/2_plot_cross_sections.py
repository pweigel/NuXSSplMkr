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

def load_cross_section(fname):
    energy, log_xs = [], []
    with open(fname, 'r') as f:
        for line in f.readlines():
            x, y = line.rstrip('\n').split(',')
            if float(y) > -60:
                energy.append(float(x))
                log_xs.append(float(y))
                
    return np.array(energy)/1e9, np.array(log_xs)
# xs_path = '/n/holylfs05/LABS/arguelles_delgado_lab/Everyone/pweigel/sandbox/src/NuXSSplMkr/reference_data/CT18A_NNLO/replica_0/cross_sections/'
# xs_path = '/n/holylfs05/LABS/arguelles_delgado_lab/Everyone/pweigel/sandbox/src/NuXSSplMkr/data/EPPS21_NLO_O16/replica_0/cross_sections/'
xs_path = '/data/submit/pweigel/cross_sections/CT18A_NNLO/replica_0/cross_sections/'

fig = plt.figure(figsize=(12,9))
ax = fig.add_subplot(111)
energy, total_p = load_cross_section(xs_path+'total_neutrino_proton_light.out')
energy, total_n = load_cross_section(xs_path+'total_neutrino_neutron_light.out')
energy, total_pc = load_cross_section(xs_path+'total_neutrino_proton_charm.out')
energy, total_nc = load_cross_section(xs_path+'total_neutrino_neutron_charm.out')
total = 0.5*(total_p + total_n + total_pc + total_nc)
ax.plot(energy, total * 1e38 / energy, linewidth=2, label=r'$\nu N$')
energy, total_p = load_cross_section(xs_path+'total_antineutrino_proton_light.out')
energy, total_n = load_cross_section(xs_path+'total_antineutrino_neutron_light.out')
energy, total_pc = load_cross_section(xs_path+'total_antineutrino_proton_charm.out')
energy, total_nc = load_cross_section(xs_path+'total_antineutrino_neutron_charm.out')
total = 0.5*(total_p + total_n + total_pc + total_nc)
ax.plot(energy, total * 1e38 / energy, linewidth=2, label=r'$\bar{\nu} N$')
ax.plot([10, 1000], [0.675, 0.675], color='k', alpha=0.33)
ax.plot([10, 1000], [0.33, 0.33], color='k', linestyle='--', alpha=0.33)
ax.set_xscale('log')
ax.set_yscale('linear')
ax.set_xlim([10, 1000])
ax.set_ylim([0, 1.3])
ax.set_xlabel('Energy [GeV]')
ax.set_ylabel(r'$\sigma^{CC}/A/E_{\nu}~[10^{-38}~\rm{cm^2}/\rm{GeV}]$')
plt.tight_layout()
plt.savefig('2_xs_test_ct18annlo_lowe.pdf')


fig = plt.figure(figsize=(12,9))
ax = fig.add_subplot(111)
energy, total_p = load_cross_section(xs_path+'total_neutrino_proton_light.out')
energy, total_n = load_cross_section(xs_path+'total_neutrino_neutron_light.out')
energy, total_pc = load_cross_section(xs_path+'total_neutrino_proton_charm.out')
energy, total_nc = load_cross_section(xs_path+'total_neutrino_neutron_charm.out')
total = 0.5*(total_p + total_n + total_pc + total_nc)
ax.plot(energy, total * 1e38, linewidth=2, label=r'$\nu N$')
energy, total_p = load_cross_section(xs_path+'total_antineutrino_proton_light.out')
energy, total_n = load_cross_section(xs_path+'total_antineutrino_neutron_light.out')
energy, total_pc = load_cross_section(xs_path+'total_antineutrino_proton_charm.out')
energy, total_nc = load_cross_section(xs_path+'total_antineutrino_neutron_charm.out')
total = 0.5*(total_p + total_n + total_pc + total_nc)
ax.plot(energy, total * 1e38, linewidth=2, label=r'$\bar{\nu} N$')
# ax.plot([10, 1000], [0.675, 0.675], color='k', alpha=0.33)
# ax.plot([10, 1000], [0.33, 0.33], color='k', linestyle='--', alpha=0.33)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim([1e1, 1e9])
# ax.set_ylim([0, 1.3])
ax.set_xlabel('Energy [GeV]')
ax.set_ylabel(r'$\sigma^{CC}/A/E_{\nu}~[10^{-38}~\rm{cm^2}/\rm{GeV}]$')
plt.tight_layout()
plt.savefig('2_xs_test_ct18annlo.pdf')