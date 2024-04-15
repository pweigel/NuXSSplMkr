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
        
    return np.array(energy)/1e9, 10**np.array(log_xs)
xs_path = '/n/holylfs05/LABS/arguelles_delgado_lab/Everyone/pweigel/sandbox/src/NuXSSplMkr/data/CT18A_NNLO/replica_0/cross_sections/'
# xs_path = '/n/holylfs05/LABS/arguelles_delgado_lab/Everyone/pweigel/sandbox/src/NuXSSplMkr/data/EPPS21_NLO_O16/replica_0/cross_sections/'

energy, total = load_cross_section(xs_path+'total_neutrino_proton_total.out')
_, light = load_cross_section(xs_path+'total_neutrino_proton_light.out')
_, charm = load_cross_section(xs_path+'total_neutrino_proton_charm.out')
_, bottom = load_cross_section(xs_path+'total_neutrino_proton_bottom.out')
energy_top, top = load_cross_section(xs_path+'total_neutrino_proton_top.out')

total_new = light + charm + bottom
total_new[(len(total)-len(top)):] += top

fig = plt.figure(figsize=(12,9))
ax = fig.add_subplot(111)
ax.plot(energy, total, linewidth=3)
ax.plot(energy, total_new, linewidth=1, color='k')
ax.plot(energy, light, linewidth=3)
ax.plot(energy, charm, linewidth=3)
ax.plot(energy, bottom, linewidth=3)
ax.plot(energy_top, top, linewidth=3)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim([1e2, 1e9])
ax.set_ylim([1e-38, 2e-32])
ax.set_xlabel(r'Energy [GeV]')
ax.set_ylabel(r'$\sigma~[\textrm{cm}^{2}]$')

plt.tight_layout()
plt.savefig('2_cross_sections.png')

fig = plt.figure(figsize=(12,9))
ax = fig.add_subplot(111)
ax.plot(energy, total/total, linewidth=3)
ax.plot(energy, light/total, linewidth=3)
ax.plot(energy, charm/total, linewidth=3)
ax.plot(energy, bottom/total, linewidth=3)
ax.plot(energy_top, top/total[(len(total)-len(top)):], linewidth=3)

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_ylim([1e-6, 1])
ax.set_xlim([1e2, 1e9])
ax.set_xlabel(r'Energy [GeV]')
ax.set_ylabel(r'$\sigma~[\textrm{cm}^{2}]$')
plt.tight_layout()
plt.savefig('2_cross_sections_fraction.png')
