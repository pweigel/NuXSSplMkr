"""
Demonstration of the phase space constraints
"""
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


E = 5e2
m_N = 0.938
s = 2 * m_N * E

x_values = np.logspace(-4, -1e-9, 1001)
W2_min = 4.0
# Q2 = s x y
# W2 = Q2 * (1. - x) / x
# Q2 = W2min * x / (1 - x)
fig = plt.figure(figsize=(12,9))
ax = fig.add_subplot(111)


ax.text(1e-2, 1.5, s=r'$Q^2 = Q_{min}^2$', fontsize=20)
ax.plot([x_values[0], x_values[-1]], [1.3, 1.3], color='k')

ax.text(5e-2, W2_min * 5e-2 / (1 - 5e-2) + 7e-2, s=r'$W^2 = W_{min}^2$', rotation=32.5, fontsize=20)
ax.plot(x_values, W2_min * x_values / (1 - x_values), color='k')

ax.text(1e-2, s * 1e-2 + 4, s=r'$Q^2 = s$', rotation=30, fontsize=20)
ax.plot(x_values, s * x_values, color='k')

ax.plot([x_values[0], x_values[-1]], [4, 4], color='gray', linestyle='--')

ax.set_xlabel(r'$x$')
ax.set_ylabel(r'$Q^2~[\textrm{GeV}]$')
ax.set_xscale('log')
ax.set_yscale('log')

ax.set_xlim([1e-3, 1])
ax.set_ylim([1e-1, 1e3])

plt.tight_layout()
plt.savefig('3_phase_space.png')
