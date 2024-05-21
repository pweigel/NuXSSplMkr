import numpy as np
from matplotlib import pyplot as plt
import matplotlib as mpl
path = '/n/holylfs05/LABS/arguelles_delgado_lab/Everyone/pweigel/sandbox/src/NuXSSplMkr/data/nCTEQ15_H_PCAC/replica_0'
intest = np.loadtxt(path+'/F1_neutrino_proton_total.grid', skiprows=2, delimiter=',')

with open(path+'/F1_neutrino_proton_total.grid', 'r') as f:
    lines = f.readlines()
    num = [int(x) for x in lines[0].rstrip('\n').split(' ')]
    bounds = [float(x) for x in lines[1].rstrip('\n').split(' ')]

_x = np.logspace(bounds[0], bounds[1], num[0])
_y = np.logspace(bounds[2], bounds[3], num[1])

fig = plt.figure()
ax = fig.add_subplot(111)
x, y = np.meshgrid(_x, _y)
ax_img = ax.pcolormesh(x, y, intest.T, norm=mpl.colors.LogNorm())

ax.set_xscale('log')
ax.set_yscale('log')

ax.set_xlim([0.01, 1e4])

plt.colorbar(ax_img, ax=ax)
plt.savefig('10_sf_test.png', dpi=300)