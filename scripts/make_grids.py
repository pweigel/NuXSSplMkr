"""
Make the grids used for SF and XS calculations
"""
import numpy as np
from matplotlib import pyplot as plt




fig = plt.figure()
ax = fig.add_subplot(111)

Npts = 50
x_values = np.logspace(-9, -1, 50)
q_values = np.logspace(1, 10, 50)
for x in x_values:
    for q in q_values:
        ax.scatter(x, q, color='k', s=1)

x_values = np.linspace(0.1, 1, 30)
q_values = np.logspace(1, 10, 50)
for x in x_values:
    for q in q_values:
        ax.scatter(x, q, color='g', s=1)

x_values = np.linspace(0.1, 1, 30)
q_values = np.linspace(1, 10, 30)
for x in x_values:
    for q in q_values:
        ax.scatter(x, q, color='m', s=1)

x_values = np.logspace(-9, -1, 50)
q_values = np.linspace(1, 10, 30)
for x in x_values:
    for q in q_values:
        ax.scatter(x, q, color='r', s=1)

ax.set_xlim([1e-9, 1])
ax.set_ylim([1e-1, 1e10])
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('x')
ax.set_ylabel(r'$Q^2$')
plt.savefig('grid_test.png', dpi=300)