"""
Check that the splines are not bumpy!
"""

import numpy as np
import photospline
from matplotlib import pyplot as plt

base_dir = '/n/holylfs05/LABS/arguelles_delgado_lab/Everyone/pweigel/sandbox/src/NuXSSplMkr/data/CT18_NLO/replica_0'
base_dir_parton = '/n/holylfs05/LABS/arguelles_delgado_lab/Everyone/pweigel/sandbox/src/NuXSSplMkr/data/CT18_NLO_parton/replica_0'
light_fit_file = '/F2_neutrino_proton_light.fits'
charm_fit_file = '/F2_neutrino_proton_charm.fits'
total_fit_file = '/F2_neutrino_proton_total.fits'

nlo_light = photospline.SplineTable(base_dir + light_fit_file)
nlo_charm = photospline.SplineTable(base_dir + charm_fit_file)
parton = photospline.SplineTable(base_dir_parton + total_fit_file)

# print(light_spline.extents)
fig = plt.figure()
ax = fig.add_subplot(111)

Q2_value = np.log10(1.69)
print(Q2_value)
x_values = np.linspace(-2, 0, 250)

f2_values = nlo_light.evaluate_simple([Q2_value, x_values]) + nlo_charm.evaluate_simple([Q2_value, x_values])
ax.plot(10**x_values, f2_values)

f2_values = parton.evaluate_simple([Q2_value, x_values])
ax.plot(10**x_values, f2_values)
ax.set_xscale('log')
print(min(f2_values))
plt.savefig('1_spline_check.png')
