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

### EM CKMT Parameters ###
# CKMT_A = 0.1502
# CKMT_B = 1.2064
# CKMT_f = 0.15

## nu F2
CKMT_A = 0.5967
CKMT_B = 2.7145
CKMT_f = 0.5962

CKMT_a = 0.2631
CKMT_b = 0.6452
CKMT_c = 3.5489
CKMT_d = 1.1170
CKMT_Delta0 = 0.07684
CKMT_AlphaR = 0.4150

### CKMT Functions ###

def n_CKMT(Q2):
    return 1.5 * (1.0 + Q2 / (Q2 + CKMT_c))
  
def Delta_CKMT(Q2):
    return CKMT_Delta0 * (1.0 + 2.0 * Q2 / (Q2 + CKMT_d))

## nu F2
CKMT_A = 0.5967
CKMT_B = 2.7145
CKMT_f = 0.5962
def CKMT_SF(x, Q2):
    _delta = Delta_CKMT(Q2)
    _n = n_CKMT(Q2)
    term1 = CKMT_A * x**(-_delta) * (1.0 - x)**(_n + 4.0)
    term2 = (Q2 / (Q2 + CKMT_a))**(1+_delta)
    
    term3 = CKMT_B * x**(1.0 - CKMT_AlphaR) * (1.0 - x)**(_n)
    term4 = (Q2 / (Q2 + CKMT_b))**(CKMT_AlphaR) * (1.0 + CKMT_f*(1.0 - x))

    return term1*term2 + term3*term4

b1,b2,b3 = 0.635,0.5747,-0.3534
lmd=0.2
def CKMT_R(x, Q2):
    big_theta = 1.0 + 12 * (Q2 / (Q2 + 1)) * (0.125**2 / (1.25**2 + x**2))
    R = b1 / np.log10(Q2/lmd**2) * big_theta + b2 / Q2 + b3 / (Q2**2 + 0.3**2)
    return R
  
def CKMT_F1(x, Q2):
    _f2 = CKMT_SF(x, Q2)
    _r = CKMT_R(x, Q2)
    return _f2 * (1 + 4*(0.938**2 * x / Q2)) / (2 * x * (_r + 1))
    
    
fig = plt.figure(figsize=(12,9))
ax = fig.add_subplot(111)

Q2 = 2.0
x_values = np.logspace(np.log10(5e-3), np.log10(5e-1), 101)

f2 = CKMT_SF(x_values, Q2)
f1 = CKMT_F1(x_values, Q2)

## nu xF3
CKMT_A = 9.3995e-3
CKMT_B = 2.4677
CKMT_f = 0.5962
def CKMT_SF(x, Q2):
    _delta = Delta_CKMT(Q2)
    _n = n_CKMT(Q2)
    term1 = CKMT_A * x**(-_delta) * (1.0 - x)**(_n + 4.0)
    term2 = (Q2 / (Q2 + CKMT_a))**(1+_delta)
    
    term3 = CKMT_B * x**(1.0 - CKMT_AlphaR) * (1.0 - x)**(_n)
    term4 = (Q2 / (Q2 + CKMT_b))**(CKMT_AlphaR) * (1.0 + CKMT_f*(1.0 - x))

    return term1*term2 + term3*term4

f3 = CKMT_SF(x_values, Q2)

ax.plot(x_values,2* x_values * f1, linewidth=3, label='2xF1')
# ax.plot(x_values, f2 / (2*x_values), linewidth=3, label='F2/2x')
ax.plot(x_values, f2, linewidth=3, label='F2')
ax.plot(x_values, f3, linewidth=3, label='F3')
ax.set_xlim([5e-3, 5e-1])
ax.set_ylim([0, 2.5])
ax.set_xscale('log')

ax.legend()

plt.savefig('9_ckmt_test.pdf')