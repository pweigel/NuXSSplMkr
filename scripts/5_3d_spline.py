import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib import rcParams
import photospline

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

def load_file(fname):
    energies = []
    dsdxdy = []
    with open(fname, 'r') as f:
        n = 0
        for line in f.readlines():
            if n == 0:
                energies = np.array([float(x) for x in line.rstrip('\n').split(',')[1:]])
            elif n == 1:
                y_values = np.array([float(x) for x in line.rstrip('\n').split(',')[1:]])
            elif n == 2:
                x_values = np.array([float(x) for x in line.rstrip('\n').split(',')[1:]])
            else:
                d = np.array([float(x) for x in line.rstrip('\n').split(',')])
                dsdxdy.append(d)
            n += 1
    return energies / 1e9, y_values, x_values, np.array(dsdxdy)


if __name__ == '__main__':
    base_path = '/n/holylfs05/LABS/arguelles_delgado_lab/Everyone/pweigel/sandbox/src/NuXSSplMkr/data/CT18A_NNLO/replica_0/cross_sections'
    filename = base_path + '/dsdxdy_neutrino_proton_light.out'
    
    fig = plt.figure(figsize=(12, 9))
    ax = fig.add_subplot(111)
    
    xs = np.loadtxt(filename, skiprows=3, delimiter=',')
    print(xs.shape)
    # 0.00546228,0.00107227
    # Let's get the E, y, x vectors
    with open(filename, 'r') as f:
        energies = np.array([float(x) for x in f.readline().rstrip('\n').split(',')[1:]])
        y_values = np.array([float(x) for x in f.readline().rstrip('\n').split(',')[1:]])
        x_values = np.array([float(x) for x in f.readline().rstrip('\n').split(',')[1:]])
    
    xs = xs.reshape(len(energies), len(y_values), len(x_values))
    print(xs.shape, energies.shape, y_values.shape, x_values.shape)
    
    spline_path = '3d_spline_light.fits'
    spline = photospline.SplineTable(spline_path)
    
    x = np.log10(x_values)[-10]
    y = np.log10(y_values)[-10]
    # energies = np.linspace(2, 9, 100)
    
    fig = plt.figure(figsize=(12, 9))
    ax = fig.add_subplot(111)
    out = spline.evaluate_simple([np.log10(energies) - 9.0, x, y])
    # print(xs[10, :, :])
    ax.scatter(energies / 1e9, xs[:, -10, -10], s=5)
    print(out[:25])
    ax.plot(energies / 1e9, 10**out)

    ax.set_xscale('log')
    ax.set_yscale('log')
    plt.savefig('5_3d_spline.pdf')
    
    E = 8.0
    y_values = np.linspace(-0.1, 0, 100)
    
    fig = plt.figure(figsize=(12, 9))
    ax = fig.add_subplot(111)

    x = np.log10(0.001)
    out = spline.evaluate_simple([E, x, y_values])
    ax.plot(10**y_values, 10**out)
    
    x = np.log10(0.005)
    out = spline.evaluate_simple([E, x, y_values])
    ax.plot(10**y_values, 10**out)
    
    x = np.log10(0.01)
    out = spline.evaluate_simple([E, x, y_values])
    ax.plot(10**y_values, 10**out)
    
    x = np.log10(0.015)
    out = spline.evaluate_simple([E, x, y_values])
    ax.plot(10**y_values, 10**out)
    
    ax.set_xscale('log')
    ax.set_yscale('log')
    
    plt.savefig('test.pdf')

    
