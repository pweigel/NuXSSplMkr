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
    base_path = '/n/holylfs05/LABS/arguelles_delgado_lab/Everyone/pweigel/sandbox/src/NuXSSplMkr/data/CSMS/cross_sections'
    filename = base_path + '/dsdxdy_nu_CC_iso_corrected.out'
        
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
    
    spline_path = base_path + '/dsdxdy_nu_CC_iso_corrected.fits'
    spline = photospline.SplineTable(spline_path)
    csms_path = '~/utils/xs_iso/dsdxdy_nu_CC_iso.fits'
    csms_spline = photospline.SplineTable(csms_path)
    nx = -11
    ny = -11
    
    x = np.log10(x_values)[nx]
    y = np.log10(y_values)[ny]
    # energies = np.linspace(2, 9, 100)
    
    fig = plt.figure(figsize=(12, 9))
    ax = fig.add_subplot(111)
    out = spline.evaluate_simple([np.log10(energies) - 9.0, x, y])
    # print(xs[10, :, :])
    ax.scatter(energies / 1e9, xs[:, nx, ny], s=5)
    print(out[:25])
    ax.plot(energies / 1e9, 10**out)

    ax.set_xscale('log')
    ax.set_yscale('log')
    plt.savefig('5_3d_spline_rc.pdf')

    fig = plt.figure(figsize=(12, 9))
    ax = fig.add_subplot(111)
    
    nE = 25
    selected_energy = np.log10(energies[nE]) - 9
    print(selected_energy)
    
    xscorr = np.zeros((len(y_values), len(x_values)))
    for i, y in enumerate(y_values):
        for j, x in enumerate(x_values):
            v = 10**spline.evaluate_simple([selected_energy, np.log10(x), np.log10(y)])
            w = 10**csms_spline.evaluate_simple([selected_energy, np.log10(x), np.log10(y)])
            print(x, y, v, w)
            xscorr[i, j] = v

    _x, _y = np.meshgrid(x_values, y_values)
    ax_img = ax.pcolormesh(_x, _y, xscorr, norm=mpl.colors.LogNorm())
    ax.set_xscale('log')
    ax.set_yscale('log')
    plt.colorbar(ax_img, ax=ax)
    plt.savefig('test.pdf')
    
    fig = plt.figure(figsize=(12, 9))
    ax = fig.add_subplot(111)
    _x, _y = np.meshgrid(x_values, y_values)
    ax_img = ax.pcolormesh(_x, _y, xs[nE, :, :], norm=mpl.colors.LogNorm())
    ax.set_xscale('log')
    ax.set_yscale('log')
    plt.colorbar(ax_img, ax=ax)
    plt.savefig('test2.pdf')
    
    # E = 3.0
    # y_values = np.linspace(-5, -1, 100)
    
    # fig = plt.figure(figsize=(12, 9))
    # ax = fig.add_subplot(111)

    # # x = np.log10(0.1)
    # # out = spline.evaluate_simple([E, x, y_values])
    # # ax.plot(10**y_values, 10**out)
    
    # # x = np.log10(0.005)
    # # out = spline.evaluate_simple([E, x, y_values])
    # # ax.plot(10**y_values, 10**out)
    
    # x = np.log10(0.1)
    # out = spline.evaluate_simple([E, x, y_values])
    # ax.plot(10**y_values, 10**out)
    
    # x = np.log10(0.15)
    # out = spline.evaluate_simple([E, x, y_values])
    # ax.plot(10**y_values, 10**out)
    
    # ax.set_xscale('linear')
    # ax.set_yscale('linear')
    
    # plt.savefig('test.pdf')
    


    
