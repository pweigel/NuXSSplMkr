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
    filename = base_path + '/dsdxdy_nu_CC_iso_corrected_logspaced.out'
        
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
    
    y_values = y_values[:-1]
    x_values = x_values[:-1]
    xs = xs[:, :-1, :-1]
    
    spline_path = base_path + '/dsdxdy_nu_CC_iso_corrected_logspaced.fits'
    spline = photospline.SplineTable(spline_path)
    csms_path = '~/utils/xs_iso/dsdxdy_nu_CC_iso.fits'
    csms_spline = photospline.SplineTable(csms_path)
    
    nE = 50
    nx = -25
    ny = -25
    
    x = np.log10(x_values)[nx]
    y = np.log10(y_values)[ny]    
    selected_energy = np.log10(energies[nE]) - 9
    print('ENERGY = ', selected_energy)
    spline_energies = np.linspace(np.log10(energies[0]), np.log10(energies[-1]), 250) - 9
    spline_x = np.linspace(np.log10(x_values[0]), np.log10(x_values[-1]), 250)
    spline_y = np.linspace(np.log10(y_values[0]), np.log10(y_values[-1]), 250)

    fig = plt.figure(figsize=(12, 9))
    ax = fig.add_subplot(111)
    out = spline.evaluate_simple([spline_energies, x, y])
    out_csms = csms_spline.evaluate_simple([spline_energies, x, y])
    ax.scatter(energies / 1e9, xs[:, nx, ny], s=5)
    ax.plot(10**spline_energies, 10**out)
    ax.plot(10**spline_energies, 10**out_csms)
    Emin = 1.0 / (2*0.938*x_values[nx]*y_values[ny])
    ax.plot([Emin, Emin], [1e-40, 1e-32], color='k', alpha=0.5)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylim([1e-40, 1e-32])
    plt.savefig('5_3d_spline_rc_E.pdf')
    
    fig = plt.figure(figsize=(12, 9))
    ax = fig.add_subplot(111)
    out = spline.evaluate_simple([selected_energy, spline_x, y])
    out_csms = csms_spline.evaluate_simple([selected_energy, spline_x, y])
    # print(out)
    ax.scatter(x_values, xs[nE, :, ny], s=5)
    ax.plot(10**spline_x, 10**out)
    ax.plot(10**spline_x, 10**out_csms)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylim([1e-40, 1e-32])
    plt.savefig('5_3d_spline_rc_x.pdf')
    
    fig = plt.figure(figsize=(12, 9))
    ax = fig.add_subplot(111)
    out = spline.evaluate_simple([selected_energy, x, spline_y])
    out_csms = csms_spline.evaluate_simple([selected_energy, x, spline_y])
    # print(out)
    ax.scatter(y_values, xs[nE, nx, :], s=5)
    ax.plot(10**spline_y, 10**out)
    ax.plot(10**spline_y, 10**out_csms)
    # Emin = 1.0 / (2*0.938*x_values[nx]*y_values[ny])
    # ax.plot([Emin, Emin], [1e-40, 1e-32], color='k', alpha=0.5)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylim([1e-40, 1e-32])
    plt.savefig('5_3d_spline_rc_y.pdf')
    # exit()
    
    
    fig = plt.figure(figsize=(12, 9))
    ax = fig.add_subplot(111)
    
    print(selected_energy)
    
    xscorr = np.zeros((len(y_values), len(x_values)))
    xs_spl = np.zeros((len(y_values), len(x_values)))
    for i, y in enumerate(y_values):
        for j, x in enumerate(x_values):
            a = spline.evaluate_simple([selected_energy, np.log10(x), np.log10(y)])
            print(y, x, a, np.log10(xs[nE, j, i]))

            v = 10**a
            w = 10**csms_spline.evaluate_simple([selected_energy, np.log10(x), np.log10(y)])
            # print(x, y, v, w)
            
            xscorr[j, i] = v
            xs_spl[j, i] = w

    _x, _y = np.meshgrid(x_values, y_values)
    ax_img = ax.pcolormesh(_x, _y, xscorr.T, norm=mpl.colors.LogNorm(vmin=1e-42, vmax=1e-37))
    ax.set_xscale('log')
    ax.set_yscale('log')
    plt.colorbar(ax_img, ax=ax)
    plt.savefig('test.pdf')
    
    fig = plt.figure(figsize=(12, 9))
    ax = fig.add_subplot(111)
    _x, _y = np.meshgrid(x_values, y_values)
    ax_img = ax.pcolormesh(_x, _y, np.abs(xs[nE, :, :].T), norm=mpl.colors.LogNorm(vmin=1e-42, vmax=1e-37))#, norm=mpl.colors.LogNorm(vmin=1e-42))
    ax.set_xscale('log')
    ax.set_yscale('log')
    plt.colorbar(ax_img, ax=ax)
    plt.savefig('test2.pdf')
    
    fig = plt.figure(figsize=(12, 9))
    ax = fig.add_subplot(111)
    _x, _y = np.meshgrid(x_values, y_values)
    ax_img = ax.pcolormesh(_x, _y, xs_spl.T, norm=mpl.colors.LogNorm(vmin=1e-42, vmax=1e-37))
    ax.set_xscale('log')
    ax.set_yscale('log')
    plt.colorbar(ax_img, ax=ax)
    plt.savefig('test3.pdf')
    
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
    


    
