import numpy as np
from matplotlib import pyplot as plt
import matplotlib as mpl
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

### STANDARD PARAMETERS ###
NQ2 = 500
Nx = 400
Ny = 100
NE = 200


  
data_folder = '../data'
plot_path = 'tests'

### Loaders ###
def load_dsdy(fname):
    E_values = np.logspace(1, 11, NE)
    y_values = np.logspace(-6, 0, Ny)
    
    xs_values = np.zeros((NE, Ny))
    with open(fname, 'r') as f:
        n = 0
        for line in f.readlines():
            _xs = float(line.rstrip('\n'))
            xs_values[n // Ny, n % Ny] = 10**_xs
            n += 1
    
    return E_values, y_values, xs_values

## SF Plots ##
def plot_2d_SF(name, projectile='neutrino', target='proton', sftype='total', 
               Q2min=1.69, Q2max=1e12, xmin=1e-9, xmax=1,
               outfile='2d_sf.pdf'):
    log_Q2_values = np.linspace(np.log10(Q2min), np.log10(Q2max), NQ2)
    log_x_values = np.linspace(np.log10(xmin), np.log10(xmax), Nx)
    Q2_values = 10**log_Q2_values
    x_values = 10**log_x_values
    
    _x, _y = np.meshgrid(Q2_values, x_values)

    F1 = photospline.SplineTable(data_folder + '/' + name + '/F1_' + projectile + '_' + target + '_' + sftype + '.fits')
    F2 = photospline.SplineTable(data_folder + '/' + name + '/F2_' + projectile + '_' + target + '_' + sftype + '.fits')
    F3 = photospline.SplineTable(data_folder + '/' + name + '/F3_' + projectile + '_' + target + '_' + sftype + '.fits')

    f1 = F1.grideval([log_Q2_values, log_x_values]).squeeze()
    f2 = F2.grideval([log_Q2_values, log_x_values]).squeeze()
    f3 = F3.grideval([log_Q2_values, log_x_values]).squeeze()

    fig = plt.figure(figsize=(12, 9*3))
    
    ax = fig.add_subplot(311)
    ax_img = ax.pcolormesh(_x, _y, _x * f1.T, norm=mpl.colors.LogNorm(), rasterized=True)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim([Q2_values[0], Q2_values[-1]])
    ax.set_ylim([x_values[0], x_values[-1]])
    ax.set_xlabel(r'$Q^{2}~[\textrm{GeV}^2]$')
    ax.set_ylabel(r'$x$')
    plt.colorbar(ax_img, ax=ax, label=r'$xF_{1}(x, Q^2)$')
    
    ax = fig.add_subplot(312)
    ax_img = ax.pcolormesh(_x, _y, f2.T, norm=mpl.colors.LogNorm(), rasterized=True)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim([Q2_values[0], Q2_values[-1]])
    ax.set_ylim([x_values[0], x_values[-1]])
    ax.set_xlabel(r'$Q^{2}~[\textrm{GeV}^2]$')
    ax.set_ylabel(r'$x$')
    plt.colorbar(ax_img, ax=ax, label=r'$F_{2}(x, Q^2)$')
    
    ax = fig.add_subplot(313)
    ax_img = ax.pcolormesh(_x, _y, np.abs(f3.T), norm=mpl.colors.LogNorm(), rasterized=True)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim([Q2_values[0], Q2_values[-1]])
    ax.set_ylim([x_values[0], x_values[-1]])
    ax.set_xticks([1e0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e10, 1e11, 1e12])
    ax.set_yticks([1e0, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9])
    ax.set_xlabel(r'$Q^{2}~[\textrm{GeV}^2]$')
    ax.set_ylabel(r'$x$')
    plt.colorbar(ax_img, ax=ax, label=r'$\vert xF_{3}(x, Q^2) \vert$')
    
    plt.tight_layout()
    plt.savefig(plot_path + '/structure_functions.pdf')
    plt.close()

def plot_1d_SF(name, Q2=[2.0, 4.0, 10.0, 100.0], 
               projectile='neutrino', target='proton', sftype='total', 
               Q2min=1.69, Q2max=1e12, xmin=1e-9, xmax=1,
               outfile='1d_sf.pdf'):
    log_Q2_values = np.log10(Q2)
    log_x_values = np.linspace(np.log10(xmin), np.log10(xmax), Nx)
    Q2_values = 10**log_Q2_values
    x_values = 10**log_x_values
    
    _x, _y = np.meshgrid(Q2_values, x_values)

    F1 = photospline.SplineTable(data_folder + '/' + name + '/F1_' + projectile + '_' + target + '_' + sftype + '.fits')
    F2 = photospline.SplineTable(data_folder + '/' + name + '/F2_' + projectile + '_' + target + '_' + sftype + '.fits')
    F3 = photospline.SplineTable(data_folder + '/' + name + '/F3_' + projectile + '_' + target + '_' + sftype + '.fits')

    f1 = F1.grideval([log_Q2_values, log_x_values]).squeeze()
    f2 = F2.grideval([log_Q2_values, log_x_values]).squeeze()
    f3 = F3.grideval([log_Q2_values, log_x_values]).squeeze()
    fL = f2 - 2*x_values*f1
    
    fig, axes = plt.subplots(ncols=2, nrows=2, figsize=(12*2, 9*2))
    
    for i, ax in enumerate(axes.flatten()):
        ax.plot(x_values, fL[i, :], label=r'$F_{L}$', color='r', linewidth=3)
        ax.plot(x_values, f2[i, :], label=r'$F_{2}$', color='g', linewidth=3)
        ax.plot(x_values, x_values * f3[i, :], label=r'$xF_{3}$', color='b', linewidth=3)
        ax.set_xscale('log')
        ax.set_xlim([1e-4, 1])
        ax.set_ylim([0, 2.5])
        ax.set_xlabel(r'$x$')
        ax.set_ylabel(r'$F_{i}(x, Q^2)$')
        ax.set_title(r'$Q^2 = {:.1f}~\textrm{{GeV}}$'.format(Q2[i]))
        ax.legend()
        
    
    plt.tight_layout()
    
    plt.savefig(plot_path + '/' + outfile)
    plt.close()

def compare_1d_SF(name1, name2, Q2=[2.0, 4.0, 10.0, 100.0], 
                  projectile1='neutrino', projectile2='neutrino', 
                  target1='proton', target2='proton', 
                  sftype1='total', sftype2='total',
                  Q2min=1.69, Q2max=1e12, xmin=1e-9, xmax=1,
                  outfile='compare_1d_sf.pdf'):
    log_Q2_values = np.log10(Q2)
    log_x_values = np.linspace(np.log10(xmin), np.log10(xmax), Nx)
    Q2_values = 10**log_Q2_values
    x_values = 10**log_x_values
    
    _x, _y = np.meshgrid(Q2_values, x_values)

    F1a = photospline.SplineTable(data_folder + '/' + name1 + '/F1_' + projectile1 + '_' + target1 + '_' + sftype1 + '.fits')
    F2a = photospline.SplineTable(data_folder + '/' + name1 + '/F2_' + projectile1 + '_' + target1 + '_' + sftype1 + '.fits')
    F3a = photospline.SplineTable(data_folder + '/' + name1 + '/F3_' + projectile1 + '_' + target1 + '_' + sftype1 + '.fits')

    f1a = F1a.grideval([log_Q2_values, log_x_values]).squeeze()
    f2a = F2a.grideval([log_Q2_values, log_x_values]).squeeze()
    f3a = F3a.grideval([log_Q2_values, log_x_values]).squeeze()
    fLa = f2a - 2*x_values*f1a
    
    F1b = photospline.SplineTable(data_folder + '/' + name2 + '/F1_' + projectile2 + '_' + target2 + '_' + sftype2 + '.fits')
    F2b = photospline.SplineTable(data_folder + '/' + name2 + '/F2_' + projectile2 + '_' + target2 + '_' + sftype2 + '.fits')
    F3b = photospline.SplineTable(data_folder + '/' + name2 + '/F3_' + projectile2 + '_' + target2 + '_' + sftype2 + '.fits')

    f1b = F1b.grideval([log_Q2_values, log_x_values]).squeeze()
    f2b = F2b.grideval([log_Q2_values, log_x_values]).squeeze()
    f3b = F3b.grideval([log_Q2_values, log_x_values]).squeeze()
    fLb = f2b - 2*x_values*f1b
    
    fig, axes = plt.subplots(ncols=2, nrows=2, figsize=(12*2, 9*2))
    
    for i, ax in enumerate(axes.flatten()):
        ax.plot(x_values, fLa[i, :], label=r'$F_{L}^{A}$', color='r', linewidth=3)
        ax.plot(x_values, f2a[i, :], label=r'$F_{2}^{A}$', color='g', linewidth=3)
        ax.plot(x_values, x_values * f3a[i, :], label=r'$xF_{3}^{A}$', color='b', linewidth=3)
        ax.plot(x_values, fLb[i, :], label=r'$F_{L}^{B}$', color='r', linewidth=3, linestyle='--')
        ax.plot(x_values, f2b[i, :], label=r'$F_{2}^{B}$', color='g', linewidth=3, linestyle='--')
        ax.plot(x_values, x_values * f3b[i, :], label=r'$xF_{3}^{B}$', color='b', linewidth=3, linestyle='--')
        ax.set_xscale('log')
        ax.set_xlim([1e-4, 1])
        ax.set_ylim([0, 2.5])
        ax.set_xlabel(r'$x$')
        ax.set_ylabel(r'$F_{i}(x, Q^2)$')
        ax.set_title(r'$Q^2 = {:.1f}~\textrm{{GeV}}$'.format(Q2[i]))
        ax.legend()
        
    plt.tight_layout()
    
    plt.savefig(plot_path + '/' + outfile)
    plt.close()

def plot_dsdy_2d(name, projectile='neutrino', target='proton', sftype='total',
                 outfile='2d_dsdy.pdf'):
    dsdy_fn = data_folder + '/' + name + '/cross_sections/dsdy_' + projectile + '_' + target + '_' + sftype + '.out'
    E, y, dsdy = load_dsdy(dsdy_fn)
    print(E.shape, y.shape, dsdy.shape)
    _x, _y = np.meshgrid(E, y)
    
    fig = plt.figure(figsize=(12, 9))
    ax = fig.add_subplot(111)
    
    ax_img = ax.pcolormesh(_x, _y, dsdy.T, norm=mpl.colors.LogNorm(vmin=1e-42), rasterized=True)
    
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim([E[0], E[-1]])
    ax.set_ylim([y[0], y[-1]])
    ax.set_xlabel(r'$E~[\textrm{GeV}]$')
    ax.set_ylabel(r'$y$')
    
    plt.colorbar(ax_img, ax=ax, label=r'$\frac{d\sigma}{dy}$')
    plt.tight_layout()
    plt.savefig(plot_path + '/' + outfile)
    plt.close()

def plot_sigma(name, 
               projectile='neutrino', target='proton', sftype='total', 
               Q2min=1.69, Q2max=1e12, xmin=1e-9, xmax=1,
               outfile='1d_sf.pdf'):

    log_Q2_values = np.log10(Q2)
    log_x_values = np.linspace(np.log10(xmin), np.log10(xmax), Nx)
    Q2_values = 10**log_Q2_values
    x_values = 10**log_x_values
    
    _x, _y = np.meshgrid(Q2_values, x_values)

    F1 = photospline.SplineTable(data_folder + '/' + name + '/F1_' + projectile + '_' + target + '_' + sftype + '.fits')
    F2 = photospline.SplineTable(data_folder + '/' + name + '/F2_' + projectile + '_' + target + '_' + sftype + '.fits')
    F3 = photospline.SplineTable(data_folder + '/' + name + '/F3_' + projectile + '_' + target + '_' + sftype + '.fits')

    f1 = F1.grideval([log_Q2_values, log_x_values]).squeeze()
    f2 = F2.grideval([log_Q2_values, log_x_values]).squeeze()
    f3 = F3.grideval([log_Q2_values, log_x_values]).squeeze()
    fL = f2 - 2*x_values*f1
    
    fig, axes = plt.subplots(ncols=2, nrows=2, figsize=(12*2, 9*2))
    
    for i, ax in enumerate(axes.flatten()):
        ax.plot(x_values, fL[i, :], label=r'$F_{L}$', color='r', linewidth=3)
        ax.plot(x_values, f2[i, :], label=r'$F_{2}$', color='g', linewidth=3)
        ax.plot(x_values, x_values * f3[i, :], label=r'$xF_{3}$', color='b', linewidth=3)
        ax.set_xscale('log')
        ax.set_xlim([1e-4, 1])
        ax.set_ylim([0, 2.5])
        ax.set_xlabel(r'$x$')
        ax.set_ylabel(r'$F_{i}(x, Q^2)$')
        ax.set_title(r'$Q^2 = {:.1f}~\textrm{{GeV}}$'.format(Q2[i]))
        ax.legend()
        
    
    plt.tight_layout()
    
    plt.savefig(plot_path + '/' + outfile)
    plt.close()
    
# plot_1d_SF('CSMS', sftype='charm'. outfile='1d_sf_charm.pdf')
# compare_1d_SF(name1='CT18A_NNLO', name2='CSMS', outfile='compare_1d_sf.pdf')
plot_2d_SF('CT18A_NNLO', Q2min=1, outfile='2d_sf.pdf', projectile='antineutrino', target='proton', sftype='total')
# plot_dsdy_2d('CSMS', projectile='neutrino', target='proton', sftype='total', outfile='2d_dsdy.pdf')
