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
  
def load_SF_grid(fname):
    with open(fname, 'r') as f:
        lines = f.readlines()
        header_1 = [int(x) for x in lines[0].rstrip('\n').split(' ')]
        NQ2 = header_1[0]
        Nx = header_1[1]
        
        header_2 = [float(x) for x in lines[1].rstrip('\n').split(' ')]
        logQ2min = header_2[0]
        logQ2max = header_2[1]
        logxmin = header_2[2]
        logxmax = header_2[3]
        
        sf_data = []
        
        for line in lines[2:]:
            sf = [float(x) for x in line.rstrip('\n').split(',')]
            sf_data.append(sf)
            
        return NQ2, Nx, logQ2min, logQ2max, logxmin, logxmax, np.array(sf_data)

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
    ax_img = ax.pcolormesh(_x, _y, _y * f1.T, norm=mpl.colors.LogNorm(), rasterized=True)
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
    ax_img = ax.pcolormesh(_x, _y, _y * f3.T, norm=mpl.colors.SymLogNorm(linthresh=1), rasterized=True)
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
    plt.savefig(plot_path + '/' + outfile)
    plt.close()
    
def plot_2d_SF_grid(name, projectile='neutrino', target='proton', sftype='total', 
               Q2min=1.69, Q2max=1e12, xmin=1e-9, xmax=1,
               outfile='2d_sf.pdf'):
    fn1 = data_folder + '/' + name + '/F1_' + projectile + '_' + target + '_' + sftype + '.grid'
    fn2 = data_folder + '/' + name + '/F2_' + projectile + '_' + target + '_' + sftype + '.grid'
    fn3 = data_folder + '/' + name + '/F3_' + projectile + '_' + target + '_' + sftype + '.grid'
    NQ2, Nx, logQ2min, logQ2max, logxmin, logxmax, f1 = load_SF_grid(fn1)
    NQ2, Nx, logQ2min, logQ2max, logxmin, logxmax, f2 = load_SF_grid(fn2)
    NQ2, Nx, logQ2min, logQ2max, logxmin, logxmax, f3 = load_SF_grid(fn3)
    
    print(NQ2, Nx, logQ2min, logQ2max, logxmin, logxmax, f1.shape)
    
    log_Q2_values = np.linspace(logQ2min,logQ2max, NQ2)
    log_x_values = np.linspace(logxmin, logxmax, Nx)
    Q2_values = 10**log_Q2_values
    x_values = 10**log_x_values
    
    _x, _y = np.meshgrid(Q2_values, x_values)

    fig = plt.figure(figsize=(12, 9*3))
    
    ax = fig.add_subplot(311)
    # ax_img = ax.pcolormesh(_x, _y, _y * f1.T, norm=mpl.colors.LogNorm(), rasterized=True)
    ax_img = ax.pcolormesh(_x, _y, _y * f1.T, norm=mpl.colors.SymLogNorm(linthresh=1), rasterized=True)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim([Q2_values[0], Q2_values[-1]])
    ax.set_ylim([x_values[0], x_values[-1]])
    ax.set_xlabel(r'$Q^{2}~[\textrm{GeV}^2]$')
    ax.set_ylabel(r'$x$')
    plt.colorbar(ax_img, ax=ax, label=r'$xF_{1}(x, Q^2)$')
    
    ax = fig.add_subplot(312)
    # ax_img = ax.pcolormesh(_x, _y, f2.T, norm=mpl.colors.LogNorm(), rasterized=True)
    ax_img = ax.pcolormesh(_x, _y, f2.T, norm=mpl.colors.SymLogNorm(linthresh=1), rasterized=True)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim([Q2_values[0], Q2_values[-1]])
    ax.set_ylim([x_values[0], x_values[-1]])
    ax.set_xlabel(r'$Q^{2}~[\textrm{GeV}^2]$')
    ax.set_ylabel(r'$x$')
    plt.colorbar(ax_img, ax=ax, label=r'$F_{2}(x, Q^2)$')
    
    print(np.max(f3), np.min(f3))
    ax = fig.add_subplot(313)
    # ax_img = ax.pcolormesh(_x, _y, _y * f3.T, norm=mpl.colors.LogNorm(), rasterized=True)
    ax_img = ax.pcolormesh(_x, _y, _y * f3.T, norm=mpl.colors.SymLogNorm(linthresh=1), rasterized=True)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim([Q2_values[0], Q2_values[-1]])
    ax.set_ylim([x_values[0], x_values[-1]])
    ax.set_xticks([1e0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e10, 1e11, 1e12])
    ax.set_yticks([1e0, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9])
    ax.set_xlabel(r'$Q^{2}~[\textrm{GeV}^2]$')
    ax.set_ylabel(r'$x$')
    plt.colorbar(ax_img, ax=ax, label=r'$xF_{3}(x, Q^2)$')
    
    plt.tight_layout()
    plt.savefig(plot_path + '/' + outfile)
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

def get_sf_sums(name, projectile, target):
    f1, f2, f3 = None, None, None
    for sftype in ['light', 'charm', 'bottom', 'top']:
        fn1 = data_folder + '/' + name + '/F1_' + projectile + '_' + target + '_' + sftype + '.grid'
        fn2 = data_folder + '/' + name + '/F2_' + projectile + '_' + target + '_' + sftype + '.grid'
        fn3 = data_folder + '/' + name + '/F3_' + projectile + '_' + target + '_' + sftype + '.grid'
        _, _, _, _, _, _, _f1 = load_SF_grid(fn1)
        _, _, _, _, _, _, _f2 = load_SF_grid(fn2)
        _, _, _, _, _, _, _f3 = load_SF_grid(fn3)
        if f1 is None:
            f1 = _f1
        else:
            f1 += _f1
        if f2 is None:
            f2 = _f2
        else:
            f2 += _f2
        if f3 is None:
            f3 = _f3
        else:
            f3 += _f3
    return f1, f2, f3
        

def plot_2d_SF_ratio_grid(name, projectile='neutrino', target='proton', sftype='charm',
               Q2min=1.69, Q2max=1e12, xmin=1e-9, xmax=1,
               outfile='2d_sf_ratio.pdf'):
    fn1 = data_folder + '/' + name + '/F1_' + projectile + '_' + target + '_' + sftype + '.grid'
    fn2 = data_folder + '/' + name + '/F2_' + projectile + '_' + target + '_' + sftype + '.grid'
    fn3 = data_folder + '/' + name + '/F3_' + projectile + '_' + target + '_' + sftype + '.grid'
    NQ2, Nx, logQ2min, logQ2max, logxmin, logxmax, f1 = load_SF_grid(fn1)
    NQ2, Nx, logQ2min, logQ2max, logxmin, logxmax, f2 = load_SF_grid(fn2)
    NQ2, Nx, logQ2min, logQ2max, logxmin, logxmax, f3 = load_SF_grid(fn3)
    
    s = 2 * (0.938) * 1e6 # 1 PeV
    test_1pev_x = np.logspace(np.log10(1/s), 0, 100)
    test_1pev_q2 = s*test_1pev_x
    #W^2 = mt^2 = Q^2 ( 1/x - 1 )
    # Q^2 = mt^2 / (1/x - 1)
    test_xval = np.logspace(-9, 0, 100)
    test_top_w2_threshold =  (173*173)/(1/(1e-9+test_xval) - 1) 
    
    s = 2 * (0.938) * 1e9 # 1 EeV
    test_1eev_x = np.logspace(np.log10(1/s), 0, 100)
    test_1eev_q2 = s*test_1eev_x
    
    tot_f1, tot_f2, tot_f3 = get_sf_sums(name, projectile, target)
    
    print(NQ2, Nx, logQ2min, logQ2max, logxmin, logxmax, f1.shape)
    
    log_Q2_values = np.linspace(logQ2min,logQ2max, NQ2)
    log_x_values = np.linspace(logxmin, logxmax, Nx)
    Q2_values = 10**log_Q2_values
    x_values = 10**log_x_values
    
    _x, _y = np.meshgrid(x_values, Q2_values)

    fig = plt.figure(figsize=(12, 9*3))
    
    ax = fig.add_subplot(311)
    # ax_img = ax.pcolormesh(_x, _y, _y * f1.T, norm=mpl.colors.LogNorm(), rasterized=True)

    ax_img = ax.pcolormesh(_x, _y, _x * f1 / tot_f1, rasterized=True)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylim([Q2_values[0], Q2_values[-1]])
    ax.set_xlim([x_values[0], x_values[-1]])
    ax.set_ylabel(r'$Q^{2}~[\textrm{GeV}^2]$')
    ax.set_xlabel(r'$x$')
    plt.colorbar(ax_img, ax=ax, label=r'$xF_{1}(x, Q^2)$')
    ax.plot(test_1pev_x, test_1pev_q2, linestyle='--', color='w', linewidth=3)
    ax.plot(test_1eev_x, test_1eev_q2, linestyle=':', color='w', linewidth=3)
    ax.plot(test_xval, test_top_w2_threshold, linestyle='-', color='w', linewidth=3)
    
    ax = fig.add_subplot(312)
    # ax_img = ax.pcolormesh(_x, _y, f2.T, norm=mpl.colors.LogNorm(), rasterized=True)
    ax_img = ax.pcolormesh(_x, _y, f2 / tot_f2, rasterized=True)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylim([Q2_values[0], Q2_values[-1]])
    ax.set_xlim([x_values[0], x_values[-1]])
    ax.set_ylabel(r'$Q^{2}~[\textrm{GeV}^2]$')
    ax.set_xlabel(r'$x$')
    plt.colorbar(ax_img, ax=ax, label=r'$F_{2}(x, Q^2)$')
    ax.plot(test_1pev_x, test_1pev_q2, linestyle='--', color='w', linewidth=3)
    ax.plot(test_1eev_x, test_1eev_q2, linestyle=':', color='w', linewidth=3)
    ax.plot(test_xval, test_top_w2_threshold, linestyle='-', color='w', linewidth=3)

    print(np.max(f3), np.min(f3))
    ax = fig.add_subplot(313)
    # ax_img = ax.pcolormesh(_x, _y, _y * f3.T, norm=mpl.colors.LogNorm(), rasterized=True)
    ax_img = ax.pcolormesh(_x, _y, _x * f3 / tot_f3, rasterized=True)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylim([Q2_values[0], Q2_values[-1]])
    ax.set_xlim([x_values[0], x_values[-1]])
    ax.set_yticks([1e0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e10, 1e11, 1e12])
    ax.set_xticks([1e0, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9])
    ax.set_ylabel(r'$Q^{2}~[\textrm{GeV}^2]$')
    ax.set_xlabel(r'$x$')
    plt.colorbar(ax_img, ax=ax, label=r'$xF_{3}(x, Q^2)$')
    
    # draw limit
    # Q = s x

    ax.plot(test_1pev_x, test_1pev_q2, linestyle='--', color='w', linewidth=3)
    ax.plot(test_1eev_x, test_1eev_q2, linestyle=':', color='w', linewidth=3)
    ax.plot(test_xval, test_top_w2_threshold, linestyle='-', color='w', linewidth=3)

    
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
# plot_2d_SF('CT18A_NNLO', Q2min=1, outfile='2d_sf.pdf', projectile='neutrino', target='proton', sftype='total')
# plot_2d_SF_grid('CSMS', outfile='grid_structure_functions_CSMS.pdf')
# plot_2d_SF_grid('CT18A_NNLO', outfile='grid_structure_functions.pdf')
# plot_2d_SF_grid('nCTEQ15_H', sftype='top', outfile='grid_structure_functions_nCT.pdf')
# plot_2d_SF_grid('CT18A_NLO', outfile='grid_structure_functions_CT18ANLO.pdf')
# plot_dsdy_2d('CSMS', projectile='neutrino', target='proton', sftype='total', outfile='2d_dsdy.pdf')

# plot_2d_SF_grid('CT18A_NNLO', sftype='light', outfile='sfs_ct18annlo_light.pdf')
# plot_2d_SF_grid('CT18A_NNLO', sftype='charm', outfile='sfs_ct18annlo_charm.pdf')
# plot_2d_SF_grid('CT18A_NNLO', sftype='bottom', outfile='sfs_ct18annlo_bottom.pdf')
# plot_2d_SF_grid('CT18A_NNLO', sftype='top', outfile='sfs_ct18annlo_top.pdf')

plot_2d_SF_ratio_grid('CT18A_NNLO', sftype='top', outfile='sfs_ct18annlo_top_ratio.pdf')