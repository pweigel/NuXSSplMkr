import numpy as np
import lhapdf
import os, sys, glob

class CrossSections():
    def __init__(self, path, pdfname, nreps, NE=100, Ny=100, A=1, Z=0.5, target_name=r'$\textrm{Isoscalar}$'):
        self.path = path
        self.pdfname = pdfname
        self.nreps = nreps
        
        self.NE = NE
        self.Ny = Ny
        
        self.A = A
        self.Z = Z
        self.N = A - Z
        self.target_name = target_name
        self.get_dsdy()
        self.get_total_xs()
        
    def get_dsdy(self):
      
        E, y, nu_p_light = self.load_replicas('neutrino', 'proton', 'light')
        E, y, nu_n_light = self.load_replicas('neutrino', 'neutron', 'light')
        E, y, nu_p_charm = self.load_replicas('neutrino', 'proton', 'charm')
        E, y, nu_n_charm = self.load_replicas('neutrino', 'neutron', 'charm')
        E, y, nubar_p_light = self.load_replicas('antineutrino', 'proton', 'light')
        E, y, nubar_n_light = self.load_replicas('antineutrino', 'neutron', 'light')
        E, y, nubar_p_charm = self.load_replicas('antineutrino', 'proton', 'charm')
        E, y, nubar_n_charm = self.load_replicas('antineutrino', 'neutron', 'charm')
      
        nu_i_light = (self.Z * nu_p_light + self.N * nu_n_light) / self.A
        nu_i_charm = (self.Z * nu_p_charm + self.N * nu_n_charm) / self.A
        nubar_i_light = (self.Z * nubar_p_light + self.N * nubar_n_light) / self.A
        nubar_i_charm = (self.Z * nubar_p_charm + self.N * nubar_n_charm) / self.A

        nu_i_total = nu_i_light + nu_i_charm
        nubar_i_total = nubar_i_light + nubar_i_charm

        nu_i_tot_cen, nu_i_tot_minus, nu_i_tot_plus = self.get_err_bands(nu_i_total)
        nubar_i_tot_cen, nubar_i_tot_minus, nubar_i_tot_plus = self.get_err_bands(nubar_i_total)

        nu_i_c_cen, nu_i_c_minus, nu_i_c_plus = self.get_err_bands(nu_i_charm)
        nubar_i_c_cen, nubar_i_c_minus, nubar_i_c_plus = self.get_err_bands(nubar_i_charm)
        
        self.dsdy = {'E': E, 'y': y}
        self.dsdy['nu'] = {'total': self.get_err_bands(nu_i_total), 'charm': self.get_err_bands(nu_i_charm)}
        self.dsdy['nubar'] = {'total': self.get_err_bands(nubar_i_total), 'charm': self.get_err_bands(nubar_i_charm)}
        
        self.mean_y = {
                        'nu': {'total': self.get_meany_err_bands(self.get_mean_y(nu_i_total, y)), 'charm': self.get_meany_err_bands(self.get_mean_y(nu_i_charm, y))}, 
                        'nubar': {'total': self.get_meany_err_bands(self.get_mean_y(nubar_i_total, y)), 'charm': self.get_meany_err_bands(self.get_mean_y(nubar_i_charm, y))}
                      }
        
        
    def get_total_xs(self):
        self.NE = 200
        E, nu_p_light = self.load_xs_replicas('neutrino', 'proton', 'light')
        E, nu_n_light = self.load_xs_replicas('neutrino', 'neutron', 'light')
        E, nu_p_charm = self.load_xs_replicas('neutrino', 'proton', 'charm')
        E, nu_n_charm = self.load_xs_replicas('neutrino', 'neutron', 'charm')
        # E, nu_p_bottom = self.load_xs_replicas('neutrino', 'proton', 'bottom')
        # E, nu_n_bottom = self.load_xs_replicas('neutrino', 'neutron', 'bottom')
        # E, nu_p_top = self.load_xs_replicas('neutrino', 'proton', 'top')
        # E, nu_n_top = self.load_xs_replicas('neutrino', 'neutron', 'top')
        E, nubar_p_light = self.load_xs_replicas('antineutrino', 'proton', 'light')
        E, nubar_n_light = self.load_xs_replicas('antineutrino', 'neutron', 'light')
        E, nubar_p_charm = self.load_xs_replicas('antineutrino', 'proton', 'charm')
        E, nubar_n_charm = self.load_xs_replicas('antineutrino', 'neutron', 'charm')
        # E, nubar_p_bottom = self.load_xs_replicas('antineutrino', 'proton', 'bottom')
        # E, nubar_n_bottom = self.load_xs_replicas('antineutrino', 'neutron', 'bottom')
        # E, nubar_p_top = self.load_xs_replicas('antineutrino', 'proton', 'top')
        # E, nubar_n_top = self.load_xs_replicas('antineutrino', 'neutron', 'top')
        
        nu_i_light = (self.Z * nu_p_light + self.N * nu_n_light) / self.A
        nu_i_charm = (self.Z * nu_p_charm + self.N * nu_n_charm) / self.A
        # nu_i_bottom = (self.Z * nu_p_bottom + self.N * nu_n_bottom) / self.A
        # nu_i_top = (self.Z * nu_p_top + self.N * nu_n_top) / self.A
        nubar_i_light = (self.Z * nubar_p_light + self.N * nubar_n_light) / self.A
        nubar_i_charm = (self.Z * nubar_p_charm + self.N * nubar_n_charm) / self.A
        # nubar_i_bottom = (self.Z * nubar_p_bottom + self.N * nubar_n_bottom) / self.A
        # nubar_i_top = (self.Z * nubar_p_top + self.N * nubar_n_top) / self.A

        nu_i_total = nu_i_light + nu_i_charm #+ nu_i_bottom + nu_i_top
        nubar_i_total = nubar_i_light + nubar_i_charm #+ nubar_i_bottom + nubar_i_top

        nu_i_tot_cen, nu_i_tot_minus, nu_i_tot_plus = self.get_xs_err_bands(nu_i_total)
        nubar_i_tot_cen, nubar_i_tot_minus, nubar_i_tot_plus = self.get_xs_err_bands(nubar_i_total)
        
        nu_i_l_cen, nu_i_l_minus, nu_i_l_plus = self.get_xs_err_bands(nu_i_light)
        nubar_i_l_cen, nubar_i_l_minus, nubar_i_l_plus = self.get_xs_err_bands(nubar_i_light)

        nu_i_c_cen, nu_i_c_minus, nu_i_c_plus = self.get_xs_err_bands(nu_i_charm)
        nubar_i_c_cen, nubar_i_c_minus, nubar_i_c_plus = self.get_xs_err_bands(nubar_i_charm)
        
        # nu_i_b_cen, nu_i_b_minus, nu_i_b_plus = self.get_xs_err_bands(nu_i_bottom)
        # nubar_i_b_cen, nubar_i_b_minus, nubar_i_b_plus = self.get_xs_err_bands(nubar_i_bottom)
      
        # nu_i_t_cen, nu_i_t_minus, nu_i_t_plus = self.get_xs_err_bands(nu_i_top)
        # nubar_i_t_cen, nubar_i_t_minus, nubar_i_t_plus = self.get_xs_err_bands(nubar_i_top)
        
        self.total_xs = {'E': E}
        self.total_xs['nu'] = {'total': self.get_xs_err_bands(nu_i_total),
                               'light': self.get_xs_err_bands(nu_i_light),
                               'charm': self.get_xs_err_bands(nu_i_charm)}
        self.total_xs['nubar'] = {'total': self.get_xs_err_bands(nubar_i_total), 
                                  'light': self.get_xs_err_bands(nubar_i_light),
                                  'charm': self.get_xs_err_bands(nubar_i_charm)}
        
        # self.total_xs['nu'] = {'total': self.get_xs_err_bands(nu_i_total), 
        #                        'charm': self.get_xs_err_bands(nu_i_charm),
        #                        'bottom': self.get_xs_err_bands(nu_i_bottom),
        #                        'top': self.get_xs_err_bands(nu_i_top)}
        # self.total_xs['nubar'] = {'total': self.get_xs_err_bands(nubar_i_total), 
        #                           'charm': self.get_xs_err_bands(nubar_i_charm),
        #                           'bottom': self.get_xs_err_bands(nubar_i_bottom),
        #                           'top': self.get_xs_err_bands(nubar_i_top)}

    def load_dsdy(self, fname, linear=False):

        energy_vals = np.logspace(1, 9, self.NE)
        if not linear:
            y_vals = np.logspace(-6, 0, self.Ny)
        else:
            y_vals = np.linspace(1e-3, 1, self.Ny)
        xs_values = np.zeros((self.NE, self.Ny))
        
        with open(fname, 'r') as f:
            n = 0
            for line in f.readlines():
                _xs = float(line.rstrip('\n'))
                xs_values[n // self.Ny, n % self.Ny] = 10**_xs
                n += 1
        return energy_vals, y_vals, xs_values

    def load_replicas(self, projectile='neutrino', target='proton', flavor='light'):
        replicas = []
        for nrep in range(self.nreps):
            fn = self.path + '/replica_{}/cross_sections/dsdy_{}_{}_{}.out'.format(nrep, projectile, target, flavor)
            E, y, xs = self.load_dsdy(fn)
            replicas.append(xs)
        
        replicas = np.array(replicas)
        return E, y, replicas

    def load_xs_replicas(self, projectile='neutrino', target='proton', flavor='light'):
        replicas = []
        for nrep in range(self.nreps):
            fn = self.path + '/replica_{}/cross_sections/total_{}_{}_{}.out'.format(nrep, projectile, target, flavor)
            E, xs = self.load_xs(fn)
            replicas.append(xs)
        
        replicas = np.array(replicas)
        return E, replicas

    def load_xs(self, fname):
        E = []
        xs = []
        with open(fname, 'r') as f:
            for l in f.readlines():
                d = [float(x) for x in l.rstrip('\n').split(',') ]
                E.append(d[0] / 1e9)
                xs.append(10**d[1])
        return np.array(E), np.array(xs)

    def get_err_bands(self, replicas, cl=68.268949):
        s = lhapdf.getPDFSet(self.pdfname)
        central, minus, plus = np.zeros((self.NE, self.Ny)), np.zeros((self.NE, self.Ny)), np.zeros((self.NE, self.Ny))
        # print(replicas.shape)

        for i in range(self.NE):
            for j in range(self.Ny):
                vals = replicas[:, i, j]
                unc = s.uncertainty(vals, cl, False)
                
                central[i, j] = unc.central
                minus[i, j] = unc.errminus
                plus[i, j] = unc.errplus
                
        return central, minus, plus

    def get_meany_err_bands(self, replicas, cl=68.268949):
        s = lhapdf.getPDFSet(self.pdfname)
        central, minus, plus = np.zeros(self.NE), np.zeros(self.NE), np.zeros(self.NE)
        for i in range(self.NE):
            vals = replicas[:, i]
            unc = s.uncertainty(vals, cl, False)
            
            central[i] = unc.central
            minus[i] = unc.errminus
            plus[i] = unc.errplus
                
        return central, minus, plus

    def get_xs_err_bands(self, replicas, cl=68.268949):
        s = lhapdf.getPDFSet(self.pdfname)
        central, minus, plus = np.zeros(self.NE), np.zeros(self.NE), np.zeros(self.NE)
        
        for i in range(self.NE):
            vals = replicas[:, i]
            unc = s.uncertainty(vals, cl, False)
            
            central[i] = unc.central
            minus[i] = unc.errminus
            plus[i] = unc.errplus
                
        return central, minus, plus

    def get_mean_y(self, replicas, y):
        means = np.zeros((self.nreps, self.Ny))
        for nrep in range(self.nreps):
            rep = replicas[nrep, :, :]
            tot=np.zeros(self.Ny)
            val=np.zeros(self.Ny)
            for i in range(self.Ny - 1):
                yavg = 0.5 * (y[i+1]+y[i])
                dy = y[i+1] - y[i]
                val = val + 0.5*(rep[:, i] + rep[:, i+1]) * yavg * dy
                tot = tot + 0.5*(rep[:, i] + rep[:, i+1]) * dy
            means[nrep, :] = val / tot
        return means

