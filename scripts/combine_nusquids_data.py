import numpy as np
import h5py

proj_map = {'neutrino': 'nu', 'antineutrino': 'nubar'}

base_path = '/n/holylfs05/LABS/arguelles_delgado_lab/Everyone/pweigel/sandbox/src/NuXSSplMkr/data/CT18ANNLO_nf6_sx_v2/replica_0/cross_sections'
currents = ['CC', 'NC']
flavors = ['light', 'charm', 'bottom', 'top']
targets = ['proton', 'neutron']
projectiles = ['neutrino', 'antineutrino']
data = {}



for target in targets:
    hf = h5py.File(f'wcg24_{target}.h5', 'w')
        
    # todo: currently these are all the same
    electron = hf.create_group('electron')
    muon = hf.create_group('muon')
    tau = hf.create_group('tau')
    
    energies = None
    zs = None
    
    for proj in projectiles:
        _p = proj_map[proj]
        for current in currents:
            data[f'dsdy_{current}_{proj}'] = None
            data[f's_{current}_{proj}'] = None
            suffix = '4'
            if current == 'NC':
                suffix = '1'
            
            for flavor in flavors:
                dsdy_fn = base_path + f'/nusquids_dsdy_{current}_{proj}_{target}_{flavor}.{suffix}.out'
                s_fn = base_path + f'/nusquids_sigma_{current}_{proj}_{target}_{flavor}.{suffix}.out'
                
                if data[f'dsdy_{current}_{proj}'] is None:
                    data[f'dsdy_{current}_{proj}'] = 10**np.loadtxt(dsdy_fn, skiprows=2).T
                else:
                    data[f'dsdy_{current}_{proj}'] += 10**np.loadtxt(dsdy_fn, skiprows=2).T
                if data[f's_{current}_{proj}'] is None:
                    data[f's_{current}_{proj}'] = np.loadtxt(s_fn, delimiter=',')[:, 1]
                else:
                    data[f's_{current}_{proj}'] += np.loadtxt(s_fn, delimiter=',')[:, 1]
            
            data[f'dsdy_{current}_{proj}'] = np.log10(data[f'dsdy_{current}_{proj}'])
            data[f's_{current}_{proj}'] = np.log10(data[f's_{current}_{proj}'])
            
            for grp in [electron, muon, tau]:
                ds = grp.create_dataset(f's_{current}_{_p}', (551,))
                ds[:] = data[f's_{current}_{proj}']
                
                ds = grp.create_dataset(f'dsdy_{current}_{_p}', (551, 551))
                ds[:, :] = data[f'dsdy_{current}_{proj}']

        if zs is None and energies is None:
            with open(dsdy_fn, 'r') as f:
                energies = [np.log10(float(x)) for x in f.readline().split(',')[1:]]
                zs = [float(x) for x in f.readline().split(',')[1:]]

            energies = np.array(energies)
            ds = hf.create_dataset('energies', (551,))
            ds[:] = energies
            
            zs = np.array(zs)
            ds = hf.create_dataset('zs', (551,))
            ds[:] = zs
        