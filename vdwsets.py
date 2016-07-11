from caflib.Tools.dataset import Dataset, Cluster
from caflib.Tools import geomlib
from caflib.Tools.geomlib import Atom, Molecule, bohr
from pathlib import Path
import json
import re
from math import exp

root = Path(__file__).parent
kcal = 627.509


def get_s22():
    ds = Dataset('S22')
    with (root/'s22/energies.json').open() as f:
        enes = json.load(f)
    for idx, row in enumerate(enes):
        cluster = Cluster()
        for i, name in enumerate(['complex'] + 2*['monomer']):
            filename = '{:02}-{}-{}.xyz'.format(idx+1, name, i)
            geom = geomlib.readfile(root/'s22/geoms'/filename)
            geomid = geom.hash()
            ds.geoms[geomid] = geom
            if name == 'complex':
                fragment = name
            else:
                fragment = 'fragment-{}'.format(i)
            cluster[fragment] = geomid
        ds[(row['system name'],)] = cluster
    return ds


def get_s66x8():
    ds = Dataset('S66x8')
    with (root/'s66x8/energies.json').open() as f:
        enes = json.load(f)
    for row in enes:
        cluster = Cluster()
        m = re.match(r'(?P<idx>\d+) (?P<label>.*) \((?P<dist>[\d.]+)\)', row['system name'])
        idx, label, dist = int(m.group('idx')), m.group('label'), float(m.group('dist'))
        for i, name in enumerate(['complex'] + 2*['monomer']):
            filename = '{:02}-{:.2f}-{}-{}.xyz'.format(
                idx, dist if name == 'complex' else 1.0, name, i
            )
            geom = geomlib.readfile(root/'s66x8/geoms'/filename)
            geomid = geom.hash()
            ds.geoms[geomid] = geom
            if name == 'complex':
                fragment = name
            else:
                fragment = 'fragment-{}'.format(i)
            cluster[fragment] = geomid
        ds[(label, dist)] = cluster
    return ds


def get_x40x10():
    ds = Dataset('X40x10')
    with (root/'x40x10/energies.json').open() as f:
        enes = json.load(f)
    for row in enes:
        cluster = Cluster()
        m = re.match(r'(?P<idx>\d+) (?P<label>.*) (?P<dist>[\d.]+)', row['system name'])
        idx, label, dist = int(m.group('idx')), m.group('label'), float(m.group('dist'))
        for i, name in enumerate(['complex'] + 2*['monomer']):
            filename = '{:02}-{:.2f}-{}-{}.xyz'.format(
                idx, dist if name == 'complex' else 1.0, name, i
            )
            geom = geomlib.readfile(root/'x40x10/geoms'/filename)
            geomid = geom.hash()
            ds.geoms[geomid] = geom
            if name == 'complex':
                fragment = name
            else:
                fragment = 'fragment-{}'.format(i)
            cluster[fragment] = geomid
        ds[(label, dist)] = cluster
    return ds


def get_s12l():
    ds = Dataset('S12L')
    for idx in range(2, 8):
        for subidx in 'ab':
            cluster = Cluster()
            for fragment in ['host', 'complex', 'monomer']:
                filename = '{}-{}-{}.xyz'.format(
                    idx, fragment, subidx if fragment != 'host' else ''
                )
                geom = geomlib.readfile(root/'s12l/geoms'/filename)
                geomid = geom.hash()
                ds.geoms[geomid] = geom
                if fragment == 'monomer':
                    fragment = 'guest'
                cluster[fragment] = geomid
            ds[(str(idx) + subidx,)] = cluster
    return ds


def get_x23():
    ds = Dataset('X23')
    for path in (root/'x23/geoms').glob('*_g.xyz'):
        name = path.stem.split('_')[0]
        cluster = Cluster()
        for fragment, geom in [
                ('molecule', Molecule(geomlib.readfile(path, 'xyzc').atoms)),
                ('crystal', geomlib.readfile(str(path).replace('_g', ''), 'xyzc'))
        ]:
            geomid = geom.hash()
            ds.geoms[geomid] = geom
            cluster[fragment] = geomid
        ds[(name,)] = cluster
    return ds


def get_l7():
    ds = Dataset('L7')
    with (root/'l7/energies.json').open() as f:
        enes = json.load(f)
    for idx, row in enumerate(enes):
        cluster = Cluster()
        for path in (root/'l7/geoms').glob('{}-*.xyz'.format(idx+1)):
            geom = geomlib.readfile(path)
            idx, fragment, i = path.stem.split('-')
            geomid = geom.hash()
            ds.geoms[geomid] = geom
            if fragment == 'monomer':
                fragment = 'fragment-{}'.format(i)
            cluster[fragment] = geomid
        ds[(row['system name'],)] = cluster
    return ds


def get_rare_gas():
    def eval_potential(R, D, Rmin, A, C6, C8, C10, alpha, beta, **kwargs):
        DRmin = D*Rmin
        Fx = exp(-(DRmin/R-1)**2) if R < DRmin else 1
        E_vdw = -Fx*(C6/R**6+C8/R**8+C10/R**10)
        E_pauli = A*exp(-alpha*R+beta*R**2)
        E = E_pauli+E_vdw
        return E

    ds = Dataset('rare-gas')
    distances = [
        2.1, 2.3, 2.5, 2.7, 2.9, 3.1, 3.3, 3.5, 3.7,
        3.9, 4.1, 4.3, 4.5, 5., 5.5, 6.5, 8., 10.
    ]
    with (root/'rare-gas/potentials.json').open() as f:
        potentials = json.load(f)
    energies = []
    for specie, potential in potentials.items():
        atom = Molecule([Atom(specie, (0, 0, 0))])
        geomid_atom = atom.hash()
        ds.geoms[geomid_atom] = atom
        for distance in distances:
            cluster = Cluster()
            cluster['atom'] = geomid_atom
            ene = eval_potential(distance/bohr, **potential)
            energies.append((specie, distance, kcal*ene))
            dimer = Molecule([Atom(specie, (0, 0, 0)), Atom(specie, (distance, 0, 0))])
            geomid = dimer.hash()
            ds.geoms[geomid] = dimer
            cluster['dimer'] = geomid
            ds[(specie, distance)] = cluster
    return ds


def get_all_datasets():
    all_ds = {}
    for get_ds in [
            get_s22, get_s66x8, get_x40x10, get_s12l, get_x23, get_l7,
            get_rare_gas
    ]:
        ds = get_ds()
        all_ds[ds.name] = ds
    return all_ds


if __name__ == '__main__':
    geomlib.settings['eq_precision'] = 3
    print(get_rare_gas())
