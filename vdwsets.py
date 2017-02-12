from caflib.Tools.dataset import Dataset, Cluster
from caflib.Tools import geomlib
from caflib.Tools.geomlib import Atom, Molecule, bohr
from csv import DictReader
from pathlib import Path
import json
import re
from math import exp, inf

root = Path(__file__).parent
kcal = 627.509
kjmol = 2625.5


def get_s22(limit=inf):
    ds = Dataset('S22')
    with (root/'s22/energies.json').open() as f:
        enes = json.load(f)
    for idx, row in enumerate(enes):
        if idx >= limit:
            continue
        cluster = Cluster(
            energies={'ref': float(row['CCSD(T) /CBS CP'])},
            intene=lambda x: x['complex']-x['fragment-1']-x['fragment-2']
        )
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


def get_s66x8(limit=inf):
    ds = Dataset('S66x8')
    with (root/'s66x8/energies.json').open() as f:
        enes = json.load(f)
    for row in enes:
        cluster = Cluster(
            energies={'ref': float(row['CCSD(T) /CBS CP'])},
            intene=lambda x: x['complex']-x['fragment-1']-x['fragment-2']
        )
        m = re.match(r'(?P<idx>\d+) (?P<label>.*) \((?P<dist>[\d.]+)\)', row['system name'])
        idx, label, dist = int(m.group('idx')), m.group('label'), float(m.group('dist'))
        if idx > limit:
            continue
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


def get_3b_69():
    ds = Dataset('3B-69')
    with (root/'3b-69/energies.csv').open() as f:
        for row in DictReader(f, delimiter=';'):
            number, system, variant = int(row['number']), row['system'], row['variant']
            if (number, variant) in [(12, 'c'), (19, 'a')]:
                continue
            cluster = Cluster(
                energies={'ref': float(row['CCSD(T)/CBS'])},
                intene=lambda x: x['ABC']-x['AB']-x['BC']-x['AC']+x['A']+x['B']+x['C']
            )
            systempath = system.replace(' ', '_').replace('/', '_')
            path = root/'3b-69/geoms'/f'{number:02}{variant}_{systempath}.xyz'
            print(path)
            cmplx = geomlib.readfile(path)
            frags = cmplx.get_fragments()
            assert len(frags) == 3 and geomlib.concat(frags) == cmplx
            geoms = {
                'ABC': cmplx,
                'AB': frags[0]+frags[1],
                'BC': frags[1]+frags[2],
                'AC': frags[0]+frags[2],
                'A': frags[0],
                'B': frags[1],
                'C': frags[2]
            }
            for fragname, geom in geoms.items():
                geomid = geom.hash()
                ds.geoms[geomid] = geom
                cluster[fragname] = geomid
            ds[(f'{system}[{variant}]',)] = cluster
    return ds


def get_x40x10(limit=inf):
    ds = Dataset('X40x10')
    with (root/'x40x10/energies.json').open() as f:
        enes = json.load(f)
    for row in enes:
        cluster = Cluster(
            energies={'ref': float(row['CCSD(T) /CBS CP'])},
            intene=lambda x: x['complex']-x['fragment-1']-x['fragment-2']
        )
        m = re.match(r'(?P<idx>\d+) (?P<label>.*) (?P<dist>[\d.]+)', row['system name'])
        idx, label, dist = int(m.group('idx')), m.group('label'), float(m.group('dist'))
        if idx > limit:
            continue
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


def get_s12l(limit=inf):
    ds = Dataset('S12L')
    with (root/'s12l/energies.csv').open() as f:
        lines = [l.strip().split(';') for l in f]
    refs = {line[0]: {
        refname: float(ene) if ene else None for refname, ene
        in zip(lines[0][1:], line[1:])
    } for line in lines[1:]}
    for idx in range(2, 8):
        if idx-1 > limit:
            continue
        for subidx in 'ab':
            cluster = Cluster(
                energies=refs[str(idx) + subidx],
                intene=lambda x: x['complex']-x['host']-x['guest']
            )
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


def get_x23(limit=inf):
    ds = Dataset('X23')
    with (root/'x23/energies.csv').open() as f:
        lines = [l.strip().split(';') for l in f]
    refs = {line[0]: float(line[2]) for line in lines[1:]}
    for i, path in enumerate(sorted((root/'x23/geoms').glob('*_g.xyz'))):
        if i >= limit:
            continue
        name = path.stem.split('_')[0]
        cluster = Cluster(energies={'ref': refs[name]/kjmol*kcal})
        for fragment, geom in [
                ('molecule', Molecule(geomlib.readfile(path, 'xyzc').atoms)),
                ('crystal', geomlib.readfile(str(path).replace('_g', ''), 'xyzc'))
        ]:
            geomid = geom.hash()
            ds.geoms[geomid] = geom
            cluster[fragment] = geomid
        n = len(ds.geoms[cluster.fragments['crystal']])//len(ds.geoms[cluster.fragments['molecule']])
        cluster._intene = lambda x, n=n: x['crystal']/n-x['molecule']
        ds[(name,)] = cluster
    return ds


def get_l7(limit=inf):
    ds = Dataset('L7')
    with (root/'l7/energies.json').open() as f:
        enes = json.load(f)
    for idx, row in enumerate(enes):
        if idx >= limit:
            continue
        components = sorted((root/'l7/geoms').glob('{}-*.xyz'.format(idx+1)))
        cluster = Cluster(
            energies={'ref': float(row['QCISD(T) /CBS CP'] or row['CCSD(T) /CBS CP'])},
            intene=lambda x, n=len(components)-1: x['complex']-sum(x['fragment-' + str(i+1)] for i in range(n))
        )
        for path in components:
            geom = geomlib.readfile(path)
            idx, fragment, i = path.stem.split('-')
            geomid = geom.hash()
            ds.geoms[geomid] = geom
            if fragment == 'monomer':
                fragment = 'fragment-{}'.format(i)
            cluster[fragment] = geomid
        ds[(row['system name'],)] = cluster
    return ds


def get_rare_gas(limit=inf):
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
        for i, distance in enumerate(distances):
            if i >= limit:
                continue
            ene = eval_potential(distance/bohr, **potential)*kcal
            cluster = Cluster(
                energies={'ref': ene},
                intene=lambda x: x['dimer']-2*x['atom']
            )
            cluster['atom'] = geomid_atom
            energies.append((specie, distance, kcal*ene))
            dimer = Molecule([Atom(specie, (0, 0, 0)), Atom(specie, (distance, 0, 0))])
            geomid = dimer.hash()
            ds.geoms[geomid] = dimer
            cluster['dimer'] = geomid
            ds[(specie, distance)] = cluster
    return ds


def get_all_datasets(include=[], exclude=[], limit=inf):
    all_ds = {}
    for get_ds in [
            get_s22, get_s66x8, get_x40x10, get_s12l, get_x23, get_l7,
            get_rare_gas, get_3b_69
    ]:
        ds = get_ds(limit=limit)
        if ds.name in exclude or (include and ds.name not in include):
            continue
        all_ds[ds.name] = ds
    return all_ds


if __name__ == '__main__':
    geomlib.settings['eq_precision'] = 3
    print(get_all_datasets())
