from caflib.Tools.dataset import Dataset, Cluster
from caflib.Tools import geomlib
from caflib.Tools.geomlib import Atom, Molecule, bohr
from csv import DictReader
from pathlib import Path
import json
import re
from math import exp, inf
import numpy as np

root = Path(__file__).parent
kcal = 627.509
kjmol = 2625.5


def intene_cp(x):
    ene_int = x['complex']-x['fragment-1']-x['fragment-2']
    try:
        ene_int_cp = x['complex']-x['fragment-1-cp']-x['fragment-2-cp']
    except KeyError:
        return ene_int
    return ene_int, ene_int_cp


def get_cp_frags(geom):
    frags = geom.get_fragments()
    frags_empty = [frag.copy() for frag in frags]
    for frag in frags_empty:
        for atom in frag:
            atom.flags['dummy'] = True
            atom.prop = atom.prop.copy()
            atom.prop['mass'] *= 0.1
    frags_cp = (frags[0] + frags_empty[1], frags[1] + frags_empty[0])
    for frag in frags_cp:
        frag.metadata['cp'] = True
    return frags_cp


def get_s22(limit=inf):
    ds = Dataset('S22')
    with (root/'s22/energies.json').open() as f:
        enes = json.load(f)
    for idx, row in enumerate(enes):
        if idx >= limit:
            continue
        cluster = Cluster(
            energies={'ref': float(row['CCSD(T) /CBS CP'])},
            intene=intene_cp
        )
        geoms = {}
        geom_names = [
            ('complex', 0, 'complex'),
            ('monomer', 1, 'fragment-1'),
            ('monomer', 2, 'fragment-2'),
        ]
        for name, i, fragment in geom_names:
            filename = '{:02}-{}-{}.xyz'.format(idx+1, name, i)
            geom = geomlib.readfile(root/'s22/geoms'/filename)
            if geom.metadata.get('comment') == '':
                geom.metadata['comment'] = None
            geom['comment'] = f'Formula: {geom!r}'
            geoms[fragment] = geom
        geoms['fragment-1-cp'], geoms['fragment-2-cp'] = get_cp_frags(geoms['complex'])
        for fragment, geom in geoms.items():
            geomid = geom.hash()
            ds.geoms[geomid] = geom
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
            intene=intene_cp
        )
        m = re.match(r'(?P<idx>\d+) (?P<label>.*) \((?P<dist>[\d.]+)\)', row['system name'])
        idx, label, dist = int(m.group('idx')), m.group('label'), float(m.group('dist'))
        if idx > limit:
            continue
        geoms = {}
        geom_names = [
            (dist, 'complex', 0, 'complex'),
            (1.0, 'monomer', 1, 'fragment-1'),
            (1.0, 'monomer', 2, 'fragment-2'),
        ]
        for dist_real, name, i, fragment in geom_names:
            filename = '{:02}-{:.2f}-{}-{}.xyz'.format(idx, dist_real, name, i)
            geom = geomlib.readfile(root/'s66x8/geoms'/filename)
            if geom.metadata.get('comment') == '':
                geom.metadata['comment'] = None
            geoms[fragment] = geom
        geoms['fragment-1-cp'], geoms['fragment-2-cp'] = get_cp_frags(geoms['complex'])
        for fragment, geom in geoms.items():
            geomid = geom.hash()
            ds.geoms[geomid] = geom
            cluster[fragment] = geomid
        ds[(label, dist)] = cluster
    return ds


def get_3b_69(limit=inf):
    ds = Dataset('3B-69')
    with (root/'3b-69/energies.csv').open() as f:
        for row in DictReader(f, delimiter=';'):
            number, system, variant = int(row['number']), row['system'], row['variant']
            if number > inf:
                continue
            if (number, variant) in [(12, 'c'), (19, 'a')]:
                continue
            cluster = Cluster(
                energies={'ref': float(row['CCSD(T)/CBS'])},
                intene=lambda x: x['ABC']-x['AB']-x['BC']-x['AC']+x['A']+x['B']+x['C']
            )
            systempath = system.replace(' ', '_').replace('/', '_')
            path = root/'3b-69/geoms'/f'{number:02}{variant}_{systempath}.xyz'
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
            ds[(f'{system}@{variant}',)] = cluster
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
            if geom.metadata.get('comment') == '':
                geom.metadata['comment'] = None
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
    def intene_cp(x):
        ene_int = x['complex']-x['host']-x['guest']
        try:
            ene_int_cp = (
                x['complex']-x['host']-x['guest'] +
                (-x['fragment-1-cp']-x['fragment-2-cp']) +
                x['fragment-1-nocp']+x['fragment-2-nocp']
            )
        except KeyError:
            return ene_int
        return ene_int, ene_int_cp

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
                intene=intene_cp
            )
            geoms = {}
            for fragment in ['host', 'complex', 'monomer']:
                filename = '{}-{}-{}.xyz'.format(
                    idx, fragment, subidx if fragment != 'host' else ''
                )
                geom = geomlib.readfile(root/'s12l/geoms'/filename)
                geoms[fragment] = geom
            geoms['fragment-1-nocp'], geoms['fragment-2-nocp'] = geoms['complex'].get_fragments()
            geoms['fragment-1-nocp'].metadata['cp'] = True
            geoms['fragment-2-nocp'].metadata['cp'] = True
            geoms['fragment-1-cp'], geoms['fragment-2-cp'] = get_cp_frags(geoms['complex'])
            for fragment, geom in geoms.items():
                geomid = geom.hash()
                ds.geoms[geomid] = geom
                if fragment == 'monomer':
                    fragment = 'guest'
                cluster[fragment] = geomid
            ds[(str(idx) + subidx,)] = cluster
    return ds


def get_sc_frags(crystal, cutoff=6):
    sc = crystal.complete_molecules().supercell((3, 3, 3))
    base = geomlib.Molecule([a for a in sc.atoms if a.flags['cell'] == (1, 1, 1)])
    frags = base.get_fragments()
    cp_shells = []
    for frag in frags:
        dists = np.sqrt(np.sum((frag.xyz[:, None]-sc.xyz[None, :])**2, 2)).min(0)
        cp_shells.append(geomlib.Molecule(
            [sc.atoms[i] for i in np.flatnonzero((0 < dists) & (dists < cutoff))]
        ))
    return frags, cp_shells


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
        geoms = {
            'molecule': Molecule(geomlib.readfile(path, 'xyzc').atoms),
            'crystal': geomlib.readfile(str(path).replace('_g', ''), 'xyzc'),
        }
        for i, (frag, cp_shell) in enumerate(zip(*get_sc_frags(geoms['crystal']))):
            i += 1
            frag.metadata['cp'] = True
            geoms[f'fragment-{i}-nocp'] = frag
            for atom in cp_shell:
                atom.flags['dummy'] = True
                atom.prop = atom.prop.copy()
                atom.prop['mass'] *= 0.1
            frag_cp = frag + cp_shell
            frag_cp.metadata['cp'] = True
            geoms[f'fragment-{i}-cp'] = frag_cp
        for fragment, geom in geoms.items():
            geomid = geom.hash()
            ds.geoms[geomid] = geom
            cluster[fragment] = geomid
        n = len(geoms['crystal'])//len(geoms['molecule'])

        def intene_cp(x, n=n):
            ene_int = x['crystal']/n-x['molecule']
            try:
                ene_int_cp = (
                    x['crystal']-n*x['molecule'] +
                    sum(-x[f'fragment-{i+1}-cp']+x[f'fragment-{i+1}-nocp'] for i in range(n))
                )/n
            except KeyError:
                return ene_int
            return ene_int, ene_int_cp

        cluster._intene = intene_cp
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
