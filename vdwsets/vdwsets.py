# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
from pkg_resources import resource_stream, resource_filename

from .dataset import Dataset, Cluster


def get_s22():
    import pandas as pd
    df = pd.read_csv(
        resource_stream(__name__, 'data/s22/energies.csv'),
        index_col='label scale'.split()
    )
    ds = Dataset('S22', df)
    for (label, scale), row in df.iterrows():
        file_lbl = label.replace(' ', '-').lower()
        fragments = {}
        for j, fragment in enumerate(['complex'] + 2*['fragment']):
            filename = 'data/s22/geoms/' \
                f'{int(row["idx"]):02}_{file_lbl}_{fragment}_{j}.xyz'
            path = resource_filename(__name__, filename)
            fragments[fragment if fragment == 'complex' else f'fragment-{j}'] = path
        ds[label, scale] = Cluster(
            fragments,
            lambda x: x['complex']-x['fragment-1']-x['fragment-2']
        )
    return ds


def get_s12l():
    import pandas as pd
    df = pd.read_csv(
        resource_stream(__name__, 'data/s12l/energies.csv'),
        index_col='label scale'.split()
    )
    ds = Dataset('S12L', df)
    for (label, scale), row in df.iterrows():
        idx, subidx = label
        fragments = {}
        for fragment in ['host', 'complex', 'guest']:
            filename = f'data/s12l/geoms/{idx}' \
                f'-{fragment if fragment != "guest" else "monomer"}' \
                f'-{subidx if fragment != "host" else ""}.xyz'
            fragments[fragment] = resource_filename(__name__, filename)
        ds[label, scale] = Cluster(
            fragments,
            lambda x: x['complex']-x['host']-x['guest']
        )
    return ds


def get_s66x8():
    import pandas as pd
    df = pd.read_csv(
        resource_stream(__name__, 'data/s66x8/energies.csv'),
        index_col='label scale'.split()
    )
    ds = Dataset('S66', df)
    for (label, scale), row in df.iterrows():
        file_lbl = label.replace(' ', '-').lower()
        fragments = {}
        for j, fragment in enumerate(['complex'] + 2*['fragment']):
            scale_lbl = scale if fragment == 'complex' else 1.
            filename = 'data/s66x8/geoms/' \
                f'{int(row["idx"]):02}-{scale_lbl:.2f}_{file_lbl}_{fragment}_{j}.xyz'
            path = resource_filename(__name__, filename)
            fragments[fragment if fragment == 'complex' else f'fragment-{j}'] = path
        ds[label, scale] = Cluster(
            fragments,
            lambda x: x['complex']-x['fragment-1']-x['fragment-2']
        )
    return ds


def get_x23():
    import pandas as pd
    df = pd.read_csv(
        resource_stream(__name__, 'data/x23/energies.csv'),
        index_col='label scale'.split()
    )
    ds = Dataset('X23', df)
    for (label, scale), row in df.iterrows():
        prefix = 'data/x23/geoms/' + row['geomname']
        fragments = {
            'molecule': resource_filename(__name__, prefix + '_g.xyz'),
            'crystal': resource_filename(__name__, prefix + '.xyzc'),
        }
        with open(fragments['molecule']) as fm, open(fragments['crystal']) as fc:
            n = int(next(fc))/int(next(fm))
        ds[label, scale] = Cluster(
            fragments,
            lambda x, n=n: x['crystal']/n-x['molecule']
        )
    return ds


# def get_3b_69():
#     ds = Dataset('3B-69')
#     with (root/'3b-69/energies.csv').open() as f:
#         for row in DictReader(f, delimiter=';'):
#             number, system, variant = int(row['number']), row['system'], row['variant']
#             if (number, variant) in [(12, 'c'), (19, 'a')]:
#                 continue
#             cluster = Cluster(
#                 energies={'ref': float(row['CCSD(T)/CBS'])},
#                 intene=lambda x: x['ABC']-x['AB']-x['BC']-x['AC']+x['A']+x['B']+x['C']
#             )
#             systempath = system.replace(' ', '_').replace('/', '_')
#             path = root/'3b-69/geoms'/f'{number:02}{variant}_{systempath}.xyz'
#             print(path)
#             cmplx = geomlib.readfile(path)
#             frags = cmplx.get_fragments()
#             assert len(frags) == 3 and geomlib.concat(frags) == cmplx
#             geoms = {
#                 'ABC': cmplx,
#                 'AB': frags[0]+frags[1],
#                 'BC': frags[1]+frags[2],
#                 'AC': frags[0]+frags[2],
#                 'A': frags[0],
#                 'B': frags[1],
#                 'C': frags[2]
#             }
#             for fragname, geom in geoms.items():
#                 geomid = geom.hash()
#                 ds.geoms[geomid] = geom
#                 cluster[fragname] = geomid
#             ds[(f'{system}[{variant}]',)] = cluster
#     return ds
#
#
# def get_x40x10(limit=inf):
#     ds = Dataset('X40x10')
#     with (root/'x40x10/energies.json').open() as f:
#         enes = json.load(f)
#     for row in enes:
#         cluster = Cluster(
#             energies={'ref': float(row['CCSD(T) /CBS CP'])},
#             intene=lambda x: x['complex']-x['fragment-1']-x['fragment-2']
#         )
#         m = re.match(r'(?P<idx>\d+) (?P<label>.*) (?P<dist>[\d.]+)', row['system name'])
#         idx, label, dist = int(m.group('idx')), m.group('label'), float(m.group('dist'))
#         if idx > limit:
#             continue
#         for i, name in enumerate(['complex'] + 2*['monomer']):
#             filename = '{:02}-{:.2f}-{}-{}.xyz'.format(
#                 idx, dist if name == 'complex' else 1.0, name, i
#             )
#             geom = geomlib.readfile(root/'x40x10/geoms'/filename)
#             geomid = geom.hash()
#             ds.geoms[geomid] = geom
#             if name == 'complex':
#                 fragment = name
#             else:
#                 fragment = 'fragment-{}'.format(i)
#             cluster[fragment] = geomid
#         ds[(label, dist)] = cluster
#     return ds
#
#
# def get_l7(limit=inf):
#     ds = Dataset('L7')
#     with (root/'l7/energies.json').open() as f:
#         enes = json.load(f)
#     for idx, row in enumerate(enes):
#         if idx >= limit:
#             continue
#         components = sorted((root/'l7/geoms').glob('{}-*.xyz'.format(idx+1)))
#         cluster = Cluster(
#             energies={'ref': float(row['QCISD(T) /CBS CP'] or row['CCSD(T) /CBS CP'])},
#             intene=lambda x, n=len(components)-1: x['complex']-sum(x['fragment-' + str(i+1)] for i in range(n))
#         )
#         for path in components:
#             geom = geomlib.readfile(path)
#             idx, fragment, i = path.stem.split('-')
#             geomid = geom.hash()
#             ds.geoms[geomid] = geom
#             if fragment == 'monomer':
#                 fragment = 'fragment-{}'.format(i)
#             cluster[fragment] = geomid
#         ds[(row['system name'],)] = cluster
#     return ds
#
#
# def get_rare_gas(limit=inf):
#     def eval_potential(R, D, Rmin, A, C6, C8, C10, alpha, beta, **kwargs):
#         DRmin = D*Rmin
#         Fx = exp(-(DRmin/R-1)**2) if R < DRmin else 1
#         E_vdw = -Fx*(C6/R**6+C8/R**8+C10/R**10)
#         E_pauli = A*exp(-alpha*R+beta*R**2)
#         E = E_pauli+E_vdw
#         return E
#
#     ds = Dataset('rare-gas')
#     distances = [
#         2.1, 2.3, 2.5, 2.7, 2.9, 3.1, 3.3, 3.5, 3.7,
#         3.9, 4.1, 4.3, 4.5, 5., 5.5, 6.5, 8., 10.
#     ]
#     with (root/'rare-gas/potentials.json').open() as f:
#         potentials = json.load(f)
#     energies = []
#     for specie, potential in potentials.items():
#         atom = Molecule([Atom(specie, (0, 0, 0))])
#         geomid_atom = atom.hash()
#         ds.geoms[geomid_atom] = atom
#         for i, distance in enumerate(distances):
#             if i >= limit:
#                 continue
#             ene = eval_potential(distance/bohr, **potential)*kcal
#             cluster = Cluster(
#                 energies={'ref': ene},
#                 intene=lambda x: x['dimer']-2*x['atom']
#             )
#             cluster['atom'] = geomid_atom
#             energies.append((specie, distance, kcal*ene))
#             dimer = Molecule([Atom(specie, (0, 0, 0)), Atom(specie, (distance, 0, 0))])
#             geomid = dimer.hash()
#             ds.geoms[geomid] = dimer
#             cluster['dimer'] = geomid
#             ds[(specie, distance)] = cluster
#     return ds
#
#
# def get_all_datasets(include=[], exclude=[], limit=inf):
#     all_ds = {}
#     for get_ds in [
#             get_s22, get_s66x8, get_x40x10, get_s12l, get_x23, get_l7,
#             get_rare_gas, get_3b_69
#     ]:
#         ds = get_ds(limit=limit)
#         if ds.name in exclude or (include and ds.name not in include):
#             continue
#         all_ds[ds.name] = ds
#     return all_ds
#
#
# if __name__ == '__main__':
#     geomlib.settings['eq_precision'] = 3
#     print(get_all_datasets())
