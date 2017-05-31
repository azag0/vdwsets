# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
from pathlib import Path
import csv
from math import inf

from .dataset import Dataset, Cluster

datadir = Path(__file__).parent/'data'
kcal = 627.509
kjmol = 2625.5


def get_s22(limit=inf):
    ds = Dataset('S22')
    with open(datadir/'s22/energies.csv') as f:
        enes = list(csv.DictReader(f, quoting=csv.QUOTE_NONNUMERIC))
    for i, row in enumerate(enes):
        if i >= limit:
            continue
        cluster = Cluster(
            energies={'ref': float(row['CCSD(T)/CBS CP'])},
            intene=lambda x: x['complex']-x['fragment-1']-x['fragment-2']
        )
        name = row['system name'].replace(' ', '-')
        for j, fragment in enumerate(['complex'] + 2*['fragment']):
            path = datadir/'s22/geoms'/f'{i+1:02}_{name}_{fragment}_{j}.xyz'
            assert path.exists()
            if fragment == 'complex':
                key = 'complex'
            else:
                key = f'fragment-{j}'
            cluster[key] = str(path)
        ds[row['system name']] = cluster
    return ds


def get_s66x8(limit=inf):
    ds = Dataset('S66x8')
    with open(datadir/'s66x8/energies.csv') as f:
        enes = list(csv.DictReader(f, quoting=csv.QUOTE_NONNUMERIC))
    for i, row in enumerate(enes):
        if i >= limit:
            continue
        cluster = Cluster(
            energies={'ref': float(row['CCSD(T)/CBS CP'])},
            intene=lambda x: x['complex']-x['fragment-1']-x['fragment-2']
        )
        name = row['system name'].replace(' ', '-')
        for j, fragment in enumerate(['complex'] + 2*['fragment']):
            scale = row["scale"] if fragment == 'complex' else 1.
            path = datadir/'s66x8/geoms'/f'{int(row["idx"]):02}-{scale:.2f}_{name}_{fragment}_{j}.xyz'
            assert path.exists()
            if fragment == 'complex':
                key = 'complex'
            else:
                key = f'fragment-{j}'
            cluster[key] = str(path)
        ds[(row['system name'], row['scale'])] = cluster
    return ds
