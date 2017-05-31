#!/usr/bin/env python3
import sys
from pathlib import Path
import geomlib
import csv


energies = list(csv.DictReader(sys.stdin, quotechar="'"))
prefix = Path(sys.argv[1])
paths = map(Path, sys.argv[2:])

geoms = []
for path in paths:
    geom = geomlib.readfile(path)
    label = path.stem.split('_')[1]
    _, label, scale = int(label[0:2]), label[2:-3], float(label[-3:])/100
    if scale == 1.0:
        frags = geom.get_fragments()
        if not (len(frags) == 2 and geomlib.concat(frags) == geom):
            print(
                f'error: {label} ({scale}) was not fragmented correctly',
                file=sys.stderr
            )
    else:
        frags = []
    geoms.append({'complex': geom, 'fragments': frags})

for row in energies:
    for key, val in row.items():
        try:
            row[key] = float(val)
        except ValueError:
            pass
    idx, *words, scale = row['system name'].split()
    row['scale'] = float(scale[1:-1])
    row.move_to_end('scale', False)
    row['idx'] = int(idx)
    row.move_to_end('idx', False)
    row['system name'] = ' '.join(words)

writer = csv.DictWriter(sys.stdout, energies[0].keys(), quoting=csv.QUOTE_NONNUMERIC)
writer.writeheader()
writer.writerows(energies)
for gs, es in zip(geoms, energies):
    idx = f'{es["idx"]:02}-{es["scale"]:.2f}'
    name = es['system name'].lower().replace(' ', '-')
    gs['complex'].write(prefix/f'{idx}_{name}_complex_0.xyz')
    for j, fragment in enumerate(gs['fragments']):
        fragment.write(prefix/f'{idx}_{name}_fragment_{j+1}.xyz')
