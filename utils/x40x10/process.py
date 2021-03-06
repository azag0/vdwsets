#!/usr/bin/env python3
import sys
from pathlib import Path
sys.path.append('..')
from caflib.Tools import geomlib
import json


energies = json.load(sys.stdin)
prefix = Path(sys.argv[1])
paths = map(Path, sys.argv[2:])

geoms = []
for path in paths:
    geom = geomlib.readfile(path)
    label = path.stem.split('_')[1]
    idx, label, scale = int(label[0:2]), label[2:-3], float(label[-3:])/100
    if scale == 1.0:
        frags = geom.get_fragments()
        if not (len(frags) == 2 and geomlib.concat(frags) == geom):
            print(
                'error: {} ({}) was not fragmented correctly'.format(label, scale),
                file=sys.stderr
            )
    else:
        frags = []
    geoms.append({'label': label,
                  'idx': idx,
                  'scale': scale,
                  'complex': geom,
                  'fragments': frags})

json.dump(energies, sys.stdout)

for row in geoms:
    idx = '{:02}-{:.2f}'.format(row['idx'], row['scale'])
    row['complex'].write(prefix/'{}-complex-0.xyz'.format(idx))
    for i, fragment in enumerate(row['fragments']):
        fragment.write(prefix/'{}-monomer-{}.xyz'.format(idx, i+1))
