#!/usr/bin/env python
import sys
from pathlib import Path
import geomlib
import json
import csv
from difflib import SequenceMatcher


energies = json.load(sys.stdin)
prefix = Path(sys.argv[1])
paths = map(Path, sys.argv[2:])


def similar(a, b):
    return SequenceMatcher(None, a, b).ratio()


geoms = []
for path in paths:
    geom = geomlib.readfile(path)
    code, label = path.stem.split('_')
    code = int(code)
    frags = geom.get_fragments()
    if code == 4107:
        frags[1].join(frags[2])
    elif code == 4109:
        frags[0].join(frags[1])
        frags[1] = frags[2]
    elif code == 4110:
        frags[1].join(frags[2])
    elif code == 4112:
        frags[0].join(frags[1])
        frags[1] = frags[2].joined(frags[3])
    frags = frags[:2]
    if not (len(frags) == 2 and geomlib.concat(frags) == geom):
        print('error: {} ({}) was not fragmented correctly'.format(label, code),
              file=sys.stderr)
    geoms.append({'label': label,
                  'code': code,
                  'complex': geom,
                  'fragments': frags})

geoms.sort(key=lambda x: x['code'])
geomlbls = [g['label'] for g in geoms]
enelbls = [row['system name'] for row in energies]
energies = [energies[l.index(max(l))] for l in
            [[similar(a, b) for a in enelbls] for b in geomlbls]]

writer = csv.DictWriter(sys.stdout, fieldnames=energies[0].keys())
writer.writeheader()
writer.writerows(energies)

for idx, row in enumerate(geoms):
    row['complex'].write(prefix/'{}-complex.xyz'.format(idx+1))
    for i in range(2):
        row['fragments'][i].write(prefix/'{}-monomer-{}.xyz'.format(idx+1, i+1))
