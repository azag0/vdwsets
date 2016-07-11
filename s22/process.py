#!/usr/bin/env python3
import sys
from pathlib import Path
sys.path.append('..')
import geomlib
import json
from difflib import SequenceMatcher


def similarity(a, b):
    return SequenceMatcher(None, a, b).ratio()


energies = json.load(sys.stdin)
prefix = Path(sys.argv[1])
paths = map(Path, sys.argv[2:])

geoms = []
for path in paths:
    geom = geomlib.readfile(path)
    code, label = path.stem.split('_')
    code = int(code)
    frags = geom.get_fragments()
    if not (len(frags) == 2 and geomlib.concat(frags) == geom):
        print(
            'error: {} ({}) was not fragmented correctly'.format(label, code),
            file=sys.stderr
        )
    geoms.append({'label': label,
                  'code': code,
                  'complex': geom,
                  'fragments': frags})

geoms.sort(key=lambda x: x['code'])
geom_labels = [g['label'] for g in geoms]
energy_labels = [row['system name'] for row in energies]
energies = [energies[l.index(max(l))] for l in [
    [similarity(a, b) for a in energy_labels] for b in geom_labels
]]

json.dump(energies, sys.stdout)
for idx, row in enumerate(geoms):
    row['complex'].write(prefix/'{:02}-complex-0.xyz'.format(idx+1))
    for i, fragment in enumerate(row['fragments']):
        fragment.write(prefix/'{:02}-monomer-{}.xyz'.format(idx+1, i+1))
