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
    if code == 4109:
        frags[0].join(frags.pop(1))
    elif code == 4110:
        frags[1].join(frags.pop(2))
    elif code == 4107:
        frags[1].join(frags.pop(2))
    elif code == 4112:
        frags[0].join(frags.pop(1))
        frags[1].join(frags.pop(2))
    if geomlib.concat(frags) != geom:
        print(
            'error: {} ({}) was not fragmented correctly'.format(label, code),
            file=sys.stderr
        )
    geoms.append({'label': label,
                  'code': code,
                  'complex': geom,
                  'fragments': frags})

geoms.sort(key=lambda x: x['code'])
geomlbls = [g['label'] for g in geoms]
enelbls = [row['system name'] for row in energies]
energies = [energies[l.index(max(l))] for l in [
    [similarity(a, b) for a in enelbls] for b in geomlbls
]]

json.dump(energies, sys.stdout)
for idx, row in enumerate(geoms):
    row['complex'].write(prefix/'{}-complex-0.xyz'.format(idx+1))
    for i, fragment in enumerate(row['fragments']):
        fragment.write(prefix/'{}-monomer-{}.xyz'.format(idx+1, i+1))
