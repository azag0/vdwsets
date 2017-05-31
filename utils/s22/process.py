#!/usr/bin/env python3
import sys
from pathlib import Path
import geomlib
import csv
from difflib import SequenceMatcher


def similarity(a, b):
    return SequenceMatcher(None, a, b).ratio()


energies = list(csv.DictReader(sys.stdin, quotechar="'"))
prefix = Path(sys.argv[1])
paths = map(Path, sys.argv[2:])

geoms = []
for path in paths:
    geom = geomlib.readfile(path)
    code, label = path.stem.split('_')
    code = int(code)
    geom['name'] = label
    frags = geom.get_fragments()
    if not (len(frags) == 2 and geomlib.concat(frags) == geom):
        print(
            f'error: {label} ({code}) was not fragmented correctly',
            file=sys.stderr
        )
        sys.exit(1)
    geoms.append({
        'label': label,
        'code': code,
        'complex': geom,
        'fragments': frags
    })

geoms.sort(key=lambda x: x['code'])
geoms.insert(9, geoms.pop(2))
geom_labels = [g['label'] for g in geoms]
energy_labels = [row['system name'] for row in energies]
energies = [energies[l.index(max(l))] for l in [
    [similarity(a, b) for a in energy_labels] for b in geom_labels
]]

for i, row in enumerate(energies):
    for key, val in row.items():
        try:
            row[key] = float(val)
        except ValueError:
            pass
    row['idx'] = i+1
    row.move_to_end('idx', False)

writer = csv.DictWriter(sys.stdout, energies[0].keys(), quoting=csv.QUOTE_NONNUMERIC)
writer.writeheader()
writer.writerows(energies)
for i, (gs, es) in enumerate(zip(geoms, energies)):
    name = es['system name'].lower().replace(' ', '-')
    gs['complex'].write(prefix/f'{i+1:02}_{name}_complex_0.xyz')
    for j, fragment in enumerate(gs['fragments']):
        fragment.write(prefix/f'{i+1:02}_{name}_fragment_{j+1}.xyz')
