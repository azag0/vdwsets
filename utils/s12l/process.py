#!/usr/bin/env python3
import json
from itertools import groupby
import numpy as np
import sys
import re
from geomlib import Molecule, Atom
from pathlib import Path
from pprint import pprint

bohr = 0.52917721092

with open('res/suppl-info.json') as f:
    pages = json.load(f)


def get_lines(pages):
    for page in pages:
        top_last = 0
        for top, tokens in groupby(page['text'], lambda tok: tok['top']):
            if top < top_last+12:
                print(page['number'], list(tokens))
                return
            top_last = top
            token_data = [tok['data'] for tok in tokens]
            if len(token_data) == 1:
                continue
            yield token_data


lines = get_lines(pages)
while next(lines)[-1] != 'COORDINATES':
    pass
next(lines)
geoms = {}
try:
    for line in lines:
        if len(line) in (2, 3):
            geom = Molecule([])
            if len(line) == 3:
                geom['charge'] = int(line[2].split('CHARGE=')[1])
            geoms[(int(line[0]), line[1].lower())] = geom
        elif len(line) == 4:
            geom.atoms.append(Atom(line[3], tuple(float(x)*bohr for x in line[0:3])))
except:
    print('Line: ', line)
    raise

pprint([
    (idx, system, geom.metadata.get('charge'), len(geom))
    for (idx, system), geom in sorted(geoms.items())
])


for (idx, system), geom in sorted(geoms.items()):
    system, subidx = re.findall(r'(host|complex|monomer)(\d)?', system)[0]
    subidx = 'ab'[int(subidx)-1] if subidx else ''
    filename = '{}-{}-{}.xyz'.format(idx, system, subidx)
    geom.write('geoms/' + filename)

