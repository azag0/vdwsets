#!/usr/bin/env python3
import sys
sys.path.append('..')
from geomlib import Molecule, Atom
from math import exp
import json

kcal = 627.509
bohr = 0.5291769999999987


def eval_potential(R, D, Rmin, A, C6, C8, C10, alpha, beta, **kwargs):
    DRmin = D*Rmin
    Fx = exp(-(DRmin/R-1)**2) if R < DRmin else 1
    E_vdw = -Fx*(C6/R**6+C8/R**8+C10/R**10)
    E_pauli = A*exp(-alpha*R+beta*R**2)
    E = E_pauli+E_vdw
    return E


distances = [
    2.1, 2.3, 2.5, 2.7, 2.9, 3.1, 3.3, 3.5, 3.7,
    3.9, 4.1, 4.3, 4.5, 5., 5.5, 6.5, 8., 10.
]
with open('res/potentials.json') as f:
    potentials = json.load(f)
energies = []
for specie, potential in potentials.items():
    Molecule([Atom(specie, (0, 0, 0))]).write(
        'geoms/{}--atom.xyz'.format(specie.lower())
    )
    for distance in distances:
        ene = eval_potential(distance/bohr, **potential)
        energies.append((specie, distance, kcal*ene))
        geom = Molecule([Atom(specie, (0, 0, 0)), Atom(specie, (distance, 0, 0))])
        geom.write('geoms/{}-{}-dimer.xyz'.format(specie.lower(), distance))
with open('energies.json', 'w') as f:
    json.dump(energies, f)
