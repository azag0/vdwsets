# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
from itertools import chain


class Dataset:
    def __init__(self, name):
        self.name = name
        self.clusters = {}

    @property
    def geoms(self):
        return set(chain.from_iterable(c.fragments.values() for c in self.clusters.values()))

    def __repr__(self):
        return (
            f'<Dataset {self.name!r} containing '
            f'{len(self.clusters)} clusters and {len(self.geoms)} structures>'
        )

    def get_task(self, ctx, taskgen):
        tasks = {
            geomid: taskgen(ctx, geom, self.name)
            for geomid, geom in self.geoms.items()
        }
        tasktree = [(
            key,
            [
                tasks[geomid] + ctx.link(fragment)
                for fragment, geomid
                in cluster.fragments.items()
                if tasks[geomid]
            ] + ctx()
        ) for key, cluster in self.clusters.items()]
        tasktree.sort(key=lambda x: x[0])
        return [
            task + ctx.link('_'.join(str(k) for k in key))
            for key, task in tasktree
        ] + ctx()

    def __setitem__(self, key, value):
        self.clusters[key] = value

    def get_int_enes(self, energies, scale=1):
        return {
            key: scale*cluster.get_int_ene(energies[key])
            for key, cluster in self.clusters.items()
            if key in energies
        }


class Cluster:
    def __init__(self, fragments=None, energies=None, intene=None):
        self.fragments = fragments or {}
        self.energies = energies or {}
        self._intene = intene

    def __repr__(self):
        return f'Cluster({self.fragments!r})'

    def __setitem__(self, key, value):
        self.fragments[key] = value

    def get_int_ene(self, energies):
        return self._intene(energies)
