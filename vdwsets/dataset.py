# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.


class Dataset:
    def __init__(self, name):
        self.name = name
        self.clusters = {}
        self._geoms = None

    def __repr__(self):
        if self._geoms is not None:
            ngeoms = len(self._geoms)
        else:
            ngeoms = len(set(
                path
                for cluster in self.clusters.values()
                for path in cluster.fragments.values()
            ))
        return f'<Dataset {self.name!r} containing ' \
            f'{len(self.clusters)} clusters and {ngeoms} structures>'

    def __setitem__(self, key, value):
        self.clusters[key] = value

    def get_int_enes(self, energies, scale=1):
        return {
            key: scale*cluster.get_int_ene(energies[key])
            for key, cluster in self.clusters.items()
            if key in energies
        }

    def load_geoms(self):
        from caflib.Tools import geomlib

        assert self._geoms is None
        path_geoms = {}
        self._geoms = {}
        for cluster in self.clusters.values():
            for name, path in cluster.fragments.items():
                geom = path_geoms.get(path)
                if geom is None:
                    geom = path_geoms.setdefault(path, geomlib.readfile(path))
                    geom = self._geoms.setdefault(geom.hash(), geom)
                cluster.fragments[name] = geom

    def generate_tasks(self, ctx, taskgen):
        self.load_geoms()
        for key, cluster in self.clusters.items():
            for fragment, geom in cluster.fragments.items():
                with ctx.cd(f'{"_".join(map(str, key))}/{fragment}'):
                    taskgen(ctx, geom)


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
