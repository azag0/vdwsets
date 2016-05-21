all: build

build: energies.json

energies.json: generate.py | res/potentials.json
	rm -rf geoms
	mkdir geoms
	./$^

clean:
	rm -rf geoms energies.json

distclean: clean
