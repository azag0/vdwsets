all: build

build: energies.json

energies.json: process.py | structures.zip energies.csv
	$(eval TMPDIR != mktemp -d)
	tar -C $(TMPDIR) -xf structures.zip
	for f in $(TMPDIR)/*.xyz; do sed -i "" '/[^[:space:]]/,$$!d' $$f; done
	mkdir -p geoms
	../csv2json.py <res/energies.csv | ./$^ geoms $(TMPDIR)/*.xyz >$@
	rm -r $(TMPDIR)

structures.zip:
	wget -O $@ http://www.begdb.com/moldown.php?id=4

energies.csv:
	$(error $@ needs to be downloaded manually at http://www.begdb.com/index.php?action=oneDataset&id=4&state=show)

web:
	open http://dx.doi.org/10.1039/b600027d

