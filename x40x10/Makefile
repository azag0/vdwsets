all: build

build: energies.json

energies.json: process.py | res/structures.zip res/energies.csv
	$(eval TMPDIR != mktemp -d)
	tar -C $(TMPDIR) -xf res/structures.zip
	for f in $(TMPDIR)/*.xyz; do sed -i "" '/[^[:space:]]/,$$!d' $$f; done
	mkdir -p geoms
	../csv2json.py <res/energies.csv | ./$^ geoms $(TMPDIR)/*.xyz >$@
	rm -r $(TMPDIR)

download: $(addprefix src/,${sources})

res/structures.zip:
	wget -O $@ http://www.begdb.com/moldown.php?id=32

res/energies.csv:
	$(error $@ needs to be downloaded manually at http://www.begdb.com/index.php?action=oneDataset&id=32&state=show)

res:
	mkdir $@

clean:
	rm -rf geoms energies.json

distclean: clean
	rm -f res/structures.zip

web:
	open http://dx.doi.org/10.1021/ct2002946
