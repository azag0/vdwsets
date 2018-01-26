ID = 26

all: build

install: build
	mkdir -p $(DATADIR)
	cp -r energies.csv geoms $(DATADIR)

build: energies.csv

energies.csv: process.py clean-csv.py | res/structures.zip res/energies.csv
	$(eval TMPDIR := $(shell mktemp -d))
	tar -C $(TMPDIR) -xf res/structures.zip
	for f in $(TMPDIR)/*.xyz; do sed -i "" '/[^[:space:]]/,$$!d' $$f; done
	mkdir -p geoms
	$(word 2,$^) <res/energies.csv | $< geoms $(TMPDIR)/*.xyz >$(TMPDIR)/$@
	mv $(TMPDIR)/$@ .
	rm -r $(TMPDIR)

res/structures.zip:
	mkdir -p $(@D)
	wget -O $@ http://www.begdb.com/moldown.php?id=$(ID)

res/energies.csv:
	$(error $@ needs to be downloaded manually at http://www.begdb.com/index.php?action=oneDataset&id=$(ID)&state=show)

clean:
	rm -rf geoms energies.csv

distclean: clean
	rm -rf res
