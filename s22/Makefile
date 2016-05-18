sources = suppl-info.txt structures.zip energies.csv article.pdf
doi = 10.1039/b600027d

all: download
	$(eval tmpdir := $(shell mktemp -d))
	tar -C ${tmpdir} -xf src/structures.zip
	for f in ${tmpdir}/*.xyz; do sed -i "" '/[^[:space:]]/,$$!d' $$f; done
	mkdir -p build/geoms
	python3 ../sanitize.py <src/energies.csv |\
		PYTHONPATH=.. python3 finalize.py build/geoms ${tmpdir}/*.xyz >build/energies.json
	rm -r ${tmpdir}

web:
	open http://dx.doi.org/${doi}

download: $(addprefix src/,${sources})

build/energies.json: src/energies.csv

src/suppl-info.txt: | src
	wget -O $@ http://www.rsc.org/suppdata/cp/b6/b600027d/geometries-s22.xyz

src/article.pdf: | src
	wget -O $@ http://pubs.rsc.org/en/content/articlepdf/2006/cp/b600027d

src/structures.zip: | src
	wget -O $@ http://www.begdb.com/moldown.php?id=4

src/energies.csv: | src
	@echo "Go to"
	@echo "http://www.begdb.com/index.php?action=oneDataset&id=4&state=show"
	@echo "click on the csv icon and download the file."
	@read -p "Enter path to csv file: " path; \
		[[ -e "$$path" ]] && cp "$$path" $@ \
		|| { echo "Error: $$path does not exist, exiting..." && exit 1; }

src:
	mkdir src

clean:
	-rm -r build

distclean: clean
	-rm -r src
