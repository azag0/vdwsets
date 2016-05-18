sources = structures.zip energies.csv article.pdf
doi = 10.1021/ct2002946

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

src/article.pdf: | src
	wget -O $@ http://pubs.acs.org/doi/pdfplus/10.1021/ct2002946

src/structures.zip: | src
	wget -O $@ http://www.begdb.com/moldown.php?id=26

src/energies.csv: | src
	@echo "Go to"
	@echo "http://www.begdb.com/index.php?action=oneDataset&id=26&state=show"
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
