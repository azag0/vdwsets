all: build

build: geoms

geoms: process.ipynb | res/suppl-info.json
	rm -rf $@
	mkdir $@
	jupyter nbconvert $^ --log-level=WARN --to=script --stdout | python3

%.json: %.pdf
	pdf2json -i -enc UTF-8 -f 10 $^ $@

res/suppl-info.pdf: | res
	wget -O $@ http://api.onlinelibrary.wiley.com/asset/v1/doi/10.1002%2Fchem.201200497/asset/supinfo%2Fchem_201200497_sm_miscellaneous_information.pdf?l=SkaBT8QEx2qC%2FkuQm2Z%2BpG9d6XDb3oxTSe%2FNfRai8PRE43YROeSuE8bc2Kf0d%2F1KrrsGKQarqx%2BE%0AXQXDXOoYgKFVxQEAF2iOzVjifLDuY8yWXSoEwFmVog%3D%3D

res:
	mkdir $@

clean:
	rm -rf geoms

distclean: clean
	rm -rf res

web:
	open http://dx.doi.org/10.1002/chem.201200497

