all: build

install: build
	mkdir -p $(DATADIR)
	cp -r geoms $(DATADIR)

build: geoms

geoms: process.py | res/suppl-info.json
	rm -rf $@ && mkdir $@
	$^

%.json: %.pdf
	pdf2json -i -enc UTF-8 -f 10 $^ $@

res/suppl-info.pdf:
	mkdir -p $(@D)
	wget -O $@ http://api.onlinelibrary.wiley.com/asset/v1/doi/10.1002%2Fchem.201200497/asset/supinfo%2Fchem_201200497_sm_miscellaneous_information.pdf?l=SkaBT8QEx2qC%2FkuQm2Z%2BpG9d6XDb3oxTSe%2FNfRai8PRE43YROeSuE8bc2Kf0d%2F1KrrsGKQarqx%2BE%0AXQXDXOoYgKFVxQEAF2iOzVjifLDuY8yWXSoEwFmVog%3D%3D

clean:
	rm -rf geoms

distclean: clean
	rm -rf res
