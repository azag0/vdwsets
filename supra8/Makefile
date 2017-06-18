all: build

build: geoms

geoms: process.ipynb | res/suppl-info.json
	rm -rf $@
	mkdir $@
	jupyter nbconvert $^ --log-level=WARN --to=script --stdout | python3

%.json: %.pdf
	pdf2json -i -enc UTF-8 -f 2 $^ $@

res/suppl-info.pdf: | res
	wget -O $@ http://www.rsc.org/suppdata/cc/c4/c4cc06722c/c4cc06722c1.pdf

res:
	mkdir $@

clean:
	rm -rf geoms

distclean: clean
	rm -rf res

web:
	open http://dx.doi.org/10.1039/C4CC06722C

