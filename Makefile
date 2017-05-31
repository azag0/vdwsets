DATASETS = s22
BLDDIR = build
DATADIR = data
UTILDIR = utils

all: $(DATASETS)

$(DATASETS): | $(BLDDIR)/geomlib.py
	mkdir -p $(BLDDIR)/$@
	cp $(UTILDIR)/$@/Makefile $(BLDDIR)/$@
	$(MAKE) -C $(BLDDIR)/$@

$(BLDDIR)/geomlib.py:
	mkdir -p $(@D)
	wget -O $@ https://raw.githubusercontent.com/azag0/caf/a8f7af756c39675156b3754834cd6082ab9159a4/caflib/Tools/$(@F)

clean:
	@for d in $(DIRS); do $(MAKE) -C $$d clean; done

distclean:
	@for d in $(DIRS); do $(MAKE) -C $$d distclean; done
	rm -f geomlib.py
