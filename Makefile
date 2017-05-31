DATASETS = s22 s66x8
BLDDIR = build
DATADIR = data
UTILDIR = utils

all: $(DATASETS)

$(DATASETS): | $(BLDDIR)/geomlib.py
	mkdir -p $(BLDDIR)/$@
	$(MAKE) -f $(abspath $(UTILDIR)/$@/Makefile) -C $(BLDDIR)/$@ install VPATH=$(abspath $(UTILDIR)):$(abspath $(UTILDIR)/$@) PYTHONPATH=$(PYTHONPATH):$(abspath $(BLDDIR)) DATADIR=$(abspath $(DATADIR))/$@

$(BLDDIR)/geomlib.py:
	mkdir -p $(@D)
	wget -O $@ https://raw.githubusercontent.com/azag0/caf/a8f7af756c39675156b3754834cd6082ab9159a4/caflib/Tools/$(@F)

clean:
	@for d in $(DATASETS); do $(MAKE) -f $(abspath $(UTILDIR)/$$d/Makefile) -C $(BLDDIR)/$$d clean; done

distclean:
	@for d in $(DATASETS); do $(MAKE) -f $(abspath $(UTILDIR)/$$d/Makefile) -C $(BLDDIR)/$$d distclean; done
	rm -f geomlib.py
	rmdir build
