DIRS := s22 s66x8 l7 s12l x40x10 supra8

.PHONY: $(DIRS)

all: $(DIRS)

$(DIRS): | geomlib.py
	@$(MAKE) -C $@

geomlib.py:
	wget https://raw.githubusercontent.com/azag0/caf/5f9ca6f8c3a430da351584d4756a0509c87276cf/caflib/Tools/$@

clean:
	for d in $(DIRS); do $(MAKE) -C $$d clean; done

distclean:
	for d in $(DIRS); do $(MAKE) -C $$d distclean; done
	rm -f geomlib.py
