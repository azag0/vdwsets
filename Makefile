DIRS := s22 s66x8 l7 s12l x40x10 supra8

.PHONY: $(DIRS)

all: $(DIRS)

$(DIRS): | geomlib.py
	@$(MAKE) -C $@

geomlib.py:
	wget https://raw.githubusercontent.com/azag0/caf/79527a6e6247d0bd96790fc72e7ebc2891bd20a7/caflib/Tools/$@

clean:
	for d in $(DIRS); do $(MAKE) -C $$d clean; done

distclean:
	for d in $(DIRS); do $(MAKE) -C $$d distclean; done
	rm -f geomlib.py
