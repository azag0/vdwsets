dirs := s22/ s66x8/

.PHONY: ${dirs}

all: ${dirs}

${dirs}: | geomlib.py
	@cd $@ && ${MAKE}
	
geomlib.py:
	wget https://raw.githubusercontent.com/azag0/caf/79527a6e6247d0bd96790fc72e7ebc2891bd20a7/caflib/Tools/geomlib.py

clean:
	for d in ${dirs}; do ${MAKE} -C $$d clean; done

distclean:
	-rm */geomlib.py
	for d in ${dirs}; do ${MAKE} -C $$d distclean; done
