#
# Makefile to build the FABM models from HZG
#

SUBDIRS += omexdia_p
SUBDIRS += maecs
SUBDIRS += benthic_pool
SUBDIRS += n2pzdq
SUBDIRS += jelly

all: objs
	$(MOVE_MODULES_COMMAND)

objs doc clean:
	set -e; for dir in $(SUBDIRS); do $(MAKE) -C $$dir $@ ; done

.PHONY: objs doc clean $(SUBDIRS)
