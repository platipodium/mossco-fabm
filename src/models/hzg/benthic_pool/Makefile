#
# Makefile to build the benthic_predator example model
#

include ../../../Rules.make

DOCSRC	= benthic_pool.F90 

OBJS    = ${LIBFABM}(benthic_pool.o)

all: objs

objs: ${OBJS}
	$(MOVE_MODULES_COMMAND)

doc:    $(DOCSRC)
	$(PROTEX) $(DOCSRC) > ../../../../doc/examples_benthic_pool.tex

clean:
	$(RM) *.o *~

#-----------------------------------------------------------------------
# Copyright (C) 2008 - Hans Burchard and Karsten Bolding (BBH)         !
#-----------------------------------------------------------------------
