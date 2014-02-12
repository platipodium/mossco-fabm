#
# Makefile to build the benthic_predator example model
#

include ../../../Rules.make

DOCSRC	= n2pzdq.F90 

OBJS    = ${LIBFABM}(n2pzdq.o)

all: objs

objs: ${OBJS}
	$(MOVE_MODULES_COMMAND)

doc:    $(DOCSRC)
	$(PROTEX) $(DOCSRC) > ../../../../doc/n2pzdq.tex

clean:
	$(RM) *.o *~

#-----------------------------------------------------------------------
# Copyright (C) 2008 - Hans Burchard and Karsten Bolding (BBH)         !
#-----------------------------------------------------------------------
