#
# Makefile to build the benthic_predator example model
#

include ../../../Rules.make

DOCSRC	= Ndepoden.F90 

OBJS    = ${LIBFABM}(Ndepoden.o)

all: objs

objs: ${OBJS}
	$(MOVE_MODULES_COMMAND)

doc:    $(DOCSRC)
	$(PROTEX) $(DOCSRC) > ../../../../doc/examples_denit_deposit.tex

clean:
	$(RM) *.o *~

#-----------------------------------------------------------------------
# Copyright (C) 2008 - Hans Burchard and Karsten Bolding (BBH)         !
#-----------------------------------------------------------------------
