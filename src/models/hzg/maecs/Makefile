#
# Makefile to build the MAECS module
#

include ../../../Rules.make

LIB	= $(LIBDIR)/libfabm$(buildtype).a

DOCSRC	=  maecs.F90

OBJS   = \
${LIB}(maecs_types.o) \
${LIB}(maecs.o) \
${LIB}(maecs_functions.o) \
${LIB}(maecs_primprod.o) \
${LIB}(maecs_grazing.o) \
${LIB}(maecs_do.o)

all: objs

objs: ${OBJS}
	$(MOVE_MODULES_COMMAND)

doc:    $(DOCSRC)
	$(PROTEX) $(DOCSRC) > ../../../../doc/fabm.tex 
	touch doc

clean:
	$(RM) *.o *~

#-----------------------------------------------------------------------
# Copyright (C) 2008 - Hans Burchard and Karsten Bolding (BBH)         !
#-----------------------------------------------------------------------
