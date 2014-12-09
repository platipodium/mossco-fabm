#
# Makefile to build the size-based jelly model
#

include ../../../Rules.make

DOCSRC	= jelly.F90 

OBJS    = ${LIBFABM}(jelly.o)

all: objs

objs: ${OBJS}
	$(MOVE_MODULES_COMMAND)

doc:    $(DOCSRC)
	$(PROTEX) $(DOCSRC) > ../../../../doc/jelly.tex

clean:
	$(RM) *.o *~