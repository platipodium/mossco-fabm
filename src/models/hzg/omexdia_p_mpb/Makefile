#
# Makefile to build the omexdia_p_mpb model
#

include ../../../Rules.make

DOCSRC	= omexdia_p_mpb.F90

OBJS    = ${LIBFABM}(omexdia_p_mpb.o)

all: objs

${OBJS}: $(FABMBASE)

objs: ${OBJS}
	$(MOVE_MODULES_COMMAND)

doc:    $(DOCSRC)
	$(PROTEX) $(DOCSRC) > ../../../../doc/examples_omexdia_p_mpb.tex

clean:
	$(RM) *.o *~

