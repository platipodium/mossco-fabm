#
# Makefile to build the MAECS module
#

include ../../../Rules.make

DATE    =   `date +%d%m%y`

LIB	= $(LIBDIR)/libfabm$(buildtype).a

DOCSRC	=  maecs.F90

OBJS   = \
${LIB}(maecs_types.o) \
${LIB}(maecs_functions.o) \
${LIB}(maecs_primprod.o) \
${LIB}(maecs_grazing.o) \
${LIB}(maecs_do.o) \
${LIB}(maecs.o) 

all: objs

objs: ${OBJS}
	$(MOVE_MODULES_COMMAND)

doc:    $(DOCSRC)
	$(PROTEX) $(DOCSRC) > ../../../../doc/fabm.tex 
	touch doc

clean:
	$(RM) *.o *~

backup:	clean
	@echo -n "Creating Archiv of MAECS .. \n"; \
	tar zcvf maecs_$(DATE).tgz *.*  Makefile ; \
	echo "Done.";
