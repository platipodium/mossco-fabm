#
# Makefile to build the MAECS module
##${LIB}(maecs_do.o) \

include ../../../Rules.make

DATE    =   `date +%d%m%y`

LIB	= $(LIBDIR)/libfabm$(buildtype).a

DOCSRC	=  maecs.F90 maecs_do.F90 maecs_functions.F90 maecs_grazing.F90 maecs_primprod.F90 maecs_types.F90

OBJS   = \
${LIB}(maecs_types.o) \
${LIB}(maecs_functions.o) \
${LIB}(maecs_primprod.o) \
${LIB}(maecs_grazing.o) \
${LIB}(maecs.o) 

all: objs

# maecs.o: maecs.F90 maecs_do.F90

objs: ${OBJS}
	$(MOVE_MODULES_COMMAND)

doc:    $(DOCSRC)
	$(MAKE) -C doc 
	
clean:  
	$(RM) *.o *~

backup:	clean
	@echo -n "Creating Archiv of MAECS .. \n"; \
	tar zcvf maecs_$(DATE).tgz *.*  Makefile ; \
	echo "Done.";
	
.PHONY: doc
