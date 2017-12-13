
DOXYGEN ?= doxygen
DOXYGEN := $(shell which $(DOXYGEN))
ifeq (x$(DOXYGEN),x)
$(error  The program "doxygen" required to build the documentation was not found. Either specify the correct location of the doxygen executable with DOXYGEN=/your/path/to/bin/doxygen or type "make help" to learn how to access the online documentation.) 
endif

DOXYGEN_CONFIG := maecs-doxyfile
#EXTRA_DIST=README $(DOXYGEN_CONFIG) doxygen.css doxygen.png

.PHONY: default doc #www

default: doc 

doc:
	#combine the problematic maecs_do.F90 with maecs.F90
	@./combine_maecs_maecsdo.sh
	
	#regular doxygen stuff
	$(DOXYGEN) $(DOXYGEN_CONFIG)
	$(MAKE) -C latex
	#build the nomenclature
	makeindex ./latex/refman.nlo -s nomencl.ist -o ./latex/refman.nls 	
	$(MAKE) -C latex
	@cp latex/refman.pdf maecs_documentation.pdf

	#restore the _F90 codes spared by the combine_ script above
	@mv ../maecs_F90 ../maecs.F90
	@mv ../maecs_do_F90 ../maecs_do.F90	

restore:
	#restore the maecs.F90 and maecs_do.F90 that might have been renamed as maecs_F90 and maecs_do_F90
	@mv ../maecs_F90 ../maecs.F90
	@mv ../maecs_do_F90 ../maecs_do.F90
	
clean:
	@rm -rf latex
	@rm -rf html
	@rm maecs_maecsdo_combined.F90
	
realclean: clean
	$(RM) *.dvi *.ps *.pdf
