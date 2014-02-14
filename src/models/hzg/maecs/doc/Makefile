
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
	$(DOXYGEN) $(DOXYGEN_CONFIG)
	$(MAKE) -C latex
	@cp latex/refman.pdf maecs_documentation.pdf

clean:
	@rm -rf latex
	@rm -f maecs_documentation.pdf

