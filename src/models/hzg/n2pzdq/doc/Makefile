
DOXYGEN ?= doxygen
DOXYGEN := $(shell which $(DOXYGEN))
ifeq (x$(DOXYGEN),x)
$(error  The program "doxygen" required to build the documentation was not found. Either specify the correct location of the doxygen executable with DOXYGEN=/your/path/to/bin/doxygen or type "make help" to learn how to access the online documentation.) 
endif

DOXYGEN_CONFIG := n2pzdq-doxyfile
#EXTRA_DIST=README $(DOXYGEN_CONFIG) doxygen.css doxygen.png

.PHONY: default doc #www

default: doc 

doc:
	$(DOXYGEN) $(DOXYGEN_CONFIG)
	$(MAKE) -C latex
	makeindex ./latex/refman.nlo -s nomencl.ist -o ./latex/refman.nls 	
	$(MAKE) -C latex
	@cp latex/refman.pdf n2pzdq_documentation.pdf

clean:
	@rm -rf latex
	@rm -rf html
	@rm -f n2pzdq_documentation.pdf

