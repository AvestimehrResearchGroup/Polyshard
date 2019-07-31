SAGE=sage
SAGES=$(wildcard *.sage)
PYS=$(SAGES:%.sage=%.py)
all: $(PYS) test readme.md

.PHONY: build
build: $(PYS)

test: $(PYS)
	@echo "Testing the files"
	$(SAGE) -c "from codinglib import *" 1>/dev/null

%.py: %.sage
	-rm $*.py $*.pyc
	$(SAGE) -preparse $*.sage
	-mv $*.sage.py $*.py

readme.md: $(SAGES) $(PYS) readme.md.pre
	$(SAGE) generate_readme.py > _doclines.tmp
	./merge_doclines > readme.md

readme.html: readme.md
	markdown readme.md > readme.html

clean:
	-rm -rf $(PYS) $(SAGES:%sage=%pyc)

