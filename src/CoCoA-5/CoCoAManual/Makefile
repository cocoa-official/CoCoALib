# Makefile for CoCoALib/src/CoCoA-5/CoCoAManual/ directory.

# COCOA_ROOT=../../..
# include $(COCOA_ROOT)/configuration/autoconf.mk
include ../version.mk
COCOA5_VERSION_DEFINES=-DCOCOA5_VER_MAJ=$(COCOA5_VER_MAJ)  -DCOCOA5_VER_MIN=$(COCOA5_VER_MIN)  -DCOCOA5_VER_MINMIN=$(COCOA5_VER_MINMIN)
CWD=src/CoCoA-5/CoCoAManual/

COCOA5_DOC_DIR=CoCoAManual
COCOA5_DOC_XML=CoCoAHelp.xml
COCOA5_DOC_AUX_DIR=aux-files
COCOA5_DOC_HTML_DIR=html
COCOA5_DOC_TEX_DIR=tex

# Default target
.PHONY: pdf-and-html-doc
pdf-and-html-doc: $(COCOA5_DOC_HTML_DIR)/TimeStamp  $(COCOA5_DOC_TEX_DIR)/TimeStamp

# Target for HTML files: uses a TimeStamp to tell whether needs to be rebuilt
$(COCOA5_DOC_HTML_DIR)/TimeStamp:  CoCoAHelp.xml  aux-files/GUI_help.xsl
	/bin/mkdir -p $(COCOA5_DOC_HTML_DIR)
	cd $(COCOA5_DOC_HTML_DIR); \
	touch WillBecomeTimeStamp; \
	java -jar  ../aux-files/saxon.jar  ../$(COCOA5_DOC_XML)  ../aux-files/GUI_help.xsl  CoCoAVersion="$(COCOA5_VER_MAJ).$(COCOA5_VER_MIN).$(COCOA5_VER_MINMIN)" ; \
	/bin/mv WillBecomeTimeStamp TimeStamp



$(COCOA5_DOC_TEX_DIR)/TimeStamp:  CoCoAHelp.xml  aux-files/TeX.xsl  aux-files/TeX-extra-files/mybook.cls
	/bin/mkdir -p $(COCOA5_DOC_TEX_DIR)
	cd $(COCOA5_DOC_TEX_DIR); \
	touch WillBecomeTimeStamp; \
	java -jar ../aux-files/saxon.jar  ../$(COCOA5_DOC_XML)  ../aux-files/TeX.xsl  CoCoAVersion="$(COCOA5_VER_MAJ).$(COCOA5_VER_MIN).$(COCOA5_VER_MINMIN)"  > CoCoAManual.tex
	/bin/cp -r $(COCOA5_DOC_AUX_DIR)/TeX-extra-files/*  $(COCOA5_DOC_TEX_DIR)/.
	@cd $(COCOA5_DOC_TEX_DIR); \
	echo "[[CoCoA5:LaTeX-->PDF]]  Pass 1 of 4"; pdflatex -halt-on-error CoCoAManual.tex < /dev/null > /dev/null  2>&1; \
	 if [ $$? -ne 0 ]; then echo "!!!LaTeX error!!! -- see file $(COCOA5_DOC_TEX_DIR)/CoCoAManual.log"; exit 1; fi; \
	echo "[[CoCoA5:LaTeX-->PDF]]  Pass 2 of 4"; pdflatex -halt-on-error CoCoAManual.tex > /dev/null  2>&1; \
	echo "[[CoCoA5:LaTeX-->PDF]]  Pass 3 of 4"; pdflatex -halt-on-error CoCoAManual.tex > /dev/null  2>&1; \
	echo "[[CoCoA5:LaTeX-->PDF]]  Final pass";  pdflatex -halt-on-error CoCoAManual.tex > /dev/null  2>&1; \
	grep -F "LaTeX Warning:" CoCoAManual.log; \
	/bin/mv WillBecomeTimeStamp TimeStamp

.PHONY:  clean-local
clean-local:
	@/bin/rm -f  texput.log

.PHONY: clean
clean:  clean-local
	@echo "Cleaned CoCoALib/$(CWD)"

.PHONY: veryclean
veryclean:  clean-local
	@/bin/rm -rf  $(COCOA5_DOC_HTML_DIR)  $(COCOA5_DOC_TEX_DIR)
	@echo "Verycleaned CoCoALib/$(CWD)"
