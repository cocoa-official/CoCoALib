# Makefile for CoCoALib/doc/ directory.

COCOA_ROOT=..
include $(COCOA_ROOT)/configuration/autoconf.mk
CWD=doc/

SRC_DIR=txt
TEX_DIR=tex
TASK_DIR=CoCoALib-tasks
HTML_DIR=html
AUX_DIR=aux-txt2tags
DOC_SRCS = $(wildcard $(SRC_DIR)/*.txt)
DOC_HTML = $(patsubst $(SRC_DIR)/%.txt,$(HTML_DIR)/%.html,$(DOC_SRCS))
TASKS_OUT = CoCoALib-tasks.html  CoCoALib-CompletedTasks.html

.SUFFIXES:

# Rule for compiling *.txt files into *.html in HTML_DIR
$(HTML_DIR)/%.html:  $(SRC_DIR)/%.txt
	@echo "Generating $*.html"
	@txt2tags -t html -o $@ -i $<  >/dev/null


.PHONY : alldoc
alldoc: examples-index  htmldoc  texdoc  tasktables

$(TEX_DIR):
	@/bin/mkdir -p $(TEX_DIR)

$(TEX_DIR)/CoCoALib.tex:  $(DOC_SRCS)  $(TEX_DIR)  $(AUX_DIR)/DocTeX.t2t
	@echo "Generating CoCoALib.tex..."
	@sed -e "s/VERSION_NUMBER/${COCOALIB_VERSION}/" -e "s/DATE_OF_CREATION/`date +%Y-%m-%d`/"  $(AUX_DIR)/DocTeX.t2t  > $(AUX_DIR)/DocTeX-VERSION.t2t
	@txt2tags -o $(TEX_DIR)/CoCoALib.tex  $(AUX_DIR)/DocTeX-VERSION.t2t   >/dev/null

$(TEX_DIR)/cocoalib-doc.sty:  $(TEX_DIR)  $(AUX_DIR)/cocoalib-doc.sty
	@/bin/cp -p $(AUX_DIR)/cocoalib-doc.sty $(TEX_DIR)


CoCoALib.pdf:  $(COCOA_ROOT)/configuration/version  $(DOC_SRCS)  $(AUX_DIR)/DocTeX.t2t  $(AUX_DIR)/cocoalib-doc.sty
	@command -v txt2tags >/dev/null 2>&1 ; \
	 if [ $$? -ne 0 ]; \
	 then \
	   echo ">>>>>>>>>>>>>>>>>>>>>> WARNING <<<<<<<<<<<<<<<<<<<<<<"; \
	   echo ">>>  Cannot build CoCoALib PDF documentation!     <<<"; \
	   echo ">>>  Install txt2tags from http://txt2tags.org/   <<<"; \
	   echo ">>>  (make sure you get version 3.0 or 2.7-DEV)   <<<"; \
	   echo ">>>>>>>>>>>>>>>>>>>> END WARNING <<<<<<<<<<<<<<<<<<<<"; \
	   echo; \
	   sleep 2; \
	   exit; \
	 fi; \
	 ( cd txt; txt2tags -q -t tex -o /dev/null bool3.txt >/dev/null 2>&1 ) ; \
	 if [ $$? -ne 0 ]; \
	 then \
	   echo ">>>>>>>>>>>>>>>>>>>>>> WARNING <<<<<<<<<<<<<<<<<<<<<<<<<"; \
	   echo ">>>  txt2tags check failed -- python not installed?  <<<"; \
	   echo ">>>  (txt2tags is needed to build doc)               <<<"; \
	   echo ">>>>>>>>>>>>>>>>>>>> END WARNING <<<<<<<<<<<<<<<<<<<<<<<"; \
	   sleep 2; \
	   exit; \
	 fi; \
	 $(MAKE)  $(TEX_DIR)/CoCoALib.tex  $(TEX_DIR)/cocoalib-doc.sty; \
	 if [ \! -f "$(TEX_DIR)/CoCoALib.tex" ]; then echo; echo "ERROR!!!  FAILED TO MAKE $(TEX_DIR)/CoCoALib.tex"; echo; exit 1; fi; \
	 echo "Doing LaTeX-->PDF"; \
	 cd $(TEX_DIR); \
	 echo "[[CoCoALib:LaTeX-->PDF]]  Pass 1 of 3"; pdflatex -halt-on-error CoCoALib.tex < /dev/null > /dev/null  2>&1; \
	 if [ $$? -ne 0 ]; then echo "!!!LaTeX error!!! -- see file tex/CoCoALib.log"; exit 1; fi; \
	 echo "[[CoCoALib:LaTeX-->PDF]]  Pass 2 of 3"; pdflatex -halt-on-error CoCoALib.tex > /dev/null  2>&1; \
	 echo "[[CoCoALib:LaTeX-->PDF]]  Final pass";  pdflatex -halt-on-error CoCoALib.tex > /dev/null  2>&1; \
	 grep -F "LaTeX Warning:" CoCoALib.log; \
	 /bin/mv CoCoALib.pdf ..

# examples/index.html is really part of the documentation...
.PHONY: examples-index
examples-index:
	@cd $(COCOA_ROOT)/examples; $(MAKE) -f Makefile-indexhtml


.PHONY: texdoc
texdoc: CoCoALib.pdf


.PHONY: htmldoc
htmldoc:  $(COCOA_ROOT)/configuration/version  $(HTML_DIR)/cocoalib-doc.css  $(HTML_DIR)/index.html  $(HTML_DIR)/index-alpha.html
	@$(MAKE) -q $(DOC_HTML); \
	 if [ $$? -ne 0 ]; \
	 then \
	   command -v txt2tags >/dev/null 2>&1 ; \
	   if [ $$? -ne 0 ]; \
	   then \
	     echo ">>>>>>>>>>>>>>>>>>>>>> WARNING <<<<<<<<<<<<<<<<<<<<<<"; \
	     echo ">>>  Cannot build CoCoALib HTML documentation!    <<<"; \
	     echo ">>>  Install txt2tags from http://txt2tags.org/   <<<"; \
	     echo ">>>  (make sure you get version 3.0 or 2.7-DEV)   <<<"; \
	     echo ">>>>>>>>>>>>>>>>>>>> END WARNING <<<<<<<<<<<<<<<<<<<<"; \
	     echo; \
	     sleep 2; \
	     exit; \
	   fi; \
	 ( cd txt; txt2tags -q -t html -o /dev/null bool3.txt >/dev/null 2>&1 ) ; \
	 if [ $$? -ne 0 ]; \
	 then \
	   echo ">>>>>>>>>>>>>>>>>>>>>> WARNING <<<<<<<<<<<<<<<<<<<<<<<<<"; \
	   echo ">>>  txt2tags check failed -- python not installed?  <<<"; \
	   echo ">>>  (txt2tags is needed to build doc)               <<<"; \
	   echo ">>>>>>>>>>>>>>>>>>>> END WARNING <<<<<<<<<<<<<<<<<<<<<<<"; \
	   sleep 2; \
	   exit; \
	 fi; \
	   $(MAKE) -s $(DOC_HTML); \
	 fi


$(HTML_DIR):
	@/bin/mkdir -p $(HTML_DIR)

$(HTML_DIR)/cocoalib-doc.css:  $(HTML_DIR)  $(AUX_DIR)/cocoalib-doc.css
	@/bin/cp -p  $(AUX_DIR)/cocoalib-doc.css  $(HTML_DIR)

$(HTML_DIR)/index-alpha.html:  $(HTML_DIR)  $(AUX_DIR)/index-alpha.html
	@/bin/rm -f $(HTML_DIR)/index-alpha.html
	@sed -e "s/VERSION_NUMBER/${COCOALIB_VERSION}/"  $(AUX_DIR)/index-alpha.html  > $(HTML_DIR)/index-alpha.html

$(HTML_DIR)/index.html:  $(HTML_DIR)  $(AUX_DIR)/index.html
	@/bin/rm -f $(HTML_DIR)/index.html
	@sed -e "s/VERSION_NUMBER/${COCOALIB_VERSION}/"  $(AUX_DIR)/index.html  > $(HTML_DIR)/index.html

.PHONY: tasktables
tasktables:
	@cd $(TASK_DIR); $(MAKE) -s


.PHONY:  clean  clean-local  clean-subdirs
clean:  clean-local  clean-subdirs
	@echo "Cleaned CoCoALib/$(CWD)"

clean-local:
	@/bin/rm -f ./*~  ./.*~  ./.\#*
	@/bin/rm -rf $(TEX_DIR)

clean-subdirs:
	@cd $(TASK_DIR);  $(MAKE) -s clean
	@cd $(SRC_DIR);   /bin/rm -f  ./*~  ./.*~  ./.\#*
	@cd $(AUX_DIR);   /bin/rm -f  ./*~  ./.*~  ./.\#*  DocTeX-VERSION.t2t


.PHONY:  veryclean  veryclean-local  veryclean-subdirs
veryclean: veryclean-local  veryclean-subdirs
	@echo "Verycleaned CoCoALib/$(CWD)"

veryclean-local:  clean-local
	@/bin/rm -f $(HTML_DIR).tgz  CoCoALib.pdf  $(TASKS_OUT)  TIMESTAMP
	@/bin/rm -rf $(HTML_DIR)

veryclean-subdirs:
	@cd $(TASK_DIR);  $(MAKE) -s veryclean
	@cd $(SRC_DIR);   /bin/rm -f  ./*~  ./.*~  ./.\#*
	@cd $(AUX_DIR);   /bin/rm -f  ./*~  ./.*~  ./.\#*  DocTeX-VERSION.t2t
